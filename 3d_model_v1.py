#!/usr/bin/env python3

from poisson_3d import m_solver_3d
import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Reduced discharge model')
parser.add_argument('-run', type=int, default=0,
                    help='Index of run to simulate')
parser.add_argument('-cycle', type=int, default=10,
                    help='Start from this cycle')
parser.add_argument('-r_scale', type=float, default=1.5,
                    help='Scale factor for the radius')
parser.add_argument('-nsteps', type=int, default=2,
                    help='How many steps to simulate')
parser.add_argument('-dt_frac', type=int, default=1,
                    help='Reduce dt in model by this factor')
parser.add_argument('-box_size', type=int, default=32,
                    help='Size of boxes in afivo')
parser.add_argument('h5file', type=str,
                    help='Input hdf5 file')
parser.add_argument('-siloname', type=str, default='output/simulation_3d',
                    help='Base filename for output Silo files')
parser.add_argument('-plot', action='store_true',
                    help='Make plot of solution and original data')
args = parser.parse_args()

E_threshold = 5e6


def get_high_field_length(phi, dz):
    E = np.abs(np.gradient(phi, dz))
    i_max = np.argmax(E)
    d_i = np.argmax(E[i_max:] < E_threshold)
    return d_i * dz


def sigma_func(x, x_head, sigma_head, radius, slope_behind):
    # Scale coordinate
    xs = (x - x_head + radius)
    cond = [xs < 0, xs < radius, xs >= radius]
    values = [sigma_head * np.clip(1 + (x-x_head+radius) *
                                   slope_behind/sigma_head, 0., 1.),
              sigma_head * (1 - (xs/radius)**3),
              np.zeros_like(xs)]
    return np.select(cond, values)


def get_radius(L_E):
    return args.r_scale * \
        np.sqrt(8.441501306469196e-08 + 0.6097078090287649 * L_E)


def get_velocity(sigma_max):
    return 307269.6823523594 + 1412173247560.5532 * sigma_max


def get_dsigma_dt(L_E):
    return -34.43873725692916 + 262462782.20047116 * L_E**2


def get_radial_weights(r_max, n_points, dz):
    r = np.linspace(0., r_max, n_points)
    dr = r[1] - r[0]
    w = np.maximum(0.0, 1 - 3*(r/r_max)**2 + 2*(r/r_max)**3)
    w_sum = np.sum(2 * np.pi * r * dr * w)
    return w/w_sum


def map_cyl_to_3d(cyl_data, r_grid):
    Nx = Ny = cyl_data.shape[0]
    Nz = cyl_data.shape[1]
    xyz_data = np.zeros([Nx//1, Ny//1, Nz//1])

    # # Assume even-sized array, cut half off
    # r_grid = r_grid[:r_grid.size//1]
    # cyl_data = cyl_data[:cyl_data.size//1]
    # i_center = r_grid.size//1
    # r_center = 0.5 * (r_grid[i_center-1] + r_grid[i_center])
    # r_xy = np.meshgrid(np.abs(r_grid-r_center), np.abs(r_grid-r_center))
    # rr = (r_xy[0]**2 + r_xy[1]**2)**0.5

    # for iz in range(Nz):
    #     xyz_data[:, :, iz] = np.interp(rr, cyl_data, r_grid)
    return xyz_data


with h5py.File(args.h5file, 'r') as h5f:
    z = np.array(h5f['z_line'])  # Could in the future differ from z_grid
    sigma_z = np.array(h5f[f'run_{args.run}/sigmaz_{args.cycle}'])
    sigma = np.array(h5f[f'run_{args.run}/sigma_{args.cycle}'])
    rhs = np.array(h5f[f'run_{args.run}/rhs_{args.cycle}'])
    phi_bc = np.array(h5f[f'run_{args.run}/voltage'][args.cycle])
    times = np.array(h5f[f'run_{args.run}/time'])
    r_grid = np.array(h5f[f'run_{args.run}/r_grid'])
    z_grid = np.array(h5f[f'run_{args.run}/z_grid'])

n_cycles = len(times)
Nx, Ny, Nz = 256, 256, 256
domain_size = [20e-3, 20e-3, 20e-3]
Lx = domain_size[0]

dt = times[1] - times[0]
dt_model = dt / args.dt_frac

rhs_3d = np.zeros((Nx, Ny, Nz))
sigma_3d = np.zeros((Nx, Ny, Nz))

m_solver_3d.set_rod_electrode([0.5*Lx, 0.5*Lx, 0.0],
                              [0.5*Lx, 0.5*Lx, 0.15*domain_size[2]], 0.5e-3)
m_solver_3d.initialize(domain_size, [Nx, Ny, Nz], args.box_size, phi_bc)
m_solver_3d.set_rhs_and_sigma(rhs_3d, sigma_3d)

# Compute initial solution
m_solver_3d.solve(0.0)
m_solver_3d.write_solution(f'{args.siloname}_{0:04d}')
z_phi, phi = m_solver_3d.get_line_potential([0.5*Lx, 0.5*Lx, z[0]],
                                            [0.5*Lx, 0.5*Lx, z[-1]], len(z))

n_streamers = 2
sigma_heads = np.zeros(2)
sigma_heads[:] = [sigma_z.max(), sigma_z.max()]
sigma_heads_prev = sigma_heads.copy()

i_head = np.argmax(np.abs(np.gradient(phi)))
r_heads = np.zeros((n_streamers, 3))
r_heads[:, :] = [[0.5*Lx, 0.5*Lx, 3.5e-3], [0.55*Lx, 0.5*Lx, 3.5e-3]]
r_heads_prev = r_heads.copy()

v_heads = np.zeros((n_streamers, 3))
v_heads[:, :] = [[0., 0., 1e6], [0., 0., 1e6]]
v_heads_prev = v_heads.copy()

radius_heads = np.zeros(2)
radius_heads_prev = radius_heads.copy()

for step in range(1, args.nsteps+1):

    for i in range(n_streamers):
        # Get location ahead of streamer
        r_ahead = r_heads_prev[i] + 1.5*radius_heads_prev[i] * \
            v_heads_prev[i]/np.linalg.norm(v_heads_prev[i])

        # Get electric field unit vector
        E_vec = m_solver_3d.get_electric_field(r_ahead)
        E_hat = E_vec / np.linalg.norm(E_vec)

        v_heads[i] = 1e6 * E_hat
        r_heads[i] = r_heads_prev[i] + v_heads[i] * dt_model
        radius_heads[i] = 0.5e-3
        sigma_heads[i] = sigma_heads_prev[i]

    m_solver_3d.update_sigma(r_heads_prev, r_heads,
                             sigma_heads_prev, sigma_heads,
                             radius_heads_prev, radius_heads, step == 1)
    m_solver_3d.solve(dt_model)
    m_solver_3d.write_solution(f'{args.siloname}_{step:04d}')

    radius_heads_prev = radius_heads.copy()
    r_heads_prev = r_heads.copy()
    v_heads_prev = v_heads.copy()
    sigma_heads_prev = sigma_heads.copy()
