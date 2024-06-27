#!/usr/bin/env python3

from poisson import m_solver
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
parser.add_argument('-siloname', type=str, default='output/simulation',
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
    xs = (x - x_head)
    cond = [xs < 0, xs < radius, xs >= radius]
    values = [sigma_head * np.clip(1 + (x-x_head) *
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
Nr, Nz = rhs.shape

dz = z_grid[1] - z_grid[0]
dr = r_grid[1] - r_grid[0]
dt = times[1] - times[0]
dt_model = dt / args.dt_frac

domain_size = [Nr * dr, Nz * dz]

sigma_z_pred = np.zeros((args.nsteps+1, Nz))
phi_z_pred = np.zeros((args.nsteps+1, Nz))
sigma_head = np.zeros((args.nsteps+1))
z_head = np.zeros((args.nsteps+1))
L_E = np.zeros((args.nsteps+1))

m_solver.set_rod_electrode([0.0, 0.0], [0.0, 0.15*domain_size[1]], 0.5e-3)
m_solver.initialize(domain_size, [Nr, Nz], args.box_size, phi_bc)
m_solver.set_rhs_and_sigma(rhs, sigma)

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}_{0:04d}')
z_phi, phi = m_solver.get_line_potential(z[0], z[-1], len(z))

phi_z_pred[0] = phi
sigma_z_pred[0] = sigma_z
# i_head = np.argmax(np.abs(np.gradient(phi)))
i_head = np.argmax(np.abs(sigma_z_pred[0]))
sigma_head[0] = sigma_z_pred[0].max()
z_head[0] = z[i_head]


for step in range(1, args.nsteps+1):
    # Get input features for model
    L_E[step-1] = get_high_field_length(phi_z_pred[step-1], dz)

    dsigma_dt = get_dsigma_dt(L_E[step-1])
    v = get_velocity(sigma_head[step-1])
    R = get_radius(sigma_head[step-1])
    # slope = dsigma_dt/v

    sigma_head[step] = sigma_head[step-1] + dsigma_dt * dt_model
    z_head[step] = z_head[step-1] + v*dt_model

    # Determine slope between new z_head and z_head - 2 * R
    i_z = np.argmax(z > z_head[step] - 2 * R)
    slope = (sigma_head[step] - sigma_z_pred[step-1][i_z]) / (2 * R)

    new_sigma_z = sigma_func(z, z_head[step], sigma_head[step], R, slope)
    sigma_z_pred[step] = np.where(z > z_head[step] - 2 * R,
                                  new_sigma_z, sigma_z_pred[step-1])

    w_r = get_radial_weights(R, 20, dz)
    sigma_z_diff = sigma_z_pred[step] - sigma_z_pred[step-1]

    m_solver.update_sigma(z[0], z[-1], sigma_z_diff, R, w_r)
    m_solver.solve(dt_model)

    z_phi, phi = m_solver.get_line_potential(z[0], z[-1], len(z))
    phi_z_pred[step] = phi

    m_solver.write_solution(f'{args.siloname}_{step:04d}')


if args.plot:
    fig, ax = plt.subplots(4, sharex=True, layout='constrained')

    nsteps_data = min(args.nsteps//args.dt_frac + 1,
                      n_cycles - args.cycle)

    with h5py.File(args.h5file, 'r') as h5f:
        sigmaz_data = []
        phiz_data = []

        for i in range(args.cycle, args.cycle+nsteps_data):
            tmp = np.array(h5f[f'run_{args.run}/sigmaz_{i}'])
            ax[0].plot(z, tmp, c='gray')
            sigmaz_data.append(tmp)

            tmp = np.array(h5f[f'run_{args.run}/phiz_{i}'])
            ax[1].plot(z, tmp, c='gray')
            ax[2].plot(z, -np.gradient(tmp, dz), c='gray')
            phiz_data.append(tmp)

    for step in range(0, args.nsteps+1, args.dt_frac):
        ax[0].plot(z, sigma_z_pred[step], ls='--')
        ax[1].plot(z, phi_z_pred[step], ls='--')
        ax[2].plot(z, -np.gradient(phi_z_pred[step], dz), ls='--')

    ax[0].set_ylabel('sigma')
    ax[1].set_ylabel('phi (V)')
    ax[2].set_ylabel('max(E)')

    ax[3].plot(z_head[:-1], L_E[:-1], label='model')

    # Compare with data
    L_E_data = [get_high_field_length(phi, dz) for phi in phiz_data]
    i_head_data = [np.argmax(sigma) for sigma in sigmaz_data]
    z_head_data = [z[i] for i in i_head_data]

    ax[3].plot(z_head_data, L_E_data, label='data')
    ax[3].set_ylabel('length(E)')
    ax[3].legend()
    ax[-1].set_xlabel('z (m)')

    plt.show()
