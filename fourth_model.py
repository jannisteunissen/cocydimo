#!/usr/bin/env python3

from poisson_2d import m_solver
import numpy as np
from numpy.linalg import norm
import h5py
import copy
import matplotlib.pyplot as plt
import argparse
import model_lib as mlib

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Reduced discharge model')
parser.add_argument('-run', type=int, default=0,
                    help='Index of run to simulate')
parser.add_argument('-r_scale', type=float, default=1.2,
                    help='Scale factor for the radius')
parser.add_argument('-nsteps', type=int, default=2,
                    help='How many steps from the data to simulate')
parser.add_argument('-dt_factor', type=float, default=1.0,
                    help='Multiply dt in model by this factor')
parser.add_argument('-grid_size', type=int, nargs=2,
                    help='Size of computational grid (Nr, Nz)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo')
parser.add_argument('h5file', type=str,
                    help='Input hdf5 file')
parser.add_argument('-alpha', type=float, default=1.0,
                    help='Exponential smoothing coefficient')
parser.add_argument('-c0_L_E_dx', type=float,
                    help='Correction factor for L_E w.r.t. grid spacing')
parser.add_argument('-siloname', type=str, default='output/run',
                    help='Base filename for output Silo files')
parser.add_argument('-plot', action='store_true',
                    help='Make plot of solution and original data')
args = parser.parse_args()

args.cycle = 0

with h5py.File(args.h5file, 'r') as h5f:
    z_data = np.array(h5f['z_line'])  # Could in the future differ from z_grid
    sigma_z = np.array(h5f[f'run_{args.run}/sigmaz_{args.cycle}'])
    phi_bc = np.array(h5f[f'run_{args.run}/voltage'][args.cycle])
    times = np.array(h5f[f'run_{args.run}/time'])
    radius = np.array(h5f[f'run_{args.run}/radius'][args.cycle])
    velocity = np.array(h5f[f'run_{args.run}/velocity'][args.cycle])
    field_rod_r1 = np.array(h5f['field_rod_r1'][args.run])
    field_rod_radius = np.array(h5f['field_rod_radius'][args.run])
    field_amplitude = np.array(h5f['field_amplitude'][args.run])

n_cycles = len(times)
Nz_data = len(z_data)
dz_data = z_data[1] - z_data[0]
domain_size = np.array([Nz_data * dz_data, Nz_data * dz_data])

if args.grid_size is None:
    args.grid_size = [Nz_data, Nz_data]

dz = domain_size[1] / args.grid_size[1]
z = np.linspace(0.5*dz, domain_size[1]-0.5*dz, args.grid_size[1])
dt = times[1] - times[0]
time = times[args.cycle]
end_time = times[args.cycle + args.nsteps]
n_steps_model = np.ceil(args.nsteps/args.dt_factor).astype(int)
dt_model = (end_time - time)/n_steps_model

print(f'Time step in data:  {dt:.3e} s')
print(f'Time step in model: {dt_model:.3e} s')
print(f'End time:           {end_time:.3e} s')

m_solver.set_rod_electrode([0.0, 0.0], field_rod_r1*domain_size,
                           field_rod_radius)
m_solver.initialize_domain(domain_size, args.grid_size, args.box_size, phi_bc)

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}{args.run}_{0:04d}', 0, time)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt('k_eff_air.txt').T
m_solver.store_k_eff(table_fld[0], table_fld[-1], table_k_eff)

# Get phi and locate Emax
_, phi = m_solver.get_var_along_line('phi', [0., z[0]], [0., 1.],
                                     z[-1] - z[0], len(z))

Emax, r_Emax = m_solver.get_max_field_location()

# Estimate initial streamer radius
L_E = mlib.get_high_field_length(z, np.abs(np.gradient(phi, dz)),
                                 args.c0_L_E_dx, dz_data, dz)
radius0 = 0.5 * args.r_scale * mlib.get_radius_v3(L_E)

# z-coordinate lies at the center of the streamer head
z0 = r_Emax[1] - radius0

streamers = [mlib.Streamer([0.0, z0], [0., 0.], radius0, sigma_z.max())]

# For plots
phi_z_pred = np.zeros((n_steps_model+1, args.grid_size[1]))
phi_z_pred[0] = phi

sigma_head = np.zeros((n_steps_model+1))
sigma_head[0] = streamers[0].sigma

z_head = np.zeros((n_steps_model+1))
z_head[0] = streamers[0].r[1] + streamers[0].R

L_E_all = np.zeros((n_steps_model+1))
L_E_all[0] = mlib.get_high_field_length(z, np.abs(np.gradient(phi, dz)),
                                        args.c0_L_E_dx, dz_data, dz)

for step in range(1, n_steps_model+1):
    streamers_prev = copy.deepcopy(streamers)

    # Get input features for model
    L_E_new = mlib.get_high_field_length(z, np.abs(np.gradient(phi, dz)),
                                         args.c0_L_E_dx, dz_data, dz)

    if step == 1:
        L_E = L_E_new
    else:
        L_E = args.alpha * L_E_new + (1 - args.alpha) * L_E

    for s in streamers:
        old_sigma = s.sigma
        old_v = s.v.copy()
        s.sigma = mlib.get_sigma_v3(L_E)
        s.v[1] = mlib.get_velocity_v3(L_E)

        dR = min(args.r_scale * mlib.get_radius_v3(L_E) - s.R,
                 norm(s.v) * dt_model)
        s.R = s.R + dR

        # Subtract dR to so that charge layer moves with velocity v
        s.r = s.r + s.v * dt_model - [0., 0.99*dR]

    sigma_head[step] = streamers[0].sigma
    z_head[step] = streamers[0].r[1] + streamers[0].R
    L_E_all[step] = L_E

    mlib.update_sigma(m_solver.update_sigma, streamers, streamers_prev,
                      time, dt_model, 1e-9, False)
    m_solver.solve(dt_model)

    time += dt_model

    _, phi = m_solver.get_var_along_line('phi', [0., z[0]], [0., 1.],
                                         z[-1] - z[0], len(z))
    phi_z_pred[step] = phi

    m_solver.write_solution(f'{args.siloname}{args.run}_{step:04d}',
                            step, time)

    streamers = [s for s in streamers if s.keep]

    if len(streamers) == 0:
        raise ValueError('All streamers gone')


if args.plot:
    fig, ax = plt.subplots(4, sharex=True, layout='constrained')

    with h5py.File(args.h5file, 'r') as h5f:
        i = args.cycle + args.nsteps
        sigmaz_data = np.array(h5f[f'run_{args.run}/sigmaz_{i}'])
        ax[0].plot(z_data, sigmaz_data, c='gray')

        phiz_data = np.array(h5f[f'run_{args.run}/phiz_{i}'])
        ax[1].plot(z_data, phiz_data, c='gray')
        ax[2].plot(z_data, -np.gradient(phiz_data, dz_data), c='gray')

    ax[0].plot(z_head, sigma_head, marker='.')

    ax[1].plot(z, phi_z_pred[n_steps_model], ls='--')
    ax[2].plot(z, -np.gradient(phi_z_pred[n_steps_model], dz), ls='--')

    ax[0].set_ylabel('sigma')
    ax[1].set_ylabel('phi (V)')
    ax[2].set_ylabel('E (V/m)')
    ax[3].plot(z_head, L_E_all, label='model')

    ax[3].set_ylabel('length(E)')
    ax[3].legend()
    ax[-1].set_xlabel('z (m)')

    plt.show()
