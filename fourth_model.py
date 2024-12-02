#!/usr/bin/env python3

from poisson_2d import m_solver
import numpy as np
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
                    help='How many steps to simulate')
parser.add_argument('-dt_frac', type=int, default=1,
                    help='Reduce dt in model by this factor')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo')
parser.add_argument('h5file', type=str,
                    help='Input hdf5 file')
parser.add_argument('-alpha', type=float, default=1.0,
                    help='Exponential smoothing coefficient')
parser.add_argument('-siloname', type=str, default='output/run',
                    help='Base filename for output Silo files')
parser.add_argument('-plot', action='store_true',
                    help='Make plot of solution and original data')
parser.add_argument('-print_error', action='store_true',
                    help='Print difference with simulation data')
args = parser.parse_args()


mlib.E_threshold = 6e6
args.cycle = 0

with h5py.File(args.h5file, 'r') as h5f:
    z = np.array(h5f['z_line'])  # Could in the future differ from z_grid
    sigma_z = np.array(h5f[f'run_{args.run}/sigmaz_{args.cycle}'])
    phi_bc = np.array(h5f[f'run_{args.run}/voltage'][args.cycle])
    times = np.array(h5f[f'run_{args.run}/time'])
    radius = np.array(h5f[f'run_{args.run}/radius'][args.cycle])
    velocity = np.array(h5f[f'run_{args.run}/velocity'][args.cycle])
    field_rod_r1 = np.array(h5f['field_rod_r1'][args.run])
    field_rod_radius = np.array(h5f['field_rod_radius'][args.run])
    field_amplitude = np.array(h5f['field_amplitude'][args.run])

n_cycles = len(times)
Nz = len(z)
dz = z[1] - z[0]
dt = times[1] - times[0]
dt_model = dt / args.dt_frac
time = times[args.cycle]

domain_ncells = [Nz, Nz]
domain_size = [Nz * dz, Nz * dz]
print(domain_ncells, domain_size)

m_solver.set_rod_electrode([0.0, 0.0], field_rod_r1*domain_size,
                           field_rod_radius)
m_solver.initialize(domain_size, domain_ncells, args.box_size, phi_bc)

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}{args.run}_{0:04d}', 0, time)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt('k_eff_air.txt').T
m_solver.store_k_eff(table_fld[0], table_fld[-1], table_k_eff)

# Locate Emax
Emax, r_Emax = m_solver.get_max_field_location()

# z-coordinate lies at the center of the streamer head
radius0 = args.r_scale * radius
z0 = r_Emax[1] - radius0

_, phi = m_solver.get_line_potential([0., z[0]], [0., z[-1]], len(z))

streamers = [mlib.Streamer([0.0, z0], [0., 0.], radius0, sigma_z.max())]

# For plots
phi_z_pred = np.zeros((args.nsteps+1, Nz))
phi_z_pred[0] = phi

sigma_head = np.zeros((args.nsteps+1))
sigma_head[0] = streamers[0].sigma

z_head = np.zeros((args.nsteps+1))
z_head[0] = streamers[0].r[1] + streamers[0].R

L_E_all = np.zeros((args.nsteps+1))
L_E_all[0] = mlib.get_high_field_length(z, np.abs(np.gradient(phi, dz)))


for step in range(1, args.nsteps+1):
    streamers_prev = copy.deepcopy(streamers)

    # Get input features for model
    L_E_new = mlib.get_high_field_length(z, np.abs(np.gradient(phi, dz)))

    if step == 1:
        L_E = L_E_new
    else:
        L_E = args.alpha * L_E_new + (1 - args.alpha) * L_E

    for s in streamers:
        old_sigma = s.sigma
        old_v = s.v.copy()
        tmp = mlib.get_sigma_v2(L_E, field_amplitude)
        s.sigma = min(s.sigma + mlib.max_dt_sigma * dt_model, tmp)

        s.R = args.r_scale * mlib.get_radius_v2(L_E, field_amplitude)
        s.v[1] = mlib.get_velocity_v2(L_E)
        s.r = s.r + 0.5 * (old_v + s.v) * dt_model

    sigma_head[step] = streamers[0].sigma
    z_head[step] = streamers[0].r[1] + streamers[0].R
    L_E_all[step] = L_E

    mlib.update_sigma(m_solver.update_sigma, streamers, streamers_prev,
                      time, dt_model, 1e-9, False)
    m_solver.solve(dt_model)

    if s.sigma < mlib.sigma_min:
        print('Streamer removed due to too low sigma')
        s.keep = False

    time += dt_model

    _, phi = m_solver.get_line_potential([0., z[0]], [0., z[-1]], len(z))
    phi_z_pred[step] = phi

    m_solver.write_solution(f'{args.siloname}{args.run}_{step:04d}',
                            step, time)

    streamers = [s for s in streamers if s.keep]

    if len(streamers) == 0:
        raise ValueError('All streamers gone')


if args.print_error:
    ix = 20

    with h5py.File(args.h5file, 'r') as h5f:
        phiz_sim = np.array(h5f[f'run_{args.run}/phiz_{ix}'])
        sigmaz_sim = np.array(h5f[f'run_{args.run}/sigmaz_{ix}'])
        iEmax_sim = np.argmax(np.abs(np.gradient(phiz_sim)))
        iEmax_model = np.argmax(np.abs(np.gradient(phi_z_pred[ix])))
        print(iEmax_model - iEmax_sim)


if args.plot:
    fig, ax = plt.subplots(4, sharex=True, layout='constrained')

    nsteps_data = min(args.nsteps//args.dt_frac+1,
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

    ax[0].plot(z_head, sigma_head, marker='.')

    for step in range(0, args.nsteps+1, args.dt_frac):
        ax[1].plot(z, phi_z_pred[step], ls='--')
        ax[2].plot(z, -np.gradient(phi_z_pred[step], dz), ls='--')

    ax[0].set_ylabel('sigma')
    ax[1].set_ylabel('phi (V)')
    ax[2].set_ylabel('max(E)')
    ax[3].plot(z_head, L_E_all, label='model')

    # Compare with data
    L_E_data = [mlib.get_high_field_length(z, np.abs(np.gradient(phi, dz)))
                for phi in phiz_data]
    i_head_data = [np.argmax(np.abs(np.gradient(phi))) for phi in phiz_data]
    z_head_data = [z[i] for i in i_head_data]

    ax[3].plot(z_head_data, L_E_data, label='data')
    ax[3].set_ylabel('length(E)')
    ax[3].legend()
    ax[-1].set_xlabel('z (m)')

    plt.show()
