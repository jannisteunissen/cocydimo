#!/usr/bin/env python3
"""Cocydimo: Conducting Cylinder Discharge Model"""

import copy
import argparse
import numpy as np
import json
import os
from numpy.linalg import norm
import model_lib as mlib
from poisson_2d import m_solver as p2d

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Cocydimo: Conducting Cylinder Discharge Model')
parser.add_argument('-r_scale', type=float, default=1.2,
                    help='Scale factor compared to electrodynamic radius')
parser.add_argument('-domain_size', type=float, nargs=2,
                    default=[30e-3, 30e-3], help='Domain size (m)')
parser.add_argument('-coarse_grid_size', type=int, nargs=2,
                    default=[16, 16],
                    help='Size of coarse grid (Nr, Nz)')
parser.add_argument('-grid_size', type=int, nargs=2, default=[256, 256],
                    help='Size of computational grid (Nr, Nz)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo (#cells)')
parser.add_argument('-rod_r0', type=float, nargs=2,
                    default=[0.0e-3, 0.0e-3],
                    help='First point of rod electrode (m)')
parser.add_argument('-rod_r1', type=float, nargs=2,
                    default=[0.0, 5.0e-3],
                    help='Second point of rod electrode (m)')
parser.add_argument('-rod_radius', type=float, default=0.75e-3,
                    help='Radius of rod electrode (m)')
parser.add_argument('-n_steps', type=int, default=10,
                    help='How many steps to simulate')
parser.add_argument('-dt', type=float, default=2.5e-10,
                    help='Time step (s)')
parser.add_argument('-dz_data', type=float, default=30e-3/256,
                    help='Grid spacing used to obtain L_E from dataset (m)')
parser.add_argument('-phi_bc', type=float, default=-4e4,
                    help='Applied potential (V)')
parser.add_argument('-alpha', type=float, default=0.5,
                    help='Exponential smoothing coefficient')
parser.add_argument('-channel_update_delay', type=float, default=1e-9,
                    help='Delay for first updating channel conductivity (s)')
parser.add_argument('-channel_no_ionization', action='store_true',
                    help='Do not increase channel conductivity, which can be '
                    'problematic near domain boundaries)')
parser.add_argument('-channel_max_sigma', type=float, default=5.0,
                    help='Limit growth of volume conductivity to this value '
                    'to prevent issues at domain boundaries [A/(m V)]')
parser.add_argument('-L_E_max', type=float, default=5e-3,
                    help='Maximum value of L_E (m)')
parser.add_argument('-L_E_min', type=float, default=1e-4,
                    help='Minimum value of L_E (m)')
parser.add_argument('-c0_L_E_dx', type=float, default=0.75,
                    help='Correction factor for L_E w.r.t. data grid spacing')
parser.add_argument('-k_eff_file', type=str, default='data/k_eff_air.txt',
                    help='File with k_eff (1/s) vs electric field (V/m)')
parser.add_argument('-poisson_rtol', type=float, default=1e-5,
                    help='Relative tolerance for Poisson solver')
parser.add_argument('-siloname', type=str, default='output/simulation_2d',
                    help='Base filename for output Silo files')
parser.add_argument('-write_eps', action='store_true',
                    help='Write epsilon (of Poisson eq.) to output')
parser.add_argument('-write_time', action='store_true',
                    help='Write time that channel was added to output')
parser.add_argument('-write_rhs', action='store_true',
                    help='Write r.h.s. of Poisson eq. to output')
parser.add_argument('-rng_seed', type=int,
                    help='Seed for the random number generator')
parser.add_argument('-memory_limit', type=float, default=0.1,
                    help='Memory limit (GB)')
parser.add_argument('-steps_per_output', type=int, default=1,
                    help='Write output every N steps')

args = parser.parse_args()

# Make sure output folder exists
os.makedirs(os.path.dirname(args.siloname), exist_ok=True)

# Save settings
fname = f'{args.siloname}.json'
with open(fname, 'w') as f:
    json.dump(args.__dict__, f, indent=2)
    print(f'Wrote settings to {fname}')

model = mlib.AirStreamerModel(c0=args.c0_L_E_dx, dz0=args.dz_data)

np.random.seed(args.rng_seed)

p2d.set_rod_electrode(args.rod_r0, args.rod_r1, args.rod_radius)
p2d.initialize_domain(args.domain_size, args.coarse_grid_size,
                      args.box_size, args.phi_bc, args.memory_limit,
                      args.write_eps, args.write_time, args.write_rhs)
p2d.use_uniform_grid(args.grid_size)

dz = p2d.get_finest_grid_spacing()
print(f'Minimum grid spacing: {dz:.2e}')

# Compute initial solution
p2d.solve(0.0, args.poisson_rtol)
p2d.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt(args.k_eff_file).T
if args.channel_no_ionization:
    table_k_eff = np.minimum(table_k_eff, 0.0)

p2d.store_k_eff(table_fld[0], table_fld[-1], table_k_eff)

# Get L_E to estimate initial streamer radius
Emax, r_Emax = p2d.get_max_field_location()
z, E, success = p2d.get_var_along_line('E_norm', [0.0, r_Emax[1]], [0., 1.0],
                                       args.L_E_max, 2*args.L_E_max/dz)
if not success:
    raise RuntimeError('Interpolation error')

L_E = model.get_L_E(z, E, dz)

# Start with a smaller radius to approximate initial phase
radius0 = 0.5 * args.r_scale * model.get_radius(L_E)

# Single streamer in the z-direction
streamers = [mlib.Streamer([0.0, r_Emax[1] - radius0],
                           [0., 1.0], radius0, 0.0)]

for step in range(1, args.n_steps+1):
    time = (step-1) * args.dt
    print(f'{step:4d} t = {time*1e9:.1f} ns')

    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:

        # Get samples of |E| ahead of the streamer to determine L_E
        r_tip = s.r + 0.5 * s.R * s.v/norm(s.v)
        z, E, success = p2d.get_var_along_line('E_norm', r_tip, s.v,
                                               args.L_E_max, 2*args.L_E_max/dz)
        if not success:
            print('Could not sample L_E, removing streamer')
            s.keep = False
            continue

        L_E_new = model.get_L_E(z, E, dz)

        if L_E_new < args.L_E_min:
            print(f'L_E too small {L_E_new:.2e}, removing streamer')
            s.keep = False
            continue

        if step == 1:
            L_E = L_E_new
        else:
            L_E = args.alpha * L_E_new + (1 - args.alpha) * L_E

        # Propagation in +z direction
        E_hat = np.array([0.0, 1.0])

        s.sigma = model.get_sigma(L_E)
        s.v = model.get_velocity(L_E) * E_hat

        dR = min(args.r_scale * model.get_radius(L_E) - s.R,
                 norm(s.v) * args.dt)
        s.R = s.R + dR
        s.r = s.r + s.v * (args.dt - 0.99 * dR/norm(s.v))

    mlib.update_sigma(p2d.update_sigma, streamers, streamers_prev,
                      time, args.dt, args.channel_update_delay, step == 1,
                      args.channel_max_sigma)
    p2d.solve(args.dt, args.poisson_rtol)

    time += args.dt

    # Write output every N steps, with N = args.steps_per_output
    if step % args.steps_per_output == 0:
        i_output = step // args.steps_per_output
        p2d.write_solution(f'{args.siloname}_{i_output:04d}', i_output, time)

    streamers = [s for s in streamers if s.keep]
    if len(streamers) == 0:
        raise ValueError('All streamers gone')
