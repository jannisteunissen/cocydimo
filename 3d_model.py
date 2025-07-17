#!/usr/bin/env python3
"""Cocydimo: Conducting Cylinder Discharge Model"""

import copy
import argparse
import numpy as np
import json
from time import perf_counter
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
import model_lib as mlib
from poisson_3d import m_solver as p3d

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Cocydimo: Conducting Cylinder Discharge Model')
parser.add_argument('-r_scale', type=float, default=1.2,
                    help='Scale factor compared to electrodynamic radius')
parser.add_argument('-domain_size', type=float, nargs=3,
                    default=[80e-3, 80e-3, 80e-3], help='Domain size (m)')
parser.add_argument('-coarse_grid_size', type=int, nargs=3,
                    default=[16, 16, 16],
                    help='Size of coarse grid (Nx, Ny, Nz)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo (#cells)')
parser.add_argument('-rod_r0', type=float, nargs=3,
                    default=[40e-3, 40e-3, 0e-3],
                    help='First point of rod electrode (m)')
parser.add_argument('-rod_r1', type=float, nargs=3,
                    default=[40e-3, 40e-3, 40e-3],
                    help='Second point of rod electrode (m)')
parser.add_argument('-rod_radius', type=float, default=0.75e-3,
                    help='Radius of rod electrode (m)')
parser.add_argument('-n_steps', type=int, default=100,
                    help='How many steps to simulate')
parser.add_argument('-n_initial_streamers', type=int, default=5,
                    help='How many streamers to start with')
parser.add_argument('-r_start', type=float, nargs=3,
                    help='Initial location of initial streamers (m)')
parser.add_argument('-dt', type=float, default=5e-10,
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
parser.add_argument('-siloname', type=str, default='output/simulation_3d',
                    help='Base filename for output Silo files')
parser.add_argument('-write_eps', action='store_true',
                    help='Write epsilon (of Poisson eq.) to output')
parser.add_argument('-write_time', action='store_true',
                    help='Write time that channel was added to output')
parser.add_argument('-write_rhs', action='store_true',
                    help='Write r.h.s. of Poisson eq. to output')
parser.add_argument('-rng_seed', type=int,
                    help='Seed for the random number generator')
parser.add_argument('-c_b', type=float, default=15.,
                    help='Branching coeff. - larger means less branching')
parser.add_argument('-L_b', type=float, default=2.0e-4,
                    help='Branching coeff. (m) - thin streamers branch less')
parser.add_argument('-branch_gamma', type=float, default=90.0,
                    help='Angle (degrees) that affects branching')
parser.add_argument('-refine_E', type=float, default=3e6,
                    help='Refine if E is above this value (V/m)')
parser.add_argument('-derefine_E', type=float, default=2e6,
                    help='Derefine if E is below this value (V/m)')
parser.add_argument('-derefine_nlevels', type=int, default=1,
                    help='Derefine at most this many levels')
parser.add_argument('-min_dx', type=float, default=1e-4,
                    help='Minimum allowed grid spacing (m)')
parser.add_argument('-max_dx', type=float, default=2e-3,
                    help='Maximum allowed grid spacing (m)')
parser.add_argument('-max_dx_electrode', type=float, default=8e-4,
                    help='Maximum allowed grid spacing around electrode (m)')
parser.add_argument('-memory_limit', type=float, default=8.0,
                    help='Memory limit (GB)')
parser.add_argument('-print_performance', action='store_true',
                    help='Show performance information')
parser.add_argument('-steps_per_output', type=int, default=1,
                    help='Write output every N steps')

args = parser.parse_args()

# Save settings
fname = f'{args.siloname}.json'
with open(fname, 'w') as f:
    json.dump(args.__dict__, f, indent=2)
    print(f'Wrote settings to {fname}')

model = mlib.AirStreamerModel(c0=args.c0_L_E_dx, dz0=args.dz_data)

np.random.seed(args.rng_seed)


def get_tau_branch(radius, velocity):
    """Get expected branching time for a streamer

    :param radius: radius (m)
    :param velocity: velocity (m/s)
    :returns: expected branching time
    """
    return args.c_b * radius/velocity * (1 + (args.L_b/radius)**2)


def find_orthogonal_unit_vector(y):
    """Randomly sample a unit vector orthogonal to y

    :param y: input vector
    :returns: unit vector orthogonal to y
    """
    orthvec = np.cross(y, np.random.uniform(-1., 1., 3))
    return orthvec / norm(orthvec)


p3d.set_rod_electrode(args.rod_r0, args.rod_r1, args.rod_radius)

p3d.initialize_domain(args.domain_size, args.coarse_grid_size,
                      args.box_size, args.phi_bc, args.memory_limit,
                      args.write_eps, args.write_time, args.write_rhs)

p3d.set_refinement(args.refine_E, args.derefine_E,
                   args.min_dx, args.max_dx,
                   args.max_dx_electrode, args.derefine_nlevels,
                   args.poisson_rtol)

dz = p3d.get_finest_grid_spacing()
print(f'Minimum grid spacing: {dz:.2e}')

# Compute initial solution
p3d.solve(0.0, args.poisson_rtol)
p3d.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt(args.k_eff_file).T
if args.channel_no_ionization:
    table_k_eff = np.minimum(table_k_eff, 0.0)

p3d.store_k_eff(table_fld[0], table_fld[-1], table_k_eff)

if args.r_start is not None:
    r_start = np.array(args.r_start)
else:
    Emax, r_start = p3d.get_max_field_location()
    print(f'Streamers start from {r_start}')

# Get L_E to estimate initial streamer radius
z, E, success = p3d.get_var_along_line('E_norm', r_start, [0., 0., 1.0],
                                       args.L_E_max, 2*args.L_E_max/dz)
if not success:
    raise RuntimeError('Interpolation error at r_start')
L_E = model.get_L_E(z, E, dz)

# Start with a smaller radius to approximate initial phase
radius0 = 0.5 * args.r_scale * model.get_radius(L_E)

# Start with multiple streamers in random directions
streamers = []
for i in range(args.n_initial_streamers):
    # Sample multiple velocity directions
    axis = find_orthogonal_unit_vector([0., 0., 1.0])
    angle = np.random.uniform(0., args.branch_gamma*np.pi/180.)
    rot = Rotation.from_rotvec(axis * angle)
    v_hat = rot.apply([0., 0., 1.0])
    streamers.append(mlib.Streamer(r_start - [0., 0., radius0],
                                   v_hat, radius0, 0.0))

wct_refinement = 0.0
wct_update_sigma = 0.0
wct_poisson = 0.0
wct_output = 0.0
t_start = perf_counter()

for step in range(1, args.n_steps+1):
    time = (step-1) * args.dt
    print(f'{step:4d} t = {time*1e9:.1f} ns n_streamers = {len(streamers)}')

    if args.print_performance:
        t_total = perf_counter() - t_start
        print(f' refinement: {1e2*wct_refinement/t_total:.2f}% '
              f'update_sigma: {1e2*wct_update_sigma/t_total:.2f}% '
              f'poisson: {1e2*wct_poisson/t_total:.2f}% '
              f'output: {1e2*wct_output/t_total:.2f}%')

    new_branches = []
    for s in streamers:
        tau_branch = get_tau_branch(s.R, norm(s.v))

        if np.random.exponential(tau_branch) < args.dt:
            s.is_branching = True
            s.branching_angle = np.random.uniform(0., args.branch_gamma *
                                                  np.pi/180.)
            s.branching_axis = find_orthogonal_unit_vector(s.v)

            new_branch = copy.copy(s)
            new_branch.branching_angle -= np.pi/2
            new_branch.n_steps = 0
            new_branches.append(new_branch)

    streamers = streamers + new_branches
    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:
        if s.is_branching:
            rot = Rotation.from_rotvec(s.branching_axis * s.branching_angle)
            s.v = rot.apply(s.v)
            s.is_branching = False

        # Take the field ahead of the streamer. This field will tend to bend
        # towards the background electric field.
        v_hat = s.v/norm(s.v)
        E_hat, success = p3d.get_field_vector_at(s.r + 1.5 * s.R * v_hat)

        if not success:
            print('Could not sample E_hat, removing streamer')
            s.keep = False
            continue

        E_hat = E_hat/norm(E_hat)

        # Get samples of |E| ahead of the streamer to determine L_E
        r_tip = s.r + 0.5 * s.R * s.v/norm(s.v)
        z, E, success = p3d.get_var_along_line('E_norm', r_tip, E_hat,
                                               args.L_E_max, 2*args.L_E_max/dz)

        if success:
            L_E_new = model.get_L_E(z, E, dz)
        else:
            # Use previous value
            L_E_new = s.L_E

        if L_E_new < args.L_E_min:
            print('L_E too small, removing streamer')
            s.keep = False
            continue

        if s.n_steps == 0:
            # At the first step, use only the new field
            L_E = L_E_new
        else:
            L_E = args.alpha * L_E_new + (1 - args.alpha) * s.L_E

        s.L_E = L_E
        s.sigma = model.get_sigma(L_E)
        s.v = model.get_velocity(L_E) * E_hat

        dR = min(args.r_scale * model.get_radius(L_E) - s.R,
                 norm(s.v) * args.dt)
        s.R = s.R + dR
        s.r = s.r + s.v * (args.dt - 0.99 * dR/norm(s.v))
        s.n_steps += 1

    t0 = perf_counter()
    n_add = p3d.adjust_refinement()
    t1 = perf_counter()
    wct_refinement += t1 - t0

    mlib.update_sigma(p3d.update_sigma, streamers, streamers_prev,
                      time, args.dt, args.channel_update_delay, step == 1,
                      args.channel_max_sigma)
    t0 = perf_counter()
    wct_update_sigma += t0 - t1

    p3d.solve(args.dt, args.poisson_rtol)
    t1 = perf_counter()
    wct_poisson += t1 - t0

    time += args.dt

    # Write output every N steps, with N = args.steps_per_output
    if step % args.steps_per_output == 0:
        i_output = step // args.steps_per_output
        p3d.write_solution(f'{args.siloname}_{i_output:04d}', i_output, time)

    t0 = perf_counter()
    wct_output += t0 - t1

    streamers = [s for s in streamers if s.keep]
    if len(streamers) == 0:
        raise ValueError('All streamers gone')
