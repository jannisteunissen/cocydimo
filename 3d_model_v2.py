#!/usr/bin/env python3

from poisson_3d import m_solver
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import copy
import argparse
from scipy.spatial.transform import Rotation
import model_lib as mlib


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Reduced discharge model')
parser.add_argument('-run', type=int, default=0,
                    help='Index of run to simulate')
parser.add_argument('-r_scale', type=float, default=1.2,
                    help='Scale factor for the radius')
parser.add_argument('-domain_size', type=float, nargs=3,
                    default=[30e-3, 30e-3, 30e-3], help='Domain size')
parser.add_argument('-rod_r0', type=float, nargs=3,
                    default=[15e-3, 15e-3, 0e-3],
                    help='First point of rod electrode')
parser.add_argument('-rod_r1', type=float, nargs=3,
                    default=[15e-3, 15e-3, 5e-3],
                    help='Second point of rod electrode')
parser.add_argument('-rod_radius', type=float, default=0.75e-3,
                    help='Radius of rod electrode')
parser.add_argument('-nsteps', type=int, default=2,
                    help='How many steps to simulate')
parser.add_argument('-grid_size', type=int, nargs=3, default=[256, 256, 256],
                    help='Size of computational grid')
parser.add_argument('-dt', type=float, default=2.5e-10,
                    help='Time step (s)')
parser.add_argument('-dz_data', type=float, default=30e-3/256,
                    help='Grid spacing used to obtain L_E from simulations')
parser.add_argument('-phi_bc', type=float, default=-4e4,
                    help='Applied potential (V)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo')
parser.add_argument('-alpha', type=float, default=0.5,
                    help='Exponential smoothing coefficient')
parser.add_argument('-L_E_max', type=float, default=5e-3,
                    help='Maximum value of L_E')
parser.add_argument('-L_E_min', type=float, default=1e-4,
                    help='Minimum value of L_E')
parser.add_argument('-c0_L_E_dx', type=float,
                    help='Correction factor for L_E w.r.t. data grid spacing')
parser.add_argument('-siloname', type=str, default='output/simulation_3d',
                    help='Base filename for output Silo files')
parser.add_argument('-rng_seed', type=int, default=8912734,
                    help='Seed for the random number generator')
parser.add_argument('-time_between_branching', type=float, default=2e-9,
                    help='Minimum time between branching events of a channel')
parser.add_argument('-c_b', type=float, default=10.,
                    help='Branching coeff. (larger means less branching)')
parser.add_argument('-L_b', type=float, default=1.0e-4,
                    help='Branching coeff. (thin streamers branch less')

args = parser.parse_args()
np.random.seed(args.rng_seed)

m_solver.set_rod_electrode(args.rod_r0, args.rod_r1, args.rod_radius)
m_solver.initialize_domain(args.domain_size, args.grid_size,
                           args.box_size, args.phi_bc, 8.0)
dz = args.domain_size[2]/args.grid_size[2]

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt('k_eff_air.txt').T
m_solver.store_k_eff(table_fld[0], table_fld[-1], table_k_eff)

# Locate Emax
Emax, r_Emax = m_solver.get_max_field_location()
r_Emax[0:2] = args.rod_r1[0:2]

# Estimate initial streamer radius
z, E = m_solver.get_var_along_line('E_norm', r_Emax, [0., 0., 1.0],
                                   args.L_E_max, 2*args.L_E_max/dz)
L_E = mlib.get_high_field_length(z, E, args.c0_L_E_dx, args.dz_data, dz)
radius0 = 0.5 * args.r_scale * mlib.get_radius_v3(L_E)

streamers = [mlib.Streamer(r_Emax - [0., 0., radius0],
                           [0., 0., 1.0], radius0, 0.0)]


def get_tau_branch(R, v):
    tau_branch = args.c_b * R/v * (1 + (args.L_b/R)**2)
    return tau_branch


def find_orthogonal_unit_vector(y):
    # Find arbitrary vector not parallel to y
    vec1 = np.cross(y, np.random.uniform(-1., 1., 3))
    return vec1 / np.linalg.norm(vec1)


L_E_data = np.zeros(args.nsteps)

for step in range(1, args.nsteps+1):
    time = (step-1) * args.dt
    print(f'{step:4d} t = {time*1e9:.1f} ns n_streamers = {len(streamers)}')

    new_branches = []
    for s in streamers:
        tau_branch = get_tau_branch(s.R, norm(s.v))

        if np.random.exponential(tau_branch) < args.dt and \
           s.time_previous_branch < time - args.time_between_branching:
            s.is_branching = True
            s.time_previous_branch = time
            s.branching_angle = np.random.uniform(0., 90.)/180. * np.pi
            s.branching_axis = find_orthogonal_unit_vector(s.v)

            new_branch = copy.copy(s)
            new_branch.branching_angle -= np.pi/2
            new_branches.append(new_branch)

    streamers = streamers + new_branches
    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:
        if s.is_branching:
            rot = Rotation.from_rotvec(s.branching_axis * s.branching_angle)
            s.v = rot.apply(s.v)
            s.is_branching = False

        # Get samples of sigma in the vicinity of the streamer
        _, sigma = m_solver.get_var_along_line('sigma', s.r, s.v, 2*s.R, 5)

        # Get samples of |E| ahead of the streamer to determine L_E
        r_tip = s.r + 0.5 * s.R * s.v/norm(s.v)
        z, E = m_solver.get_var_along_line('E_norm', r_tip, s.v,
                                           args.L_E_max, 2*args.L_E_max/dz)
        L_E_new = mlib.get_high_field_length(z, E, args.c0_L_E_dx,
                                             args.dz_data, dz)

        if L_E_new < args.L_E_min:
            print('L_E too small, removing streamer')
            s.keep = False
            continue

        if step == 1:
            L_E = L_E_new
        else:
            L_E = args.alpha * L_E_new + (1 - args.alpha) * L_E

        L_E_data[step-1] = L_E

        # Take the field ahead of the streamer. This field will tend to bend
        # towards the background electric field.
        v_hat = s.v/norm(s.v)
        E_hat = m_solver.get_field_vector_at(s.r + 1.5 * s.R * v_hat)
        E_hat = E_hat/norm(E_hat)

        s.sigma = mlib.get_sigma_v3(L_E)
        s.v = mlib.get_velocity_v3(L_E) * E_hat

        if sigma.min() > 0.2 * sigma.max() and sigma.max() > 1e-2:
            # No dip in sigma ahead of the channel, so a streamer encounter
            s.keep = False

        dR = min(args.r_scale * mlib.get_radius_v3(L_E) - s.R,
                 norm(s.v) * args.dt)
        s.R = s.R + dR
        s.r = s.r + s.v * (args.dt - 0.99 * dR/norm(s.v))

    n_add = m_solver.adjust_refinement()
    mlib.update_sigma(m_solver.update_sigma, streamers, streamers_prev,
                      time, args.dt, 1e-9, step == 1)
    m_solver.solve(args.dt)

    time += args.dt
    m_solver.write_solution(f'{args.siloname}_{step:04d}', step, time)

    streamers = [s for s in streamers if s.keep]
    if len(streamers) == 0:
        raise ValueError('All streamers gone')


plt.plot(L_E_data)
plt.show()
