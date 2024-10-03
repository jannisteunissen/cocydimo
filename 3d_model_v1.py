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
parser.add_argument('-r_scale', type=float, default=1.25,
                    help='Scale factor for the radius')
parser.add_argument('-nsteps', type=int, default=2,
                    help='How many steps to simulate')
parser.add_argument('-grid_size', type=int, nargs=3, default=[128, 128, 128],
                    help='Size of computational grid')
parser.add_argument('-dt', type=float, default=2e-10,
                    help='Time step (s)')
parser.add_argument('-phi_bc', type=float, default=-4e4,
                    help='Applied potential (V)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo')
parser.add_argument('-siloname', type=str, default='output/simulation_3d',
                    help='Base filename for output Silo files')
args = parser.parse_args()

E_threshold = 5e6

Nx, Ny, Nz = args.grid_size
Lx = 20e-3
domain_size = [Lx, Lx, Lx]

m_solver.set_rod_electrode([0.5*Lx, 0.5*Lx, 0.0],
                           [0.5*Lx, 0.5*Lx, 0.15*Lx], 0.5e-3)
m_solver.initialize(domain_size, [Nx, Ny, Nz], args.box_size, args.phi_bc)

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt('k_eff_air.txt').T
m_solver.store_k_eff(table_fld[0], table_fld[-1], 0*table_k_eff)

z_phi, phi = m_solver.get_line_potential([0.5*Lx, 0.5*Lx, 0.],
                                         [0.5*Lx, 0.5*Lx, Lx], Nz)

# Locate Emax
Emax, r_Emax = m_solver.get_max_field_location()

streamers = [mlib.Streamer([0.5*Lx, 0.5*Lx, 3.4e-3],
                           [0., 0., 1.0], 1.0e-4, 0.0)]

time = 0.

for step in range(1, args.nsteps+1):

    # if step == 8 or step == 20:
    #     streamers[0].is_branching = True
    #     streamers.append(copy.copy(streamers[0]))
    #     streamers[0].branching_angle = 30/180. * np.pi
    #     streamers[-1].branching_angle = -30/180. * np.pi

    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:
        # Get samples of sigma in the vicinity of the streamer
        _, sigma = m_solver.get_var_along_line('sigma', s.r, s.v, 2*s.R, 5)

        # Get samples of |E| ahead of the streamer to determine L_E
        r_tip = s.r + s.R * s.v/norm(s.v)
        z, E = m_solver.get_var_along_line('E_norm', r_tip, s.v, 2e-3, 100)
        L_E = mlib.get_high_field_length(z, E)

        # Take the field ahead of the streamer. This field will tend to bend
        # towards the background electric field.
        v_hat = s.v/norm(s.v)
        E_hat = m_solver.get_field_vector_at(s.r + 2 * s.R * v_hat)
        E_hat = E_hat/norm(E_hat)

        old_sigma = s.sigma
        old_v = s.v.copy()
        s.sigma = max(0.5 * s.sigma, min(s.sigma + mlib.max_dt_sigma * args.dt,
                                         mlib.get_sigma(L_E)))

        if s.is_branching:
            rot = Rotation.from_euler('y', s.branching_angle)
            E_hat = rot.apply(E_hat)
            s.is_branching = False
        elif sigma.min() > 0.2 * sigma.max() and sigma.max() > 1e-2:
            # No dip in sigma ahead of the channel, so a streamer encounter
            s.keep = False
        elif s.sigma < mlib.sigma_min:
            # TODO: adjust for initial streamer growth
            print('Streamer removed due to too low sigma')
            s.keep = False

        dR = mlib.get_radius(s.sigma, args.r_scale) - s.R
        s.R = s.R + dR
        s.v = mlib.get_velocity(s.sigma) * E_hat
        s.r = s.r + 0.5 * (old_v + s.v) * args.dt

        print(L_E, s)

    n_add = m_solver.adjust_refinement()
    mlib.update_sigma(m_solver.update_sigma, streamers, streamers_prev,
                      time, args.dt, 1e-9, step == 1)
    m_solver.solve(args.dt)

    time += args.dt
    m_solver.write_solution(f'{args.siloname}_{step:04d}', step, time)

    streamers = [s for s in streamers if s.keep]
    if len(streamers) == 0:
        raise ValueError('All streamers gone')
