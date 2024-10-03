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
phi_bc = -4e4

m_solver.set_rod_electrode([0.5*Lx, 0.5*Lx, 0.0],
                           [0.5*Lx, 0.5*Lx, 0.15*Lx], 0.5e-3)
m_solver.initialize(domain_size, [Nx, Ny, Nz], args.box_size, phi_bc)

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = np.loadtxt('k_eff_air.txt').T
m_solver.store_k_eff(table_fld[0], table_fld[-1], table_k_eff)

z_phi, phi = m_solver.get_line_potential([0.5*Lx, 0.5*Lx, 0.],
                                         [0.5*Lx, 0.5*Lx, Lx], Nz)

# Locate Emax
Emax, r_Emax = m_solver.get_max_field_location()

streamers = [mlib.Streamer([0.5*Lx, 0.5*Lx, r_Emax[2]],
                           [0., 0., 1e6], 0.5e-3, 1e-6)]

time = 0.

for step in range(1, args.nsteps+1):

    if step == 8 or step == 20:
        streamers[0].is_branching = True
        streamers.append(copy.copy(streamers[0]))
        streamers[0].branching_angle = 30/180. * np.pi
        streamers[-1].branching_angle = -30/180. * np.pi

    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:
        # Get samples of the electric field and sigma in the direction of the
        # previous streamer velocity
        E_trace, sigma = m_solver.get_head_trace(s.r, s.v, 2*s.R, 5)

        # Take the field at the farthest sample. This field will tend to bend
        # towards the background electric field.
        E_hat = E_trace[:, -1]
        E_hat = E_hat / norm(E_hat)

        if s.is_branching:
            rot = Rotation.from_euler('y', s.branching_angle)
            E_hat = rot.apply(E_hat)
            s.is_branching = False
        elif sigma.min() > 0.2 * sigma.max():
            # No dip in sigma ahead of the channel, so a streamer encounter
            s.keep = False

        s.v = 1e6 * E_hat
        s.r = s.r + s.v * args.dt
        s.R = 0.5e-3
        s.sigma = s.sigma

        print(s)

    n_add = m_solver.adjust_refinement()
    mlib.update_sigma(m_solver.update_sigma, streamers, streamers_prev,
                      time, args.dt, 1e-9, step == 1)
    m_solver.solve(args.dt)

    time += args.dt
    m_solver.write_solution(f'{args.siloname}_{step:04d}', step, time)

    streamers = [s for s in streamers if s.keep]
    if len(streamers) == 0:
        raise ValueError('All streamers gone')
