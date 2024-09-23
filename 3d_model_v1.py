#!/usr/bin/env python3

from poisson_3d import m_solver_3d
import numpy as np
from numpy.linalg import norm
import h5py
import matplotlib.pyplot as plt
import copy
import argparse
from scipy.spatial.transform import Rotation

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


def get_high_field_length(phi, dz):
    E = np.abs(np.gradient(phi, dz))
    i_max = np.argmax(E)
    d_i = np.argmax(E[i_max:] < E_threshold)
    return d_i * dz


def get_radius(L_E):
    return args.r_scale * \
        np.sqrt(8.441501306469196e-08 + 0.6097078090287649 * L_E)


def get_velocity(sigma_max):
    return 307269.6823523594 + 1412173247560.5532 * sigma_max


def get_dsigma_dt(L_E):
    return -34.43873725692916 + 262462782.20047116 * L_E**2


Nx, Ny, Nz = args.grid_size
Lx = 20e-3
domain_size = [Lx, Lx, Lx]
phi_bc = -4e4

m_solver_3d.set_rod_electrode([0.5*Lx, 0.5*Lx, 0.0],
                              [0.5*Lx, 0.5*Lx, 0.15*Lx], 0.5e-3)
m_solver_3d.initialize(domain_size, [Nx, Ny, Nz], args.box_size, phi_bc)

# Compute initial solution
m_solver_3d.solve(0.0)
m_solver_3d.write_solution(f'{args.siloname}_{0:04d}')
z_phi, phi = m_solver_3d.get_line_potential([0.5*Lx, 0.5*Lx, 0.],
                                            [0.5*Lx, 0.5*Lx, Lx], Nz)


class Streamer(object):

    def __init__(self, r, v, R, sigma):
        self.r = np.array(r)
        self.v = np.array(v)
        self.R = R
        self.sigma = sigma
        self.keep = True
        self.is_branching = False
        self.branching_angle = None


def update_sigma(streamers_t1, streamers_t0, first_step):
    n = len(streamers_t1)

    if len(streamers_t0) != n:
        raise ValueError('Same number of streamers required')

    r = np.zeros((n, 3))
    r_prev = np.zeros((n, 3))
    sigma = np.zeros(n)
    sigma_prev = np.zeros(n)
    radius = np.zeros(n)
    radius_prev = np.zeros(n)

    for i in range(n):
        r[i] = streamers_t1[i].r
        r_prev[i] = streamers_t0[i].r
        sigma[i] = streamers_t1[i].sigma
        sigma_prev[i] = streamers_t0[i].sigma
        radius[i] = streamers_t1[i].R
        radius_prev[i] = streamers_t0[i].R

    m_solver_3d.update_sigma(r_prev, r, sigma_prev, sigma,
                             radius_prev, radius, first_step)


streamers = [Streamer([0.5*Lx, 0.5*Lx, 3.5e-3],
                      [0., 0., 1e6],
                      0.5e-3, 1e-6)]
             # Streamer([0.6*Lx, 0.5*Lx, 6.5e-3],
             #          [0., 0., 1e6],
             #          0.5e-3, 1e-6)]


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
        E_trace, sigma = m_solver_3d.get_head_trace(s.r, s.v, 2*s.R, 5)

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

    n_add = m_solver_3d.adjust_refinement()
    update_sigma(streamers, streamers_prev, step == 1)
    m_solver_3d.solve(args.dt)
    m_solver_3d.write_solution(f'{args.siloname}_{step:04d}')

    streamers = [s for s in streamers if s.keep]
    if len(streamers) == 0:
        raise ValueError('All streamers gone')
