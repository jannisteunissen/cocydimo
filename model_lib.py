#!/usr/bin/env python3

import numpy as np

E_threshold = 5e6               # Used for fitting data
max_dt_sigma = 100e3            # Maximum time derivative of sigma


class Streamer(object):
    def __init__(self, r, v, R, sigma):
        self.ndim = len(r)
        self.r = np.array(r)
        self.v = np.array(v)
        self.R = R
        self.sigma = sigma
        self.keep = True
        self.is_branching = False
        self.branching_angle = None
        self.branching_axis = None
        self.time_previous_branch = -1e100

    def __repr__(self):
        with np.printoptions(formatter={'float': lambda x: format(x, '.2E')}):
            r = f'Streamer(r = {self.r}, v = {self.v}, ' + \
                f'sigma = {self.sigma:.2e}, R = {self.R:.2e})'
        return r


def update_sigma(method, streamers_t1, streamers_t0, time, dt,
                 channel_delay, first_step):
    n = len(streamers_t1)
    ndim = streamers_t1[0].ndim

    if len(streamers_t0) != n:
        raise ValueError('Same number of streamers required')

    r = np.zeros((n, ndim))
    r_prev = np.zeros((n, ndim))
    sigma = np.zeros(n)
    sigma_prev = np.zeros(n)
    radius = np.zeros(n)
    radius_prev = np.zeros(n)

    for i in range(n):
        r[i] = streamers_t1[i].r
        r_prev[i] = streamers_t0[i].r
        sigma[i] = streamers_t1[i].sigma

        if first_step:
            # Use same value for both points
            sigma_prev[i] = streamers_t1[i].sigma
        else:
            sigma_prev[i] = streamers_t0[i].sigma

        radius[i] = streamers_t1[i].R
        radius_prev[i] = streamers_t0[i].R

    method(r_prev, r, sigma_prev, sigma, radius_prev, radius, time,
           dt, channel_delay, first_step)


def get_high_field_length(z, E, c0=None, dz0=None, dz=None):
    i_max = np.argmax(E)
    i_diff = np.argmax(E[i_max:] < E_threshold)

    if i_diff == 0:
        # Threshold was not reached (or at domain boundary)
        return None

    i_threshold = i_max + i_diff
    L_E = abs(z[i_threshold] - z[i_max])

    if c0 is not None:
        L_E += c0 * (dz - dz0)

    return L_E


def get_radius_v3(L_E):
    R = np.where(L_E < 1e-3,
                 2.897e-05 + 1.229 * L_E,
                 2.897e-05 + 1.229 * 1e-3 + 6.271e-01 * (L_E - 1e-3))
    return R


def get_velocity_v3(L_E):
    return 1.78e+09 * L_E


def get_sigma_v3(L_E):
    sigma = np.where(L_E < 1e-3, 1e-8 + 1.397 * L_E**2,
                     1e-8 + 1.397 * 1e-6 + (L_E - 1e-3) * 2 * 1.397 * 1e-3)
    return sigma
