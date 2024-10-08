#!/usr/bin/env python3

import numpy as np

E_threshold = 5e6               # Used in fits
max_dt_sigma = 1e3              # Maximum time derivative of sigma
sigma_min = 5e-9                # Minimum sigma


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


def get_high_field_length(z, E):
    i_max = np.argmax(E)
    i_diff = np.argmax(E[i_max:] < E_threshold)

    if i_diff == 0:
        raise ValueError('The threshold field was not reached')

    i_threshold = i_max + i_diff
    return abs(z[i_threshold] - z[i_max])


def get_radius(sigma, r_scale):
    tmp = np.maximum(0., 6.78662043384393e-08 + 0.6133314336312898 * sigma)
    return r_scale * np.sqrt(tmp)


def get_velocity(sigma):
    return 318891.6649006611 + 1417943826064.345 * sigma


def get_sigma(L_E):
    return -3.9234629547877727e-07 + 0.0015436863032232415 * L_E
