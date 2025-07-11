#!/usr/bin/env python3

import numpy as np


class Streamer():
    """
    A class to represent a streamer discharge in a simulation.

    Attributes:
        ndim: Dimension of the space (length of r).
        r: Head position of the streamer (numpy array).
        v: Velocity of the streamer (numpy array).
        R: Radius of the streamer.
        sigma: Line conductivity of the streamer.
        keep: A flag to indicate whether to keep the streamer active.
        is_branching: A flag to indicate whether the streamer will branch.
        branching_angle: Angle at which the streamer branches, if applicable.
        branching_axis: Axis around which the streamer will branch, if applicable.
    """

    def __init__(self, r, v, R, sigma):
        self.ndim = len(r)
        self.r = np.array(r)
        self.v = np.array(v)
        self.R = R
        self.sigma = sigma
        self.L_E = 0.0
        self.n_steps = 0
        self.keep = True
        self.is_branching = False
        self.branching_angle = None
        self.branching_axis = None

    def __repr__(self):
        """Return a string representation of the Streamer instance."""
        with np.printoptions(formatter={'float': lambda x: format(x, '.2E')}):
            r = f'Streamer(r = {self.r}, v = {self.v}, ' + \
                f'sigma = {self.sigma:.2e}, R = {self.R:.2e}, ' + \
                f'steps = {self.n_steps})'
        return r


class AirStreamerModel():
    """A class to model the properties of positive streamer discharges in air
    based on the length of the high-field region (L_E).

    Attributes
    ----------
    E_threshold : float
        A threshold value used for fitting data, set to 5e6.
    c0 : float
        A correction factor of order unity.
    dz0 : float
        The grid spacing used when fitting data.

    Methods
    -------
    get_radius(L_E)
        Calculates the radius of the streamer based on L_E.

    get_velocity(L_E)
        Calculates the velocity of the streamer based on L_E.

    get_sigma(L_E)
        Calculates the line conductivity (sigma) of the streamer based on L_E.

    get_L_E(z, E, dz=None) Calculates the length of the high-field region
        (L_E) based on the electric field data.

    """

    E_threshold = 5e6  # Used for fitting data

    def __init__(self, c0=0.0, dz0=0.0):
        self.c0 = c0
        self.dz0 = dz0

    @staticmethod
    def get_radius(L_E):
        """Calculate the radius of the streamer.

        Parameters
        ----------
        L_E : float or ndarray
            The length of the high-field region ahead of the streamer.

        Returns
        -------
        float or ndarray
            The calculated radius of the streamer.

        """
        R = np.where(L_E < 1e-3,
                     2.897e-05 + 1.229 * L_E,
                     2.897e-05 + 1.229 * 1e-3 + 6.271e-01 * (L_E - 1e-3))
        return R

    @staticmethod
    def get_velocity(L_E):
        """Calculate the velocity of the streamer.

        Parameters
        ----------
        L_E : float or ndarray
            The length of the high-field region ahead of the streamer.

        Returns
        -------
        float or ndarray
            The calculated velocity of the streamer.

        """
        return 1.78e+09 * L_E

    @staticmethod
    def get_sigma(L_E):
        """Calculate the line conductivity (sigma) of the streamer.

        Parameters
        ----------
        L_E : float or ndarray
            The length of the high-field region ahead of the streamer.

        Returns
        -------
        float or ndarray
            The calculated line conductivity of the streamer.

        """
        sigma = np.where(L_E < 1e-3,
                         1e-8 + 1.397 * L_E**2,
                         1e-8 + 1.397 * 1e-6 + (L_E - 1e-3) * 2 * 1.397 * 1e-3)
        return sigma

    def get_L_E(self, z, E, dz=None):
        """Calculate the length of the high-field region (L_E) based on the
        electric field data.

        Parameters
        ----------
        z : ndarray
            The spatial coordinates.
        E : ndarray
            The electric field profile (in the forward direction)
        dz : float, optional
            The actual grid spacing.

        Returns
        -------
        float
            The length of the high-field region (L_E).

        """
        # Locate maximum of E
        i_max = np.argmax(E)

        # Determine distance between maximum and threshold.
        i_diff = np.argmax(E[i_max:] < self.E_threshold)

        if i_diff == 0:
            # Threshold was not reached (or at domain boundary)
            return 0.0

        i_threshold = i_max + i_diff
        # Convert to length
        L_E = abs(z[i_threshold] - z[i_max])

        # Apply correction for finite grid spacing
        if dz is not None:
            L_E += self.c0 * (dz - self.dz0)

        return L_E


def update_sigma(method, streamers_t1, streamers_t0, time, dt,
                 channel_delay, first_step):
    """
    Update the conductivity for streamers on a grid based on current and previous time steps.

    Parameters:
        method: Function to call for updating the conductivity.
        streamers_t1: List of streamers at the current time step.
        streamers_t0: List of streamers at the previous time step.
        time: Current simulation time.
        dt: Time increment since the last update.
        channel_delay: Delay time for updating conductivity within a channel,
                       accounting for processes such as attachment.
        first_step: Boolean indicating if this is the first time step of the simulation.

    Raises:
        ValueError: If the number of streamers at the current and previous time steps do not match.
    """
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
