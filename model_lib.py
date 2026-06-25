#!/usr/bin/env python3

import numpy as np

K_BOLTZMANN = 1.380649e-23  # J/K


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
        A threshold value used for fitting data, set to 5e6 V/m.
    N0 : float
        Gas number density at 300 K and 1 bar
    c0 : float
        A correction factor for L_E on coarser grids than dz0
    c1 : float
        A correction factor for L_E on finer grids than dz0
    dz0 : float
        The grid spacing used when fitting data.

    Methods
    -------
    get_radius(L_E, N)
        Calculates the radius of the streamer based on L_E.

    get_velocity(L_E, N)
        Calculates the velocity of the streamer based on L_E.

    get_sigma(L_E, N)
        Calculates the line conductivity (sigma) of the streamer based on L_E.

    get_L_E(z, E, N, dz=None) Calculates the length of the high-field region
        (L_E) based on the electric field data.

    """

    E_threshold = 5e6  # Used for fitting data
    N0 = 1e5 / (300 * 1.380649e-23)  # N0 at 300 K and 1 bar

    def __init__(self, c0=0.0, c1=0.0, dz0=0.0):
        self.c0 = c0
        self.c1 = c1
        self.dz0 = dz0

    def get_radius(self, L_E, N):
        """Calculate the radius of the streamer.

        Parameters
        ----------
        L_E : float or ndarray
            The length of the high-field region ahead of the streamer.
        N : float or ndarray
            The gas number density

        Returns
        -------
        float or ndarray
            The calculated radius of the streamer.

        """

        L_EN = L_E * N/self.N0
        R = self.N0/N * np.where(L_EN < 1e-3,
                                 2.897e-05 + 1.229 * L_EN,
                                 2.897e-05 + 1.229 * 1e-3 +
                                 6.271e-01 * (L_EN - 1e-3))
        return R

    def get_velocity(self, L_E, N):
        """Calculate the velocity of the streamer.

        Parameters
        ----------
        L_E : float or ndarray
            The length of the high-field region ahead of the streamer.
        N : float or ndarray
            The gas number density

        Returns
        -------
        float or ndarray
            The calculated velocity of the streamer.

        """
        L_EN = L_E * N/self.N0
        return 1.78e+09 * L_EN

    def get_sigma(self, L_E, N):
        """Calculate the line conductivity (sigma) of the streamer.

        Parameters
        ----------
        L_E : float or ndarray
            The length of the high-field region ahead of the streamer.
        N : float or ndarray
            The gas number density

        Returns
        -------
        float or ndarray
            The calculated line conductivity of the streamer.

        """
        L_EN = L_E * N/self.N0
        sigma = (self.N0/N)**2 * \
            np.where(L_EN < 1e-3,
                     1e-8 + 1.397 * L_EN**2,
                     1e-8 + 1.397 * 1e-6 + (L_EN - 1e-3) * 2 * 1.397 * 1e-3)
        return sigma

    def get_L_E(self, z, E, N, dz=None, prev=None):
        """Calculate the length of the high-field region (L_E) based on the
        electric field data.

        Parameters
        ----------
        z : ndarray
            The spatial coordinates.
        E : ndarray
            The electric field profile (in the forward direction)
        N : ndarray
            the gas number density
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
        threshold = self.E_threshold * (N/self.N0)
        i_diff = np.argmax(E[i_max:] < threshold)

        if i_diff == 0:
            if E[i_max:].min() >= threshold:
                # The field does not drop below the threshold
                return prev
            else:
                # Threshold was not reached
                return 0.0

        i_threshold = i_max + i_diff
        # Convert to length
        L_E = abs(z[i_threshold] - z[i_max])

        # Apply correction for finite grid spacing
        if dz is not None:
            if dz > self.dz0:
                L_E += self.c0 * (dz - self.dz0)
            else:
                L_E += self.c1 * (dz - self.dz0)

        return L_E


def update_sigma(ndim, method, streamers_t1, streamers_t0, time, dt,
                 channel_delay, first_step):
    """
    Update the conductivity for streamers on a grid based on current and previous time steps.

    Parameters:
        ndim: number of spatial dimensions
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

    if len(streamers_t0) != n:
        raise ValueError('Same number of streamers required')

    # We cannot easily pass zero-sized arrays, so allocate with at least one
    # element. Note that n is passed explicitly to method (last argument).
    m = max(n, 1)
    r = np.zeros((m, ndim))
    r_prev = np.zeros((m, ndim))
    sigma = np.zeros(m)
    sigma_prev = np.zeros(m)
    radius = np.zeros(m)
    radius_prev = np.zeros(m)

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
           dt, channel_delay, first_step, n)


def read_swarm_file(filename):
    """Parse a swarm-parameter file.

    Returns a dict with:
        'gas_composition' : {species: fraction, ...}
        'mobility', 'alpha', 'eta', 'three_body' : (E_N_Td, values) arrays
    """
    header_map = {
        'gas_composition':        'gas_composition',
        'Mobility':               'mobility',
        'Townsend ioniz. coef.':  'alpha',
        'Townsend attach. coef.': 'eta',
        'Three-body attachment':  'three_body',
    }

    data = {}
    current = None
    rows = []

    def store():
        if current is None or not rows:
            return
        if current == 'gas_composition':
            data[current] = {p[0]: float(p[1]) for p in (r.split() for r in rows)}
        else:
            arr = np.array([[float(x) for x in r.split()] for r in rows])
            data[current] = (arr[:, 0], arr[:, 1])

    with open(filename) as f:
        for line in (l.strip() for l in f):
            if not line or line.startswith('-'):
                continue
            new = next((name for key, name in header_map.items()
                        if line.startswith(key)), None)
            if new is not None:        # header line -> start a new section
                store()
                current, rows = new, []
            else:                      # data line
                rows.append(line)
    store()
    return data


def effective_ionization_rate(filename, gas_temperature, pressure_bar):
    """Compute the effective ionization rate coefficient [1/s].

    k_eff = v_drift * (alpha - eta - k3 * N_O2)
    """
    data = read_swarm_file(filename)

    # Number density (ideal gas law) and O2 density
    N = pressure_bar * 1e5 / (K_BOLTZMANN * gas_temperature)  # [1/m^3]
    N_O2 = data['gas_composition'].get('O2', 0.0) * N      # [1/m^3]

    # Reduced field in Townsend (1 Td = 1e-21 V m^2)
    E_over_N_Td = data['alpha'][0]

    interp = lambda key: np.interp(E_over_N_Td, *data[key])
    muN   = interp('mobility')
    alpha = interp('alpha') * N           # ionization     [1/m]
    eta   = interp('eta') * N             # attachment     [1/m]
    k3    = interp('three_body')          # 3-body rate    [m^6/s]

    v_drift = muN * E_over_N_Td * 1e-21  # [m/s]
    k_eff = v_drift * (alpha - eta) - k3 * N_O2**2

    return E_over_N_Td, k_eff
