# Module including the ODE solver based on the code in: "https://github.com/MD-CWI/streamer-head-ode/tree/main".
# The ODE solver is described in Bouwman et al PSST, 2025 (https://doi.org/10.1088/1361-6595/adaf53).


#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from scipy.signal import savgol_filter
from scipy.optimize import root_scalar, root
import scipy.constants as c
import pickle
import sys
import argparse

def ODE_model(Q=0.5, v=5e4, R=None, E_max=None, fixed_I_ph=None, 
              E_bg=4.5e5, ne_bc=1e10, T=300.,
              p=1.0, frac_O2=0.2, p_q=40e-3, photoi_eff=0.075, L_factor=1.0,
              E_threshold=1.0, rtol=1e-5, neg_streamer=False,
              use_diffusion=False,
              alpha_table="input/reduced_alpha_phelps_jannis.txt",
              eta_table="input/reduced_eta_phelps_jannis.txt",
              mu_table="input/reduced_mu_phelps_jannis.txt",
              dif_table="input/reduced_dif_Phelps.txt"):

    # Gas number density from ideal gas law
    N0 = p * 1.0e5 / (c.k * T)
    
    # Conversion factor for Townsend to SI units
    Td_to_SI = 1e-21 * N0
    
    # Minimal z-coordinate to consider
    z_min = 1e-5
    
    # Indices of variables
    i_ne, i_ni, i_q = 0, 1, 2
    if use_diffusion:
        i_w = 3
    
    # Load transport data files
    TD_alpha = np.loadtxt(alpha_table, skiprows=2).T
    TD_eta = np.loadtxt(eta_table, skiprows=2).T
    TD_mu = np.loadtxt(mu_table, skiprows=2).T
    TD_D = np.loadtxt(dif_table, skiprows=2).T
    
    # Construct functions to interpolate transport data
    f_alpha = interp1d(Td_to_SI * TD_alpha[0], TD_alpha[1]*N0,
                       fill_value='extrapolate')
    f_eta = interp1d(Td_to_SI * TD_eta[0], TD_eta[1]*N0,
                     fill_value='extrapolate')
    f_mu = interp1d(Td_to_SI * TD_mu[0], TD_mu[1]/N0,
                    fill_value='extrapolate')
    f_dif = interp1d(Td_to_SI * TD_D[0], TD_D[1]/N0,
                     fill_value='extrapolate')
    
    # Apply filter when computing dmu/dE
    TD_dmu_dE = np.gradient(TD_mu[1]/N0, Td_to_SI * TD_mu[0])
    TD_dmu_dE = savgol_filter(TD_dmu_dE, 15, 2)
    f_dmu_dE = interp1d(Td_to_SI * TD_mu[0], TD_dmu_dE,
                        bounds_error=False,
                        fill_value=(TD_dmu_dE[0], TD_dmu_dE[-1]))
    
    # Calculate diffusion coefficients
    TD_dD_dE = np.gradient(TD_D[1]/N0, Td_to_SI * TD_D[0])
    TD_dD_dE = savgol_filter(TD_dD_dE, 15, 2)
    
    f_dD_dE = interp1d(Td_to_SI * TD_D[0], TD_dD_dE,
                       bounds_error=False,
                       fill_value=(TD_dD_dE[0], TD_dD_dE[-1]))
    
    # Determine critical field
    sol = root_scalar(lambda E: f_alpha(E) - f_eta(E), bracket=[0., 1e7])
    E_breakdown = sol.root
    
    # For numerical reasons, it is convenient if numbers are around unity.
    # Therefore we use this scaling factor. I_ph represents the number of photons
    # produced per unit time.
    I_ph_fac = 1e15
    
    # Scaling parameters for root-finding algorithm
    R0_scale = 1e-3
    E_max_scale = 1e7
    
    # We add
    v_scale = 1e5
    Q_scale = 1e-3
    
    # Initial condition
    if neg_streamer:
        if E_max is not None:
            E_max = -E_max
        E_bg = -E_bg
        Q = -Q
    else:
        E_bg = E_bg
        Q = Q
    
    if use_diffusion:
        y0 = np.zeros(4)
        y0[[i_ne, i_ni, i_q, i_w]] = [ne_bc, ne_bc, Q, 0.0]
    else:
        y0 = np.zeros(3)
        y0[[i_ne, i_ni, i_q]] = [ne_bc, ne_bc, Q]
    
    def get_E(z, q, E_bg):
        """Get electric field"""
        return E_bg + q/z**2 
    
    
    def get_dE_dz(z, q, dq_dz, E_bg):
        """Spatial derivative of electric field"""
        return dq_dz/z**2 - 2 * q/z**3
    
    
    def get_z_crit(Q, E_bg):
        """Solve for location of breakdown field, but avoid negative charge"""
        Q_nonneg = np.maximum(Q, 0.)
        if R is not None and E_max is None:
            return R
        else:
            return np.sqrt(Q_nonneg / (E_breakdown - E_bg))
    
    
    def f_absorption(z):
        """Absorption function of Zheleznyak model"""
        xmin = 2.6e3 * frac_O2 * p
        xmax = 150e3 * frac_O2 * p
        z = np.maximum(z, 1e-12)    # prevent negative values
        return (np.exp(-xmin * z) - np.exp(-xmax * z))/(z * np.log(xmax/xmin))
    
    
    def get_S_ph(z, I_ph, R):
        """Get photoionization source term, assuming that all photons come from a
        point source located at z=R/2
    
        I_ph is the number of photons produced per unit time:
        I_ph = p_q/(p + p_q) * photoi_eff * n_ch * pi * R**2 * v,
        where n_ch * pi * R**2 * v is the ionization produced per unit time
        """
        zstar = np.maximum(z - 0.5 * R, z_min)
        S_ph = abs(I_ph) * I_ph_fac * f_absorption(zstar) / (4 * np.pi * zstar**2)
        return S_ph
    
    
    def ODE_model_rhs(z, y, I_ph, R):
        """Right-hand side of the ODE model"""
        dy = np.zeros_like(y)
        if use_diffusion:
            ne, ni, q, w = y
        else:
            ne, ni, q = y
    
        if neg_streamer:
            rho_eps = (-ni + ne) * c.e / c.epsilon_0
        else:
            rho_eps = (ni - ne) * c.e / c.epsilon_0
        
        E = get_E(z, q, E_bg)
        E_abs = np.abs(E)
        
        alpha, eta = f_alpha(E_abs), f_eta(E_abs)
        mu, dmu_dE = f_mu(E_abs), f_dmu_dE(E_abs)
        D, dD_dE   = f_dif(E_abs), f_dD_dE(E_abs)
        
        S_ph = get_S_ph(z, I_ph, R) if I_ph > 0 else 0.
        src = (alpha - eta) * mu * E_abs * ne + S_ph
    
        dq_dz = rho_eps * z**2
        dE_dz = get_dE_dz(z, q, dq_dz, E_bg)
        dmu_dz = dmu_dE * dE_dz
        dD_dz = dD_dE * dE_dz
    
        # d/dz [ne, ni, q]
        if neg_streamer and not use_diffusion:
            # Negative streamers
            dy[i_ne] = -1 / (v - mu * E) * (src - dmu_dz * E * ne -
                                                 mu * rho_eps * ne) 
            dy[i_ni] = -src / v
            dy[i_q] = dq_dz
        elif not neg_streamer and not use_diffusion:
            # Postive streamers
            dy[i_ne] = -1 / (v + mu * E) * (src + dmu_dz * E * ne +
                                                 mu * rho_eps * ne)
            dy[i_ni] = -src / v 
            dy[i_q] = dq_dz
        elif use_diffusion:
            #"Include diffusion"
            if neg_streamer:
                dy[i_w] = 1/D * ((v + dD_dz - mu * E) * w - ne * E * dmu_dz - 
                                                        mu * ne * rho_eps + src)
            else:
                dy[i_w] = 1/D * ((v + dD_dz + mu * E) * w + ne * E * dmu_dz +
                                                        mu * ne * rho_eps + src)
            dy[i_ne] = w
            dy[i_ni] = -src / v
            dy[i_q] = dq_dz
        else:
            sys.exit("Je moet wel wat kiezen")
        return dy
    
    def get_radius(sol):
        """Determine the radius for an ODE solution"""
        field = get_E(sol.t, sol.y[i_q], E_bg)
        return sol.t[np.argmax(field)]
    
    
    def get_E_max(sol):
        """Determine E_max for an ODE solution"""
        field = get_E(sol.t, sol.y[i_q], E_bg)
        return field.max()
    
    
    def stop_condition(z, y, *other_args):
        """Stop integrating (towards z = 0) when this condition is met"""
        return get_E(z, y[i_q], E_bg) - E_threshold * E_breakdown
    
    
    stop_condition.terminal = True
    stop_condition.direction = -1.
    
    
    def solve_fixed_I_ph(y0, I_ph, R=R):
        """Integrate ODE using I_ph and R for the photoionization source term"""
        y0[i_q] = abs(y0[i_q])
        z_crit = get_z_crit(y0[i_q], E_bg)
        z_max = max(L_factor * z_crit, 2 * z_min)
    
        # We need high accuracy here to later do root finding on solutions
        sol = solve_ivp(ODE_model_rhs, [z_max, z_min],
                        y0, dense_output=True,
                        events=stop_condition, rtol=rtol,
                        first_step=1e-6,
                        args=(I_ph, R), method='Radau')
        return sol
    
    
    def residual_photoi(x0, y0):
        """Return difference between guess and updated guess for I_ph and R"""
        I_ph_guess, R_guess = x0
        sol = solve_fixed_I_ph(y0, I_ph_guess, R_guess)
        n_ch = sol.y[i_ni].max()
        R = get_radius(sol)
        I_ph = p_q / (1 + p_q) * photoi_eff * \
            n_ch * np.pi * R**2 * v / I_ph_fac
        return (I_ph - I_ph_guess)/max(I_ph, 1.0), (R - R_guess)/R0_scale
    
    
    def residual_R(Q_guess, y0, I_ph):
        """Return relative difference between obtained radius and goal radius"""
        y0[i_q] = Q_guess
        sol = solve_fixed_I_ph(y0, I_ph)
        R_obtained = get_radius(sol)
        return (R_obtained - R)/R0_scale
    
    
    def residual_E_max(Q_guess, y0, I_ph, R):
        """Return relative difference between obtained E_max and goal E_max"""
        y0[i_q] = Q_guess
        sol = solve_fixed_I_ph(y0, I_ph, R)
        E_max = get_E_max(sol)
        return (E_max - E_max)/E_max_scale
    
    
    def residual_photoi_R(x0, y0):
        """Return relative differences for I_ph and radius"""
        I_ph_guess = x0[0]
        y0[i_q] = x0[1]
        sol = solve_fixed_I_ph(y0, I_ph_guess)
        n_ch = sol.y[i_ni].max()
        R_obtained = get_radius(sol)
        I_ph = p_q / (1 + p_q) * photoi_eff * \
            n_ch * np.pi * R**2 * v / I_ph_fac
        return (I_ph - I_ph_guess)/max(I_ph, 1.0), (R_obtained - R)/R0_scale
    
    
    def residual_photoi_E_max(x0, y0):
        """Return relative differences for I_ph and E_max"""
        I_ph_guess, R_guess, y0[i_q] = x0
        # print(x0)
        sol = solve_fixed_I_ph(y0, I_ph_guess, R_guess)
        # print("done")
        n_ch = sol.y[i_ni].max()
        E_max = get_E_max(sol)
        R = get_radius(sol)
        I_ph = p_q / (1 + p_q) * photoi_eff * \
            n_ch * np.pi * R**2 * v / I_ph_fac
        residual = [(I_ph - I_ph_guess)/max(I_ph, 1.0), (R - R_guess)/R0_scale,
                    (E_max - E_max)/E_max_scale]
        return residual
    
    def residual_v_Q(x, y0):
        """
        Solve for velocity and Q
        while matching R and E_max.
        """
        
        # Update parameters
        v = x[0] * v_scale
        y0[i_q] = x[1] * Q_scale
        I_ph = fixed_I_ph
    
        # Solve ODE
        sol = solve_fixed_I_ph(y0, I_ph)
        R_obtained = get_radius(sol)
        E_max = get_E_max(sol)
        residual = [(R_obtained - R) / R0_scale, (E_max - E_max) / E_max_scale]
        return residual
    
    def check_convergence(sol):
            conv_fac = 1e-2
            if np.abs(sol.fun).max() > conv_fac:
                sys.stderr.write(str(sol) + '\n')
                raise ValueError('No convergence, adjust numerical parameters') 
    
    if R is not None and E_max is not None:
        print("Solving for v")
        x0 = (v / v_scale, Q / Q_scale)
        tmp = root( residual_v_Q, x0=x0, args=(y0,), method='lm', tol=rtol,
            options={'eps': 1e-2, 'factor': 0.1, 'maxiter': 500, 'xtol': 1e-4, 'ftol': 1e-4})
    
        # Recover physical variables
        v = tmp.x[0] * v_scale
        y0[i_q] = tmp.x[1] * Q_scale
    
        I_ph = fixed_I_ph
        sol = solve_fixed_I_ph(y0, I_ph)
    
    elif R is not None:
        # Solve iterative problem to find solution with given radius
        if fixed_I_ph is not None:
            # Fixed photoionization source, determine Q
            a, b = 30e1*Q, 9e-1*Q 
            #print(residual_R(a, y0, fixed_I_ph), residual_R(b, y0, fixed_I_ph))
            tmp = root_scalar(residual_R, args=(y0, fixed_I_ph),
                              x0=y0[i_q], bracket=[a,b],
                              rtol=rtol)
            I_ph = fixed_I_ph
            y0[i_q] = tmp.root
            sol = solve_fixed_I_ph(y0, I_ph)
        else:
            Q_guesses = np.array([1.0, 0.2, 5, 0.04, 25.]) * Q
            for Q in Q_guesses:
                y0[i_q] = Q
                tmp = root(residual_photoi_R, x0=(1., y0[i_q]), args=(y0, ),
                           method='lm', tol=rtol, options={'eps': 5e-3})
                if np.abs(tmp.fun).max() < 1e-2:
                    break
            check_convergence(tmp)
    
            I_ph = tmp.x[0]
            y0[i_q] = tmp.x[1]
            sol = solve_fixed_I_ph(y0, I_ph)
    elif E_max is not None:
        if fixed_I_ph is None:
            Q_guesses = np.array([1.0, 0.2, 5, 0.04, 25.]) * Q
            for Q in Q_guesses:
                y0[i_q] = Q 
                tmp = root(residual_photoi_E_max, x0=(1., 5e-6, y0[i_q]),
                           args=(y0, ), method='lm', tol=rtol,
                           options={'eps': 5e-3})
                if np.abs(tmp.fun).max() < 1e-2:
                    break
            check_convergence(tmp)
    
            I_ph, R, y0[i_q] = tmp.x
            sol = solve_fixed_I_ph(y0, I_ph, R)
        elif fixed_I_ph == 0.0:
            # Fixed photoionization source, determine Q
            a, b = -9.6e2*Q, 9.6e2*Q
            #print(residual_E_max(a, y0, fixed_I_ph, 0.), residual_E_max(b, y0, fixed_I_ph, 0.))
            tmp = root_scalar(residual_E_max, args=(y0, fixed_I_ph, 0.0),
                              x0=y0[i_q], bracket=[a,b],
                              rtol=rtol)
            I_ph = fixed_I_ph
            y0[i_q] = tmp.root
            sol = solve_fixed_I_ph(y0, I_ph)
        else:
            raise ValueError('fixed_I_ph > 0 not supported')
    else:
        if fixed_I_ph is not None:
            # Solve for fixed Q and I_ph
            I_ph = fixed_I_ph
            sol = solve_fixed_I_ph(y0, I_ph)
        else:
            I_ph, R = 1.0, 5e-4
            tmp = root(residual_photoi, x0=(I_ph, R), args=(y0, ),
                       method='lm', tol=rtol, options={'eps': 5e-3})
            I_ph, R = tmp.x
            sol = solve_fixed_I_ph(y0, I_ph, R)
    
    sol['E'] = get_E(sol.t, sol.y[i_q], E_bg)
    sol['q'] = sol.y[i_q] * 4 * np.pi * c.epsilon_0
    sol['z'] = sol.t
    sol['ne'] = sol.y[i_ne]
    sol['ni'] = sol.y[i_ni]
    field = get_E(sol.t, sol.y[i_q], E_bg)
    E_max = field.max()
    R = get_radius(sol)
    
    sol['R'] = R
    sol['Emax'] = E_max
    return sol
