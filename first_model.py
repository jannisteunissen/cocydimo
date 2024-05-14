#!/usr/bin/env python3

from poisson import m_solver
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Reduced discharge model')
parser.add_argument('-run', type=int, default=0,
                    help='Index of run to simulate')
parser.add_argument('-cycle', type=int, default=10,
                    help='Start from this cycle')
parser.add_argument('-r_scale', type=float, default=1.5,
                    help='Scale factor for the radius')
parser.add_argument('-nsteps', type=int, default=2,
                    help='How many steps to simulate')
parser.add_argument('-dt_frac', type=int, default=1,
                    help='Reduce dt in model by this factor')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo')
parser.add_argument('line_data', type=str,
                    help='Input numpy file with line data')
parser.add_argument('rz_data', type=str,
                    help='Input numpy file with 2D data')
parser.add_argument('-siloname', type=str, default='simulation',
                    help='Base filename for output Silo files')
args = parser.parse_args()


def get_high_field_length(phi, dz, E_threshold):
    E = np.abs(np.gradient(phi, dz))
    i_max = np.argmax(E)
    d_i = np.argmax(E[i_max:] < E_threshold)
    return d_i * dz


# Parameters to describe the 'tip' of sigma_z
c = np.array([-1.05254158, -0.63331224,  0.15173764, -0.06770565, 1.0])
x_root = 0.87882323
x0 = 0.0015625
s0 = 6.013819664423325e-07

# Coefficients for linear fit dsigma/dt
L_fit_coeff = np.array([106058644.01730129, -8.410521689589444])

# Coefficients for linear fit velocity
v_fit_coeff = np.array([2766342719656.4473, 298201.45554197486])

# Coefficients for quadratic fit radius
R_fit_coeff = np.array([1.55, 4.18])


def sigma_func(x, x_head, sigma_head, slope_behind):
    # Scale coordinate
    xs = (x - x_head) / (x0 * (sigma_head/s0)**0.5)
    cond = [xs < 0, xs < x_root, xs >= x0]
    values = [sigma_head * np.clip(1 + (x-x_head) *
                                   slope_behind/sigma_head, 0., 1.),
              np.polyval(c, xs) * sigma_head,
              np.zeros_like(xs)]
    return np.select(cond, values)


def get_radius(sigma_head):
    return args.r_scale * (R_fit_coeff[0] * sigma_head -
                           R_fit_coeff[1] * sigma_head**2)**0.5


def get_velocity(sigma_head):
    return v_fit_coeff[0] * sigma_head + v_fit_coeff[1]


def get_dsigma_dt(high_field_length):
    return L_fit_coeff[0] * high_field_length**2 + L_fit_coeff[1]


def get_radial_weights(r_max, n_points, dz):
    r = np.linspace(0., r_max, n_points)
    dr = r[1] - r[0]
    w = np.maximum(0.0, 1 - 3*(r/r_max)**2 + 2*(r/r_max)**3)
    w_sum = np.sum(2 * np.pi * r * dr * w)
    return w/w_sum


line_data = np.load(args.line_data)
phi_rz_data = np.load(args.rz_data)

line_mask = (line_data['run_index'] == args.run)
rz_mask = (phi_rz_data['run_index'] == args.run)

n_cycles = line_mask.sum()
sigma_z = line_data['sigma_z'][line_mask]
phi_z = line_data['phi_z'][line_mask]
cycles = line_data['cycles'][line_mask]
times = line_data['times'][line_mask]
z = line_data['z']
phi_rz = phi_rz_data['phi'][rz_mask]
rhs_rz = phi_rz_data['rhs'][rz_mask]
sigma_rz = phi_rz_data['sigma'][rz_mask]

Nr = phi_rz.shape[1]
Nz = phi_rz.shape[2]

# This is a pretty good approximation of the applied voltage
phi_bc = phi_rz[args.cycle, -1, -1] + phi_rz[args.cycle, -1, 0]
dz = z[1] - z[0]
dt = times[args.cycle+1] - times[args.cycle]
z0 = int(3.5e-3 / dz)
dt_model = dt / args.dt_frac

m_solver.initialize([Nr * dz, Nz * dz], [Nr, Nz], args.box_size, phi_bc)
m_solver.set_rhs_and_sigma(rhs_rz[args.cycle], sigma_rz[args.cycle])

sigma_z_pred = np.zeros((args.nsteps+1, Nz))
phi_z_pred = np.zeros((args.nsteps+1, Nz))
sigma_head = np.zeros((args.nsteps+1))
z_sigma_head = np.zeros((args.nsteps+1))
z_Emax = np.zeros((args.nsteps+1))
L_E = np.zeros((args.nsteps+1))

sigma_z_pred[0] = sigma_z[args.cycle]
phi_z_pred[0] = phi_z[args.cycle]
sigma_head[0] = sigma_z_pred[0, z0:].max()
i_sigma_head = np.argmax(sigma_z_pred[0, z0:]) + z0
z_sigma_head[0] = z[i_sigma_head]
z_Emax[0] = z[np.argmax(np.abs(np.gradient(phi_z_pred[0])))]

for step in range(1, args.nsteps+1):
    # Get input features for model
    L_E[step-1] = get_high_field_length(phi_z_pred[step-1], dz, 5.2e6)

    dsigma_dt = get_dsigma_dt(L_E[step-1])
    sigma_head[step] = sigma_head[step-1] + dsigma_dt * dt_model

    v = get_velocity(sigma_head[step-1])
    R = get_radius(sigma_head[step-1])
    slope = dsigma_dt/v

    z_sigma_head[step] = z_sigma_head[step-1] + v*dt_model
    new_sigma_z = sigma_func(z, z_sigma_head[step], sigma_head[step], slope)
    sigma_z_pred[step] = np.where(z > z_sigma_head[step-1],
                                  new_sigma_z, sigma_z_pred[step-1])
    sigma_z_diff = sigma_z_pred[step] - sigma_z_pred[step-1]

    w_r = get_radial_weights(R, 20, dz)
    m_solver.update_sigma(z[0], z[-1], sigma_z_diff, R, w_r)
    m_solver.solve(dt_model)
    m_solver.write_solution(f'{args.siloname}_{step:04d}')
    z_phi, phi = m_solver.get_line_potential(z[0], z[-1], len(z))

    phi_z_pred[step] = phi
    z_Emax[step] = z[np.argmax(np.abs(np.gradient(phi_z_pred[step])))]

fig, ax = plt.subplots(4, sharex=True, layout='constrained')

nsteps_data = min(args.nsteps//args.dt_frac + 1,
                  len(sigma_z) - args.cycle)

for i in range(nsteps_data):
    ax[0].plot(z, sigma_z[args.cycle+i], c='gray')
    ax[1].plot(z, phi_z[args.cycle+i], c='gray')
    ax[2].plot(z, -np.gradient(phi_z[args.cycle+i], dz), c='gray')

for step in range(0, args.nsteps+1, args.dt_frac):
    ax[0].plot(z, sigma_z_pred[step], ls='--')
    ax[1].plot(z, phi_z_pred[step], ls='--')
    ax[2].plot(z, -np.gradient(phi_z_pred[step], dz), ls='--')

ax[0].set_ylabel('sigma')
ax[1].set_ylabel('phi (V)')
ax[2].set_ylabel('max(E)')

ax[3].plot(z_Emax[:-1], L_E[:-1], label='model')

# Compare with data
L_E_data = [get_high_field_length(phi, dz, 5.2e6)
            for phi in phi_z[args.cycle:args.cycle+nsteps_data]]
i_z_head_data = [np.argmax(np.abs(np.gradient(phi, dz)))
                 for phi in phi_z[args.cycle:args.cycle+nsteps_data]]
z_head_data = [z[i] for i in i_z_head_data]
ax[3].plot(z_head_data, L_E_data, label='data')
ax[3].set_ylabel('length(E)')
ax[3].legend()
ax[-1].set_xlabel('z (m)')

plt.show()
