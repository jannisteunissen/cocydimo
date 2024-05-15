#!/usr/bin/env python3

from poisson import m_solver
import numpy as np
import h5py
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
parser.add_argument('h5file', type=str,
                    help='Input hdf5 file')
parser.add_argument('-siloname', type=str, default='output/simulation',
                    help='Base filename for output Silo files')
parser.add_argument('-plot', action='store_true',
                    help='Make plot of solution and original data')
args = parser.parse_args()


# Polynomial coefficients to describe the 'tip' of sigma_z
c = np.array([-1.05254158, -0.63331224,  0.15173764, -0.06770565, 1.0])
x_root = 0.87882323
x0 = 0.0015625
s0 = 6.013819664423325e-07

# Coefficients for linear fit of dsigma/dt
E_threshold = 5.2e6
L_fit_coeff = np.array([106058644.01730129, -8.410521689589444])

# Coefficients for linear fit of velocity
v_fit_coeff = np.array([2766342719656.4473, 298201.45554197486])

# Coefficients for quadratic-type fit of radius
R_fit_coeff = np.array([1.55, 4.18])


def get_high_field_length(phi, dz):
    E = np.abs(np.gradient(phi, dz))
    i_max = np.argmax(E)
    d_i = np.argmax(E[i_max:] < E_threshold)
    return d_i * dz


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


with h5py.File(args.h5file, 'r') as h5f:
    z = np.array(h5f['z_line'])  # Could in the future differ from z_grid
    sigma_z = np.array(h5f[f'run_{args.run}/sigmaz_{args.cycle}'])
    sigma = np.array(h5f[f'run_{args.run}/sigma_{args.cycle}'])
    rhs = np.array(h5f[f'run_{args.run}/rhs_{args.cycle}'])
    phi_bc = np.array(h5f[f'run_{args.run}/voltage'][args.cycle])
    times = np.array(h5f[f'run_{args.run}/time'])
    r_grid = np.array(h5f[f'run_{args.run}/r_grid'])
    z_grid = np.array(h5f[f'run_{args.run}/z_grid'])

n_cycles = len(times)
Nr, Nz = rhs.shape

dz = z_grid[1] - z_grid[0]
dr = r_grid[1] - r_grid[0]
dt = times[1] - times[0]
z0 = int(3.5e-3 / dz)           # To skip initial seed
dt_model = dt / args.dt_frac

sigma_z_pred = np.zeros((args.nsteps+1, Nz))
phi_z_pred = np.zeros((args.nsteps+1, Nz))
sigma_head = np.zeros((args.nsteps+1))
z_sigma_head = np.zeros((args.nsteps+1))
L_E = np.zeros((args.nsteps+1))

sigma_z_pred[0] = sigma_z
sigma_head[0] = sigma_z_pred[0, z0:].max()
i_sigma_head = np.argmax(sigma_z_pred[0, z0:]) + z0
z_sigma_head[0] = z[i_sigma_head]

domain_size = [Nr * dr, Nz * dz]
m_solver.initialize(domain_size, [Nr, Nz], args.box_size, phi_bc)
m_solver.set_rhs_and_sigma(rhs, sigma)

# Compute initial solution
m_solver.solve(0.0)
m_solver.write_solution(f'{args.siloname}_{0:04d}')
z_phi, phi = m_solver.get_line_potential(z[0], z[-1], len(z))
phi_z_pred[0] = phi

for step in range(1, args.nsteps+1):
    # Get input features for model
    L_E[step-1] = get_high_field_length(phi_z_pred[step-1], dz)

    dsigma_dt = get_dsigma_dt(L_E[step-1])
    v = get_velocity(sigma_head[step-1])
    R = get_radius(sigma_head[step-1])
    slope = dsigma_dt/v

    sigma_head[step] = sigma_head[step-1] + dsigma_dt * dt_model
    z_sigma_head[step] = z_sigma_head[step-1] + v*dt_model
    new_sigma_z = sigma_func(z, z_sigma_head[step], sigma_head[step], slope)
    sigma_z_pred[step] = np.where(z > z_sigma_head[step-1],
                                  new_sigma_z, sigma_z_pred[step-1])

    w_r = get_radial_weights(R, 20, dz)
    sigma_z_diff = sigma_z_pred[step] - sigma_z_pred[step-1]

    m_solver.update_sigma(z[0], z[-1], sigma_z_diff, R, w_r)
    m_solver.solve(dt_model)

    z_phi, phi = m_solver.get_line_potential(z[0], z[-1], len(z))
    phi_z_pred[step] = phi

    m_solver.write_solution(f'{args.siloname}_{step:04d}')


if args.plot:
    fig, ax = plt.subplots(4, sharex=True, layout='constrained')

    nsteps_data = min(args.nsteps//args.dt_frac + 1,
                      n_cycles - args.cycle)

    with h5py.File(args.h5file, 'r') as h5f:
        sigmaz_data = []
        phiz_data = []

        for i in range(args.cycle, args.cycle+nsteps_data):
            tmp = np.array(h5f[f'run_{args.run}/sigmaz_{i}'])
            ax[0].plot(z, tmp, c='gray')
            sigmaz_data.append(tmp)

            tmp = np.array(h5f[f'run_{args.run}/phiz_{i}'])
            ax[1].plot(z, phi, c='gray')
            ax[2].plot(z, -np.gradient(tmp, dz), c='gray')
            phiz_data.append(tmp)

    for step in range(0, args.nsteps+1, args.dt_frac):
        ax[0].plot(z, sigma_z_pred[step], ls='--')
        ax[1].plot(z, phi_z_pred[step], ls='--')
        ax[2].plot(z, -np.gradient(phi_z_pred[step], dz), ls='--')

    ax[0].set_ylabel('sigma')
    ax[1].set_ylabel('phi (V)')
    ax[2].set_ylabel('max(E)')

    ax[3].plot(z_sigma_head[:-1], L_E[:-1], label='model')

    # Compare with data
    L_E_data = [get_high_field_length(phi, dz) for phi in phiz_data]
    i_head_data = [np.argmax(sigma[z0:]) + z0 for sigma in sigmaz_data]
    z_head_data = [z[i] for i in i_head_data]

    ax[3].plot(z_head_data, L_E_data, label='data')
    ax[3].set_ylabel('length(E)')
    ax[3].legend()
    ax[-1].set_xlabel('z (m)')

    plt.show()
