#!/usr/bin/env python3
"""Cocydimo: Conducting Cylinder Discharge Model"""

import copy
import argparse
import numpy as np
import json
import os
from time import perf_counter
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
import model_lib as mlib
from poisson_3d import m_solver as p3d

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Cocydimo: Conducting Cylinder Discharge Model')
parser.add_argument('-r_scale', type=float, default=1.2,
                    help='Scale factor compared to electrodynamic radius')
parser.add_argument('-domain_size', type=float, nargs=3,
                    default=[80e-3, 80e-3, 80e-3], help='Domain size (m)')
parser.add_argument('-coarse_grid_size', type=int, nargs=3,
                    default=[16, 16, 16],
                    help='Size of coarse grid (Nx, Ny, Nz)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo (#cells)')
parser.add_argument('-rod_r0', type=float, nargs=3,
                    default=[40e-3, 40e-3, 0e-3],
                    help='First point of rod electrode (m)')
parser.add_argument('-rod_r1', type=float, nargs=3,
                    default=[40e-3, 40e-3, 40e-3],
                    help='Second point of rod electrode (m)')
parser.add_argument('-rod_radius', type=float, default=0.75e-3,
                    help='Radius of rod electrode (m)')
parser.add_argument('-n_steps', type=int, default=100,
                    help='How many steps to simulate')
parser.add_argument('-n_initial_streamers', type=int, default=5,
                    help='How many streamers to start with')
parser.add_argument('-r_start', type=float, nargs=3,
                    help='Initial location of initial streamers (m)')
parser.add_argument('-dt', type=float, default=5e-10,
                    help='Time step (s)')
parser.add_argument('-dt_factor', type=float, default=1.0,
                    help='Increase dt by this factor when streamers are gone')
parser.add_argument('-dz_data', type=float, default=30e-3/256,
                    help='Grid spacing used to obtain L_E from dataset (m)')
parser.add_argument('-phi_bc', type=float, default=-4e4,
                    help='Applied potential (V)')
parser.add_argument('-phi_factor_vs_time', type=str,
                    help='File with factor for applied voltage vs time')
parser.add_argument('-use_circuit', action='store_true',
                    help='Use R-C circuit for voltage source')
parser.add_argument('-capacitance', type=float, default=1e-9,
                    help='Capacitance C of the voltage source (farad)')
parser.add_argument('-resistance', type=float, default=300.0,
                    help='Internal resistance R of the voltage source (Ohm)')
parser.add_argument('-alpha', type=float, default=0.5,
                    help='Exponential smoothing coefficient')
parser.add_argument('-channel_update_delay', type=float, default=1e-9,
                    help='Delay for first updating channel conductivity (s).'
                    'Should be larger than the dielectric relaxation time.')
parser.add_argument('-channel_no_ionization', action='store_true',
                    help='Do not increase channel conductivity, which can be '
                    'problematic near domain boundaries)')
parser.add_argument('-channel_max_sigma', type=float, default=5.0,
                    help='Limit growth of volume conductivity to this value '
                    'to prevent issues at domain boundaries [A/(m V)]')
parser.add_argument('-channel_min_sigma', type=float, default=1e-9,
                    help='Keep electron condtivity above this value to allow '
                    'it to grow again [A/(m V)]')
parser.add_argument('-mu_electron', type=float, default=0.04,
                    help='Effective electron mobility at 1 bar, 300 K. '
                    'mu_e/mu_i is used to update ion conductivity (m2/(V s))')
parser.add_argument('-mu_ion', type=float, default=2e-4,
                    help='Effective ion mobility at 1 bar, 300 K')
parser.add_argument('-k_ion_recombination', type=float, default=1e-13,
                    help='Ion-ion recombination rate constant (m^3/s)')
parser.add_argument('-L_E_max', type=float, default=5e-3,
                    help='Maximum value of L_E (m)')
parser.add_argument('-L_E_min', type=float, default=1e-4,
                    help='Minimum value of L_E (m)')
parser.add_argument('-c0_L_E_dx', type=float, default=0.75,
                    help='Correction factor for L_E w.r.t. data grid spacing')
parser.add_argument('-c1_L_E_dx', type=float, default=0.0,
                    help='Correction factor for L_E when dx < dx_data')
parser.add_argument('-transport_data_file', type=str,
                    default='data/TD_N2_0.8_Phelps_O2_0.2_Phelps.txt',
                    help='Transport data file')
parser.add_argument('-k_eff_num_points', type=int, default=200,
                    help='Number of points to use internally for k_eff_table')
parser.add_argument('-poisson_rtol', type=float, default=1e-3,
                    help='Relative tolerance for Poisson solver')
parser.add_argument('-poisson_atol', type=float, default=1.0,
                    help='Absolute tolerance for Poisson solver')
parser.add_argument('-siloname', type=str, default='output/simulation_3d',
                    help='Base filename for output Silo files')
parser.add_argument('-write_eps', action='store_true',
                    help='Write epsilon (of Poisson eq.) to output')
parser.add_argument('-write_time', action='store_true',
                    help='Write time that channel was added to output')
parser.add_argument('-write_rhs', action='store_true',
                    help='Write r.h.s. of Poisson eq. to output')
parser.add_argument('-rng_seed', type=int,
                    help='Seed for the random number generator')
parser.add_argument('-c_b', type=float, default=15.,
                    help='Branching coeff. - larger means less branching')
parser.add_argument('-L_b', type=float, default=2.0e-4,
                    help='Branching coeff. (m) - thin streamers branch less')
parser.add_argument('-branch_gamma', type=float, default=90.0,
                    help='Angle (degrees) that affects branching')
parser.add_argument('-refine_E', type=float, default=3e6,
                    help='Refine if E is above this value (V/m)')
parser.add_argument('-derefine_E', type=float, default=2e6,
                    help='Derefine if E is below this value (V/m)')
parser.add_argument('-derefine_nlevels', type=int, default=1,
                    help='Derefine at most this many levels')
parser.add_argument('-min_dx', type=float, default=1e-4,
                    help='Minimum allowed grid spacing (m)')
parser.add_argument('-max_dx', type=float, default=2e-3,
                    help='Maximum allowed grid spacing (m)')
parser.add_argument('-max_dx_electrode', type=float, default=8e-4,
                    help='Maximum allowed grid spacing around electrode (m)')
parser.add_argument('-refine_max_dist_head', type=float, default=5e-3,
                    help='Max. distance from active head for refinement (m)')
parser.add_argument('-memory_limit', type=float, default=8.0,
                    help='Memory limit (GB)')
parser.add_argument('-print_performance', action='store_true',
                    help='Show performance information')
parser.add_argument('-steps_per_output', type=int, default=1,
                    help='Write output every N steps')
parser.add_argument('-gas_dynamics', action='store_true',
                    help='Simulate gas dynamics')
parser.add_argument('-pressure', type=float, default=1.0,
                    help='Gas pressure (bar)')
parser.add_argument('-temperature', type=float, default=300.0,
                    help='Gas temperature (Kelvin)')
parser.add_argument('-mean_molecular_weight', type=float, default=28.97,
                    help='Mean molecular weight of gas molecules (Dalton)')
parser.add_argument('-gas_gamma', type=float, default=1.4,
                    help='Gas adiabatic index')
parser.add_argument('-gas_fast_heat_factor', type=float, default=1.0,
                    help='Fraction of Joule heating that is immediately'
                    'converted to gas heating')
parser.add_argument('-gas_slow_heat_factor', type=float, default=0.0,
                    help='Fraction of Joule heating that is slowly converted'
                    'to gas heating')
parser.add_argument('-gas_slow_heat_timescale', type=float, default=20.0e-6,
                    help='Time scale for slow heating (s)')
parser.add_argument('-verbose', type=int, default=0,
                    help='How verbose the code is (> 0 shows more info)')

args = parser.parse_args()

# Make sure output folder exists
os.makedirs(os.path.dirname(args.siloname), exist_ok=True)

# Save settings
fname = f'{args.siloname}.json'
with open(fname, 'w') as f:
    json.dump(args.__dict__, f, indent=2)
    print(f'Wrote settings to {fname}')

model = mlib.AirStreamerModel(c0=args.c0_L_E_dx, c1=args.c1_L_E_dx,
                              dz0=args.dz_data)

np.random.seed(args.rng_seed)


def get_tau_branch(radius, velocity):
    """Get expected branching time for a streamer

    :param radius: radius (m)
    :param velocity: velocity (m/s)
    :returns: expected branching time
    """
    return args.c_b * radius/velocity * (1 + (args.L_b/radius)**2)


def find_orthogonal_unit_vector(y):
    """Randomly sample a unit vector orthogonal to y

    :param y: input vector
    :returns: unit vector orthogonal to y
    """
    orthvec = np.cross(y, np.random.uniform(-1., 1., 3))
    return orthvec / norm(orthvec)


def set_voltage(time, V0, voltage_table):
    voltage = V0
    if voltage_table is not None:
        voltage *= np.interp(time, voltage_table[0], voltage_table[1])
    p3d.set_voltage(voltage)
    print(time, voltage)


if args.phi_factor_vs_time is not None:
    voltage_table = np.loadtxt(args.phi_factor_vs_time).T
else:
    voltage_table = None


p3d.store_parameters(args.channel_min_sigma,
                     args.channel_max_sigma,
                     args.mu_electron,
                     args.mu_ion,
                     args.k_ion_recombination,
                     args.resistance,
                     args.capacitance,
                     args.verbose)

p3d.set_rod_electrode(args.rod_r0, args.rod_r1, args.rod_radius)

p3d.initialize_domain(args.domain_size, args.coarse_grid_size,
                      args.box_size, args.phi_bc, args.memory_limit,
                      args.write_eps, args.write_time, args.write_rhs,
                      args.gas_dynamics)

p3d.set_refinement(args.refine_E, args.derefine_E,
                   args.min_dx, args.max_dx,
                   args.max_dx_electrode, args.derefine_nlevels,
                   args.refine_max_dist_head,
                   args.poisson_rtol, args.poisson_atol)

# Initial gas density
N0 = 1e5 * args.pressure / (args.temperature * 1.380649e-23)

if args.gas_dynamics:
    p3d.set_gas(args.pressure, args.temperature, args.mean_molecular_weight,
                args.gas_gamma, args.gas_fast_heat_factor,
                args.gas_slow_heat_factor, args.gas_slow_heat_timescale)

dz = p3d.get_finest_grid_spacing()
print(f'Minimum grid spacing: {dz:.2e}')

# Compute initial solution
set_voltage(0.0, args.phi_bc, voltage_table)
p3d.solve(0.0, args.poisson_rtol, args.poisson_atol)
p3d.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = mlib.effective_ionization_rate(
    args.transport_data_file, args.temperature, args.pressure)

# Ensure table has uniform spacing
x = np.linspace(table_fld[0], table_fld[-1], args.k_eff_num_points)
y = np.interp(x, table_fld, table_k_eff)

if args.channel_no_ionization:
    y = np.minimum(y, 0.0)

p3d.store_k_eff(x[0], x[-1], y)

if args.r_start is not None:
    r_start = np.array(args.r_start)
else:
    Emax, r_start = p3d.get_max_field_location()
    print(f'Streamers start from {r_start}')

# Get L_E to estimate initial streamer radius
z, E, success = p3d.get_var_along_line('E_norm', r_start, [0., 0., 1.0],
                                       args.L_E_max, 2*args.L_E_max/dz)
if not success:
    raise RuntimeError('Interpolation error at r_start')
L_E = model.get_L_E(z, E, N0, dz, prev=args.L_E_max)

# Start with a smaller radius to approximate initial phase
radius0 = 0.5 * args.r_scale * model.get_radius(L_E, N0)

# Start with multiple streamers in random directions
streamers = []
for i in range(args.n_initial_streamers):
    # Sample multiple velocity directions
    axis = find_orthogonal_unit_vector([0., 0., 1.0])
    angle = np.random.uniform(0., args.branch_gamma*np.pi/180.)
    rot = Rotation.from_rotvec(axis * angle)
    v_hat = rot.apply([0., 0., 1.0])
    streamers.append(mlib.Streamer(r_start - [0., 0., radius0],
                                   v_hat, radius0, 0.0))

wct_refinement = 0.0
wct_update_sigma = 0.0
wct_poisson = 0.0
wct_output = 0.0
time = 0.0
V_cap = args.phi_bc
V_gap = args.phi_bc
step = 0

# Write current to a file
f_current = open(f'{args.siloname}_current.txt', 'w')
f_current.write('# time(s) J_tot J_displ\n')

J_tot, J_displ, G_eff = p3d.compute_current(time)
f_current.write(f'{time:.6e} {J_tot:.6e} {J_displ:.6e}\n')
f_current.flush()

t_start = perf_counter()

for step in range(1, args.n_steps+1):
    print(f'{step:4d} t = {time*1e9:.1f} ns n_streamers = {len(streamers)}')

    # Potentially increase dt when there are no more streamers
    if len(streamers) == 0:
        dt = args.dt * args.dt_factor
    else:
        dt = args.dt

    if args.print_performance:
        t_total = perf_counter() - t_start
        print(f' refinement: {1e2*wct_refinement/t_total:.2f}% '
              f'update_sigma: {1e2*wct_update_sigma/t_total:.2f}% '
              f'poisson: {1e2*wct_poisson/t_total:.2f}% '
              f'output: {1e2*wct_output/t_total:.2f}%')

    new_branches = []
    for s in streamers:
        tau_branch = get_tau_branch(s.R, norm(s.v))

        if np.random.exponential(tau_branch) < dt:
            s.is_branching = True
            s.branching_angle = np.random.uniform(0., args.branch_gamma *
                                                  np.pi/180.)
            s.branching_axis = find_orthogonal_unit_vector(s.v)

            new_branch = copy.copy(s)
            new_branch.branching_angle -= np.pi/2
            new_branch.n_steps = 0
            new_branches.append(new_branch)

    streamers = streamers + new_branches
    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:
        if s.is_branching:
            rot = Rotation.from_rotvec(s.branching_axis * s.branching_angle)
            s.v = rot.apply(s.v)
            s.is_branching = False

        # Take the field ahead of the streamer. This field will tend to bend
        # towards the background electric field.
        v_hat = s.v/norm(s.v)
        E_hat, success = p3d.get_field_vector_at(s.r + 1.5 * s.R * v_hat)

        if not success:
            print('Could not sample E_hat, removing streamer')
            s.keep = False
            continue

        E_hat = E_hat/norm(E_hat)

        # Get samples of |E| ahead of the streamer to determine L_E
        r_tip = s.r + 0.5 * s.R * s.v/norm(s.v)
        z, E, success = p3d.get_var_along_line('E_norm', r_tip, E_hat,
                                               args.L_E_max, 2*args.L_E_max/dz)

        if success:
            L_E_new = model.get_L_E(z, E, N0, dz, prev=s.L_E)
        else:
            # Use previous value
            L_E_new = s.L_E

        if L_E_new < args.L_E_min:
            print('L_E too small, removing streamer')
            s.keep = False
            continue

        if s.n_steps == 0:
            # At the first step, use only the new field
            L_E = L_E_new
        else:
            L_E = args.alpha * L_E_new + (1 - args.alpha) * s.L_E

        s.L_E = L_E
        s.sigma = model.get_sigma(L_E, N0)
        s.v = model.get_velocity(L_E, N0) * E_hat

        dR = min(args.r_scale * model.get_radius(L_E, N0) - s.R,
                 norm(s.v) * dt)
        s.R = s.R + dR
        s.r = s.r + s.v * (dt - 0.99 * dR/norm(s.v))
        s.n_steps += 1

    t0 = perf_counter()
    n_add = p3d.adjust_refinement()
    t1 = perf_counter()
    wct_refinement += t1 - t0

    if args.gas_dynamics:
        p3d.update_gas(dt)

    mlib.update_sigma(3, p3d.update_sigma, streamers, streamers_prev,
                      time, dt, args.channel_update_delay, step == 1)
    t0 = perf_counter()
    wct_update_sigma += t0 - t1

    time += dt
    set_voltage(time, args.phi_bc, voltage_table)
    p3d.solve(dt, args.poisson_rtol, args.poisson_atol)
    t1 = perf_counter()
    wct_poisson += t1 - t0

    J_tot, J_displ, G_eff = p3d.compute_current(time)

    if args.use_circuit:
        V_cap, V_gap = p3d.update_voltage_rc(dt, G_eff)

    f_current.write(f'{time:.6e} {J_tot:.6e} {J_displ:.6e} '
                    f'{V_gap:.6e} {V_cap:.6e}\n')
    f_current.flush()

    # Write output every N steps, with N = args.steps_per_output
    if step % args.steps_per_output == 0:
        i_output = step // args.steps_per_output
        p3d.write_solution(f'{args.siloname}_{i_output:04d}', i_output, time)

    t0 = perf_counter()
    wct_output += t0 - t1

    streamers = [s for s in streamers if s.keep]
