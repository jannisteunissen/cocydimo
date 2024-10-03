#!/usr/bin/env python3

import numpy as np

from poisson_3d import m_solver as m_solver_3d
from poisson_2d import m_solver as m_solver_2d

domain_size = [1.0, 1.0]
grid_size = [256, 256]
box_size = 256
phi_bc = 1.0

rhs = np.zeros(grid_size)
sigma = np.zeros(grid_size)

# Define electrode
rod_r0 = np.array([0.0, 0.0]) * domain_size
rod_r1 = np.array([0.0, 0.15]) * domain_size
rod_radius = 0.1

m_solver_2d.set_rod_electrode(rod_r0, rod_r1, rod_radius)
m_solver_2d.initialize(domain_size, grid_size, box_size, phi_bc)

m_solver_2d.set_rhs_and_sigma(rhs, sigma)
m_solver_2d.solve(0.0)
m_solver_2d.write_solution('output/test_poisson_solver_2d', 0, 0.)

domain_size = [1.0, 1.0, 1.0]
grid_size = [128, 128, 128]
box_size = 32
phi_bc = 1.0

# Define electrode
rod_r0 = np.array([0.5, 0.5, 0.0]) * domain_size
rod_r1 = np.array([0.5, 0.5, 0.15]) * domain_size
rod_radius = 0.1

m_solver_3d.set_rod_electrode(rod_r0, rod_r1, rod_radius)
m_solver_3d.initialize(domain_size, grid_size, box_size, phi_bc)

m_solver_3d.solve(0.0)
m_solver_3d.write_solution('output/test_poisson_solver_3d', 0, 0.)
