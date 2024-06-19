program test_poisson_solver
  use m_solver
  implicit none

  real(dp) :: domain_size(2) = [1.0_dp, 1.0_dp]
  integer :: grid_size(2) = [256, 256]
  integer :: box_size = 8

  rod_r0 = [0.0_dp, 0.0_dp] * domain_size
  rod_r1 = [0.0_dp, 0.15_dp] * domain_size
  rod_radius = 0.1_dp
  phi_bc = 1.0_dp

  call set_rod_electrode(rod_r0, rod_r1, rod_radius)
  call initialize(domain_size, grid_size, box_size, phi_bc)
  call solve(0.0_dp)
  call write_solution('output/test_poisson_solver_fortran')

end program test_poisson_solver
