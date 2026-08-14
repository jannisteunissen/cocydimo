module m_solver_lib
#include "cpp_macros.h"

  use m_af_all

  implicit none
  public

  real(dp), parameter :: eps0 = 8.8541878128e-12_dp ! permitivity of vacuum (SI)
  real(dp), parameter :: elem_charge = 1.602176634e-19_dp
  real(dp), parameter :: Da = 1.66053906892e-27_dp ! Dalton (kg)
  real(dp), parameter :: k_b = 1.380649e-23_dp ! Boltzmann constant (J/K)
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: SI_to_Townsend = 1e21_dp ! Convert V/m to Townsend
  real(dp), parameter, public :: Townsend_to_SI = 1e-21_dp ! Convert Townsend to V/m

  type(af_t) :: tree
  type(mg_t) :: mg
  type(mg_t) :: mg_lpl
  integer    :: i_sigma_tot
  integer    :: i_sigma_e
  integer    :: i_sigma_i
  integer    :: i_E_norm
  integer    :: i_E_vec
  integer    :: i_lsf
  integer    :: i_time

  ! How verbose the code is
  integer :: verbose = 0

  ! For simple R-C circuit
  real(dp) :: C_gap = 0.0_dp             ! Geometric gap capacitance [F]
  real(dp) :: applied_voltage = 0.0_dp   ! Gap voltage [V]
  real(dp) :: capacitor_voltage = 0.0_dp ! Capacitor voltage [V]
  real(dp) :: rc_resistance = 0.0_dp ! RC resistance [Ohm]
  real(dp) :: rc_capacitance = 0.0_dp ! RC capacitance [farad]

  ! Maximum electron conductivity. Relevant in regions where the field remains
  ! above the critical field.
  real(dp) :: max_sigma = 5.0_dp

  ! Minimum electron conductivity. This allows the electron conductivity to
  ! grow again at a later time.
  real(dp) :: min_sigma = 1e-9_dp

  ! Effective electron mobility, only used to update ion conductivity
  real(dp) :: mu_electron = 0.04_dp

  ! Effective ion mobility
  real(dp) :: mu_ion = 2e-4_dp

  ! Ion-ion recombination rate constant (m^3/s)
  real(dp) :: k_ion_rec = 1e-13_dp

  ! For gas dynamics
  real(dp) :: gas_gamma = 1.4_dp
  real(dp) :: gas_inv_gamma_m1 = 1/(1.4_dp - 1)

  ! Mean molecular weight of gas
  real(dp) :: gas_inv_molecular_weight = 0 ! will be set later

  ! Inverse of initial gas number density
  real(dp) :: gas_inv_N0 = 300_dp * k_b / 1e5_dp

  ! Number of gas variables
  integer, parameter :: n_gas_vars = 2 + fndims

  ! Density variable (relative index)
  integer, parameter :: i_gas_rho = 1
  ! Index offset for momentum (relative index)
#if fndims == 2
  integer, parameter :: i_gas_mom(fndims) = [2, 3]
#elif fndims == 3
  integer, parameter :: i_gas_mom(fndims) = [2, 3, 4]
#endif
  ! Energy variable (relative index)
  integer, parameter :: i_gas_e = 2 + fndims

  ! Indices of temporal variables in tree data structure
  integer :: i_gas_vars(n_gas_vars) = -1

  ! Indices of fluxes in tree data structure
  integer :: i_gas_fluxes(n_gas_vars) = -1

  ! Index of slow gas heating variable in tree data structure
  integer :: i_gas_slow_heat = -1

  ! Fraction of Joule heating that is immediately converted to gas heating
  real(dp) :: gas_fast_heating_factor = 1.0_dp

  ! Fraction of Joule heating that is slowly converted to gas heating
  real(dp) :: gas_slow_heating_factor = 0.0_dp

  ! Time scale for slow heating, related to V-T relaxation (s)
  real(dp) :: gas_slow_heating_timescale = 20e-6_dp

#if fndims == 2
  ! Names of variables
  character(len=10), parameter :: gas_var_names(n_gas_vars) = [character(len=10) :: &
       "gas_rho", "gas_momx", "gas_momy", "gas_e"]
#elif fndims == 3
  ! Names of variables
  character(len=10), parameter :: gas_var_names(n_gas_vars) = [character(len=10) :: &
       "gas_rho", "gas_momx", "gas_momy", "gas_momz", "gas_e"]
#endif

  ! Electrode parameters
  real(dp) :: rod_r0(fndims), rod_r1(fndims), rod_radius = 0.0_dp

  ! Table with k_eff for updating channel conductivity
  real(dp), allocatable :: k_eff_table(:)
  integer               :: k_eff_table_n_points
  real(dp)              :: k_eff_table_x_min
  real(dp)              :: k_eff_table_inv_fac

  ! For mesh refinement
  real(dp) :: refine_min_dx            = -1.0_dp ! Minimum grid spacing (m)
  real(dp) :: refine_max_dx            = -1.0_dp ! Maximum grid spacing (m)
  real(dp) :: refine_field_threshold   = -1.0_dp ! Refinement field threshold (V/m)
  real(dp) :: derefine_field_threshold = -1.0_dp ! Derefinement field threshold (V/m)
  real(dp) :: derefine_dx              = -1.0_dp ! Only derefine up to this value (m)
  real(dp) :: refine_electrode_max_dx  = -1.0_dp ! Maximum dx around electrode (m)

contains

  subroutine refinement_criterion(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: max_field, dx
    integer                 :: nc

    nc = box%n_cell
    max_field = maxval(box%cc(DTIMES(1:nc), i_E_norm))
    dx = minval(box%dr)

    if (max_field > refine_field_threshold .and. dx > 2 * refine_min_dx) then
       cell_flags = af_do_ref
    else if (iand(box%tag, mg_lsf_box) > 0 .and. &
         dx > refine_electrode_max_dx) then
       cell_flags = af_do_ref
    else if (max_field < derefine_field_threshold .and. dx < derefine_dx) then
       cell_flags = af_rm_ref
    else
       cell_flags = af_keep_ref
    end if

    ! Ensure dx <= refine_max_dx
    if (dx > refine_max_dx) cell_flags = af_do_ref

  end subroutine refinement_criterion

  subroutine set_epsilon_from_sigma(box, dt_vec)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt_vec(:)
    integer                    :: nc

    nc = box%n_cell
    box%cc(DTIMES(1:nc), tree%mg_i_eps) = 1 + (dt_vec(1)/eps0) * &
         box%cc(DTIMES(1:nc), i_sigma_tot)
  end subroutine set_epsilon_from_sigma

  subroutine compute_rhs(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer                   :: lvl, i, id

    call mg_use(tree, mg)

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call mg%box_op(tree%boxes(id), mg%i_rhs, mg)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine compute_rhs

  ! This routine sets boundary conditions for a box
  subroutine sides_bc(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(fndims, box%n_cell**(fndims-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(fndims-1))
    integer, intent(out)    :: bc_type

    if (nb == 2 * fndims - 1) then
       bc_type = af_bc_dirichlet
       bc_val = 0.0_dp
    else if (nb == 2 * fndims) then
       bc_type = af_bc_dirichlet
       bc_val = applied_voltage
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine sides_bc

  ! Compute distance from a line
  pure subroutine dist_vec_line(r, r0, r1, n_dim, dist_vec, dist_line, frac)
    integer, intent(in)   :: n_dim
    real(dp), intent(in)  :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp), intent(out) :: dist_vec(n_dim) !< Distance vector from line segment
    real(dp), intent(out) :: dist_line       !< Distance from line if it were infinite
    real(dp), intent(out) :: frac            !< Fraction [0,1] along line
    real(dp)              :: line_len2
    real(dp), parameter   :: eps = 1e-100_dp

    ! Distance to infinite line
    line_len2 = sum((r1 - r0)**2)
    if (line_len2 > eps) then
       frac = sum((r - r0) * (r1 - r0))/line_len2
    else
       frac = 0.0_dp
    end if

    dist_vec = r - r0 - frac * (r1 - r0)
    dist_line = norm2(dist_vec)

    ! Adjust for finite line
    if (frac <= 0.0_dp) then
       dist_vec = r - r0
    else if (frac >= 1.0_dp) then
       dist_vec = r - r1
    end if
  end subroutine dist_vec_line

  ! Level set function for rod electrode
  real(dp) function rod_lsf(r)
    real(dp), intent(in) :: r(fndims)
    rod_lsf = get_dist_line(r, rod_r0, rod_r1, fndims) - rod_radius
  end function rod_lsf

  subroutine set_lsf_box(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(fndims)

    nc = box%n_cell
    do KJI_DO(0, nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = rod_lsf(rr)
    end do; CLOSE_DO
  end subroutine set_lsf_box

  function get_dist_line(r, r0, r1, n_dim) result(dist)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp)             :: dist, dist_vec(n_dim), dist_line, frac
    call dist_vec_line(r, r0, r1, n_dim, dist_vec, dist_line, frac)
    dist = norm2(dist_vec)
  end function get_dist_line

  !> Boundary condition for radial momentum flux (in axisymmetric coordinates)
  subroutine bc_radial_momentum(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(fndims, box%n_cell**(fndims-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(fndims-1))
    integer, intent(out)    :: bc_type

    if (nb == af_neighb_lowx) then
       ! This will ensure the radial momentum is zero on the axis (by having
       ! ghost values with opposite sign)
       bc_type = af_bc_dirichlet
       bc_val  = 0.0_dp
    else
       bc_type = af_bc_neumann
       bc_val  = 0.0_dp
    end if
  end subroutine bc_radial_momentum

  !> Set initial condition for the gas
  subroutine set_initial_condition_gas(box, arguments)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: arguments(:)
    integer                    :: IJK, nc

    nc = box%n_cell

    ! Initialize Euler variables: density, momentum, energy
    do KJI_DO(0, nc+1)
       box%cc(IJK, i_gas_vars(i_gas_rho)) = arguments(1)
       box%cc(IJK, i_gas_vars(i_gas_mom)) = arguments(2:2+fndims-1)
       box%cc(IJK, i_gas_vars(i_gas_e)) = arguments(2+fndims)
    end do; CLOSE_DO
  end subroutine set_initial_condition_gas

  subroutine add_gas_source_terms(box, dt_vec)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt_vec(:)
    integer                    :: IJK, nc
    real(dp)                   :: dt, J_dot_E
    real(dp)                   :: E_vt_release

    dt = dt_vec(1)
    nc = box%n_cell

    do KJI_DO(1, nc)
       ! Joule heating term is sigma * E**2
       J_dot_E = box%cc(IJK, i_sigma_tot) * box%cc(IJK, i_E_norm)**2 * dt

       ! How much energy is released from slow heating
       E_vt_release = box%cc(IJK, i_gas_slow_heat)/gas_slow_heating_timescale * dt

       box%cc(IJK, i_gas_slow_heat) = box%cc(IJK, i_gas_slow_heat) + &
            gas_slow_heating_factor * J_dot_E - E_vt_release

       box%cc(IJK, i_gas_vars(i_gas_e)) = box%cc(IJK, i_gas_vars(i_gas_e)) + &
            gas_fast_heating_factor * J_dot_E + E_vt_release
    end do; CLOSE_DO

  end subroutine add_gas_source_terms

  subroutine gas_forward_euler(tree, dt, dt_stiff, dt_lim, time, s_deriv, n_prev, &
       s_prev, w_prev, s_out, i_step, n_steps)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    real(dp), intent(in)      :: dt_stiff       !< Time step for stiff terms
    real(dp), intent(inout)   :: dt_lim
    real(dp), intent(in)      :: time
    integer, intent(in)       :: s_deriv
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out
    integer, intent(in)       :: i_step, n_steps
    real(dp)                  :: dummy_dt(0)

    call flux_generic_tree(tree, n_gas_vars, i_gas_vars, s_deriv, i_gas_fluxes, dt_lim, &
         max_wavespeed, get_fluxes, flux_dummy_modify, flux_dummy_line_modify, &
         to_primitive, to_conservative, af_limiter_vanleer_t)

    if (tree%coord_t == af_cyl) then
       call flux_update_densities(tree, dt, n_gas_vars, i_gas_vars, n_gas_vars, &
            i_gas_vars, i_gas_fluxes, s_deriv, n_prev, s_prev, w_prev, s_out, &
            add_geometric_source, 0, dummy_dt)
    else
       call flux_update_densities(tree, dt, n_gas_vars, i_gas_vars, n_gas_vars, &
            i_gas_vars, i_gas_fluxes, s_deriv, n_prev, s_prev, w_prev, s_out, &
            flux_dummy_source, 0, dummy_dt)
    end if

  end subroutine gas_forward_euler

  !> Convert gas variables to primitive form
  subroutine to_primitive(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
    integer                 :: i

    do i = 1, fndims
       u(:, i_gas_mom(i)) = u(:, i_gas_mom(i)) / u(:, i_gas_rho)
    end do
    u(:, i_gas_e) = (gas_gamma-1.0_dp) * (u(:, i_gas_e) - &
         0.5_dp*u(:, i_gas_rho)* sum(u(:, i_gas_mom(:))**2, dim=2))
  end subroutine to_primitive

  !> Convert gas variables to conservative form
  subroutine to_conservative(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
    real(dp)                :: kin_en(n_values)
    real(dp)                :: inv_fac
    integer                 :: i

    ! Compute kinetic energy (0.5 * rho * velocity^2)
    kin_en = 0.5_dp * u(:, i_gas_rho) * sum(u(:, i_gas_mom(:))**2, dim=2)

    ! Compute energy from pressure and kinetic energy
    inv_fac = 1/(gas_gamma - 1.0_dp)
    u(:, i_gas_e) = u(:, i_gas_e) * inv_fac + kin_en

    ! Compute momentum from density and velocity components
    do i = 1, fndims
       u(:, i_gas_mom(i)) = u(:, i_gas_rho) * u(:, i_gas_mom(i))
    end do
  end subroutine to_conservative

  !> Estimate of maximum wavespeed in gas
  subroutine max_wavespeed(n_values, n_var, flux_dim, u, w)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var) !< Primitive variables
    real(dp), intent(out) :: w(n_values) !< Maximum speed
    real(dp)              :: sound_speeds(n_values)

    sound_speeds = sqrt(gas_gamma * u(:, i_gas_e) / u(:, i_gas_rho))
    w = sound_speeds + abs(u(:, i_gas_mom(flux_dim)))
  end subroutine max_wavespeed

  !> Compute fluxes from primitive variables
  subroutine get_fluxes(n_values, n_var, flux_dim, u, flux, box, line_ix, s_deriv)
    integer, intent(in)     :: n_values !< Number of cell faces
    integer, intent(in)     :: n_var    !< Number of variables
    integer, intent(in)     :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(n_values, n_var)
    real(dp), intent(out)   :: flux(n_values, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: line_ix(fndims-1)
    integer, intent(in)     :: s_deriv
    real(dp)                :: E(n_values), inv_fac
    integer                 :: i

    ! Compute left and right flux for conservative variables from the primitive
    ! reconstructed values.

    ! Density flux
    flux(:, i_gas_rho) = u(:, i_gas_rho) * u(:, i_gas_mom(flux_dim))

    ! Momentum flux
    do i = 1, fndims
       flux(:,  i_gas_mom(i)) = u(:, i_gas_rho) * &
            u(:, i_gas_mom(i)) * u(:, i_gas_mom(flux_dim))
    end do

    ! Add pressure term
    flux(:, i_gas_mom(flux_dim)) = flux(:, i_gas_mom(flux_dim)) + u(:, i_gas_e)

    ! Compute energy
    inv_fac = 1/(gas_gamma-1.0_dp)
    E = u(:, i_gas_e) * inv_fac + 0.5_dp * u(:, i_gas_rho) * &
         sum(u(:, i_gas_mom(:))**2, dim=2)

    ! Energy flux
    flux(:, i_gas_e) = u(:, i_gas_mom(flux_dim)) * (E + u(:, i_gas_e))

  end subroutine get_fluxes

  !> Geometric source term for axisymmetric simulations of gas dynamics
  subroutine add_geometric_source(box, dt, n_vars, i_cc, s_deriv, s_out, &
       n_dt, dt_lim, mask)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: n_vars
    integer, intent(in)        :: i_cc(n_vars)
    integer, intent(in)        :: s_deriv
    integer, intent(in)        :: s_out
    logical, intent(in)        :: mask(DTIMES(box%n_cell))
    integer, intent(in)        :: n_dt
    real(dp), intent(inout)    :: dt_lim(n_dt)

#if fndims == 2
    real(dp)                   :: pressure(DTIMES(box%n_cell))
    real(dp)                   :: inv_radius
    integer                    :: nc, i

    nc = box%n_cell
    pressure = get_pressure(box, s_deriv)

    do i = 1, nc
       inv_radius = 1/af_cyl_radius_cc(box, i)
       where (mask(i, :))
          box%cc(i, 1:nc, i_cc(i_gas_mom(1))+s_out) = &
               box%cc(i, 1:nc, i_cc(i_gas_mom(1))+s_out) + dt * &
               pressure(i, :) * inv_radius
       end where
    end do
#endif
  end subroutine add_geometric_source

#if fndims == 2
  pure function get_pressure(box, s_in) result(pressure)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: s_in
    real(dp)                :: pressure(DTIMES(box%n_cell))
    integer                 :: nc

    nc = box%n_cell
    pressure = (gas_gamma-1.0_dp) * (&
         box%cc(DTIMES(1:nc), i_gas_vars(i_gas_e)+s_in) - 0.5_dp * &
         sum(box%cc(DTIMES(1:nc), i_gas_vars(i_gas_mom)+s_in)**2, dim=fndims+1) / &
         box%cc(DTIMES(1:nc), i_gas_vars(i_gas_rho)+s_in))
  end function get_pressure
#endif

  !> Write primitive variables to output: velocities, pressure, temperature
  subroutine write_gas_primitive(box, new_vars, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: n_var
    integer                 :: i
    real(dp)                :: new_vars(DTIMES(0:box%n_cell+1), n_var)
    real(dp)                :: inv_rho(DTIMES(0:box%n_cell+1))

    if (n_var /= fndims + 3) error stop "Invalid value for n_var"

    inv_rho = 1/box%cc(DTIMES(:), i_gas_vars(i_gas_rho))

    ! Velocity
    do i = 1, fndims
       new_vars(DTIMES(:), i) = box%cc(DTIMES(:), i_gas_vars(i_gas_mom(i))) * inv_rho
    end do

    ! Pressure = (gamma - 1) * (gas_e - kinetic energy)
    new_vars(DTIMES(:), fndims+1) = (gas_gamma-1.0_dp) * &
         (box%cc(DTIMES(:), i_gas_vars(i_gas_e)) - &
         0.5_dp * inv_rho * &
         sum(box%cc(DTIMES(:), i_gas_vars(i_gas_mom(:)))**2, dim=fndims+1))

    ! Gas number density N = rho / molecular weight
    new_vars(DTIMES(:), fndims+2) = box%cc(DTIMES(:), i_gas_vars(i_gas_rho)) * &
         gas_inv_molecular_weight

    ! Temperature = P / (k_b * N)
    new_vars(DTIMES(:), fndims+3) = new_vars(DTIMES(:), fndims+1) / &
         (k_b * new_vars(DTIMES(:), fndims+2))

  end subroutine write_gas_primitive

  !> Compute total field energy in Joule, defined as the volume integral over
  !> 1/2 * epsilon * E^2
  subroutine compute_field_energy(field_energy)
    real(dp), intent(out) :: field_energy

    call af_reduction(tree, field_energy_box, reduce_sum, 0.0_dp, field_energy)
  end subroutine compute_field_energy

  !> Compute the integral over J.E over the whole domain
  subroutine compute_JdotE_integral(JdotE)
    real(dp), intent(out) :: JdotE

    call af_reduction(tree, sum_JdotE_box, reduce_sum, 0.0_dp, JdotE)
  end subroutine compute_JdotE_integral

  !> Get the electrostatic field energy in a box
  real(dp) function field_energy_box(box)
    type(box_t), intent(in) :: box
#if fndims == 2
    integer                 :: i
    real(dp), parameter     :: twopi = 2 * acos(-1.0_dp)
#endif
    real(dp)                :: w(DTIMES(box%n_cell))
    integer                 :: nc

    nc = box%n_cell
    w = 0.5_dp * eps0 * product(box%dr)

#if fndims == 2
    if (box%coord_t == af_cyl) then
       ! Weight by 2 * pi * r
       do i = 1, nc
          w(i, :) = w(i, :) * twopi * af_cyl_radius_cc(box, i)
       end do
    end if
#endif

    field_energy_box = sum(w * box%cc(DTIMES(1:nc), i_E_norm)**2)
  end function field_energy_box

  !> Get the sum of J dot E in a box
  real(dp) function sum_JdotE_box(box)
    type(box_t), intent(in) :: box
#if fndims == 2
    integer                 :: i
    real(dp), parameter     :: twopi = 2 * acos(-1.0_dp)
#endif
    real(dp)                :: w(DTIMES(box%n_cell))
    integer                 :: nc

    nc = box%n_cell
    w = product(box%dr)

#if fndims == 2
    if (box%coord_t == af_cyl) then
       ! Weight by 2 * pi * r
       do i = 1, nc
          w(i, :) = w(i, :) * twopi * af_cyl_radius_cc(box, i)
       end do
    end if
#endif

    sum_JdotE_box = sum(w * box%cc(DTIMES(1:nc), i_sigma_tot) * &
         box%cc(DTIMES(1:nc), i_E_norm)**2)
  end function sum_JdotE_box

  real(dp) function reduce_sum(a, b)
    real(dp), intent(in) :: a, b
    reduce_sum = a + b
  end function reduce_sum

end module m_solver_lib
