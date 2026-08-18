module m_solver
#include "cpp_macros.h"

  use m_solver_lib

  implicit none
  private

  public :: initialize_domain
  public :: use_uniform_grid
  public :: set_refinement
  public :: adjust_refinement
  public :: set_rod_electrode
  public :: set_voltage
  public :: update_voltage_rc
  public :: update_sigma
  public :: store_k_eff
  public :: store_parameters
  public :: get_finest_grid_spacing
  public :: get_max_field_location
  public :: get_field_vector_at
  public :: get_var_along_line
  public :: solve
  public :: write_solution
  public :: set_gas
  public :: update_gas
  public :: get_gas_number_density_at
  public :: compute_current

contains

  ! Initialize the computational domain
  subroutine initialize_domain(domain_len, coarse_grid_size, box_size, &
       voltage, mem_limit_gb, write_eps, write_time, write_rhs, &
       gas_dynamics)
    real(dp), intent(in) :: domain_len(fndims)       ! Domain size (m)
    integer, intent(in)  :: coarse_grid_size(fndims) ! Coarse grid size
    integer, intent(in)  :: box_size                 ! Size of grid boxes
    real(dp), intent(in) :: voltage                  ! Applied voltage (V)
    real(dp), intent(in) :: mem_limit_gb             ! Memory limit (GB)
    logical, intent(in)  :: write_eps                ! Write epsilon to output
    logical, intent(in)  :: write_time               ! Write time to output
    logical, intent(in)  :: write_rhs                ! Write rhs to output
    logical, intent(in)  :: gas_dynamics             ! Simulate gas dynamics
    integer              :: coord_t, n

    if (verbose > 0) print *, "log: initialize_domain()"

    coord_t = af_xyz
    if (fndims == 2) coord_t = af_cyl

    applied_voltage = voltage
    capacitor_voltage = voltage

    call af_add_cc_variable(tree, "phi", ix=mg%i_phi)
    call af_add_cc_variable(tree, "rhs", ix=mg%i_rhs, write_out=write_rhs)
    call af_add_cc_variable(tree, "tmp", ix=mg%i_tmp, write_out=.false.)
    call af_add_cc_variable(tree, "eps", ix=tree%mg_i_eps, write_out=write_eps)
    call af_add_cc_variable(tree, "sigma_tot", ix=i_sigma_tot)
    call af_add_cc_variable(tree, "sigma_e", ix=i_sigma_e)
    call af_add_cc_variable(tree, "sigma_i", ix=i_sigma_i)
    call af_add_cc_variable(tree, "electric_fld", ix=i_E_norm)
    call af_add_cc_variable(tree, "time", ix=i_time, write_out=write_time)
    call af_add_fc_variable(tree, "E_vec", ix=i_E_vec)

    call af_set_cc_methods(tree, tree%mg_i_eps, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_sigma_tot, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_sigma_e, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_sigma_i, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_E_norm, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_time, af_bc_neumann_zero)

    if (gas_dynamics) then
       call af_add_cc_variable(tree, "slow_heat", ix=i_gas_slow_heat)

       do n = 1, n_gas_vars
          call af_add_cc_variable(tree, gas_var_names(n), ix=i_gas_vars(n), &
               n_copies=2, write_out=.false.)
          call af_add_fc_variable(tree, "flux", ix=i_gas_fluxes(n))

          if (coord_t == af_cyl .and. n == i_gas_mom(1)) then
             call af_set_cc_methods(tree, i_gas_vars(n), bc_radial_momentum)
          else
             call af_set_cc_methods(tree, i_gas_vars(n), af_bc_neumann_zero)
          end if
       end do
    end if

    if (rod_radius > 0) then
       call af_add_cc_variable(tree, "lsf", ix=i_lsf)

       mg%lsf_boundary_value = 0.0_dp ! Electrode is grounded
       mg%lsf => rod_lsf
       mg%lsf_dist => mg_lsf_dist_gss
       mg%lsf_length_scale = rod_radius

       tree%mg_i_lsf = i_lsf
       call af_set_cc_methods(tree, i_lsf, funcval=set_lsf_box)
    end if

    call af_init(tree, box_size, domain_len, coarse_grid_size, &
         coord=coord_t, mem_limit_gb=mem_limit_gb)

    mg%sides_bc => sides_bc ! Method for boundary conditions

    ! Create a copy of the operator but without the variable coefficient
    mg_lpl = mg
    mg_lpl%operator_mask = mg_normal_box + mg_lsf_box

    ! Estimate initial gap capacitance, used in RC model
    if (coord_t == af_cyl) then
       C_gap = eps0 * pi * domain_len(1)**2 / domain_len(fndims)
    else
       C_gap = eps0 * product(domain_len(1:2)) / domain_len(fndims)
    end if

    if (verbose > 0) print *, "log: initialize_domain() done"

  end subroutine initialize_domain

  !> Set initial state for gas
  subroutine set_gas(pressure, temperature, mean_molecular_weight, gamma, &
       f_fast_heat, f_slow_heat, tau_slow_heat)
    real(dp), intent(in) :: pressure ! in bar
    real(dp), intent(in) :: temperature ! in Kelvin
    real(dp), intent(in) :: mean_molecular_weight ! in Dalton
    real(dp), intent(in) :: gamma ! Adiabatic index
    real(dp), intent(in) :: f_fast_heat ! Fast heating factor
    real(dp), intent(in) :: f_slow_heat ! Slow heating factor
    real(dp), intent(in) :: tau_slow_heat ! Slow heating time scale (s)

    real(dp)            :: N0, rho, momentum(fndims), energy

    if (.not. allocated(tree%boxes)) &
         error stop "Call initialize_domain before set_gas"
    if (i_gas_vars(1) == -1) &
         error stop "Gas was not initialized when calling initialize_domain"

    ! Ideal gas law (approximation)
    N0 = 1e5_dp * pressure / (k_b * temperature)

    gas_gamma = gamma
    gas_inv_gamma_m1 = 1/(gamma - 1)
    gas_fast_heating_factor    = f_fast_heat
    gas_slow_heating_factor    = f_slow_heat
    gas_slow_heating_timescale = tau_slow_heat
    gas_inv_molecular_weight   = 1/(mean_molecular_weight * Da)
    gas_inv_N0                 = 1/N0

    ! Set initial gas density
    rho = N0 * mean_molecular_weight * Da

    ! Initial momentum
    momentum = 0.0_dp

    ! Initial energy
    energy = pressure * 1e5_dp / (gamma - 1)

    ! Set initial density, momentum and energy in domain
    call af_loop_box_arg(tree, set_initial_condition_gas, [rho, momentum, energy])
  end subroutine set_gas

  subroutine set_voltage(voltage)
    real(dp), intent(in) :: voltage

    applied_voltage = voltage
  end subroutine set_voltage

  !> Get gas number density at some location
  subroutine get_gas_number_density_at(r, N, success)
    real(dp), intent(in)  :: r(fndims)
    real(dp), intent(out) :: N
    logical, intent(out)  :: success
    real(dp)              :: tmp(1)

    tmp = af_interp1(tree, r, [i_gas_vars(i_gas_rho)], success)

    if (success) then
       N = tmp(1) * gas_inv_molecular_weight
    else
       N = 0.0_dp
    end if
  end subroutine get_gas_number_density_at

  !> Add source terms from the discharge in the Euler equations
  subroutine update_gas(dt, max_dt)
    real(dp), intent(in)  :: dt
    real(dp), intent(out) :: max_dt     ! Limit on dt due to gas dynamics
    real(dp)              :: time_dummy ! Unused

    time_dummy = 0.0_dp
    call af_loop_box_arg(tree, add_gas_source_terms, [dt], .true.)
    call af_advance(tree, dt, max_dt, time_dummy, i_gas_vars, &
         af_heuns_method, gas_forward_euler)
  end subroutine update_gas

  ! Perform uniform initial refinement of the domain
  subroutine use_uniform_grid(uniform_grid_size)
    integer, intent(in) :: uniform_grid_size(fndims)
    integer             :: max_lvl

    ! Refine uniformly
    max_lvl = nint(log(uniform_grid_size(1) / real(tree%coarse_grid_size(1), dp)) / &
         log(2.0_dp)) + 1

    if (any(tree%coarse_grid_size * 2**(max_lvl-1) /= uniform_grid_size)) &
         error stop "Incompatible grid size"
    call af_refine_up_to_lvl(tree, max_lvl)
  end subroutine use_uniform_grid

  ! Set initial grid refinement
  subroutine set_refinement(refine_field, derefine_field, &
       min_dx, max_dx, electrode_max_dx, derefine_levels, max_head_dist, rtol, atol)
    real(dp), intent(in) :: refine_field     ! Refine when the field is above this value
    real(dp), intent(in) :: derefine_field   ! Derefine when the field is below this value
    real(dp), intent(in) :: min_dx           ! Minimum allowed grid spacing
    real(dp), intent(in) :: max_dx           ! Maximum allowed grid spacing
    real(dp), intent(in) :: electrode_max_dx ! Maximum grid spacing around electrode
    integer, intent(in)  :: derefine_levels  ! How many levels can be derefined
    real(dp), intent(in) :: max_head_dist    ! Max. distance from head for refinement
    real(dp), intent(in) :: rtol             ! Relative tolerance for Poisson solver
    real(dp), intent(in) :: atol             ! Absolute tolerance for Poisson solver
    integer              :: n_add, n, n_its
    real(dp)             :: residu

    if (verbose > 0) print *, "log: set_refinement()"

    refine_field_threshold   = refine_field
    derefine_field_threshold = derefine_field
    refine_min_dx            = min_dx
    refine_max_dx            = max_dx
    refine_electrode_max_dx  = electrode_max_dx
    ! Finest dx is between min_dx and 2*min_dx; only derefine when dx < derefine_dx
    derefine_dx              = min_dx * 2**derefine_levels
    refine_max_distance_head = max_head_dist

    call solve(0.0_dp, rtol, atol, n_its, residu)

    do n = 1, 20
       if (verbose > 0) print *, "log: iteration", n
       call adjust_refinement(n_add)
       if (n_add == 0) exit

       ! Reset r.h.s. since we are not advancing in time
       call af_tree_clear_cc(tree, mg%i_rhs)
       call solve(0.0_dp, rtol, atol, n_its, residu)
    end do

    if (verbose > 0) print *, "log: set_refinement() done"
  end subroutine set_refinement

  ! Update the refinement of the mesh. Changes the local refinement by at most
  ! one level at a time, so multiple calls might be needed.
  subroutine adjust_refinement(n_add)
    integer, intent(out) :: n_add
    type(ref_info_t) :: refine_info

    if (verbose > 0) print *, "log: adjust_refinement()"

    ! Restrict species, for the ghost cells near refinement boundaries
    call af_restrict_tree(tree, [i_sigma_e, i_sigma_i])
    call af_gc_tree(tree, [i_sigma_e, i_sigma_i])

    call af_adjust_refinement(tree, refinement_criterion, refine_info, 0)
    n_add = refine_info%n_add

    if (verbose > 0) print *, "log: adjust_refinement() done"
  end subroutine adjust_refinement

  ! Specify geometry of rod electrode
  subroutine set_rod_electrode(r0, r1, radius)
    real(dp), intent(in) :: r0(fndims), r1(fndims), radius

    if (allocated(tree%boxes)) error stop "Set electrode before initialization"
    rod_r0 = r0
    rod_r1 = r1
    rod_radius = radius
  end subroutine set_rod_electrode

  ! Update sigma (conductivity)
  subroutine update_sigma(n_in, r0c, r1c, sigma0, sigma1, radius0, radius1, &
       t, dt, channel_delay, first_step, n_streamers)
    integer, intent(in)  :: n_in
    real(dp), intent(in) :: r0c(n_in, fndims), r1c(n_in, fndims)
    real(dp), intent(in) :: sigma0(n_in), sigma1(n_in)
    real(dp), intent(in) :: radius0(n_in), radius1(n_in)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: channel_delay
    logical, intent(in)  :: first_step
    !> Limit sigma to this value when updating channel conductivity
    integer, intent(in)  :: n_streamers
    integer              :: lvl, n, id, IJK, nc, ix, jx
    real(dp)             :: r(fndims), dist_vec(fndims), r_dist, frac, fld_Td
    real(dp)             :: k_eff, dsigma, box_rmax(fndims), length, radius
    real(dp)             :: mu_rel, ion_fac
    real(dp)             :: r_min(fndims, n_streamers)
    real(dp)             :: r_max(fndims, n_streamers)
    real(dp)             :: r0(fndims, n_in), r1(fndims, n_in)
    integer              :: n_in_box, ix_in_box(n_streamers)

    if (verbose > 0) print *, "log: update_sigma()"

    ! Store transposed arrays for better memory access
    r0(:, :) = transpose(r0c(1:n_streamers, :))
    r1(:, :) = transpose(r1c(1:n_streamers, :))

    ! Store streamer information for refinement
    if (n_streamers > max_streamers) error stop "Increase max_streamers"
    global_n_streamers = n_streamers
    global_r_heads(:, 1:n_streamers) = r1
    global_time = t

    nc = tree%n_cell
    ion_fac = elem_charge * mu_ion
    mu_rel = mu_ion / mu_electron

    ! Determine the extent of channels, with some margin
    do ix = 1, n_streamers
       length = norm2(r1(:, ix) - r0(:, ix))
       radius = max(radius0(ix), radius1(ix))
       r_min(:, ix) = min(r0(:, ix), r1(:, ix)) - radius - 0.5_dp * length
       r_max(:, ix) = max(r0(:, ix), r1(:, ix)) + radius + 0.5_dp * length
    end do

    if (.not. allocated(k_eff_table)) error stop "Call store_k_eff first"

    !$omp parallel private(lvl, n, id, IJK, r, dist_vec, r_dist, &
    !$omp &frac, ix, k_eff, dsigma, box_rmax, n_in_box, ix_in_box, jx, fld_Td)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)

          associate (box => tree%boxes(id))

            ! Determine which channels can possibly lie in the box
            n_in_box = 0
            box_rmax = box%r_min + nc * box%dr
            do ix = 1, n_streamers
               if (all(r_max(:, ix) >= box%r_min .and. &
                    r_min(:, ix) <= box_rmax)) then
                  n_in_box = n_in_box + 1
                  ix_in_box(n_in_box) = ix
               end if
            end do

            do KJI_DO(1, nc)
               if (rod_radius > 0) then
                  if (box%cc(IJK, i_lsf) < 0.0_dp) cycle
               end if

               r = af_r_cc(box, [IJK])

               do jx = 1, n_in_box
                  ix = ix_in_box(jx)
                  call dist_vec_line(r, r0(:, ix), r1(:, ix), fndims, &
                       dist_vec, r_dist, frac)

                  ! Exclude semi-sphere of previous point
                  if (norm2(dist_vec) <= radius1(ix) .and. (first_step .or. &
                       (frac >= 0 .and. norm2(r0(:, ix) - r) > radius0(ix)))) then

                     call get_sigma_profile(r_dist, radius0(ix), radius1(ix), frac, &
                          sigma0(ix), sigma1(ix), dsigma)

                     box%cc(IJK, i_sigma_e) = box%cc(IJK, i_sigma_e) + dsigma
                     box%cc(IJK, i_sigma_i) = box%cc(IJK, i_sigma_i) + dsigma * mu_rel
                     box%cc(IJK, i_time) = t
                  end if
               end do

               ! Update channel electron and ion conductivity, but only where
               ! the channel has already existed for some time
               if (box%cc(IJK, i_time) < t - channel_delay .and. &
                    box%cc(IJK, i_sigma_e) > 0.0_dp) then
                  if (i_gas_vars(1) /= -1) then
                     fld_Td = box%cc(IJK, i_E_norm) / (box%cc(IJK, i_gas_vars(i_gas_rho)) * &
                          gas_inv_molecular_weight) * SI_to_Townsend
                  else
                     fld_Td = box%cc(IJK, i_E_norm) * gas_inv_N0 * SI_to_Townsend
                  end if

                  ! TODO: include scaling of k_eff with gas density. It would
                  ! also be good to split into ionization, attachment and
                  ! recombination (to better update ion conductivity)
                  call get_k_eff(fld_Td, k_eff)

                  ! Electron conductivity change, using analytic expression
                  ! for integral. Limit growth factor per time step to prevent
                  ! instabilities.
                  dsigma = min(2.0_dp, (exp(dt * k_eff) - 1.0_dp)) * box%cc(IJK, i_sigma_e)

                  ! Limit conductivity to at most max_sigma. Relevant in
                  ! regions where the field remains above the critical field.
                  dsigma = min(dsigma, max_sigma - box%cc(IJK, i_sigma_e))

                  ! Limit conductivity to at least min_sigma. This allows the
                  ! electron conductivity to grow at a later time.
                  dsigma = max(dsigma, min_sigma - box%cc(IJK, i_sigma_e))

                  box%cc(IJK, i_sigma_e) = box%cc(IJK, i_sigma_e) + dsigma

                  ! Ion conductivity change. Any change in electron
                  ! conductivity produces net ions. TODO: we could separate
                  ! out attachment and impact ionization.
                  box%cc(IJK, i_sigma_i) = box%cc(IJK, i_sigma_i) + abs(dsigma) * mu_rel

                  ! Ion recombination
                  ! n_i(t+dt) = 1/(dt * k_ion_rec + 1/n_i(t))
                  ! n_i(t) = sigma_i / (e * mu_ion) = sigma_i / ion_fac
                  box%cc(IJK, i_sigma_i) = ion_fac/(dt * k_ion_rec + &
                       ion_fac/box%cc(IJK, i_sigma_i))
               end if

               ! Set total sigma
               box%cc(IJK, i_sigma_tot) = box%cc(IJK, i_sigma_e) + box%cc(IJK, i_sigma_i)

            end do; CLOSE_DO
          end associate
       end do
       !$omp end do
    end do
    !$omp end parallel

    if (verbose > 0) print *, "log: update_sigma() done"

  end subroutine update_sigma

  subroutine get_sigma_profile(r_dist, radius0, radius1, z_frac, s0, s1, dsigma)
    real(dp), intent(in)  :: r_dist ! Radial distance
    real(dp), intent(in)  :: radius0, radius1 ! Radius at s0 and s1
    real(dp), intent(in)  :: z_frac ! Relative z-coordinate, 0.0 at s0, 1.0 at s1
    real(dp), intent(in)  :: s0, s1 ! Value at start and end point
    real(dp), intent(out) :: dsigma
    real(dp)              :: radius, fac_r, frac_bnd
    real(dp), parameter   :: pi = acos(-1.0_dp)

    frac_bnd = max(min(z_frac, 1.0_dp), 0.0_dp)
    radius = radius1

    ! 1-(r/R)^2 profile in radial direction
    fac_r = 2 * max(0.0_dp, 1 - (r_dist/radius)**2) / (pi * radius**2)

    dsigma = fac_r * (frac_bnd * s1 + (1-frac_bnd) * s0)
  end subroutine get_sigma_profile

  ! Linearly interpolate tabulated data for effective ionization rate
  subroutine get_k_eff(fld_Td, k_eff)
    real(dp), intent(in)  :: fld_Td
    real(dp), intent(out) :: k_eff
    real(dp)              :: frac, low_frac
    integer               :: low_ix

    frac = (fld_Td - k_eff_table_x_min) * k_eff_table_inv_fac

    ! Check bounds
    if (frac <= 0) then
       low_ix   = 1
       low_frac = 1
    else if (frac >= k_eff_table_n_points - 1) then
       low_ix   = k_eff_table_n_points - 1
       low_frac = 0
    else
       low_ix   = ceiling(frac)
       low_frac = low_ix - frac
    end if

    k_eff = low_frac * k_eff_table(low_ix) + &
         (1-low_frac) * k_eff_table(low_ix+1)

  end subroutine get_k_eff

  ! Store tabulated data for effective ionization rate
  subroutine store_k_eff(Td_min, Td_max, n_points, k_eff)
    real(dp), intent(in) :: Td_min, Td_max
    integer, intent(in)  :: n_points
    real(dp), intent(in) :: k_eff(n_points)

    allocate(k_eff_table(n_points))
    k_eff_table(:) = k_eff
    k_eff_table_n_points = n_points
    k_eff_table_x_min = Td_min
    k_eff_table_inv_fac = (n_points - 1)/(Td_max - Td_min)
  end subroutine store_k_eff

  ! Store parameters for the model
  subroutine store_parameters(min_sigma_arg, max_sigma_arg, mu_electron_arg, &
       mu_ion_arg, k_ion_rec_arg, resistance, capacitance, verbose_arg)
    real(dp), intent(in) :: min_sigma_arg
    real(dp), intent(in) :: max_sigma_arg
    real(dp), intent(in) :: mu_electron_arg
    real(dp), intent(in) :: mu_ion_arg
    real(dp), intent(in) :: k_ion_rec_arg
    real(dp), intent(in) :: resistance
    real(dp), intent(in) :: capacitance
    integer, intent(in)  :: verbose_arg

    min_sigma = min_sigma_arg
    max_sigma = max_sigma_arg
    mu_electron = mu_electron_arg
    mu_ion = mu_ion_arg
    k_ion_rec = k_ion_rec_arg

    rc_resistance = resistance
    rc_capacitance = capacitance

    verbose = verbose_arg
  end subroutine store_parameters

  ! Get the finest grid spacing of the mesh
  subroutine get_finest_grid_spacing(dx_min)
    real(dp), intent(out) :: dx_min
    dx_min = af_min_dr(tree)
  end subroutine get_finest_grid_spacing

  ! Get the maximum electric field and its location
  subroutine get_max_field_location(E_max_vec, r)
    real(dp), intent(out) :: E_max_vec(fndims)
    real(dp), intent(out) :: r(fndims)
    real(dp)              :: E_max_norm
    type(af_loc_t)        :: loc
    logical               :: success

    call af_tree_max_cc(tree, i_E_norm, E_max_norm, loc)
    r = af_r_loc(tree, loc)
    call get_field_vector_at(r, E_max_vec, success)
  end subroutine get_max_field_location

  ! Get the electric field vector at a location
  subroutine get_field_vector_at(r, E_vec, success)
    real(dp), intent(in)  :: r(fndims)
    real(dp), intent(out) :: E_vec(fndims)
    logical, intent(out)  :: success

    E_vec = af_interp1_fc(tree, r, i_E_vec, success)
  end subroutine get_field_vector_at

  ! Trace the value of a variable along a line
  subroutine get_var_along_line(varname, r0, direction, length, n_steps, &
       z_line, line, success)
    character(len=*), intent(in) :: varname
    real(dp), intent(in)         :: r0(fndims), direction(fndims), length
    integer, intent(in)          :: n_steps
    real(dp), intent(out)        :: z_line(n_steps)
    real(dp), intent(out)        :: line(n_steps)
    logical, intent(out)         :: success
    real(dp)                     :: r(fndims), dr(fndims)
    integer                      :: n, i_var

    if (n_steps <= 1) error stop "n_steps should be at least 2"

    select case (varname)
    case ('sigma')
       i_var = i_sigma_tot
    case ('phi')
       i_var = mg%i_phi
    case ('E_norm')
       i_var = i_E_norm
    case default
       error stop 'Unknown variable'
    end select

    r = r0
    dr = length * direction / (norm2(direction) * (n_steps - 1))
    success = .false.

    do n = 1, n_steps
       line(n:n) = af_interp1(tree, r, [i_var], success)
       if (.not. success) return
       z_line(n) = (n-1) * norm2(dr)
       r = r + dr
    end do
  end subroutine get_var_along_line

  ! Compute new potential for a given time step using the current sigma
  subroutine solve(dt, rtol, atol, n_iterations, residu)
    real(dp), intent(in)  :: dt
    real(dp), intent(in)  :: rtol
    real(dp), intent(in)  :: atol
    integer, intent(out)  :: n_iterations
    real(dp), intent(out) :: residu
    integer, parameter    :: max_iterations = 100
    real(dp)              :: prev_residu, max_rhs, initial_residu
    logical               :: converged

    if (verbose > 0) print *, "log: solve() - n_cells = ", &
         af_num_leaves_used(tree) * real(tree%n_cell)**3
    call af_loop_box_arg(tree, set_epsilon_from_sigma, [dt], leaves_only=.true.)
    call af_restrict_tree(tree, [tree%mg_i_eps])
    call af_gc_tree(tree, [tree%mg_i_eps], corners=.false.)

    if (.not. mg%initialized) then
       call mg_init(tree, mg)
       call mg_init(tree, mg_lpl)
    else
       call mg_update_operator_stencil(tree, mg, .false., .true.)
    end if

    call af_tree_maxabs_cc(tree, mg%i_rhs, max_rhs)
    prev_residu = huge(1.0_dp)
    initial_residu = huge(1.0_dp)
    converged = .false.

    do n_iterations = 1, max_iterations
       call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=.true.)
       if (verbose > 1) print *, "log: mg_fas_fmg done"

       call af_tree_maxabs_cc(tree, mg%i_tmp, residu)
       if (n_iterations == 1) initial_residu = residu

       if (verbose > 0) print *, "log: iteration = ", &
            n_iterations, "residu = ", residu

       converged = (residu < max(atol, rtol * max_rhs))
       if (converged .or. residu > 1e2_dp * initial_residu) exit
       prev_residu = residu
    end do

    if (.not. converged) then
       print *, "Multigrid residual:     ", residu
       print *, "Multigrid n_iterations: ", n_iterations
       error stop "the multigrid solve did not converge, reduce dt?"
    end if

    ! Compute new rhs with standard Laplace operator
    call compute_rhs(tree, mg_lpl)

    ! Compute electric field with standard Laplace operator
    call mg_compute_phi_gradient(tree, mg_lpl, i_E_vec, -1.0_dp, i_E_norm)
    call af_gc_tree(tree, [i_E_norm])
    if (verbose > 0) print *, "log: solve() done"

  end subroutine solve

  ! Write a silo file
  subroutine write_solution(fname, i_cycle, time)
    character(len=*), intent(in) :: fname
    integer, intent(in)          :: i_cycle
    real(dp), intent(in)         :: time
    character(len=20)            :: gas_primitive_names(fndims+3)

    ! Ensure valid ghost cells
    call af_restrict_tree(tree, [i_sigma_e, i_sigma_i, i_sigma_tot])
    call af_gc_tree(tree, [i_sigma_e, i_sigma_i, i_sigma_tot])

    if (i_gas_vars(1) == -1) then
       ! No gas output
       call af_write_silo(tree, trim(fname), i_cycle, time)
    else
       ! Store extra variables for gas
       gas_primitive_names(1) = "gas_vx"
       gas_primitive_names(2) = "gas_vy"
#if fndims == 3
       gas_primitive_names(3) = "gas_vz"
#endif
       gas_primitive_names(fndims+1) = "gas_p"
       gas_primitive_names(fndims+2) = "gas_N"
       gas_primitive_names(fndims+3) = "gas_T"
       call af_write_silo(tree, trim(fname), i_cycle, time, &
            add_vars=write_gas_primitive, add_names=gas_primitive_names)
    end if
  end subroutine write_solution

  !> Estimate electric current according to Sato's equation V*I = sum(J.E),
  !> where J includes both the conduction current and the displacement
  !> current, see 10.1088/0022-3727/32/5/005.
  !> The latter is computed through the field energy
  subroutine compute_current(time, J_tot, J_displ, gap_conductance)
    real(dp), intent(in)  :: time
    real(dp), intent(out) :: J_tot   ! Total current
    real(dp), intent(out) :: J_displ ! Displacement current
    real(dp), intent(out) :: gap_conductance ! Effective gap conductance
    real(dp)              :: energy_deriv, JdotE_integral, new_field_energy
    logical, save         :: first_call        = .true.
    real(dp), save        :: prev_time         = 0.0_dp
    real(dp), save        :: prev_field_energy = 0.0_dp

    call compute_field_energy(new_field_energy)
    call compute_JdotE_integral(JdotE_integral)

    if (first_call) then
       first_call = .false.
       energy_deriv = 0.0_dp
    else
       ! Time derivative of field energy
       energy_deriv = (new_field_energy - prev_field_energy)/(time - prev_time)
    end if

    if (abs(applied_voltage) > 0.0_dp) then
       J_displ = energy_deriv/applied_voltage
       J_tot = J_displ + JdotE_integral/applied_voltage
       gap_conductance = JdotE_integral/applied_voltage**2
    else
       J_displ = 0.0_dp
       J_tot = 0.0_dp
       gap_conductance = 0.0_dp
    end if

    prev_time = time
    prev_field_energy = new_field_energy
  end subroutine compute_current

  !> Update voltage according to simple RC circuit model. Uses a fixed
  !> geometric C_gap and a conductance G for the conduction current.
  !>
  !> Topology:
  !>     V_C --[C]-- gnd ,  V_C --[R]-- V_gap --( C_gap || G )-- gnd
  !>
  !>   C   dV_C/dt   = -(V_C - V_gap)/R
  !>   C_gap dV_gap/dt = (V_C - V_gap)/R - G*V_gap
  subroutine update_voltage_rc(dt, gap_conductance, V_cap, V_gap)
    real(dp), intent(in)  :: dt, gap_conductance
    real(dp), intent(out) :: V_cap, V_gap
    real(dp)              :: Vc, Vg, a, b, g, det

    Vc = capacitor_voltage
    Vg = applied_voltage

    ! Use a backward Euler approach, which is stable when dt > RC and when dt
    ! > C_gap/gap_conductance
    a = dt/(rc_resistance*rc_capacitance)
    b = dt/(rc_resistance*C_gap)
    g = dt*gap_conductance/C_gap

    ! Solve the 2x2 implicit system:
    !   (1+a) Vc_new -  a     Vg_new = Vc
    !   -b    Vc_new + (1+b+g)Vg_new = Vg
    det = (1.0_dp+a)*(1.0_dp+b+g) - a*b

    V_cap = ((1.0_dp+b+g)*Vc + a*Vg) / det
    V_gap = ((1.0_dp+a)*Vg   + b*Vc) / det
    capacitor_voltage = V_cap
    applied_voltage = V_gap
  end subroutine update_voltage_rc

end module m_solver
