module m_solver
#include "cpp_macros.h"

  use m_solver_lib

  implicit none

contains

  ! Initialize the computational domain
  subroutine initialize_domain(domain_len, coarse_grid_size, box_size, &
       voltage, mem_limit_gb)
    real(dp), intent(in) :: domain_len(fndims)       ! Domain size (m)
    integer, intent(in)  :: coarse_grid_size(fndims) ! Coarse grid size
    integer, intent(in)  :: box_size                 ! Size of grid boxes
    real(dp), intent(in) :: voltage                  ! Applied voltage (V)
    real(dp), intent(in) :: mem_limit_gb             ! Memory limit (GB)
    integer              :: coord_t

    coord_t = af_xyz
    if (fndims == 2) coord_t = af_cyl

    applied_voltage = voltage

    call af_add_cc_variable(tree, "phi", ix=mg%i_phi)
    call af_add_cc_variable(tree, "rhs", ix=mg%i_rhs)
    call af_add_cc_variable(tree, "tmp", ix=mg%i_tmp)
    call af_add_cc_variable(tree, "eps", ix=tree%mg_i_eps)
    call af_add_cc_variable(tree, "sigma", ix=i_sigma)
    call af_add_cc_variable(tree, "dsigma", ix=i_dsigma)
    call af_add_cc_variable(tree, "phi", ix=mg%i_phi)
    call af_add_cc_variable(tree, "electric_fld", ix=i_E_norm)
    call af_add_cc_variable(tree, "time", ix=i_time)
    call af_add_fc_variable(tree, "E_vec", ix=i_E_vec)

    call af_set_cc_methods(tree, tree%mg_i_eps, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_sigma, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_E_norm, af_bc_neumann_zero)
    call af_set_cc_methods(tree, i_time, af_bc_neumann_zero)

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
  end subroutine initialize_domain

  ! Perform uniform initial refinement of the domain
  subroutine use_uniform_grid(uniform_grid_size)
    integer, intent(in) :: uniform_grid_size(fndims)
    integer             :: max_lvl

    ! Refine uniformly
    max_lvl = nint(log(uniform_grid_size(1) / real(tree%n_cell, dp)) / &
         log(2.0_dp)) + 1

    if (any(tree%n_cell * 2**(max_lvl-1) /= uniform_grid_size)) &
         error stop "Incompatible grid size"
    call af_refine_up_to_lvl(tree, max_lvl)
  end subroutine use_uniform_grid

  ! Set parameters for grid refinement
  subroutine set_refinement(refine_field, derefine_field, &
       min_dx, max_dx, electrode_max_dx, derefine_levels)
    real(dp), intent(in) :: refine_field     ! Refine when the field is above this value
    real(dp), intent(in) :: derefine_field   ! Derefine when the field is below this value
    real(dp), intent(in) :: min_dx           ! Minimum allowed grid spacing
    real(dp), intent(in) :: max_dx           ! Maximum allowed grid spacing
    real(dp), intent(in) :: electrode_max_dx ! Maximum grid spacing around electrode
    integer, intent(in)  :: derefine_levels  ! How many levels can be derefined
    integer              :: n_add, n, n_its
    real(dp)             :: residu

    refine_field_threshold   = refine_field
    derefine_field_threshold = derefine_field
    refine_min_dx            = min_dx
    refine_max_dx            = max_dx
    refine_electrode_max_dx  = electrode_max_dx
    ! Finest dx is between min_dx and 2*min_dx; only derefine when dx < derefine_dx
    derefine_dx              = min_dx * 2**derefine_levels

    call solve(0.0_dp, n_its, residu)

    do n = 1, 20
       call adjust_refinement(n_add)
       if (n_add == 0) exit
       call solve(0.0_dp, n_its, residu)
    end do
  end subroutine set_refinement

  ! Update the refinement of the mesh. Changes the local refinement by at most
  ! one level at a time, so multiple calls might be needed.
  subroutine adjust_refinement(n_add)
    integer, intent(out) :: n_add
    type(ref_info_t) :: refine_info

    ! Restrict species, for the ghost cells near refinement boundaries
    call af_restrict_tree(tree, [i_sigma])
    call af_gc_tree(tree, [i_sigma])

    call af_adjust_refinement(tree, refinement_criterion, refine_info, 0)
    n_add = refine_info%n_add
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
  subroutine update_sigma(n_streamers, r0, r1, sigma0, sigma1, radius0, radius1, &
       t, dt, channel_delay, first_step)
    integer, intent(in)  :: n_streamers
    real(dp), intent(in) :: r0(n_streamers, fndims), r1(n_streamers, fndims)
    real(dp), intent(in) :: sigma0(n_streamers), sigma1(n_streamers)
    real(dp), intent(in) :: radius0(n_streamers), radius1(n_streamers)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: channel_delay
    logical, intent(in)  :: first_step
    integer              :: lvl, n, id, IJK, nc, ix
    real(dp)             :: r(fndims), dist_vec(fndims), r_dist, frac
    real(dp)             :: k_eff, dsigma

    nc = tree%n_cell

    if (.not. allocated(k_eff_table)) error stop "Call store_k_eff first"

    !$omp parallel private(lvl, n, id, IJK, r, dist_vec, r_dist, &
    !$omp &frac, ix, k_eff, dsigma)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)

          associate (box => tree%boxes(id))
            do KJI_DO(1, nc)
               box%cc(IJK, i_dsigma) = 0.0_dp

               if (box%cc(IJK, i_lsf) < 0.0_dp) cycle

               r = af_r_cc(box, [IJK])

               do ix = 1, n_streamers
                  call dist_vec_line(r, r0(ix, :), r1(ix, :), &
                       fndims, dist_vec, r_dist, frac)

                  ! Exclude semi-sphere of previous point
                  if (norm2(dist_vec) <= radius1(ix) .and. (first_step .or. &
                       (frac >= 0 .and. norm2(r0(ix, :) - r) > radius0(ix)))) then

                     call get_sigma_profile(r_dist, radius0(ix), radius1(ix), frac, &
                          sigma0(ix), sigma1(ix), dsigma)

                     box%cc(IJK, i_dsigma) = box%cc(IJK, i_dsigma) + dsigma
                     box%cc(IJK, i_time) = t
                  end if
               end do

               ! Update channel conductivity, but only where the channel
               ! has already existed for some time
               if (box%cc(IJK, i_time) < t - channel_delay) then
                  call get_k_eff(box%cc(IJK, i_E_norm), k_eff)

                  ! Limit increase to a factor 2 per time step
                  box%cc(IJK, i_dsigma) = min(1.0_dp, exp(dt * k_eff) - 1.0_dp) * &
                       box%cc(IJK, i_sigma)
               end if

            end do; CLOSE_DO
          end associate
       end do
       !$omp end do
    end do
    !$omp end parallel

    ! Add the change in sigma
    call af_tree_apply(tree, i_sigma, i_dsigma, '+')

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
  subroutine get_k_eff(fld, k_eff)
    real(dp), intent(in)  :: fld
    real(dp), intent(out) :: k_eff
    real(dp)              :: frac, low_frac
    integer               :: low_ix

    frac = (fld - k_eff_table_x_min) * k_eff_table_inv_fac

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
  subroutine store_k_eff(E_min, E_max, n_points, k_eff)
    real(dp), intent(in) :: E_min, E_max
    integer, intent(in)  :: n_points
    real(dp), intent(in) :: k_eff(n_points)

    allocate(k_eff_table(n_points))
    k_eff_table(:) = k_eff
    k_eff_table_n_points = n_points
    k_eff_table_x_min = E_min
    k_eff_table_inv_fac = (n_points - 1)/(E_max - E_min)
  end subroutine store_k_eff

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

    select case (varname)
    case ('sigma')
       i_var = i_sigma
    case ('phi')
       i_var = mg%i_phi
    case ('E_norm')
       i_var = i_E_norm
    case default
       error stop 'Unknown variable'
    end select

    r = r0
    dr = length * direction / (norm2(direction) * (n_steps - 1))

    do n = 1, n_steps
       line(n:n) = af_interp1(tree, r, [i_var], success)
       if (.not. success) return
       z_line(n) = (n-1) * norm2(dr)
       r = r + dr
    end do
  end subroutine get_var_along_line

  ! Compute new potential for a given time step using the current sigma
  subroutine solve(dt, n_iterations, residu)
    real(dp), intent(in)  :: dt
    integer, intent(out)  :: n_iterations
    real(dp), intent(out) :: residu
    integer, parameter    :: max_iterations = 100
    real(dp)              :: prev_residu, max_rhs

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

    do n_iterations = 1, max_iterations
       call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=.true.)
       call af_tree_maxabs_cc(tree, mg%i_tmp, residu)
       if (residu < 1e-5_dp * max_rhs .or. residu > 0.5 * prev_residu) exit
       prev_residu = residu
    end do

    ! Compute new rhs with standard Laplace operator
    call compute_rhs(tree, mg_lpl)

    ! Compute electric field with standard Laplace operator
    call mg_compute_phi_gradient(tree, mg_lpl, i_E_vec, -1.0_dp, i_E_norm)
    call af_gc_tree(tree, [i_E_norm])

  end subroutine solve

  ! Write a silo file
  subroutine write_solution(fname, i_cycle, time)
    character(len=*), intent(in) :: fname
    integer, intent(in)          :: i_cycle
    real(dp), intent(in)         :: time

    ! Ensure valid ghost cells
    call af_restrict_tree(tree, [i_sigma])
    call af_gc_tree(tree, [i_sigma])

    call af_write_silo(tree, trim(fname), i_cycle, time)
  end subroutine write_solution

end module m_solver
