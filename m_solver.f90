module m_solver
#include "cpp_macros.h"

  use m_solver_lib

  implicit none

contains

  ! Set the domain size, the finest grid spacing (if the mesh was uniform) and
  ! other domain properties
  subroutine initialize_domain(domain_len, grid_size, box_size, applied_voltage)
    real(dp), intent(in) :: domain_len(fndims)
    integer, intent(in)  :: grid_size(fndims)
    integer, intent(in)  :: box_size
    real(dp), intent(in) :: applied_voltage
    integer              :: n, n_add, max_lvl, coord_t, n_its
    real(dp)             :: residu

    coord_t = af_xyz
    if (fndims == 2) coord_t = af_cyl

    phi_bc = applied_voltage
    uniform_grid_size = grid_size
    min_dr = minval(domain_len / uniform_grid_size)

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

    call af_init(tree, box_size, domain_len, &
         [box_size, box_size, box_size], coord=coord_t, &
         mem_limit_gb=2.0_dp)

    mg%sides_bc => sides_bc ! Method for boundary conditions

    ! Create a copy of the operator but without the variable coefficient
    mg_lpl = mg
    mg_lpl%operator_mask = mg_normal_box + mg_lsf_box

    if (fndims == 2) then
       ! Refine uniformly so that we can use simulation data initially
       max_lvl = nint(log(grid_size(1) / real(box_size, dp))/log(2.0_dp)) + 1
       if (any(box_size * 2**(max_lvl-1) /= grid_size)) &
            error stop "Incompatible grid size"
       call af_refine_up_to_lvl(tree, max_lvl)
    else
       ! Use grid refinement, which requires the electric field
       call solve(0.0_dp, n_its, residu)

       do n = 1, 10
          call adjust_refinement(n_add)
          if (n_add == 0) exit
          print *, "Initial refinement, step", n
          call solve(0.0_dp, n_its, residu)
       end do
    end if

  end subroutine initialize_domain

  ! Set the electric field thresholds for grid refinement
  subroutine set_refinement_thresholds(refine_field, derefine_field)
    real(dp), intent(in) :: refine_field ! Refine when the field is above this value
    real(dp), intent(in) :: derefine_field ! Derefine when the field is below this value

    refine_threshold = refine_field
    derefine_threshold = derefine_field
  end subroutine set_refinement_thresholds

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

  ! Get the maximum electric field and its location
  subroutine get_max_field_location(E_max_vec, r)
    real(dp), intent(out) :: E_max_vec(fndims)
    real(dp), intent(out) :: r(fndims)
    real(dp)              :: E_max_norm
    type(af_loc_t)        :: loc

    call af_tree_max_cc(tree, i_E_norm, E_max_norm, loc)
    r = af_r_loc(tree, loc)
    call get_field_vector_at(r, E_max_vec)
  end subroutine get_max_field_location

  ! Get the electric field vector at a location
  subroutine get_field_vector_at(r, E_vec)
    real(dp), intent(in)  :: r(fndims)
    real(dp), intent(out) :: E_vec(fndims)
    logical               :: success

    E_vec = af_interp1_fc(tree, r, i_E_vec, success)
    if (.not. success) error stop "Interpolation error"
  end subroutine get_field_vector_at

  ! Trace the value of a variable along a line
  subroutine get_var_along_line(varname, r0, direction, length, n_steps, &
       z_line, line)
    character(len=*), intent(in) :: varname
    real(dp), intent(in)         :: r0(fndims), direction(fndims), length
    integer, intent(in)          :: n_steps
    real(dp), intent(out)        :: z_line(n_steps)
    real(dp), intent(out)        :: line(n_steps)
    real(dp)                     :: r(fndims), dr(fndims)
    integer                      :: n, i_var
    logical                      :: success

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
       z_line(n) = (n-1) * norm2(dr)
       if (.not. success) error stop "Interpolation error"
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
