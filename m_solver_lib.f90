module m_solver_lib
#include "cpp_macros.h"

  use m_af_all

  implicit none
  public

  real(dp), parameter :: eps0 = 8.8541878128e-12_dp ! permitivity of vacuum (SI)
  real(dp), parameter :: pi = acos(-1.0_dp)

  type(af_t) :: tree
  type(mg_t) :: mg
  type(mg_t) :: mg_lpl
  real(dp)   :: applied_voltage
  integer    :: i_sigma
  integer    :: i_E_norm
  integer    :: i_E_vec
  integer    :: i_dsigma
  integer    :: i_lsf
  integer    :: i_time

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
         box%cc(DTIMES(1:nc), i_sigma)
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

end module m_solver_lib
