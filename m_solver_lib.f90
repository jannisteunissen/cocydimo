module m_solver_lib
#include "cpp_macros.h"

  use m_af_all

  implicit none

  real(dp), parameter :: eps0 = 8.8541878128e-12_dp ! permitivity of vacuum (SI)
  real(dp), parameter :: pi = acos(-1.0_dp)

  type(af_t) :: tree
  type(mg_t) :: mg
  type(mg_t) :: mg_lpl
  real(dp)   :: phi_bc
  integer    :: i_sigma
  integer    :: i_E_norm
  integer    :: i_E_vec
  integer    :: i_dsigma
  integer    :: i_lsf
  integer    :: i_time
  integer    :: uniform_grid_size(fndims)
  real(dp)   :: min_dr

  ! Electrode parameters
  real(dp) :: rod_r0(fndims), rod_r1(fndims), rod_radius = 0.0_dp

  ! Table with k_eff for updating channel conductivity
  real(dp), allocatable :: k_eff_table(:)
  integer               :: k_eff_table_n_points
  real(dp)              :: k_eff_table_x_min
  real(dp)              :: k_eff_table_inv_fac

  ! To store initial conditions
  real(dp), allocatable :: rhs_input(DTIMES(:)), sigma_input(DTIMES(:))

contains

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc, ix_offset(fndims), ix(fndims)

    nc = box%n_cell
    box%cc(DTIMES(:), mg%i_phi) = 0
    ix_offset = (box%ix - 1) * nc

    do KJI_DO(1, nc)
       ix = ix_offset + [IJK]
       box%cc(IJK, i_sigma)  = sigma_input(DINDEX(ix))
       box%cc(IJK, mg%i_rhs) = rhs_input(DINDEX(ix))
    end do; CLOSE_DO
  end subroutine set_init_cond

  subroutine refinement_criterion(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    integer                 :: nc

    nc = box%n_cell

    if (minval(box%dr) > 1.9_dp * min_dr .and. ( &
         maxval(box%cc(DTIMES(1:nc), i_E_norm)) > 3e6_dp .or. &
         iand(box%tag, mg_lsf_box) > 0)) then
       cell_flags = af_do_ref
    else if (maxval(box%cc(DTIMES(1:nc), i_E_norm)) < 2e6_dp) then
       cell_flags = af_rm_ref
    else
       cell_flags = af_keep_ref
    end if
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
       bc_val = phi_bc
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
