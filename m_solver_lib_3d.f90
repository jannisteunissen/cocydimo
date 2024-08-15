module m_solver_lib_3d
  use m_af_all

  implicit none

  real(dp), parameter :: eps0 = 8.8541878176d-12 ! permitivity of vacuum (SI)
  real(dp), parameter :: pi = acos(-1.0_dp)

  type(af_t) :: tree
  type(mg_t) :: mg
  type(mg_t) :: mg_lpl
  real(dp)   :: phi_bc
  integer    :: i_sigma
  integer    :: i_dsigma
  integer    :: i_lsf
  integer    :: i_E_norm
  integer    :: i_E_vec
  integer    :: uniform_grid_size(3)

  ! Electrode parameters
  real(dp) :: rod_r0(3), rod_r1(3), rod_radius = 0.0_dp

  real(dp), allocatable :: rhs_input(:, :, :), sigma_input(:, :, :)

contains

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                    :: i, j, k, nc, ix_offset(3), k1, k2, k3

    nc = box%n_cell
    box%cc(:, :, :, mg%i_phi) = 0
    ix_offset = (box%ix - 1) * nc

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             k1 = ix_offset(1) + i
             k2 = ix_offset(2) + j
             k3 = ix_offset(3) + k

             box%cc(i, j, k, i_sigma)  = sigma_input(k1, k2, k3)
             box%cc(i, j, k, mg%i_rhs) = rhs_input(k1, k2, k3)
          end do
       end do
    end do
  end subroutine set_init_cond

  subroutine set_epsilon_from_sigma(box, dt_vec)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt_vec(:)
    integer                    :: nc

    nc = box%n_cell
    box%cc(1:nc, 1:nc, 1:nc, tree%mg_i_eps) = 1 + (dt_vec(1)/eps0) * &
         box%cc(1:nc, 1:nc, 1:nc, i_sigma)
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
    real(dp), intent(in)    :: coords(3, box%n_cell**2)
    real(dp), intent(out)   :: bc_val(box%n_cell**2)
    integer, intent(out)    :: bc_type

    if (nb == af_neighb_lowz) then
       bc_type = af_bc_dirichlet
       bc_val = 0.0_dp
    else if (nb == af_neighb_highz) then
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
    real(dp), intent(out) :: dist_line !< Distance from line if it were infinite
    real(dp), intent(out) :: frac !< Fraction [0,1] along line
    real(dp)              :: line_len2

    ! Distance to infinite line
    line_len2 = sum((r1 - r0)**2)
    frac = sum((r - r0) * (r1 - r0))/line_len2
    dist_vec = r - r0 - frac * (r1 - r0)
    dist_line = norm2(dist_vec)

    ! Adjust for finite line
    if (frac <= 0.0_dp) then
       frac = 0.0_dp
       dist_vec = r - r0
    else if (frac >= 1.0_dp) then
       frac = 1.0_dp
       dist_vec = r - r1
    end if
  end subroutine dist_vec_line

  ! Level set function for rod electrode
  real(dp) function rod_lsf(r)
    real(dp), intent(in) :: r(3)
    rod_lsf = get_dist_line(r, rod_r0, rod_r1, 3) - rod_radius
  end function rod_lsf

  subroutine set_lsf_box(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: i, j, k, nc
    real(dp)                   :: rr(3)

    nc = box%n_cell
    do k = 0, nc+1
       do j = 0, nc+1
          do i = 0, nc+1
             rr = af_r_cc(box, [i, j, k])
             box%cc(i, j, k, iv) = rod_lsf(rr)
          end do
       end do
    end do
  end subroutine set_lsf_box

  function get_dist_line(r, r0, r1, n_dim) result(dist)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp)             :: dist, dist_vec(n_dim), dist_line, frac
    call dist_vec_line(r, r0, r1, n_dim, dist_vec, dist_line, frac)
    dist = norm2(dist_vec)
  end function get_dist_line

end module m_solver_lib_3d
