module spk_finites
    !> Contains calculation of finite differences for selected quantities
    USE, INTRINSIC :: IEEE_ARITHMETIC
    use spk_input_lattice, only: N_atom, latt_vecs
    use spk_precision_const, only: xp
    use spk_rotations, only: rotate_spinconfiguration, rotate_vector_between_tangent_spaces
    use spk_tangentspace, only: vec_to_tangentspace
    use spk_first_derivatives, only: gradient, gradient_subsystem
    use spk_tangentspace, only: create_tangentspace_basis
    use spk_vector_lina, only: norm
    use spk_bilinear_shells, only: calc_shell_distances_a1a2, calc_shell_Nmembers_a1a2
    use spk_bilinear_neighbors, only: calc_a1a2_indexshift_allbds, calc_a1a2_neightable_allbds
    use spk_atomic_lattice, only: atom_positions
    use spk_modeconstruction, only: translationmode, rotationmode
    use spk_casting, only: dp_e_to_str
    use spk_logging, only: log
    use spk_abbreviations_const, only: C_WARNING
    implicit none

    interface fd_hessvec_forward
        module procedure fd_hessvec_forward_richardson_order1
        module procedure fd_hessvec_forward_richardson_order1_subsys
        module procedure fd_hessvec_forward_orderN
        module procedure fd_hessvec_forward_orderN_subsys
    end interface fd_hessvec_forward

    public :: findiff_shapeoperator, findiff_basis, findiff_jacobian_xy, &
            & findiff_jacobian_preconditioner, fd_hessvec_forward
    private
    contains
        subroutine findiff_shapeoperator(spin,grad,spin_next,grad_next,step,diff_shape)
            !> Calculates the finite difference of two adjacent shape operators Gamma:
            !> Diff_Gamma = 1/step * (Gamma_next - Gamma)
            !> With the shape operator being: Gamma_i = dot (spin_i, grad_i)
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom), spin_next(3*N_atom)    !< Two "adjacent" spin configuration
            real(kind=xp), intent(in) :: grad(3*N_atom), grad_next(3*N_atom)    !< Gradients of "adjacent" spin config.
            real(kind=xp), intent(in) :: step                                   !< distance
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: diff_shape(N_atom)                    !< finite difference of shape operator
            !================================Local Variable=============================================================
            integer :: i_atom                                                   !< iterator for atom
            real(kind=xp) :: shape_i, shape_next_i                              !< shape entries for i-th atom

            diff_shape = 0.0d0
            do i_atom=1,N_atom
                shape_i = (spin(3*i_atom-2) * grad(3*i_atom-2) + spin(3*i_atom-1) * grad(3*i_atom-1) &
                        & + spin(3*i_atom) * grad(3*i_atom))
                shape_next_i = (spin_next(3*i_atom-2) * grad_next(3*i_atom-2) &
                        & + spin_next(3*i_atom-1) * grad_next(3*i_atom-1) + spin_next(3*i_atom) * grad_next(3*i_atom))
                diff_shape(i_atom) = (1.0d0 / step) * (shape_next_i-shape_i)
            end do
        end subroutine findiff_shapeoperator

        subroutine findiff_basis(basis, basis_next, step, diff_basis)
            !> Calculates the finite difference for a 3x2xN_atom tangent space basis.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: basis(3,2,N_atom), basis_next(3,2,N_atom)  !< basis of two adjacent configs.
            real(kind=xp), intent(in) :: step                                       !< distance
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: diff_basis(3,2,N_atom)                    !< finite difference for the basis
            !================================Local Variable=============================================================
            integer :: i_atom                                                       !< iterator for atom
            diff_basis = 0.0d0
            do i_atom=1,N_atom
                diff_basis(:,:,i_atom) = (1.0d0 / step) * (basis_next(:,:,i_atom)-basis(:,:,i_atom))
            end do
        end subroutine findiff_basis

        subroutine findiff_jacobian_preconditioner(n_shells, n_neigh, neigh_table, neigh_mask)
            !> Preconditioner routine for calling the jacobian finite difference approximation
            !================================Input Variable=============================================================
            integer, intent(in) :: n_shells                         !< How many neighbor shells to be considered for the finite
                                                                    !< difference approximation
            !================================Output Variable============================================================
            integer, intent(out) :: n_neigh                         !< number of neighbors for all shells
            integer, allocatable, intent(out) :: neigh_table(:,:)   !< neighbor table
            integer, allocatable, intent(out) :: neigh_mask(:,:)    !< accounting for periodic boundaries
            !================================Local Variable=============================================================
            real(kind=xp) :: shell_distances(n_shells)              !< shell distances in the a1 a2 plane
            integer :: shell_members(n_shells)                      !< number of members for these shells
            integer, allocatable :: index_shifts(:,:)               !< index shifts for the neighbors

            call calc_shell_distances_a1a2(n_shells, shell_distances)
            call calc_shell_Nmembers_a1a2(n_shells, shell_distances, shell_members)
            call calc_a1a2_indexshift_allbds(n_shells, shell_distances, shell_members, index_shifts)
            n_neigh = size(index_shifts, 1)
            allocate(neigh_table(N_atom,n_neigh))
            allocate(neigh_mask(N_atom, n_neigh))
            call calc_a1a2_neightable_allbds(n_neigh, index_shifts, neigh_table, neigh_mask)
        end subroutine findiff_jacobian_preconditioner

        subroutine findiff_jacobian_xy(spin, n_neigh, neigh_table, neigh_mask, J)
            !> Calculates the jacobian matrix of the spin configuration J(3,N) which is the 3 x N matrix which gives the
            !> derivative of each spin with respect to the positional coordinates x,y,z. So if we have a spin i the
            !> Jacobian would be            d_x s_x  d_y s_x d_z s_x
            !> J(3,3*i-2:3*i) = J_i(3,3) =  d_x s_y  d_y s_y d_z s_y.
            !>                              d_x s_z  d_y s_z d_z s_z
            !> These derivatives will be calculated by a finite difference sheme based on spins at discrete points. All
            !> derivatives in z direction (d_z) are set to zero (J_i(3,:)=0.0). Due to the general restrictions of
            !> Spinaker the a1 and a2 lattice vectors are always in the xy-plane (without loss of generality for the
            !> lattice model). Also only one atom per atomic species is allowed per unit cell. Therefore all interaction
            !> class groups like (j,j) for have the same index shift. However, a translation mode might exist although
            !> a certain (j,j) interaction class group is not defined. Therefore, we can't use the neighbor table as
            !> only neighbors for existing interaction class groups are included. To calculate the partial derivatives in x,y
            !> direction the following procedure is applied.
            !> For each spin a given number of
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom)             !< Current spin configuration
            integer, intent(in) :: n_neigh                          !< Number of neighbors for all shells
            integer, allocatable, intent(in) :: neigh_table(:,:)    !< Neighbor table
            integer, allocatable, intent(in) :: neigh_mask(:,:)     !< Accounting for periodic boundaries
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: J(3,3*N_atom)             !< Jacobian Matrix for each spin
            !================================Local Variable=============================================================
            integer :: i_atom                                       !< Iterator for atoms
            integer :: j_neigh                                      !< Iterator for neighbors
            integer :: n_neigh_counter                              !< Neighbor counter
                                                                    !< (can be reduced by not per. bound.)
            integer :: j_idx_neigh                                  !< Index of the neighbor atom in the atom list
            real(kind=xp) :: current_pos(3)                         !< Position of current atom
            real(kind=xp) :: neigh_pos(3)                           !< Position of the current neighbor
            real(kind=xp) :: dist                                   !< Distance between atom and neighbor
            real(kind=xp) :: dsx_dx, dsy_dx, dsz_dx, dsx_dy, &      !< Partial derivatives
                    & dsy_dy, dsz_dy
            real(kind=xp) :: J_i(3,3)                               !< Jacobian for the ith spin
            J=0.0d0
            do i_atom=1,N_atom
                J_i = 0.0d0
                current_pos = atom_positions(3*i_atom-2:3*i_atom)
                dsx_dx = 0.0
                dsx_dy = 0.0
                dsy_dx = 0.0
                dsy_dy = 0.0
                dsz_dx = 0.0
                dsz_dy = 0.0
                n_neigh_counter = n_neigh
                do j_neigh=1,n_neigh
                    if (neigh_mask(i_atom,j_neigh)==0) then
                        n_neigh_counter = n_neigh_counter - 1
                        cycle
                    end if
                    j_idx_neigh = neigh_table(i_atom, j_neigh)
                    neigh_pos = atom_positions(3*j_idx_neigh-2:3*j_idx_neigh)
                    dist = norm(neigh_pos-current_pos)
                    dsx_dx = dsx_dx + (spin(3*j_idx_neigh-2) - spin(3*i_atom-2)) / dist**2 * (neigh_pos(1)-current_pos(1))
                    dsy_dx = dsy_dx + (spin(3*j_idx_neigh-1) - spin(3*i_atom-1)) / dist**2 * (neigh_pos(1)-current_pos(1))
                    dsz_dx = dsz_dx + (spin(3*j_idx_neigh) - spin(3*i_atom)) / dist**2 * (neigh_pos(1)-current_pos(1))
                    dsx_dy = dsx_dy + (spin(3*j_idx_neigh-2) - spin(3*i_atom-2)) / dist**2 * (neigh_pos(2)-current_pos(2))
                    dsy_dy = dsy_dy + (spin(3*j_idx_neigh-1) - spin(3*i_atom-1)) / dist**2 * (neigh_pos(2)-current_pos(2))
                    dsz_dy = dsz_dy + (spin(3*j_idx_neigh) - spin(3*i_atom)) / dist**2 * (neigh_pos(2)-current_pos(2))
                end do
                ! Derivatives in z-direction are considered 0.
                J_i(1,1) = dsx_dx / n_neigh_counter
                J_i(1,2) = dsx_dy / n_neigh_counter
                J_i(2,1) = dsy_dx / n_neigh_counter
                J_i(2,2) = dsy_dy / n_neigh_counter
                J_i(3,1) = dsz_dx / n_neigh_counter
                J_i(3,2) = dsz_dy / n_neigh_counter
                J(:,3*i_atom-2:3*i_atom) = J_i
            end do
        end subroutine findiff_jacobian_xy

        subroutine fd_hessvec_forward_richardson_order1(spin, vector, error_allowed,step_guess, step_used, diff_hessianvec,&
                & n_error_rescale)
            !> Computes product of the hessian with some vector based on finite differences of the gradients. Pearlmutter
            !> trick. For finite difference method richardson-extrapolation is used. This allows for providing a desired
            !> accuracy (by "error_allowed" input argument). The finite difference scheme used is of order 1 (forward).
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom)                     !< Spin configuration
            real(kind=xp), intent(in) :: step_guess                         !< Finite difference step (will be rescaled
                                                                            !< If error is above error tolerance)
            real(kind=xp), intent(in) :: vector(3*N_atom)                   !< Input vector
            real(kind=xp), intent(in) :: error_allowed                      !< Allowed error for finite difference
                                                                            !< Approximation.
            integer, intent(in) :: n_error_rescale                          !< Number of attempts to rescale the step
                                                                            !< In order to fullfill the allowed error
                                                                            !< Criterium
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: diff_hessianvec(3*N_atom)         !< Output approximation of hessian * vec
            real(kind=xp), intent(out) :: step_used                         !< Used finite difference step size.
            !================================Local Variable=============================================================
            real(kind=xp) :: step_half                                      !< Half Step. For Richardson Extrapolation
            real(kind=xp) :: fd_step(3*N_atom), fd_step_half(3*N_atom)      !< Finite Difference Approximation for the
                                                                            !< the step and the half step.
            real(kind=xp) :: current_step                                   !< Current used step size
            real(kind=xp) :: spin_rotated(3*N_atom)                         !< Spin configuration rotated along step
            real(kind=xp) :: current_gradient(3*N_atom)                     !< Gradient
            real(kind=xp) :: gradient_tspace_notrotated(3*N_atom)           !< Gradient of not rotated Spin Conifg.
            real(kind=xp) :: current_gradient_tspace(3*N_atom)              !< Gradient in tangent space
            real(kind=xp) :: current_gradient_tspace_initial(3*N_atom)      !< Gradient in tangent space of NOT rotated
                                                                            !< Spin Texture
            real(kind=xp) :: back_displacement(3*N_atom)                    !< Negative Step to rotate back to initial
                                                                            !< Tangent frame
            real(kind=xp) :: fd_richardson(3*N_atom)                        !< Finite difference of directional der.
                                                                            !< By Richardson-Extrapolation
            real(kind=xp) :: fd_richardson_previous(3*N_atom)               !< Approximation of previous iteration
            real(kind=xp) :: error_previous                                 !< Error of previous iteration
            real(kind=xp) :: error_richardson                               !< Error of the approximation
            real(kind=xp) :: step_previous                                  !< Step size of previous iteration
            real(kind=xp) :: error_rescale_factor                           !< Prefactor to do careful error rescaling
            integer :: i_rescale_attempt                                    !< Iterator for step size rescale attempts

            diff_hessianvec = 0.0d0
            error_richardson = 1.0d0
            step_used = step_guess
            error_rescale_factor = 0.9d0
            call gradient(spin, current_gradient)
            call vec_to_tangentspace(spin, current_gradient, gradient_tspace_notrotated)
            outer: do i_rescale_attempt=1,n_error_rescale
                step_half = step_used / 2.0
                ! Calculate the finite difference approximation of the step-displacement
                spin_rotated = spin
                current_step = step_used
                call rotate_spinconfiguration(spin_rotated, vector, current_step)
                call gradient(spin_rotated, current_gradient)
                call vec_to_tangentspace(spin_rotated, current_gradient, current_gradient_tspace)
                back_displacement = -1.0d0 * current_step * vector
                call rotate_vector_between_tangent_spaces(spin_rotated, back_displacement, current_gradient_tspace, &
                        & current_gradient_tspace_initial)
                fd_step = (current_gradient_tspace_initial - gradient_tspace_notrotated) / current_step
                ! Calculate the finite difference approximation of the half step-displacement
                spin_rotated = spin
                current_step = step_half
                call rotate_spinconfiguration(spin_rotated, vector, current_step)
                call gradient(spin_rotated, current_gradient)
                call vec_to_tangentspace(spin_rotated, current_gradient, current_gradient_tspace)
                back_displacement = -1.0d0 * current_step * vector
                call rotate_vector_between_tangent_spaces(spin_rotated, back_displacement, current_gradient_tspace, &
                        & current_gradient_tspace_initial)
                fd_step_half = (current_gradient_tspace_initial - gradient_tspace_notrotated) / current_step
                ! Calculate the Richardson Extrapolation
                fd_richardson_previous = fd_richardson
                fd_richardson = 2.0 * fd_step_half - fd_step
                if (step_used<=1.0d-10) then
                    call log("Step below allowed threshold: "//trim(adjustl(dp_e_to_str(current_step))), &
                            & prefix = C_WARNING, idt_level=2)
                    call log("Exit finite difference with current approx.",prefix=C_WARNING, idt_level=2)
                    exit outer
                end if
                error_previous = error_richardson
                error_richardson = norm(fd_richardson-fd_step_half)
                if (IEEE_IS_NAN(error_richardson)) then
                    call log("Richardson error is NaN. Will exit",prefix=C_WARNING,idt_level=3)
                    stop
                    spin_rotated = spin
                    current_step = step_guess
                    call rotate_spinconfiguration(spin_rotated, vector, current_step)
                    call gradient(spin_rotated, current_gradient)
                    call vec_to_tangentspace(spin_rotated, current_gradient, current_gradient_tspace)
                    back_displacement = -1.0d0 * current_step * vector
                    call rotate_vector_between_tangent_spaces(spin_rotated, back_displacement, current_gradient_tspace, &
                            & current_gradient_tspace_initial)
                    ! Just assign the simple FD approximation without extrapolation in this case
                    fd_richardson = (current_gradient_tspace_initial - gradient_tspace_notrotated) / current_step
                    exit outer
                end if
                if (error_richardson<=error_allowed) then
                    exit outer
                end if
                if (error_richardson>=error_previous) then
                    fd_richardson = fd_richardson_previous
                    error_richardson = error_previous
                    step_previous = step_used ! Is that correct? Shouldnt it be step_used = step_previous?
                    call log("Minimum error which could be reached was: "//trim(adjustl(dp_e_to_str(error_richardson))), &
                            & prefix = C_WARNING, idt_level=2)
                    exit outer
                end if
                if (i_rescale_attempt == n_error_rescale) then
                    call log("Minimum error which could be reached was: "//trim(adjustl(dp_e_to_str(error_richardson))), &
                            & prefix = C_WARNING, idt_level=2)
                    exit outer
                else
                    step_previous = step_used
                    step_used = error_rescale_factor * step_used * sqrt(error_allowed / error_richardson)
                end if
            end do outer
            diff_hessianvec = fd_richardson

        end subroutine fd_hessvec_forward_richardson_order1

        subroutine fd_hessvec_forward_richardson_order1_subsys(spin, vector,subsys_indices, error_allowed,step_guess, &
                & step_used, diff_hessianvec,n_error_rescale)
            !> Computes product of the hessian with some vector based on finite differences of the gradients. Pearlmutter
            !> trick. For finite difference method richardson-extrapolation is used. This allows for providing a desired
            !> accuracy (by "error_allowed" input argument). The finite difference scheme used is of order 1 (forward).
            !> All quantities are considered for subsystem
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(:)                            !< Spin configuration
            real(kind=xp), intent(in) :: step_guess                         !< Finite difference step (will be rescaled
                                                                            !< If error is above error tolerance)
            integer, intent(in) :: subsys_indices(:)                        !< Indices of sub-system
            real(kind=xp), intent(in) :: vector(:)                          !< Input vector (3N)
            real(kind=xp), intent(in) :: error_allowed                      !< Allowed error for finite difference
                                                                            !< Approximation.
            integer, intent(in) :: n_error_rescale                          !< Number of attempts to rescale the step
                                                                            !< In order to fullfill the allowed error
                                                                            !< Criterium
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: diff_hessianvec(:)                !< Output approximation of hessian * vec (3N)
            real(kind=xp), intent(out) :: step_used                         !< Used finite difference step size.
            !================================Local Variable=============================================================
            integer :: N                                                    !< Number of atoms in sub-system
            real(kind=xp) :: step_half                                      !< Half Step. For Richardson Extrapolation
            real(kind=xp), allocatable :: fd_step(:), fd_step_half(:)       !< Finite Difference Approximation for the
                                                                            !< the step and the half step.
            real(kind=xp) :: current_step                                   !< Current used step size
            real(kind=xp), allocatable :: spin_rotated(:)                   !< Spin configuration rotated along step
            real(kind=xp), allocatable :: current_gradient(:)               !< Gradient
            real(kind=xp), allocatable :: gradient_tspace_notrotated(:)     !< Gradient of not rotated Spin Conifg.
            real(kind=xp), allocatable :: current_gradient_tspace(:)        !< Gradient in tangent space
            real(kind=xp), allocatable :: current_gradient_tspace_initial(:)!< Gradient in tangent space of NOT rotated
                                                                            !< Spin Texture
            real(kind=xp), allocatable :: back_displacement(:)              !< Negative Step to rotate back to initial
                                                                            !< Tangent frame
            real(kind=xp), allocatable :: fd_richardson(:)                  !< Finite difference of directional der.
                                                                            !< By Richardson-Extrapolation
            real(kind=xp), allocatable :: fd_richardson_previous(:)         !< Approximation of previous iteration
            real(kind=xp) :: error_previous                                 !< Error of previous iteration
            real(kind=xp) :: error_richardson                               !< Error of the approximation
            real(kind=xp) :: step_previous                                  !< Step size of previous iteration
            real(kind=xp) :: error_rescale_factor                           !< Prefactor to do careful error rescaling
            integer :: i_rescale_attempt                                    !< Iterator for step size rescale attempts

            N = size(subsys_indices)
            allocate(fd_step(3*N),fd_step_half(3*N),current_gradient(3*N),gradient_tspace_notrotated(3*N), &
                    & current_gradient_tspace(3*N),current_gradient_tspace_initial(3*N),back_displacement(3*N), &
                    & fd_richardson(3*N),fd_richardson_previous(3*N))
            allocate(spin_rotated(3*N_atom))

            diff_hessianvec = 0.0d0
            error_richardson = 1.0d0
            step_used = step_guess
            error_rescale_factor = 0.9d0
            call gradient_subsystem(spin,subsys_indices,current_gradient)
            call vec_to_tangentspace(spin,subsys_indices,current_gradient, gradient_tspace_notrotated)
            outer: do i_rescale_attempt=1,n_error_rescale
                step_half = step_used / 2.0
                ! Calculate the finite difference approximation of the step-displacement
                spin_rotated = spin
                current_step = step_used
                if (current_step<=1.0d-10) then
                    call log("Step below allowed threshold: "//trim(adjustl(dp_e_to_str(current_step))), &
                            & prefix = C_WARNING, idt_level=2)
                    call log("Exit finite difference with current approx.",prefix=C_WARNING, idt_level=2)
                    exit
                end if
                call rotate_spinconfiguration(spin_rotated,subsys_indices, vector, current_step)
                call gradient_subsystem(spin_rotated,subsys_indices,current_gradient)
                call vec_to_tangentspace(spin_rotated,subsys_indices,current_gradient, current_gradient_tspace)
                back_displacement = -1.0d0 * current_step * vector
                call rotate_vector_between_tangent_spaces(spin_rotated,subsys_indices,back_displacement, &
                        & current_gradient_tspace, current_gradient_tspace_initial)
                fd_step = (current_gradient_tspace_initial - gradient_tspace_notrotated) / current_step
                ! Calculate the finite difference approximation of the half step-displacement
                spin_rotated = spin
                current_step = step_half
                call rotate_spinconfiguration(spin_rotated,subsys_indices, vector, current_step)
                call gradient_subsystem(spin_rotated,subsys_indices,current_gradient)
                call vec_to_tangentspace(spin_rotated,subsys_indices, current_gradient, current_gradient_tspace)
                back_displacement = -1.0d0 * current_step * vector
                call rotate_vector_between_tangent_spaces(spin_rotated,subsys_indices, back_displacement, &
                        & current_gradient_tspace, current_gradient_tspace_initial)
                fd_step_half = (current_gradient_tspace_initial - gradient_tspace_notrotated) / current_step
                ! Calculate the Richardson Extrapolation
                fd_richardson_previous = fd_richardson
                fd_richardson = 2.0 * fd_step_half - fd_step
                error_previous = error_richardson
                error_richardson = norm(fd_richardson-fd_step_half)
                if (error_richardson<=error_allowed) then
                    exit outer
                end if
                if (error_richardson>=error_previous) then
                    fd_richardson = fd_richardson_previous
                    error_richardson = error_previous
                    call log("Minimum error which could be reached was: "//trim(adjustl(dp_e_to_str(error_richardson))), &
                            & prefix = C_WARNING, idt_level=2)
                    exit outer
                end if
                if (i_rescale_attempt == n_error_rescale) then
                    call log("Minimum error which could be reached was: "//trim(adjustl(dp_e_to_str(error_richardson))), &
                            & prefix = C_WARNING, idt_level=2)
                    exit outer
                else
                    step_previous = step_used
                    step_used = error_rescale_factor * step_used * sqrt(error_allowed / error_richardson)
                end if
            end do outer
            diff_hessianvec = fd_richardson
        end subroutine fd_hessvec_forward_richardson_order1_subsys

        subroutine fd_hessvec_forward_orderN(spin, step, order, vector, diff_hessianvec)
            !> Computes product of the hessian with some vector based on finite differences of the gradients. Pearlmutter
            !> trick. The finite difference scheme is forward with a variable number of stencil points.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom)                     !< Spin configuration
            real(kind=xp), intent(in) :: step                               !< Finite difference step
            integer, intent(in) :: order                                    !< Order of the scheme (stencil points)
            real(kind=xp), intent(in) :: vector(3*N_atom)                   !< Input vector
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: diff_hessianvec(3*N_atom)         !< Output approximation of hessian * vec
            !================================Local Variable=============================================================
            real(kind=xp), allocatable :: coefficients(:)                   !< Coefficients
            real(kind=xp) :: spin_rotated(3*N_atom)
            real(kind=xp) :: current_step
            real(kind=xp) :: current_gradient(3*N_atom)
            real(kind=xp) :: current_gradient_tspace(3*N_atom)
            real(kind=xp) :: current_gradient_tspace_initial(3*N_atom)
            real(kind=xp) :: back_displacement(3*N_atom)
            real(kind=xp) :: vector_tspace(3*N_atom)
            integer :: i_order                                              !< Current order of the sheme

            diff_hessianvec = 0.0d0
            vector_tspace = vector
            allocate(coefficients(order+1))
            !Finite difference coefficients for forward sheme. See: wikipedia: Finite difference coefficient
            select case(order)
                case(1)
                    coefficients=(/-1.0d0, 1.0d0/)
                case(2)
                    coefficients=(/-3.0d0/2.0d0, 2.0d0, -1.0d0/2.0d0/)
                case(3)
                    coefficients=(/-11.0d0/6.0d0, 3.0d0, -3.0d0/2.0d0, 1.0d0/3.0d0/)
                case(4)
                    coefficients=(/-25.0d0/12.0d0, 4.0d0, -3.0d0, 4.0d0/3.0d0, -1.0d0/4.0d0/)
                case(5)
                    coefficients=(/-137.0d0/60.0d0, 5.0d0, -5.0d0, 10.0d0/3.0d0, -5.0d0/4.0d0, 1.0d0/5.0d0/)
                case(6)
                    coefficients = (/-49.0d0/20.0d0, 6.0d0, -15.0d0/2.d0, 20.0d0 /3.0d0, -15.0d0/4.0d0, 6.0d0/5.0d0, &
                            &-1.0d0/6.0d0/)
            end select

            call gradient(spin, current_gradient)
            call vec_to_tangentspace(spin, current_gradient, current_gradient_tspace)
            do i_order=1,order + 1
                spin_rotated = spin
                current_step = step * (i_order - 1)
                call rotate_spinconfiguration(spin_rotated, vector_tspace, current_step)
                call gradient(spin_rotated, current_gradient)
                call vec_to_tangentspace(spin_rotated, current_gradient, current_gradient_tspace)
                back_displacement = -1.0d0 * current_step * vector_tspace
                call rotate_vector_between_tangent_spaces(spin_rotated, back_displacement, current_gradient_tspace, &
                        & current_gradient_tspace_initial)
                diff_hessianvec = diff_hessianvec + current_gradient_tspace_initial * coefficients(i_order)
            end do
            diff_hessianvec = diff_hessianvec / (step*1)

        end subroutine fd_hessvec_forward_orderN

        subroutine fd_hessvec_forward_orderN_subsys(spin, step,subsys_indices, order, vector, diff_hessianvec)
            !> Computes product of the hessian with some vector based on finite differences of the gradients. Pearlmutter
            !> trick. The finite difference scheme is forward with a variable number of stencil points.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(:)                            !< Spin configuration
            real(kind=xp), intent(in) :: step                               !< Finite difference step
            integer, intent(in) :: order                                    !< Order of the scheme (stencil points)
            real(kind=xp), intent(in) :: vector(:)                          !< Input vector (3N)
            integer, intent(in) :: subsys_indices(:)                        !< Indices of atoms in sub-system
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: diff_hessianvec(:)                !< Output approximation of hessian * vec (3N)
            !================================Local Variable=============================================================
            integer :: N                                                    !< Number of atoms in subsystem
            real(kind=xp), allocatable :: coefficients(:)                   !< Coefficients
            real(kind=xp), allocatable :: spin_rotated(:)
            real(kind=xp) :: current_step
            real(kind=xp), allocatable :: current_gradient(:)
            real(kind=xp), allocatable :: current_gradient_tspace(:)
            real(kind=xp), allocatable :: current_gradient_tspace_initial(:)
            real(kind=xp), allocatable :: back_displacement(:)
            real(kind=xp), allocatable :: vector_tspace(:)
            integer :: i_order                                              !< Current order of the sheme

            diff_hessianvec = 0.0d0
            N = size(subsys_indices)
            allocate(coefficients(order+1))
            allocate(spin_rotated(3*N_atom))
            allocate(current_gradient(3*N))
            allocate(current_gradient_tspace(3*N))
            allocate(current_gradient_tspace_initial(3*N))
            allocate(back_displacement(3*N))
            allocate(vector_tspace(3*N))
            vector_tspace = vector

            !Finite difference coefficients for forward sheme. See: wikipedia: Finite difference coefficient
            select case(order)
                case(1)
                    coefficients=(/-1.0d0, 1.0d0/)
                case(2)
                    coefficients=(/-3.0d0/2.0d0, 2.0d0, -1.0d0/2.0d0/)
                case(3)
                    coefficients=(/-11.0d0/6.0d0, 3.0d0, -3.0d0/2.0d0, 1.0d0/3.0d0/)
                case(4)
                    coefficients=(/-25.0d0/12.0d0, 4.0d0, -3.0d0, 4.0d0/3.0d0, -1.0d0/4.0d0/)
                case(5)
                    coefficients=(/-137.0d0/60.0d0, 5.0d0, -5.0d0, 10.0d0/3.0d0, -5.0d0/4.0d0, 1.0d0/5.0d0/)
                case(6)
                    coefficients = (/-49.0d0/20.0d0, 6.0d0, -15.0d0/2.d0, 20.0d0 /3.0d0, -15.0d0/4.0d0, 6.0d0/5.0d0, &
                            &-1.0d0/6.0d0/)
            end select

            call gradient_subsystem(spin, subsys_indices, current_gradient)
            call vec_to_tangentspace(spin,subsys_indices, current_gradient, current_gradient_tspace)
            do i_order=1,order + 1
                spin_rotated = spin
                current_step = step * (i_order - 1)
                call rotate_spinconfiguration(spin_rotated,subsys_indices, vector_tspace, current_step)
                call gradient_subsystem(spin_rotated, subsys_indices, current_gradient)
                call vec_to_tangentspace(spin_rotated,subsys_indices, current_gradient, current_gradient_tspace)
                back_displacement = -1.0d0 * current_step * vector_tspace
                call rotate_vector_between_tangent_spaces(spin_rotated,subsys_indices, back_displacement, &
                        & current_gradient_tspace, current_gradient_tspace_initial)
                diff_hessianvec = diff_hessianvec + current_gradient_tspace_initial * coefficients(i_order)
            end do
            diff_hessianvec = diff_hessianvec / (step*1)

        end subroutine fd_hessvec_forward_orderN_subsys

end module spk_finites
