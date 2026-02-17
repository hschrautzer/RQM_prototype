module spk_lbfgs
    !> limited-memory Broyden–Fletcher–Goldfarb–Shanno (LBFGS). The method is described in:
    !> Ivanov, Aleksei V., Valery M. Uzdin, and Hannes Jónsson.
    !> "Fast and robust algorithm for energy minimization of spin systems applied in an analysis of high temperature
    !> spin configurations in terms of skyrmion density." Computer Physics Communications 260 (2021): 107749.
    !> It can be used in different optimization problems like minimum mode following or simple minimization.
    use spk_precision_const, only: xp
    use spk_input_lattice, only: N_atom
    use spk_vector_lina, only: norm
    use spk_rotations, only: rotate_vector_between_tangent_spaces, rotate_spinconfiguration
    implicit none
    private
    interface lbfgs_step
        !> Performs lbfgs integration. Interface decides whether vpo is applied to the configuration
        !> space build which represents the combined configuration space of a series of configurations (images) or to
        !> a single spin configuration.
        module procedure lbfgs_step_single_configuration
        module procedure lbfgs_step_single_configuration_subsystem
        module procedure lbfgs_step_path
        module procedure lbfgs_step_noncurved
    end interface lbfgs_step

    public :: lbfgs_step
    contains
        subroutine lbfgs_step_single_configuration(spin, iteration, lbfgs_memory, force, force_previous, step_previous, &
                & steplength_previous, rho, gamma, d, y, theta_max, new_step, stepl, i_rotatespin)
            !================================Input Variable=============================================================
            integer, intent(in) :: iteration                        !< current iteration
                                                                    !< (needed to access correct memory quantity)
            integer, intent(in) :: lbfgs_memory                     !< number of memory steps
            real(kind=xp), intent(in) :: force(3*N_atom)            !< force of current configuration
                                                                    !< (in tangent space of curr. configuration)
            real(kind=xp), intent(in) :: force_previous(3*N_atom)   !< force of previous configuration
                                                                    !< (in tangent space of prev. configuration)
            real(kind=xp), intent(in) :: step_previous(3*N_atom)    !< previous step in 3N tangent space of
                                                                    !< previous configuration
            real(kind=xp), intent(in) :: steplength_previous        !< steplength (scaling factor) of previous step
            real(kind=xp), intent(in) :: theta_max
            logical, optional, intent(in) :: i_rotatespin           !< whether spin configuration is rotated.
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: new_step(3*N_atom)        !< new step
            real(kind=xp), intent(out) :: stepl                     !< steplength
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: rho(lbfgs_memory)       !< memorized rho
            real(kind=xp), intent(inout) :: gamma(lbfgs_memory)     !< memorized rho
            real(kind=xp), intent(inout) :: d(3*N_atom,lbfgs_memory)!< memorized difference in spin configurations in 2N
            real(kind=xp), intent(inout) :: y(3*N_atom,lbfgs_memory)!< memorized difference in gradient in 2N
            real(kind=xp), intent(inout) :: spin(3*N_atom)          !< Inout spin configuration
            !================================Local Variable=============================================================
            integer :: k                                            !< local copy of iteration count
            integer :: l,n                                          !< iterator through memory quantities
            integer :: i,j                                          !< iterator
            real(kind=xp) :: step_previous_dumy(3*N_atom)           !< local copy of previous step in 3N tangent space
                                                                    !< of previous iter multiplied by previous step length
            real(kind=xp) :: step_previous_r(3*N_atom)              !< local copy of previous step in 3N tangent space
                                                                    !< of current configuration
            real(kind=xp) :: force_previous_r(3*N_atom)             !< local copy of previous step in 3N tangent space
                                                                    !< of current configuration
            real(kind=xp) :: dumy_mem(3*N_atom)                     !<dummy memory quantity vector
            real(kind=xp) :: res, theta_rms
            real(kind=xp) :: q(3*N_atom), dummy_step(3*N_atom)
            logical :: local_i_spin_rotation

            if (present(i_rotatespin)) then
                local_i_spin_rotation = i_rotatespin
            else
                local_i_spin_rotation = .True.
            end if

            k = iteration - 1
    444     n = modulo(k, lbfgs_memory) + 1
            if (k==0) then
                new_step = force
                rho=0.0d0
                d=0.0d0
                y=0.0d0
                gamma=0.0d0
            else
                step_previous_dumy = steplength_previous * step_previous
                ! Rotate the previous step and the perp gradient (which are defined in the tangent space of the spin
                ! of the current spin configuration) towards the current spin configuration tangent space
                call rotate_vector_between_tangent_spaces(spin, step_previous_dumy, step_previous_dumy, step_previous_r)
                call rotate_vector_between_tangent_spaces(spin, step_previous_dumy, force_previous, force_previous_r)
                d(:,n)=step_previous_r
                y(:,n)=force_previous_r - force !switch representation between force and gradient: (*-1.0d0)
                do l=1,lbfgs_memory
                    if (l==n) then
                        cycle
                    end if
                    call rotate_vector_between_tangent_spaces(spin, step_previous_dumy,d(:,l),dumy_mem)
                    d(:,l)=dumy_mem
                    call rotate_vector_between_tangent_spaces(spin, step_previous_dumy,y(:,l),dumy_mem)
                    y(:,l)=dumy_mem
                end do
                res = 0.0d0
                do i=1,3*N_atom
                    res = res + y(i,n) * d(i,n)
                end do
                rho(n)=1.0d0 / res

                if (rho(n)<0.0d0) then
                    k=0
                    goto 444
                end if
                q = -1.0d0 * force
                do l=lbfgs_memory,1,-1
                    j = modulo(l+n,lbfgs_memory) + 1
                    res = 0.0d0
                    do i=1,3*N_atom
                        res = res + q(i) * d(i,j)
                    end do
                    gamma(j)=rho(j) * res
                    q(:)=q(:) - gamma(j) * y(:,j)
                end do

                res = 0.0d0
                do i=1,3*N_atom
                    res = res + y(i,n) * y(i,n)
                end do

                dummy_step(:) = q(:) / (rho(n) * res)
                do l=1,lbfgs_memory
                    if (k<lbfgs_memory) then
                        j = l
                    else
                        j=modulo(l+n,lbfgs_memory)+1
                    end if
                    res = 0.0d0
                    do i=1,3*N_atom
                        res = res + y(i,j) * dummy_step(i)
                    end do
                    dummy_step(:)=dummy_step(:) + d(:,j) * (gamma(j) - rho(j) * (res))
                end do
                new_step = -1.0d0 * dummy_step
            end if

            !==================ROTATION AND SCALING=====================================================================
            !-Calculate root means square of displacement angles
            theta_rms=0
            do i=1,3*N_atom
                theta_rms = theta_rms + new_step(i) * new_step(i)
            end do
            theta_rms=sqrt(theta_rms)
            !-If rms angle is bigger than the maximum rotation angle defined by user scale steplength
            if ((theta_max/theta_rms)<1) then
                stepl=(theta_max / theta_rms)
            else
                stepl=1.0d0
            end if

            if (local_i_spin_rotation) call rotate_spinconfiguration(spin, new_step, stepl)

        end subroutine lbfgs_step_single_configuration

        subroutine lbfgs_step_single_configuration_subsystem(spin, subsys_indices, iteration, lbfgs_memory, force, &
                & force_previous, step_previous, steplength_previous, rho, gamma, d, y, theta_max, new_step, stepl, &
                & i_rotatespin)
            !================================Input Variable=============================================================
            integer, intent(in) :: iteration                        !< current iteration
                                                                    !< (needed to access correct memory quantity)
            integer, intent(in) :: lbfgs_memory                     !< number of memory steps
            integer, intent(in) :: subsys_indices(:)                !< Indices of the atoms in the subsystem
            real(kind=xp), intent(in) :: force(:)                   !< Force of current configuration (subsystem)
                                                                    !< (in tangent space of curr. configuration)
            real(kind=xp), intent(in) :: force_previous(:)          !< Force of previous configuration (subsystem)
                                                                    !< (in tangent space of prev. configuration)
            real(kind=xp), intent(in) :: step_previous(:)           !< previous step in 3N tangent space of
                                                                    !< previous configuration (subsystem)
            real(kind=xp), intent(in) :: steplength_previous        !< steplength (scaling factor) of previous step
            real(kind=xp), intent(in) :: theta_max
            logical, optional, intent(in) :: i_rotatespin           !< whether spin configuration is rotated.
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: new_step(:)               !< new step (subsystem)
            real(kind=xp), intent(out) :: stepl                     !< steplength
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: rho(:)                  !< memorized rho (length equals memory)
            real(kind=xp), intent(inout) :: gamma(:)                !< memorized gamma (length equals memory)
            real(kind=xp), intent(inout) :: d(:,:)                  !< mem. diff. in configs in 2N (subsystem)
            real(kind=xp), intent(inout) :: y(:,:)                  !< mem. diff. in gradient in 2N (subsystem)
            real(kind=xp), intent(inout) :: spin(:)                 !< Inout spin configuration (full lattice)
            !================================Local Variable=============================================================
            integer :: k                                            !< Local copy of iteration count
            integer :: l,n                                          !< Iterator through memory quantities
            integer :: i,j                                          !< Iterator
            real(kind=xp), allocatable :: step_previous_dumy(:)     !< Local copy of previous step in 3N tangent space
                                                                    !< of previous iter multiplied by previous step length
            real(kind=xp), allocatable :: step_previous_r(:)        !< Local copy of previous step in 3N tangent space
                                                                    !< of current configuration
            real(kind=xp), allocatable :: force_previous_r(:)       !< Local copy of previous step in 3N tangent space
                                                                    !< of current configuration
            real(kind=xp), allocatable :: dumy_mem(:)               !< Dummy memory quantity vector
            real(kind=xp) :: res, theta_rms
            real(kind=xp), allocatable :: q(:), dummy_step(:)
            logical :: local_i_spin_rotation
            integer :: N_sub                                        !< Number of atoms in the subsystem

            if (present(i_rotatespin)) then
                local_i_spin_rotation = i_rotatespin
            else
                local_i_spin_rotation = .True.
            end if

            N_sub = size(subsys_indices)
            allocate(step_previous_dumy(3*N_sub),step_previous_r(3*N_sub),force_previous_r(3*N_sub),dumy_mem(3*N_sub),&
             & q(3*N_sub), dummy_step(3*N_sub))

            k = iteration - 1
    555     n = modulo(k, lbfgs_memory) + 1
            if (k==0) then
                new_step = force
                rho=0.0d0
                d=0.0d0
                y=0.0d0
                gamma=0.0d0
            else
                step_previous_dumy = steplength_previous * step_previous
                ! Rotate the previous step and the perp gradient (which are defined in the tangent space of the spin
                ! of the current spin configuration) towards the current spin configuration tangent space
                call rotate_vector_between_tangent_spaces(spin, subsys_indices, step_previous_dumy, &
                            & step_previous_dumy, step_previous_r)
                call rotate_vector_between_tangent_spaces(spin, subsys_indices, step_previous_dumy, &
                            & force_previous, force_previous_r)
                d(:,n)=step_previous_r
                y(:,n)=force_previous_r - force !switch representation between force and gradient: (*-1.0d0)
                !Now we need to rotate all the other memory vectors to the current coordinate frame (tangent space)
                do l=1,lbfgs_memory
                    if (l==n) then
                        cycle
                    end if
                    call rotate_vector_between_tangent_spaces(spin, subsys_indices, step_previous_dumy, &
                                & d(:,l),dumy_mem)
                    d(:,l)=dumy_mem
                    call rotate_vector_between_tangent_spaces(spin, subsys_indices, step_previous_dumy, &
                                & y(:,l),dumy_mem)
                    y(:,l)=dumy_mem
                end do
                res = 0.0d0
                do i=1,3*N_sub
                    res = res + y(i,n) * d(i,n)
                end do
                rho(n)=1.0d0 / res

                if (rho(n)<0.0d0) then
                    k=0
                    goto 555
                end if
                q = -1.0d0 * force
                do l=lbfgs_memory,1,-1
                    j = modulo(l+n,lbfgs_memory) + 1
                    res = 0.0d0
                    do i=1,3*N_sub
                        res = res + q(i) * d(i,j)
                    end do
                    gamma(j)=rho(j) * res
                    q(:)=q(:) - gamma(j) * y(:,j)
                end do

                res = 0.0d0
                do i=1,3*N_sub
                    res = res + y(i,n) * y(i,n)
                end do

                dummy_step(:) = q(:) / (rho(n) * res)
                do l=1,lbfgs_memory
                    if (k<lbfgs_memory) then
                        j = l
                    else
                        j=modulo(l+n,lbfgs_memory)+1
                    end if
                    res = 0.0d0
                    do i=1,3*N_sub
                        res = res + y(i,j) * dummy_step(i)
                    end do
                    dummy_step(:)=dummy_step(:) + d(:,j) * (gamma(j) - rho(j) * (res))
                end do
                new_step = -1.0d0 * dummy_step
            end if

            !==================ROTATION AND SCALING=====================================================================
            !-Calculate root means square of displacement angles
            theta_rms=0
            do i=1,3*N_sub
                theta_rms = theta_rms + new_step(i) * new_step(i)
            end do
            theta_rms=sqrt(theta_rms)
            !-If rms angle is bigger than the maximum rotation angle defined by user scale steplength
            if ((theta_max/theta_rms)<1) then
                stepl=(theta_max / theta_rms)
            else
                stepl=1.0d0
            end if

            if (local_i_spin_rotation) call rotate_spinconfiguration(spin,subsys_indices,new_step*stepl)

        end subroutine lbfgs_step_single_configuration_subsystem

        subroutine lbfgs_step_noncurved(iteration, lbfgs_memory, force, force_previous, step_previous, &
                & steplength_previous, rho, gamma, d, y, theta_max, new_step, stepl)
            !> Version of LBFGS without rotations
            !================================Input Variable=============================================================
            integer, intent(in) :: iteration                        !< Current iteration of the optimization algorithm
                                                                    !< (needed to access correct memory quantity)
            integer, intent(in) :: lbfgs_memory                     !< Number of steps saved in Memory
            real(kind=xp), intent(in) :: force(:)                   !< Force of current configuration
                                                                    !< (in tangent space of curr. configuration)
            real(kind=xp), intent(in) :: force_previous(:)          !< Force of configuration in previous iteration
                                                                    !< (in tangent space of prev. configuration)
            real(kind=xp), intent(in) :: step_previous(:)           !< Step of previous iteration in 3N tangent space of
                                                                    !< previous configuration
            real(kind=xp), intent(in) :: steplength_previous        !< Steplength (scaling factor) of step of prev. iter.
            real(kind=xp), intent(in) :: theta_max                  !< Maximum Rotation angle of individual component
                                                                    !< Alternative to Line Search
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: new_step(:)               !< New step (Displacement)
            real(kind=xp), intent(out) :: stepl                     !< Steplength
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: rho(lbfgs_memory)       !< memorized rho
            real(kind=xp), intent(inout) :: gamma(lbfgs_memory)     !< memorized rho
            real(kind=xp), intent(inout) :: d(:,:)                  !< memorized difference in spin configurations in 2N
            real(kind=xp), intent(inout) :: y(:,:)                  !< memorized difference in gradient in 2N
            !================================Local Variable=============================================================
            integer :: la                                           !< Length of the Arrays
            integer :: k                                            !< Local copy of iteration count
            integer :: l,n                                          !< Iterator through memory quantities
            integer :: j                                            !< Iterators
            real(kind=xp), allocatable :: step_previous_dumy(:)     !< Local copy of previous step in 3N tangent space
                                                                    !< of previous iter multiplied by previous step length
            real(kind=xp) :: res, theta_rms
            real(kind=xp), allocatable :: q(:), dummy_step(:)

            la = size(force)
            allocate(step_previous_dumy(la),q(la),dummy_step(la))

            k = iteration - 1
    555     n = modulo(k, lbfgs_memory) + 1
            if (k==0) then
                new_step = force
                rho=0.0d0
                d=0.0d0
                y=0.0d0
                gamma=0.0d0
            else
                step_previous_dumy = steplength_previous * step_previous
                ! Rotate the previous step and the perp gradient (which are defined in the tangent space of the spin
                ! of the current spin configuration) towards the current spin configuration tangent space

                d(:,n)=step_previous_dumy
                y(:,n)=force_previous - force !switch representation between force and gradient: (*-1.0d0)

                res = DOT_PRODUCT(y(:,n),d(:,n))
                rho(n)=1.0d0 / res

                if (rho(n)<0.0d0) then
                    k=0
                    goto 555
                end if
                q = -1.0d0 * force

                do l=lbfgs_memory,1,-1
                    j = modulo(l+n,lbfgs_memory) + 1
                    res = DOT_PRODUCT(q,d(:,j))
                    gamma(j)=rho(j) * res
                    q(:)=q(:) - gamma(j) * y(:,j)
                end do

                res = DOT_PRODUCT(y(:,n),y(:,n))

                dummy_step(:) = q(:) / (rho(n) * res)

                do l=1,lbfgs_memory
                    if (k<lbfgs_memory) then
                        j = l
                    else
                        j=modulo(l+n,lbfgs_memory)+1
                    end if
                    res = 0.0d0
                    res = DOT_PRODUCT(y(:,j),dummy_step)
                    dummy_step(:)=dummy_step(:) + d(:,j) * (gamma(j) - rho(j) * (res))
                end do
                new_step = -1.0d0 * dummy_step
            end if

            !==================ROTATION AND SCALING=====================================================================
            !-Calculate root means square of displacement angles
            theta_rms=DOT_PRODUCT(new_step,new_step)
            theta_rms=sqrt(theta_rms)
            !-If rms angle is bigger than the maximum rotation angle defined by user scale steplength
            if ((theta_max/theta_rms)<1) then
                stepl=(theta_max / theta_rms)
            else
                stepl=1.0d0
            end if
        end subroutine lbfgs_step_noncurved

        subroutine lbfgs_step_path(n_image,path, iteration, lbfgs_memory, forces, forces_previous, steps_previous, &
                & steplength_previous, rho, gamma, ds, ys, theta_max, new_steps, stepl)
            !================================Input Variable=============================================================
            integer, intent(in) :: iteration                        !< current iteration
                                                                    !< (needed to access correct memory quantity)
            integer, intent(in) :: n_image                          !< number of images forming the path
            integer, intent(in) :: lbfgs_memory                     !< number of memory steps
            real(kind=xp), intent(in) :: forces(3*N_atom,n_image)   !< forces of all images
            real(kind=xp), intent(in) :: forces_previous(3*N_atom,n_image)
                                                                    !< force of previous configuration
            real(kind=xp), intent(in) :: steps_previous(3*N_atom,n_image)
                                                                    !< previous step in 3N tangent space of
                                                                    !< previous configuration
            real(kind=xp), intent(in) :: steplength_previous        !< steplength (scaling factor) of previous step
            real(kind=xp), intent(in) :: theta_max
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: new_steps(3*N_atom,n_image)
                                                                    !< new step
            real(kind=xp), intent(out) :: stepl                     !< steplength
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: rho(lbfgs_memory)       !< memorized rho
            real(kind=xp), intent(inout) :: gamma(lbfgs_memory)     !< memorized rho
            real(kind=xp), intent(inout) :: ds(3*N_atom,lbfgs_memory,n_image)
                                                                    !< memorized difference in spin configurations in 2N
            real(kind=xp), intent(inout) :: ys(3*N_atom,lbfgs_memory,n_image)
                                                                    !< memorized difference in gradient in 2N
            real(kind=xp), intent(inout) :: path(3*N_atom,n_image)  !< Inout spin configuration
            !================================Local Variable=============================================================
            integer :: k                                            !< local copy of iteration count
            integer :: l,n                                          !< iterator through memory quantities
            integer :: i,j                                          !< iterator
            integer :: i_image                                      !< iterator for images
            real(kind=xp) :: step_previous_dumy(3*N_atom)           !< local copy of previous step in 3N tangent space
                                                                    !< of previous iter multiplied by previous step length
            real(kind=xp) :: step_previous_r(3*N_atom)              !< local copy of previous step in 3N tangent space
                                                                    !< of current configuration
            real(kind=xp) :: force_previous_r(3*N_atom)             !< local copy of previous step in 3N tangent space
                                                                    !< of current configuration
            real(kind=xp) :: dumy_mem(3*N_atom)                     !<dummy memory quantity vector
            real(kind=xp) :: res, theta_rms
            real(kind=xp) :: qs(3*N_atom,n_image)

            k = iteration - 1
    666     n = modulo(k, lbfgs_memory) + 1
            if (k==0) then
                new_steps = forces
                rho=0.0d0
                ds=0.0d0
                ys=0.0d0
                gamma=0.0d0
            else
                !$OMP PARALLEL DO PRIVATE(i_image,step_previous_dumy,step_previous_r,force_previous_r,l,dumy_mem) SHARED(n)
                do i_image=1,n_image
                    step_previous_dumy = steplength_previous * steps_previous(:,i_image)
                    ! Rotate the previous step and the perp gradient (which are defined in the tangent space of the spin
                    ! of the current spin configuration) towards the current spin configuration tangent space
                    call rotate_vector_between_tangent_spaces(path(:,i_image), step_previous_dumy, step_previous_dumy, &
                            & step_previous_r)
                    call rotate_vector_between_tangent_spaces(path(:,i_image), step_previous_dumy, &
                            & forces_previous(:,i_image), force_previous_r)
                    ds(:,n,i_image) = step_previous_r
                    ys(:,n,i_image) = force_previous_r - forces(:,i_image)
                    !Now we need to rotate all the other memory vectors to the current coordinate frame (tangent space)
                    do l=1,lbfgs_memory
                        if (l==n) then
                            cycle
                        end if
                        call rotate_vector_between_tangent_spaces(path(:,i_image), step_previous_dumy,ds(:,l,i_image),dumy_mem)
                        ds(:,l,i_image)=dumy_mem
                        call rotate_vector_between_tangent_spaces(path(:,i_image), step_previous_dumy,ys(:,l,i_image),dumy_mem)
                        ys(:,l,i_image)=dumy_mem
                    end do
                end do
                !$OMP END PARALLEL DO
                !Calculate dot products seriell in order to achieve global memory quantities for the whole path
                res = 0.0d0
                do i_image=1,n_image
                    do i=1,3*N_atom
                        res = res + ys(i,n,i_image) * ds(i,n,i_image)
                    end do
                end do
                rho(n)=1.0d0 / res

                if (rho(n)<0.0d0) then
                    k=0
                    goto 666
                end if

                do i_image=1,n_image
                    qs(:,i_image) = -1.0d0 * forces(:,i_image)
                end do

                do l=lbfgs_memory,1,-1
                    j = modulo(l+n,lbfgs_memory) + 1
                    res = 0.0d0
                    do i_image=1,n_image
                        do i=1,3*N_atom
                            res = res + qs(i,i_image) * ds(i,j,i_image)
                        end do
                    end do
                    gamma(j)=rho(j) * res
                    do i_image=1,n_image
                        qs(:,i_image)=qs(:,i_image) - gamma(j) * ys(:,j,i_image)
                    end do
                end do

                res = 0.0d0
                do i_image=1,n_image
                    do i=1,3*N_atom
                        res = res + ys(i,n,i_image) * ys(i,n,i_image)
                    end do
                end do

                do i_image=1,n_image
                    new_steps(:,i_image) = qs(:,i_image) / (rho(n) * res)
                end do


                do l=1,lbfgs_memory
                    if (k<lbfgs_memory) then
                        j = l
                    else
                        j=modulo(l+n,lbfgs_memory)+1
                    end if
                    res = 0.0d0
                    do i_image=1,n_image
                        do i=1,3*N_atom
                            res = res + ys(i,j,i_image) * new_steps(i,i_image)
                        end do
                    end do
                    do i_image=1,n_image
                        new_steps(:,i_image)=new_steps(:,i_image) + ds(:,j,i_image) * (gamma(j) - rho(j) * (res))
                    end do
                end do
                new_steps = -1.0d0 * new_steps
            end if

            !==================ROTATION AND SCALING=====================================================================
            !-Calculate root means square of displacement angles
            theta_rms=0
            do i_image=1,n_image
                do i=1,3*N_atom
                    theta_rms = theta_rms + new_steps(i,i_image) * new_steps(i,i_image)
                end do
            end do
            theta_rms=sqrt(theta_rms)

            !-If rms angle is bigger than the maximum rotation angle defined by user scale steplength
            if ((theta_max/theta_rms)<1) then
                stepl=(theta_max / theta_rms)
            else
                stepl=1.0d0
            end if


            !$OMP PARALLEL DO PRIVATE(i_image) SHARED(stepl)
            do i_image=1,n_image
                call rotate_spinconfiguration(path(:,i_image), new_steps(:,i_image), stepl)
            end do
            !$OMP END PARALLEL DO

        end subroutine lbfgs_step_path

end module spk_lbfgs
