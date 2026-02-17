module spk_rqm
    !> Module for Rayleigh Quotient Minimization
    use spk_eigensolverparams, only: rqm_params
    use spk_precision_const, only: xp
    use spk_logging, only: log
    use spk_algoconst, only: INTEGRATOR_EULER, INTEGRATOR_VPO, INTEGRATOR_LBFGS, FD_SCHEME_FORWARD, FD_SCHEME_CENTRAL, &
            & FD_SCHEME_BACKWARD
    use spk_abbreviations_const, only: C_NOTCODED_ERROR
    use spk_tangentspace, only: vec_to_tangentspace, vec_to_embeddingspace, to_tangentspace_grassmann_manifold
    use spk_finites, only: fd_hessvec_forward
    use spk_retraction, only: GM_retract_qr, GM_retract_exp
    use spk_vectortransport, only: GM_paralleltransport_calcTM, GM_paralleltransport_applyTM
    use spk_input_lattice, only: N_atom
    use spk_euler, only: euler_step
    use spk_vector_lina, only: norm
    use spk_vpo, only: vpo_step, vpo_update_velocity
    use spk_products, only: innerproduct_frobenius
    use spk_lbfgs, only: lbfgs_step
    use lapack95, only: syev
    use omp_lib
    use spk_casting, only: int_to_str
    implicit none
    private
    interface rqm
        module procedure rqm_grassmann_performance
        module procedure rqm_grassmann_information
        module procedure rqm_grassmann_subsystem_performance
        module procedure rqm_grassmann_subsystem_information
    end interface rqm
    public :: rqm
    contains
        subroutine rqm_grassmann_performance(spin,basis,p,v_ini,v_fin,rqm_settings,evals,rqm_iterations)
            !> This routine performs an optimization routine on the Grassmann Manifold with the objective function
            !> tr X^T H X, where H is the Hessian of the spin configuration and X is a 2N x p matrix representative
            !> for an equivalence class of matrices on the Stiefel manifold corresponding to one point on the
            !> Grassmannian. The subspace spanned by the columns of X_min is the invariant subspace belonging to
            !> the p smallest eigenvalues. This version of the routine only computes the rayleigh ritz procedure once
            !> at the end.
            !===================================Input Variable======================================================
            real(kind=xp), intent(in) :: spin(:)                                    !< Spin Configuration
            real(kind=xp), intent(in) :: basis(:,:,:)                               !< Basis of Tangent Space of Spin
            integer, intent(in) :: p                                                !< Dimension of the subspace
            real(kind=xp), intent(in) :: v_ini(:,:)                                 !< Start Vector for RQM (3N)
            type(rqm_params), intent(in) :: rqm_settings                            !< Parameter Container for RQM
            !================================Output Variable========================================================
            real(kind=xp), intent(out) :: v_fin(:,:)                                !< Result Vector of RQM (3N)
            real(kind=xp), intent(out) :: evals(:)                                  !< Result of the two lowest evals
            integer, optional, intent(out) :: rqm_iterations                        !< Number RQM iterations
            !================================Local Variable=========================================================
            real(kind=xp), allocatable :: X(:,:),X_previous(:,:)                    !< Points on Grassmannian
            integer :: iter                                                         !< Iteration Counter
            integer :: p_iter, p_iter2                                              !< Iteration for Subspace Dim.
            real(kind=xp), allocatable :: dummyvec3N(:),dummyvec3N_2(:),dummyvec2N(:)!< Dummy Vectors
            real(kind=xp) :: dummyval                                               !< Dummy Value
            real(kind=xp) :: rq_matrix(p,p)                                         !< Matrix X^T * H * X
            real(kind=xp), allocatable :: HX_findiff(:)                             !< Finite difference approximation
                                                                                        !< of the matrix-vector product
            real(kind=xp), allocatable :: rq_gradient(:)                            !< Gradient of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_gradient_2d(:,:)                       !< Matrix representation
            real(kind=xp), allocatable :: rq_force(:)                               !< Force of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_step(:)                                !< Step for updating the RQ
            real(kind=xp), allocatable :: rq_step_2d(:,:)                           !< Step for updating the X quantity
            real(kind=xp) :: rq_gradient_norm                                       !< Norm of the Ralyeigh Gradient G
                                                                                    !< sqrt(tr(G^T G))
            real(kind=xp) :: rq
            real(kind=xp) :: rayleighritz_evec(p,p)                                 !< Eigenvalues and Eigenvector of
                                                                                        !< Rayleigh Ritz Procedure
            real(kind=xp), allocatable :: U(:,:), S(:), VT(:,:)                     !< thin SVD of rq_step_2d
            real(kind=xp), allocatable :: T(:,:)                                    !< Parallel transport matrix
            !------------------------------FD Information---------------------------------------------------------------
            real(kind=xp) :: fd_step_used(p)                                        !< Finite Difference Step used (only
                                                                                    !< needed in case of Richardson-FD)
            !------------------------------Algorithm Depending Variables------------------------------------------------
            ! Depending on which algorithm is chosen one might need different local variable, they will be only
            ! allocated if needed. In this case all the solver quantities have a larger dimension depending on the size
            ! of the invariant subspace
            ! .................................................VPO......................................................
            real(kind=xp), allocatable :: vpo_velocity(:)                           !< Velocity in Velocity-Projection-Algo.
            real(kind=xp), allocatable :: vpo_b(:)                                  !< Local vpo quantity
            ! .................................................LBFGS....................................................
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_d
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_y
            real(kind=xp), allocatable, dimension(:) :: lbfgs_rho
            real(kind=xp), allocatable, dimension(:) :: lbfgs_gamma
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_force
            real(kind=xp), allocatable, dimension(:) :: lbfgs_current_step
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_step
            real(kind=xp) :: lbfgs_previous_steplength, lbfgs_steplength
            integer :: i_m
            ! Linesearch Variables
            real(kind=xp) :: leftside, rightside, alpha, rho, c
            real(kind=xp), allocatable :: lbfgs_Hstep(:,:),lbfgs_linesearchdummy(:,:)
            !.............................Local Variable Allocation.....................................................
            ! These local variables are made allocatable, since they can reach sizes beyond stack limit (heap arrays)
            allocate(dummyvec2N(2*N_atom))                      ! Dummy 2N vector
            allocate(dummyvec3N(3*N_atom))                      ! Dummy 3N vector
            allocate(dummyvec3N_2(3*N_atom))                    ! Dummy 3N vector
            allocate(X(2*N_atom,p),X_previous(2*N_atom,p))      ! Optimization quantity
            allocate(HX_findiff(2*N_atom*p))                    ! Store the 2N x p Matrix as 1D array
            allocate(rq_gradient(2*N_atom*p))                   ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_gradient_2d(2*N_atom,p))                ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_force(2*N_atom*p))                      ! Forcematrix (2N x p) of the RQ as 1D
            allocate(rq_step(2*N_atom*p))                       ! Stepmatrix (2N x p) of the RQ as 1D
            allocate(rq_step_2d(2*N_atom,p))                    ! Stepmatrix (2N x p) of the RQ as 2D
            allocate(T(2*N_atom,p),U(2*N_atom,p),VT(p,p),S(p))  ! Transport matrix and thin SVD of rq_step_2d
            !.............................Prepare Solver............................................................
            select case(rqm_settings%rqm_solver%solver)
                case(INTEGRATOR_VPO)
                    allocate(vpo_velocity(2*N_atom*p),vpo_b(2*N_atom*p))
                    vpo_velocity = 0.0d0
                    vpo_b = 0.0d0
                case(INTEGRATOR_LBFGS)
                    allocate(lbfgs_d(2*N_atom*p,rqm_settings%rqm_solver%lbfgs_memory),lbfgs_y(2*N_atom*p,&
                            & rqm_settings%rqm_solver%lbfgs_memory), lbfgs_rho(rqm_settings%rqm_solver%lbfgs_memory), &
                            & lbfgs_gamma(rqm_settings%rqm_solver%lbfgs_memory))
                    allocate(lbfgs_previous_force(2*N_atom*p), lbfgs_current_step(2*N_atom*p),&
                            & lbfgs_previous_step(2*N_atom*p))
                    allocate(lbfgs_Hstep(2*N_atom,p),lbfgs_linesearchdummy(2*N_atom,p))
                    lbfgs_steplength = 1.0d0
                    lbfgs_previous_force = 0.0d0
                    rho = 0.5d0
                    c = 0.0001d0
            end select
            !.............................Iteration Zero (Preparation)..................................................
            if (present(rqm_iterations)) rqm_iterations = 0
            v_fin = 0.0d0
            rq = 0.0d0
            X = 0.0d0
            evals = 0.0d0
            ! Initialize 2N x p matrix X
            do p_iter=1,p
                call vec_to_tangentspace(basis,v_ini(:,p_iter),dummyvec2N)
                X(:,p_iter) = dummyvec2N / norm(dummyvec2N)
            end do

            call GM_retract_qr(X,"G")

            ! Compute the Hessian-vector product H * X (in R^{2Nxp}), which is needed for computed the Gradient
            ! of the Grassmanian as well as evaluating the objective function
            if (rqm_settings%fd_settings%i_richardson) then
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                    & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                    & rqm_settings%fd_settings%richardson_iter)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            fd_step_used(p) = dummyval
                            HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                end select
            else
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, &
                                    & rqm_settings%fd_settings%order, dummyvec3N, dummyvec3N_2)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    end select
            end if

            ! Compute Rayleigh Matrix X^T H X
            rq_matrix = 0.0d0
            do p_iter=1,p
                do p_iter2=1,p
                    ! col idx, row idx
                    rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N_atom*(p_iter2-1)+1:2*N_atom*p_iter2))
                end do
            end do
            ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
            rq = 0.0d0
            do p_iter=1,p
                rq = rq + rq_matrix(p_iter,p_iter)
            end do

            ! The gradient is given by 2 H X - 2 X (X^T H X) and is a 2N x p Matrix
            rq_gradient = HX_findiff
            do p_iter=1,p
                do p_iter2=1,p
                    rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) &
                            & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                end do
            end do
            rq_gradient = rq_gradient * 2.0d0

            ! Change representation of Gradient to a 2N x p matrix:
            do p_iter=1,p
                rq_gradient_2d(:,p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
            end do

            ! Compute the norm of the Gradient which is given by the Frobenius Norm:
            ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
            rq_gradient_norm = 0.0d0
            do p_iter=1,p
                rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter),&
                        & rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter))
            end do
            rq_gradient_norm = sqrt(rq_gradient_norm)
            ! For our solvers we need the force representation
            rq_force = -1.0d0 * rq_gradient

            do iter=1,rqm_settings%rqm_solver%n_steps
                !----------------------------Check convergence crit-----------------------------------------------------
                if (rq_gradient_norm<=rqm_settings%rqm_solver%conv_crit) then
                    !write(*,*) 'Rayleigh optimization converged at step: ', trim(adjustl(int_to_str(iter)))
                    if (present(rqm_iterations)) rqm_iterations = iter
                    exit
                end if
                if (iter == rqm_settings%rqm_solver%n_steps) then
                    if (present(rqm_iterations)) rqm_iterations = iter
                    write(*,*) "RQM exceeded iterations"
                end if
                !----------------------------Solver---------------------------------------------------------------------
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_EULER)
                        call euler_step(rq_force,rqm_settings%rqm_solver%euler_dt,rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
                        end do
                    case(INTEGRATOR_VPO)
                        call vpo_step(vpo_velocity,rq_force, vpo_b, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt, rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
                        end do
                    case(INTEGRATOR_LBFGS)
                        ! Calculate the step and rotate the spin configuration
                        call lbfgs_step(iter, rqm_settings%rqm_solver%lbfgs_memory,  rq_force, lbfgs_previous_force, &
                                & lbfgs_previous_step, lbfgs_previous_steplength, lbfgs_rho, lbfgs_gamma, lbfgs_d, lbfgs_y, &
                                & rqm_settings%rqm_solver%lbfgs_theta_max, rq_step, lbfgs_steplength)
                        lbfgs_previous_force =  rq_force
                        lbfgs_previous_steplength = lbfgs_steplength
                        lbfgs_previous_step = rq_step
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) * lbfgs_steplength
                        end do
                        ! Simple Armijo Linesearch
                        if (.False.) then
                            alpha = 1.0d0
                            rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                            rightside = rightside *c*alpha
                            !Left side computation
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis, rq_step_2d(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                        & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                        & rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                lbfgs_Hstep(:,p_iter) = dummyvec2N
                            end do
                            lbfgs_linesearchdummy = 2.0*alpha*X+rq_step_2d*alpha**2.0
                            call GM_retract_qr(lbfgs_linesearchdummy,"G")
                            leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                            i_m = 0
                            do while(leftside>rightside)
                                i_m = i_m + 1
                                if (i_m==10) then
                                    exit
                                end if
                                write(*,*) leftside, " > " ,rightside
                                alpha = alpha * rho
                                rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                                rightside = rightside *c*alpha
                                !Left side computation
                                do p_iter=1,p
                                    call vec_to_embeddingspace(basis, rq_step_2d(:,p_iter), dummyvec3N)
                                    call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                            & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                            & rqm_settings%fd_settings%richardson_iter)
                                    call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                    lbfgs_Hstep(:,p_iter) = dummyvec2N
                                end do
                                lbfgs_linesearchdummy = 2.0*alpha*X+alpha**2*rq_step_2d
                                call GM_retract_qr(lbfgs_linesearchdummy,"G")
                                leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                                write(*,*) "alpha: ", alpha
                            end do
                            write(*,*) "finished LS:", alpha
                            lbfgs_previous_steplength = alpha
                            lbfgs_steplength = alpha
                            rq_step_2d = rq_step_2d * lbfgs_steplength
                        end if
                    end select

                ! Save the current iterate X_k
                X_previous = X
                ! Compute the next iterate X_k+1 by moving along the geodesic defined by the step and the exp. map
                call GM_retract_exp(X,rq_step_2d,1.0d0,U,S,VT,"Y")
                !X = X_previous + rq_step_2d

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                    case(INTEGRATOR_LBFGS)
                        if(.True.)then
                            call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                            call GM_paralleltransport_applyTM(lbfgs_previous_step,T,U)
                            call GM_paralleltransport_applyTM(lbfgs_previous_force,T,U)
                            do i_m=1,rqm_settings%rqm_solver%lbfgs_memory
                                call GM_paralleltransport_applyTM(lbfgs_d(:,i_m),T,U)
                                call GM_paralleltransport_applyTM(lbfgs_y(:,i_m),T,U)
                            end do
                        end if
                end select

                ! Retraction on the Grassmannian
                call GM_retract_qr(X,"G")

                ! Project Transport Back to retraction
                !select case(rqm_settings%rqm_solver%solver)
                !    case(INTEGRATOR_VPO)
                !    case(INTEGRATOR_LBFGS)
                !        if(.False.)then
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_step,X)
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_force,X)
                !            do i_m=1,solversettings%lbfgs_memory
                !                call to_tangentspace_grassmann_manifold(lbfgs_d(:,i_m),X)
                !                call to_tangentspace_grassmann_manifold(lbfgs_y(:,i_m),X)
                !            end do
                !        end if
                !end select
                !---------------------------Quantities for next iteration-----------------------------------------------
                if (rqm_settings%fd_settings%i_richardson) then
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                        & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                        & rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                fd_step_used(p) = dummyval
                                HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                    end select
                else
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, &
                                        & rqm_settings%fd_settings%order, dummyvec3N, &
                                        & dummyvec3N_2)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        end select
                end if

                ! Compute Rayleigh Matrix X^T H X
                rq_matrix = 0.0d0
                do p_iter=1,p
                    do p_iter2=1,p
                        ! col idx, row idx
                        rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N_atom*(p_iter2-1)+1:2*N_atom*p_iter2))
                    end do
                end do
                ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
                rq = 0.0d0
                do p_iter=1,p
                    rq = rq + rq_matrix(p_iter,p_iter)
                end do
                ! The gradient is given by H X - X (X^T H X) and is a 3N x p Matrix
                rq_gradient = HX_findiff
                do p_iter=1,p
                    do p_iter2=1,p
                        rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) &
                                & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                    end do
                end do
                rq_gradient = rq_gradient * 2.0d0

                ! Change representation of Gradient:
                do p_iter=1,p
                    rq_gradient_2d(:,p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
                end do
                ! Compute the norm of the Gradient which is given by the Frobenius Norm:
                ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
                rq_gradient_norm = 0.0d0
                do p_iter=1,p
                    rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter),&
                            & rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter))
                end do
                ! For our solvers we need the force representation
                rq_force = -1.0d0 * rq_gradient

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        !Update the velocity based on the force for the next iteration and the rotated b
                        call vpo_update_velocity(rq_force, vpo_b, vpo_velocity, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt)
                end select
            end do

            rayleighritz_evec = rq_matrix
            call syev(rayleighritz_evec,evals,jobz="V")
            v_fin = 0.0d0
            do p_iter=1,p
                dummyvec2N = 0.0d0
                do p_iter2=1,p
                    dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                end do
                call vec_to_embeddingspace(basis, dummyvec2N, v_fin(:,p_iter))
            end do
        end subroutine rqm_grassmann_performance

        subroutine rqm_grassmann_information(interface_char,spin,basis,p,v_ini,v_fin,rqm_settings,evals,v_fin_exact, &
                & info_out_file, rqm_iterations)
            !> This routine performs an optimization routine on the Grassmann Manifold with the objective function
            !> tr X^T H X, where H is the Hessian of the spin configuration and X is a 2N x p matrix representative
            !> for an equivalence class of matrices on the Stiefel manifold corresponding to one point on the
            !> Grassmannian. The subspace spanned by the columns of X_min is the invariant subspace belonging to
            !> the p smallest eigenvalues. This routine performs the Rayleigh Ritz procedure in every iteration and
            !> compares the eigenvectors with the correct ones.
            !===================================Input Variable======================================================
            real(kind=xp), intent(in) :: spin(:)                                    !< Spin Configuration
            real(kind=xp), intent(in) :: basis(:,:,:)                               !< Basis of Tangent Space of Spin
            integer, intent(in) :: p                                                !< Dimension of the subspace
            real(kind=xp), intent(in) :: v_ini(:,:)                                 !< Start Vector for RQM (3N)
            type(rqm_params), intent(in) :: rqm_settings                            !< Parameter Container for RQM
            real(kind=xp), optional, intent(in) :: v_fin_exact(:,:)                 !< Correct eigenvectors
            character(len=*), optional, intent(in) :: info_out_file                 !< Info Out: Filename
            character, intent(in) :: interface_char                                 !< Character for interface disting.
            !================================Output Variable========================================================
            real(kind=xp), intent(out) :: v_fin(:,:)                                !< Result Vector of RQM (3N)
            real(kind=xp), intent(out) :: evals(:)                                  !< Result of the two lowest evals
            integer, optional, intent(out) :: rqm_iterations                        !< Number RQM iterations
            !================================Local Variable=========================================================
            real(kind=xp), allocatable :: X(:,:),X_previous(:,:)                    !< Points on Grassmannian
            character(len=:), allocatable :: info_out_file_local                    !< Local Copy of optional param
            logical :: i_info_out_local
            integer :: io                                                           !< Io-Stream for output file
            integer :: iter                                                         !< Iteration Counter
            integer :: p_iter, p_iter2                                              !< Iteration for Subspace Dim.
            real(kind=xp), allocatable :: dummyvec3N(:),dummyvec3N_2(:),dummyvec2N(:)!< Dummy Vectors
            real(kind=xp) :: dummyval                                               !< Dummy Value
            real(kind=xp) :: rq_matrix(p,p)                                         !< Matrix X^T * H * X
            real(kind=xp), allocatable :: HX_findiff(:)                             !< Finite difference approximation
                                                                                        !< of the matrix-vector product
            real(kind=xp), allocatable :: rq_gradient(:)                            !< Gradient of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_gradient_2d(:,:)                       !< Matrix representation
            real(kind=xp), allocatable :: rq_force(:)                               !< Force of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_step(:)                                !< Step for updating the RQ
            real(kind=xp), allocatable :: rq_step_2d(:,:)                           !< Step for updating the X quantity
            real(kind=xp) :: rq_gradient_norm                                       !< Norm of the Ralyeigh Gradient G
                                                                                    !< sqrt(tr(G^T G))
            real(kind=xp) :: rq                                                     !< Rayleigh Quotient
            real(kind=xp) :: rayleighritz_evec(p,p)                                 !< Eigenvalues and Eigenvector of
                                                                                    !< Rayleigh Ritz Procedure
            real(kind=xp) :: comp_w_exact(p)                                        !< Dot product between the current
                                                                                    !< evec and the correct eigenvec.
            real(kind=xp), allocatable :: U(:,:), S(:), VT(:,:)                     !< thin SVD of rq_step_2d
            real(kind=xp), allocatable :: T(:,:)                                    !< Parallel transport matrix
            !------------------------------FD Information---------------------------------------------------------------
            real(kind=xp) :: fd_step_used(p)                                        !< Finite Difference Step used (only
                                                                                    !< needed in case of Richardson-FD)
            !------------------------------Algorithm Depending Variables------------------------------------------------
            ! Depending on which algorithm is chosen one might need different local variable, they will be only
            ! allocated if needed. In this case all the solver quantities have a larger dimension depending on the size
            ! of the invariant subspace
            ! .................................................VPO......................................................
            real(kind=xp), allocatable :: vpo_velocity(:)                           !< Velocity in Velocity-Projection-Algo.
            real(kind=xp), allocatable :: vpo_b(:)                                  !< Local vpo quantity
            ! .................................................LBFGS....................................................
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_d
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_y
            real(kind=xp), allocatable, dimension(:) :: lbfgs_rho
            real(kind=xp), allocatable, dimension(:) :: lbfgs_gamma
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_force
            real(kind=xp), allocatable, dimension(:) :: lbfgs_current_step
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_step
            real(kind=xp) :: lbfgs_previous_steplength, lbfgs_steplength
            integer :: i_m
            ! Linesearch Variables
            real(kind=xp) :: leftside, rightside, alpha, rho, c
            real(kind=xp), allocatable :: lbfgs_Hstep(:,:),lbfgs_linesearchdummy(:,:)
            ! Timing variables
            real(kind=xp) :: fintime, initime
            !.............................Local Variable Allocation.....................................................
            ! These local variables are made allocatable, since they can reach sizes beyond stack limit (heap arrays)
            allocate(dummyvec2N(2*N_atom))                      ! Dummy 2N vector
            allocate(dummyvec3N(3*N_atom))                      ! Dummy 3N vector
            allocate(dummyvec3N_2(3*N_atom))                    ! Dummy 3N vector
            allocate(X(2*N_atom,p),X_previous(2*N_atom,p))      ! Optimization quantity
            allocate(HX_findiff(2*N_atom*p))                    ! Store the 2N x p Matrix as 1D array
            allocate(rq_gradient(2*N_atom*p))                   ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_gradient_2d(2*N_atom,p))                ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_force(2*N_atom*p))                      ! Forcematrix (2N x p) of the RQ as 1D
            allocate(rq_step(2*N_atom*p))                       ! Stepmatrix (2N x p) of the RQ as 1D
            allocate(rq_step_2d(2*N_atom,p))                    ! Stepmatrix (2N x p) of the RQ as 2D
            allocate(T(2*N_atom,p),U(2*N_atom,p),VT(p,p),S(p))  ! Transport matrix and thin SVD of rq_step_2d
            !.............................Weird Fortran Default Argument Setting........................................
            if (present(info_out_file)) then
                info_out_file_local = info_out_file
                i_info_out_local = .True.
            else
                info_out_file_local = "info_rqm.csv"
            end if
            !.............................Prepare Solver............................................................
            select case(rqm_settings%rqm_solver%solver)
                case(INTEGRATOR_VPO)
                    allocate(vpo_velocity(2*N_atom*p),vpo_b(2*N_atom*p))
                    vpo_velocity = 0.0d0
                    vpo_b = 0.0d0
                case(INTEGRATOR_LBFGS)
                    allocate(lbfgs_d(2*N_atom*p,rqm_settings%rqm_solver%lbfgs_memory),lbfgs_y(2*N_atom*p,&
                            & rqm_settings%rqm_solver%lbfgs_memory), lbfgs_rho(rqm_settings%rqm_solver%lbfgs_memory), &
                            & lbfgs_gamma(rqm_settings%rqm_solver%lbfgs_memory))
                    allocate(lbfgs_previous_force(2*N_atom*p), lbfgs_current_step(2*N_atom*p),&
                            & lbfgs_previous_step(2*N_atom*p))
                    allocate(lbfgs_Hstep(2*N_atom,p),lbfgs_linesearchdummy(2*N_atom,p))
                    lbfgs_steplength = 1.0d0
                    lbfgs_previous_force = 0.0d0
                    rho = 0.5d0
                    c = 0.0001d0
            end select
            !.............................Prepare Info File.............................................................
            if (i_info_out_local) then
                open(newunit=io,file=trim(adjustl(info_out_file_local)),action='write', form='formatted')
                write(io,'((A10,3x),4(A23,3x))',advance="no") 'iteration', 'rq', 'norm_grad', "time_r_grad"
                do p_iter=1,p
                    write(io,"(1(A23,3x))",advance="no") "eig" // trim(adjustl(int_to_str(p_iter)))
                end do
                if (present(v_fin_exact)) then
                    do p_iter=1,p
                        write(io,"(1(A23,3x))",advance="no") "comp_w_correct_" // trim(adjustl(int_to_str(p_iter)))
                    end do
                end if
                if (rqm_settings%fd_settings%i_richardson) then
                    do p_iter=1,p
                        write(io,'(1(A23,3x))',advance="no") "fd_step_" // trim(adjustl(int_to_str(p_iter)))
                    end do
                end if
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        write(io,'((A23,3x))') 'velocity_norm'
                    case(INTEGRATOR_LBFGS)
                        write(io,'((A23,3x))') 'steplength'
                    case(INTEGRATOR_EULER)
                        write(io,'((A23,3x))') 'eulerdummy'
                end select
            end if
            !.............................Iteration Zero (Preparation)..................................................
            if (present(rqm_iterations)) rqm_iterations = 0
            v_fin = 0.0d0
            rq = 0.0d0
            X = 0.0d0
            evals = 0.0d0
            ! Initialize 2N x p matrix X
            do p_iter=1,p
                call vec_to_tangentspace(basis,v_ini(:,p_iter),dummyvec2N)
                X(:,p_iter) = dummyvec2N / norm(dummyvec2N)
            end do

            call GM_retract_qr(X,"G")

            ! Compute the Hessian-vector product H * X (in R^{2Nxp}), which is needed for computed the Gradient
            ! of the Grassmanian as well as evaluating the objective function
            initime = omp_get_wtime()
            if (rqm_settings%fd_settings%i_richardson) then
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                    & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                    & rqm_settings%fd_settings%richardson_iter)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            fd_step_used(p) = dummyval
                            HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                end select
            else
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, &
                                    & rqm_settings%fd_settings%order, dummyvec3N, dummyvec3N_2)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    end select
            end if

            ! Compute Rayleigh Matrix X^T H X
            rq_matrix = 0.0d0
            do p_iter=1,p
                do p_iter2=1,p
                    ! col idx, row idx
                    rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N_atom*(p_iter2-1)+1:2*N_atom*p_iter2))
                end do
            end do
            ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
            rq = 0.0d0
            do p_iter=1,p
                rq = rq + rq_matrix(p_iter,p_iter)
            end do

            ! The gradient is given by 2 H X - 2 X (X^T H X) and is a 2N x p Matrix
            rq_gradient = HX_findiff
            do p_iter=1,p
                do p_iter2=1,p
                    rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) &
                            & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                end do
            end do
            rq_gradient = rq_gradient * 2.0d0

            ! Change representation of Gradient to a 2N x p matrix:
            do p_iter=1,p
                rq_gradient_2d(:,p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
            end do

            fintime = omp_get_wtime()
            ! Compute the norm of the Gradient which is given by the Frobenius Norm:
            ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
            rq_gradient_norm = 0.0d0
            do p_iter=1,p
                rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter),&
                        & rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter))
            end do
            rq_gradient_norm = sqrt(rq_gradient_norm)
            ! For our solvers we need the force representation
            rq_force = -1.0d0 * rq_gradient

            if (present(info_out_file)) then
                ! If output information is requested compute the Rayleigh Ritz procedure in every iteration
                rayleighritz_evec = rq_matrix
                call syev(rayleighritz_evec,evals,jobz="V")
                v_fin = 0.0d0
                do p_iter=1,p
                    dummyvec2N = 0.0d0
                    do p_iter2=1,p
                        dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                    end do
                    call vec_to_embeddingspace(basis, dummyvec2N, v_fin(:,p_iter))
                end do
                if (present(v_fin_exact)) then
                    do p_iter=1,p
                        comp_w_exact(p_iter) = DOT_PRODUCT(v_fin(:,p_iter),v_fin_exact(:,p_iter))
                    end do
                end if
            end if

            do iter=1,rqm_settings%rqm_solver%n_steps
                !-------------------------------Write outputs-----------------------------------------------------------
                if (i_info_out_local) then
                    write(io,'((I10,3x),5(E23.12E3,3x))',advance="no") iter, rq, rq_gradient_norm,fintime-initime
                    do p_iter=1,p
                        write(io,"(1(E23.12E3,3x))",advance="no") evals(p_iter)
                    end do
                    if (present(v_fin_exact)) then
                        do p_iter=1,p
                            write(io,"(1(E23.12E3,3x))",advance="no") comp_w_exact(p_iter)
                        end do
                    end if
                    if (rqm_settings%fd_settings%i_richardson) then
                        do p_iter=1,p
                            write(io,'(1(E23.12E3,3x))',advance="no") fd_step_used(p_iter)
                        end do
                    end if
                    select case(rqm_settings%rqm_solver%solver)
                        case(INTEGRATOR_VPO)
                            write(io,'(1(E23.12E3,3x))') norm(vpo_velocity)
                        case(INTEGRATOR_LBFGS)
                            write(io,'(1(E23.12E3,3x))') lbfgs_steplength
                        case(INTEGRATOR_EULER)
                            write(io,'(1(E23.12E3,3x))') 0.0
                    end select
                end if
                !----------------------------Check convergence crit-----------------------------------------------------
                if (rq_gradient_norm<=rqm_settings%rqm_solver%conv_crit) then
                    !write(*,*) 'Rayleigh optimization converged at step: ', trim(adjustl(int_to_str(iter)))
                    if (present(rqm_iterations)) rqm_iterations = iter
                    exit
                end if
                if (iter == rqm_settings%rqm_solver%n_steps) then
                    if (present(rqm_iterations)) rqm_iterations = iter
                    write(*,*) "RQM exceeded iterations"
                end if
                !----------------------------Solver---------------------------------------------------------------------
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_EULER)
                        call euler_step(rq_force,rqm_settings%rqm_solver%euler_dt,rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
                        end do
                    case(INTEGRATOR_VPO)
                        call vpo_step(vpo_velocity,rq_force, vpo_b, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt, rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
                        end do
                    case(INTEGRATOR_LBFGS)
                        ! Calculate the step and rotate the spin configuration
                        call lbfgs_step(iter, rqm_settings%rqm_solver%lbfgs_memory,  rq_force, lbfgs_previous_force, &
                                & lbfgs_previous_step, lbfgs_previous_steplength, lbfgs_rho, lbfgs_gamma, lbfgs_d, lbfgs_y, &
                                & rqm_settings%rqm_solver%lbfgs_theta_max, rq_step, lbfgs_steplength)
                        lbfgs_previous_force =  rq_force
                        lbfgs_previous_steplength = lbfgs_steplength
                        lbfgs_previous_step = rq_step
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) * lbfgs_steplength
                        end do
                        ! Simple Armijo Linesearch
                        if (.False.) then
                            alpha = 1.0d0
                            rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                            rightside = rightside *c*alpha
                            !Left side computation
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis, rq_step_2d(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                        & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                        & rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                lbfgs_Hstep(:,p_iter) = dummyvec2N
                            end do
                            lbfgs_linesearchdummy = 2.0*alpha*X+rq_step_2d*alpha**2.0
                            call GM_retract_qr(lbfgs_linesearchdummy,"G")
                            leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                            i_m = 0
                            do while(leftside>rightside)
                                i_m = i_m + 1
                                if (i_m==10) then
                                    exit
                                end if
                                write(*,*) leftside, " > " ,rightside
                                alpha = alpha * rho
                                rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                                rightside = rightside *c*alpha
                                !Left side computation
                                do p_iter=1,p
                                    call vec_to_embeddingspace(basis, rq_step_2d(:,p_iter), dummyvec3N)
                                    call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                            & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                            & rqm_settings%fd_settings%richardson_iter)
                                    call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                    lbfgs_Hstep(:,p_iter) = dummyvec2N
                                end do
                                lbfgs_linesearchdummy = 2.0*alpha*X+alpha**2*rq_step_2d
                                call GM_retract_qr(lbfgs_linesearchdummy,"G")
                                leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                                write(*,*) "alpha: ", alpha
                            end do
                            write(*,*) "finished LS:", alpha
                            lbfgs_previous_steplength = alpha
                            lbfgs_steplength = alpha
                            rq_step_2d = rq_step_2d * lbfgs_steplength
                        end if
                    end select

                ! Save the current iterate X_k
                X_previous = X
                ! Compute the next iterate X_k+1 by moving along the geodesic defined by the step and the exp. map
                call GM_retract_exp(X,rq_step_2d,1.0d0,U,S,VT,"Y")
                !X = X_previous + rq_step_2d

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                    case(INTEGRATOR_LBFGS)
                        if(.True.)then
                            call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                            call GM_paralleltransport_applyTM(lbfgs_previous_step,T,U)
                            call GM_paralleltransport_applyTM(lbfgs_previous_force,T,U)
                            do i_m=1,rqm_settings%rqm_solver%lbfgs_memory
                                call GM_paralleltransport_applyTM(lbfgs_d(:,i_m),T,U)
                                call GM_paralleltransport_applyTM(lbfgs_y(:,i_m),T,U)
                            end do
                        end if
                end select

                ! Retraction on the Grassmannian
                call GM_retract_qr(X,"G")

                ! Project Transport Back to retraction
                !select case(rqm_settings%rqm_solver%solver)
                !    case(INTEGRATOR_VPO)
                !    case(INTEGRATOR_LBFGS)
                !        if(.False.)then
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_step,X)
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_force,X)
                !            do i_m=1,solversettings%lbfgs_memory
                !                call to_tangentspace_grassmann_manifold(lbfgs_d(:,i_m),X)
                !                call to_tangentspace_grassmann_manifold(lbfgs_y(:,i_m),X)
                !            end do
                !        end if
                !end select
                !---------------------------Quantities for next iteration-----------------------------------------------
                initime = omp_get_wtime()
                if (rqm_settings%fd_settings%i_richardson) then
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N, rqm_settings%fd_settings%richardson_error, &
                                        & rqm_settings%fd_settings%step,  dummyval, dummyvec3N_2, &
                                        & rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                fd_step_used(p) = dummyval
                                HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                    end select
                else
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, &
                                        & rqm_settings%fd_settings%order, dummyvec3N, dummyvec3N_2)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                HX_findiff(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        end select
                end if

                ! Compute Rayleigh Matrix X^T H X
                rq_matrix = 0.0d0
                do p_iter=1,p
                    do p_iter2=1,p
                        ! col idx, row idx
                        rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N_atom*(p_iter2-1)+1:2*N_atom*p_iter2))
                    end do
                end do
                ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
                rq = 0.0d0
                do p_iter=1,p
                    rq = rq + rq_matrix(p_iter,p_iter)
                end do
                ! The gradient is given by H X - X (X^T H X) and is a 3N x p Matrix
                rq_gradient = HX_findiff
                do p_iter=1,p
                    do p_iter2=1,p
                        rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter) &
                                & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                    end do
                end do
                rq_gradient = rq_gradient * 2.0d0

                ! Change representation of Gradient:
                do p_iter=1,p
                    rq_gradient_2d(:,p_iter) = rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter)
                end do
                fintime = omp_get_wtime()
                ! Compute the norm of the Gradient which is given by the Frobenius Norm:
                ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
                rq_gradient_norm = 0.0d0
                do p_iter=1,p
                    rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter),&
                            & rq_gradient(2*N_atom*(p_iter-1)+1:2*N_atom*p_iter))
                end do
                ! For our solvers we need the force representation
                rq_force = -1.0d0 * rq_gradient

                if (present(info_out_file)) then
                    ! If output information is requested compute the Rayleigh Ritz procedure in every iteration
                    rayleighritz_evec = rq_matrix
                    call syev(rayleighritz_evec,evals,jobz="V")
                    v_fin = 0.0d0
                    do p_iter=1,p
                        dummyvec2N = 0.0d0
                        do p_iter2=1,p
                            dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                        end do
                        call vec_to_embeddingspace(basis, dummyvec2N, v_fin(:,p_iter))
                    end do
                    if (present(v_fin_exact)) then
                        do p_iter=1,p
                            comp_w_exact(p_iter) = DOT_PRODUCT(v_fin(:,p_iter),v_fin_exact(:,p_iter))
                        end do
                    end if
                end if

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        !Update the velocity based on the force for the next iteration and the rotated b
                        call vpo_update_velocity(rq_force, vpo_b, vpo_velocity, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt)
                end select
            end do

            rayleighritz_evec = rq_matrix
            call syev(rayleighritz_evec,evals,jobz="V")
            v_fin = 0.0d0
            do p_iter=1,p
                dummyvec2N = 0.0d0
                do p_iter2=1,p
                    dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                end do
                call vec_to_embeddingspace(basis, dummyvec2N, v_fin(:,p_iter))
            end do

            if (i_info_out_local) then
                write(io,'((I10,3x),3(E23.12E3,3x))',advance="no") iter, rq, rq_gradient_norm
                do p_iter=1,p
                    write(io,"(1(E23.12E3,3x))",advance="no") evals(p_iter)
                end do
                if (rqm_settings%fd_settings%i_richardson) then
                    do p_iter=1,p
                        write(io,'(1(E23.12E3,3x))',advance="no") fd_step_used(p_iter)
                    end do
                end if
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        write(io,'(1(E23.12E3,3x))') norm(vpo_velocity)
                    case(INTEGRATOR_LBFGS)
                        write(io,'(1(E23.12E3,3x))') lbfgs_steplength
                    case(INTEGRATOR_EULER)
                        write(io,'(1(E23.12E3,3x))') 0.0
                end select
            end if

            if (i_info_out_local) close(io)
        end subroutine rqm_grassmann_information

        subroutine rqm_grassmann_subsystem_performance(spin,basis,subsys_indices,p,v_ini,v_fin,rqm_settings,evals,rqm_iterations)
            !> This routine performs an optimization routine on the Grassmann Manifold with the objective function
            !> tr X^T H X, where H is the Hessian of the spin configuration and X is a 2N x p matrix representative
            !> for an equivalence class of matrices on the Stiefel manifold corresponding to one point on the
            !> Grassmannian. The subspace spanned by the columns of X_min is the invariant subspace belonging to
            !> the p smallest eigenvalues. This version of the routine only computes the rayleigh ritz procedure once
            !> at the end. N is here the number of spins in the subsystem
            !===================================Input Variable======================================================
            real(kind=xp), intent(in) :: spin(:)                                    !< Spin Configuration
            real(kind=xp), intent(in) :: basis(:,:,:)                               !< Basis of Tangent Space of Spin
            integer, intent(in) :: subsys_indices(:)                                !< Indices of atoms in subsystem
            integer, intent(in) :: p                                                !< Dimension of the subspace
            real(kind=xp), intent(in) :: v_ini(:,:)                                 !< Start Vector for RQM (3N)
            type(rqm_params), intent(in) :: rqm_settings                            !< Parameter Container for RQM
            !================================Output Variable========================================================
            real(kind=xp), intent(out) :: v_fin(:,:)                                !< Result Vector of RQM (3N)
            real(kind=xp), intent(out) :: evals(:)                                  !< Result of the two lowest evals
            integer, optional, intent(out) :: rqm_iterations                        !< Number RQM iterations
            !================================Local Variable=========================================================
            integer :: N                                                            !< Number of atoms
            real(kind=xp), allocatable :: X(:,:),X_previous(:,:)                    !< Points on Grassmannian
            integer :: iter                                                         !< Iteration Counter
            integer :: p_iter, p_iter2                                              !< Iteration for Subspace Dim.
            real(kind=xp), allocatable :: dummyvec3N(:),dummyvec3N_2(:),dummyvec2N(:)!< Dummy Vectors
            real(kind=xp) :: dummyval                                               !< Dummy Value
            real(kind=xp) :: rq_matrix(p,p)                                         !< Matrix X^T * H * X
            real(kind=xp), allocatable :: HX_findiff(:)                             !< Finite difference approximation
                                                                                        !< of the matrix-vector product
            real(kind=xp), allocatable :: rq_gradient(:)                            !< Gradient of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_gradient_2d(:,:)                       !< Matrix representation
            real(kind=xp), allocatable :: rq_force(:)                               !< Force of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_step(:)                                !< Step for updating the RQ
            real(kind=xp), allocatable :: rq_step_2d(:,:)                           !< Step for updating the X quantity
            real(kind=xp) :: rq_gradient_norm                                       !< Norm of the Ralyeigh Gradient G
                                                                                    !< sqrt(tr(G^T G))
            real(kind=xp) :: rq
            real(kind=xp) :: rayleighritz_evec(p,p)                                 !< Eigenvalues and Eigenvector of
                                                                                        !< Rayleigh Ritz Procedure
            real(kind=xp), allocatable :: U(:,:), S(:), VT(:,:)                     !< thin SVD of rq_step_2d
            real(kind=xp), allocatable :: T(:,:)                                    !< Parallel transport matrix
            !------------------------------FD Information---------------------------------------------------------------
            real(kind=xp) :: fd_step_used(p)                                        !< Finite Difference Step used (only
                                                                                    !< needed in case of Richardson-FD)
            !------------------------------Algorithm Depending Variables------------------------------------------------
            ! Depending on which algorithm is chosen one might need different local variable, they will be only
            ! allocated if needed. In this case all the solver quantities have a larger dimension depending on the size
            ! of the invariant subspace
            ! .................................................VPO......................................................
            real(kind=xp), allocatable :: vpo_velocity(:)                           !< Velocity in Velocity-Projection-Algo.
            real(kind=xp), allocatable :: vpo_b(:)                                  !< Local vpo quantity
            ! .................................................LBFGS....................................................
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_d
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_y
            real(kind=xp), allocatable, dimension(:) :: lbfgs_rho
            real(kind=xp), allocatable, dimension(:) :: lbfgs_gamma
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_force
            real(kind=xp), allocatable, dimension(:) :: lbfgs_current_step
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_step
            real(kind=xp) :: lbfgs_previous_steplength, lbfgs_steplength
            integer :: i_m
            ! Linesearch Variables
            real(kind=xp) :: leftside, rightside, alpha, rho, c
            real(kind=xp), allocatable :: lbfgs_Hstep(:,:),lbfgs_linesearchdummy(:,:)
            !.............................Local Variable Allocation.....................................................
            N = size(subsys_indices)
            ! These local variables are made allocatable, since they can reach sizes beyond stack limit (heap arrays)
            allocate(dummyvec2N(2*N))                       ! Dummy 2N vector
            allocate(dummyvec3N(3*N))                       ! Dummy 3N vector
            allocate(dummyvec3N_2(3*N))                     ! Dummy 3N vector
            allocate(X(2*N,p),X_previous(2*N,p))            ! Optimization quantity
            allocate(HX_findiff(2*N*p))                     ! Store the 2N x p Matrix as 1D array
            allocate(rq_gradient(2*N*p))                    ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_gradient_2d(2*N,p))                 ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_force(2*N*p))                       ! Forcematrix (2N x p) of the RQ as 1D
            allocate(rq_step(2*N*p))                        ! Stepmatrix (2N x p) of the RQ as 1D
            allocate(rq_step_2d(2*N,p))                     ! Stepmatrix (2N x p) of the RQ as 2D
            allocate(T(2*N,p),U(2*N,p),VT(p,p),S(p))        ! Transport matrix and thin SVD of rq_step_2d
            !.............................Prepare Solver............................................................
            select case(rqm_settings%rqm_solver%solver)
                case(INTEGRATOR_VPO)
                    allocate(vpo_velocity(2*N*p),vpo_b(2*N*p))
                    vpo_velocity = 0.0d0
                    vpo_b = 0.0d0
                case(INTEGRATOR_LBFGS)
                    allocate(lbfgs_d(2*N*p,rqm_settings%rqm_solver%lbfgs_memory),lbfgs_y(2*N*p,&
                            & rqm_settings%rqm_solver%lbfgs_memory), lbfgs_rho(rqm_settings%rqm_solver%lbfgs_memory), &
                            & lbfgs_gamma(rqm_settings%rqm_solver%lbfgs_memory))
                    allocate(lbfgs_previous_force(2*N*p), lbfgs_current_step(2*N*p),lbfgs_previous_step(2*N*p))
                    allocate(lbfgs_Hstep(2*N,p),lbfgs_linesearchdummy(2*N,p))
                    lbfgs_steplength = 1.0d0
                    lbfgs_previous_force = 0.0d0
                    rho = 0.5d0
                    c = 0.0001d0
            end select

            !.............................Iteration Zero (Preparation)..................................................
            if (present(rqm_iterations)) rqm_iterations = 0
            v_fin = 0.0d0
            rq = 0.0d0
            X = 0.0d0
            evals = 0.0d0
            ! Initialize 2N x p matrix X
            do p_iter=1,p
                call vec_to_tangentspace(basis,v_ini(:,p_iter),dummyvec2N)
                X(:,p_iter) = dummyvec2N / norm(dummyvec2N)
            end do

            call GM_retract_qr(X,"G")

            ! Compute the Hessian-vector product H * X (in R^{2Nxp}), which is needed for computed the Gradient
            ! of the Grassmanian as well as evaluating the objective function
            if (rqm_settings%fd_settings%i_richardson) then
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis,N,X(:,p_iter),dummyvec3N)
                            call fd_hessvec_forward(spin,dummyvec3N,subsys_indices, &
                                    & rqm_settings%fd_settings%richardson_error,rqm_settings%fd_settings%step,dummyval,&
                                    & dummyvec3N_2,rqm_settings%fd_settings%richardson_iter)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            fd_step_used(p) = dummyval
                            HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                end select
            else
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, rqm_settings%fd_settings%step,subsys_indices, &
                                    & rqm_settings%fd_settings%order, dummyvec3N, dummyvec3N_2)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    end select
            end if


            ! Compute Rayleigh Matrix X^T H X
            rq_matrix = 0.0d0
            do p_iter=1,p
                do p_iter2=1,p
                    ! col idx, row idx
                    rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N*(p_iter2-1)+1:2*N*p_iter2))
                end do
            end do
            ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
            rq = 0.0d0
            do p_iter=1,p
                rq = rq + rq_matrix(p_iter,p_iter)
            end do

            ! The gradient is given by 2 H X - 2 X (X^T H X) and is a 2N x p Matrix
            rq_gradient = HX_findiff
            do p_iter=1,p
                do p_iter2=1,p
                    rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) &
                            & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                end do
            end do
            rq_gradient = rq_gradient * 2.0d0

            ! Change representation of Gradient to a 2N x p matrix:
            do p_iter=1,p
                rq_gradient_2d(:,p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter)
            end do

            ! Compute the norm of the Gradient which is given by the Frobenius Norm:
            ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
            rq_gradient_norm = 0.0d0
            do p_iter=1,p
                rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter),&
                        & rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter))
            end do
            rq_gradient_norm = sqrt(rq_gradient_norm)
            ! For our solvers we need the force representation
            rq_force = -1.0d0 * rq_gradient

            do iter=1,rqm_settings%rqm_solver%n_steps
                !----------------------------Check convergence crit-----------------------------------------------------
                if (rq_gradient_norm<=rqm_settings%rqm_solver%conv_crit) then
                    !write(*,*) 'Rayleigh optimization converged at step: ', trim(adjustl(int_to_str(iter)))
                    if (present(rqm_iterations)) rqm_iterations = iter
                    exit
                end if
                if (iter == rqm_settings%rqm_solver%n_steps) then
                    if (present(rqm_iterations)) rqm_iterations = iter
                    write(*,*) "RQM exceeded iterations"
                end if
                !----------------------------Solver---------------------------------------------------------------------
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_EULER)
                        call euler_step(rq_force,rqm_settings%rqm_solver%euler_dt,rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N*(p_iter-1)+1:2*N*p_iter)
                        end do
                    case(INTEGRATOR_VPO)
                        call vpo_step(vpo_velocity,rq_force, vpo_b, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt, rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N*(p_iter-1)+1:2*N*p_iter)
                        end do
                    case(INTEGRATOR_LBFGS)
                        ! Calculate the step and rotate the spin configuration
                        call lbfgs_step(iter, rqm_settings%rqm_solver%lbfgs_memory,  rq_force, lbfgs_previous_force, &
                                & lbfgs_previous_step, lbfgs_previous_steplength, lbfgs_rho, lbfgs_gamma, lbfgs_d, lbfgs_y, &
                                & rqm_settings%rqm_solver%lbfgs_theta_max, rq_step, lbfgs_steplength)
                        lbfgs_previous_force =  rq_force
                        lbfgs_previous_steplength = lbfgs_steplength
                        lbfgs_previous_step = rq_step
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N*(p_iter-1)+1:2*N*p_iter) * lbfgs_steplength
                        end do
                        ! Simple Armijo Linesearch
                        if (.False.) then
                            alpha = 1.0d0
                            rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                            rightside = rightside *c*alpha
                            !Left side computation
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis,N, rq_step_2d(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N,subsys_indices,&
                                        & rqm_settings%fd_settings%richardson_error, rqm_settings%fd_settings%step,  &
                                        & dummyval, dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                lbfgs_Hstep(:,p_iter) = dummyvec2N
                            end do
                            lbfgs_linesearchdummy = 2.0*alpha*X+rq_step_2d*alpha**2.0
                            call GM_retract_qr(lbfgs_linesearchdummy,"G")
                            leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                            i_m = 0
                            do while(leftside>rightside)
                                i_m = i_m + 1
                                if (i_m==10) then
                                    exit
                                end if
                                write(*,*) leftside, " > " ,rightside
                                alpha = alpha * rho
                                rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                                rightside = rightside *c*alpha
                                !Left side computation
                                do p_iter=1,p
                                    call vec_to_embeddingspace(basis,N, rq_step_2d(:,p_iter), dummyvec3N)
                                    call fd_hessvec_forward(spin, dummyvec3N,subsys_indices, &
                                            & rqm_settings%fd_settings%richardson_error,rqm_settings%fd_settings%step,  &
                                            & dummyval, dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                                    call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                    lbfgs_Hstep(:,p_iter) = dummyvec2N
                                end do
                                lbfgs_linesearchdummy = 2.0*alpha*X+alpha**2*rq_step_2d
                                call GM_retract_qr(lbfgs_linesearchdummy,"G")
                                leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                                write(*,*) "alpha: ", alpha
                            end do
                            write(*,*) "finished LS:", alpha
                            lbfgs_previous_steplength = alpha
                            lbfgs_steplength = alpha
                            rq_step_2d = rq_step_2d * lbfgs_steplength
                        end if
                    end select

                ! Save the current iterate X_k
                X_previous = X
                ! Compute the next iterate X_k+1 by moving along the geodesic defined by the step and the exp. map
                call GM_retract_exp(X,rq_step_2d,1.0d0,U,S,VT,"Y")
                !X = X_previous + rq_step_2d

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                    case(INTEGRATOR_LBFGS)
                        if(.True.)then
                            call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                            call GM_paralleltransport_applyTM(lbfgs_previous_step,T,U)
                            call GM_paralleltransport_applyTM(lbfgs_previous_force,T,U)
                            do i_m=1,rqm_settings%rqm_solver%lbfgs_memory
                                call GM_paralleltransport_applyTM(lbfgs_d(:,i_m),T,U)
                                call GM_paralleltransport_applyTM(lbfgs_y(:,i_m),T,U)
                            end do
                        end if
                end select

                ! Retraction on the Grassmannian
                call GM_retract_qr(X,"G")

                ! Project Transport Back to retraction
                !select case(rqm_settings%rqm_solver%solver)
                !    case(INTEGRATOR_VPO)
                !    case(INTEGRATOR_LBFGS)
                !        if(.False.)then
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_step,X)
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_force,X)
                !            do i_m=1,solversettings%lbfgs_memory
                !                call to_tangentspace_grassmann_manifold(lbfgs_d(:,i_m),X)
                !                call to_tangentspace_grassmann_manifold(lbfgs_y(:,i_m),X)
                !            end do
                !        end if
                !end select
                !---------------------------Quantities for next iteration-----------------------------------------------
                if (rqm_settings%fd_settings%i_richardson) then
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N,subsys_indices, &
                                        & rqm_settings%fd_settings%richardson_error, rqm_settings%fd_settings%step,  &
                                        & dummyval, dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                fd_step_used(p) = dummyval
                                HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                    end select
                else
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, subsys_indices, &
                                        & rqm_settings%fd_settings%order, dummyvec3N, &
                                        & dummyvec3N_2)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        end select
                end if

                ! Compute Rayleigh Matrix X^T H X
                rq_matrix = 0.0d0
                do p_iter=1,p
                    do p_iter2=1,p
                        ! col idx, row idx
                        rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N*(p_iter2-1)+1:2*N*p_iter2))
                    end do
                end do
                ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
                rq = 0.0d0
                do p_iter=1,p
                    rq = rq + rq_matrix(p_iter,p_iter)
                end do
                ! The gradient is given by H X - X (X^T H X) and is a 3N x p Matrix
                rq_gradient = HX_findiff
                do p_iter=1,p
                    do p_iter2=1,p
                        rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) &
                                & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                    end do
                end do
                rq_gradient = rq_gradient * 2.0d0

                ! Change representation of Gradient:
                do p_iter=1,p
                    rq_gradient_2d(:,p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter)
                end do
                ! Compute the norm of the Gradient which is given by the Frobenius Norm:
                ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
                rq_gradient_norm = 0.0d0
                do p_iter=1,p
                    rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter),&
                            & rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter))
                end do
                ! For our solvers we need the force representation
                rq_force = -1.0d0 * rq_gradient

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        !Update the velocity based on the force for the next iteration and the rotated b
                        call vpo_update_velocity(rq_force, vpo_b, vpo_velocity, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt)
                end select
            end do

            rayleighritz_evec = rq_matrix
            call syev(rayleighritz_evec,evals,jobz="V")
            v_fin = 0.0d0
            do p_iter=1,p
                dummyvec2N = 0.0d0
                do p_iter2=1,p
                    dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                end do
                call vec_to_embeddingspace(basis,N, dummyvec2N, v_fin(:,p_iter))
            end do
        end subroutine rqm_grassmann_subsystem_performance

        subroutine rqm_grassmann_subsystem_information(interface_char,spin,basis,subsys_indices,p,v_ini,v_fin, &
                & rqm_settings,evals,v_fin_exact, info_out_file, rqm_iterations)
            !> This routine performs an optimization routine on the Grassmann Manifold with the objective function
            !> tr X^T H X, where H is the Hessian of the spin configuration and X is a 2N x p matrix representative
            !> for an equivalence class of matrices on the Stiefel manifold corresponding to one point on the
            !> Grassmannian. The subspace spanned by the columns of X_min is the invariant subspace belonging to
            !> the p smallest eigenvalues. This routine performs the Rayleigh Ritz procedure in every iteration and
            !> compares the eigenvectors with the correct ones.
            !===================================Input Variable======================================================
            real(kind=xp), intent(in) :: spin(:)                                    !< Spin Configuration
            real(kind=xp), intent(in) :: basis(:,:,:)                               !< Basis of Tangent Space of Spin
            integer, intent(in) :: subsys_indices(:)                                !< Indices Atoms sub-system
            integer, intent(in) :: p                                                !< Dimension of the subspace
            real(kind=xp), intent(in) :: v_ini(:,:)                                 !< Start Vector for RQM (3N)
            type(rqm_params), intent(in) :: rqm_settings                            !< Parameter Container for RQM
            real(kind=xp), optional, intent(in) :: v_fin_exact(:,:)                 !< Correct eigenvectors
            character(len=*), optional, intent(in) :: info_out_file                 !< Info Out: Filename
            character, intent(in) :: interface_char                                 !< Character for interface disting.
            !================================Output Variable========================================================
            real(kind=xp), intent(out) :: v_fin(:,:)                                !< Result Vector of RQM (3N)
            real(kind=xp), intent(out) :: evals(:)                                  !< Result of the two lowest evals
            integer, optional, intent(out) :: rqm_iterations                        !< Number RQM iterations
            !================================Local Variable=========================================================
            integer :: N                                                            !< Number atoms sub-system
            real(kind=xp), allocatable :: X(:,:),X_previous(:,:)                    !< Points on Grassmannian
            character(len=:), allocatable :: info_out_file_local                    !< Local Copy of optional param
            logical :: i_info_out_local
            integer :: io                                                           !< Io-Stream for output file
            integer :: iter                                                         !< Iteration Counter
            integer :: p_iter, p_iter2                                              !< Iteration for Subspace Dim.
            real(kind=xp), allocatable :: dummyvec3N(:),dummyvec3N_2(:),dummyvec2N(:)!< Dummy Vectors
            real(kind=xp) :: dummyval                                               !< Dummy Value
            real(kind=xp) :: rq_matrix(p,p)                                         !< Matrix X^T * H * X
            real(kind=xp), allocatable :: HX_findiff(:)                             !< Finite difference approximation
                                                                                        !< of the matrix-vector product
            real(kind=xp), allocatable :: rq_gradient(:)                            !< Gradient of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_gradient_2d(:,:)                       !< Matrix representation
            real(kind=xp), allocatable :: rq_force(:)                               !< Force of Rayleigh Quotient
            real(kind=xp), allocatable :: rq_step(:)                                !< Step for updating the RQ
            real(kind=xp), allocatable :: rq_step_2d(:,:)                           !< Step for updating the X quantity
            real(kind=xp) :: rq_gradient_norm                                       !< Norm of the Ralyeigh Gradient G
                                                                                    !< sqrt(tr(G^T G))
            real(kind=xp) :: rq                                                     !< Rayleigh Quotient
            real(kind=xp) :: rayleighritz_evec(p,p)                                 !< Eigenvalues and Eigenvector of
                                                                                    !< Rayleigh Ritz Procedure
            real(kind=xp) :: comp_w_exact(p)                                        !< Dot product between the current
                                                                                    !< evec and the correct eigenvec.
            real(kind=xp), allocatable :: U(:,:), S(:), VT(:,:)                     !< thin SVD of rq_step_2d
            real(kind=xp), allocatable :: T(:,:)                                    !< Parallel transport matrix
            !------------------------------FD Information---------------------------------------------------------------
            real(kind=xp) :: fd_step_used(p)                                        !< Finite Difference Step used (only
                                                                                    !< needed in case of Richardson-FD)
            !------------------------------Algorithm Depending Variables------------------------------------------------
            ! Depending on which algorithm is chosen one might need different local variable, they will be only
            ! allocated if needed. In this case all the solver quantities have a larger dimension depending on the size
            ! of the invariant subspace
            ! .................................................VPO......................................................
            real(kind=xp), allocatable :: vpo_velocity(:)                           !< Velocity in Velocity-Projection-Algo.
            real(kind=xp), allocatable :: vpo_b(:)                                  !< Local vpo quantity
            ! .................................................LBFGS....................................................
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_d
            real(kind=xp), allocatable, dimension(:,:) :: lbfgs_y
            real(kind=xp), allocatable, dimension(:) :: lbfgs_rho
            real(kind=xp), allocatable, dimension(:) :: lbfgs_gamma
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_force
            real(kind=xp), allocatable, dimension(:) :: lbfgs_current_step
            real(kind=xp), allocatable, dimension(:) :: lbfgs_previous_step
            real(kind=xp) :: lbfgs_previous_steplength, lbfgs_steplength
            integer :: i_m
            ! Linesearch Variables
            real(kind=xp) :: leftside, rightside, alpha, rho, c
            real(kind=xp), allocatable :: lbfgs_Hstep(:,:),lbfgs_linesearchdummy(:,:)
            !.............................Local Variable Allocation.....................................................
            N=size(subsys_indices)
            ! These local variables are made allocatable, since they can reach sizes beyond stack limit (heap arrays)
            allocate(dummyvec2N(2*N))                       ! Dummy 2N vector
            allocate(dummyvec3N(3*N))                       ! Dummy 3N vector
            allocate(dummyvec3N_2(3*N))                     ! Dummy 3N vector
            allocate(X(2*N,p),X_previous(2*N,p))            ! Optimization quantity
            allocate(HX_findiff(2*N*p))                     ! Store the 2N x p Matrix as 1D array
            allocate(rq_gradient(2*N*p))                    ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_gradient_2d(2*N,p))                 ! Gradientmatrix (2N x p) of the RQ as 1D
            allocate(rq_force(2*N*p))                       ! Forcematrix (2N x p) of the RQ as 1D
            allocate(rq_step(2*N*p))                        ! Stepmatrix (2N x p) of the RQ as 1D
            allocate(rq_step_2d(2*N,p))                     ! Stepmatrix (2N x p) of the RQ as 2D
            allocate(T(2*N,p),U(2*N,p),VT(p,p),S(p))  !     Transport matrix and thin SVD of rq_step_2d
            !.............................Weird Fortran Default Argument Setting........................................
            if (present(info_out_file)) then
                info_out_file_local = info_out_file
                i_info_out_local = .True.
            else
                info_out_file_local = "info_rqm.csv"
            end if
            !.............................Prepare Solver............................................................
            select case(rqm_settings%rqm_solver%solver)
                case(INTEGRATOR_VPO)
                    allocate(vpo_velocity(2*N*p),vpo_b(2*N*p))
                    vpo_velocity = 0.0d0
                    vpo_b = 0.0d0
                case(INTEGRATOR_LBFGS)
                    allocate(lbfgs_d(2*N*p,rqm_settings%rqm_solver%lbfgs_memory),lbfgs_y(2*N*p,&
                            & rqm_settings%rqm_solver%lbfgs_memory), lbfgs_rho(rqm_settings%rqm_solver%lbfgs_memory), &
                            & lbfgs_gamma(rqm_settings%rqm_solver%lbfgs_memory))
                    allocate(lbfgs_previous_force(2*N*p), lbfgs_current_step(2*N*p),&
                            & lbfgs_previous_step(2*N*p))
                    allocate(lbfgs_Hstep(2*N,p),lbfgs_linesearchdummy(2*N,p))
                    lbfgs_steplength = 1.0d0
                    lbfgs_previous_force = 0.0d0
                    rho = 0.5d0
                    c = 0.0001d0
            end select
            !.............................Prepare Info File.............................................................
            if (i_info_out_local) then
                open(newunit=io,file=trim(adjustl(info_out_file_local)),action='write', form='formatted')
                write(io,'((A10,3x),3(A23,3x))',advance="no") 'iteration', 'rq', 'norm_grad'
                do p_iter=1,p
                    write(io,"(1(A23,3x))",advance="no") "eig" // trim(adjustl(int_to_str(p_iter)))
                end do
                if (present(v_fin_exact)) then
                    do p_iter=1,p
                        write(io,"(1(A23,3x))",advance="no") "comp_w_correct_" // trim(adjustl(int_to_str(p_iter)))
                    end do
                end if
                if (rqm_settings%fd_settings%i_richardson) then
                    do p_iter=1,p
                        write(io,'(1(A23,3x))',advance="no") "fd_step_" // trim(adjustl(int_to_str(p_iter)))
                    end do
                end if
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        write(io,'((A23,3x))') 'velocity_norm'
                    case(INTEGRATOR_LBFGS)
                        write(io,'((A23,3x))') 'steplength'
                    case(INTEGRATOR_EULER)
                        write(io,'((A23,3x))') 'eulerdummy'
                end select
            end if
            !.............................Iteration Zero (Preparation)..................................................
            if (present(rqm_iterations)) rqm_iterations = 0
            v_fin = 0.0d0
            rq = 0.0d0
            X = 0.0d0
            evals = 0.0d0
            ! Initialize 2N x p matrix X
            do p_iter=1,p
                call vec_to_tangentspace(basis,v_ini(:,p_iter),dummyvec2N)
                X(:,p_iter) = dummyvec2N / norm(dummyvec2N)
            end do

            call GM_retract_qr(X,"G")

            ! Compute the Hessian-vector product H * X (in R^{2Nxp}), which is needed for computed the Gradient
            ! of the Grassmanian as well as evaluating the objective function
            if (rqm_settings%fd_settings%i_richardson) then
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, dummyvec3N,subsys_indices, &
                                    & rqm_settings%fd_settings%richardson_error,rqm_settings%fd_settings%step,dummyval,&
                                    & dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            fd_step_used(p) = dummyval
                            HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                end select
            else
                select case(rqm_settings%fd_settings%scheme)
                    case(FD_SCHEME_FORWARD)
                        do p_iter=1,p
                            call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                            call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, subsys_indices, &
                                    & rqm_settings%fd_settings%order, dummyvec3N, dummyvec3N_2)
                            call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                            HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                        end do
                    case(FD_SCHEME_BACKWARD)
                        call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    case(FD_SCHEME_CENTRAL)
                        call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                        stop
                    end select
            end if

            ! Compute Rayleigh Matrix X^T H X
            rq_matrix = 0.0d0
            do p_iter=1,p
                do p_iter2=1,p
                    ! col idx, row idx
                    rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N*(p_iter2-1)+1:2*N*p_iter2))
                end do
            end do
            ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
            rq = 0.0d0
            do p_iter=1,p
                rq = rq + rq_matrix(p_iter,p_iter)
            end do

            ! The gradient is given by 2 H X - 2 X (X^T H X) and is a 2N x p Matrix
            rq_gradient = HX_findiff
            do p_iter=1,p
                do p_iter2=1,p
                    rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) &
                            & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                end do
            end do
            rq_gradient = rq_gradient * 2.0d0

            ! Change representation of Gradient to a 2N x p matrix:
            do p_iter=1,p
                rq_gradient_2d(:,p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter)
            end do

            ! Compute the norm of the Gradient which is given by the Frobenius Norm:
            ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
            rq_gradient_norm = 0.0d0
            do p_iter=1,p
                rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter),&
                        & rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter))
            end do
            rq_gradient_norm = sqrt(rq_gradient_norm)
            ! For our solvers we need the force representation
            rq_force = -1.0d0 * rq_gradient

            if (present(info_out_file)) then
                ! If output information is requested compute the Rayleigh Ritz procedure in every iteration
                rayleighritz_evec = rq_matrix
                call syev(rayleighritz_evec,evals,jobz="V")
                v_fin = 0.0d0
                do p_iter=1,p
                    dummyvec2N = 0.0d0
                    do p_iter2=1,p
                        dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                    end do
                    call vec_to_embeddingspace(basis, dummyvec2N, v_fin(:,p_iter))
                end do
                if (present(v_fin_exact)) then
                    do p_iter=1,p
                        comp_w_exact(p_iter) = DOT_PRODUCT(v_fin(:,p_iter),v_fin_exact(:,p_iter))
                    end do
                end if
            end if

            do iter=1,rqm_settings%rqm_solver%n_steps
                !-------------------------------Write outputs-----------------------------------------------------------
                if (i_info_out_local) then
                    write(io,'((I10,3x),3(E23.12E3,3x))',advance="no") iter, rq, rq_gradient_norm
                    do p_iter=1,p
                        write(io,"(1(E23.12E3,3x))",advance="no") evals(p_iter)
                    end do
                    if (present(v_fin_exact)) then
                        do p_iter=1,p
                            write(io,"(1(E23.12E3,3x))",advance="no") comp_w_exact(p_iter)
                        end do
                    end if
                    if (rqm_settings%fd_settings%i_richardson) then
                        do p_iter=1,p
                            write(io,'(1(E23.12E3,3x))',advance="no") fd_step_used(p_iter)
                        end do
                    end if
                    select case(rqm_settings%rqm_solver%solver)
                        case(INTEGRATOR_VPO)
                            write(io,'(1(E23.12E3,3x))') norm(vpo_velocity)
                        case(INTEGRATOR_LBFGS)
                            write(io,'(1(E23.12E3,3x))') lbfgs_steplength
                        case(INTEGRATOR_EULER)
                            write(io,'(1(E23.12E3,3x))') 0.0
                    end select
                end if
                !----------------------------Check convergence crit-----------------------------------------------------
                if (rq_gradient_norm<=rqm_settings%rqm_solver%conv_crit) then
                    !write(*,*) 'Rayleigh optimization converged at step: ', trim(adjustl(int_to_str(iter)))
                    if (present(rqm_iterations)) rqm_iterations = iter
                    exit
                end if
                if (iter == rqm_settings%rqm_solver%n_steps) then
                    if (present(rqm_iterations)) rqm_iterations = iter
                    write(*,*) "RQM exceeded iterations"
                end if
                !----------------------------Solver---------------------------------------------------------------------
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_EULER)
                        call euler_step(rq_force,rqm_settings%rqm_solver%euler_dt,rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N*(p_iter-1)+1:2*N*p_iter)
                        end do
                    case(INTEGRATOR_VPO)
                        call vpo_step(vpo_velocity,rq_force, vpo_b, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt, rq_step)
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N*(p_iter-1)+1:2*N*p_iter)
                        end do
                    case(INTEGRATOR_LBFGS)
                        ! Calculate the step and rotate the spin configuration
                        call lbfgs_step(iter, rqm_settings%rqm_solver%lbfgs_memory,  rq_force, lbfgs_previous_force, &
                                & lbfgs_previous_step, lbfgs_previous_steplength, lbfgs_rho, lbfgs_gamma, lbfgs_d, lbfgs_y, &
                                & rqm_settings%rqm_solver%lbfgs_theta_max, rq_step, lbfgs_steplength)
                        lbfgs_previous_force =  rq_force
                        lbfgs_previous_steplength = lbfgs_steplength
                        lbfgs_previous_step = rq_step
                        ! Change representation of Step:
                        do p_iter=1,p
                            rq_step_2d(:,p_iter) = rq_step(2*N*(p_iter-1)+1:2*N*p_iter) * lbfgs_steplength
                        end do
                        ! Simple Armijo Linesearch
                        if (.False.) then
                            alpha = 1.0d0
                            rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                            rightside = rightside *c*alpha
                            !Left side computation
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis,N, rq_step_2d(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N,subsys_indices, &
                                        & rqm_settings%fd_settings%richardson_error,rqm_settings%fd_settings%step,  &
                                        & dummyval, dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                lbfgs_Hstep(:,p_iter) = dummyvec2N
                            end do
                            lbfgs_linesearchdummy = 2.0*alpha*X+rq_step_2d*alpha**2.0
                            call GM_retract_qr(lbfgs_linesearchdummy,"G")
                            leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                            i_m = 0
                            do while(leftside>rightside)
                                i_m = i_m + 1
                                if (i_m==10) then
                                    exit
                                end if
                                write(*,*) leftside, " > " ,rightside
                                alpha = alpha * rho
                                rightside = innerproduct_frobenius(rq_step_2d,rq_gradient_2d)
                                rightside = rightside *c*alpha
                                !Left side computation
                                do p_iter=1,p
                                    call vec_to_embeddingspace(basis,N, rq_step_2d(:,p_iter), dummyvec3N)
                                    call fd_hessvec_forward(spin, dummyvec3N, subsys_indices, &
                                            & rqm_settings%fd_settings%richardson_error, rqm_settings%fd_settings%step,&
                                            & dummyval, dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                                    call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                    lbfgs_Hstep(:,p_iter) = dummyvec2N
                                end do
                                lbfgs_linesearchdummy = 2.0*alpha*X+alpha**2*rq_step_2d
                                call GM_retract_qr(lbfgs_linesearchdummy,"G")
                                leftside = innerproduct_frobenius(lbfgs_linesearchdummy,lbfgs_Hstep)
                                write(*,*) "alpha: ", alpha
                            end do
                            write(*,*) "finished LS:", alpha
                            lbfgs_previous_steplength = alpha
                            lbfgs_steplength = alpha
                            rq_step_2d = rq_step_2d * lbfgs_steplength
                        end if
                    end select

                ! Save the current iterate X_k
                X_previous = X
                ! Compute the next iterate X_k+1 by moving along the geodesic defined by the step and the exp. map
                call GM_retract_exp(X,rq_step_2d,1.0d0,U,S,VT,"Y")
                !X = X_previous + rq_step_2d

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                    case(INTEGRATOR_LBFGS)
                        if(.True.)then
                            call GM_paralleltransport_calcTM(X_previous,1.0d0,U,S,VT,T)
                            call GM_paralleltransport_applyTM(lbfgs_previous_step,T,U)
                            call GM_paralleltransport_applyTM(lbfgs_previous_force,T,U)
                            do i_m=1,rqm_settings%rqm_solver%lbfgs_memory
                                call GM_paralleltransport_applyTM(lbfgs_d(:,i_m),T,U)
                                call GM_paralleltransport_applyTM(lbfgs_y(:,i_m),T,U)
                            end do
                        end if
                end select

                ! Retraction on the Grassmannian
                call GM_retract_qr(X,"G")

                ! Project Transport Back to retraction
                !select case(rqm_settings%rqm_solver%solver)
                !    case(INTEGRATOR_VPO)
                !    case(INTEGRATOR_LBFGS)
                !        if(.False.)then
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_step,X)
                !            call to_tangentspace_grassmann_manifold(lbfgs_previous_force,X)
                !            do i_m=1,solversettings%lbfgs_memory
                !                call to_tangentspace_grassmann_manifold(lbfgs_d(:,i_m),X)
                !                call to_tangentspace_grassmann_manifold(lbfgs_y(:,i_m),X)
                !            end do
                !        end if
                !end select
                !---------------------------Quantities for next iteration-----------------------------------------------
                if (rqm_settings%fd_settings%i_richardson) then
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, dummyvec3N,subsys_indices, &
                                        & rqm_settings%fd_settings%richardson_error, rqm_settings%fd_settings%step, &
                                        & dummyval, dummyvec3N_2, rqm_settings%fd_settings%richardson_iter)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                fd_step_used(p) = dummyval
                                HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                    end select
                else
                    select case(rqm_settings%fd_settings%scheme)
                        case(FD_SCHEME_FORWARD)
                            do p_iter=1,p
                                call vec_to_embeddingspace(basis,N, X(:,p_iter), dummyvec3N)
                                call fd_hessvec_forward(spin, rqm_settings%fd_settings%step, subsys_indices, &
                                        & rqm_settings%fd_settings%order, dummyvec3N, dummyvec3N_2)
                                call vec_to_tangentspace(basis,dummyvec3N_2,dummyvec2N)
                                HX_findiff(2*N*(p_iter-1)+1:2*N*p_iter) = dummyvec2N
                            end do
                        case(FD_SCHEME_BACKWARD)
                            call log("Backwards finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        case(FD_SCHEME_CENTRAL)
                            call log("Central finite difference not coded.", prefix=C_NOTCODED_ERROR, idt_level=3)
                            stop
                        end select
                end if

                ! Compute Rayleigh Matrix X^T H X
                rq_matrix = 0.0d0
                do p_iter=1,p
                    do p_iter2=1,p
                        ! col idx, row idx
                        rq_matrix(p_iter,p_iter2) = DOT_PRODUCT(X(:,p_iter),HX_findiff(2*N*(p_iter2-1)+1:2*N*p_iter2))
                    end do
                end do
                ! The Rayleigh Quotient is the trace of the subspace matrix X^T H X
                rq = 0.0d0
                do p_iter=1,p
                    rq = rq + rq_matrix(p_iter,p_iter)
                end do
                ! The gradient is given by H X - X (X^T H X) and is a 3N x p Matrix
                rq_gradient = HX_findiff
                do p_iter=1,p
                    do p_iter2=1,p
                        rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter) &
                                & - rq_matrix(p_iter2,p_iter) * X(:,p_iter2)
                    end do
                end do
                rq_gradient = rq_gradient * 2.0d0

                ! Change representation of Gradient:
                do p_iter=1,p
                    rq_gradient_2d(:,p_iter) = rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter)
                end do
                ! Compute the norm of the Gradient which is given by the Frobenius Norm:
                ! sqrt(tr(G^T*G)) = sqrt(g_1^2 + ... + g_p^2)
                rq_gradient_norm = 0.0d0
                do p_iter=1,p
                    rq_gradient_norm = rq_gradient_norm + DOT_PRODUCT(rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter),&
                            & rq_gradient(2*N*(p_iter-1)+1:2*N*p_iter))
                end do
                ! For our solvers we need the force representation
                rq_force = -1.0d0 * rq_gradient

                if (present(info_out_file)) then
                    ! If output information is requested compute the Rayleigh Ritz procedure in every iteration
                    rayleighritz_evec = rq_matrix
                    call syev(rayleighritz_evec,evals,jobz="V")
                    v_fin = 0.0d0
                    do p_iter=1,p
                        dummyvec2N = 0.0d0
                        do p_iter2=1,p
                            dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                        end do
                        call vec_to_embeddingspace(basis, dummyvec2N, v_fin(:,p_iter))
                    end do
                    if (present(v_fin_exact)) then
                        do p_iter=1,p
                            comp_w_exact(p_iter) = DOT_PRODUCT(v_fin(:,p_iter),v_fin_exact(:,p_iter))
                        end do
                    end if
                end if

                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        !Update the velocity based on the force for the next iteration and the rotated b
                        call vpo_update_velocity(rq_force, vpo_b, vpo_velocity, rqm_settings%rqm_solver%vpo_mass, &
                                & rqm_settings%rqm_solver%vpo_dt)
                end select
            end do

            rayleighritz_evec = rq_matrix
            call syev(rayleighritz_evec,evals,jobz="V")
            v_fin = 0.0d0
            do p_iter=1,p
                dummyvec2N = 0.0d0
                do p_iter2=1,p
                    dummyvec2N = dummyvec2N + X(:,p_iter2) * rayleighritz_evec(p_iter2,p_iter)
                end do
                call vec_to_embeddingspace(basis,N, dummyvec2N, v_fin(:,p_iter))
            end do

            if (i_info_out_local) then
                write(io,'((I10,3x),3(E23.12E3,3x))',advance="no") iter, rq, rq_gradient_norm
                do p_iter=1,p
                    write(io,"(1(E23.12E3,3x))",advance="no") evals(p_iter)
                end do
                if (present(v_fin_exact)) then
                    do p_iter=1,p
                        write(io,"(1(E23.12E3,3x))",advance="no") comp_w_exact(p_iter)
                    end do
                end if
                if (rqm_settings%fd_settings%i_richardson) then
                    do p_iter=1,p
                        write(io,'(1(E23.12E3,3x))',advance="no") fd_step_used(p_iter)
                    end do
                end if
                select case(rqm_settings%rqm_solver%solver)
                    case(INTEGRATOR_VPO)
                        write(io,'(1(E23.12E3,3x))') norm(vpo_velocity)
                    case(INTEGRATOR_LBFGS)
                        write(io,'(1(E23.12E3,3x))') lbfgs_steplength
                    case(INTEGRATOR_EULER)
                        write(io,'(1(E23.12E3,3x))') 0.0
                end select
            end if

            if (i_info_out_local) close(io)
        end subroutine rqm_grassmann_subsystem_information

end module spk_rqm