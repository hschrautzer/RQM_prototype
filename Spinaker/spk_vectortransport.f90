module spk_vectortransport
    !> Compute vector transports for optimization algorithms Riemannian Manifolds, which is about moving a vector from
    !> the tangent space of iterate X_k to the tangent space of X_k+1. The two iterates are connected by a Linesearch
    !> X_k+1 = X_k + t_k * DX_k.
    !> The vector transports depend on the specific manifold. The prefixes describe the manifold:
    !> GM: Grassmann Manifold
    !> The parallel transport is the "exact" way to do that. One alternative are differentiated retractions.
    use spk_precision_const, only: xp
    implicit none
    private
    interface GM_paralleltransport_applyTM
        module procedure GM_paralleltransport_applyTM_1D
        module procedure GM_paralleltransport_applyTM_2D
    end interface GM_paralleltransport_applyTM

    public :: GM_paralleltransport_calcTM, GM_paralleltransport_applyTM
    contains
        subroutine GM_paralleltransport_calcTM(X,t,U,S,VT,TM)
            !> The formula for parallel transport of a vector b on the tangent space of a point X on the Grassmannian
            !> along the geodesic defined by the direction DX=USV^T with steplength t is:
            !> b' = b - [X*V*sin(S*t)+U*cos(1-S*t)]*U^T*b
            !> b' = b - TM*U^T*b
            !> This routine calculates the matrix TM. The application of this to a vector b can be done by use of the
            !> routine: GM_paralleltransport_applyTM
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: U(:,:),S(:),VT(:,:)            !< Singular value decomposition of DX
            real(kind=xp) :: t                                          !< Step length
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: X(:,:)                      !< Representative of equivalence class of n by
                                                                        !< p matrices on the Stiefel Manifold
                                                                        !< corresponding to the same point on the Grass-
                                                                        !< mannian.
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: TM(:,:)                       !< Transport matrix and matrix U from DX = USV^T
            !===================================Local Variable==========================================================
            real(kind=xp), allocatable :: dummyvec(:)
            real(kind=xp), allocatable :: XVT(:,:)
            integer :: n, p, shapeX(2)
            integer :: i_p, i_p2

            shapeX = shape(X)
            n = shapeX(1)
            p = shapeX(2)
            allocate(dummyvec(n))
            allocate(XVT(n,p))
            TM = 0.0d0
            XVT = 0.0d0
            !do i_p=1,p
            !    do i_p2=1,p
            !        XVT(:,i_p) = XVT(:,i_p) + X(:,i_p2) * VT(i_p,i_p2)
            !    end do
            !end do
            !do i_p=1,p
            !    TM(:,i_p) = U(:,i_p) * (1.0d0-cos(S(i_p)*t)) + XVT(:,i_p) * sin(S(i_p)*t)
            !end do
            do i_p=1,p
                dummyvec = 0.0d0
                do i_p2=1,p
                    ! We want to comput X * V (therefore we have to choose VT(i_n,i_n2) for V(i_n2,i_n))
                    dummyvec = dummyvec + X(:,i_p2) * VT(i_p,i_p2)
                end do
                TM(:,i_p) = U(:,i_p) * (1.0d0-cos(S(i_p)*t)) + dummyvec * sin(S(i_p)*t)
            end do
            deallocate(dummyvec,XVT)
        end subroutine GM_paralleltransport_calcTM

        subroutine GM_paralleltransport_applyTM_2D(b,T,U)
            !> Parallel transport of a vector b from the tangent space of the Grassmanian Gnp with a given transport
            !> matrix T and matrix U from the thin SVD of the geodesic direction.
            !> b' = b - T*(U^T*b)
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: U(:,:)                     !< (n x p) matrix from the thin SVD of the geodesic
            real(kind=xp), intent(in) :: T(:,:)                     !< (n x p) transport matrix
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: b(:,:)                  !< Vector (n x p) in tangent space of Gnp at point X
                                                                    !< which is implicitly included in the transport
                                                                    !< matrix
            !===================================Local Variable==========================================================
            integer :: n, p, shapeB(2)
            integer :: i_p, i_p2
            real(kind=xp), allocatable :: UTB(:,:)                  !< Product U^T B
            real(kind=xp), allocatable :: bnew(:,:)
            shapeB = shape(b)
            n = shapeB(1)
            p = shapeB(2)
            allocate(UTB(p,p),bnew(n,p))
            bnew = 0.0d0
            do i_p=1,p
                do i_p2=1,p
                    UTB(i_p,i_p2) = DOT_PRODUCT(U(:,i_p),b(:,i_p2))
                end do
            end do
            do i_p=1,p
                do i_p2=1,p
                    bnew(:,i_p) = b(:,i_p) - T(:,i_p2) * UTB(i_p2,i_p)
                end do
            end do
            b = bnew
            deallocate(UTB,bnew)
        end subroutine GM_paralleltransport_applyTM_2D

        subroutine GM_paralleltransport_applyTM_1D(b,T,U)
            !> Parallel transport of a vector B from the tangent space of the Grassmanian Gnp with a given transport
            !> matrix T and matrix U from the thin SVD of the geodesic direction.
            !> Btrans = b - T*b(U^T*b)
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: U(:,:)                     !< (n x p) matrix from the thin SVD of the geodesic
            real(kind=xp), intent(in) :: T(:,:)                     !< (n x p) transport matrix
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: b(:)                    !< Vector (n*p) in tangent space of Gnp at point X
                                                                    !< which is implicitly included in the transport
                                                                    !< matrix
            !===================================Local Variable==========================================================
            integer :: n, p, shapeU(2)
            integer :: i_p, i_p2
            real(kind=xp), allocatable :: UTB(:,:)                  !< Product U^T B
            real(kind=xp), allocatable :: dummyvec(:)
            real(kind=xp), allocatable :: bnew(:)
            shapeU = shape(U)
            n = shapeU(1)
            p = shapeU(2)
            allocate(UTB(p,p),dummyvec(n),bnew(n*p))
            do i_p=1,p
                do i_p2=1,p
                    UTB(i_p,i_p2) = DOT_PRODUCT(U(:,i_p),b(n*(i_p2-1)+1:n*i_p2))
                end do
            end do
            bnew = 0.0d0
            do i_p=1,p
                dummyvec = 0.0d0
                do i_p2=1,p
                    dummyvec = dummyvec + T(:,i_p2) * UTB(i_p2,i_p)
                end do
                bnew(n*(i_p-1)+1:n*i_p) = b(n*(i_p-1)+1:n*i_p) - dummyvec
            end do
            b = bnew
            deallocate(UTB,dummyvec,bnew)
        end subroutine GM_paralleltransport_applyTM_1D

end module spk_vectortransport