module spk_retraction
    !> Module contains reatraction routines for optimizations on Riemannian Manifolds, which is a generalization of the
    !> Linesearch procedure X_k+1 = X_k + t_k * DX_k. On Riemmannian Manifolds the exact way to perform this update
    !> is the exponential map, which might have closed forms depending on the specific manifolds. The other way is to
    !> think of the current iterate and the search direction eta_k, which is in the tangent space of X_k, of being em-
    !> bedded in the euclidean vectorspace. Then we can perform the addition in the "usual" euclidean way and afterwards
    !> "project" back to the manifold. These types of retraction can be viewed as approximations to the exponential map.
    !> The retractions depend on the specific manifold. The prefixes describe the manifold:
    !> GM: Grassmann Manifold
    USE, INTRINSIC :: IEEE_ARITHMETIC
    use spk_precision_const, only: xp
    use spk_abbreviations_const, only: C_VALUE_ERROR
    use spk_qr, only: qr_gramschmidt, qr_householder, qr_decomposition
    use spk_singularvalues, only: svd
    use spk_vector_lina, only: norm
    use spk_logging, only: log
    implicit None
    private

    public :: GM_retract_qr, GM_retract_polar, GM_retract_exp
    contains

        subroutine GM_retract_qr(X,method)
            !> Performs a Retraction on the Grassmannian of a new Iterate X (X = X_k+1) by using QR-decomposition.
            !===================================Input Variable==========================================================
            character(len=1), intent(in) :: method                      !< Method used for QR decomposition:
                                                                        !< G: Gram-Schmidt
                                                                        !< H: Householder Transformation
                                                                        !< M: Wrapper to MKL Lapack QR decomposition
                                                                        !< If any character is provided the Lapack impl.
                                                                        !< is chosen.
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: X(:,:)                      !< Representative of equivalence class of n by
                                                                        !< p matrices on the Stiefel Manifold
                                                                        !< corresponding to the same point on the Grass-
                                                                        !< mannian.
            !===================================Local Variable==========================================================
            real(kind=xp), allocatable :: X_new(:,:)
            integer :: n,p,shapeX(2),i
            shapeX = shape(X)
            n = shapeX(1)
            p = shapeX(2)
            allocate(X_new(n,p))
            if (method=="G") then
                call qr_gramschmidt(X,X_new)
            elseif (method=="H") then
                call qr_householder(X,X_new)
            else
                call qr_decomposition(X,X_new)
            end if
            X = X_new
            do i=1,p
                if (IEEE_IS_NAN(norm(X_new(:,i)))) then
                    call log("NaN vector during QR-decomp. Duplicate vector?", prefix=C_VALUE_ERROR,idt_level=2)
                    stop
                end if
            end do
            deallocate(X_new)
        end subroutine GM_retract_qr

        subroutine GM_retract_polar(X)
            !> Performs a Retraction on the Grassmannian of a new Iterate X (X = X_k+1) by using Polar-decomposition,
            !> which means performing a (thin/compact) Singular Value Decomposition and using X=UV^T as the new iterate
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: X(:,:)                      !< Representative of equivalence class of n by
                                                                        !< p matrices on the Stiefel Manifold
                                                                        !< corresponding to the same point on the Grass-
                                                                        !< mannian.
            !===================================Local Variable==========================================================
            integer :: n, p, shapeX(2), i_p, i_p2
            real(kind=xp), allocatable :: U(:,:), S(:), VT(:,:)
            shapeX = shape(X)
            n = shapeX(1)
            p = shapeX(2)
            allocate(U(n,p),S(p),VT(p,p))
            call svd(X,U,S,VT)
            X = 0.0d0
            do i_p=1,p
                do i_p2=1,p
                    X(:,i_p) = X(:,i_p) + U(:,i_p2) * VT(i_p,i_p2)
                end do
            end do
            deallocate(U,S,VT)
        end subroutine GM_retract_polar
    
        subroutine GM_retract_exp(X,DX,t,U,S,VT,postmultiplication)
            !> Calculates the exponential map for the Grassmann Manifold and move along the geodesic defined by it.
            !> The formula is given by
            !> X' = X*V*cos(S*t)*V^T+U*sin(S*t)*V^T
            !> where DX = U * S * V^T is the thin Singular Value Decomposition of the direction of the geodesic defined
            !> by DX in the tangent space of X.
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: t                              !< Step Length for translating along geodesic
            real(kind=xp), intent(in) :: DX(:,:)                        !< Direction of the geodesic, DX is in the
                                                                        !< tangent space of X
            character(len=1), intent(in) :: postmultiplication          !< Post Multiplication with VT (orthonormal)
                                                                        !< does not change the point on the Grassmann
                                                                        !< Manifold only its representative of the
                                                                        !< equivalence class. So unless if this arg. is
                                                                        !< set to "Y" the post-multiplication wont be
                                                                        !< performed.
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: U(:,:), S(:), VT(:,:)         !< Thin SVD of DX
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: X(:,:)                      !< Representative of equivalence class of n by
                                                                        !< p matrices on the Stiefel Manifold
                                                                        !< corresponding to the same point on the Grass-
                                                                        !< mannian.
            !===================================Local Variable==========================================================
            real(kind=xp), allocatable :: B(:,:)
            real(kind=xp), allocatable :: XV(:,:)
            integer :: n, p, shapeX(2)
            integer :: i_p, i_p2
            shapeX = shape(X)
            n = shapeX(1)
            p = shapeX(2)
            allocate(B(n,p))
            allocate(XV(n,p))
            ! First compute the compact singular value decomposition of the direction defining the transport (DX):
            ! Left singular vector U, singular values S as well as right singular vectors V will be allocated within SVD
            U=0.0d0
            S=0.0d0
            VT=0.0d0
            call svd(DX,U,S,VT)
            B = 0.0d0
            XV = 0.0d0
            do i_p=1,p
                do i_p2=1,p
                    !B(:,i_p) = B(:,i_p) + X(:,i_p2) * VT(i_p,i_p2) * cos(S(i_p2)*t) + U(:,i_p2) * sin(S(i_p2)*t)
                    XV(:,i_p) = XV(:,i_p) + X(:,i_p2) * VT(i_p,i_p2)
                end do
            end do
            do i_p=1,p
                B(:,i_p) = B(:,i_p) + XV(:,i_p) * cos(S(i_p)*t)+U(:,i_p)*sin(S(i_p)*t)
            end do

            if (postmultiplication=="Y") then
                X = 0.0d0
                do i_p=1,p
                    do i_p2=1,p
                        X(:,i_p) = X(:,i_p) + B(:,i_p2) * VT(i_p2,i_p)
                    end do
                end do
            else
                X = B
            end if

        end subroutine GM_retract_exp

end module spk_retraction