module spk_qr
    !> Module contains routines for QR-decompositions.
    use spk_precision_const, only: xp
    use spk_vector_lina, only: norm
    use spk_logging, only: log
    use spk_abbreviations_const, only: C_PRECISION_ERROR
    implicit none
    external dgeqrf, dorgqr
    private
    interface qr_decomposition
        module procedure qr_decomposition
    end interface qr_decomposition

    public :: qr_gramschmidt, qr_householder, qr_decomposition
    contains
        subroutine qr_decomposition(A,Q)
            !> Performs QR Decomposition of rectangular m x n matrix A using geqrf (no pivoting)
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: A(:,:)                         !< Input matrix to be factored
            !===================================Output Variable=========================================================
            real(kind=xp), allocatable, intent(out) :: Q(:,:)           !< Q Factor of QR decomposition
            !===================================Local Variable==========================================================
            real(kind=xp), allocatable :: A_local(:,:)                  !< Local copy of A
            integer :: lda                                              !< Leading dimension of A
            real(kind=xp), allocatable :: tau(:)                        !< Contains scalars that define elementary
                                                                        !< reflectors for the matrix Q in its decomp.
                                                                        !< in a product of elementary reflectors.
            real(kind=xp), allocatable :: work(:)                       !< work array
            integer :: lwork                                            !< workspace query
            integer :: info
            integer :: m,n,shapeA(2)
            shapeA=shape(A)
            m = shapeA(1)
            n = shapeA(2)

            lda = m
            allocate(A_local(m,n),tau(n))
            A_local = A
            ! Workspace query
            allocate(work(1))
            lwork = -1
            call dgeqrf(m, n, A_local, lda, tau, work, -1, info)
            lwork = max(int(work(1)),1)
            deallocate(work)
            allocate(work(lwork))
            call dgeqrf(m, n, A_local, lda, tau, work, lwork, info)
            deallocate(work)
            allocate(work(1))
            call dorgqr(m, n, n, A_local, lda, tau, work, -1, info)
            lwork= max(1,int(work(1)))
            deallocate(work)
            allocate(work(lwork))
            call dorgqr(m, n, n, A_local, lda, tau, work, lwork, info)
            deallocate(work)
            Q=A_local
        end subroutine qr_decomposition

        subroutine qr_householder(A,Q)
            !> Performs QR using Householder Transformation (not really tested)
            !========================================Input Variable=====================================================
            real(kind=xp), intent(in) :: A(:,:)                         !< Tall skinny matrix
            !========================================Output Variable====================================================
            real(kind=xp), intent(out) :: Q(:,:)                        !< Q-Factor of QR Decomposition
            !========================================Local Variable=====================================================
            real(kind=xp),allocatable :: A_local(:,:)
            real(kind=xp),allocatable :: R(:,:)
            integer :: i,j,k
            real(kind=xp) :: alpha, normx
            real(kind=xp),allocatable :: v(:)
            integer :: m, n, shapeA(2)
            shapeA=shape(A)
            m = shapeA(1)
            n = shapeA(2)
            allocate(A_local(m,n),R(m,n),v(m))
            A_local = A
            Q = 0.0d0
            R = 0.0d0

            do k = 1, n
                ! Compute the norm of the k-th column
                normx = 0.0d0
                do i = k, m
                    normx = normx + A_local(i, k)**2
                end do
                normx = sqrt(normx)

                ! Calculate the sign of the first element
                alpha = -sign(normx, A_local(k, k))

                ! Construct the Householder vector
                v = 0.0d0
                v(k) = A_local(k, k) - alpha
                do i = k + 1, m
                    v(i) = A_local(i, k)
                end do

                ! Normalize the Householder vector
                v = v / sqrt(sum(v**2))

                ! Update A and accumulate Q
                do j = k, n
                    R(k, j) = dot_product(v(k:m), A_local(k:m, j))
                    A_local(k:m, j) = A_local(k:m, j) - 2.0d0 * v(k:m) * R(k, j)
                end do

                do i = k + 1, m
                    Q(i, k) = v(i)
                end do
            end do

            deallocate(A_local,R,v)
        end subroutine qr_householder

        subroutine qr_gramschmidt(A,Q)
            !> Will compute the orthonormalized matrix Q based on the matrix A using basic gram-schmidt procedure. This
            !> can be numerically unstable and it is recommended to use Householder-Reflections
            !========================================Input Variable=====================================================
            real(kind=xp), intent(in) :: A(:,:)                         !< Input Matrix A of dimensions n x p
            !========================================Output Variable====================================================
            real(kind=xp), intent(out) :: Q(:,:)                        !< Output matrix Q (n x p)
            !========================================Local Variable=====================================================
            integer :: shapeA(2), n, p
            integer :: i,j
            real(kind=xp), allocatable :: r(:), s(:)
            shapeA = shape(A)
            n = shapeA(1)
            p = shapeA(2)
            allocate(r(n),s(n))
            Q=0.0d0
            do i=1,p
                r = 0.0d0
                do j=1,i-1
                    r = r + DOT_PRODUCT(Q(:,j),A(:,i)) * Q(:,j)
                end do
                s = A(:,i) - r
                Q(:,i) = s / norm(s)
            end do
            deallocate(r,s)
        end subroutine qr_gramschmidt

end module spk_qr