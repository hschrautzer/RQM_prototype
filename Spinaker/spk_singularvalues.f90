module spk_singularvalues
    !> This module contains routines for computing Singular value decompositions. Mostly wrappers to Lapack Routines
    use spk_vector_lina, only: norm
    use spk_precision_const, only: xp
    use lapack95, only: gesvd
    implicit none
    private
    interface svd
        !> Interface for singular value decompositions
        module procedure svd_rUSVT                  !< Returns U, S and VT
    end interface svd
    public :: svd
    contains
        subroutine svd_rUSVT(A,U,S,VT)
            !> Computes the compact singular value decomposition of a rectangular Matrix A of size n x p with n>=p and
            !> return left singular vector matrix U (n,p) and right singular vector matrix V (p,p) as well as the
            !> singular values in an array in decreasing order.
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: A(:,:)                         !< n by p matrix A
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: U(:,:), S(:), VT(:,:)         !< See documentation of gesvd
            !===================================Local Variable==========================================================
            integer :: m, n, shapeA(2)
            real(kind=xp), allocatable :: l_A(:,:)
            shapeA = shape(A)
            m = shapeA(1)
            n = shapeA(2)
            allocate(l_A(m,n))
            l_A = A
            ! Call Fortran 95 Interface (no work arrays needed)
            call gesvd(l_A,S,U,VT)
        end subroutine svd_rUSVT
    
end module spk_singularvalues