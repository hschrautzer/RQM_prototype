module spk_finiteparams
    !> Contains derived types summarizing settings for finite difference calculations
    use spk_precision_const, only: xp
    use spk_algoconst, only: FD_SCHEME_FORWARD
    implicit none
    private
    type, public :: finiteparams
        integer :: scheme = FD_SCHEME_FORWARD       !< Forward, Central, Backward
        real(kind=xp) :: step = 1.0d-6              !< Denominator of final difference
        integer :: order = 1                        !< Stencil Points
        logical :: i_richardson = .False.           !< If Richardson Extrapolation is selected
        integer :: richardson_iter = 10             !< Maximum number of Richardson Extrapolation Iterations
        real(kind=xp) :: richardson_error = 1.0d-8  !< Maximum allowed deviation of Richardson Extrapolation
    end type finiteparams


end module spk_finiteparams