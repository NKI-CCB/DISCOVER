module ieee_constants

  use types, only: dp

  implicit none

  real(dp), parameter :: INFINITY = transfer(Z'7ff0000000000000', 1.0d0)
  real(dp), parameter :: NEGATIVE_INFINITY = transfer(Z'fff0000000000000', 1.0d0)
  real(dp), parameter :: NAN = transfer(Z'7ff8000000000000', 1.0d0)

  interface isneginf
     module procedure isneginf_dp
  end interface isneginf

contains

  pure function isneginf_dp(x)
    logical :: isneginf_dp
    real(dp), intent(in) :: x

    isneginf_dp = x == NEGATIVE_INFINITY
  end function isneginf_dp

end module ieee_constants
