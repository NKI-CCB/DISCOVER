module logspace

  use types, only: dp

  implicit none

contains
  
  impure elemental function log1p(x) result (y)
    real(dp), intent(in) :: x
    real(dp) :: y
    real(dp), volatile :: u, v
  
    u = 1.0d0 + x
    if (u == 0) then
       y = log(u)
    else
       v = u - 1.0d0
       y = log(u) - (v - x) / u    ! cancels errors with IEEE arithmetic
    end if
  end function log1p


  function logadd(x, y) result (sum)
    use ieee_constants, only: isneginf

    real(dp), intent(in) :: x, y
    real(dp) :: sum

    if (isneginf(x)) then
       sum = y
    else if (isneginf(y)) then
       sum = x
    else
       sum = max(x, y) + log1p(exp(-abs(x - y)))
    end if
  end function logadd

end module logspace
