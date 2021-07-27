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
    real(dp), intent(in) :: x, y
    real(dp) :: sum

    real(dp), parameter :: LOG2 = log(2.0_dp)
    real(dp) :: x_minus_y

    if (x == y) then
       sum = x + LOG2
    else
       x_minus_y = x - y
       if (x_minus_y > 0) then
          sum = x + log1p(exp(-x_minus_y))
       else if (x_minus_y <= 0) then
          sum = y + log1p(exp(x_minus_y))
       else
          sum = x_minus_y
       end if
    end if
  end function logadd

end module logspace
