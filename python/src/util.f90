module util

  use types, only: dp

  implicit none

  interface cumsum
     module procedure cumsum_dp
     module procedure cumsum_i
  end interface cumsum

contains

  function bincount(x)
    integer, dimension(:), intent(in) :: x
    integer, dimension(:), allocatable :: bincount

    integer :: i, k

    allocate(bincount(maxval(x)))
    bincount = 0

    do i = 1, size(x)
       k = x(i)
       bincount(k) = bincount(k) + 1
    end do
  end function bincount


  function cumsum_i(x) result (cumsum)
    integer, dimension(:), intent(in) :: x
    integer, dimension(size(x)) :: cumsum

    integer :: i

    cumsum(1) = x(1)
    do i = 2, size(x)
       cumsum(i) = cumsum(i - 1) + x(i)
    end do
  end function cumsum_i


  function cumsum_dp(x) result (cumsum)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: cumsum

    integer :: i

    cumsum(1) = x(1)
    do i = 2, size(x)
       cumsum(i) = cumsum(i - 1) + x(i)
    end do
  end function cumsum_dp


  subroutine cummin(x)
    real(dp), dimension(:), intent(inout) :: x

    integer :: i

    do i = size(x) - 1, 1, -1
       x(i) = min(x(i), x(i + 1))
    end do
  end subroutine cummin


  subroutine unique(x, values, inverse)
    use m_uniinv, only: uniinv

    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), allocatable, intent(out) :: values
    integer, dimension(size(x)), intent(out) :: inverse

    integer :: i

    call uniinv(x, inverse)
    allocate(values(maxval(inverse)))
    
    do i = 1, size(inverse)
       values(inverse(i)) = x(i)
    end do
  end subroutine unique


  function exponentialSearch(values, key, offset) result (position)
    real(dp), dimension(:), intent(in) :: values
    real(dp), intent(in) :: key
    integer, intent(in) :: offset
    integer :: position

    integer :: lbound, ubound, mid, window

    lbound = offset
    window = min(1, size(values) - offset)  !! assert offset <= size(values) and key <= maxval(values)

    do while (values(lbound + window) < key)
       lbound = lbound + window
       window = min(window * 2, size(values) - lbound)
    end do

    ubound = lbound + window

    do while (lbound <= ubound)
       mid = lbound + (ubound - lbound) / 2
       if (values(mid) >= key) then
          ubound = mid - 1
       else
          lbound = mid + 1
       end if
    end do

    position = lbound
  end function exponentialSearch


  function insertionPoints(keys, values)
    real(dp), dimension(:), intent(in) :: keys
    real(dp), dimension(:), intent(in) :: values
    integer, dimension(size(keys)) :: insertionPoints

    integer :: i, offset

    offset = 1

    do i = 1, size(keys)
       if (keys(i) <= values(size(values))) then
          insertionPoints(i) = exponentialSearch(values, keys(i), offset)
          offset = insertionPoints(i)
       else
          insertionPoints(i) = -1
       end if
    end do
  end function insertionPoints

end module util
