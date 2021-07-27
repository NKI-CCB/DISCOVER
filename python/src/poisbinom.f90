module poisbinom

  use types, only: dp

  implicit none

contains

  function cdf(p, x)
    use logspace, only: logadd, log1p

    real(dp), dimension(:), intent(in) :: p
    integer, intent(in) :: x
    real(dp) :: cdf

    real(dp), dimension(size(p)) :: logp, lognotp
    real(dp), dimension(2, size(p)) :: memory

    integer :: j, k, n, prev, curr, tmp

    n = size(p)

    if (x > n) then
      cdf = 1
    else if (x < 0) then
      cdf = 0
    else
      logp = log(p)
      lognotp = log1p(-p)

      memory(1, 1) = lognotp(1)
      do j = 2, size(p)
        memory(1, j) = lognotp(j) + memory(1, j - 1)
      end do

      prev = 1
      curr = 2

      cdf = memory(1, n)

      do k = 1, x
         if (k > 1) then
            memory(curr, k) = logp(k) + memory(prev, k - 1)
         else
            memory(curr, k) = logp(k)
         end if

         do j = k + 1, n
            memory(curr, j) = logadd(lognotp(j) + memory(curr, j - 1), &
                                     logp(j) + memory(prev, j - 1))
         end do

         cdf = logadd(cdf, memory(curr, n))
         if (cdf >= 0) then
            cdf = 0
            exit
         end if

         tmp = curr
         curr = prev
         prev = tmp
      end do

      cdf = exp(cdf)
    end if
  end function cdf


  function fullCdf(p, x) result (cdf)
    use logspace, only: logadd, log1p

    real(dp), dimension(:), intent(in) :: p
    integer, intent(in) :: x
    real(dp), dimension(0:x) :: cdf

    real(dp), dimension(size(p)) :: logp, lognotp
    real(dp), dimension(2, size(p)) :: memory

    integer :: j, k, n, prev, curr, tmp

    n = size(p)

    if (x > n) then
      cdf = 1.0
    else
      logp = log(p)
      lognotp = log1p(-p)

      memory(1, 1) = lognotp(1)
      do j = 2, size(p)
        memory(1, j) = lognotp(j) + memory(1, j - 1)
      end do

      prev = 1
      curr = 2

      cdf(0) = memory(1, n)

      do k = 1, x
         if (k > 1) then
            memory(curr, k) = logp(k) + memory(prev, k - 1)
         else
            memory(curr, k) = logp(k)
         end if

         do j = k + 1, n
            memory(curr, j) = logadd(lognotp(j) + memory(curr, j - 1), &
                                     logp(j) + memory(prev, j - 1))
         end do

         cdf(k) = logadd(cdf(k - 1), memory(curr, n))
         if (cdf(k) >= 0) then
            cdf(k:) = 0
            exit
         end if

         tmp = curr
         curr = prev
         prev = tmp
      end do

      cdf = exp(cdf)
    end if
  end function fullCdf

end module poisbinom
