module fdr

  use types, only: dp

  implicit none

contains

  subroutine coordsFromLinear(n, k, i, j)
    integer, intent(in) :: n, k
    integer, intent(out) :: i, j

    i = int(n - (sqrt(real(-8 * k + (2 * n - 1)**2 + 8, dp)) - 1) / 2)
    j = i * (i + 1) / 2 - i * n + k + n
  end subroutine coordsFromLinear


  function poisbinomTest(data, bg, lowerTail) result(p)
    use poisbinom, only: cdf

    integer, dimension(:, :), intent(in) :: data
    real(dp), dimension(size(data, 1), size(data, 2)), intent(in) :: bg
    logical, intent(in) :: lowerTail

    real(dp), dimension(size(data, 1) * (size(data, 1) - 1) / 2) :: p
    
    integer :: i, j, k, n

    n = size(data, 1)

    !$omp parallel do private(i, j)
    do k = 1, n * (n - 1) / 2
       call coordsFromLinear(n, k, i, j)    ! k = j + n * (i - 1) - i * (i + 1) / 2
       if (lowerTail) then
          p(k) = cdf(bg(i, :) * bg(j, :), dot_product(data(i, :), data(j, :)))
       else
          p(k) = 1 - cdf(bg(i, :) * bg(j, :), dot_product(data(i, :), data(j, :)) - 1)
       end if
    end do
    !$omp end parallel do
  end function poisbinomTest


  function expectedP(cdf)
    real(dp), dimension(0:), intent(in) :: cdf
    real(dp) :: expectedP

    integer :: i

    expectedP = cdf(0)**2
    do i = 1, ubound(cdf, 1)
       expectedP = expectedP + cdf(i) * (cdf(i) - cdf(i - 1))
    end do
  end function expectedP


  function mutex(data, bg, lowerTail, estimate_fdr, q, pi0) result (p)
    use util, only: bincount, cummin, cumsum, unique

    integer, dimension(:, :), intent(in) :: data
    real(dp), dimension(size(data, 1), size(data, 2)), intent(in) :: bg
    logical, intent(in) :: lowerTail
    logical, intent(in) :: estimate_fdr
    real(dp), intent(out) :: pi0

    real(dp), dimension(size(data, 1) * (size(data, 1) - 1) / 2) :: p

    real(dp), dimension(size(p)), intent(out) :: q
    
    real(dp), dimension(:), allocatable :: sortedLevels
    integer, dimension(size(p)) :: inv
    integer, dimension(:), allocatable :: counts
    real(dp), dimension(:), allocatable :: qUnique

    real(dp) :: expectedP0

    integer :: i, j, n, t

    p = poisbinomTest(data, bg, lowerTail)

    if (estimate_fdr) then
       call unique(p, sortedLevels, inv)
       allocate(qUnique(size(sortedLevels)))

       qUnique = 0
       expectedP0 = 0

       n = size(data, 1)

       !$omp parallel do private(i, j) reduction(+:qUnique,expectedP0)
       do t = 1, n * (n - 1) / 2
          call coordsFromLinear(n, t, i, j)
          call updateQ(bg(i, :) * bg(j, :), lowerTail, sortedLevels, qUnique, expectedP0)
       end do
       !$omp end parallel do

       counts = cumsum(bincount(inv))
       qUnique = cumsum(qUnique) / counts
       where (qUnique > 1) qUnique = 1
       call cummin(qUnique)

       do i = 1, size(inv)
          q(i) = qUnique(inv(i))
       end do

       print *, "expected P:", expectedP0 / size(p)
       print *, "mean P:", sum(p) / size(p)
       pi0 = min(1.0_dp, sum(p) / expectedP0)
       print *, "pi0:", pi0
    end if
  end function mutex


  function computeP(chrom1, p1, chrom2, p2, lowerTail) result(p)
    use poisbinom

    integer, dimension(:, :), intent(in) :: chrom1
    real(dp), dimension(size(chrom1, 1), size(chrom1, 2)), intent(in) :: p1
    integer, dimension(:, :), intent(in) :: chrom2
    real(dp), dimension(size(chrom2, 1), size(chrom2, 2)), intent(in) :: p2
    logical, intent(in) :: lowerTail

    real(dp), dimension(size(chrom1, 1), size(chrom2, 1)) :: p

    integer :: i, j

    !$omp parallel do collapse(2)
    do i = 1, size(chrom1, 1)
       do j = 1, size(chrom2, 1)
          if (lowerTail) then
             p(i, j) = cdf(p1(i, :) * p2(j, :), dot_product(chrom1(i, :), chrom2(j, :)))
          else
             p(i, j) = max(1 - cdf(p1(i, :) * p2(j, :), dot_product(chrom1(i, :), chrom2(j, :)) - 1), 0.0_dp)
          end if
       end do
    end do
    !$omp end parallel do
  end function computeP


  subroutine updateQ(p, lowerTail, sortedLevels, qUnique, expectedP0)
    use poisbinom, only: fullCdf
    use util, only: insertionPoints

    real(dp), dimension(:), intent(in) :: p
    logical, intent(in) :: lowerTail
    real(dp), dimension(:), intent(in) :: sortedLevels
    real(dp), dimension(size(sortedLevels)), intent(inout) :: qUnique
    real(dp), intent(inout) :: expectedP0

    real(dp), dimension(0:size(p)) :: cdf
    integer, dimension(0:size(p)) :: k
    integer :: l

    if (lowerTail) then
       cdf = fullCdf(p, size(p))
    else
       cdf = fullCdf(1 - p, size(p))
    end if

    expectedP0 = expectedP0 + expectedP(cdf)

    k = insertionPoints(cdf, sortedLevels)
    qUnique(k(0)) = qUnique(k(0)) + cdf(0)

    do l = 1, ubound(k, 1)
       if (k(l) == -1) then
          exit
       else
          qUnique(k(l)) = qUnique(k(l)) + cdf(l) - cdf(l - 1)
       end if
    end do
  end subroutine updateQ


  subroutine updateMultiQ(bg1, bg2, lowerTail, sortedLevels, qUnique, expectedP0)
    real(dp), dimension(:, :), intent(in) :: bg1
    real(dp), dimension(:, :), intent(in) :: bg2
    logical, intent(in) :: lowerTail
    real(dp), dimension(:), intent(in) :: sortedLevels
    real(dp), dimension(size(sortedLevels)), intent(inout) :: qUnique
    real(dp), intent(inout) :: expectedP0

    integer :: i, j

    !$omp parallel do collapse(2) reduction(+:qUnique,expectedP0)
    do i = 1, size(bg1, 1)
       do j = 1, size(bg2, 1)
          call updateQ(bg1(i, :) * bg2(j, :), lowerTail, sortedLevels, qUnique, expectedP0)
       end do
    end do
    !$omp end parallel do
  end subroutine updateMultiQ


  function analyseBlockStructure(data, bg, lowerTail, blockLengths, estimate_fdr, qMatrix, pi0) result (pMatrix)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use util, only: bincount, cummin, cumsum, unique

    integer, dimension(:, :), intent(in) :: data
    real(dp), dimension(size(data, 1), size(data, 2)), intent(in) :: bg
    logical, intent(in) :: lowerTail
    integer, dimension(:), intent(in) :: blockLengths
    logical, intent(in) :: estimate_fdr
    real(dp), dimension(size(data, 1), size(data, 1)), intent(out) :: qMatrix
    real(dp), intent(out) :: pi0
    real(dp), dimension(size(data, 1), size(data, 1)) :: pMatrix

    real(dp), dimension(:), allocatable :: p
    integer, dimension(size(blockLengths)) :: offsets, flatOffsets
    integer :: i, start1, end1, start2, k1, k2

    real(dp), dimension(:), allocatable :: sortedLevels
    integer, dimension(:), allocatable :: inv
    integer, dimension(:), allocatable :: counts
    real(dp), dimension(:), allocatable :: qUnique
    real(dp), dimension(:), allocatable :: q
    real(dp) :: expectedP0


    ! size(blockLengths) should be > 1

    !offsets(1) = 0
    !offsets(2:) = cumsum(blockLengths)(2:)
    offsets = cumsum(blockLengths)
    offsets(2:) = offsets(:size(offsets - 1)) + 1
    offsets(1) = 1

    flatOffsets(1) = 1
    flatOffsets(2:) = cumsum([ (blockLengths(i) * sum(blockLengths(i + 1:)), i = 1, size(blockLengths) - 1) ]) + 1

    allocate(p(sum([ (blockLengths(i) * sum(blockLengths(i+1:)), i = 1, size(blockLengths) - 1) ])))

    do i = 1, size(blockLengths) - 1
       start1 = offsets(i)
       end1 = offsets(i + 1) - 1
       start2 = offsets(i + 1)
       k1 = flatOffsets(i)
       k2 = flatOffsets(i + 1) - 1
       p(k1:k2) = reshape(computeP(data(start1:end1, :), bg(start1:end1, :), data(start2:, :), bg(start2:, :), lowerTail), [k2 - k1 + 1])
    end do

    if (estimate_fdr) then
       allocate(inv(size(p)))
       call unique(p, sortedLevels, inv)
       allocate(qUnique(size(sortedLevels)))

       qUnique = 0
       expectedP0 = 0
       do i = 1, size(blockLengths) - 1
          start1 = offsets(i)
          end1 = offsets(i + 1) - 1
          start2 = offsets(i + 1)
          call updateMultiQ(bg(start1:end1, :), bg(start2:, :), lowerTail, sortedLevels, qUnique, expectedP0)
       end do

       counts = cumsum(bincount(inv))
       qUnique = cumsum(qUnique) / counts
       where (qUnique > 1) qUnique = 1
       call cummin(qUnique)
       allocate(q(size(p)))

       do i = 1, size(inv)
          q(i) = qUnique(inv(i))
       end do

       print *, "expected P:", expectedP0 / size(p)
       print *, "mean P:", sum(p) / size(p)
       pi0 = min(1.0_dp, sum(p) / expectedP0)
       print *, "pi0:", pi0
    end if
    
    pMatrix = ieee_value(1.0_dp, ieee_quiet_nan)
    qMatrix = ieee_value(1.0_dp, ieee_quiet_nan)
    do i = 1, size(blockLengths) - 1
       start1 = offsets(i)
       end1 = offsets(i + 1) - 1
       start2 = offsets(i + 1)
       k1 = flatOffsets(i)
       k2 = flatOffsets(i + 1) - 1
       pMatrix(start1:end1, start2:) = reshape(p(k1:k2), [end1 - start1 + 1, size(data, 1) - start2 + 1])
       if (estimate_fdr) then
          qMatrix(start1:end1, start2:) = reshape(q(k1:k2), [end1 - start1 + 1, size(data, 1) - start2 + 1])
       end if
    end do
  end function analyseBlockStructure

end module fdr
