subroutine estimateBackground(numRowValues, rowValues, rowWeights, numColValues, colValues, colWeights, mu)
  use types, only: dp
  use maxent, only: fit

  integer, intent(in) :: numRowValues
  integer, dimension(numRowValues), intent(in) :: rowValues
  integer, dimension(numRowValues), intent(in) :: rowWeights
  integer, intent(in) :: numColValues
  integer, dimension(numColValues), intent(in) :: colValues
  integer, dimension(numColValues), intent(in) :: colWeights
  real(dp), dimension(numRowValues + numColValues), intent(out) :: mu

  call fit(rowValues, rowWeights, colValues, colWeights, mu)
end subroutine estimateBackground


subroutine analyseBlockStructure(numFeatures, numSamples, data, bg, lowerTail, numBlocks, blockLengths, estimate_fdr, pMatrix, qMatrix, pi0)
  use types, only: dp
  use fdr, only: analyseBlockStructure_ => analyseBlockStructure

  integer, intent(in) :: numFeatures
  integer, intent(in) :: numSamples
  integer, dimension(numFeatures, numSamples), intent(in) :: data
  real(dp), dimension(numFeatures, numSamples), intent(in) :: bg
  logical, intent(in) :: lowerTail
  integer, intent(in) :: numBlocks
  integer, dimension(numBlocks), intent(in) :: blockLengths
  logical, intent(in) :: estimate_fdr
  real(dp), dimension(numFeatures, numFeatures) :: pMatrix
  real(dp), dimension(numFeatures, numFeatures), intent(out) :: qMatrix
  real(dp), intent(out) :: pi0

  pMatrix = analyseBlockStructure_(data, bg, lowerTail, blockLengths, estimate_fdr, qMatrix, pi0)
end subroutine analyseBlockStructure


subroutine mutex(numFeatures, numSamples, data, bg, lowerTail, estimate_fdr, p, q, pi0)
  use types, only: dp
  use fdr, only: mutex_ => mutex

  integer, intent(in) :: numFeatures
  integer, intent(in) :: numSamples
  integer, dimension(numFeatures, numSamples), intent(in) :: data
  real(dp), dimension(numFeatures, numSamples), intent(in) :: bg
  logical, intent(in) :: lowerTail
  logical, intent(in) :: estimate_fdr
  real(dp), dimension(numFeatures * (numFeatures - 1) / 2), intent(out) :: p
  real(dp), dimension(numFeatures * (numFeatures - 1) / 2), intent(out) :: q
  real(dp), intent(out) :: pi0

  p = mutex_(data, bg, lowerTail, estimate_fdr, q, pi0)
end subroutine mutex


subroutine ppoisbinom(q, size, probs, p)
  use poisbinom, only: cdf
  use types, only: dp

  integer, intent(in) :: q
  integer, intent(in) :: size
  real(dp), dimension(size), intent(in) :: probs
  real(dp), intent(out) :: p

  p = cdf(probs, q)
end subroutine ppoisbinom


subroutine colLogSumExps(x, numRows, numCols, logSums)
  use types, only: dp

  real(dp), dimension(numRows, numCols), intent(in) :: x
  integer, intent(in) :: numRows, numCols
  real(dp), dimension(numCols), intent(out) :: logSums

  integer :: j

  if (numRows == 1) then
     logSums = x(1, :)
  else
     do j = 1, numCols
        logSums(j) = logSumExp(x(:, j))
     end do
  end if

contains

  function logSumExp(x)
    use types, only: dp
    use logspace, only: log1p

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: logSumExp

    real(dp) :: x_max
    integer :: i_max

    integer :: i
    real(dp) :: sumExp

    i_max = maxloc(x, 1)
    x_max = x(i_max)

    sumExp = 0
    do i = 1, size(x)
       if (i /= i_max) then
          sumExp = sumExp + exp(x(i) - x_max)
       end if
    end do

    logSumExp = x_max + log1p(sumExp)
  end function logSumExp

end subroutine colLogSumExps
