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


subroutine analyseBlockStructure(numFeatures, numSamples, data, bg, lowerTail, numBlocks, blockLengths, pMatrix, qMatrix, pi0)
  use types, only: dp
  use fdr, only: analyseBlockStructure_ => analyseBlockStructure

  integer, intent(in) :: numFeatures
  integer, intent(in) :: numSamples
  integer, dimension(numFeatures, numSamples), intent(in) :: data
  real(dp), dimension(numFeatures, numSamples), intent(in) :: bg
  logical, intent(in) :: lowerTail
  integer, intent(in) :: numBlocks
  integer, dimension(numBlocks), intent(in) :: blockLengths
  real(dp), dimension(numFeatures, numFeatures) :: pMatrix
  real(dp), dimension(numFeatures, numFeatures), intent(out) :: qMatrix
  real(dp), intent(out) :: pi0

  pMatrix = analyseBlockStructure_(data, bg, lowerTail, blockLengths, qMatrix, pi0)
end subroutine analyseBlockStructure


subroutine mutex(numFeatures, numSamples, data, bg, lowerTail, p, q, pi0)
  use types, only: dp
  use fdr, only: mutex_ => mutex

  integer, intent(in) :: numFeatures
  integer, intent(in) :: numSamples
  integer, dimension(numFeatures, numSamples), intent(in) :: data
  real(dp), dimension(numFeatures, numSamples), intent(in) :: bg
  logical, intent(in) :: lowerTail
  real(dp), dimension(numFeatures * (numFeatures - 1) / 2), intent(out) :: p
  real(dp), dimension(numFeatures * (numFeatures - 1) / 2), intent(out) :: q
  real(dp), intent(out) :: pi0

  p = mutex_(data, bg, lowerTail, q, pi0)
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
