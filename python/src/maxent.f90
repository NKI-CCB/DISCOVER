module maxent

  use types, only: dp

  implicit none

contains

  function outer(x, y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(size(x), size(y)) :: outer

    integer :: i, j

    do concurrent (i = 1:size(x), j = 1:size(y))
       outer(i, j) = x(i) * y(j)
    end do
  end function outer 

  subroutine fit(rowValues, rowWeights, colValues, colWeights, mu)
    integer, dimension(:), intent(in) :: rowValues
    integer, dimension(size(rowValues)), intent(in) :: rowWeights
    integer, dimension(:), intent(in) :: colValues
    integer, dimension(size(colValues)), intent(in) :: colWeights
    real(dp), dimension(size(rowValues) + size(colValues)), intent(out), target :: mu

    character(len=60) :: task, csave
    logical :: lsave(4)
    integer :: isave(44)
    real(dp) :: dsave(29)
    integer :: n, m
    real(dp) :: f
    real(dp), dimension(size(rowValues) + size(colValues)), target :: g
    real(dp), dimension(size(rowValues) + size(colValues)) :: l, u
    integer, dimension(size(rowValues) + size(colValues)) :: nbd
    real(dp), allocatable  :: wa(:)
    integer,  allocatable  :: iwa(:)

    real(dp), dimension(size(rowValues)) :: row_violations
    real(dp), dimension(size(colValues)) :: col_violations
    real(dp), dimension(size(rowValues), size(colValues)) :: P, eA
    integer :: nrows, ncols

    real(dp), dimension(:), pointer :: mu_rows
    real(dp), dimension(:), pointer :: mu_cols

    real(dp), dimension(:), pointer :: g_rows
    real(dp), dimension(:), pointer :: g_cols
    
    nrows = size(rowValues)
    ncols = size(colValues)

    mu_rows(1:nrows) => mu(:nrows)
    mu_cols(1:ncols) => mu(nrows+1:)

    g_rows(1:nrows) => g(:nrows)
    g_cols(1:ncols) => g(nrows+1:)

    mu = 0
    task = 'START'
    n = size(mu)
    m = 10
    nbd = 0

    allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))
    allocate(iwa(3*n))

    do while (task(1:2) .eq. 'FG' .or. task .eq. 'NEW_X' .or. task .eq. 'START')
       !call setulb(n, m, mu, l, u, nbd, f, g, 1.d+7, 1.d-5, wa, &
       call setulb(n, m, mu, l, u, nbd, f, g, 1.d+1, 1.d-5, wa, &
            iwa, task, 0, csave, lsave, &
            isave, dsave)

       if (task(1:2) .eq. 'FG') then
          eA = outer(exp(mu_rows), exp(mu_cols))
          P = 1.0d0 / (eA + 1)

          row_violations = rowWeights * (matmul(P, colWeights) - rowValues)
          col_violations = colWeights * (matmul(rowWeights, P) - colValues)
          f = -1 * (dot_product(mu_rows, row_violations) &
               + dot_product(mu_cols, col_violations) &
               + dot_product(rowWeights, matmul(log(eA) * eA / (eA + 1) - log(eA + 1), colWeights)))

          g_rows = -1 * row_violations
          g_cols = -1 * col_violations
       end if
    end do

    deallocate(wa)
    deallocate(iwa)
    
  end subroutine fit

end module maxent
