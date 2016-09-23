module maxent

  use types, only: dp

  implicit none

contains

  function scaleMatrix(matrix, rowScaling, colScaling) result (scaled)
    real(dp), dimension(:, :), intent(in) :: matrix
    integer, dimension(size(matrix, 1)), intent(in) :: rowScaling
    integer, dimension(size(matrix, 2)), intent(in) :: colScaling
    real(dp), dimension(size(matrix, 1), size(matrix, 2)) :: scaled

    integer :: i

    forall (i = 1:size(matrix, 1))
       scaled(i, :) = matrix(i, :) / colScaling
    end forall

    forall (i = 1:size(matrix, 2))
       scaled(:, i) = scaled(:, i) / rowScaling
    end forall
  end function scaleMatrix

  function multiplyMatrix(matrix, rowScaling, colScaling) result (scaled)
    real(dp), dimension(:, :), intent(in) :: matrix
    integer, dimension(size(matrix, 1)), intent(in) :: rowScaling
    integer, dimension(size(matrix, 2)), intent(in) :: colScaling
    real(dp), dimension(size(matrix, 1), size(matrix, 2)) :: scaled

    integer :: i

    forall (i = 1:size(matrix, 1))
       scaled(i, :) = matrix(i, :) * colScaling
    end forall

    forall (i = 1:size(matrix, 2))
       scaled(:, i) = scaled(:, i) * rowScaling
    end forall
  end function multiplyMatrix

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

    real(dp), dimension(size(rowValues), size(colValues)) :: P, eA
    integer :: nrows, ncols
    integer, dimension(2) :: rowShape, colShape

    real(dp), dimension(:, :), pointer :: mu_rows
    real(dp), dimension(:, :), pointer :: mu_cols

    real(dp), dimension(:), pointer :: g_rows
    real(dp), dimension(:), pointer :: g_cols
    
    nrows = size(rowValues)
    ncols = size(colValues)
    rowShape = [1, ncols]
    colShape = [nrows, 1]

    mu_rows(1:nrows, 1:1) => mu(:nrows)
    mu_cols(1:1, 1:ncols) => mu(nrows+1:)

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
          eA = matmul(exp(mu_rows / reshape(rowWeights, colShape)), &
                      exp(mu_cols / reshape(colWeights, rowShape)))
          P = 1.0d0 / (eA + 1)

          f = -1 * (dot_product(mu_rows(:, 1), matmul(P, colWeights) - rowValues) &
               + dot_product(mu_cols(1, :), matmul(rowWeights, P) - colValues) &
               + sum(multiplyMatrix(log(eA) * eA / (eA + 1) - log(eA + 1), rowWeights, colWeights)))

          g_rows = -1 * (matmul(P, colWeights) - rowValues)
          g_cols = -1 * (matmul(rowWeights, P) - colValues)
       end if       
    end do

    deallocate(wa)
    deallocate(iwa)
    
  end subroutine fit

end module maxent
