program test_spectral_integration
    use iso_fortran_env, only: real64
    use lgl_module
    implicit none
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: real64
            integer, intent(in) :: n, nrhs, lda, ldb
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    
    integer, parameter :: nz = 9
    real(real64) :: z_nodes(nz), lgl_deriv_matrix(nz,nz)
    real(real64) :: rhs(nz), solution(nz), exact(nz)
    real(real64) :: matrix(nz,nz), temp_matrix(nz,nz)
    integer :: ipiv(nz), info
    integer :: k
    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    real(real64) :: max_error
    
    write(*,'(A)') '🧪 Testing Spectral Integration with LAPACK'
    write(*,'(A)') '==========================================='
    
    ! Generate LGL nodes and differentiation matrix
    call legendre_gauss_lobatto(nz, z_nodes)
    call differentiation_matrix(nz, z_nodes, lgl_deriv_matrix)
    
    ! Test problem: ∂w/∂z = -sin(πz) with w(±1) = 0
    ! Exact solution: w(z) = (cos(πz) - cos(π))/π = (cos(πz) + 1)/π
    
    ! Set up system matrix (differentiation matrix with boundary conditions)
    matrix = lgl_deriv_matrix
    matrix(1,:) = 0.0_real64; matrix(1,1) = 1.0_real64    ! w(-1) = 0
    matrix(nz,:) = 0.0_real64; matrix(nz,nz) = 1.0_real64  ! w(+1) = 0
    
    ! Set up right-hand side
    do k = 1, nz
        rhs(k) = -sin(pi * z_nodes(k))  ! ∂w/∂z = -sin(πz)
        exact(k) = (cos(pi * z_nodes(k)) + 1.0_real64) / pi
    end do
    rhs(1) = 0.0_real64   ! Boundary condition w(-1) = 0
    rhs(nz) = 0.0_real64  ! Boundary condition w(+1) = 0
    
    ! Solve the system
    temp_matrix = matrix
    solution = rhs
    
    call dgesv(nz, 1, temp_matrix, nz, ipiv, solution, nz, info)
    
    if (info == 0) then
        write(*,'(A)') '✅ LAPACK solver succeeded'
        
        max_error = maxval(abs(solution - exact))
        write(*,'(A,ES12.5)') 'Maximum error: ', max_error
        
        if (max_error < 1.0e-12_real64) then
            write(*,'(A)') '✅ Spectral integration HIGHLY ACCURATE'
        else if (max_error < 1.0e-8_real64) then
            write(*,'(A)') '✅ Spectral integration GOOD'
        else
            write(*,'(A)') '❌ Spectral integration NEEDS IMPROVEMENT'
        end if
        
        write(*,'(A)') 'Solution comparison:'
        do k = 1, nz
            write(*,'(A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                '  z(', k, '): computed = ', solution(k), ', exact = ', exact(k), &
                ', error = ', abs(solution(k) - exact(k))
        end do
    else
        write(*,'(A,I0)') '❌ LAPACK solver failed with info = ', info
    end if

end program test_spectral_integration
