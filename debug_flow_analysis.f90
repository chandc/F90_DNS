program debug_flow_analysis
    use lgl_module
    use iso_fortran_env, only: wp => real64
    implicit none
    
    ! Grid parameters
    integer, parameter :: nz = 33
    real(wp), parameter :: re = 100.0_wp
    real(wp), parameter :: ybar = 2.0_wp
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    
    ! Arrays
    real(wp), allocatable :: zpts(:), zwts(:)
    real(wp), allocatable :: d1(:,:), d2(:,:)
    real(wp), allocatable :: u_profile(:)
    real(wp), allocatable :: numerical_second_deriv(:)
    real(wp) :: forcing_x, second_deriv, analytical_forcing, z_scale
    integer :: k
    
    ! Allocate arrays
    allocate(zpts(nz), zwts(nz), d1(nz,nz), d2(nz,nz))
    allocate(u_profile(nz), numerical_second_deriv(nz))
    
    ! Setup LGL quadrature
    call lgl_nodes_weights(nz, zpts, zwts)
    call differentiation_matrix(nz, zpts, d1)
    d2 = matmul(d1, d1)
    
    ! Initialize Poiseuille profile
    do k = 1, nz
        u_profile(k) = 1.5_wp * (1.0_wp - zpts(k)**2)
    end do
    
    ! Calculate forcing as in the code
    forcing_x = 3.0_wp / re
    
    ! Calculate analytical second derivative of U(z) = 1.5*(1-z^2)
    ! d²U/dz² = d/dz(1.5*(-2z)) = 1.5*(-2) = -3.0
    analytical_forcing = -(-3.0_wp) / re  ! Required forcing = -ν * d²U/dz²
    
    write(*,'(A)') ' =========================================='
    write(*,'(A)') '   FLOW ANALYSIS DEBUG'
    write(*,'(A)') ' =========================================='
    write(*,'(A,F10.6)') ' Reynolds number: ', re
    write(*,'(A,F10.6)') ' Forcing used in code: ', forcing_x
    write(*,'(A,F10.6)') ' Analytical forcing needed: ', analytical_forcing
    write(*,'(A,F10.6)') ' Ratio (should be 1.0): ', forcing_x / analytical_forcing
    write(*,'(A)') ' '
    write(*,'(A)') ' Initial velocity profile (first/last 5 points):'
    do k = 1, min(5, nz)
        write(*,'(A,I2,A,F8.4,A,F10.6)') ' k=', k, ', z=', zpts(k), ', u=', u_profile(k)
    end do
    write(*,'(A)') ' ...'
    do k = max(nz-4, 1), nz
        write(*,'(A,I2,A,F8.4,A,F10.6)') ' k=', k, ', z=', zpts(k), ', u=', u_profile(k)
    end do
    write(*,'(A)') ' '
    
    ! Check numerical second derivative
    z_scale = 2.0_wp / ybar
    numerical_second_deriv = matmul(d2, u_profile) * z_scale**2
    
    write(*,'(A)') ' Numerical second derivative (middle points):'
    do k = nz/2-2, nz/2+2
        write(*,'(A,I2,A,F8.4,A,F10.6,A)') ' k=', k, ', z=', zpts(k), &
            ', d2u/dz2=', numerical_second_deriv(k), ' (analytical: -3.0)'
    end do
    
    write(*,'(A)') ' =========================================='
    
    deallocate(zpts, zwts, d1, d2, u_profile, numerical_second_deriv)
    
end program debug_flow_analysis
