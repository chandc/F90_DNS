!===============================================================================
! FFTW3_2D_EXAMPLE.F90 - Native FFTW3 2D Real-to-Complex Transforms
!
! PURPOSE: Demonstrates proper usage of FFTW3 native 2D functions for
!          real-to-complex forward and complex-to-real backward transforms
!
! AUTHOR: GitHub Copilot
! DATE: August 8, 2025
!===============================================================================

program fftw3_2d_example
    use iso_fortran_env, only: wp => real64
    use, intrinsic :: iso_c_binding
    implicit none
    
    ! FFTW3 interface declarations
    include 'fftw3.f03'
    
    ! Grid parameters
    integer, parameter :: nx = 128, ny = 32
    integer, parameter :: nxhp = nx/2 + 1  ! Complex dimension after r2c transform
    
    ! Arrays for real and complex data
    real(wp), allocatable :: real_data(:,:)
    complex(wp), allocatable :: complex_data(:,:)
    real(wp), allocatable :: reconstructed_data(:,:)
    
    ! FFTW plans
    type(C_PTR) :: plan_forward, plan_backward
    
    ! Test parameters
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(wp), parameter :: Lx = 2.0_wp * pi, Ly = 2.0_wp * pi
    real(wp) :: x, y, error
    integer :: i, j
    
    write(*,'(A)') '============================================'
    write(*,'(A)') ' FFTW3 2D Real-to-Complex Transform Demo'
    write(*,'(A)') '============================================'
    write(*,'(A,I0,A,I0)') ' Grid: ', nx, ' x ', ny
    write(*,'(A,I0,A,I0)') ' Complex grid: ', nxhp, ' x ', ny
    
    ! Allocate arrays
    allocate(real_data(nx, ny))
    allocate(complex_data(nxhp, ny))
    allocate(reconstructed_data(nx, ny))
    
    ! Initialize test data with a known analytical function
    write(*,'(A)') ' Creating test data: f(x,y) = sin(2πx/Lx) * cos(3πy/Ly)'
    do j = 1, ny
        y = real(j-1, wp) * Ly / real(ny, wp)
        do i = 1, nx
            x = real(i-1, wp) * Lx / real(nx, wp)
            real_data(i,j) = sin(2.0_wp * pi * x / Lx) * cos(3.0_wp * pi * y / Ly)
        end do
    end do
    
    !---------------------------------------------------------------------------
    ! CREATE FFTW3 PLANS
    !---------------------------------------------------------------------------
    write(*,'(A)') ' Creating FFTW3 plans...'
    
    ! Forward plan: real(nx,ny) -> complex(nx/2+1,ny)
    ! FFTW_R2C: Real-to-complex transform
    ! FFTW_ESTIMATE: Quick plan creation (use FFTW_MEASURE for production)
    plan_forward = fftw_plan_dft_r2c_2d(ny, nx, real_data, complex_data, FFTW_ESTIMATE)
    
    ! Backward plan: complex(nx/2+1,ny) -> real(nx,ny)
    ! FFTW_C2R: Complex-to-real transform
    plan_backward = fftw_plan_dft_c2r_2d(ny, nx, complex_data, reconstructed_data, FFTW_ESTIMATE)
    
    if (.not. c_associated(plan_forward)) then
        stop 'Error: Failed to create forward FFTW plan'
    end if
    if (.not. c_associated(plan_backward)) then
        stop 'Error: Failed to create backward FFTW plan'
    end if
    
    write(*,'(A)') ' Plans created successfully!'
    
    !---------------------------------------------------------------------------
    ! FORWARD TRANSFORM: real -> complex
    !---------------------------------------------------------------------------
    write(*,'(A)') ' Executing forward transform (real -> complex)...'
    
    call fftw_execute_dft_r2c(plan_forward, real_data, complex_data)
    
    ! Display some spectral coefficients
    write(*,'(A)') ' Sample spectral coefficients:'
    write(*,'(A,2ES12.4)') '   F(0,0) = ', real(complex_data(1,1)), aimag(complex_data(1,1))
    write(*,'(A,2ES12.4)') '   F(1,0) = ', real(complex_data(2,1)), aimag(complex_data(2,1))
    write(*,'(A,2ES12.4)') '   F(0,1) = ', real(complex_data(1,2)), aimag(complex_data(1,2))
    
    !---------------------------------------------------------------------------
    ! BACKWARD TRANSFORM: complex -> real
    !---------------------------------------------------------------------------
    write(*,'(A)') ' Executing backward transform (complex -> real)...'
    
    call fftw_execute_dft_c2r(plan_backward, complex_data, reconstructed_data)
    
    ! IMPORTANT: FFTW C2R transforms are unnormalized
    ! Must divide by total number of points (nx * ny)
    reconstructed_data = reconstructed_data / real(nx * ny, wp)
    
    !---------------------------------------------------------------------------
    ! VERIFY ROUND-TRIP ACCURACY
    !---------------------------------------------------------------------------
    error = maxval(abs(real_data - reconstructed_data))
    write(*,'(A,ES12.4)') ' Round-trip error (max absolute): ', error
    
    if (error < 1.0e-12_wp) then
        write(*,'(A)') ' ✓ Round-trip test PASSED!'
    else
        write(*,'(A)') ' ✗ Round-trip test FAILED!'
    end if
    
    !---------------------------------------------------------------------------
    ! DEMONSTRATION: SPECTRAL DERIVATIVE
    !---------------------------------------------------------------------------
    write(*,'(A)') ''
    write(*,'(A)') ' Demonstrating spectral derivative calculation...'
    
    call demonstrate_spectral_derivative()
    
    !---------------------------------------------------------------------------
    ! CLEANUP
    !---------------------------------------------------------------------------
    write(*,'(A)') ' Cleaning up...'
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
    call fftw_cleanup()
    
    deallocate(real_data, complex_data, reconstructed_data)
    
    write(*,'(A)') ' Demo completed successfully!'
    write(*,'(A)') '============================================'

contains

    !---------------------------------------------------------------------------
    ! SUBROUTINE: demonstrate_spectral_derivative
    !
    ! Shows how to compute derivatives in spectral space using FFTW3
    !---------------------------------------------------------------------------
    subroutine demonstrate_spectral_derivative()
        implicit none
        
        real(wp), allocatable :: test_func(:,:), analytical_dx(:,:), numerical_dx(:,:)
        complex(wp), allocatable :: func_hat(:,:), dfdx_hat(:,:)
        real(wp), allocatable :: kx_values(:)
        type(C_PTR) :: plan_fwd, plan_bwd
        real(wp) :: kx, dx_error
        integer :: i, j
        
        write(*,'(A)') '   Computing ∂f/∂x where f(x,y) = sin(4πx/Lx) * cos(2πy/Ly)'
        
        ! Allocate arrays
        allocate(test_func(nx, ny), analytical_dx(nx, ny), numerical_dx(nx, ny))
        allocate(func_hat(nxhp, ny), dfdx_hat(nxhp, ny))
        allocate(kx_values(nxhp))
        
        ! Create test function and analytical derivative
        do j = 1, ny
            y = real(j-1, wp) * Ly / real(ny, wp)
            do i = 1, nx
                x = real(i-1, wp) * Lx / real(nx, wp)
                test_func(i,j) = sin(4.0_wp * pi * x / Lx) * cos(2.0_wp * pi * y / Ly)
                analytical_dx(i,j) = (4.0_wp * pi / Lx) * cos(4.0_wp * pi * x / Lx) * cos(2.0_wp * pi * y / Ly)
            end do
        end do
        
        ! Setup wavenumber array for x-direction
        do i = 1, nxhp
            kx_values(i) = 2.0_wp * pi * real(i-1, wp) / Lx
        end do
        
        ! Create plans for derivative calculation
        plan_fwd = fftw_plan_dft_r2c_2d(ny, nx, test_func, func_hat, FFTW_ESTIMATE)
        plan_bwd = fftw_plan_dft_c2r_2d(ny, nx, dfdx_hat, numerical_dx, FFTW_ESTIMATE)
        
        ! Forward transform
        call fftw_execute_dft_r2c(plan_fwd, test_func, func_hat)
        
        ! Compute derivative in spectral space: ∂f/∂x ↔ i*kx * f_hat
        do j = 1, ny
            do i = 1, nxhp
                kx = kx_values(i)
                dfdx_hat(i,j) = cmplx(0.0_wp, kx, kind=wp) * func_hat(i,j)
            end do
        end do
        
        ! Backward transform
        call fftw_execute_dft_c2r(plan_bwd, dfdx_hat, numerical_dx)
        
        ! Normalize (FFTW C2R is unnormalized)
        numerical_dx = numerical_dx / real(nx * ny, wp)
        
        ! Check accuracy
        dx_error = maxval(abs(analytical_dx - numerical_dx))
        write(*,'(A,ES12.4)') '   Spectral derivative error: ', dx_error
        
        if (dx_error < 1.0e-12_wp) then
            write(*,'(A)') '   ✓ Spectral derivative test PASSED!'
        else
            write(*,'(A)') '   ✗ Spectral derivative test FAILED!'
        end if
        
        ! Cleanup
        call fftw_destroy_plan(plan_fwd)
        call fftw_destroy_plan(plan_bwd)
        deallocate(test_func, analytical_dx, numerical_dx)
        deallocate(func_hat, dfdx_hat, kx_values)
        
    end subroutine demonstrate_spectral_derivative

end program fftw3_2d_example
