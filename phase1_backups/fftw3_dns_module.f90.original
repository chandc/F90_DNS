!===============================================================================
! FFTW3_DNS_MODULE.F90 - Native FFTW3 Module for DNS Solver
!
! PURPOSE: Provides drop-in replacement for fft_module_2d using native FFTW3
!          2D real-to-complex transforms with proper dimension handling
!
! AUTHOR: GitHub Copilot  
! DATE: August 8, 2025
!===============================================================================

module fftw3_dns_module
    use iso_fortran_env, only: wp => real64
    use, intrinsic :: iso_c_binding
    implicit none
    
    include 'fftw3.f03'
    
    ! FFT plan storage
    type :: fftw3_dns_plans
        type(C_PTR) :: plan_forward   ! real(nx,ny) -> complex(nxhp,ny)
        type(C_PTR) :: plan_backward  ! complex(nxhp,ny) -> real(nx,ny)
        integer :: nx, ny, nxhp
        logical :: initialized = .false.
    end type fftw3_dns_plans
    
contains

    !---------------------------------------------------------------------------
    ! SUBROUTINE: setup_fftw3_plans_dns
    !
    ! Creates FFTW3 plans for the DNS solver FFT operations
    !---------------------------------------------------------------------------
    subroutine setup_fftw3_plans_dns(plans, nx, ny, sample_real, sample_complex)
        implicit none
        type(fftw3_dns_plans), intent(out) :: plans
        integer, intent(in) :: nx, ny
        real(wp), intent(inout) :: sample_real(nx, ny)
        complex(wp), intent(inout) :: sample_complex(nx/2+1, ny)
        
        plans%nx = nx
        plans%ny = ny
        plans%nxhp = nx/2 + 1
        
        write(*,'(A,I0,A,I0,A)') ' Setting up FFTW3 plans for ', nx, 'x', ny, ' grid'
        write(*,'(A,I0,A,I0)') ' Complex storage: ', plans%nxhp, 'x', ny
        
        ! Create forward plan: real(nx,ny) -> complex(nx/2+1,ny)
        ! Note: FFTW uses row-major ordering, so dimensions are swapped
        plans%plan_forward = fftw_plan_dft_r2c_2d(ny, nx, sample_real, sample_complex, FFTW_MEASURE)
        
        ! Create backward plan: complex(nx/2+1,ny) -> real(nx,ny)
        plans%plan_backward = fftw_plan_dft_c2r_2d(ny, nx, sample_complex, sample_real, FFTW_MEASURE)
        
        if (.not. c_associated(plans%plan_forward)) then
            stop 'Error: Failed to create FFTW3 forward plan'
        end if
        if (.not. c_associated(plans%plan_backward)) then
            stop 'Error: Failed to create FFTW3 backward plan'
        end if
        
        plans%initialized = .true.
        write(*,'(A)') ' FFTW3 plans created successfully'
        
    end subroutine setup_fftw3_plans_dns

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fftw3_forward_2d_dns
    !
    ! Performs 2D real-to-complex forward transform
    ! Replaces your current fft_forward_2d routine
    !---------------------------------------------------------------------------
    subroutine fftw3_forward_2d_dns(plans, real_field, complex_field)
        implicit none
        type(fftw3_dns_plans), intent(in) :: plans
        real(wp), intent(inout) :: real_field(plans%nx, plans%ny)
        complex(wp), intent(out) :: complex_field(plans%nxhp, plans%ny)
        
        if (.not. plans%initialized) then
            stop 'Error: FFTW3 plans not initialized'
        end if
        
        ! Execute the forward transform
        call fftw_execute_dft_r2c(plans%plan_forward, real_field, complex_field)
        
    end subroutine fftw3_forward_2d_dns

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fftw3_backward_2d_dns
    !
    ! Performs 2D complex-to-real backward transform  
    ! Replaces your current fft_backward_2d routine
    !---------------------------------------------------------------------------
    subroutine fftw3_backward_2d_dns(plans, complex_field, real_field)
        implicit none
        type(fftw3_dns_plans), intent(in) :: plans
        complex(wp), intent(inout) :: complex_field(plans%nxhp, plans%ny)
        real(wp), intent(out) :: real_field(plans%nx, plans%ny)
        
        if (.not. plans%initialized) then
            stop 'Error: FFTW3 plans not initialized'
        end if
        
        ! Execute the backward transform
        call fftw_execute_dft_c2r(plans%plan_backward, complex_field, real_field)
        
        ! CRITICAL: FFTW C2R transforms are unnormalized
        real_field = real_field / real(plans%nx * plans%ny, wp)
        
    end subroutine fftw3_backward_2d_dns

    !---------------------------------------------------------------------------
    ! SUBROUTINE: destroy_fftw3_plans_dns
    !
    ! Clean up FFTW3 plans and memory
    !---------------------------------------------------------------------------
    subroutine destroy_fftw3_plans_dns(plans)
        implicit none
        type(fftw3_dns_plans), intent(inout) :: plans
        
        if (plans%initialized) then
            call fftw_destroy_plan(plans%plan_forward)
            call fftw_destroy_plan(plans%plan_backward)
            call fftw_cleanup()
            plans%initialized = .false.
            write(*,'(A)') ' FFTW3 plans destroyed'
        end if
        
    end subroutine destroy_fftw3_plans_dns

end module fftw3_dns_module
