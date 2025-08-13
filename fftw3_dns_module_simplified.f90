! Simplified FFTW module for testing perturbation module compilation
module fftw3_dns_module_simplified
    use iso_fortran_env, only: real64
    implicit none
    
    integer, parameter :: wp = real64
    private :: wp
    
    type :: fftw3_dns_plans
        ! Placeholder type for compilation testing
        integer :: dummy
    end type fftw3_dns_plans
    
contains

subroutine fftw3_forward_2d_dns(plans, input_real, output_complex)
    implicit none
    type(fftw3_dns_plans), intent(in) :: plans
    real(wp), intent(in) :: input_real(:,:,:)
    complex(wp), intent(out) :: output_complex(:,:,:)
    
    ! Placeholder - just copy real part for testing
    output_complex = cmplx(input_real, 0.0_wp, wp)
    write(*,*) 'PLACEHOLDER: fftw3_forward_2d_dns called'
end subroutine fftw3_forward_2d_dns

subroutine fftw3_backward_2d_dns(plans, input_complex, output_real)
    implicit none
    type(fftw3_dns_plans), intent(in) :: plans
    complex(wp), intent(in) :: input_complex(:,:,:)
    real(wp), intent(out) :: output_real(:,:,:)
    
    ! Placeholder - just copy real part for testing
    output_real = real(input_complex, wp)
    write(*,*) 'PLACEHOLDER: fftw3_backward_2d_dns called'
end subroutine fftw3_backward_2d_dns

end module fftw3_dns_module_simplified
