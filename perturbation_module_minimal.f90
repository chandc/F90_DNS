module perturbation_module_minimal
    use iso_fortran_env, only: wp => real64
    use lgl_module         
    use fftw3_dns_module_simplified   
    implicit none
    
    private
    public :: test_minimal_perturbation

contains

subroutine test_minimal_perturbation(nx, ny, nz)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(wp) :: test_val
    
    test_val = pi * real(nx, wp)
    write(*,*) 'Minimal test: ', test_val, nx, ny, nz
    
end subroutine test_minimal_perturbation

end module perturbation_module_minimal
