module test_perturbation_module
    implicit none
    
    integer, parameter :: wp = selected_real_kind(15, 307)
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    
    ! Module-level variables for monitoring
    real(wp), save :: energy_prev = 0.0_wp, time_prev = 0.0_wp
    logical, save :: first_call = .true.
    
    private
    public :: test_basic_functionality

contains

subroutine test_basic_functionality()
    implicit none
    write(*,*) 'Perturbation module basic structure OK'
end subroutine test_basic_functionality

end module test_perturbation_module
