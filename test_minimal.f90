module test_minimal
    use iso_fortran_env, only: wp => real64
    implicit none
    
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(wp), parameter :: small_val = 1.0e-12_wp
    
contains

subroutine test_sub()
    implicit none
    real(wp) :: x
    x = pi * small_val
    write(*,*) 'Test: ', x
end subroutine test_sub

end module test_minimal
