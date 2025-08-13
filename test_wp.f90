program test_wp
    use iso_fortran_env, only: wp => real64
    implicit none
    
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(wp) :: test_val
    
    test_val = pi * 2.0_wp
    write(*,*) 'wp test: ', test_val
    
end program test_wp
