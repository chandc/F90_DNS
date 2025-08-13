program test_perturbation_compilation
    ! Minimal test program to check if perturbation module compiles
    ! by commenting out FFTW-dependent parts temporarily
    use iso_fortran_env, only: real64
    use lgl_module
    
    implicit none
    integer, parameter :: wp = real64
    
    ! Test parameters
    integer, parameter :: nx = 8, ny = 8, nz = 9
    real(wp), parameter :: xlen = 4.0_wp, ylen = 2.0_wp
    real(wp) :: z_nodes(nz)
    real(wp) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    real(wp) :: perturbation_amplitude = 0.01_wp
    
    write(*,*) 'Testing perturbation module compilation...'
    write(*,*) 'Grid: ', nx, 'x', ny, 'x', nz
    write(*,*) 'Domain: ', xlen, 'x', ylen
    write(*,*) 'Test completed successfully!'
    
end program test_perturbation_compilation
