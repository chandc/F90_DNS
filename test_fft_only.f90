program test_fft_only
    use fft_module_2d
    implicit none
    
    integer, parameter :: nx = 4, ny = 4
    integer, parameter :: nkx = nx/2+1
    
    type(fft_plans) :: plans
    real(8) :: sample_real(nx, ny), data_real(nx, ny)
    complex(8) :: sample_complex(nkx, ny), data_complex(nkx, ny)
    
    write(*,*) 'Testing FFT module...'
    
    ! Initialize data
    data_real = 1.0d0
    
    call setup_fft_plans(plans, nx, ny, sample_real, sample_complex)
    call fft_forward_2d(plans, data_real, data_complex)
    call fft_backward_2d(plans, data_complex, data_real)
    call destroy_fft_plans(plans)
    
    write(*,*) 'FFT test completed! Result:', data_real(1,1)
    
end program test_fft_only
