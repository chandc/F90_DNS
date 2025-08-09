! =============================================================================
! OPTIMIZED FFT MODULE USING FFTW3
! =============================================================================
module fft_module_2d
    use iso_fortran_env, only: real64
    use iso_c_binding
    implicit none
    
    private
    integer, parameter :: wp = real64
    
    ! Public interfaces
    public :: fft_plans, setup_fft_plans, destroy_fft_plans
    public :: fft_forward_2d, fft_backward_2d
    public :: initialize_fft_threading, cleanup_fft_threading
    
    ! FFTW3 constants
    integer(c_int), parameter :: FFTW_FORWARD = -1
    integer(c_int), parameter :: FFTW_BACKWARD = 1
    integer(c_int), parameter :: FFTW_ESTIMATE = 64
    integer(c_int), parameter :: FFTW_MEASURE = 0
    
    ! FFTW3 interfaces
    interface
        ! Double precision real-to-complex 1D FFT
        function fftw_plan_dft_r2c_1d(n, in, out, flags) bind(C, name='fftw_plan_dft_r2c_1d')
            import :: c_ptr, c_int, c_double, c_double_complex
            integer(c_int), value :: n, flags
            real(c_double), intent(inout) :: in(*)
            complex(c_double_complex), intent(inout) :: out(*)
            type(c_ptr) :: fftw_plan_dft_r2c_1d
        end function
        
        ! Double precision complex-to-real 1D FFT
        function fftw_plan_dft_c2r_1d(n, in, out, flags) bind(C, name='fftw_plan_dft_c2r_1d')
            import :: c_ptr, c_int, c_double, c_double_complex
            integer(c_int), value :: n, flags
            complex(c_double_complex), intent(inout) :: in(*)
            real(c_double), intent(inout) :: out(*)
            type(c_ptr) :: fftw_plan_dft_c2r_1d
        end function
        
        ! Double precision real-to-complex FFT
        function fftw_plan_dft_r2c_2d(nx, ny, in, out, flags) bind(C, name='fftw_plan_dft_r2c_2d')
            import :: c_ptr, c_int, c_double, c_double_complex
            integer(c_int), value :: nx, ny, flags
            real(c_double), intent(inout) :: in(*)
            complex(c_double_complex), intent(inout) :: out(*)
            type(c_ptr) :: fftw_plan_dft_r2c_2d
        end function
        
        ! Double precision complex-to-real FFT
        function fftw_plan_dft_c2r_2d(nx, ny, in, out, flags) bind(C, name='fftw_plan_dft_c2r_2d')
            import :: c_ptr, c_int, c_double, c_double_complex
            integer(c_int), value :: nx, ny, flags
            complex(c_double_complex), intent(inout) :: in(*)
            real(c_double), intent(inout) :: out(*)
            type(c_ptr) :: fftw_plan_dft_c2r_2d
        end function
        
        ! Execute real-to-complex FFT
        subroutine fftw_execute_dft_r2c(plan, in, out) bind(C, name='fftw_execute_dft_r2c')
            import :: c_ptr, c_double, c_double_complex
            type(c_ptr), value :: plan
            real(c_double), intent(inout) :: in(*)
            complex(c_double_complex), intent(inout) :: out(*)
        end subroutine
        
        ! Execute complex-to-real FFT
        subroutine fftw_execute_dft_c2r(plan, in, out) bind(C, name='fftw_execute_dft_c2r')
            import :: c_ptr, c_double, c_double_complex
            type(c_ptr), value :: plan
            complex(c_double_complex), intent(inout) :: in(*)
            real(c_double), intent(inout) :: out(*)
        end subroutine
        
        ! Destroy plan
        subroutine fftw_destroy_plan(plan) bind(C, name='fftw_destroy_plan')
            import :: c_ptr
            type(c_ptr), value :: plan
        end subroutine
        
        ! Initialize threads (disabled)
        ! function fftw_init_threads() bind(C, name='fftw_init_threads')
        !     import :: c_int
        !     integer(c_int) :: fftw_init_threads
        ! end function
        
        ! Set number of threads (disabled)
        ! subroutine fftw_plan_with_nthreads(nthreads) bind(C, name='fftw_plan_with_nthreads')
        !     import :: c_int
        !     integer(c_int), value :: nthreads
        ! end subroutine
        
        ! Cleanup threads (disabled)
        ! subroutine fftw_cleanup_threads() bind(C, name='fftw_cleanup_threads')
        ! end subroutine
    end interface
    
    ! FFT plan storage
    type :: fft_plans
        type(c_ptr) :: forward_plan = c_null_ptr
        type(c_ptr) :: backward_plan = c_null_ptr
        integer :: nx, ny
        logical :: initialized = .false.
    end type fft_plans
    
contains

    ! =========================================================================
    ! INITIALIZE FFTW3 WITH THREADING SUPPORT
    ! =========================================================================
    subroutine initialize_fft_threading(num_threads)
        integer, intent(in), optional :: num_threads
        integer :: nthreads
        
        ! Set number of threads (for informational purposes only)
        if (present(num_threads)) then
            nthreads = num_threads
        else
            nthreads = 1  ! Single-threaded
        endif
        
        write(*,'(A,I0,A)') 'FFTW initialized (single-threaded), requested ', nthreads, ' threads'
        
    end subroutine initialize_fft_threading

    ! =========================================================================
    ! CLEANUP FFTW THREADING
    ! =========================================================================
    subroutine cleanup_fft_threading()
        ! No cleanup needed for single-threaded version
        write(*,*) 'FFTW cleanup completed'
    end subroutine cleanup_fft_threading    ! =========================================================================
    ! SETUP FFT PLANS FOR POISSON SOLVER
    ! =========================================================================
    subroutine setup_fft_plans(plans, nx, ny, sample_real, sample_complex)
        type(fft_plans), intent(inout) :: plans
        integer, intent(in) :: nx, ny
        real(wp), intent(inout), target :: sample_real(nx, ny)
        complex(wp), intent(inout), target :: sample_complex(nx/2+1, ny)
        
        ! Store dimensions
        plans%nx = nx
        plans%ny = ny
        
        ! Create forward plan (real to complex) - 1D FFT in x-direction only
        ! We'll apply this plan to each z-level separately
        plans%forward_plan = fftw_plan_dft_r2c_1d(nx, sample_real(:,1), sample_complex(:,1), FFTW_MEASURE)
        
        ! Create backward plan (complex to real) - 1D FFT in x-direction only  
        plans%backward_plan = fftw_plan_dft_c2r_1d(nx, sample_complex(:,1), sample_real(:,1), FFTW_MEASURE)
        
        plans%initialized = .true.
        write(*,'(A,I0,A,I0,A)') 'FFT plans created for ', nx, ' x ', ny, ' grid'
        
    end subroutine setup_fft_plans
    
    ! =========================================================================
    ! DESTROY FFT PLANS
    ! =========================================================================
    subroutine destroy_fft_plans(plans)
        type(fft_plans), intent(inout) :: plans
        
        if (plans%initialized) then
            if (c_associated(plans%forward_plan)) call fftw_destroy_plan(plans%forward_plan)
            if (c_associated(plans%backward_plan)) call fftw_destroy_plan(plans%backward_plan)
            plans%initialized = .false.
        endif
        
    end subroutine destroy_fft_plans
    
    ! =========================================================================
    ! FORWARD FFT (Real to Complex) - 1D FFTs in x-direction for each z-level
    ! =========================================================================
    subroutine fft_forward_2d(plans, real_data, complex_data)
        type(fft_plans), intent(in) :: plans
        real(wp), intent(inout), target :: real_data(:,:)
        complex(wp), intent(inout), target :: complex_data(:,:)
        integer :: k
        
        if (.not. plans%initialized) then
            write(*,*) 'Error: FFT plans not initialized'
            return
        endif
        
        ! Apply 1D FFT to each z-level separately
        do k = 1, plans%ny
            call fftw_execute_dft_r2c(plans%forward_plan, real_data(:,k), complex_data(:,k))
        end do
        
    end subroutine fft_forward_2d
    
    ! =========================================================================
    ! BACKWARD FFT (Complex to Real)
    ! =========================================================================
    ! BACKWARD FFT (Complex to Real) - 1D FFTs in x-direction for each z-level
    ! =========================================================================
    subroutine fft_backward_2d(plans, complex_data, real_data)
        type(fft_plans), intent(in) :: plans
        complex(wp), intent(inout), target :: complex_data(:,:)
        real(wp), intent(inout), target :: real_data(:,:)
        integer :: k
        real(wp) :: norm_factor
        
        if (.not. plans%initialized) then
            write(*,*) 'Error: FFT plans not initialized'
            return
        endif
        
        ! Apply 1D inverse FFT to each z-level separately
        do k = 1, plans%ny
            call fftw_execute_dft_c2r(plans%backward_plan, complex_data(:,k), real_data(:,k))
        end do
        
        ! Normalize (FFTW doesn't normalize by default) - only by nx since we're doing 1D FFTs
        norm_factor = 1.0_wp / plans%nx
        real_data = real_data * norm_factor
        
    end subroutine fft_backward_2d
    
    ! =========================================================================
    ! OPTIMIZED FFT FOR 4D ARRAYS (FOR POISSON SOLVER)
    ! =========================================================================
    subroutine fft_forward_4d(plans, real_4d, complex_4d)
        type(fft_plans), intent(in) :: plans
        real(wp), intent(inout) :: real_4d(:,:,:,:)
        complex(wp), intent(inout) :: complex_4d(:,:,:,:)
        
        integer :: k, ne, nk, nelem
        
        nk = size(real_4d, 3)
        nelem = size(real_4d, 4)
        
        !$omp parallel do private(ne, k) collapse(2)
        do ne = 1, nelem
            do k = 1, nk
                call fft_forward_2d(plans, real_4d(:,:,k,ne), complex_4d(:,:,k,ne))
            end do
        end do
        !$omp end parallel do
        
    end subroutine fft_forward_4d
    
    ! =========================================================================
    ! OPTIMIZED INVERSE FFT FOR 4D ARRAYS
    ! =========================================================================
    subroutine fft_backward_4d(plans, complex_4d, real_4d)
        type(fft_plans), intent(in) :: plans
        complex(wp), intent(inout) :: complex_4d(:,:,:,:)
        real(wp), intent(inout) :: real_4d(:,:,:,:)
        
        integer :: k, ne, nk, nelem
        
        nk = size(real_4d, 3)
        nelem = size(real_4d, 4)
        
        !$omp parallel do private(ne, k) collapse(2)
        do ne = 1, nelem
            do k = 1, nk
                call fft_backward_2d(plans, complex_4d(:,:,k,ne), real_4d(:,:,k,ne))
            end do
        end do
        !$omp end parallel do
        
    end subroutine fft_backward_4d
    
end module fft_module_2d
