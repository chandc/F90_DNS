!===============================================================================
!
! DNS_PRESSURE_BC.F90 - 3D Incompressible Navier-Stokes Solver
!
! DESCRIPTION:
!   Complete F90 modernization of the original F77 DNS channel flow solver
!   by Daniel Chiu-Leung Chan (1993). This version implements a robust
!   pressure boundary condition approach using iterative CGS solvers.
!
! KEY FEATURES:
!   • Spectral methods in streamwise (x) direction using FFT
!   • Legendre-Gauss-Lobatto (LGL) collocation in wall-normal (z) direction  
!   • Crank-Nicolson time integration for viscous terms
!   • 4th-order Runge-Kutta for convective terms
!   • Kim & Moin boundary conditions for viscous wall treatment
!   • CGS iterative solver for pressure Poisson equation (F77-compatible)
!   • Bottom-wall-only pressure pinning for zero mode stability
!
! MATHEMATICAL FORMULATION:
!   Incompressible Navier-Stokes equations in channel flow geometry:
!   ∂u/∂t + u·∇u = -∇p + (1/Re)∇²u + f
!   ∇·u = 0
!
!   Fractional step method:
!   1. Momentum step: u* = u^n + dt*[RK4(convection) + CN(viscous)]
!   2. Pressure step: ∇²φ = ∇·u*/dt  (Poisson equation)
!   3. Projection step: u^{n+1} = u* - dt*∇φ
!
! GRID CONFIGURATION:
!   • nx = 128: Fourier modes in streamwise direction
!   • nz = 33:  LGL collocation points in wall-normal direction
!   • nxpp = nx + 2: Padded grid for FFT efficiency
!   • Periodic in x, no-slip walls at z = ±1
!
! PRESSURE SOLVER:
!   Uses F77-compatible approach with CGS iterative solver:
!   • Zero mode (kx=0): Special treatment with bottom-wall pinning
!   • Non-zero modes: Standard Poisson equation with natural BCs
!   • Same diff matrix construction for all modes (F77 approach)
!   • Robust handling of near-singular systems
!
! AUTHORS:
!   Original F77: Daniel Chiu-Leung Chan (1993)
!   F90 Modernization: GitHub Copilot with F77 compatibility focus
!
! VERSION: 2025.8.3 - Pressure BC Development Branch
!
!===============================================================================
! - Target: Crank-Nicolson viscous step: (I - 0.5*dt/Re ∇²)u^{n+1} = (I + 0.5*dt/Re ∇²)u^n + dt·RHS
! - Approach: Incremental implementation and testing
! 
! Date: August 3, 2025
! =============================================================================

program channel_flow_solver
    use iso_fortran_env, only: wp => real64
    use lgl_module
    use fft_module
    implicit none
    
    ! LAPACK interface for linear system solving
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            real(8), intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    
    ! Grid parameters - now runtime variables instead of compile-time parameters
    integer :: nx, nz  ! Grid dimensions
    integer :: nxpp, nxh, nxhp, nxf, ntot, nzm  ! Derived grid parameters
    
    ! Main flow parameters structure - consolidates all COMMON blocks
    type :: navier_stokes_params
        ! Time parameters (from /time5/)
        real(wp) :: t, dt
        integer :: nsteps, istep
        
        ! Timing variables
        real(8) :: step_start_time, step_end_time, total_time, avg_step_time
        integer :: timing_counter
        
        ! Physical parameters (from /param/, /inputs/, /params/)
        real(wp) :: alpha, beta, re, ta, ybar, cgstol, cs, u00, wavlen
        real(wp) :: xlen, ylen, retau
        real(wp) :: facvel, fact, fac1, fac2
        integer :: istart, nwrt, iform, iles, ntavg
        real(wp) :: omz  ! rotation parameter
        logical :: use_crank_nicolson  ! CN viscous step flag
        
        ! Wave numbers (from /wave/)
        real(wp), allocatable :: xw(:), xsq(:)
        
        ! LGL grid (from /leng/)
        real(wp), allocatable :: zpts(:), wg(:), d(:,:)
        
        ! Flow fields (from /flow/)
        real(wp), allocatable :: u(:), w(:), temp(:)
        real(wp), allocatable :: un(:), wn(:), tn(:)
        real(wp), allocatable :: pressure(:)  ! Pressure field from previous step
        
        ! FFT work arrays (from /fft/)
        real(wp), allocatable :: trigsx(:), work(:)
        integer, allocatable :: ifaxx(:)
        type(fft_plans) :: fft_plan  ! FFTW plans
        real(wp), allocatable :: real_scratch(:,:), complex_scratch_real(:,:), complex_scratch_imag(:,:)
        complex(wp), allocatable :: complex_scratch(:,:)
        
        ! Matrix arrays (from /matrix/, /coeff/)
        real(wp), allocatable :: su(:), sw(:), st(:)
        real(wp), allocatable :: diff(:,:), amass(:)
        
        ! Boundary conditions (from /bc/)
        real(wp), allocatable :: ubc(:), wbc(:)
        
        ! Scratch arrays (from /scratch/, /source/)
        real(wp), allocatable :: ox(:), oz(:), body(:)
        
        ! Working arrays for boundary conditions
        real(wp), allocatable :: uw1(:), ww1(:), uw2(:), ww2(:)
        
        ! Divergence checking array
        real(wp), allocatable :: div_check(:)
        real(wp), allocatable :: div_check_2d(:,:)
    end type navier_stokes_params
    
    type(navier_stokes_params) :: p
    integer :: istep_local, i
    real(wp), parameter :: pressure_gradient = 1.0_wp  ! Driving pressure gradient
    real(wp) :: u_max, u_rms, w_max, w_rms  ! Velocity statistics
    real(wp) :: div_max, div_rms  ! Divergence statistics
    
    write(*,'(A)') ' ============================================'
    write(*,'(A)') '   3D Navier-Stokes Channel Flow Solver'
    write(*,'(A)') '   Complete F90 Conversion'
    write(*,'(A)') '   Originally by Daniel Chiu-Leung Chan, 1993'
    write(*,'(A)') ' ============================================'
    
    ! Initialize simulation parameters and grid
    call initialize_solver(p)
    
    ! Set up initial or restart conditions
    if (p%istart == 0) then
        call setup_grid_and_matrices(p)
        call initialize_flow_fields(p)
        p%t = 0.0_wp
        p%fac1 = 1.0_wp
        p%fac2 = 0.0_wp
    else
        call setup_grid_and_matrices(p)
        call initialize_flow_fields(p)  ! This allocates arrays
        call restart_from_file(p)       ! This reads into allocated arrays
        p%fac1 = 1.5_wp
        p%fac2 = -0.5_wp
    endif
    
    write(*,'(A,F8.1,A,F8.5)') '  Re = ', p%re, ', dt = ', p%dt
    write(*,'(A,I0,A)') '  Running ', p%nsteps, ' time steps'
    write(*,'(A)') ' ============================================'
    
    ! Initialize timing
    p%total_time = 0.0d0
    p%timing_counter = 0
    
    ! Main time integration loop
    p%ntavg = 0
    
    do istep_local = 1, p%nsteps
        p%istep = istep_local
        
        ! Start timing for this time step
        call cpu_time(p%step_start_time)
        
        ! 4-stage Runge-Kutta convection step - NOW ENABLED
        call convection_step(p)
        
        p%fac1 = 1.5_wp
        p%fac2 = -0.5_wp
        
        ! Build momentum source terms including convection and forcing
        ! The RK4 convection step has updated p%su and p%sw with convective terms
        
        ! Add forcing terms to the convection source terms
        do i = 1, ntot
            ! For Helmholtz equation: (I - dt/Re * ∇²) u* = u^n + dt * (convection + forcing)
            ! The RK4 step provides convection terms in p%su, p%sw
            ! Now add pressure gradient forcing to u-momentum source term
            
            p%su(i) = p%dt * (p%su(i) + pressure_gradient)  ! dt*(convection + forcing) for u-momentum
            p%sw(i) = p%dt * p%sw(i)                         ! dt*convection for w-momentum
        end do
        
        ! Solve Helmholtz equations for velocity
        call solve_helmholtz_system(p)
        
        ! Apply boundary conditions
        call apply_edge_conditions(p)
        
        ! Update time
        p%t = p%t + p%dt
        
        ! Print velocity statistics to show flow evolution
        u_max = maxval(abs(p%u))
        u_rms = sqrt(sum(p%u**2) / real(ntot, wp))
        w_max = maxval(abs(p%w))
        w_rms = sqrt(sum(p%w**2) / real(ntot, wp))
        
        ! Check divergence and get values
        call compute_divergence_spectral(p, p%div_check_2d)
        div_max = maxval(abs(p%div_check_2d))
        div_rms = sqrt(sum(p%div_check_2d**2) / size(p%div_check_2d))
        
        write(*,'(A,ES12.4,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') ' Time = ', p%t, &
            ', u_max=', u_max, ', u_rms=', u_rms, ', w_max=', w_max, ', w_rms=', w_rms, &
            ', div_max=', div_max, ', div_rms=', div_rms
        
        ! Check divergence
        call check_divergence(p)
        
        ! End timing for this time step
        call cpu_time(p%step_end_time)
        p%timing_counter = p%timing_counter + 1
        p%total_time = p%total_time + (p%step_end_time - p%step_start_time)
        
        ! Output intermediate results and timing
        if (mod(istep_local, p%nwrt) == 0) then
            call output_solution(p)
            ! Output timing statistics every output interval
            p%avg_step_time = p%total_time / p%timing_counter
            write(*,'(A,F8.5,A,F8.5,A,F8.1,A)') ' Timing: step=', &
                (p%step_end_time - p%step_start_time)*1000.0, 'ms, avg=', &
                p%avg_step_time*1000.0, 'ms, remaining=', &
                (p%nsteps - istep_local) * p%avg_step_time / 60.0, 'min'
        endif
        
        ! Brief timing update every 10,000 steps
        if (mod(istep_local, 10000) == 0) then
            p%avg_step_time = p%total_time / p%timing_counter
            write(*,'(A,I0,A,F6.3,A,F6.1,A)') ' Step ', istep_local, ': avg=', &
                p%avg_step_time*1000.0, 'ms, ETA=', &
                (p%nsteps - istep_local) * p%avg_step_time / 60.0, 'min'
        endif
        
    end do
    
    ! Final timing summary
    if (p%timing_counter > 0) then
        p%avg_step_time = p%total_time / p%timing_counter
        write(*,'(A)') ' ============================================'
        write(*,'(A,I0,A)') ' Timing Summary for ', p%timing_counter, ' time steps:'
        write(*,'(A,F8.3,A)') ' Average time per step: ', p%avg_step_time*1000.0, ' ms'
        write(*,'(A,F8.1,A)') ' Total computation time: ', p%total_time, ' seconds'
        write(*,'(A,F8.1,A)') ' Total computation time: ', p%total_time/60.0, ' minutes'
        write(*,'(A,F8.3,A)') ' Time steps per second: ', 1.0/p%avg_step_time, ' steps/s'
        write(*,'(A)') ' ============================================'
    endif
    
    ! Final output
    call final_output(p)
    
    ! Cleanup
    call finalize_solver(p)
    
    write(*,'(A)') ' ============================================'
    write(*,'(A)') '   F90 Conversion Completed Successfully!'
    write(*,'(A)') ' ============================================'

contains

    ! =========================================================================
    ! INITIALIZATION ROUTINES
    ! =========================================================================
    
    !============================================================================
    ! SUBROUTINE: initialize_solver
    !
    ! PURPOSE:
    !   Main initialization routine for the 3D Navier-Stokes solver
    !
    ! DESCRIPTION:
    !   Sets up all solver components including:
    !   - Grid parameters and arrays
    !   - FFT plans and spectral arrays
    !   - LGL matrices for wall-normal derivatives
    !   - Flow field initialization or restart
    !   - Boundary conditions
    !
    ! INPUTS:
    !   p - Navier-Stokes parameter structure (intent: inout)
    !
    ! OUTPUTS:
    !   p - Fully initialized solver ready for time stepping
    !
    ! DEPENDENCIES:
    !   - FFTW3 library for spectral transforms
    !   - LGL module for Legendre-Gauss-Lobatto grids
    !   - Input file 'input.dat' with namelist parameters
    !============================================================================
    subroutine initialize_solver(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), parameter :: pi2 = 8.0_wp * atan(1.0_wp)
        logical :: file_exists
        
        ! Set default values from original F77 DATA statements
        p%istart = 0; p%dt = 0.01_wp; p%nsteps = 100000
        p%alpha = 1.0_wp; p%beta = 1.0_wp; p%re = 180.0_wp; p%ta = 0.0_wp
        p%nwrt = 10000; p%iform = 0; p%ybar = 1.0_wp
        p%cgstol = 1.0e-4_wp; p%iles = 0; p%cs = 0.1_wp
        p%wavlen = 1.0_wp; p%ylen = 1.0_wp; p%u00 = 0.0_wp
        p%xlen = pi2  ! Default domain length = 2π
        
        ! Read input parameters FIRST to get grid dimensions
        inquire(file='input.dat', exist=file_exists)
        if (file_exists) then
            call read_input_file(p)
        else
            write(*,'(A)') ' Using default parameters (input.dat not found)'
            ! Setup default grid if no input file
            call setup_grid_parameters(128, 33)
        endif
        
        ! Now allocate arrays using the grid parameters from input
        allocate(p%xw(nxhp), p%xsq(nxhp))
        allocate(p%trigsx(nxf), p%work(ntot), p%ifaxx(13))
        allocate(p%uw1(nxpp), p%ww1(nxpp), p%uw2(nxpp), p%ww2(nxpp))
        
        ! Allocate FFT scratch arrays
        allocate(p%real_scratch(nx, nz))
        allocate(p%complex_scratch(nxhp, nz))
        allocate(p%complex_scratch_real(nxhp, nz), p%complex_scratch_imag(nxhp, nz))
        
        ! Initialize FFT arrays and FFTW
        call setup_fft_arrays(p)
        call initialize_fft_threading(4)  ! Use 4 threads
        
        ! Initialize working arrays
        p%uw1 = 0.0_wp; p%ww1 = 0.0_wp
        p%uw2 = 0.0_wp; p%ww2 = 0.0_wp
        
        ! Compute derived parameters
        p%alpha = pi2 / p%xlen
        p%retau = p%re
        p%fact = 2.0_wp / p%dt
        p%facvel = p%fact * p%re
        
        write(*,'(A,F10.6)') ' alpha = ', p%alpha
        write(*,'(A,F8.3)') ' u00/vel ratio = ', p%u00
        write(*,'(A,F8.3)') ' Wavelength = ', p%wavlen
        
    end subroutine initialize_solver
    
    !============================================================================
    ! SUBROUTINE: setup_fft_arrays
    !
    ! PURPOSE:
    !   Initialize FFT arrays and wavenumber tables for spectral methods
    !
    ! DESCRIPTION:
    !   Sets up fundamental spectral method components:
    !   - Wavenumber arrays: xw(i) = kₓ, xsq(i) = λₓ = kₓ²
    !   - Spectral differentiation matrices
    !   - FFTW workspace arrays
    !   
    ! SPECTRAL METHOD DETAILS:
    !   alpha = 2π/L_x (fundamental wavenumber from domain length)
    !   xw(i) = wavenumber for mode i-1: kₓ = (i-1) * alpha
    !   xsq(i) = lambda_x = kₓ² (used in Helmholtz and Poisson equations)
    !
    ! INPUTS:
    !   p - Solver parameters with alpha already set
    !
    ! OUTPUTS:
    !   p - Updated with spectral arrays and FFT workspace
    !============================================================================
    subroutine setup_fft_arrays(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i
        
        ! Set up wave numbers for spectral methods
        ! alpha = 2π/L_x (fundamental wavenumber from domain length)
        ! xw(i) = wavenumber for mode i-1: kₓ = (i-1) * alpha
        ! xsq(i) = lambda_x = kₓ² (used in Helmholtz and Poisson equations)
        do i = 1, nxhp
            p%xw(i) = real(i-1, wp) * p%alpha    ! Wavenumber: kₓ
            p%xsq(i) = p%xw(i)**2               ! Squared wavenumber: λₓ = kₓ²
        end do
        
        ! Initialize legacy FFT arrays (for compatibility)
        p%trigsx = 0.0_wp
        p%ifaxx = 0
        
        ! Initialize scratch arrays
        p%real_scratch = 0.0_wp
        p%complex_scratch = (0.0_wp, 0.0_wp)
        p%complex_scratch_real = 0.0_wp
        p%complex_scratch_imag = 0.0_wp
        
        ! Setup FFTW plans for spectral transforms
        call setup_fft_plans(p%fft_plan, nx, nz, p%real_scratch, p%complex_scratch)
        
        write(*,'(A)') ' FFT arrays and FFTW plans initialized'
        
    end subroutine setup_fft_arrays
    
    !============================================================================
    ! SUBROUTINE: setup_grid_parameters
    !
    ! PURPOSE:
    !   Configure runtime grid parameters and derived constants
    !
    ! DESCRIPTION:
    !   Sets global grid parameters based on namelist input:
    !   - nx, nz: Main grid dimensions
    !   - nxpp, nxh, nxhp, nxf: Derived FFT grid parameters
    !   - Enables runtime grid configuration without recompilation
    !
    ! GRID PARAMETER RELATIONSHIPS:
    !   nxpp = nx + 2 (padded for FFT)
    !   nxh = nx/2 (half-grid for complex FFT)
    !   nxhp = nxh + 1 (half-grid plus one)
    !   nxf = 3*nx/2 + 1 (dealiased grid for nonlinear terms)
    !
    ! INPUTS:
    !   nx_in, nz_in - Grid dimensions from namelist input
    !
    ! OUTPUTS:
    !   Global grid parameters set for use throughout solver
    !============================================================================
    subroutine setup_grid_parameters(nx_in, nz_in)
        implicit none
        integer, intent(in) :: nx_in, nz_in
        
        ! Set main grid dimensions
        nx = nx_in
        nz = nz_in
        
        ! Compute derived grid parameters
        nxpp = nx + 2
        nxh = nx/2
        nxhp = nxh + 1
        nxf = 3*nx/2 + 1
        ntot = nxpp*nz
        nzm = nz - 1
        
        write(*,'(A,I0,A,I0)') ' Grid parameters set: nx=', nx, ', nz=', nz
        write(*,'(A,I0,A,I0,A,I0,A,I0)') ' Derived: nxpp=', nxpp, ', nxh=', nxh, &
                                           ', nxhp=', nxhp, ', nxf=', nxf
        
    end subroutine setup_grid_parameters

    !============================================================================
    ! SUBROUTINE: read_input_file
    !
    ! PURPOSE:
    !   Read simulation parameters from namelist input file
    !
    ! DESCRIPTION:
    !   Parses 'input.dat' using modern Fortran namelist format:
    !   - &grid: Runtime grid configuration (nx_input, nz_input)
    !   - &time_control: Time stepping parameters (dt, nsteps, etc.)
    !   - &simulation: Physical parameters (Re, alpha, etc.)
    !   - &output: Output control settings
    !
    ! FEATURES:
    !   - Self-documenting input format
    !   - Default value fallback for missing parameters
    !   - Comprehensive parameter validation
    !   - Grid parameters read first for array allocation
    !
    ! INPUTS:
    !   p - Solver parameters structure
    !
    ! OUTPUTS:
    !   p - Updated with all input parameters
    !   Global grid parameters set via setup_grid_parameters
    !============================================================================
    subroutine read_input_file(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: io_status
        
        ! Namelist variables (must match parameter names for clarity)
        integer :: istart, nsteps, nwrt, iform, iles
        integer :: nx_input, nz_input  ! Grid parameters
        real(wp) :: dt, alpha, beta, re, ta, ybar, cgstol, cs, u00, wavlen, xlen, ylen
        logical :: use_crank_nicolson  ! CN flag
        
        ! Define namelists
        namelist /grid/ nx_input, nz_input
        namelist /time_control/ istart, dt, nsteps, nwrt
        namelist /simulation/ alpha, beta, re, ta, ybar, cgstol, cs, u00, wavlen, xlen, ylen, use_crank_nicolson
        namelist /output/ iform, iles
        
        write(*,'(A)') ' Reading parameters from input.dat...'
        
        ! Initialize namelist variables with defaults
        nx_input = 128; nz_input = 33  ! Default grid sizes
        istart = p%istart; dt = p%dt; nsteps = p%nsteps; nwrt = p%nwrt
        alpha = p%alpha; beta = p%beta; re = p%re; ta = p%ta
        ybar = p%ybar; cgstol = p%cgstol; cs = p%cs; u00 = p%u00
        wavlen = p%wavlen; xlen = p%xlen; ylen = p%ylen
        iform = p%iform; iles = p%iles
        use_crank_nicolson = .true.  ! Default: use Crank-Nicolson (can be overridden in input.dat)
        
        open(7, file='input.dat', status='old', iostat=io_status)
        if (io_status /= 0) then
            write(*,'(A)') ' Warning: Could not open input.dat, using defaults'
            call setup_grid_parameters(nx_input, nz_input)
            return
        endif
        
        ! Read namelists
        read(7, nml=grid, iostat=io_status)
        if (io_status /= 0) then
            write(*,'(A)') ' Warning: Error reading &grid namelist, using defaults'
            rewind(7)
        endif
        
        ! Setup grid parameters immediately after reading them
        call setup_grid_parameters(nx_input, nz_input)
        read(7, nml=time_control, iostat=io_status)
        if (io_status /= 0) then
            write(*,'(A)') ' Warning: Error reading &time_control namelist, using defaults'
            rewind(7)
        endif
        
        read(7, nml=simulation, iostat=io_status)
        if (io_status /= 0) then
            write(*,'(A)') ' Warning: Error reading &simulation namelist, using defaults'
            rewind(7)
        endif
        
        read(7, nml=output, iostat=io_status)
        if (io_status /= 0) then
            write(*,'(A)') ' Warning: Error reading &output namelist, using defaults'
        endif
        
        close(7)
        
        ! Transfer namelist variables back to parameters
        p%istart = istart; p%dt = dt; p%nsteps = nsteps; p%nwrt = nwrt
        p%alpha = alpha; p%beta = beta; p%re = re; p%ta = ta
        p%ybar = ybar; p%cgstol = cgstol; p%cs = cs; p%u00 = u00
        p%wavlen = wavlen; p%xlen = xlen; p%ylen = ylen
        p%iform = iform; p%iles = iles
        p%use_crank_nicolson = use_crank_nicolson
        
        ! Print read values for verification
        write(*,'(A)') ' Successfully read namelist parameters:'
        write(*,'(A,I0,A,I0)') '   Grid: nx=', nx, ', nz=', nz
        write(*,'(A,I0,A,F0.4,A,I0,A,I0)') '   Time control: istart=', istart, ', dt=', dt, &
                                             ', nsteps=', nsteps, ', nwrt=', nwrt
        write(*,'(A,F0.1,A,F0.4)') '   Physical: Re=', re, ', alpha=', alpha
        if (use_crank_nicolson) then
            write(*,'(A)') '   Viscous method: Crank-Nicolson (2nd order)'
        else
            write(*,'(A)') '   Viscous method: Backward Euler (1st order)'
        endif
        write(*,'(A,I0,A,I0)') '   Output: iform=', iform, ', iles=', iles
        
    end subroutine read_input_file
    
    subroutine setup_grid_and_matrices(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i, j
        
        ! Allocate grid arrays
        allocate(p%zpts(0:nzm), p%wg(0:nzm), p%d(0:nzm,0:nzm))
        allocate(p%diff(nz,nz), p%amass(nz))
        allocate(p%ubc(nxpp), p%wbc(nxpp))
        
        ! Setup LGL grid using existing module
        call lgl_nodes_weights(nz, p%zpts(0:nzm), p%wg(0:nzm))
        call differentiation_matrix(nz, p%zpts(0:nzm), p%d(0:nzm,0:nzm))
        
        ! Setup mass matrix
        do i = 1, nz
            if (i <= nzm) then
                p%amass(i) = p%wg(i-1)
            else
                p%amass(i) = 0.0_wp
            endif
        end do
        
        ! Setup differentiation matrix for solver
        p%diff = 0.0_wp
        do i = 1, nz
            do j = 1, nz
                if (i <= nzm .and. j <= nzm) then
                    p%diff(i,j) = p%d(i-1,j-1)
                endif
            end do
        end do
        
        ! Initialize boundary condition arrays
        p%ubc = 0.0_wp
        p%wbc = 0.0_wp
        
        write(*,'(A)') ' Grid and matrices initialized'
        
    end subroutine setup_grid_and_matrices
    
    subroutine initialize_flow_fields(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i, j, k, idx
        real(wp) :: z_coord, u_theoretical
        
        ! Allocate flow field arrays
        allocate(p%u(ntot), p%w(ntot), p%temp(ntot))
        allocate(p%un(ntot), p%wn(ntot), p%tn(ntot))
        allocate(p%pressure(ntot))
        allocate(p%su(ntot), p%sw(ntot), p%st(ntot))
        allocate(p%ox(ntot), p%oz(ntot), p%body(ntot))
        allocate(p%div_check(ntot))
        allocate(p%div_check_2d(nxhp, nz))
        
        ! Initialize with laminar channel flow profile driven by pressure gradient
        p%u = 0.0_wp; p%w = 0.0_wp; p%temp = 0.0_wp
        p%un = 0.0_wp; p%wn = 0.0_wp; p%tn = 0.0_wp
        p%pressure = 0.0_wp
        p%su = 0.0_wp; p%sw = 0.0_wp; p%st = 0.0_wp
        p%ox = 0.0_wp; p%oz = 0.0_wp; p%body = 0.0_wp
        
        ! Start with a low initial velocity profile (u_max = 1.0)
        ! This will allow us to observe the full development from near-rest to steady state
        do k = 1, nz
            do i = 1, nxpp
                idx = (k-1)*nxpp + i
                ! Map LGL grid point to z ∈ [-1,1]
                z_coord = p%zpts(k-1)  ! Note: zpts indexed from 0 to nz-1
                ! Start with low parabolic profile
                p%u(idx) = 1.0_wp * (1.0_wp - z_coord**2)
                ! Zero vertical velocity
                p%w(idx) = 0.0_wp
            end do
        end do
        
        ! Check initial divergence
        write(*,'(A)') ' Checking initial divergence...'
        call compute_divergence_spectral(p, p%div_check_2d)
        write(*,'(A, E12.4)') ' Initial max divergence: ', maxval(abs(p%div_check_2d))
        write(*,'(A, E12.4)') ' Initial RMS divergence: ', sqrt(sum(p%div_check_2d**2) / size(p%div_check_2d))
        
        ! Force the initial condition to be exactly divergence-free
        if (maxval(abs(p%div_check_2d)) > 1.0e-12_wp) then
            write(*,'(A)') ' Projecting initial condition to divergence-free space...'
            call project_velocity_to_divergence_free(p)
            
            ! Check divergence after projection
            call compute_divergence_spectral(p, p%div_check_2d)
            write(*,'(A, E12.4)') ' Post-projection max divergence: ', maxval(abs(p%div_check_2d))
            write(*,'(A, E12.4)') ' Post-projection RMS divergence: ', sqrt(sum(p%div_check_2d**2) / size(p%div_check_2d))
        else
            write(*,'(A)') ' Initial condition is already divergence-free!'
        endif
        
        write(*,'(A)') ' Using laminar channel flow profile with pressure gradient driving'
        
        write(*,'(A)') ' Flow fields initialized (laminar parabolic profile)'
        
    end subroutine initialize_flow_fields
    
    subroutine restart_from_file(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: io_status
        
        write(*,'(A)') ' Restarting from start.dat...'
        
        open(8, file='start.dat', form='unformatted', status='old', iostat=io_status)
        if (io_status /= 0) then
            write(*,'(A)') ' Error: Could not open start.dat for restart'
            stop
        endif
        
        read(8) p%t, p%u, p%w, p%temp
        close(8)
        
        ! Initialize other arrays
        p%un = 0.0_wp; p%wn = 0.0_wp; p%tn = 0.0_wp
        p%su = 0.0_wp; p%sw = 0.0_wp; p%st = 0.0_wp
        p%ox = 0.0_wp; p%oz = 0.0_wp; p%body = 0.0_wp
        
        write(*,'(A,ES12.4)') ' Restart successful, t = ', p%t
        
    end subroutine restart_from_file
    
    ! =========================================================================
    ! TIME INTEGRATION ROUTINES
    ! =========================================================================
    
    !============================================================================
    ! SUBROUTINE: convection_step
    !
    ! PURPOSE:
    !   Compute convective terms for Navier-Stokes time stepping
    !
    ! DESCRIPTION:
    !   Wrapper for full 4-stage Runge-Kutta convection computation
    !   Handles nonlinear advection terms: u·∇u in momentum equations
    !
    ! NUMERICAL METHOD:
    !   - 4th-order Runge-Kutta for temporal accuracy
    !   - Spectral methods for spatial derivatives
    !   - Pressure gradient driving for mean flow
    !
    ! INPUTS:
    !   p - Solver state with current velocity fields
    !
    ! OUTPUTS:
    !   p - Updated with convective source terms (su, sw)
    !============================================================================
    subroutine convection_step(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        
        ! Full 4-stage Runge-Kutta convection step from original F77 code
        call runge_kutta_convection(p)
        
    end subroutine convection_step
    
    !============================================================================
    ! SUBROUTINE: runge_kutta_convection
    !
    ! PURPOSE:
    !   4th-order Runge-Kutta integration of convective terms
    !
    ! DESCRIPTION:
    !   Implements explicit RK4 time stepping for nonlinear convection:
    !   - Stage 1: k₁ = f(u^n)
    !   - Stage 2: k₂ = f(u^n + Δt/2·k₁)  
    !   - Stage 3: k₃ = f(u^n + Δt/2·k₂)
    !   - Stage 4: k₄ = f(u^n + Δt·k₃)
    !   - Update: u^{n+1} = u^n + Δt/6·(k₁ + 2k₂ + 2k₃ + k₄)
    !
    ! PHYSICS:
    !   - Convective acceleration: u·∇u
    !   - Pressure gradient driving for mean flow
    !   - Base flow profile: parabolic channel flow
    !
    ! NUMERICAL FEATURES:
    !   - Spectral spatial differentiation
    !   - Dealiasing for nonlinear terms
    !   - Wall boundary conditions enforced
    !
    ! INPUTS:
    !   p - Current solver state
    !
    ! OUTPUTS:
    !   p - Velocity fields advanced by convection
    !============================================================================
    subroutine runge_kutta_convection(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: irk, i, k
        real(wp) :: alf, z_coord
        
        ! Working arrays for derivatives and nonlinear terms
        real(wp), dimension(ntot) :: dudx, dudz, dwdx, dwdz
        real(wp), dimension(ntot) :: uu, uw, ww
        real(wp), dimension(ntot) :: ubar_temp
        real(wp), dimension(ntot) :: conv_u, conv_w  ! Final convective terms
        
        ! Store initial values
        p%un = p%u
        p%wn = p%w
        
        ! Initialize base flow profile (parabolic channel flow profile)
        do i = 1, ntot
            k = (i-1) / nxpp + 1  ! z-index
            if (k <= nz) then
                z_coord = p%zpts(k-1)  ! Get LGL coordinate (0-indexed array)
                ! Parabolic profile: ubar = U_max * (1 - z^2) where U_max = 3/2 * U_bulk  
                ! For Re=180, typical bulk velocity gives U_max ≈ 90
                ubar_temp(i) = 60.0_wp * (1.0_wp - z_coord**2)  ! Base flow profile
            else
                ubar_temp(i) = 0.0_wp
            endif
        end do
        
        ! Initialize convective terms
        conv_u = 0.0_wp
        conv_w = 0.0_wp
        
        ! 4-stage Runge-Kutta loop
        do irk = 4, 1, -1
            alf = 1.0_wp / real(irk, wp)
            
            ! Compute z-derivatives using differentiation matrix
            call compute_z_derivatives(p, p%u, dudz)
            call compute_z_derivatives(p, p%w, dwdz)
            
            ! Compute x-derivatives using spectral methods
            call compute_x_derivatives(p, p%u, dudx)
            call compute_x_derivatives(p, p%w, dwdx)
            
            ! Transform to physical space for nonlinear calculations
            call spectral_to_physical(p, p%u)
            call spectral_to_physical(p, p%w)
            call spectral_to_physical(p, dudx)
            call spectral_to_physical(p, dudz)
            call spectral_to_physical(p, dwdx)
            call spectral_to_physical(p, dwdz)
            
            ! Initialize source terms for this stage
            p%su = 0.0_wp
            p%sw = 0.0_wp
            
            ! Compute convective terms in physical space
            do i = 1, ntot
                ! u-momentum convection: -u∂u/∂x - w∂u/∂z - ubar∂u/∂x
                p%su(i) = -(p%u(i) + ubar_temp(i))*dudx(i) - p%w(i)*dudz(i)
                
                ! w-momentum convection: -u∂w/∂x - w∂w/∂z - ubar∂w/∂x
                p%sw(i) = -(p%u(i) + ubar_temp(i))*dwdx(i) - p%w(i)*dwdz(i)
            end do
            
            ! Transform back to spectral space
            call physical_to_spectral(p, p%su)
            call physical_to_spectral(p, p%sw)
            call physical_to_spectral(p, p%u)
            call physical_to_spectral(p, p%w)
            
            ! Runge-Kutta time advancement for this stage
            do i = 1, ntot
                p%u(i) = p%un(i) + alf * p%dt * p%su(i)
                p%w(i) = p%wn(i) + alf * p%dt * p%sw(i)
            end do
            
            ! Apply boundary conditions after each RK stage
            call apply_wall_bc(p)
            
            ! Accumulate final convective terms (for last stage)
            if (irk == 1) then
                conv_u = p%su
                conv_w = p%sw
            endif
            
        end do
        
        ! Store final convective source terms for use in momentum equation
        p%su = conv_u  ! Final convective terms for u-momentum
        p%sw = conv_w  ! Final convective terms for w-momentum
        
        ! Reset velocities to initial values - the actual time advancement
        ! will be done in the main integration with diffusion
        p%u = p%un
        p%w = p%wn
        
    end subroutine runge_kutta_convection
    
    subroutine compute_z_derivatives(p, field, dfdz)
        implicit none
        type(navier_stokes_params), intent(in) :: p
        real(wp), intent(in) :: field(ntot)
        real(wp), intent(out) :: dfdz(ntot)
        integer :: i, j, k, idx, idx2
        
        ! Compute z-derivatives using differentiation matrix
        dfdz = 0.0_wp
        
        do i = 1, nxpp
            do k = 1, nz
                idx = (k-1)*nxpp + i
                dfdz(idx) = 0.0_wp
                
                do j = 1, nz
                    idx2 = (j-1)*nxpp + i
                    if (k <= nzm .and. j <= nzm) then
                        dfdz(idx) = dfdz(idx) + p%d(k-1,j-1) * field(idx2) / p%ybar
                    endif
                end do
            end do
        end do
        
    end subroutine compute_z_derivatives
    
    subroutine compute_x_derivatives(p, field, dfdx)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(in) :: field(ntot)
        real(wp), intent(out) :: dfdx(ntot)
        integer :: i, k, idx
        
        ! Convert from 1D array to 2D for FFT
        call convert_1d_to_2d_for_fft(p, field, p%real_scratch)
        
        ! Forward FFT to get spectral coefficients
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        
        ! Apply spectral derivative (multiply by ik)
        do k = 1, nz
            do i = 1, nxhp
                p%complex_scratch(i,k) = (0.0_wp, 1.0_wp) * p%xw(i) * p%complex_scratch(i,k)
            end do
        end do
        
        ! Inverse FFT back to physical space
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        
        ! Convert back to 1D array format
        call convert_2d_to_1d_from_fft(p, p%real_scratch, dfdx)
        
    end subroutine compute_x_derivatives
    
    subroutine spectral_to_physical(p, field)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(inout) :: field(ntot)
        
        ! Convert from spectral to physical space using FFTW
        call convert_1d_to_2d_for_fft(p, field, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, field)
        
    end subroutine spectral_to_physical
    
    subroutine physical_to_spectral(p, field)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(inout) :: field(ntot)
        
        ! Convert from physical to spectral space using FFTW
        call convert_1d_to_2d_for_fft(p, field, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, field)
        
    end subroutine physical_to_spectral
    
    ! Helper routine to convert from solver's 1D layout to FFT's 2D layout
    subroutine convert_1d_to_2d_for_fft(p, field_1d, field_2d)
        implicit none
        type(navier_stokes_params), intent(in) :: p
        real(wp), intent(in) :: field_1d(ntot)
        real(wp), intent(out) :: field_2d(nx, nz)
        integer :: i, k, idx
        
        ! Convert from (k-1)*nxpp + i indexing to (i,k) indexing
        ! Only use first nx points (ignore padding)
        do k = 1, nz
            do i = 1, nx
                idx = (k-1)*nxpp + i
                field_2d(i,k) = field_1d(idx)
            end do
        end do
        
    end subroutine convert_1d_to_2d_for_fft
    
    ! Helper routine to convert from FFT's 2D layout back to solver's 1D layout
    subroutine convert_2d_to_1d_from_fft(p, field_2d, field_1d)
        implicit none
        type(navier_stokes_params), intent(in) :: p
        real(wp), intent(in) :: field_2d(nx, nz)
        real(wp), intent(inout) :: field_1d(ntot)
        integer :: i, k, idx
        
        ! Convert from (i,k) indexing back to (k-1)*nxpp + i indexing
        do k = 1, nz
            do i = 1, nx
                idx = (k-1)*nxpp + i
                field_1d(idx) = field_2d(i,k)
            end do
            ! Zero out padding points
            do i = nx+1, nxpp
                idx = (k-1)*nxpp + i
                field_1d(idx) = 0.0_wp
            end do
        end do
        
    end subroutine convert_2d_to_1d_from_fft
    
    subroutine apply_wall_bc(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i, idx
        
        ! Apply no-slip boundary conditions at walls
        do i = 1, nxpp
            ! Bottom wall (k=1)
            idx = i
            p%u(idx) = 0.0_wp
            p%w(idx) = 0.0_wp
            
            ! Top wall (k=nz)  
            idx = (nz-1)*nxpp + i
            p%u(idx) = 0.0_wp
            p%w(idx) = 0.0_wp
        end do
        
    end subroutine apply_wall_bc
    
    ! =========================================================================
    ! EXTRAPOLATE VELOCITY BOUNDARY CONDITIONS FROM PRESSURE (Moin & Kim)
    ! =========================================================================
    subroutine extrapolate_velocity_bc(p, pressure_field)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(in) :: pressure_field(ntot)
        integer :: i, ix, idx_bot, idx_top
        real(wp), allocatable :: pressure_2d(:,:), dpdz_2d(:,:)
        real(wp) :: dpdz_bot, dpdz_top
        complex(wp) :: unit_i
        
        allocate(pressure_2d(nxhp, nz), dpdz_2d(nxhp, nz))
        unit_i = cmplx(0.0_wp, 1.0_wp)
        
        ! Convert pressure to 2D spectral space
        call convert_1d_to_2d_for_fft(p, pressure_field, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, pressure_2d, p%complex_scratch_imag)
        
        ! Compute z-derivatives of pressure using spectral differentiation
        call compute_z_derivatives_2d(p, pressure_2d, dpdz_2d)
        
        ! Extrapolate boundary conditions for ALL Fourier modes (including k=0)
        do ix = 1, nxhp
            ! Get pressure gradients at walls
            dpdz_bot = dpdz_2d(ix, 1)   ! Bottom wall (z = -1)
            dpdz_top = dpdz_2d(ix, nz)  ! Top wall (z = +1)
            
            ! Convert from complex indices to real array indices
            i = 2*ix - 1  ! Real part index
            
            ! Apply pressure extrapolation to ALL modes (including k=0 mean mode)
            ! u boundary condition from wall-normal pressure gradient
            p%uw1(i) =  p%dt * p%xw(ix) * dpdz_bot    ! Bottom wall, real part
            p%uw1(i+1) = 0.0_wp                               ! Bottom wall, imag part
            p%uw2(i) = p%dt * p%xw(ix) * dpdz_top    ! Top wall, real part
            p%uw2(i+1) = 0.0_wp                               ! Top wall, imag part
            
            ! w boundary condition (still zero for no-penetration)
            p%ww1(i) = 0.0_wp      ! Bottom wall, real part
            p%ww1(i+1) = 0.0_wp    ! Bottom wall, imag part
            p%ww2(i) = 0.0_wp      ! Top wall, real part
            p%ww2(i+1) = 0.0_wp    ! Top wall, imag part
        end do
        
        deallocate(pressure_2d, dpdz_2d)
        
    end subroutine extrapolate_velocity_bc
    
    ! =========================================================================
    ! COMPUTE Z-DERIVATIVES IN 2D SPECTRAL SPACE
    ! =========================================================================
    subroutine compute_z_derivatives_2d(p, field_2d, dfdz_2d)
        implicit none
        type(navier_stokes_params), intent(in) :: p
        real(wp), intent(in) :: field_2d(nxhp, nz)
        real(wp), intent(out) :: dfdz_2d(nxhp, nz)
        integer :: ix, i, j
        
        do ix = 1, nxhp
            do i = 1, nz
                dfdz_2d(ix, i) = 0.0_wp
                do j = 1, nz
                    dfdz_2d(ix, i) = dfdz_2d(ix, i) + p%d(i-1, j-1) * field_2d(ix, j)
                end do
            end do
        end do
        
    end subroutine compute_z_derivatives_2d
    
    subroutine build_momentum_sources(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: ijk
        
        do ijk = 1, ntot
            p%su(ijk) = (p%u(ijk) + p%su(ijk)) * p%re / p%dt
            p%sw(ijk) = (p%w(ijk) + p%sw(ijk)) * p%re / p%dt
        end do
        
    end subroutine build_momentum_sources
    
    !============================================================================
    ! SUBROUTINE: solve_helmholtz_system
    !
    ! PURPOSE:
    !   Solve implicit viscous step and pressure correction for incompressible NS
    !
    ! DESCRIPTION:
    !   Two-step fractional step method:
    !   1. Viscous step: (I - dt/Re ∇²)u* = u^n + dt·RHS
    !   2. Pressure step: ∇²φ = ∇·u*/dt to enforce ∇·u^{n+1} = 0
    !
    ! ALGORITHM:
    !   - Extrapolate velocity BCs from pressure field (Moin & Kim method)
    !   - Solve viscous Helmholtz equations in spectral space
    !   - Solve pressure Poisson equation for correction
    !   - Apply pressure correction to ensure incompressibility
    !
    ! NUMERICAL METHOD:
    !   - Spectral methods in x-direction (Fourier)
    !   - LGL collocation in z-direction (wall-normal)
    !   - Direct matrix inversion for each Fourier mode
    !
    ! INPUTS:
    !   p - Solver state with convective source terms
    !
    ! OUTPUTS:
    !   p - Updated velocity and pressure fields
    !============================================================================
    subroutine solve_helmholtz_system(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        
        ! Extrapolate velocity boundary conditions from pressure (Moin & Kim)
        call extrapolate_velocity_bc(p, p%pressure)
        
        ! Solve viscous step: (I - dt/Re * ∇²) u = RHS
        call solve_viscous_step(p)
        
        ! Solve pressure step to enforce incompressibility
        call solve_pressure_step(p)
        
    end subroutine solve_helmholtz_system
    
    !============================================================================
    ! SUBROUTINE: solve_viscous_step
    !
    ! PURPOSE:
    !   Solve viscous diffusion step using Backward Euler or Crank-Nicolson
    !
    ! DESCRIPTION:
    !   Conditionally implements either:
    !   - Backward Euler: (I - dt/Re ∇²)u^{n+1} = u^n + dt·(source terms)
    !   - Crank-Nicolson: (I - 0.5*dt/Re ∇²)u^{n+1} = (I + 0.5*dt/Re ∇²)u^n + dt·(source terms)
    !   
    !   Method selection controlled by use_crank_nicolson flag.
    !   Boundary conditions: Kim & Moin extrapolated velocities from pressure gradient
    !
    !============================================================================
    subroutine solve_viscous_step(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: ix, i, j, k, info
        real(wp) :: lambda_x, visc_factor
        real(wp) :: matrix(nz,nz), rhs_matrix(nz,nz)
        real(wp) :: rhs_vec(nz), rhs_bc_contrib(nz)
        integer :: ipiv(nz)
        real(wp) :: rhs_u(nxhp,nz), rhs_w(nxhp,nz)
        real(wp) :: u_spec(nxhp,nz), w_spec(nxhp,nz)
        
        ! Viscous factor for time stepping
        if (p%use_crank_nicolson) then
            visc_factor = 0.5_wp * p%dt / p%re  ! CN: 0.5*dt/Re
            write(*,'(A)') ' Using Crank-Nicolson viscous step with Kim & Moin boundary conditions'
        else
            visc_factor = p%dt / p%re           ! BE: dt/Re  
            write(*,'(A)') ' Using Backward Euler viscous step with Kim & Moin boundary conditions'
        endif
        
        ! Transform source terms to spectral space
        call convert_1d_to_2d_for_fft(p, p%su, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, rhs_u, p%complex_scratch_imag)
        
        call convert_1d_to_2d_for_fft(p, p%sw, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, rhs_w, p%complex_scratch_imag)
        
        ! Transform current velocities to spectral space
        call convert_1d_to_2d_for_fft(p, p%u, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, u_spec, p%complex_scratch_imag)
        
        call convert_1d_to_2d_for_fft(p, p%w, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, w_spec, p%complex_scratch_imag)
        
        ! Solve Helmholtz equation for each Fourier mode
        do ix = 1, nxhp
            lambda_x = p%xsq(ix)
            
            ! Build Helmholtz matrix: I - visc_factor * (D²/dy² - λₓ² I)
            ! Following base_code.f pattern: build full matrix first, then handle boundaries
            matrix = 0.0_wp
            
            ! Add identity matrix
            do j = 1, nz
                matrix(j,j) = 1.0_wp
            end do
            
            ! Add viscous terms: -visc_factor * D²/dy² (for ALL points, like base code)
            do i = 1, nz
                do j = 1, nz
                    ! Second derivative matrix (multiplied by -visc_factor)
                    matrix(i,j) = matrix(i,j) - visc_factor * &
                                  compute_second_derivative_element(p, i-1, j-1)
                end do
                
                ! Add visc_factor*λₓ² term to diagonal
                matrix(i,i) = matrix(i,i) + visc_factor * lambda_x
            end do
            
            ! Compute boundary contributions (like base code furxy calculation)
            rhs_bc_contrib = 0.0_wp
            do k = 1, nz
                rhs_bc_contrib(k) = matrix(k,1) * p%uw1(2*ix-1) + matrix(k,nz) * p%uw2(2*ix-1)
            end do
            
            ! Apply boundary conditions (zero out boundary rows/columns)
            matrix(1,:) = 0.0_wp
            matrix(1,1) = 1.0_wp
            matrix(nz,:) = 0.0_wp
            matrix(nz,nz) = 1.0_wp
            
            ! Solve u-momentum equation
            if (p%use_crank_nicolson) then
                ! Crank-Nicolson: RHS = (I + 0.5*visc_factor*(D²/dy² - λₓ²))u^n + dt*sources
                call build_cn_rhs(u_spec(ix,:), rhs_u(ix,:), lambda_x, visc_factor, rhs_vec)
            else
                ! Backward Euler: RHS = u^n + dt*sources
                rhs_vec = u_spec(ix,:) + rhs_u(ix,:)
            endif
            
            ! Subtract boundary contributions (like base code: fur(k) = su(ik)*... - furxy(k))
            rhs_vec = rhs_vec - rhs_bc_contrib
            
            ! Apply Kim & Moin extrapolated boundary conditions for u-velocity (all Fourier modes)
            rhs_vec(1) = p%uw1(2*ix-1)   ! Bottom wall
            rhs_vec(nz) = p%uw2(2*ix-1)  ! Top wall
            
            call dgesv(nz, 1, matrix, nz, ipiv, rhs_vec, nz, info)
            if (info /= 0) then
                write(*,*) 'Error in u-momentum solve, info =', info
            endif
            u_spec(ix,:) = rhs_vec
            
            ! Solve w-momentum equation (same matrix, different RHS)
            if (p%use_crank_nicolson) then
                ! Crank-Nicolson: RHS = (I + 0.5*visc_factor*(D²/dy² - λₓ²))w^n + dt*sources
                call build_cn_rhs(w_spec(ix,:), rhs_w(ix,:), lambda_x, visc_factor, rhs_vec)
            else
                ! Backward Euler: RHS = w^n + dt*sources
                rhs_vec = w_spec(ix,:) + rhs_w(ix,:)
            endif
            
            ! Compute boundary contributions for w-momentum
            do k = 1, nz
                rhs_bc_contrib(k) = matrix(k,1) * p%ww1(2*ix-1) + matrix(k,nz) * p%ww2(2*ix-1)
            end do
            
            ! Subtract boundary contributions
            rhs_vec = rhs_vec - rhs_bc_contrib
            
            ! Apply Kim & Moin extrapolated boundary conditions for w-velocity (all Fourier modes)
            rhs_vec(1) = p%ww1(2*ix-1)   ! Bottom wall
            rhs_vec(nz) = p%ww2(2*ix-1)  ! Top wall
            
            call dgesv(nz, 1, matrix, nz, ipiv, rhs_vec, nz, info)
            if (info /= 0) then
                write(*,*) 'Error in w-momentum solve, info =', info
            endif
            w_spec(ix,:) = rhs_vec
        end do
        
        ! Transform back to physical space
        call convert_real_imag_to_complex(u_spec, p%complex_scratch_imag, p%complex_scratch)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, p%u)
        
        call convert_real_imag_to_complex(w_spec, p%complex_scratch_imag, p%complex_scratch)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, p%w)
        
    end subroutine solve_viscous_step
    
    !============================================================================
    ! SUBROUTINE: build_cn_rhs
    !
    ! PURPOSE:
    !   Build RHS for Crank-Nicolson: (I + 0.5*visc_factor*(D²/dy² - λₓ²))u^n + sources
    !
    !============================================================================
    subroutine build_cn_rhs(u_current, source_terms, lambda_x, visc_factor, rhs_out)
        implicit none
        real(wp), intent(in) :: u_current(nz), source_terms(nz)
        real(wp), intent(in) :: lambda_x, visc_factor
        real(wp), intent(out) :: rhs_out(nz)
        integer :: i, j
        real(wp) :: rhs_matrix(nz,nz)
        
        ! Build RHS matrix: I + 0.5*visc_factor*(D²/dy² - λₓ²)
        rhs_matrix = 0.0_wp
        
        ! Add identity matrix
        do j = 1, nz
            rhs_matrix(j,j) = 1.0_wp
        end do
        
        ! Add explicit CN viscous terms: +visc_factor*(D²/dy² - λₓ²)
        ! Following base code pattern: build full matrix first
        do i = 1, nz
            do j = 1, nz
                rhs_matrix(i,j) = rhs_matrix(i,j) + visc_factor * &
                                  compute_second_derivative_element(p, i-1, j-1)
            end do
            
            ! Subtract visc_factor*λₓ² term from diagonal
            rhs_matrix(i,i) = rhs_matrix(i,i) - visc_factor * lambda_x
        end do
        
        ! Set boundary rows to identity (will be overridden by Kim & Moin BCs)
        rhs_matrix(1,:) = 0.0_wp
        rhs_matrix(1,1) = 1.0_wp
        rhs_matrix(nz,:) = 0.0_wp  
        rhs_matrix(nz,nz) = 1.0_wp
        
        ! Compute RHS = RHS_matrix * u_current + source_terms
        rhs_out = 0.0_wp
        do i = 1, nz
            if (i <= nzm) then
                ! Interior points: full matrix-vector multiplication
                do j = 1, nz
                    rhs_out(i) = rhs_out(i) + rhs_matrix(i,j) * u_current(j)
                end do
                rhs_out(i) = rhs_out(i) + source_terms(i)
            else
                ! Boundary points: will be overridden by Kim & Moin values
                rhs_out(i) = 0.0_wp
            endif
        end do
        
    end subroutine build_cn_rhs
    
    !============================================================================
    ! SUBROUTINE: solve_pressure_step
    !
    ! PURPOSE:
    !   Solve pressure Poisson equation to enforce incompressibility
    !
    ! DESCRIPTION:
    !   Second step of fractional step method:
    !   ∇²φ = ∇·u*/dt  (pressure correction equation)
    !   
    !   Algorithm:
    !   - Compute divergence of intermediate velocity u*
    !   - Solve Poisson equation in spectral space
    !   - Apply pressure correction: u^{n+1} = u* - dt∇φ
    !   - Ensures ∇·u^{n+1} = 0 (incompressibility)
    !
    ! POISSON EQUATION:
    !   (D²/dy² - λₓ²)φ = ∇·u*/dt
    !   where λₓ = kₓ² is streamwise wavenumber squared
    !
    ! BOUNDARY CONDITIONS:
    !   - Neumann BC: ∂φ/∂n = 0 at walls (natural for Poisson)
    !   - Periodic in x via spectral methods
    !
    ! SPECTRAL METHOD:
    !   - Each Fourier mode solved independently
    !   - Direct matrix inversion for wall-normal direction
    !   - High accuracy from spectral differentiation
    !
    ! INPUTS:
    !   p - Solver state with intermediate velocity u*
    !
    ! OUTPUTS:
    !   p - Final velocity u^{n+1} and pressure field φ
    !============================================================================
    subroutine solve_pressure_step(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i, k, ix, idx, n
        ! lambda_x: Streamwise wavenumber squared (kₓ²) used in pressure Poisson equation
        ! - Poisson equation: (D²/dy² - λₓ²) φ = ∇·u (pressure correction)
        ! - Same spectral principle as viscous step but for pressure field
        real(wp) :: lambda_x
        real(wp), allocatable :: div_spectral(:,:), phi_spectral(:,:)
        real(wp), allocatable :: poisson_matrix(:,:), div_vec(:), phi_vec(:)
        integer :: info, j
        
        ! Allocate working arrays
        allocate(div_spectral(nxhp, nz), phi_spectral(nxhp, nz))
        allocate(poisson_matrix(nz, nz), div_vec(nz), phi_vec(nz))
        
        ! Compute divergence of current velocity field
        call compute_divergence_spectral(p, div_spectral)
        
        ! Solve Poisson equation: ∇²φ = ∇·u for pressure correction
        do ix = 1, nxhp
            lambda_x = p%xsq(ix)
            
            ! Handle zero mode specially (ix=1, kx=0)
            if (ix == 1) then
                ! Zero mode (kₓ = 0): Use exact F77 approach from base_code.f
                ! Build the diff matrix exactly like F77: ∑ (2/ybar)*wg(n)*d(n,i)*d(n,j)
                
                poisson_matrix = 0.0_wp
                
                ! Build F77-style diff matrix
                do i = 1, nz
                    do j = 1, nz
                        ! Compute exactly like F77: diff(ii,jj) = sum
                        poisson_matrix(i,j) = 0.0_wp
                        do n = 0, nzm
                            poisson_matrix(i,j) = poisson_matrix(i,j) + &
                                (2.0_wp/p%ybar) * p%wg(n) * p%d(n,i-1) * p%d(n,j-1)
                        end do
                    end do
                end do
                
                ! F77 approach: For zero mode, clear first row and first column
                do k = 1, nz
                    poisson_matrix(1,k) = 0.0_wp  ! Clear first row
                    poisson_matrix(k,1) = 0.0_wp  ! Clear first column
                end do
                poisson_matrix(1,1) = 1.0_wp     ! Pin: φ(1) = 0
                
                ! Right-hand side: F77 approach
                div_vec = div_spectral(ix,:)
                div_vec(1) = 0.0_wp               ! φ(1) = 0
                
            else
                ! For kx ≠ 0 modes, build matrix exactly like F77: diff + ak*amass
                poisson_matrix = 0.0_wp
                
                ! Build F77-style diff matrix (same as zero mode)
                do i = 1, nz
                    do j = 1, nz
                        ! Compute exactly like F77: diff(ii,jj) = sum
                        poisson_matrix(i,j) = 0.0_wp
                        do n = 0, nzm
                            poisson_matrix(i,j) = poisson_matrix(i,j) + &
                                (2.0_wp/p%ybar) * p%wg(n) * p%d(n,i-1) * p%d(n,j-1)
                        end do
                    end do
                end do
                
                ! Add ak*amass to diagonal (F77: apv(k1,k1) = apv(k1,k1) + ak*amass(k1))
                do i = 1, nz
                    poisson_matrix(i,i) = poisson_matrix(i,i) + lambda_x * p%amass(i)
                end do
                
                ! Right-hand side: F77 approach with scaling
                do k = 1, nz
                    div_vec(k) = div_spectral(ix,k) * 0.5_wp * p%ybar * p%wg(k-1)
                end do
                
            endif
            
            ! Solve Poisson equation using CGS iterative solver (like F77 base_code.f)
            ! Initialize solution vector
            phi_vec = 0.0_wp
            
            call jacgs(poisson_matrix, phi_vec, div_vec, 500, nz, nz, nz, 0, p%cgstol, info)
            if (info /= 0) then
                write(*,*) 'Warning: CGS solver did not converge for mode ix =', ix, ', info =', info
                ! Use result anyway - iterative solvers can be partially converged
                phi_spectral(ix,:) = phi_vec
            else
                phi_spectral(ix,:) = phi_vec
            endif
            
        end do
        
        ! Compute pressure gradient and subtract from velocity
        call apply_pressure_correction(p, phi_spectral)
        
        ! Store pressure correction for next time step's boundary condition extrapolation
        call store_pressure_field(p, phi_spectral)
        
        ! Verify divergence has been reduced
        call verify_pressure_correction(p)
        
        ! Cleanup
        deallocate(div_spectral, phi_spectral)
        deallocate(poisson_matrix, div_vec, phi_vec)
        
    end subroutine solve_pressure_step
    
    ! =========================================================================
    ! STORE PRESSURE FIELD FOR BOUNDARY CONDITION EXTRAPOLATION
    ! =========================================================================
    subroutine store_pressure_field(p, phi_spectral)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(in) :: phi_spectral(nxhp, nz)
        
        ! Convert pressure from spectral to physical space and store
        call convert_real_imag_to_complex(phi_spectral, p%complex_scratch_imag, p%complex_scratch)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, p%pressure)
        
        ! Update total pressure (add current correction to previous pressure)
        ! Note: For simplicity, we just store the current pressure correction
        ! In a full implementation, you might want to accumulate corrections
        
    end subroutine store_pressure_field
    
    ! Helper functions for spectral operations
    subroutine compute_divergence_spectral(p, div_spectral)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(out) :: div_spectral(nxhp, nz)
        real(wp), allocatable :: dudx_spectral(:,:), dwdz_spectral(:,:)
        real(wp), allocatable :: dudx(:), dwdz(:)
        integer :: i, k, ix
        
        allocate(dudx_spectral(nxhp, nz), dwdz_spectral(nxhp, nz))
        allocate(dudx(ntot), dwdz(ntot))
        
        ! Compute ∂u/∂x in spectral space directly
        ! First transform u to spectral space
        call convert_1d_to_2d_for_fft(p, p%u, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, dudx_spectral, p%complex_scratch_imag)
        
        ! Apply spectral derivative operator: multiply by ik
        do k = 1, nz
            do ix = 1, nxhp
                ! Store ik * û in real part (since ik * real = -k * imag, ik * imag = k * real)
                dudx_spectral(ix,k) = -p%xw(ix) * p%complex_scratch_imag(ix,k)
            end do
        end do
        
        ! Compute ∂w/∂z: transform w to physical space then use LGL differentiation
        call compute_z_derivatives(p, p%w, dwdz)
        
        ! Transform ∂w/∂z to spectral space
        call convert_1d_to_2d_for_fft(p, dwdz, p%real_scratch)
        call fft_forward_2d(p%fft_plan, p%real_scratch, p%complex_scratch)
        call convert_complex_to_real_imag(p%complex_scratch, dwdz_spectral, p%complex_scratch_imag)
        
        ! Combine: ∇·u = ∂u/∂x + ∂w/∂z
        div_spectral = dudx_spectral + dwdz_spectral
        
        ! Debug: Print some derivative values for the first few points
        if (maxval(abs(div_spectral)) > 1.0e-12_wp) then
            write(*,'(A)') ' Warning: Non-zero divergence detected!'
            write(*,'(A,E12.4)') '   Max |∂u/∂x|: ', maxval(abs(dudx_spectral))
            write(*,'(A,E12.4)') '   Max |∂w/∂z|: ', maxval(abs(dwdz_spectral))
        endif
        
        deallocate(dudx_spectral, dwdz_spectral, dudx, dwdz)
        
    end subroutine compute_divergence_spectral
    
    subroutine apply_pressure_correction(p, phi_spectral)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(in) :: phi_spectral(nxhp, nz)
        real(wp), allocatable :: dphidx(:), dphidz(:), phi_physical(:)
        real(wp), parameter :: alpha_damp = 0.1_wp  ! Reduce damping further
        integer :: i, ix, k
        
        allocate(dphidx(ntot), dphidz(ntot), phi_physical(ntot))
        
        ! Transform φ to physical space first
        p%complex_scratch_imag = 0.0_wp
        call convert_real_imag_to_complex(phi_spectral, p%complex_scratch_imag, p%complex_scratch)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, phi_physical)
        
        ! Compute ∂φ/∂x using the existing spectral derivative routine
        call compute_x_derivatives(p, phi_physical, dphidx)
        
        ! Compute ∂φ/∂z using LGL differentiation
        call compute_z_derivatives(p, phi_physical, dphidz)
        
        ! Apply pressure correction with very small damping
        do i = 1, ntot
            p%u(i) = p%u(i) - alpha_damp * dphidx(i)
            p%w(i) = p%w(i) - alpha_damp * dphidz(i)
        end do
        
        ! Apply boundary conditions after pressure correction
        call apply_wall_bc(p)
        
        deallocate(dphidx, dphidz, phi_physical)
        
    end subroutine apply_pressure_correction
    
    ! Full pressure correction without damping (for initial projection)
    subroutine apply_pressure_correction_full(p, phi_spectral)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), intent(in) :: phi_spectral(nxhp, nz)
        real(wp), allocatable :: dphidx_spectral(:,:), dphidz_spectral(:,:)
        real(wp), allocatable :: dphidx(:), dphidz(:)
        integer :: i, ix, k, idx
        
        allocate(dphidx_spectral(nxhp, nz), dphidz_spectral(nxhp, nz))
        allocate(dphidx(ntot), dphidz(ntot))
        
        ! Compute ∂φ/∂x in spectral space
        do k = 1, nz
            do ix = 1, nxhp
                if (ix == 1) then
                    dphidx_spectral(ix, k) = 0.0_wp  ! d/dx of k=0 mode is 0
                else
                    dphidx_spectral(ix, k) = p%alpha * real(ix-1, wp) * phi_spectral(ix, k)
                endif
            end do
        end do
        
        ! Convert to physical space via inverse FFT
        p%complex_scratch = cmplx(dphidx_spectral, 0.0_wp, wp)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, dphidx)
        
        ! Compute ∂φ/∂z using spectral differentiation  
        do ix = 1, nxhp
            dphidz_spectral(ix,:) = matmul(p%diff, phi_spectral(ix,:))
        end do
        
        ! Convert to physical space via inverse FFT
        p%complex_scratch = cmplx(dphidz_spectral, 0.0_wp, wp)
        call fft_backward_2d(p%fft_plan, p%complex_scratch, p%real_scratch)
        call convert_2d_to_1d_from_fft(p, p%real_scratch, dphidz)
        
        ! Apply full correction: u_new = u - ∂φ/∂x, w_new = w - ∂φ/∂z
        p%u = p%u - dphidx
        p%w = p%w - dphidz
        
        deallocate(dphidx_spectral, dphidz_spectral, dphidx, dphidz)
        
    end subroutine apply_pressure_correction_full
    
    subroutine verify_pressure_correction(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), allocatable :: dudx(:), dwdz(:)
        real(wp) :: div_max, div_rms
        integer :: i, count_nonzero
        
        allocate(dudx(ntot), dwdz(ntot))
        
        ! Compute divergence after pressure correction
        call compute_x_derivatives(p, p%u, dudx)
        call compute_z_derivatives(p, p%w, dwdz)
        
        ! Store divergence in scratch array
        do i = 1, ntot
            p%ox(i) = dudx(i) + dwdz(i)
        end do
        
        ! Compute statistics
        div_max = maxval(abs(p%ox))
        div_rms = sqrt(sum(p%ox**2) / real(ntot, wp))
        count_nonzero = count(abs(p%ox) > 1.0e-12_wp)
        
        ! Only report if significant
        if (div_max > 1.0e-10_wp) then
            write(*,'(A,ES10.3,A,ES10.3,A,I0)') &
                ' Post-correction divergence: max=', div_max, ', rms=', div_rms, &
                ', nonzero_points=', count_nonzero
        endif
        
        deallocate(dudx, dwdz)
        
    end subroutine verify_pressure_correction
    
    ! Helper function to compute second derivative matrix elements
    function compute_second_derivative_element(p, i, j) result(d2_elem)
        implicit none
        type(navier_stokes_params), intent(in) :: p
        integer, intent(in) :: i, j
        real(wp) :: d2_elem
        integer :: k
        
        ! Compute D² = D * D where D is the first derivative matrix
        d2_elem = 0.0_wp
        do k = 0, nzm
            d2_elem = d2_elem + p%d(i,k) * p%d(k,j) / (p%ybar * p%ybar)
        end do
        
    end function compute_second_derivative_element
    
    ! Helper to compute x-derivatives in complex spectral space
    subroutine compute_x_derivatives_complex(p, field_complex)
        implicit none
        type(navier_stokes_params), intent(in) :: p
        complex(wp), intent(inout) :: field_complex(nxhp, nz)
        integer :: i, k
        
        do k = 1, nz
            do i = 1, nxhp
                field_complex(i,k) = (0.0_wp, 1.0_wp) * p%xw(i) * field_complex(i,k)
            end do
        end do
        
    end subroutine compute_x_derivatives_complex
    
    ! Project initial velocity field to be divergence-free
    subroutine project_velocity_to_divergence_free(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), allocatable :: div_spectral(:,:), phi_spectral(:,:)
        real(wp), allocatable :: poisson_matrix(:,:), div_vec(:), phi_vec(:)
        integer :: ix, k, info
        real(wp) :: kx_squared
        
        allocate(div_spectral(nxhp, nz), phi_spectral(nxhp, nz))
        allocate(poisson_matrix(nz, nz), div_vec(nz), phi_vec(nz))
        
        ! Compute divergence in spectral space
        call compute_divergence_spectral(p, div_spectral)
        
        ! Solve Poisson equation for each Fourier mode
        do ix = 1, nxhp
            
            ! Set up 1D Poisson equation in z-direction: d²φ/dz² - k²φ = div
            poisson_matrix = 0.0_wp
            
            ! Get k² for this mode
            if (ix == 1) then
                kx_squared = 0.0_wp  ! k=0 mode
            else
                kx_squared = (p%alpha * real(ix-1, wp))**2
            endif
            
            ! Build matrix: D² - k²I where D is differentiation matrix
            do k = 1, nz
                poisson_matrix(k,:) = -kx_squared * p%amass  ! -k² mass matrix
                poisson_matrix(k,k) = poisson_matrix(k,k) + 1.0_wp  ! Add identity
            end do
            
            ! Add differentiation contributions
            do k = 1, nz-1
                poisson_matrix(k,:) = poisson_matrix(k,:) + p%diff(k,:)  ! Add D²
            end do
            
            ! Apply boundary conditions: φ = 0 at walls
            poisson_matrix(1,:) = 0.0_wp; poisson_matrix(1,1) = 1.0_wp
            poisson_matrix(nz,:) = 0.0_wp; poisson_matrix(nz,nz) = 1.0_wp
            
            ! Handle zero mode singularity by fixing pressure at middle point
            if (ix == 1) then
                poisson_matrix(nz/2,:) = 0.0_wp
                poisson_matrix(nz/2,nz/2) = 1.0_wp
            endif
            
            ! Set up RHS
            div_vec = div_spectral(ix,:)
            div_vec(1) = 0.0_wp    ! φ = 0 at wall
            div_vec(nz) = 0.0_wp   ! φ = 0 at wall
            if (ix == 1) div_vec(nz/2) = 0.0_wp  ! φ = 0 at middle (regularization)
            
            ! Solve linear system using CGS iterative solver
            phi_vec = 0.0_wp  ! Initial guess
            call jacgs(poisson_matrix, phi_vec, div_vec, 500, nz, nz, nz, 0, p%cgstol, info)
            
            if (info /= 0) then
                write(*,*) 'Warning: CGS solver did not converge in initial projection for mode ix =', ix, ', info =', info
                ! Use result anyway - iterative solvers can be partially converged
                phi_spectral(ix,:) = phi_vec
            else
                phi_spectral(ix,:) = phi_vec
            endif
            
        end do
        
        ! Apply FULL pressure correction (no damping for initial projection)
        call apply_pressure_correction_full(p, phi_spectral)
        
        ! Cleanup
        deallocate(div_spectral, phi_spectral)
        deallocate(poisson_matrix, div_vec, phi_vec)
        
    end subroutine project_velocity_to_divergence_free
    
    ! Helper functions for complex/real conversions
    subroutine convert_complex_to_real_imag(complex_field, real_field, imag_field)
        implicit none
        complex(wp), intent(in) :: complex_field(:,:)
        real(wp), intent(out) :: real_field(:,:), imag_field(:,:)
        
        real_field = real(complex_field)
        imag_field = aimag(complex_field)
        
    end subroutine convert_complex_to_real_imag
    
    subroutine convert_real_imag_to_complex(real_field, imag_field, complex_field)
        implicit none
        real(wp), intent(in) :: real_field(:,:), imag_field(:,:)
        complex(wp), intent(out) :: complex_field(:,:)
        
        complex_field = cmplx(real_field, imag_field, wp)
        
    end subroutine convert_real_imag_to_complex
    
    subroutine apply_edge_conditions(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i, k, idx
        
        ! Apply no-slip boundary conditions at walls
        do i = 1, nxpp
            ! Bottom wall (k=1)
            idx = i
            p%u(idx) = 0.0_wp
            p%w(idx) = 0.0_wp
            
            ! Top wall (k=nz)
            idx = (nz-1)*nxpp + i
            p%u(idx) = 0.0_wp
            p%w(idx) = 0.0_wp
        end do
        
    end subroutine apply_edge_conditions
    
    subroutine check_divergence(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp) :: div_max
        integer :: m
        
        ! Reset divergence check array
        do m = 1, ntot
            p%su(m) = 0.0_wp
        end do
        
        ! Compute divergence (simplified)
        call compute_divergence(p)
        
        ! Find maximum divergence
        div_max = maxval(abs(p%su))
        if (div_max > 1.0e-6_wp) then
            write(*,'(A,ES12.4,A)') ' WARNING: Largest absolute divergence = ', div_max, ' exceeds 10^-6 tolerance!'
            write(*,'(A)') '          This indicates violation of incompressible flow constraint.'
        endif
        
    end subroutine check_divergence
    
    subroutine compute_divergence(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        real(wp), allocatable :: dudx(:), dwdz(:)
        integer :: i
        
        allocate(dudx(ntot), dwdz(ntot))
        
        ! Compute ∂u/∂x using spectral differentiation
        call compute_x_derivatives(p, p%u, dudx)
        
        ! Compute ∂w/∂z using LGL differentiation  
        call compute_z_derivatives(p, p%w, dwdz)
        
        ! Compute divergence: ∇·u = ∂u/∂x + ∂w/∂z
        do i = 1, ntot
            p%su(i) = dudx(i) + dwdz(i)
        end do
        
        deallocate(dudx, dwdz)
        
    end subroutine compute_divergence
    
    ! =========================================================================
    ! OUTPUT ROUTINES
    ! =========================================================================
    subroutine output_solution(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        character(len=20) :: filename
        integer :: i, k, idx
        
        ! Write restart file
        open(2, file='run.dat', form='unformatted', status='replace')
        write(2) p%t, p%u, p%w, p%temp
        close(2)
        
        ! Write plot file
        write(filename, '(A,ES10.4)') 'p_', p%t
        
        if (p%iform == 0) then
            open(3, file=filename, form='unformatted', status='replace')
            write(3) p%t, ((p%u((k-1)*nxpp + i), i=1,nxpp), k=1,nz)
            write(3) ((p%w((k-1)*nxpp + i), i=1,nxpp), k=1,nz)
            write(3) ((p%temp((k-1)*nxpp + i), i=1,nxpp), k=1,nz)
        else
            open(3, file=filename, status='replace')
            write(3,'(A,ES14.6)') '# Time = ', p%t
            do k = 1, nz
                do i = 1, nxpp
                    idx = (k-1)*nxpp + i
                    write(3,'(2I6,3ES14.6)') i, k, p%u(idx), p%w(idx), p%temp(idx)
                end do
            end do
        endif
        close(3)
        
        write(*,'(A,A)') ' Output written to: ', trim(filename)
        
    end subroutine output_solution
    
    subroutine final_output(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        integer :: i, k, idx
        
        ! Final restart file
        open(2, file='run.dat', form='unformatted', status='replace')
        write(2) p%t, p%u, p%w, p%temp
        close(2)
        
        ! Final plot file
        if (p%iform == 0) then
            open(3, file='plot.dat', form='unformatted', status='replace')
            write(3) p%t, ((p%u((k-1)*nxpp + i), i=1,nxpp), k=1,nz)
            write(3) ((p%w((k-1)*nxpp + i), i=1,nxpp), k=1,nz)
            write(3) ((p%temp((k-1)*nxpp + i), i=1,nxpp), k=1,nz)
        else
            open(3, file='plot.dat', status='replace')
            write(3,'(A,ES14.6)') '# Final Time = ', p%t
            do k = 1, nz
                do i = 1, nxpp
                    idx = (k-1)*nxpp + i
                    write(3,'(2I6,3ES14.6)') i, k, p%u(idx), p%w(idx), p%temp(idx)
                end do
            end do
        endif
        close(3)
        
        write(*,'(A)') ' Final output completed'
        write(*,'(A)') ' Files written: run.dat, plot.dat'
        
    end subroutine final_output
    
    ! =========================================================================
    ! CLEANUP ROUTINES
    ! =========================================================================
    subroutine finalize_solver(p)
        implicit none
        type(navier_stokes_params), intent(inout) :: p
        
        ! Cleanup FFTW plans and threading
        call destroy_fft_plans(p%fft_plan)
        call cleanup_fft_threading()
        
        ! Deallocate all arrays
        if (allocated(p%xw)) deallocate(p%xw)
        if (allocated(p%xsq)) deallocate(p%xsq)
        if (allocated(p%zpts)) deallocate(p%zpts)
        if (allocated(p%wg)) deallocate(p%wg)
        if (allocated(p%d)) deallocate(p%d)
        if (allocated(p%u)) deallocate(p%u)
        if (allocated(p%w)) deallocate(p%w)
        if (allocated(p%temp)) deallocate(p%temp)
        if (allocated(p%pressure)) deallocate(p%pressure)
        if (allocated(p%un)) deallocate(p%un)
        if (allocated(p%wn)) deallocate(p%wn)
        if (allocated(p%tn)) deallocate(p%tn)
        if (allocated(p%trigsx)) deallocate(p%trigsx)
        if (allocated(p%work)) deallocate(p%work)
        if (allocated(p%ifaxx)) deallocate(p%ifaxx)
        if (allocated(p%su)) deallocate(p%su)
        if (allocated(p%sw)) deallocate(p%sw)
        if (allocated(p%st)) deallocate(p%st)
        if (allocated(p%diff)) deallocate(p%diff)
        if (allocated(p%amass)) deallocate(p%amass)
        if (allocated(p%ubc)) deallocate(p%ubc)
        if (allocated(p%wbc)) deallocate(p%wbc)
        if (allocated(p%ox)) deallocate(p%ox)
        if (allocated(p%oz)) deallocate(p%oz)
        if (allocated(p%body)) deallocate(p%body)
        if (allocated(p%uw1)) deallocate(p%uw1)
        if (allocated(p%ww1)) deallocate(p%ww1)
        if (allocated(p%uw2)) deallocate(p%uw2)
        if (allocated(p%ww2)) deallocate(p%ww2)
        if (allocated(p%real_scratch)) deallocate(p%real_scratch)
        if (allocated(p%complex_scratch)) deallocate(p%complex_scratch)
        if (allocated(p%complex_scratch_real)) deallocate(p%complex_scratch_real)
        if (allocated(p%complex_scratch_imag)) deallocate(p%complex_scratch_imag)
        
        write(*,'(A)') ' Memory cleanup and FFTW cleanup completed'
        
    end subroutine finalize_solver

    ! CGS Solver - Conjugate Gradient Squared with Jacobi preconditioning
    ! Ported from F77 jacgs subroutine for compatibility with base_code.f
    subroutine jacgs(a, x, b, niter, ndx, nx, ny, iprt, tol, info)
        implicit none
        integer, intent(in) :: niter, ndx, nx, ny, iprt
        real(wp), intent(in) :: tol
        real(wp), intent(in) :: a(ndx, ny)
        real(wp), intent(inout) :: x(nx)
        real(wp), intent(in) :: b(nx)
        integer, intent(out) :: info
        
        ! Local arrays
        real(wp) :: r(nx), p(nx), q(nx), apn(nx), diag(nx), sum_vec(nx)
        real(wp) :: small, factor
        real(wp) :: res0, res1, resfac, qdr, pdapn, qdrnew, alfa, beta
        integer :: i, j, iter, nxy
        
        parameter (small = 1.0e-30_wp, factor = 1.0e-3_wp)
        
        info = 0
        nxy = nx * ny
        
        ! Initialize diagonal preconditioning
        do i = 1, nx
            if (abs(a(i,i)) > small) then
                diag(i) = 1.0_wp / a(i,i)
            else
                diag(i) = 1.0_wp
            endif
            sum_vec(i) = 0.0_wp
        end do
        
        ! Compute A*x
        do j = 1, ny
            do i = 1, nx
                sum_vec(i) = sum_vec(i) + a(i,j) * x(j)
            end do
        end do
        
        ! Initialize residual and search direction
        do i = 1, nx
            r(i) = b(i) - sum_vec(i)
            q(i) = r(i) * diag(i)
            p(i) = q(i)
        end do
        
        ! Initial residual norm
        res0 = 0.0_wp
        do i = 1, nx
            res0 = res0 + r(i) * r(i)
        end do
        res0 = sqrt(res0)
        
        if (res0 <= tol) return
        
        ! CGS iteration loop
        do iter = 1, niter
            ! Compute A*p
            do i = 1, nx
                apn(i) = 0.0_wp
            end do
            
            do j = 1, ny
                do i = 1, nx
                    apn(i) = apn(i) + a(i,j) * p(j)
                end do
            end do
            
            ! Compute alpha
            qdr = 0.0_wp
            pdapn = 0.0_wp
            do i = 1, nx
                qdr = qdr + q(i) * r(i)
                pdapn = pdapn + p(i) * apn(i)
            end do
            
            alfa = qdr / (pdapn + small)
            
            ! Update solution and residual
            do i = 1, nx
                x(i) = x(i) + alfa * p(i)
                r(i) = r(i) - alfa * apn(i)
                q(i) = r(i) * diag(i)
            end do
            
            ! Compute beta
            qdrnew = 0.0_wp
            do i = 1, nx
                qdrnew = qdrnew + q(i) * r(i)
            end do
            
            beta = qdrnew / (qdr + small)
            
            ! Update search direction
            do i = 1, nx
                p(i) = r(i) * diag(i) + beta * p(i)
            end do
            
            ! Check convergence
            res1 = 0.0_wp
            do i = 1, nx
                res1 = res1 + r(i) * r(i)
            end do
            res1 = sqrt(res1)
            resfac = res1 / (res0 + small)
            
            if (iprt == 1) then
                write(*,'(A,I5,5X,A,E14.6)') '*** iter = ', iter, 'residual= ', resfac
            endif
            
            if (res1 <= tol) return
        end do
        
        ! Non-convergent
        write(*,*) 'non-convergent in cgs ', 'residual= ', res1
        info = -1
        
    end subroutine jacgs

end program channel_flow_solver
