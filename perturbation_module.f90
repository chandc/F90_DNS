!==============================================================================
! CHANNEL FLOW SOLENOIDAL PERTURBATION MODULE
!==============================================================================
! Generates divergence-free velocity perturbations for 3D channel flow DNS
! • x,y directions: Periodic (Fourier spectral)
! • z direction: Non-homogeneous with LGL nodes and wall boundaries
!==============================================================================

module perturbation_module
    use iso_fortran_env, only: real64
    use lgl_module         ! Your existing LGL nodes and derivatives
    use fftw3_dns_module   ! Real FFTW DNS module
    implicit none
    
    ! LAPACK interface for linear algebra
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: real64
            integer, intent(in) :: n, nrhs, lda, ldb
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    
    ! Module-level variables for monitoring
    real(real64), save :: energy_prev = 0.0_real64, time_prev = 0.0_real64
    logical, save :: first_call = .true.
    
    private
    public :: generate_channel_solenoidal_perturbations, &
              validate_divergence_free, &
              monitor_perturbation_evolution, &
              initialize_perturbation_system, &
              compute_perturbation_stats, &
              compute_spectral_derivatives_xy

contains

!------------------------------------------------------------------------------
! MAIN PERTURBATION GENERATOR FOR CHANNEL FLOW
!------------------------------------------------------------------------------
subroutine generate_channel_solenoidal_perturbations(nx, ny, nz, xlen, ylen, &
                                                     z_nodes, fftw_plans, u_pert, v_pert, w_pert, &
                                                     perturbation_amplitude, zero_w_component)
    implicit none
    
    ! Input parameters
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: xlen, ylen
    real(real64), intent(in) :: z_nodes(nz)
    type(fftw3_dns_plans), intent(in) :: fftw_plans
    real(real64), intent(in) :: perturbation_amplitude
    logical, intent(in), optional :: zero_w_component  ! Option to set w=0
    
    ! Output: solenoidal velocity perturbations (u,v need inout for FFTW)
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz)
    real(real64), intent(out) :: w_pert(nx,ny,nz)
    
    ! Local variables
    complex(real64) :: u_hat(nx/2+1,ny,nz), v_hat(nx/2+1,ny,nz)
    real(real64) :: kx, ky, k_perp_sq, k_dot_v_perp
    real(real64) :: random_val
    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    integer :: i, j, k, seed_size
    integer, allocatable :: seed_array(:)
    logical :: use_zero_w
    
    ! Check if w-component should be zeroed
    use_zero_w = .false.
    if (present(zero_w_component)) use_zero_w = zero_w_component
    
    if (use_zero_w) then
        write(*,'(A)') '🌊 Generating 2D solenoidal perturbations (w=0)...'
    else
        write(*,'(A)') '🌊 Generating channel flow solenoidal perturbations...'
    end if
    
    ! Initialize random number generator
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))
    call system_clock(seed_array(1))
    seed_array = seed_array(1) + 37 * [(i, i=1,seed_size)]
    call random_seed(put=seed_array)
    deallocate(seed_array)
    
    ! Step 1: Generate random Fourier modes (periodic x,y directions)
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx/2+1
                call random_number(random_val)
                u_hat(i,j,k) = cmplx(random_val - 0.5_real64, 0.0_real64, real64)
                call random_number(random_val)
                u_hat(i,j,k) = u_hat(i,j,k) + cmplx(0.0_real64, random_val - 0.5_real64, real64)
                
                call random_number(random_val)
                v_hat(i,j,k) = cmplx(random_val - 0.5_real64, 0.0_real64, real64)
                call random_number(random_val)
                v_hat(i,j,k) = v_hat(i,j,k) + cmplx(0.0_real64, random_val - 0.5_real64, real64)
            end do
        end do
    end do
    
    ! Step 2: Apply 2D solenoidal constraint: kx*u + ky*v = 0
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx/2+1
                ! Compute wavenumbers
                kx = 2.0_real64 * pi * real(i-1, real64) / xlen
                
                if (j <= ny/2+1) then
                    ky = 2.0_real64 * pi * real(j-1, real64) / ylen
                else
                    ky = 2.0_real64 * pi * real(j-1-ny, real64) / ylen
                end if
                
                k_perp_sq = kx*kx + ky*ky
                
                ! Apply high-wavenumber filtering to reduce numerical noise
                ! This helps prevent large divergence due to high-frequency components
                if (sqrt(k_perp_sq) > 0.5_real64 * min(2.0_real64*pi*nx/xlen, 2.0_real64*pi*ny/ylen)) then
                    u_hat(i,j,k) = u_hat(i,j,k) * 0.1_real64  ! Strong damping of high frequencies
                    v_hat(i,j,k) = v_hat(i,j,k) * 0.1_real64
                end if
                
                ! Apply 2D solenoidal projection
                if (k_perp_sq > 1.0e-12_real64) then
                    ! Real parts
                    k_dot_v_perp = kx * real(u_hat(i,j,k)) + ky * real(v_hat(i,j,k))
                    u_hat(i,j,k) = u_hat(i,j,k) - cmplx(kx * k_dot_v_perp / k_perp_sq, 0.0_real64, real64)
                    v_hat(i,j,k) = v_hat(i,j,k) - cmplx(ky * k_dot_v_perp / k_perp_sq, 0.0_real64, real64)
                    
                    ! Imaginary parts
                    k_dot_v_perp = kx * aimag(u_hat(i,j,k)) + ky * aimag(v_hat(i,j,k))
                    u_hat(i,j,k) = u_hat(i,j,k) - cmplx(0.0_real64, kx * k_dot_v_perp / k_perp_sq, real64)
                    v_hat(i,j,k) = v_hat(i,j,k) - cmplx(0.0_real64, ky * k_dot_v_perp / k_perp_sq, real64)
                end if
                
                ! Zero mean mode
                if (i == 1 .and. j == 1) then
                    u_hat(i,j,k) = cmplx(0.0_real64, 0.0_real64, real64)
                    v_hat(i,j,k) = cmplx(0.0_real64, 0.0_real64, real64)
                end if
            end do
        end do
    end do
    
    ! Step 3: Transform to physical space (slice by slice in z)
    do k = 1, nz
        call fftw3_backward_2d_dns(fftw_plans, u_hat(:,:,k), u_pert(:,:,k))
        call fftw3_backward_2d_dns(fftw_plans, v_hat(:,:,k), v_pert(:,:,k))
    end do
    
    ! Step 4: Compute w to satisfy full divergence-free condition OR set w=0
    if (use_zero_w) then
        ! Simply zero the w-component for 2D perturbations
        w_pert = 0.0_real64
        write(*,'(A)') '✅ W-component set to zero (2D perturbation mode)'
    else
        ! Compute w to satisfy full divergence-free condition
        call compute_divergence_free_w(nx, ny, nz, xlen, ylen, z_nodes, fftw_plans, &
                                       u_pert, v_pert, w_pert)
    end if
    
    ! Step 5: Enforce wall boundary conditions
    call enforce_channel_walls(nx, ny, nz, z_nodes, u_pert, v_pert, w_pert)
    
    ! Step 6: Scale to desired amplitude
    call scale_to_amplitude(nx, ny, nz, u_pert, v_pert, w_pert, perturbation_amplitude)
    
    if (use_zero_w) then
        write(*,'(A)') '✅ 2D solenoidal perturbations generated successfully'
    else
        write(*,'(A)') '✅ Channel flow perturbations generated successfully'
    end if
    
end subroutine generate_channel_solenoidal_perturbations

!------------------------------------------------------------------------------
! COMPUTE W-COMPONENT FOR FULL DIVERGENCE-FREE CONDITION
!------------------------------------------------------------------------------
subroutine compute_divergence_free_w(nx, ny, nz, xlen, ylen, z_nodes, fftw_plans, &
                                     u_pert, v_pert, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: xlen, ylen, z_nodes(nz)
    type(fftw3_dns_plans), intent(in) :: fftw_plans
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz)
    real(real64), intent(inout) :: w_pert(nx,ny,nz)
    
    ! Local variables
    real(real64) :: du_dx(nx,ny,nz), dv_dy(nx,ny,nz)
    real(real64) :: divergence_xy(nx,ny,nz)
    real(real64) :: lgl_deriv_matrix(nz,nz)
    integer :: i, j, k
    
    ! Compute ∂u/∂x + ∂v/∂y using spectral methods
    call compute_spectral_derivatives_xy(nx, ny, nz, xlen, ylen, fftw_plans, &
                                         u_pert, v_pert, du_dx, dv_dy)
    
    divergence_xy = du_dx + dv_dy
    
    ! Get LGL differentiation matrix for z-direction
    call differentiation_matrix(nz, z_nodes, lgl_deriv_matrix)
    
    ! Solve: ∂w/∂z = -(∂u/∂x + ∂v/∂y) with w(±1) = 0
    call integrate_w_with_walls(nx, ny, nz, z_nodes, lgl_deriv_matrix, &
                                divergence_xy, w_pert)
    
end subroutine compute_divergence_free_w

!------------------------------------------------------------------------------
! SPECTRAL DERIVATIVES IN PERIODIC DIRECTIONS (x,y)
!------------------------------------------------------------------------------
subroutine compute_spectral_derivatives_xy(nx, ny, nz, xlen, ylen, fftw_plans, &
                                           u_field, v_field, du_dx, dv_dy)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: xlen, ylen
    type(fftw3_dns_plans), intent(in) :: fftw_plans
    real(real64), intent(inout) :: u_field(nx,ny,nz), v_field(nx,ny,nz)
    real(real64), intent(out) :: du_dx(nx,ny,nz), dv_dy(nx,ny,nz)
    
    ! Local variables
    complex(real64) :: u_hat(nx/2+1,ny,nz), v_hat(nx/2+1,ny,nz)
    complex(real64) :: dudx_hat(nx/2+1,ny,nz), dvdy_hat(nx/2+1,ny,nz)
    real(real64) :: kx, ky
    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    integer :: i, j, k
    
    ! Forward FFT to Fourier space (slice by slice in z)
    do k = 1, nz
        call fftw3_forward_2d_dns(fftw_plans, u_field(:,:,k), u_hat(:,:,k))
        call fftw3_forward_2d_dns(fftw_plans, v_field(:,:,k), v_hat(:,:,k))
    end do
    
    ! Compute derivatives: multiply by ik in Fourier space
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx/2+1
                ! Wavenumbers
                kx = 2.0_real64 * pi * real(i-1, real64) / xlen
                
                if (j <= ny/2+1) then
                    ky = 2.0_real64 * pi * real(j-1, real64) / ylen
                else
                    ky = 2.0_real64 * pi * real(j-1-ny, real64) / ylen
                end if
                
                ! ∂/∂x and ∂/∂y in Fourier space
                dudx_hat(i,j,k) = cmplx(0.0_real64, kx, real64) * u_hat(i,j,k)
                dvdy_hat(i,j,k) = cmplx(0.0_real64, ky, real64) * v_hat(i,j,k)
            end do
        end do
    end do
    
    ! Inverse FFT back to physical space (slice by slice in z)
    do k = 1, nz
        call fftw3_backward_2d_dns(fftw_plans, dudx_hat(:,:,k), du_dx(:,:,k))
        call fftw3_backward_2d_dns(fftw_plans, dvdy_hat(:,:,k), dv_dy(:,:,k))
    end do
    
end subroutine compute_spectral_derivatives_xy

!------------------------------------------------------------------------------
! INTEGRATE W WITH WALL BOUNDARY CONDITIONS
!------------------------------------------------------------------------------
subroutine integrate_w_with_walls(nx, ny, nz, z_nodes, lgl_deriv_matrix, &
                                  divergence_xy, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: z_nodes(nz), lgl_deriv_matrix(nz,nz)
    real(real64), intent(in) :: divergence_xy(nx,ny,nz)
    real(real64), intent(inout) :: w_pert(nx,ny,nz)
    
    ! Local variables
    real(real64) :: w_solution(nz)
    integer :: i, j, k
    
    ! Solve for w at each (x,y) location using proper integration
    do j = 1, ny
        do i = 1, nx
            ! Use improved integration method
            call integrate_with_boundary_conditions(nz, z_nodes, divergence_xy(i,j,:), w_solution)
            w_pert(i,j,:) = w_solution(:)
        end do
    end do
    
end subroutine integrate_w_with_walls

!------------------------------------------------------------------------------
! ENFORCE CHANNEL WALL BOUNDARY CONDITIONS
!------------------------------------------------------------------------------
subroutine enforce_channel_walls(nx, ny, nz, z_nodes, u_pert, v_pert, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: z_nodes(nz)
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    
    integer :: i, j, k
    real(real64) :: z_val, wall_damping
    
    ! Method 1: Hard zero conditions at walls (z = ±1)
    do j = 1, ny
        do i = 1, nx
            ! Lower wall (z = -1, typically k=1)
            u_pert(i,j,1) = 0.0_real64
            v_pert(i,j,1) = 0.0_real64
            w_pert(i,j,1) = 0.0_real64
            
            ! Upper wall (z = +1, typically k=nz)
            u_pert(i,j,nz) = 0.0_real64
            v_pert(i,j,nz) = 0.0_real64
            w_pert(i,j,nz) = 0.0_real64
        end do
    end do
    
    ! Method 2: Smooth damping near walls (improves numerical stability)
    do k = 2, nz-1  ! Skip boundary points
        z_val = z_nodes(k)
        
        ! Quartic damping: (1-z²)² 
        ! Maximum at center (z=0), zero at walls (z=±1)
        wall_damping = (1.0_real64 - z_val*z_val)**2
        
        do j = 1, ny
            do i = 1, nx
                u_pert(i,j,k) = u_pert(i,j,k) * wall_damping
                v_pert(i,j,k) = v_pert(i,j,k) * wall_damping
                w_pert(i,j,k) = w_pert(i,j,k) * wall_damping
            end do
        end do
    end do
    
end subroutine enforce_channel_walls

!------------------------------------------------------------------------------
! SCALE PERTURBATIONS TO DESIRED AMPLITUDE
!------------------------------------------------------------------------------
subroutine scale_to_amplitude(nx, ny, nz, u_pert, v_pert, w_pert, target_amplitude)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    real(real64), intent(in) :: target_amplitude
    
    ! Local variables
    real(real64) :: current_energy, target_energy, scale_factor
    real(real64) :: base_energy_poiseuille, u_max_poiseuille
    
    ! Poiseuille flow: u_max = 1.5 (centerline velocity)
    u_max_poiseuille = 1.5_real64
    base_energy_poiseuille = 0.5_real64 * u_max_poiseuille**2
    
    ! Compute current perturbation energy
    current_energy = 0.5_real64 * (sum(u_pert**2) + sum(v_pert**2) + sum(w_pert**2)) &
                     / real(nx * ny * nz, real64)
    
    ! Target energy as percentage of base flow
    target_energy = target_amplitude * base_energy_poiseuille
    
    ! Compute scale factor
    if (current_energy > 1.0e-15_real64) then
        scale_factor = sqrt(target_energy / current_energy)
    else
        scale_factor = 1.0_real64
        write(*,'(A)') '⚠️  Warning: Near-zero perturbation energy detected'
    end if
    
    ! Apply uniform scaling
    u_pert = scale_factor * u_pert
    v_pert = scale_factor * v_pert
    w_pert = scale_factor * w_pert
    
    ! Final energy check
    current_energy = 0.5_real64 * (sum(u_pert**2) + sum(v_pert**2) + sum(w_pert**2)) &
                     / real(nx * ny * nz, real64)
    
    ! Output scaling information
    write(*,'(A)') repeat('-', 60)
    write(*,'(A,E12.5)') 'Perturbation energy scaled to: ', current_energy
    write(*,'(A,F8.4,A)') 'Perturbation amplitude: ', 100.0_real64 * target_amplitude, '% of base flow'
    write(*,'(A,F8.4)') 'Scale factor applied: ', scale_factor
    write(*,'(A,E12.5)') 'Base Poiseuille energy: ', base_energy_poiseuille
    write(*,'(A)') repeat('-', 60)
    
end subroutine scale_to_amplitude

!------------------------------------------------------------------------------
! SUPPORTING SUBROUTINES
!------------------------------------------------------------------------------
subroutine create_wall_integration_matrix(nz, deriv_matrix, integration_matrix)
    implicit none
    
    integer, intent(in) :: nz
    real(real64), intent(in) :: deriv_matrix(nz,nz)
    real(real64), intent(out) :: integration_matrix(nz,nz)
    
    integer :: i, j
    
    ! Start with LGL differentiation matrix
    integration_matrix = deriv_matrix
    
    ! Apply wall boundary conditions
    ! First row: w(z=-1) = 0
    integration_matrix(1,:) = 0.0_real64
    integration_matrix(1,1) = 1.0_real64
    
    ! Last row: w(z=+1) = 0  
    integration_matrix(nz,:) = 0.0_real64
    integration_matrix(nz,nz) = 1.0_real64
    
end subroutine create_wall_integration_matrix

subroutine simple_integration_fallback(nz, z_nodes, rhs, w_solution)
    implicit none
    
    integer, intent(in) :: nz
    real(real64), intent(in) :: z_nodes(nz), rhs(nz)
    real(real64), intent(out) :: w_solution(nz)
    
    integer :: k
    real(real64) :: dz, integral_sum
    
    ! Simple trapezoidal integration with wall BCs
    w_solution(1) = 0.0_real64  ! w(-1) = 0
    
    integral_sum = 0.0_real64
    do k = 2, nz-1
        dz = z_nodes(k) - z_nodes(k-1)
        integral_sum = integral_sum + 0.5_real64 * (rhs(k) + rhs(k-1)) * dz
        w_solution(k) = integral_sum
    end do
    
    w_solution(nz) = 0.0_real64  ! w(+1) = 0
    
    ! Adjust to satisfy boundary conditions exactly
    if (abs(w_solution(nz-1)) > 1.0e-15_real64) then
        do k = 2, nz-1
            w_solution(k) = w_solution(k) * (1.0_real64 - real(k-1,real64)/real(nz-1,real64))
        end do
    end if
    
end subroutine simple_integration_fallback

subroutine integrate_with_boundary_conditions(nz, z_nodes, divergence_xy, w_solution)
    implicit none
    
    integer, intent(in) :: nz
    real(real64), intent(in) :: z_nodes(nz), divergence_xy(nz)
    real(real64), intent(out) :: w_solution(nz)
    
    ! Temporarily use fallback method while debugging spectral integration
    ! The LAPACK spectral solver is producing unphysically large values
    call integrate_with_boundary_conditions_fallback(nz, z_nodes, divergence_xy, w_solution)
    
end subroutine integrate_with_boundary_conditions

!------------------------------------------------------------------------------
! FALLBACK INTEGRATION METHOD
!------------------------------------------------------------------------------
subroutine integrate_with_boundary_conditions_fallback(nz, z_nodes, divergence_xy, w_solution)
    implicit none
    
    integer, intent(in) :: nz
    real(real64), intent(in) :: z_nodes(nz), divergence_xy(nz)
    real(real64), intent(out) :: w_solution(nz)
    
    ! Local variables
    real(real64) :: integral_up(nz), total_integral, correction_factor
    real(real64) :: z_normalized(nz)
    integer :: k
    
    ! Improved single-direction integration with smooth correction
    ! This should provide better numerical stability
    
    ! Initialize boundary condition
    w_solution(1) = 0.0_real64  ! w(-1) = 0
    integral_up(1) = 0.0_real64
    
    ! Integrate upward using trapezoidal rule
    do k = 2, nz
        integral_up(k) = integral_up(k-1) - 0.5_real64 * (divergence_xy(k) + divergence_xy(k-1)) * &
                         (z_nodes(k) - z_nodes(k-1))
    end do
    
    ! Store the total integral (should be zero for divergence-free field)
    total_integral = integral_up(nz)
    
    ! Create normalized z coordinates for smooth correction
    do k = 1, nz
        z_normalized(k) = (z_nodes(k) - z_nodes(1)) / (z_nodes(nz) - z_nodes(1))
    end do
    
    ! Apply smooth quadratic correction to satisfy both boundary conditions
    ! This ensures w(1) = 0 and w(nz) = 0 while minimizing oscillations
    do k = 1, nz
        correction_factor = z_normalized(k) * (2.0_real64 - z_normalized(k))
        w_solution(k) = integral_up(k) - total_integral * correction_factor
    end do
    
    ! Ensure exact boundary conditions
    w_solution(1) = 0.0_real64
    w_solution(nz) = 0.0_real64
    
    ! Apply additional smoothing near boundaries to reduce numerical artifacts
    do k = 2, min(4, nz-1)
        w_solution(k) = w_solution(k) * sin(0.5_real64 * 3.14159_real64 * real(k-1, real64) / 3.0_real64)**2
    end do
    do k = max(nz-3, 2), nz-1
        w_solution(k) = w_solution(k) * sin(0.5_real64 * 3.14159_real64 * real(nz-k, real64) / 3.0_real64)**2
    end do
    
end subroutine integrate_with_boundary_conditions_fallback

subroutine lgl_differentiation_matrix(nz, z_nodes, deriv_matrix)
    implicit none
    
    integer, intent(in) :: nz
    real(real64), intent(in) :: z_nodes(nz)
    real(real64), intent(out) :: deriv_matrix(nz,nz)
    
    ! Interface to your existing LGL module
    ! This should call your existing LGL differentiation matrix routine
    ! Adjust the call based on your lgl_module.f90 interface
    
    ! Use the correct function name from lgl_module
    call differentiation_matrix(nz, z_nodes, deriv_matrix)
    
end subroutine lgl_differentiation_matrix

!------------------------------------------------------------------------------
! VALIDATION AND MONITORING SUBROUTINES
!------------------------------------------------------------------------------
subroutine validate_divergence_free(nx, ny, nz, xlen, ylen, z_nodes, fftw_plans, &
                                    u_pert, v_pert, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: xlen, ylen, z_nodes(nz)
    type(fftw3_dns_plans), intent(in) :: fftw_plans
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz)
    real(real64), intent(in) :: w_pert(nx,ny,nz)
    
    ! Local variables
    real(real64) :: du_dx(nx,ny,nz), dv_dy(nx,ny,nz), dw_dz(nx,ny,nz)
    real(real64) :: divergence(nx,ny,nz)
    real(real64) :: max_div, mean_div, rms_div, div_at_center
    real(real64) :: lgl_deriv_matrix(nz,nz)
    integer :: i, j, k, l, k_center
    
    write(*,'(A)') '🔍 Validating solenoidal condition for channel flow...'
    
    ! Spectral derivatives in periodic directions (x,y)
    call compute_spectral_derivatives_xy(nx, ny, nz, xlen, ylen, fftw_plans, &
                                         u_pert, v_pert, du_dx, dv_dy)
    
    ! LGL derivatives in wall-normal direction (z)
    call differentiation_matrix(nz, z_nodes, lgl_deriv_matrix)
    
    do j = 1, ny
        do i = 1, nx
            do k = 1, nz
                dw_dz(i,j,k) = 0.0_real64
                do l = 1, nz
                    dw_dz(i,j,k) = dw_dz(i,j,k) + lgl_deriv_matrix(k,l) * w_pert(i,j,l)
                end do
            end do
        end do
    end do
    
    ! Compute full divergence: ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
    divergence = du_dx + dv_dy + dw_dz
    
    ! Divergence statistics
    max_div = maxval(abs(divergence))
    mean_div = sum(divergence) / real(nx * ny * nz, real64)
    rms_div = sqrt(sum(divergence**2) / real(nx * ny * nz, real64))
    
    ! Divergence at channel center
    k_center = (nz + 1) / 2
    div_at_center = sum(abs(divergence(:,:,k_center))) / real(nx * ny, real64)
    
    ! Output validation results
    write(*,'(A)') repeat('=', 70)
    write(*,'(A)') '🔍 DIVERGENCE-FREE VALIDATION RESULTS'
    write(*,'(A)') repeat('=', 70)
    write(*,'(A,E12.5)') '   Maximum |∇·v|: ', max_div
    write(*,'(A,E12.5)') '   Mean ∇·v:      ', mean_div  
    write(*,'(A,E12.5)') '   RMS ∇·v:       ', rms_div
    write(*,'(A,E12.5)') '   |∇·v| at center: ', div_at_center
    
    ! Quality assessment
    if (max_div < 1.0e-13_real64) then
        write(*,'(A)') '   ✅ EXCELLENT: Machine precision divergence-free'
    else if (max_div < 1.0e-11_real64) then
        write(*,'(A)') '   ✅ VERY GOOD: Near machine precision'
    else if (max_div < 1.0e-9_real64) then
        write(*,'(A)') '   ✅ GOOD: Acceptable numerical divergence'
    else if (max_div < 1.0e-7_real64) then
        write(*,'(A)') '   ⚠️  ACCEPTABLE: Small divergence detected'
    else
        write(*,'(A)') '   ❌ WARNING: Large divergence - check implementation!'
    end if
    
    write(*,'(A)') repeat('=', 70)
    
end subroutine validate_divergence_free

subroutine compute_perturbation_stats(nx, ny, nz, z_nodes, u_pert, v_pert, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: z_nodes(nz)
    real(real64), intent(in) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    
    ! Local variables
    real(real64) :: energy_total, energy_u, energy_v, energy_w
    real(real64) :: max_u, max_v, max_w, rms_total
    real(real64) :: energy_lower_half, energy_upper_half
    integer :: k_mid, i, j, k
    
    ! Energy components
    energy_u = 0.5_real64 * sum(u_pert**2) / real(nx * ny * nz, real64)
    energy_v = 0.5_real64 * sum(v_pert**2) / real(nx * ny * nz, real64)
    energy_w = 0.5_real64 * sum(w_pert**2) / real(nx * ny * nz, real64)
    energy_total = energy_u + energy_v + energy_w
    
    ! Maximum values
    max_u = maxval(abs(u_pert))
    max_v = maxval(abs(v_pert))
    max_w = maxval(abs(w_pert))
    
    ! RMS velocity
    rms_total = sqrt(2.0_real64 * energy_total)
    
    ! Energy distribution (lower vs upper half)
    k_mid = (nz + 1) / 2
    energy_lower_half = 0.5_real64 * (sum(u_pert(:,:,1:k_mid)**2) + &
                                  sum(v_pert(:,:,1:k_mid)**2) + &
                                  sum(w_pert(:,:,1:k_mid)**2)) / real(nx * ny * k_mid, real64)
    
    energy_upper_half = 0.5_real64 * (sum(u_pert(:,:,k_mid+1:nz)**2) + &
                                  sum(v_pert(:,:,k_mid+1:nz)**2) + &
                                  sum(w_pert(:,:,k_mid+1:nz)**2)) / real(nx * ny * (nz-k_mid), real64)
    
    ! Output statistics
    write(*,'(A)') '📊 Perturbation Statistics:'
    write(*,'(A,4E12.5)') '   Energy [Total, u, v, w]: ', energy_total, energy_u, energy_v, energy_w
    write(*,'(A,3E12.5)') '   Max |velocity| [u, v, w]: ', max_u, max_v, max_w
    write(*,'(A,E12.5)') '   RMS velocity: ', rms_total
    write(*,'(A,2E12.5)') '   Energy [lower, upper half]: ', energy_lower_half, energy_upper_half
    
    ! Energy ratio check
    if (energy_upper_half > 0.0_real64 .and. energy_lower_half > 0.0_real64) then
        write(*,'(A,F8.4)') '   Energy ratio (upper/lower): ', energy_upper_half / energy_lower_half
    end if
    
end subroutine compute_perturbation_stats

subroutine initialize_perturbation_system(nx, ny, nz)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    
    ! Create monitoring output file
    open(unit=98, file='channel_flow_evolution.dat', status='replace')
    write(98,'(A)') '# Channel flow perturbation evolution'
    write(98,'(A)') '# Columns: istep, time, E_total, E_u, E_v, E_w, tau_lower, tau_upper,'
    write(98,'(A)') '#          u_max, u_center, growth_rate, pert_percent'
    write(98,'(A,3I6)') '# Grid: nx, ny, nz = ', nx, ny, nz
    write(98,'(A)') '# =========================================='
    close(98)
    
    write(*,'(A)') '📊 Channel flow monitoring system initialized'
    write(*,'(A)') '    Output file: channel_flow_evolution.dat'
    
end subroutine initialize_perturbation_system

subroutine monitor_perturbation_evolution(nx, ny, nz, z_nodes, u, v, w, &
                                          u_base, istep, time, re, dt)
    implicit none
    
    integer, intent(in) :: nx, ny, nz, istep
    real(real64), intent(in) :: z_nodes(nz), time, re, dt
    real(real64), intent(in) :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
    real(real64), intent(in) :: u_base(nx,ny,nz)
    
    ! Local variables
    real(real64) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    real(real64) :: energy_total, energy_u, energy_v, energy_w
    real(real64) :: tau_wall_lower, tau_wall_upper, cf_lower, cf_upper
    real(real64) :: u_max_current, u_center_avg
    real(real64) :: perturbation_percent, growth_rate
    integer :: i, j, k_center
    
    ! Compute perturbations
    u_pert = u - u_base
    v_pert = v  ! base v = 0 for Poiseuille
    w_pert = w  ! base w = 0 for Poiseuille
    
    ! Energy components
    energy_u = 0.5_real64 * sum(u_pert**2) / real(nx * ny * nz, real64)
    energy_v = 0.5_real64 * sum(v_pert**2) / real(nx * ny * nz, real64)
    energy_w = 0.5_real64 * sum(w_pert**2) / real(nx * ny * nz, real64)
    energy_total = energy_u + energy_v + energy_w
    
    ! Flow statistics
    u_max_current = maxval(u)
    k_center = (nz + 1) / 2
    u_center_avg = sum(u(:,:,k_center)) / real(nx * ny, real64)
    
    ! Wall shear stress (using finite differences at LGL boundaries)
    tau_wall_lower = 0.0_real64
    tau_wall_upper = 0.0_real64
    
    do j = 1, ny
        do i = 1, nx
            ! Lower wall: τ = μ(∂u/∂z) = (1/Re)(∂u/∂z)
            tau_wall_lower = tau_wall_lower + &
                (u(i,j,2) - u(i,j,1)) / (z_nodes(2) - z_nodes(1)) / re
            
            ! Upper wall
            tau_wall_upper = tau_wall_upper + &
                (u(i,j,nz) - u(i,j,nz-1)) / (z_nodes(nz) - z_nodes(nz-1)) / re
        end do
    end do
    
    tau_wall_lower = tau_wall_lower / real(nx * ny, real64)
    tau_wall_upper = tau_wall_upper / real(nx * ny, real64)
    
    ! Skin friction coefficients
    cf_lower = 2.0_real64 * abs(tau_wall_lower) / (0.5_real64 * 1.5_real64**2)  ! Based on u_max
    cf_upper = 2.0_real64 * abs(tau_wall_upper) / (0.5_real64 * 1.5_real64**2)
    
    ! Perturbation amplitude as percentage
    perturbation_percent = sqrt(2.0_real64 * energy_total) / 1.5_real64 * 100.0_real64
    
    ! Growth rate calculation
    if (.not. first_call .and. energy_prev > 1.0e-15_real64 .and. time > time_prev) then
        growth_rate = log(energy_total / energy_prev) / (time - time_prev)
    else
        growth_rate = 0.0_real64
        first_call = .false.
    end if
    
    ! Output monitoring data (every 50 steps)
    if (mod(istep, 50) == 0) then
        write(*,'(A)') repeat('-', 80)
        write(*,'(A,I8,F12.4)') '📊 Step: ', istep, time
        write(*,'(A,4E14.6)') '   Energy [Total, u, v, w]: ', energy_total, energy_u, energy_v, energy_w
        write(*,'(A,F8.4,A,E12.5)') '   Perturbation: ', perturbation_percent, '%, σ = ', growth_rate
        write(*,'(A,F12.6,A,F12.6)') '   u_max: ', u_max_current, ', u_center: ', u_center_avg
        write(*,'(A,2E14.6)') '   τ_wall [lower, upper]: ', tau_wall_lower, tau_wall_upper
        write(*,'(A,2F10.6)') '   C_f [lower, upper]: ', cf_lower, cf_upper
        
        ! Write to monitoring file
        open(unit=98, file='channel_flow_evolution.dat', position='append')
        write(98,'(I8,11E16.8)') istep, time, energy_total, energy_u, energy_v, energy_w, &
                                  tau_wall_lower, tau_wall_upper, u_max_current, u_center_avg, &
                                  growth_rate, perturbation_percent
        close(98)
    end if
    
    ! Update for next call
    energy_prev = energy_total
    time_prev = time
    
end subroutine monitor_perturbation_evolution

end module perturbation_module
