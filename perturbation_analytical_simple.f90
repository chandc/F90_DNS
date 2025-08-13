!==============================================================================
! SIMPLE ANALYTICAL PERTURBATION MODULE
!==============================================================================
! Generates exactly divergence-free analytical velocity perturbations
! Based on stream function approach for guaranteed solver stability
!==============================================================================

module perturbation_analytical_simple
    use iso_fortran_env, only: real64
    implicit none
    
    private
    public :: generate_analytical_perturbations, &
              validate_analytical_divergence

contains

!------------------------------------------------------------------------------
! SIMPLE ANALYTICAL PERTURBATION GENERATOR
!------------------------------------------------------------------------------
! This implements exactly divergence-free perturbations using analytical
! expressions. Perfect for solver stability testing and initial implementation.
!------------------------------------------------------------------------------
subroutine generate_analytical_perturbations(nx, ny, nz, xlen, ylen, &
                                            z_nodes, u_pert, v_pert, w_pert, &
                                            perturbation_amplitude, mode_x, mode_y)
    implicit none
    
    ! Input parameters
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: xlen, ylen
    real(real64), intent(in) :: z_nodes(nz)
    real(real64), intent(in) :: perturbation_amplitude
    integer, intent(in), optional :: mode_x, mode_y  ! Wavenumber modes
    
    ! Output: exactly divergence-free velocity perturbations
    real(real64), intent(out) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    
    ! Local variables
    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    real(real64) :: x, y, z
    real(real64) :: kx, ky
    real(real64) :: psi, phi  ! Stream functions
    real(real64) :: z_factor, wall_damping, z_profile
    integer :: i, j, k
    integer :: mx, my  ! Local mode numbers
    
    ! Default to simple low modes for stability
    mx = 1
    my = 1
    if (present(mode_x)) mx = mode_x
    if (present(mode_y)) my = mode_y
    
    ! Wavenumbers
    kx = 2.0_real64 * pi * real(mx, real64) / xlen
    ky = 2.0_real64 * pi * real(my, real64) / ylen
    
    write(*,'(A)') '🎯 Generating simple analytical perturbations...'
    write(*,'(A,I0,A,I0)') '   Modes: mx=', mx, ', my=', my
    write(*,'(A,F8.4)') '   Amplitude: ', perturbation_amplitude
    
    ! Initialize arrays
    u_pert = 0.0_real64
    v_pert = 0.0_real64
    w_pert = 0.0_real64
    
    ! Generate exactly divergence-free perturbations using stream function approach
    do k = 1, nz
        z = z_nodes(k)
        
        ! Wall damping factor (goes to zero at walls z=±1)
        wall_damping = 1.0_real64 - z*z  ! Simple parabolic damping
        
        do j = 1, ny
            y = real(j-1, real64) * ylen / real(ny, real64)
            do i = 1, nx
                x = real(i-1, real64) * xlen / real(nx, real64)
                
                ! RETURN TO BEST METHOD: Stream function approach for perfect 2D divergence-free
                psi = sin(kx * x) * sin(ky * y) * wall_damping
                
                ! Velocity from stream function: u = -∂ψ/∂y, v = ∂ψ/∂x
                u_pert(i,j,k) = -ky * sin(kx * x) * cos(ky * y) * wall_damping
                v_pert(i,j,k) =  kx * cos(kx * x) * sin(ky * y) * wall_damping
                
                ! Set w = 0 for exactly divergence-free 2D perturbations
                w_pert(i,j,k) = 0.0_real64
                
            end do
        end do
    end do
    
    ! Scale to desired amplitude
    call scale_perturbations(nx, ny, nz, u_pert, v_pert, w_pert, perturbation_amplitude)
    
    ! Ensure wall boundary conditions
    call enforce_wall_boundaries(nx, ny, nz, u_pert, v_pert, w_pert)
    
    write(*,'(A)') '✅ Analytical perturbations generated successfully'
    
    ! Validate divergence
    call validate_analytical_divergence(nx, ny, nz, xlen, ylen, z_nodes, &
                                        u_pert, v_pert, w_pert)
    
end subroutine generate_analytical_perturbations

!------------------------------------------------------------------------------
! SCALE PERTURBATIONS TO DESIRED AMPLITUDE
!------------------------------------------------------------------------------
subroutine scale_perturbations(nx, ny, nz, u_pert, v_pert, w_pert, target_amplitude)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    real(real64), intent(in) :: target_amplitude
    
    real(real64) :: current_max, scale_factor
    
    ! Find maximum velocity magnitude
    current_max = maxval(sqrt(u_pert**2 + v_pert**2 + w_pert**2))
    
    if (current_max > 1.0e-12_real64) then
        scale_factor = target_amplitude / current_max
        u_pert = u_pert * scale_factor
        v_pert = v_pert * scale_factor
        w_pert = w_pert * scale_factor
        
        write(*,'(A,ES12.4,A,ES12.4)') '   Scaled from max=', current_max, ' to max=', target_amplitude
    else
        write(*,'(A)') '   Warning: Zero perturbation field detected'
    end if
    
end subroutine scale_perturbations

!------------------------------------------------------------------------------
! ENFORCE WALL BOUNDARY CONDITIONS
!------------------------------------------------------------------------------
subroutine enforce_wall_boundaries(nx, ny, nz, u_pert, v_pert, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(inout) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    
    ! Bottom wall (k=1) and top wall (k=nz): all velocities = 0
    u_pert(:,:,1) = 0.0_real64
    v_pert(:,:,1) = 0.0_real64
    w_pert(:,:,1) = 0.0_real64
    
    u_pert(:,:,nz) = 0.0_real64
    v_pert(:,:,nz) = 0.0_real64
    w_pert(:,:,nz) = 0.0_real64
    
    write(*,'(A)') '   Wall boundary conditions enforced'
    
end subroutine enforce_wall_boundaries

!------------------------------------------------------------------------------
! VALIDATE DIVERGENCE OF ANALYTICAL PERTURBATIONS
!------------------------------------------------------------------------------
subroutine validate_analytical_divergence(nx, ny, nz, xlen, ylen, z_nodes, &
                                          u_pert, v_pert, w_pert)
    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real(real64), intent(in) :: xlen, ylen
    real(real64), intent(in) :: z_nodes(nz)
    real(real64), intent(in) :: u_pert(nx,ny,nz), v_pert(nx,ny,nz), w_pert(nx,ny,nz)
    
    real(real64) :: dx, dy
    real(real64) :: dudx, dvdy, dwdz, divergence
    real(real64) :: max_div, rms_div, total_div_sq
    integer :: i, j, k, count
    
    dx = xlen / real(nx, real64)
    dy = ylen / real(ny, real64)
    
    max_div = 0.0_real64
    total_div_sq = 0.0_real64
    count = 0
    
    ! Check divergence at interior points
    do k = 2, nz-1
        do j = 2, ny-1
            do i = 2, nx-1
                
                ! Compute derivatives using finite differences
                dudx = (u_pert(i+1,j,k) - u_pert(i-1,j,k)) / (2.0_real64 * dx)
                dvdy = (v_pert(i,j+1,k) - v_pert(i,j-1,k)) / (2.0_real64 * dy)
                
                ! Z-derivative using LGL spacing
                dwdz = (w_pert(i,j,k+1) - w_pert(i,j,k-1)) / (z_nodes(k+1) - z_nodes(k-1))
                
                divergence = dudx + dvdy + dwdz
                
                max_div = max(max_div, abs(divergence))
                total_div_sq = total_div_sq + divergence**2
                count = count + 1
                
            end do
        end do
    end do
    
    if (count > 0) then
        rms_div = sqrt(total_div_sq / real(count, real64))
        write(*,'(A)') '📊 Analytical Perturbation Divergence Check:'
        write(*,'(A,ES12.4)') '   Maximum divergence: ', max_div
        write(*,'(A,ES12.4)') '   RMS divergence:     ', rms_div
        
        if (max_div < 1.0e-12_real64) then
            write(*,'(A)') '✅ Perfect divergence-free property (machine precision)'
        else if (max_div < 1.0e-10_real64) then
            write(*,'(A)') '✅ Excellent divergence-free property'
        else if (max_div < 1.0e-6_real64) then
            write(*,'(A)') '✅ Good divergence-free property'
        else
            write(*,'(A)') '⚠️  Warning: High divergence detected - numerical instability likely!'
        end if
    end if
    
end subroutine validate_analytical_divergence

end module perturbation_analytical_simple
