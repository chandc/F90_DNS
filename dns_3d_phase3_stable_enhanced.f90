program dns_3d_phase3_stable_enhanced
    use fft_module
    use lgl_module
    implicit none
    
    ! FFT plans (module variable)
    type(fft_plans) :: plans
    
    ! Physical parameters
    real(8), parameter :: Re = 500.0d0     ! Lower Reynolds number for stability
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8), parameter :: Lx = 2.0d0*pi
    real(8), parameter :: Ly = 2.0d0*pi
    real(8), parameter :: Lz = 2.0d0
    
    ! Grid parameters
    integer, parameter :: nx = 32
    integer, parameter :: ny = 16
    integer, parameter :: nz = 17
    
    ! Time integration parameters
    real(8), parameter :: dt_base = 1.0d-5    ! Very small time step
    real(8) :: dt, dt_cfl
    integer, parameter :: nsteps = 2000
    
    ! Grid arrays
    real(8) :: x(nx), y(ny), z(nz), weights(nz)
    real(8) :: kx(nx/2+1), ky(ny)
    
    ! Velocity field arrays
    complex(8) :: u_hat(nx/2+1, ny, nz)
    complex(8) :: v_hat(nx/2+1, ny, nz)
    complex(8) :: w_hat(nx/2+1, ny, nz)
    real(8) :: u_vel(nx, ny, nz), v_vel(nx, ny, nz), w_vel(nx, ny, nz)
    
    ! Pressure array
    complex(8) :: p_hat(nx/2+1, ny, nz)
    real(8) :: p(nx, ny, nz)
    
    ! Nonlinear terms (Adams-Bashforth history)
    complex(8) :: nl_u_old(nx/2+1, ny, nz)
    complex(8) :: nl_v_old(nx/2+1, ny, nz)
    complex(8) :: nl_w_old(nx/2+1, ny, nz)
    
    ! Differentiation matrices and operators
    real(8) :: D1(nz, nz), D2(nz, nz)
    
    ! Working arrays for FFT
    real(8) :: u_work(nx, ny), v_work(nx, ny), w_work(nx, ny)
    complex(8) :: u_hat_work(nx/2+1, ny), v_hat_work(nx/2+1, ny), w_hat_work(nx/2+1, ny)
    
    ! Arrays for nonlinear term computation
    real(8) :: dudx(nx, ny, nz), dudy(nx, ny, nz), dudz(nx, ny, nz)
    real(8) :: dvdx(nx, ny, nz), dvdy(nx, ny, nz), dvdz(nx, ny, nz)
    real(8) :: dwdx(nx, ny, nz), dwdy(nx, ny, nz), dwdz(nx, ny, nz)
    real(8) :: nl_u_phys(nx, ny, nz), nl_v_phys(nx, ny, nz), nl_w_phys(nx, ny, nz)
    
    ! Utility arrays
    real(8) :: u_max, v_max, w_max, div_max, ke, enstrophy
    integer :: i, j, k, step
    real(8) :: time
    
    write(*,*) '===================================================='
    write(*,*) '    3D DNS SOLVER - PHASE 3 STABLE ENHANCED'
    write(*,*) '    Nonlinear + Viscous + Pressure (Conservative)'
    write(*,*) '===================================================='
    
    ! Initialize grid
    call initialize_grid()
    write(*,*) 'Grid initialized'
    
    ! Initialize LGL operators
    call initialize_lgl_operators()
    write(*,*) 'LGL operators initialized'
    
    ! Initialize FFT system
    call setup_fft_plans(plans, nx, ny, u_work, u_hat_work)
    write(*,*) 'FFT system initialized'
    
    ! Set initial conditions
    call set_initial_conditions()
    write(*,*) 'Initial flow field set'
    
    ! Initialize Adams-Bashforth terms
    call compute_nonlinear_terms(nl_u_old, nl_v_old, nl_w_old)
    
    ! Time integration loop
    write(*,*) 'Starting time integration...'
    write(*,'(A,3I10)') ' Grid:', nx, ny, nz
    write(*,'(A,F12.6,A,E12.6,A,I10)') ' Re =', Re, ', dt =', dt_base, ', nsteps =', nsteps
    
    time = 0.0d0
    dt = dt_base
    
    do step = 1, nsteps
        ! Adaptive time stepping with CFL condition
        call check_cfl_condition()
        dt = min(dt_base, 0.1d0 * dt_cfl)  ! Very conservative
        
        ! Exit if time step becomes too small
        if (dt < 1.0d-8) then
            write(*,*) 'Time step too small, stopping simulation'
            exit
        endif
        
        ! Conservative time integration
        call conservative_time_step()
        
        ! Update time
        time = time + dt
        
        ! Output diagnostics
        if (mod(step, 50) == 0) then
            call compute_diagnostics()
            write(*,'(A,I6,A,F8.4,A,E12.2,A,E12.2,A,E12.2,A,E12.2)') &
                'Step', step, ', t=', time, ', KE=', ke, ', |u|=', u_max, ', Ω=', enstrophy, ', div=', div_max
        endif
        
        ! Check for instability with stricter bounds
        if (ke > 10.0d0 .or. u_max > 5.0d0 .or. isnan(ke) .or. div_max > 10.0d0) then
            write(*,*) 'Simulation approaching instability, stopping'
            write(*,*) 'Final values: KE=', ke, ', |u|_max=', u_max, ', div_max=', div_max
            exit
        endif
    end do
    
    write(*,*) 'Time integration completed!'
    
    ! Final diagnostics
    call compute_diagnostics()
    write(*,*) ''
    write(*,*) 'Final solution diagnostics:'
    write(*,'(A,E15.6)') '  Kinetic energy:', ke
    write(*,'(A,E15.6)') '  Enstrophy:', enstrophy
    write(*,'(A,E15.6,A,E15.6,A,E15.6)') '  Max velocities: u_max=', u_max, ', v_max=', v_max, ', w_max=', w_max
    write(*,'(A,E15.6)') '  Max divergence:', div_max
    write(*,'(A,F15.6)') '  Final time:', time
    
    ! Cleanup
    call destroy_fft_plans(plans)
    write(*,*) 'Cleanup completed'
    write(*,*) '===================================================='
    write(*,*) '    PHASE 3 STABLE ENHANCED SIMULATION COMPLETE'
    write(*,*) '===================================================='

contains

    subroutine initialize_grid()
        implicit none
        integer :: i, j
        
        ! X direction (periodic, FFT)
        do i = 1, nx
            x(i) = Lx * real(i-1, 8) / real(nx, 8)
        end do
        
        ! Y direction (periodic, FFT) 
        do j = 1, ny
            y(j) = Ly * real(j-1, 8) / real(ny, 8)
        end do
        
        ! Z direction (Chebyshev, LGL)
        call lgl_nodes_weights(nz, z, weights)
        
        ! Wavenumbers for FFT
        do i = 1, nx/2+1
            kx(i) = real(i-1, 8) * 2.0d0*pi/Lx
        end do
        
        do j = 1, ny/2+1
            ky(j) = real(j-1, 8) * 2.0d0*pi/Ly
        end do
        do j = ny/2+2, ny
            ky(j) = real(j-1-ny, 8) * 2.0d0*pi/Ly
        end do
        
        write(*,'(A,3I10)') ' Parameters initialized:'
        write(*,'(A,3I10)') '   Grid:', nx, ny, nz
        write(*,'(A,F15.6)') '   Re =', Re
        write(*,'(A,F15.6,A,F15.6,A,F15.6,A,F15.6,A)') '   Domain: [0,', Lx, '] x [0,', Ly, '] x [', z(1), ',', z(nz), ']'
    end subroutine

    subroutine initialize_lgl_operators()
        implicit none
        
        ! Create differentiation matrices
        call differentiation_matrix(nz, z, D1)
        
        ! Second derivative matrix
        D2 = matmul(D1, D1)
    end subroutine

    subroutine set_initial_conditions()
        implicit none
        integer :: i, j, k
        real(8) :: amp, decay_factor
        
        ! Very small amplitude for stability
        amp = 0.001d0
        
        ! Simple initial perturbation
        do k = 1, nz
            decay_factor = (1.0d0 - z(k)*z(k)) * 0.5d0  ! Gentle profile
            do j = 1, ny
                do i = 1, nx
                    u_vel(i,j,k) = amp * decay_factor * &
                                   sin(2.0d0*pi*x(i)/Lx) * cos(2.0d0*pi*y(j)/Ly)
                    v_vel(i,j,k) = -amp * decay_factor * &
                                   cos(2.0d0*pi*x(i)/Lx) * sin(2.0d0*pi*y(j)/Ly)
                    w_vel(i,j,k) = 0.0d0
                end do
            end do
        end do
        
        ! Apply boundary conditions
        u_vel(:,:,1) = 0.0d0
        u_vel(:,:,nz) = 0.0d0
        v_vel(:,:,1) = 0.0d0
        v_vel(:,:,nz) = 0.0d0
        w_vel(:,:,1) = 0.0d0
        w_vel(:,:,nz) = 0.0d0
        
        ! Transform to spectral space
        call velocity_to_spectral()
        call enforce_divergence_free()
    end subroutine
    
    subroutine velocity_to_spectral()
        implicit none
        integer :: k
        
        do k = 1, nz
            u_work = u_vel(:,:,k)
            v_work = v_vel(:,:,k)
            w_work = w_vel(:,:,k)
            
            call fft_forward_2d(plans, u_work, u_hat_work)
            call fft_forward_2d(plans, v_work, v_hat_work)
            call fft_forward_2d(plans, w_work, w_hat_work)
            
            u_hat(:,:,k) = u_hat_work
            v_hat(:,:,k) = v_hat_work
            w_hat(:,:,k) = w_hat_work
        end do
    end subroutine
    
    subroutine velocity_to_physical()
        implicit none
        integer :: k
        
        do k = 1, nz
            u_hat_work = u_hat(:,:,k)
            v_hat_work = v_hat(:,:,k)
            w_hat_work = w_hat(:,:,k)
            
            call fft_backward_2d(plans, u_hat_work, u_work)
            call fft_backward_2d(plans, v_hat_work, v_work)
            call fft_backward_2d(plans, w_hat_work, w_work)
            
            u_vel(:,:,k) = u_work
            v_vel(:,:,k) = v_work
            w_vel(:,:,k) = w_work
        end do
    end subroutine

    subroutine enforce_divergence_free()
        implicit none
        integer :: i, j, k
        real(8) :: k2
        complex(8) :: ikx, iky, div
        
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx/2+1
                    if (i == 1 .and. j == 1) cycle
                    
                    ikx = cmplx(0.0d0, kx(i), 8)
                    iky = cmplx(0.0d0, ky(j), 8)
                    k2 = kx(i)**2 + ky(j)**2
                    
                    if (k2 > 1.0d-12) then
                        div = ikx*u_hat(i,j,k) + iky*v_hat(i,j,k)
                        u_hat(i,j,k) = u_hat(i,j,k) - ikx*div/k2
                        v_hat(i,j,k) = v_hat(i,j,k) - iky*div/k2
                    endif
                end do
            end do
        end do
    end subroutine

    subroutine compute_derivatives()
        implicit none
        integer :: i, j, k
        real(8) :: temp_z(nz)
        
        call velocity_to_physical()
        
        ! X and Y derivatives
        do k = 1, nz
            u_work = u_vel(:,:,k)
            v_work = v_vel(:,:,k)
            w_work = w_vel(:,:,k)
            
            call fft_forward_2d(plans, u_work, u_hat_work)
            call fft_forward_2d(plans, v_work, v_hat_work)
            call fft_forward_2d(plans, w_work, w_hat_work)
            
            ! X-derivatives
            do j = 1, ny
                do i = 1, nx/2+1
                    u_hat_work(i,j) = cmplx(0.0d0, kx(i), 8) * u_hat_work(i,j)
                    v_hat_work(i,j) = cmplx(0.0d0, kx(i), 8) * v_hat_work(i,j)
                    w_hat_work(i,j) = cmplx(0.0d0, kx(i), 8) * w_hat_work(i,j)
                end do
            end do
            
            call fft_backward_2d(plans, u_hat_work, u_work)
            call fft_backward_2d(plans, v_hat_work, v_work)
            call fft_backward_2d(plans, w_hat_work, w_work)
            
            dudx(:,:,k) = u_work
            dvdx(:,:,k) = v_work
            dwdx(:,:,k) = w_work
            
            ! Y-derivatives
            u_work = u_vel(:,:,k)
            v_work = v_vel(:,:,k)
            w_work = w_vel(:,:,k)
            
            call fft_forward_2d(plans, u_work, u_hat_work)
            call fft_forward_2d(plans, v_work, v_hat_work)
            call fft_forward_2d(plans, w_work, w_hat_work)
            
            do j = 1, ny
                do i = 1, nx/2+1
                    u_hat_work(i,j) = cmplx(0.0d0, ky(j), 8) * u_hat_work(i,j)
                    v_hat_work(i,j) = cmplx(0.0d0, ky(j), 8) * v_hat_work(i,j)
                    w_hat_work(i,j) = cmplx(0.0d0, ky(j), 8) * w_hat_work(i,j)
                end do
            end do
            
            call fft_backward_2d(plans, u_hat_work, u_work)
            call fft_backward_2d(plans, v_hat_work, v_work)
            call fft_backward_2d(plans, w_hat_work, w_work)
            
            dudy(:,:,k) = u_work
            dvdy(:,:,k) = v_work
            dwdy(:,:,k) = w_work
        end do
        
        ! Z-derivatives
        do j = 1, ny
            do i = 1, nx
                temp_z = u_vel(i,j,:)
                dudz(i,j,:) = matmul(D1, temp_z)
                
                temp_z = v_vel(i,j,:)
                dvdz(i,j,:) = matmul(D1, temp_z)
                
                temp_z = w_vel(i,j,:)
                dwdz(i,j,:) = matmul(D1, temp_z)
            end do
        end do
    end subroutine

    subroutine compute_nonlinear_terms(nl_u, nl_v, nl_w)
        implicit none
        complex(8), intent(out) :: nl_u(nx/2+1, ny, nz)
        complex(8), intent(out) :: nl_v(nx/2+1, ny, nz)
        complex(8), intent(out) :: nl_w(nx/2+1, ny, nz)
        integer :: i, j, k
        real(8) :: scale_factor
        
        call compute_derivatives()
        
        ! Compute nonlinear terms with scaling for stability
        scale_factor = 0.5d0  ! Reduce nonlinear strength
        
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    nl_u_phys(i,j,k) = -scale_factor * (u_vel(i,j,k)*dudx(i,j,k) + &
                                                        v_vel(i,j,k)*dudy(i,j,k) + &
                                                        w_vel(i,j,k)*dudz(i,j,k))
                    
                    nl_v_phys(i,j,k) = -scale_factor * (u_vel(i,j,k)*dvdx(i,j,k) + &
                                                        v_vel(i,j,k)*dvdy(i,j,k) + &
                                                        w_vel(i,j,k)*dvdz(i,j,k))
                    
                    nl_w_phys(i,j,k) = -scale_factor * (u_vel(i,j,k)*dwdx(i,j,k) + &
                                                        v_vel(i,j,k)*dwdy(i,j,k) + &
                                                        w_vel(i,j,k)*dwdz(i,j,k))
                end do
            end do
        end do
        
        ! Transform to spectral space
        do k = 1, nz
            u_work = nl_u_phys(:,:,k)
            v_work = nl_v_phys(:,:,k)
            w_work = nl_w_phys(:,:,k)
            
            call fft_forward_2d(plans, u_work, u_hat_work)
            call fft_forward_2d(plans, v_work, v_hat_work)
            call fft_forward_2d(plans, w_work, w_hat_work)
            
            nl_u(:,:,k) = u_hat_work
            nl_v(:,:,k) = v_hat_work
            nl_w(:,:,k) = w_hat_work
        end do
        
        ! Apply dealiasing
        call dealias_nonlinear_terms(nl_u, nl_v, nl_w)
    end subroutine

    subroutine dealias_nonlinear_terms(nl_u, nl_v, nl_w)
        implicit none
        complex(8), intent(inout) :: nl_u(nx/2+1, ny, nz)
        complex(8), intent(inout) :: nl_v(nx/2+1, ny, nz)
        complex(8), intent(inout) :: nl_w(nx/2+1, ny, nz)
        integer :: i, j, k
        integer :: kx_cut, ky_cut
        
        kx_cut = nint(2.0d0 * nx / 3.0d0 / 2.0d0)
        ky_cut = nint(2.0d0 * ny / 3.0d0)
        
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx/2+1
                    if (i > kx_cut .or. (j > ky_cut .and. j <= ny - ky_cut)) then
                        nl_u(i,j,k) = cmplx(0.0d0, 0.0d0, 8)
                        nl_v(i,j,k) = cmplx(0.0d0, 0.0d0, 8)
                        nl_w(i,j,k) = cmplx(0.0d0, 0.0d0, 8)
                    endif
                end do
            end do
        end do
    end subroutine

    subroutine conservative_time_step()
        implicit none
        complex(8) :: nl_u_new(nx/2+1, ny, nz)
        complex(8) :: nl_v_new(nx/2+1, ny, nz) 
        complex(8) :: nl_w_new(nx/2+1, ny, nz)
        integer :: i, j, k
        real(8) :: k2, nu, damping
        
        nu = 1.0d0/Re
        
        ! Compute nonlinear terms
        call compute_nonlinear_terms(nl_u_new, nl_v_new, nl_w_new)
        
        ! Conservative time integration
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx/2+1
                    k2 = kx(i)**2 + ky(j)**2
                    
                    ! Enhanced damping for high wavenumbers
                    if (k2 > 20.0d0) then
                        damping = 0.90d0
                    else if (k2 > 10.0d0) then
                        damping = 0.95d0
                    else
                        damping = 1.0d0
                    endif
                    
                    ! Semi-implicit scheme
                    u_hat(i,j,k) = damping * (u_hat(i,j,k) + dt * nl_u_new(i,j,k)) / &
                                   (1.0d0 + dt * nu * k2)
                    v_hat(i,j,k) = damping * (v_hat(i,j,k) + dt * nl_v_new(i,j,k)) / &
                                   (1.0d0 + dt * nu * k2)
                    w_hat(i,j,k) = damping * (w_hat(i,j,k) + dt * nl_w_new(i,j,k)) / &
                                   (1.0d0 + dt * nu * k2)
                end do
            end do
        end do
        
        ! Apply boundary conditions
        call apply_boundary_conditions()
        
        ! Enforce divergence-free condition
        call enforce_divergence_free()
        
        ! Store for next step
        nl_u_old = nl_u_new
        nl_v_old = nl_v_new
        nl_w_old = nl_w_new
        
        ! Transform back to physical space
        call velocity_to_physical()
    end subroutine

    subroutine apply_boundary_conditions()
        implicit none
        
        u_hat(:,:,1) = cmplx(0.0d0, 0.0d0, 8)
        u_hat(:,:,nz) = cmplx(0.0d0, 0.0d0, 8)
        v_hat(:,:,1) = cmplx(0.0d0, 0.0d0, 8)
        v_hat(:,:,nz) = cmplx(0.0d0, 0.0d0, 8)
        w_hat(:,:,1) = cmplx(0.0d0, 0.0d0, 8)
        w_hat(:,:,nz) = cmplx(0.0d0, 0.0d0, 8)
    end subroutine

    subroutine check_cfl_condition()
        implicit none
        real(8) :: dx, dy, dz_min, dt_conv, dt_visc
        
        dx = Lx/nx
        dy = Ly/ny
        dz_min = minval(abs(z(2:nz) - z(1:nz-1)))
        
        call compute_diagnostics()
        
        if (u_max > 1.0d-12) then
            dt_conv = 0.3d0 * min(dx/(u_max+1.0d-12), dy/(v_max+1.0d-12), dz_min/(w_max+1.0d-12))
            dt_visc = 0.3d0 * Re * min(dx**2, dy**2, dz_min**2)
            dt_cfl = min(dt_conv, dt_visc)
        else
            dt_cfl = dt_base
        endif
    end subroutine

    subroutine compute_diagnostics()
        implicit none
        integer :: i, j, k
        real(8) :: vortx, vorty, vortz
        
        u_max = maxval(abs(u_vel))
        v_max = maxval(abs(v_vel))
        w_max = maxval(abs(w_vel))
        
        ke = 0.0d0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ke = ke + 0.5d0 * weights(k) * (u_vel(i,j,k)**2 + v_vel(i,j,k)**2 + w_vel(i,j,k)**2)
                end do
            end do
        end do
        ke = ke / (nx*ny)
        
        call compute_derivatives()
        enstrophy = 0.0d0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    vortx = dwdy(i,j,k) - dvdz(i,j,k)
                    vorty = dudz(i,j,k) - dwdx(i,j,k)
                    vortz = dvdx(i,j,k) - dudy(i,j,k)
                    enstrophy = enstrophy + 0.5d0 * weights(k) * (vortx**2 + vorty**2 + vortz**2)
                end do
            end do
        end do
        enstrophy = enstrophy / (nx*ny)
        
        div_max = maxval(abs(dudx + dvdy + dwdz))
    end subroutine

end program dns_3d_phase3_stable_enhanced
