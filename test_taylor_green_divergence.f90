program test_taylor_green_divergence
    use lgl_module
    use fft_module  ! Use the FFT module
    implicit none

    integer, parameter :: Nx = 32, Ny = 16, Nz = 17
    integer, parameter :: Nkx = Nx/2+1
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8), parameter :: Lx = 2.0d0*pi, Ly = 2.0d0*pi, Lz = 2.0d0

    ! Velocity arrays
    real(8) :: u_phys(Nx, Ny, Nz), v_phys(Nx, Ny, Nz), w_phys(Nx, Ny, Nz)
    complex(8) :: u_spec(Nkx, Ny, Nz), v_spec(Nkx, Ny, Nz), w_spec(Nkx, Ny, Nz)
    
    ! Derivative arrays
    real(8) :: dudx(Nx, Ny, Nz), dvdy(Nx, Ny, Nz), dwdz(Nx, Ny, Nz)
    complex(8) :: dudx_spec(Nkx, Ny, Nz), dvdy_spec(Nkx, Ny, Nz)
    
    ! Grid
    real(8) :: x(Nx), y(Ny), z(Nz), zp(Nz), D(Nz, Nz)
    complex(8) :: ikx(Nkx), iky(Ny)
    
    ! FFT setup
    type(fft_plans) :: plans
    real(8) :: sample_real(Nx, Ny)
    complex(8) :: sample_complex(Nkx, Ny)
    
    integer :: i, j, k
    real(8) :: dx, dy, max_div, analytical_div, numerical_div, temp_z(Nz)
    real(8) :: amp, kx_mode, ky_mode, decay_factor
    
    ! Initialize
    call setup_fft_plans(plans, Nx, Ny, sample_real, sample_complex)
    
    ! Setup grid
    dx = Lx / real(Nx, 8)
    dy = Ly / real(Ny, 8)
    
    do i = 1, Nx
        x(i) = real(i-1, 8) * dx
    end do
    
    do j = 1, Ny
        y(j) = real(j-1, 8) * dy
    end do
    
    call lgl_nodes_weights(Nz, z, zp)
    z = Lz/2.0d0 * z
    call differentiation_matrix(Nz, z, D)
    
    ! Setup wavenumbers
    do i = 1, Nkx
        ikx(i) = cmplx(0.0d0, 2.0d0*pi/Lx * real(i-1, 8), 8)
    end do
    
    do j = 1, Ny/2+1
        iky(j) = cmplx(0.0d0, 2.0d0*pi/Ly * real(j-1, 8), 8)
    end do
    
    do j = Ny/2+2, Ny
        iky(j) = cmplx(0.0d0, 2.0d0*pi/Ly * real(j-1-Ny, 8), 8)
    end do
    
    ! Create analytical Taylor-Green vortex (3D divergence-free extension)
    amp = 0.1d0
    kx_mode = 2.0d0*pi/Lx
    ky_mode = 2.0d0*pi/Ly
    
    write(*,*) '=== Taylor-Green Vortex Divergence Test ==='
    write(*,*) 'Mode wavenumbers: kx =', kx_mode, ', ky =', ky_mode
    
    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                ! 3D Taylor-Green vortex (exactly divergence-free)
                ! Based on: u = A*sin(x)*cos(y), v = -A*cos(x)*sin(y), w = 0 (2D)
                ! Extended to 3D with z-dependence that maintains ∇·u = 0
                u_phys(i,j,k) = amp * sin(kx_mode*x(i)) * cos(ky_mode*y(j)) * cos(pi*z(k)/Lz)
                v_phys(i,j,k) = -amp * cos(kx_mode*x(i)) * sin(ky_mode*y(j)) * cos(pi*z(k)/Lz)
                ! For ∇·u = 0: ∂w/∂z = -(∂u/∂x + ∂v/∂y) = 0 for this case
                w_phys(i,j,k) = 0.0d0  ! Keep w=0 for exact divergence-free condition
            end do
        end do
    end do
    
    ! Transform to spectral space
    do k = 1, Nz
        call fft_forward_2d(plans, u_phys(:,:,k), u_spec(:,:,k))
        call fft_forward_2d(plans, v_phys(:,:,k), v_spec(:,:,k))
        call fft_forward_2d(plans, w_phys(:,:,k), w_spec(:,:,k))
    end do
    
    ! Compute spectral derivatives
    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nkx
                dudx_spec(i,j,k) = ikx(i) * u_spec(i,j,k)
                dvdy_spec(i,j,k) = iky(j) * v_spec(i,j,k)
            end do
        end do
    end do
    
    ! Transform derivatives back to physical space
    do k = 1, Nz
        call fft_backward_2d(plans, dudx_spec(:,:,k), dudx(:,:,k))
        call fft_backward_2d(plans, dvdy_spec(:,:,k), dvdy(:,:,k))
    end do
    
    ! Compute z-derivatives using LGL differentiation
    do j = 1, Ny
        do i = 1, Nx
            temp_z = w_phys(i,j,:)
            dwdz(i,j,:) = matmul(D, temp_z)
        end do
    end do
    
    ! Check divergence
    max_div = 0.0d0
    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                numerical_div = dudx(i,j,k) + dvdy(i,j,k) + dwdz(i,j,k)
                max_div = max(max_div, abs(numerical_div))
                
                ! Analytical divergence should be zero for Taylor-Green
                analytical_div = 0.0d0
                
                if (i == 1 .and. j == 1 .and. k == 1) then
                    write(*,'(A,I0,A,I0,A,I0,A)') 'Point (', i, ',', j, ',', k, '):'
                    write(*,'(A,3ES12.4)') '  Velocity:', u_phys(i,j,k), v_phys(i,j,k), w_phys(i,j,k)
                    write(*,'(A,3ES12.4)') '  Derivatives:', dudx(i,j,k), dvdy(i,j,k), dwdz(i,j,k)
                    write(*,'(A,ES12.4)') '  Numerical div:', numerical_div
                    write(*,'(A,ES12.4)') '  Analytical div:', analytical_div
                end if
            end do
        end do
    end do
    
    write(*,'(A,ES12.4)') 'Maximum divergence error:', max_div
    
    if (max_div < 1e-12) then
        write(*,*) 'PASS: Taylor-Green vortex is properly divergence-free!'
    else
        write(*,*) 'FAIL: Large divergence indicates derivative error!'
    end if
    
    call destroy_fft_plans(plans)

end program test_taylor_green_divergence
