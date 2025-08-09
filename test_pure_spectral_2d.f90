program test_pure_spectral_2d
    use fft_module
    implicit none

    ! Pure 2D test - ignore z direction completely
    integer, parameter :: Nx = 16, Ny = 8
    integer, parameter :: Nkx = Nx/2+1
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8), parameter :: Lx = 2.0d0*pi, Ly = 2.0d0*pi
    real(8), parameter :: dt = 1.0d-5, Re = 500.0d0
    integer, parameter :: nsteps = 200

    real(8) :: u_phys(Nx, Ny), v_phys(Nx, Ny)
    complex(8) :: u_spec(Nkx, Ny), v_spec(Nkx, Ny)
    
    real(8) :: x(Nx), y(Ny)
    complex(8) :: ikx(Nkx), iky(Ny)
    type(fft_plans) :: plans
    real(8) :: sample_real(Nx, Ny)
    complex(8) :: sample_complex(Nkx, Ny)
    
    integer :: i, j, nt
    real(8) :: dx, dy, energy, energy_initial, max_div
    
    call setup_fft_plans(plans, Nx, Ny, sample_real, sample_complex)
    
    dx = Lx / real(Nx, 8)
    dy = Ly / real(Ny, 8)
    
    do i = 1, Nx
        x(i) = real(i-1, 8) * dx
    end do
    
    do j = 1, Ny
        y(j) = real(j-1, 8) * dy
    end do
    
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
    
    ! Create exactly divergence-free field using stream function
    ! ψ = 0.1*cos(x)*cos(y) => u = ∂ψ/∂y = -0.1*cos(x)*sin(y), v = -∂ψ/∂x = 0.1*sin(x)*cos(y)
    do j = 1, Ny
        do i = 1, Nx
            u_phys(i,j) = -0.1d0 * cos(x(i)) * sin(y(j))
            v_phys(i,j) = 0.1d0 * sin(x(i)) * cos(y(j))
        end do
    end do
    
    ! Check initial divergence
    call check_divergence_2d()
    write(*,*) 'Initial divergence:', max_div
    
    ! Transform to spectral
    call fft_forward_2d(plans, u_phys, u_spec)
    call fft_forward_2d(plans, v_phys, v_spec)
    
    energy_initial = compute_energy_2d()
    
    write(*,*) '=== 2D DNS without pressure ==='
    write(*,'(A)') 'Step   Energy Ratio    Max|Div|'
    
    do nt = 1, nsteps
        
        ! Only viscous diffusion
        do j = 1, Ny
            do i = 1, Nkx
                u_spec(i,j) = u_spec(i,j) * (1.0d0 - dt * (1.0d0/Re) * (real(ikx(i))**2 + real(iky(j))**2))
                v_spec(i,j) = v_spec(i,j) * (1.0d0 - dt * (1.0d0/Re) * (real(ikx(i))**2 + real(iky(j))**2))
            end do
        end do
        
        if (mod(nt, 50) == 0) then
            call fft_backward_2d(plans, u_spec, u_phys)
            call fft_backward_2d(plans, v_spec, v_phys)
            energy = compute_energy_2d()
            call check_divergence_2d()
            write(*,'(I4,2ES14.6)') nt, energy/energy_initial, max_div
        end if
        
    end do
    
    call destroy_fft_plans(plans)

contains

    function compute_energy_2d() result(energy)
        real(8) :: energy
        energy = 0.0d0
        do j = 1, Ny
            do i = 1, Nx
                energy = energy + u_phys(i,j)**2 + v_phys(i,j)**2
            end do
        end do
        energy = 0.5d0 * energy * dx * dy
    end function compute_energy_2d
    
    subroutine check_divergence_2d()
        real(8) :: dudx_phys(Nx, Ny), dvdy_phys(Nx, Ny), div_phys(Nx, Ny)
        complex(8) :: dudx_spec(Nkx, Ny), dvdy_spec(Nkx, Ny)
        
        ! Compute spectral derivatives
        do j = 1, Ny
            do i = 1, Nkx
                dudx_spec(i,j) = ikx(i) * u_spec(i,j)
                dvdy_spec(i,j) = iky(j) * v_spec(i,j)
            end do
        end do
        
        ! Transform to physical
        call fft_backward_2d(plans, dudx_spec, dudx_phys)
        call fft_backward_2d(plans, dvdy_spec, dvdy_phys)
        
        div_phys = dudx_phys + dvdy_phys
        max_div = maxval(abs(div_phys))
        
    end subroutine check_divergence_2d

end program test_pure_spectral_2d
