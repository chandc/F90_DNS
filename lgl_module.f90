! =============================================================================
! LEGENDRE-GAUSS-LOBATTO (LGL) SPECTRAL ELEMENT MODULE
! =============================================================================
! This module provides:
! - LGL quadrature nodes and weights
! - Legendre polynomials and their derivatives
! - Interpolation and differentiation matrices
! - Mass and stiffness matrices for spectral elements
! =============================================================================

module lgl_module
    use iso_fortran_env, only: wp => real64
    implicit none
    
    private
    public :: lgl_nodes_weights, legendre_poly, legendre_derivative
    public :: interpolation_matrix, differentiation_matrix
    public :: mass_matrix, stiffness_matrix
    public :: lgl_quadrature_1d, lgl_quadrature_3d
    public :: spectral_interpolate, spectral_differentiate
    
    real(wp), parameter :: pi = 3.14159265358979323846_wp
    real(wp), parameter :: eps = 1.0e-15_wp
    
contains

    ! =========================================================================
    ! LGL NODES AND WEIGHTS COMPUTATION
    ! =========================================================================
    subroutine lgl_nodes_weights(n, nodes, weights)
        implicit none
        integer, intent(in) :: n  ! Number of LGL points (polynomial order + 1)
        real(wp), intent(out) :: nodes(n), weights(n)
        
        integer :: i, j, max_iter
        real(wp) :: x, x_old, p, dp, d2p
        real(wp) :: tolerance
        
        tolerance = 1.0e-15_wp
        max_iter = 100
        
        if (n == 1) then
            nodes(1) = 0.0_wp
            weights(1) = 2.0_wp
            return
        elseif (n == 2) then
            nodes(1) = -1.0_wp
            nodes(2) = 1.0_wp
            weights(1) = 1.0_wp
            weights(2) = 1.0_wp
            return
        endif
        
        ! End points
        nodes(1) = -1.0_wp
        nodes(n) = 1.0_wp
        weights(1) = 2.0_wp / (real(n*(n-1), wp))
        weights(n) = weights(1)
        
        ! Interior points using Newton-Raphson iteration
        do i = 2, n-1
            ! Initial guess using Chebyshev nodes
            x = -cos(pi * real(i-1, wp) / real(n-1, wp))
            
            do j = 1, max_iter
                x_old = x
                call legendre_poly_and_deriv(n-1, x, p, dp, d2p)
                
                ! Newton-Raphson update for LGL condition: P'_n(x) = 0
                x = x_old - dp / d2p
                
                if (abs(x - x_old) < tolerance) exit
            end do
            
            nodes(i) = x
            call legendre_poly(n-1, x, p)
            weights(i) = 2.0_wp / (real(n*(n-1), wp) * p**2)
        end do
        
    end subroutine lgl_nodes_weights
    
    ! =========================================================================
    ! LEGENDRE POLYNOMIAL EVALUATION
    ! =========================================================================
    subroutine legendre_poly(n, x, p)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: x
        real(wp), intent(out) :: p
        
        integer :: k
        real(wp) :: p0, p1, p2
        
        if (n == 0) then
            p = 1.0_wp
        elseif (n == 1) then
            p = x
        else
            p0 = 1.0_wp
            p1 = x
            do k = 2, n
                p2 = (real(2*k-1, wp) * x * p1 - real(k-1, wp) * p0) / real(k, wp)
                p0 = p1
                p1 = p2
            end do
            p = p1
        endif
        
    end subroutine legendre_poly
    
    ! =========================================================================
    ! LEGENDRE POLYNOMIAL AND ITS DERIVATIVES
    ! =========================================================================
    subroutine legendre_poly_and_deriv(n, x, p, dp, d2p)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: x
        real(wp), intent(out) :: p, dp, d2p
        
        integer :: k
        real(wp) :: p0, p1, p2, dp0, dp1, dp2, d2p0, d2p1, d2p2
        
        if (n == 0) then
            p = 1.0_wp
            dp = 0.0_wp
            d2p = 0.0_wp
        elseif (n == 1) then
            p = x
            dp = 1.0_wp
            d2p = 0.0_wp
        else
            ! Initialize
            p0 = 1.0_wp; dp0 = 0.0_wp; d2p0 = 0.0_wp
            p1 = x; dp1 = 1.0_wp; d2p1 = 0.0_wp
            
            do k = 2, n
                ! Recurrence for polynomial
                p2 = (real(2*k-1, wp) * x * p1 - real(k-1, wp) * p0) / real(k, wp)
                
                ! Recurrence for first derivative
                dp2 = (real(2*k-1, wp) * (x * dp1 + p1) - real(k-1, wp) * dp0) / real(k, wp)
                
                ! Recurrence for second derivative
                d2p2 = (real(2*k-1, wp) * (x * d2p1 + 2.0_wp * dp1) - real(k-1, wp) * d2p0) / real(k, wp)
                
                p0 = p1; p1 = p2
                dp0 = dp1; dp1 = dp2
                d2p0 = d2p1; d2p1 = d2p2
            end do
            
            p = p1
            dp = dp1
            d2p = d2p1
        endif
        
    end subroutine legendre_poly_and_deriv
    
    ! =========================================================================
    ! LEGENDRE DERIVATIVE AT POINT
    ! =========================================================================
    subroutine legendre_derivative(n, x, dp)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: x
        real(wp), intent(out) :: dp
        
        real(wp) :: p, d2p
        call legendre_poly_and_deriv(n, x, p, dp, d2p)
        
    end subroutine legendre_derivative
    
    ! =========================================================================
    ! INTERPOLATION MATRIX
    ! =========================================================================
    subroutine interpolation_matrix(n_old, nodes_old, n_new, nodes_new, interp_mat)
        implicit none
        integer, intent(in) :: n_old, n_new
        real(wp), intent(in) :: nodes_old(n_old), nodes_new(n_new)
        real(wp), intent(out) :: interp_mat(n_new, n_old)
        
        integer :: i, j, k
        real(wp) :: num, den, x_target
        
        do i = 1, n_new
            x_target = nodes_new(i)
            do j = 1, n_old
                ! Lagrange interpolation basis function
                num = 1.0_wp
                den = 1.0_wp
                
                do k = 1, n_old
                    if (k /= j) then
                        num = num * (x_target - nodes_old(k))
                        den = den * (nodes_old(j) - nodes_old(k))
                    endif
                end do
                
                if (abs(den) > eps) then
                    interp_mat(i, j) = num / den
                else
                    interp_mat(i, j) = 0.0_wp
                endif
            end do
        end do
        
    end subroutine interpolation_matrix
    
    ! =========================================================================
    ! DIFFERENTIATION MATRIX FOR LGL NODES (Standard Lagrange Formula)
    ! =========================================================================
    subroutine differentiation_matrix(n, nodes, diff_mat)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: nodes(n)
        real(wp), intent(out) :: diff_mat(n, n)
        
        integer :: i, j, k
        real(wp) :: prod1, prod2
        
        diff_mat = 0.0_wp
        
        ! Special case for n = 1
        if (n == 1) then
            diff_mat(1, 1) = 0.0_wp
            return
        endif
        
        ! Standard differentiation matrix for arbitrary nodes using Lagrange basis
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    ! Diagonal elements: sum of 1/(xi - xk) for k != i
                    diff_mat(i, i) = 0.0_wp
                    do k = 1, n
                        if (k /= i) then
                            diff_mat(i, i) = diff_mat(i, i) + 1.0_wp / (nodes(i) - nodes(k))
                        endif
                    end do
                else
                    ! Off-diagonal elements
                    prod1 = 1.0_wp
                    prod2 = 1.0_wp
                    
                    ! Compute products for Lagrange basis derivative
                    do k = 1, n
                        if (k /= i .and. k /= j) then
                            prod1 = prod1 * (nodes(i) - nodes(k))
                        endif
                        if (k /= j) then
                            prod2 = prod2 * (nodes(j) - nodes(k))
                        endif
                    end do
                    
                    if (abs(prod2) > eps) then
                        diff_mat(i, j) = prod1 / prod2
                    else
                        diff_mat(i, j) = 0.0_wp
                    endif
                endif
            end do
        end do
        
    end subroutine differentiation_matrix
    
    ! =========================================================================
    ! MASS MATRIX (for L2 inner products)
    ! =========================================================================
    subroutine mass_matrix(n, nodes, weights, mass_mat)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: nodes(n), weights(n)
        real(wp), intent(out) :: mass_mat(n, n)
        
        integer :: i, j
        
        mass_mat = 0.0_wp
        
        ! For LGL quadrature, mass matrix is diagonal
        do i = 1, n
            mass_mat(i, i) = weights(i)
        end do
        
    end subroutine mass_matrix
    
    ! =========================================================================
    ! STIFFNESS MATRIX (for derivatives)
    ! =========================================================================
    subroutine stiffness_matrix(n, nodes, weights, diff_mat, stiff_mat)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: nodes(n), weights(n), diff_mat(n, n)
        real(wp), intent(out) :: stiff_mat(n, n)
        
        integer :: i, j, k
        
        ! Stiffness matrix: S_ij = ∫ dψ_i/dx * dψ_j/dx dx
        ! Using LGL quadrature: S_ij = Σ_k w_k * D_ki * D_kj
        
        do i = 1, n
            do j = 1, n
                stiff_mat(i, j) = 0.0_wp
                do k = 1, n
                    stiff_mat(i, j) = stiff_mat(i, j) + weights(k) * diff_mat(k, i) * diff_mat(k, j)
                end do
            end do
        end do
        
    end subroutine stiffness_matrix
    
    ! =========================================================================
    ! 1D LGL QUADRATURE
    ! =========================================================================
    function lgl_quadrature_1d(n, nodes, weights, func_values) result(integral)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: nodes(n), weights(n), func_values(n)
        real(wp) :: integral
        
        integer :: i
        
        integral = 0.0_wp
        do i = 1, n
            integral = integral + weights(i) * func_values(i)
        end do
        
    end function lgl_quadrature_1d
    
    ! =========================================================================
    ! 3D LGL QUADRATURE
    ! =========================================================================
    function lgl_quadrature_3d(nx, ny, nz, wx, wy, wz, func_values) result(integral)
        implicit none
        integer, intent(in) :: nx, ny, nz
        real(wp), intent(in) :: wx(nx), wy(ny), wz(nz)
        real(wp), intent(in) :: func_values(nx, ny, nz)
        real(wp) :: integral
        
        integer :: i, j, k
        
        integral = 0.0_wp
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    integral = integral + wx(i) * wy(j) * wz(k) * func_values(i, j, k)
                end do
            end do
        end do
        
    end function lgl_quadrature_3d
    
    ! =========================================================================
    ! SPECTRAL INTERPOLATION
    ! =========================================================================
    subroutine spectral_interpolate(n_old, values_old, n_new, nodes_old, nodes_new, values_new)
        implicit none
        integer, intent(in) :: n_old, n_new
        real(wp), intent(in) :: values_old(n_old), nodes_old(n_old), nodes_new(n_new)
        real(wp), intent(out) :: values_new(n_new)
        
        real(wp) :: interp_mat(n_new, n_old)
        integer :: i, j
        
        call interpolation_matrix(n_old, nodes_old, n_new, nodes_new, interp_mat)
        
        do i = 1, n_new
            values_new(i) = 0.0_wp
            do j = 1, n_old
                values_new(i) = values_new(i) + interp_mat(i, j) * values_old(j)
            end do
        end do
        
    end subroutine spectral_interpolate
    
    ! =========================================================================
    ! SPECTRAL DIFFERENTIATION
    ! =========================================================================
    subroutine spectral_differentiate(n, nodes, values, derivatives)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: nodes(n), values(n)
        real(wp), intent(out) :: derivatives(n)
        
        real(wp) :: diff_mat(n, n)
        integer :: i, j
        
        call differentiation_matrix(n, nodes, diff_mat)
        
        do i = 1, n
            derivatives(i) = 0.0_wp
            do j = 1, n
                derivatives(i) = derivatives(i) + diff_mat(i, j) * values(j)
            end do
        end do
        
    end subroutine spectral_differentiate
    
end module lgl_module
