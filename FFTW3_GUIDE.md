# FFTW3 2D Real-to-Complex Transform Guide

## Overview

This guide demonstrates how to use the native FFTW3 library functions for 2D real-to-complex and complex-to-real transforms in Fortran. FFTW3 is the gold standard for FFT computations and provides highly optimized routines.

## Key Concepts

### Real-to-Complex (R2C) Transform
- **Input**: Real array of size `(nx, ny)`
- **Output**: Complex array of size `(nx/2+1, ny)`
- **Efficiency**: Exploits Hermitian symmetry of real data
- **Storage**: Only stores non-redundant half of the spectrum

### Complex-to-Real (C2R) Transform
- **Input**: Complex array of size `(nx/2+1, ny)`
- **Output**: Real array of size `(nx, ny)`
- **Normalization**: FFTW C2R transforms are **unnormalized** - must divide by `nx*ny`
- **Requirement**: Input must have proper Hermitian symmetry

## Core FFTW3 Functions

### 1. Plan Creation

```fortran
! Forward transform: real(nx,ny) -> complex(nx/2+1,ny)
plan_forward = fftw_plan_dft_r2c_2d(ny, nx, real_data, complex_data, FFTW_ESTIMATE)

! Backward transform: complex(nx/2+1,ny) -> real(nx,ny)
plan_backward = fftw_plan_dft_c2r_2d(ny, nx, complex_data, real_data, FFTW_ESTIMATE)
```

**Important Notes:**
- FFTW uses **row-major** ordering: `fftw_plan_dft_r2c_2d(ny, nx, ...)` for Fortran arrays declared as `array(nx,ny)`
- `FFTW_ESTIMATE`: Fast plan creation (use `FFTW_MEASURE` for production code)
- Plans can be reused for multiple transforms with the same dimensions

### 2. Transform Execution

```fortran
! Execute forward transform
call fftw_execute_dft_r2c(plan_forward, real_data, complex_data)

! Execute backward transform
call fftw_execute_dft_c2r(plan_backward, complex_data, real_data)

! CRITICAL: Normalize after C2R transform
real_data = real_data / real(nx * ny, wp)
```

### 3. Memory Management

```fortran
! Destroy plans when done
call fftw_destroy_plan(plan_forward)
call fftw_destroy_plan(plan_backward)

! Clean up FFTW internal data structures
call fftw_cleanup()
```

## Spectral Derivatives

One of the most powerful applications is computing derivatives in spectral space:

```fortran
! For ∂f/∂x: multiply by i*kx in spectral space
do j = 1, ny
    do i = 1, nxhp
        kx = 2.0_wp * pi * real(i-1, wp) / Lx
        dfdx_hat(i,j) = cmplx(0.0_wp, kx, kind=wp) * f_hat(i,j)
    end do
end do
```

## Array Dimensions and Indexing

### Real Array: `real_data(nx, ny)`
- `i = 1, 2, ..., nx` corresponds to `x = 0, Lx/nx, 2*Lx/nx, ..., (nx-1)*Lx/nx`
- `j = 1, 2, ..., ny` corresponds to `y = 0, Ly/ny, 2*Ly/ny, ..., (ny-1)*Ly/ny`

### Complex Array: `complex_data(nxhp, ny)` where `nxhp = nx/2 + 1`
- `i = 1, 2, ..., nxhp` corresponds to `kx = 0, 2π/Lx, 4π/Lx, ..., π*nx/Lx`
- `j = 1, 2, ..., ny` corresponds to `ky = 0, 2π/Ly, 4π/Ly, ..., 2π*(ny-1)/Ly`

## Compilation Requirements

### Required Libraries
- `libfftw3`: Core FFTW3 library
- `libm`: Math library (usually automatic)

### Fortran Interface
Include the FFTW3 Fortran interface:
```fortran
use, intrinsic :: iso_c_binding
include 'fftw3.f03'
```

### Compiler Command
```bash
gfortran -O3 -o program program.f90 -lfftw3 -lm
```

## Installation

### macOS (Homebrew)
```bash
brew install fftw
```

### Ubuntu/Debian
```bash
sudo apt-get install libfftw3-dev
```

### CentOS/RHEL
```bash
sudo yum install fftw3-devel
```

## Performance Tips

1. **Use FFTW_MEASURE for production**: More expensive plan creation but faster execution
2. **Reuse plans**: Create once, execute many times
3. **Memory alignment**: Use `fftw_malloc()` for optimal performance (advanced)
4. **Thread safety**: FFTW3 can be compiled with OpenMP support

## Common Pitfalls

1. **Dimension ordering**: FFTW uses row-major; Fortran uses column-major
2. **Normalization**: Always divide C2R results by `nx*ny`
3. **Hermitian symmetry**: C2R input must be properly symmetric
4. **Memory management**: Always destroy plans and call `fftw_cleanup()`

## Integration with DNS Solver

For your DNS solver, this approach would replace your current 2D FFT module:

```fortran
! In your setup routine
plan_forward = fftw_plan_dft_r2c_2d(ny, nx, sample_real, sample_complex, FFTW_MEASURE)
plan_backward = fftw_plan_dft_c2r_2d(ny, nx, sample_complex, sample_real, FFTW_MEASURE)

! In your derivative routine
call fftw_execute_dft_r2c(plan_forward, physical_field, spectral_field)
! Apply derivative operators...
call fftw_execute_dft_c2r(plan_backward, spectral_field, derivative_field)
derivative_field = derivative_field / real(nx * ny, wp)
```

This approach eliminates the dimension mismatch issues in your current code by directly using FFTW3's native real-to-complex interface.
