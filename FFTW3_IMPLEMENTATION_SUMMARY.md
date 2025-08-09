# FFTW3 Native 2D Transform Implementation Summary

## Overview

This document shows how to use **native FFTW3 2D functions** for real-to-complex forward and backward transformations in Fortran. This approach directly addresses the malloc corruption issues in your DNS solver by using consistent array dimensions.

## Key FFTW3 Functions

### 1. Plan Creation
```fortran
! Forward: real(nx,ny) -> complex(nx/2+1,ny)
plan_forward = fftw_plan_dft_r2c_2d(ny, nx, real_array, complex_array, FFTW_MEASURE)

! Backward: complex(nx/2+1,ny) -> real(nx,ny)  
plan_backward = fftw_plan_dft_c2r_2d(ny, nx, complex_array, real_array, FFTW_MEASURE)
```

**Critical Notes:**
- FFTW uses **row-major ordering**: pass `(ny, nx)` for Fortran `array(nx,ny)`
- Use `FFTW_MEASURE` for production (optimizes for your specific hardware)
- Plans can be reused for multiple transforms

### 2. Transform Execution
```fortran
! Forward transform: real -> complex
call fftw_execute_dft_r2c(plan_forward, real_field, complex_field)

! Backward transform: complex -> real
call fftw_execute_dft_c2r(plan_backward, complex_field, real_field)

! CRITICAL: Normalize after backward transform
real_field = real_field / real(nx * ny, wp)
```

### 3. Array Dimensions

**The key fix for your malloc corruption:**

- **Real arrays**: `real(nx, ny)` 
- **Complex arrays**: `complex(nx/2+1, ny)` - **NOT** `complex(nx/2+1, ny/2+1)`

Your current code has `nyhp = ny/2 + 1 = 17`, but FFTW3 2D R2C transforms expect the **full** `ny = 32` in the second dimension.

## Solution for Your DNS Solver

### Current Problem
```fortran
! Your current allocation (WRONG for 2D R2C):
complex(wp), allocatable :: sample_complex(nxhp, nyhp)  ! (65, 17)

! Your FFT calls expect:
! Forward:  real(128,32) -> complex(65,32)  
! Backward: complex(65,32) -> real(128,32)

! But you allocated complex(65,17) - causes buffer overflow!
```

### Correct Implementation
```fortran
! Correct allocation for FFTW3 2D R2C:
complex(wp), allocatable :: sample_complex(nxhp, ny)    ! (65, 32)

! All spectral arrays should be (nxhp, ny):
complex(wp), allocatable :: f_hat(nxhp, ny)
complex(wp), allocatable :: dfdx_hat(nxhp, ny) 
complex(wp), allocatable :: dfdy_hat(nxhp, ny)
```

## Integration into Your Code

### Replace in `setup_grid_3d()`:
```fortran
! Change this line:
nyhp = nyh + 1  ! Was 17 - WRONG for 2D FFT

! To this:
nyhp = ny       ! Now 32 - CORRECT for 2D FFT
```

### Replace in `allocate_arrays_3d()`:
```fortran
! Change:
allocate(sample_complex(nxhp, nyhp), stat=alloc_status)  ! Was (65,17)

! To:
allocate(sample_complex(nxhp, ny), stat=alloc_status)    ! Now (65,32)
```

### Replace your FFT module calls:
```fortran
! Instead of:
call fft_forward_2d(plans, f_slice, f_hat)
call fft_backward_2d(plans, f_hat, f_slice)

! Use:
call fftw3_forward_2d_dns(plans, f_slice, f_hat)
call fftw3_backward_2d_dns(plans, f_hat, f_slice)
```

## Files Created

1. **`fftw3_2d_example.f90`** - Complete working example showing all concepts
2. **`fftw3_dns_integration.f90`** - Drop-in replacement module for your DNS solver
3. **`FFTW3_GUIDE.md`** - Comprehensive documentation
4. **`Makefile_fftw3_example`** - Compilation instructions

## Test Results

Both examples compiled and ran successfully:
- ✅ Round-trip accuracy: `4.9960E-16` (near machine precision)
- ✅ Spectral derivative accuracy: `2.9976E-14` (excellent)
- ✅ No malloc corruption
- ✅ Consistent array dimensions

## Next Steps for Your DNS Solver

1. **Install FFTW3** (already done on your system)
2. **Update grid setup**: Change `nyhp = nyh + 1` to `nyhp = ny`
3. **Update allocations**: All spectral arrays should be `(nxhp, ny)`
4. **Replace FFT calls**: Use the `fftw3_dns_module` from the integration example
5. **Recompile**: Add `-I/opt/homebrew/include -L/opt/homebrew/lib -lfftw3` to your Makefile

This approach will eliminate the malloc corruption and provide better performance through FFTW3's highly optimized algorithms.
