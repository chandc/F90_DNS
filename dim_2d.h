c===============================================================================
c DIM_2D.H - Grid Dimension Parameters for DNS Channel Flow Solver
c
c DESCRIPTION:
c   Defines fundamental grid dimensions and derived parameters for the
c   3D incompressible Navier-Stokes DNS solver. Compatible with both
c   F77 and F90 versions of the code.
c
c GRID PARAMETERS:
c   nx   = 128  : Number of Fourier modes in streamwise direction
c   nz   = 33   : Number of LGL collocation points in wall-normal direction
c
c DERIVED PARAMETERS:
c   nxpp = 130  : Padded grid for FFT efficiency (nx + 2)
c   nxh  = 64   : Half the streamwise modes (nx/2)
c   nxhp = 65   : Half modes plus one (nxh + 1)
c   nxf  = 193  : FFT working array size (3*nx/2 + 1)
c   ntot = 4290 : Total grid points (nxpp * nz)
c   nzm  = 32   : Wall-normal index bound (nz - 1)
c
c MEMORY LAYOUT:
c   Data stored in Fortran column-major order:
c   - First index:  wall-normal direction (z)
c   - Second index: streamwise direction (x)
c   - Physical domain: x ∈ [0, 2π], z ∈ [-1, +1]
c
c AUTHORS:
c   Original: Daniel Chiu-Leung Chan (1993)
c   Updated for pressure BC development (2025)
c
c VERSION: 2025.8.3
c===============================================================================
c-------------------------------------------------------------------
      parameter   ( nx = 128,  nz  = 33)
      parameter   ( nxpp = nx + 2,     nxh=nx/2, nxhp = nxh+1,
     &              nxf = 3*nx/2 + 1,  
     &              ntot = nxpp*nz,  nzm = nz - 1 )
c-------------------------------------------------------------------
