# DNS Channel Flow Solver - Advanced Flow Control

## Overview

This repository contains modern Fortran 90 implementations of incompressible Navier-Stokes DNS (Direct Numerical Simulation) solvers for channel flow. The code provides **both 2D and 3D versions** and is a complete modernization of the original F77 solver developed by Daniel Chiu-Leung Chan (1993), featuring **advanced flow control systems**, enhanced pressure boundary condition treatment, and comprehensive analysis tools.

## Available Solvers

### 2D Channel Flow Solver (`DNS_pressure_BC_2D.f90`)
- **2D Navier-Stokes equations**: Streamwise (u) and wall-normal (w) velocity components
- **Grid**: nx×nz with spectral methods in x-direction, LGL collocation in z-direction
- **No spanwise variation**: ∂/∂y = 0 (simplified 2D formulation)
- **Ideal for**: Parameter studies, method validation, rapid prototyping
- **Perfect divergence-free**: Machine precision incompressibility (div ~10⁻¹⁶)

### 3D Channel Flow Solver (`DNS_pressure_BC_3D.f90`)
- **Full 3D Navier-Stokes equations**: All three velocity components (u,v,w)
- **Grid**: nx×ny×nz with spectral methods in x,y-directions, LGL collocation in z-direction
- **Complete 3D physics**: Captures all turbulent structures and instabilities
- **Ideal for**: Full turbulence simulations, transition studies, 3D flow phenomena

## Key Features

### Advanced Flow Control System (Both 2D and 3D)

The solvers feature a sophisticated dual-mode flow control system for precise channel flow management:

#### Method 1: Constant Pressure Gradient
- **Fixed driving force**: ∂p/∂x = constant (user-specified)
- **Natural flow development**: Bulk velocity varies with flow evolution
- **Ideal for**: Fundamental studies, transition analysis, parameter sweeps
- **Simple and robust**: No feedback control complexity

#### Method 2: Constant Bulk Velocity (PI Control)
- **Target flow rate**: U_bulk = constant (user-specified) 
- **Automated control**: Pressure gradient adjusts dynamically to maintain flow rate
- **PI controller algorithm**: dp/dt = Kp×error + Ki×∫error dt
- **Anti-windup protection**: Prevents integrator saturation during transients
- **Safety bounds**: Pressure gradient automatically limited to [0.1, 10.0] range
- **Real-time feedback**: Live monitoring of control performance and statistics
- **Spectral accuracy**: LGL quadrature ensures precise bulk velocity calculation
- **Ideal for**: Industrial applications, flow rate studies, controlled experiments

#### Control System Features
- **Seamless switching**: Change control method via input parameter
- **Robust stability**: Tested across wide range of Reynolds numbers
- **Live diagnostics**: Real-time control error and performance metrics
- **Parameter flexibility**: Adjustable PI gains and update frequency

### Numerical Methods
- **Spectral Methods**: Fourier decomposition in streamwise (and spanwise for 3D) directions
- **LGL Collocation**: Legendre-Gauss-Lobatto points in wall-normal direction
- **Time Integration**: Crank-Nicolson for viscous terms, 4th-order Runge-Kutta for convection
- **Pressure Solver**: CGS iterative method with F77 compatibility
- **Boundary Conditions**: Kim & Moin approach for viscous wall treatment
- **Fractional Step Method**: Proper dt scaling in all pressure correction steps

### Pressure Equation Innovation
- **Bottom-wall-only pressure pinning** for zero mode (kₓ=0) stability
- **Unified matrix construction** for all Fourier modes (F77 approach)
- **CGS iterative solver** robust to near-singular systems
- **No artificial Dirichlet conditions** - uses natural boundary conditions

### Grid Configuration
#### 2D Solver:
- **nx = 128**: Fourier modes in streamwise direction
- **nz = 33**: LGL collocation points in wall-normal direction  
- **Domain**: x ∈ [0, 2π], z ∈ [-1, +1] (channel half-height = 1)

#### 3D Solver:
- **nx = 64**: Fourier modes in streamwise direction
- **ny = 32**: Fourier modes in spanwise direction
- **nz = 33**: LGL collocation points in wall-normal direction  
- **Domain**: x ∈ [0, 2π], y ∈ [0, 2π], z ∈ [-1, +1]

### Reynolds Number
- **Re = 180** (2D) / **Re = 100** (3D): Based on channel half-height

## Features

### Advanced Flow Control System
- **Two Control Methods**: Constant pressure gradient or constant bulk velocity
- **PI Controller**: Automated pressure gradient adjustment for flow rate control
- **Spectral Integration**: LGL quadrature for accurate bulk velocity calculation
- **Real-time Feedback**: Live monitoring of flow control performance
- **Safety Features**: Anti-windup protection and gradient bounds

### Numerical Methods
- **High-order accuracy**: Spectral methods in all directions
- **Efficient FFT**: FFTW3 library with OpenMP parallelization
- **Stable time stepping**: Crank-Nicolson viscous step, explicit convection
- **Projection method**: Pressure correction for incompressibility
- **Advanced boundary conditions**: Kim & Moin extrapolated velocity BCs

### Analysis and Visualization Tools
- **Jupyter Notebooks**: Comprehensive turbulent flow analysis
- **Velocity Profile Analysis**: Detailed spatial flow field examination
- **Reynolds-Bulk Velocity Correlations**: Laminar and turbulent flow relationships
- **Real-time Monitoring**: Live simulation progress tracking

### Performance Features
- **OpenMP parallelization**: Multi-threaded FFTW and linear algebra
- **Optimized compilation**: -O3 optimization with double precision
- **Memory efficient**: Spectral transformations minimize storage
- **Background execution**: Automated monitoring and logging

## Requirements

### Dependencies
- **Fortran compiler**: gfortran (GCC) or Intel Fortran
- **FFTW3**: Fast Fourier Transform library with OpenMP support
- **LAPACK/BLAS**: Linear algebra libraries
- **Make**: Build system

### macOS Installation (Homebrew)
```bash
brew install gcc fftw lapack
```

### Ubuntu/Debian Installation
```bash
sudo apt-get install gfortran libfftw3-dev liblapack-dev libblas-dev
```

## Quick Start

### 1. Clone and Build
```bash
git clone https://github.com/chandc/F90_DNS.git
cd F90_DNS

# Build 2D solver
make -f Makefile_2D_pressure_BC

# Build 3D solver  
make -f Makefile_3D_pressure_BC

# Or build both (legacy)
make all
```

### 2. Run Simulations
```bash
# Run 2D solver
./dns_pressure_bc

# Run 3D solver
./dns_3d_pressure_bc

# Interactive run (legacy)
make run

# Background run with monitoring
make run-bg
make monitor

# Run with specific thread count
make run-omp4
```

### 3. Check Progress
```bash
make status
```

## Input Parameters

Configure your simulation by editing the appropriate input file:

### Basic Configuration (`input.dat`)
```fortran
&grid
  nx_input = 128,    ! Streamwise grid points
  nz_input = 33      ! Wall-normal grid points
/

&time_control
  istart = 0,        ! 0=new run, 1=restart
  dt = 0.01,         ! Time step
  nsteps = 50000,    ! Number of steps
  nwrt = 2500        ! Output frequency
/

&simulation
  alpha = 1.0,       ! Streamwise wavenumber
  re = 180.0,        ! Reynolds number
  use_crank_nicolson = .true.
/

&flow_control
  flow_control_method = 1,           ! 1=const pressure, 2=const volume
  target_pressure_gradient = 1.0,   ! Target dP/dx (method 1)
  target_bulk_velocity = 1.0,       ! Target U_bulk (method 2)
  controller_gain = 0.1,             ! PI controller gain (method 2)
  controller_update_freq = 10        ! Update frequency (method 2)
/
```

### Pre-configured Examples
- **`input_constant_pressure.dat`**: Constant pressure gradient flow control
- **`input_constant_volume.dat`**: Constant bulk velocity with PI controller  
- **`input_nondimensional.dat`**: Non-dimensional parameter setup
- **`input_high_resolution.dat`**: Higher resolution simulations
- **Additional examples**: Various specialized configurations for different studies

## File Structure

```
F90_DNS/
├── DNS_pressure_BC.f90        # Main DNS solver with flow control
├── lgl_module.f90             # LGL grid and differentiation
├── fft_module.f90             # FFTW interface
├── Makefile                   # Advanced build system with monitoring
├── FLOW_CONTROL_README.md     # Detailed flow control documentation
├── input*.dat                 # Various simulation configurations
├── start.dat                  # Initial/restart conditions  
├── run.dat                    # Current simulation state
├── simulation.log             # Runtime output log
├── turbulent_channel_flow_correlations.ipynb  # Flow analysis notebook
└── velocity_profile_analysis.ipynb            # Velocity field analysis
```

### Flow Control Input Files
- `input_constant_pressure.dat` - Method 1: Fixed pressure gradient
- `input_constant_volume.dat` - Method 2: PI-controlled flow rate
- `input_nondimensional.dat` - Non-dimensional parameters
- `input_high_resolution.dat` - High-resolution studies
- Additional specialized configurations for various flow conditions

## Makefile Targets

### Basic Operations
- `make all` - Compile the solver (default)
- `make clean` - Remove generated files
- `make help` - Show all available targets

### Running Simulations
- `make run` - Interactive execution
- `make run-bg` - Background execution with logging
- `make monitor` - Monitor background simulation (30s intervals)
- `make status` - Check simulation status
- `make stop` - Stop background simulation

### Flow Control Examples
- `make run` with `input_constant_pressure.dat` - Fixed pressure gradient
- `make run` with `input_constant_volume.dat` - PI-controlled flow rate
- `make run` with `input_nondimensional.dat` - Non-dimensional setup

### Parallel Execution
- `make run-omp2` - Run with 2 OpenMP threads
- `make run-omp4` - Run with 4 OpenMP threads  
- `make run-omp8` - Run with 8 OpenMP threads

### Development
- `make debug` - Debug build with bounds checking
- `make test` - Quick compilation test
- `make check-deps` - Verify dependencies

## Flow Control System

### Method 1: Constant Pressure Gradient
- **Use case**: Fixed driving force applications
- **Implementation**: Applies constant ∂p/∂x at each time step
- **Configuration**: Set `flow_control_method = 1` and `target_pressure_gradient`
- **Advantages**: Simple, predictable, computationally efficient

### Method 2: Constant Bulk Velocity (PI Controller)
- **Use case**: Flow rate-controlled applications
- **Implementation**: PI controller adjusts pressure gradient to maintain target bulk velocity
- **Configuration**: Set `flow_control_method = 2`, `target_bulk_velocity`, `controller_gain`, `controller_update_freq`
- **Features**:
  - **Proportional control**: Immediate response to velocity error
  - **Integral control**: Eliminates steady-state error
  - **Anti-windup protection**: Prevents controller saturation
  - **Safety bounds**: Pressure gradient limited to [0.1, 10.0]

### Bulk Velocity Calculation
- **LGL Integration**: Uses Legendre-Gauss-Lobatto quadrature for spectral accuracy
- **Formula**: `U_bulk = (1/H) ∫[0 to H] u(y) dy`
- **Real-time monitoring**: Live feedback during simulation

### Controller Parameters
- **`controller_gain`**: Start with 0.05-0.1, adjust based on stability
- **`controller_update_freq`**: Typical range 5-20 time steps
- **`target_bulk_velocity`**: Should be physically reasonable for given Reynolds number

For detailed flow control documentation, see **`FLOW_CONTROL_README.md`**.

## Algorithm Details

### Fractional Step Method
1. **Viscous Step**: `(I - 0.5*dt/Re ∇²)u* = (I + 0.5*dt/Re ∇²)u^n + dt·RHS`
2. **Pressure Step**: `∇²φ = ∇·u*/dt`
3. **Velocity Correction**: `u^{n+1} = u* - dt∇φ`

### Boundary Conditions

- **Velocity**: No-slip at walls (`u = w = 0`)
- **Pressure**: Dirichlet (`φ = 0`) at walls for numerical stability
- **Periodicity**: Streamwise and spanwise directions

### Kim & Moin Intermediate Velocity Boundary Conditions

The solver implements the **Kim & Moin (1985) approach** for handling boundary conditions in the fractional step method, which avoids numerical inconsistencies that arise when applying no-slip conditions directly to the intermediate velocity field.

#### Traditional Problem

In standard fractional step methods, applying no-slip boundary conditions (`u* = 0`) to the intermediate velocity creates an inconsistency:

- The intermediate velocity `u*` from the viscous step doesn't satisfy the divergence-free condition
- Direct application of `u* = 0` at walls can lead to numerical instabilities
- The pressure gradient `∇φ` may not properly satisfy the boundary conditions

#### Kim & Moin Solution

The **Kim & Moin method** resolves this by:

1. **Intermediate Velocity Step**: Solve the viscous step with **extrapolated boundary conditions**:

   ```fortran
   u*|wall = extrapolated from interior points (not zero)
   ```

2. **Pressure Boundary Conditions**: Apply appropriate pressure boundary conditions that ensure the final velocity satisfies no-slip:

   ```fortran
   ∂φ/∂n|wall = (1/dt) * u*|wall
   ```

3. **Final Velocity Correction**: The pressure correction automatically enforces no-slip:

   ```fortran
   u^{n+1}|wall = u*|wall - dt(∂φ/∂n)|wall = 0
   ```

#### Implementation Details

- **Wall-normal derivative**: Uses LGL differentiation matrix for accurate `∂φ/∂n` calculation
- **Extrapolation**: Quadratic extrapolation from interior LGL points to wall boundaries
- **Pressure pinning**: Bottom wall pressure pinned for kₓ=0 mode stability
- **CGS solver**: Robust iterative solution of the pressure Poisson equation

#### Advantages

- **Numerical stability**: Eliminates artificial pressure boundary layer
- **Accuracy**: Maintains spectral accuracy near walls
- **Physical consistency**: Proper treatment of viscous wall effects
- **Robust convergence**: Works reliably for wide range of Reynolds numbers

This implementation in the DNS solver shows excellent agreement with theory:

- **Wall shear stress error**: < 1% compared to analytical Poiseuille flow
- **Divergence-free**: Machine precision enforcement (`div_max ≈ 0`)
- **Energy conservation**: Maintains numerical stability over long integrations

#### Reference

Kim, J., & Moin, P. (1985). *Application of a fractional-step method to incompressible Navier-Stokes equations*. Journal of Computational Physics, 59(2), 308-323.

### Spectral Methods

- **LGL Collocation**: Legendre-Gauss-Lobatto points for wall-normal direction
- **Fourier Transform**: FFTW3 for periodic directions
- **Matrix-free**: Efficient spectral differentiation operators

## Analysis Tools

### Jupyter Notebooks

#### 1. Turbulent Channel Flow Correlations (`turbulent_channel_flow_correlations.ipynb`)
- **Bulk velocity vs Reynolds number relationships** for laminar and turbulent regimes
- **Flow control system design guidelines**
- **Target velocity estimation** for unknown flow profiles  
- **Transition analysis** from laminar to turbulent flow
- **Engineering correlations** for practical applications

#### 2. Velocity Profile Analysis (`velocity_profile_analysis.ipynb`)
- **Spatial flow field examination** at different x-locations
- **Spectral method validation** with theoretical profiles
- **LGL grid performance analysis**
- **Boundary condition verification**
- **Real-time data processing** from simulation output

### Features
- **Publication-quality plots** with matplotlib
- **Interactive parameter exploration**
- **Theoretical comparison** with analytical solutions
- **Data export capabilities** for further analysis
- **Integration with simulation output** formats

## Performance Characteristics

### Typical Performance (MacBook Pro M3 Max)
- **Grid**: 128×33×64 (streamwise×wall-normal×spanwise)
- **Speed**: ~2-3 time units per wall-clock second
- **Memory**: ~50-100 MB RAM usage
- **Scaling**: Near-linear with OpenMP threads

### Accuracy
- **Spectral convergence**: Exponential error reduction with grid refinement
- **Energy conservation**: Machine precision conservation in inviscid limit
- **Divergence-free**: Incompressibility enforced to machine precision

## Original Development

Originally developed by **Daniel Chiu-Leung Chan** in 1993, this modernized F90 version maintains the core algorithmic innovations while adding:
- **Advanced flow control system** with PI controller
- **Modern Fortran 90+ features** and modular design
- **FFTW3 integration** with OpenMP parallelization
- **Comprehensive analysis tools** (Jupyter notebooks)
- **Automated build system** with background monitoring
## Documentation

### Core References
1. Kim, J., & Moin, P. (1985). Application of a fractional-step method to incompressible Navier-Stokes equations. *Journal of Computational Physics*, 59(2), 308-323.

2. Gottlieb, D., & Orszag, S. A. (1977). *Numerical analysis of spectral methods: theory and applications*. SIAM.

3. Canuto, C., Hussaini, M. Y., Quarteroni, A., & Zang, T. A. (2007). *Spectral methods: evolution to complex geometries and applications to fluid dynamics*. Springer.

### Flow Control Documentation
- **`FLOW_CONTROL_README.md`**: Comprehensive flow control system documentation
- **`turbulent_channel_flow_correlations.ipynb`**: Turbulent flow analysis and correlations
- **`velocity_profile_analysis.ipynb`**: Velocity field analysis tools

### Implementation Notes
- **Flow control methods**: Two approaches (constant pressure/volume) with unified interface
- **PI controller design**: Proportional-integral control with anti-windup protection
- **LGL integration**: Spectral-accurate bulk velocity calculation
- **Real-time monitoring**: Live flow statistics and controller feedback

## License

[Add your chosen license here]

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Contact

[Add your contact information here]

---

*This DNS solver combines classical spectral methods with modern flow control systems, providing a powerful tool for channel flow research and applications.*
