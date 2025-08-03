# =============================================================================
# Makefile for 3D Navier-Stokes Channel Flow DNS Solver
# Originally developed by Daniel Chiu-Leung Chan, 1993
# Modernized F90 version with FFTW and OpenMP support
# =============================================================================

# Compiler settings
FC = gfortran
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -fopenmp -Wall -Wextra
DEBUG_FFLAGS = -g -O0 -fdefault-real-8 -fdefault-double-8 -fbounds-check -fbacktrace -fopenmp -Wall -Wextra

# FFTW libraries (macOS Homebrew paths)
FFTW_INCLUDE = -I/opt/homebrew/include
FFTW_LIBDIR = -L/opt/homebrew/lib
FFTW_LIBS = $(FFTW_LIBDIR) -lfftw3 -lfftw3_omp -llapack -lblas -lm

# Source files and dependencies
LGL_MODULE = lgl_module.f90
FFT_MODULE = fft_module.f90
MAIN_SOURCE = base_code_complete.f90

# Object files
LGL_OBJ = lgl_module.o
FFT_OBJ = fft_module.o
MAIN_OBJ = base_code_complete.o
OBJECTS = $(LGL_OBJ) $(FFT_OBJ) $(MAIN_OBJ)

# Executable name
TARGET = dns_solver

# Module files (auto-generated)
MOD_FILES = lgl_module.mod fft_module.mod

# Default target
all: $(TARGET)
	@echo "=============================================="
	@echo "  DNS Solver compiled successfully!"
	@echo "  Executable: ./$(TARGET)"
	@echo "  Run with: make run"
	@echo "=============================================="

# Build executable (link all objects)
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(FFTW_LIBS)

# Module dependencies and compilation rules
$(LGL_OBJ): $(LGL_MODULE)
	$(FC) $(FFLAGS) $(FFTW_INCLUDE) -c $<

$(FFT_OBJ): $(FFT_MODULE)
	$(FC) $(FFLAGS) $(FFTW_INCLUDE) -c $<

$(MAIN_OBJ): $(MAIN_SOURCE) $(LGL_OBJ) $(FFT_OBJ)
	$(FC) $(FFLAGS) $(FFTW_INCLUDE) -c $<

# Run the DNS solver
run: $(TARGET)
	@echo "Running 3D Navier-Stokes Channel Flow Solver..."
	./$(TARGET)

# Debug build with full error checking
debug: FFLAGS = $(DEBUG_FFLAGS)
debug: clean $(TARGET)
	@echo "Debug version compiled with bounds checking and backtrace"

# Clean up generated files
clean:
	rm -f $(TARGET) *.o *.mod *.out

# Create sample input file if it doesn't exist
input.dat:
	@echo "Creating sample input.dat..."
	@echo "&inputs" > input.dat
	@echo " istart=0, dt=0.01, nsteps=100, nwrt=10, iform=0, iles=0" >> input.dat
	@echo "/" >> input.dat
	@echo "&params" >> input.dat
	@echo " xlen=6.28318, ylen=1.0, re=180.0, ta=0.0, ybar=1.0," >> input.dat
	@echo " cgstol=1.0e-4, cs=0.1, u00=0.0, wavlen=1.0" >> input.dat
	@echo "/" >> input.dat

# Test compilation and basic run
test: $(TARGET) input.dat
	@echo "Testing DNS solver compilation and basic functionality..."
	./$(TARGET)

# Parallel run options
run-omp2:
	@echo "Running with 2 OpenMP threads..."
	OMP_NUM_THREADS=2 ./$(TARGET)

run-omp4:
	@echo "Running with 4 OpenMP threads..."
	OMP_NUM_THREADS=4 ./$(TARGET)

run-omp8:
	@echo "Running with 8 OpenMP threads..."
	OMP_NUM_THREADS=8 ./$(TARGET)

# Run simulation in background with logging
run-bg: $(TARGET)
	@echo "Starting DNS solver in background with logging..."
	@echo "Output will be saved to simulation.log"
	@echo "Use 'make monitor' to watch progress"
	@echo "Use 'make stop' to stop the simulation"
	./$(TARGET) > simulation.log 2>&1 &
	@echo "Simulation started with PID: $$!"

# Monitor simulation progress automatically
monitor:
	@echo "Monitoring simulation progress (Ctrl+C to stop monitoring)..."
	@echo "=========================================================="
	@while ps aux | grep -v grep | grep dns_solver > /dev/null; do \
		if [ -f simulation.log ]; then \
			echo "=== $$(date) ==="; \
			tail -3 simulation.log | grep "Time =" | tail -1; \
			echo ""; \
		fi; \
		sleep 30; \
	done; \
	echo "Simulation completed or stopped."

# Stop running simulation
stop:
	@echo "Stopping DNS solver simulation..."
	@pkill -f dns_solver || echo "No running simulation found"

# Check simulation status
status:
	@if ps aux | grep -v grep | grep dns_solver > /dev/null; then \
		echo "✓ Simulation is running"; \
		if [ -f simulation.log ]; then \
			echo "Latest progress:"; \
			tail -3 simulation.log | grep "Time =" | tail -1; \
		fi; \
	else \
		echo "✗ No simulation running"; \
	fi

# Check dependencies
check-deps:
	@echo "Checking compilation dependencies..."
	@which gfortran > /dev/null && echo "✓ gfortran found" || echo "✗ gfortran not found"
	@pkg-config --exists fftw3 && echo "✓ FFTW3 found" || echo "✗ FFTW3 not found"
	@echo "FFTW include path: $(FFTW_INCLUDE)"
	@echo "FFTW library path: $(FFTW_LIBDIR)"

# Display help
help:
	@echo "Available targets:"
	@echo "  all        - Compile the DNS solver (default)"
	@echo "  run        - Compile and run the solver interactively"
	@echo "  run-bg     - Run solver in background with logging"
	@echo "  monitor    - Monitor background simulation progress (30s intervals)"
	@echo "  status     - Check if simulation is running and show progress"
	@echo "  stop       - Stop background simulation"
	@echo "  debug      - Compile with debug flags and error checking"
	@echo "  test       - Compile, create input.dat, and test run"
	@echo "  clean      - Remove all generated files"
	@echo "  run-omp2   - Run with 2 OpenMP threads"
	@echo "  run-omp4   - Run with 4 OpenMP threads"
	@echo "  run-omp8   - Run with 8 OpenMP threads"
	@echo "  input.dat  - Create sample input file"
	@echo "  check-deps - Check for required dependencies"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Background simulation workflow:"
	@echo "  1. make run-bg     # Start simulation in background"
	@echo "  2. make monitor    # Watch progress automatically"
	@echo "  3. make stop       # Stop when needed"
	@echo ""
	@echo "Source files: $(LGL_MODULE) $(FFT_MODULE) $(MAIN_SOURCE)"
	@echo "Target executable: $(TARGET)"

.PHONY: all run run-bg monitor stop status debug test clean run-omp2 run-omp4 run-omp8 check-deps help
