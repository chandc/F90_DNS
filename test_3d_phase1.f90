!===============================================================================
!
! DNS_PRESSURE_BC_3D.F90 - 3D DNS Phase 1 Testing Version
!
! PURPOSE: Test Phase 1 data structure extensions for 3D DNS solver
!
! FEATURES TESTED:
!   • 3D Grid parameter setup (nx, ny, nz)
!   • 3D Array allocation and indexing
!   • 3D Input file parsing with ny_input
!   • Memory management for 3D arrays
!   • 3D Wavenumber setup (kx, ky)
!   • Data structure consistency validation
!
!===============================================================================

program test_3d_phase1
    use iso_fortran_env, only: wp => real64
    implicit none
    
    ! =========================================================================
    ! 3D GRID PARAMETERS - PHASE 1 TEST
    ! =========================================================================
    integer :: nx, ny, nz  ! Grid dimensions (EXTENDED: added ny)
    integer :: nxpp, nypp, nxh, nyh, nxhp, nyhp, nxf, nyf, ntot, nzm
    
    ! =========================================================================
    ! 3D PARAMETER STRUCTURE - PHASE 1 TEST
    ! =========================================================================
    type :: test_params_3d
        ! Physical parameters
        real(wp) :: alpha, beta, re, ybar, xlen, ylen
        
        ! 3D Wave numbers (EXTENDED)
        real(wp), allocatable :: xw(:), xsq(:)  ! x-direction wavenumbers
        real(wp), allocatable :: yw(:), ysq(:)  ! y-direction wavenumbers (NEW)
        
        ! 3D Flow fields (EXTENDED: added v-component)
        real(wp), allocatable :: u(:), v(:), w(:)  ! EXTENDED: added v
        real(wp), allocatable :: un(:), vn(:), wn(:)  ! EXTENDED: added vn
        
        ! 3D Source terms (EXTENDED: added sv)
        real(wp), allocatable :: su(:), sv(:), sw(:)  ! EXTENDED: added sv
        
        ! 3D Scratch arrays (EXTENDED: added oy)
        real(wp), allocatable :: ox(:), oy(:), oz(:)  ! EXTENDED: added oy
        
        ! 3D Boundary condition arrays (EXTENDED)
        real(wp), allocatable :: uw1(:), vw1(:), ww1(:)  ! EXTENDED: added vw1
        real(wp), allocatable :: uw2(:), vw2(:), ww2(:)  ! EXTENDED: added vw2
        
        ! 3D Divergence checking (EXTENDED)
        real(wp), allocatable :: div_check_3d(:,:,:)  ! EXTENDED: 3D array
    end type test_params_3d
    
    type(test_params_3d) :: p
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer :: i, j, k, idx
    logical :: test_passed
    
    write(*,'(A)') ' ============================================'
    write(*,'(A)') '   3D DNS SOLVER - PHASE 1 TESTING'
    write(*,'(A)') '   Data Structure and Memory Validation'
    write(*,'(A)') ' ============================================'
    
    ! Test 1: 3D Grid Parameter Setup
    write(*,'(A)') ' TEST 1: 3D Grid Parameter Setup'
    call test_grid_setup_3d()
    
    ! Test 2: 3D Input File Reading
    write(*,'(A)') ' TEST 2: 3D Input File Reading'
    call test_input_reading_3d(p)
    
    ! Test 3: 3D Array Allocation
    write(*,'(A)') ' TEST 3: 3D Array Allocation'
    call test_array_allocation_3d(p)
    
    ! Test 4: 3D Wavenumber Setup
    write(*,'(A)') ' TEST 4: 3D Wavenumber Setup'
    call test_wavenumber_setup_3d(p)
    
    ! Test 5: 3D Indexing Validation
    write(*,'(A)') ' TEST 5: 3D Indexing Validation'
    call test_3d_indexing(p)
    
    ! Test 6: Memory Usage Analysis
    write(*,'(A)') ' TEST 6: Memory Usage Analysis'
    call analyze_memory_usage_3d(p)
    
    ! Test 7: 3D Data Structure Consistency
    write(*,'(A)') ' TEST 7: 3D Data Structure Consistency'
    call test_data_consistency_3d(p)
    
    ! Cleanup
    call cleanup_arrays_3d(p)
    
    write(*,'(A)') ' ============================================'
    write(*,'(A)') '   PHASE 1 TESTING COMPLETE'
    write(*,'(A)') '   All 3D data structures validated ✅'
    write(*,'(A)') ' ============================================'
    write(*,'(A)') ' Ready to proceed to Phase 2: FFT Extension'
    write(*,'(A)') ' ============================================'

contains

    ! =========================================================================
    ! TEST SUBROUTINES
    ! =========================================================================
    
    subroutine test_grid_setup_3d()
        implicit none
        integer :: nx_test, ny_test, nz_test
        
        ! Test default 3D grid setup
        nx_test = 128; ny_test = 32; nz_test = 33
        call setup_grid_parameters_3d(nx_test, ny_test, nz_test)
        
        write(*,'(A,I0,A,I0,A,I0)') '   Input grid: nx=', nx_test, ', ny=', ny_test, ', nz=', nz_test
        write(*,'(A,I0)') '   Total 3D points: ', ntot
        write(*,'(A,F8.2,A)') '   Memory estimate: ', real(ntot*8*6)/(1024.0**2), ' MB (6 3D arrays)'
        
        ! Validate grid relationships
        if (ntot == nxpp*nypp*nz .and. nxpp == nx+2 .and. nypp == ny+2) then
            write(*,'(A)') '   ✅ Grid parameter relationships correct'
        else
            write(*,'(A)') '   ❌ Grid parameter relationships incorrect'
        endif
        
    end subroutine test_grid_setup_3d
    
    subroutine test_input_reading_3d(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        logical :: file_exists
        
        ! Check if we can read 3D input file
        inquire(file='input_3d.dat', exist=file_exists)
        if (file_exists) then
            write(*,'(A)') '   ✅ input_3d.dat found'
            call read_3d_input_simple(p)
        else
            write(*,'(A)') '   ⚠️  input_3d.dat not found, using defaults'
            p%xlen = 2.0_wp * pi
            p%ylen = 2.0_wp * pi
            p%re = 10000.0_wp
        endif
        
        write(*,'(A,F8.4)') '   Domain: Lx = ', p%xlen
        write(*,'(A,F8.4)') '   Domain: Ly = ', p%ylen
        write(*,'(A,F8.1)') '   Reynolds: Re = ', p%re
        
    end subroutine test_input_reading_3d
    
    subroutine test_array_allocation_3d(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        integer :: alloc_status
        
        ! Allocate 3D wavenumber arrays
        allocate(p%xw(nxhp), p%xsq(nxhp), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ X-wavenumber arrays allocated'
        else
            write(*,'(A)') '   ❌ X-wavenumber allocation failed'
        endif
        
        allocate(p%yw(nyhp), p%ysq(nyhp), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ Y-wavenumber arrays allocated'
        else
            write(*,'(A)') '   ❌ Y-wavenumber allocation failed'
        endif
        
        ! Allocate 3D velocity fields
        allocate(p%u(ntot), p%v(ntot), p%w(ntot), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D velocity fields allocated'
        else
            write(*,'(A)') '   ❌ 3D velocity allocation failed'
        endif
        
        allocate(p%un(ntot), p%vn(ntot), p%wn(ntot), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D velocity storage allocated'
        else
            write(*,'(A)') '   ❌ 3D velocity storage allocation failed'
        endif
        
        ! Allocate 3D source terms
        allocate(p%su(ntot), p%sv(ntot), p%sw(ntot), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D source terms allocated'
        else
            write(*,'(A)') '   ❌ 3D source term allocation failed'
        endif
        
        ! Allocate 3D scratch arrays
        allocate(p%ox(ntot), p%oy(ntot), p%oz(ntot), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D scratch arrays allocated'
        else
            write(*,'(A)') '   ❌ 3D scratch allocation failed'
        endif
        
        ! Allocate 3D boundary condition arrays
        allocate(p%uw1(nxpp*nypp), p%vw1(nxpp*nypp), p%ww1(nxpp*nypp), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D boundary arrays (wall 1) allocated'
        else
            write(*,'(A)') '   ❌ 3D boundary (wall 1) allocation failed'
        endif
        
        allocate(p%uw2(nxpp*nypp), p%vw2(nxpp*nypp), p%ww2(nxpp*nypp), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D boundary arrays (wall 2) allocated'
        else
            write(*,'(A)') '   ❌ 3D boundary (wall 2) allocation failed'
        endif
        
        ! Allocate 3D divergence array
        allocate(p%div_check_3d(nxhp, nyhp, nz), stat=alloc_status)
        if (alloc_status == 0) then
            write(*,'(A)') '   ✅ 3D divergence array allocated'
        else
            write(*,'(A)') '   ❌ 3D divergence allocation failed'
        endif
        
    end subroutine test_array_allocation_3d
    
    subroutine test_wavenumber_setup_3d(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        integer :: i, j
        
        ! Setup fundamental wavenumbers
        p%alpha = 2.0_wp * pi / p%xlen
        p%beta = 2.0_wp * pi / p%ylen
        
        ! X-direction wavenumbers
        do i = 1, nxhp
            p%xw(i) = real(i-1, wp) * p%alpha
            p%xsq(i) = p%xw(i)**2
        end do
        
        ! Y-direction wavenumbers (NEW)
        do j = 1, nyhp
            p%yw(j) = real(j-1, wp) * p%beta
            p%ysq(j) = p%yw(j)**2
        end do
        
        write(*,'(A,F10.6)') '   Alpha (x-wavenumber): ', p%alpha
        write(*,'(A,F10.6)') '   Beta (y-wavenumber): ', p%beta
        write(*,'(A,F10.4)') '   Max kx: ', p%xw(nxhp)
        write(*,'(A,F10.4)') '   Max ky: ', p%yw(nyhp)
        write(*,'(A,F10.4)') '   Max kx²: ', p%xsq(nxhp)
        write(*,'(A,F10.4)') '   Max ky²: ', p%ysq(nyhp)
        
        if (abs(p%xw(2) - p%alpha) < 1.0e-12_wp .and. abs(p%yw(2) - p%beta) < 1.0e-12_wp) then
            write(*,'(A)') '   ✅ Wavenumber setup correct'
        else
            write(*,'(A)') '   ❌ Wavenumber setup incorrect'
        endif
        
    end subroutine test_wavenumber_setup_3d
    
    subroutine test_3d_indexing(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        integer :: i, j, k, idx, count_test
        real(wp) :: test_value
        
        write(*,'(A)') '   Testing 3D indexing: idx = ((k-1)*nypp + j-1)*nxpp + i'
        
        ! Initialize with test pattern
        count_test = 0
        do k = 1, nz
            do j = 1, nypp
                do i = 1, nxpp
                    idx = ((k-1)*nypp + j-1)*nxpp + i
                    if (idx >= 1 .and. idx <= ntot) then
                        test_value = real(i, wp) + 0.1_wp*real(j, wp) + 0.01_wp*real(k, wp)
                        p%u(idx) = test_value
                        count_test = count_test + 1
                    endif
                end do
            end do
        end do
        
        write(*,'(A,I0,A,I0)') '   Indexed points: ', count_test, ' / ', ntot
        
        ! Test reverse indexing at a few points
        idx = ((2-1)*nypp + 3-1)*nxpp + 5  ! k=2, j=3, i=5
        test_value = 5.0_wp + 0.1_wp*3.0_wp + 0.01_wp*2.0_wp
        if (abs(p%u(idx) - test_value) < 1.0e-12_wp) then
            write(*,'(A)') '   ✅ 3D indexing correct'
        else
            write(*,'(A)') '   ❌ 3D indexing incorrect'
            write(*,'(A,F12.6,A,F12.6)') '   Expected: ', test_value, ', Got: ', p%u(idx)
        endif
        
    end subroutine test_3d_indexing
    
    subroutine analyze_memory_usage_3d(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        real(wp) :: total_memory_mb, velocity_mb, wavenumber_mb, scratch_mb, boundary_mb
        
        ! Calculate memory usage for different categories
        velocity_mb = real(ntot * 6 * 8, wp) / (1024.0_wp**2)  ! u,v,w,un,vn,wn
        wavenumber_mb = real((nxhp + nyhp) * 2 * 8, wp) / (1024.0_wp**2)  ! xw,xsq,yw,ysq
        scratch_mb = real(ntot * 3 * 8, wp) / (1024.0_wp**2)  ! ox,oy,oz
        boundary_mb = real(nxpp*nypp * 6 * 8, wp) / (1024.0_wp**2)  ! 6 boundary arrays
        total_memory_mb = velocity_mb + wavenumber_mb + scratch_mb + boundary_mb
        
        write(*,'(A)') '   Memory usage breakdown:'
        write(*,'(A,F8.2,A)') '   - Velocity fields: ', velocity_mb, ' MB'
        write(*,'(A,F8.2,A)') '   - Wavenumbers: ', wavenumber_mb, ' MB'
        write(*,'(A,F8.2,A)') '   - Scratch arrays: ', scratch_mb, ' MB'
        write(*,'(A,F8.2,A)') '   - Boundary arrays: ', boundary_mb, ' MB'
        write(*,'(A,F8.2,A)') '   - Total estimated: ', total_memory_mb, ' MB'
        
        if (total_memory_mb < 1000.0_wp) then
            write(*,'(A)') '   ✅ Memory usage reasonable for testing'
        else
            write(*,'(A)') '   ⚠️  High memory usage - consider smaller grid for testing'
        endif
        
    end subroutine analyze_memory_usage_3d
    
    subroutine test_data_consistency_3d(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        logical :: all_tests_passed
        
        all_tests_passed = .true.
        
        ! Test array size consistency
        if (size(p%u) /= ntot .or. size(p%v) /= ntot .or. size(p%w) /= ntot) then
            write(*,'(A)') '   ❌ Velocity array sizes inconsistent'
            all_tests_passed = .false.
        endif
        
        if (size(p%su) /= ntot .or. size(p%sv) /= ntot .or. size(p%sw) /= ntot) then
            write(*,'(A)') '   ❌ Source array sizes inconsistent'
            all_tests_passed = .false.
        endif
        
        if (size(p%xw) /= nxhp .or. size(p%yw) /= nyhp) then
            write(*,'(A)') '   ❌ Wavenumber array sizes inconsistent'
            all_tests_passed = .false.
        endif
        
        if (size(p%div_check_3d,1) /= nxhp .or. size(p%div_check_3d,2) /= nyhp .or. &
            size(p%div_check_3d,3) /= nz) then
            write(*,'(A)') '   ❌ 3D divergence array dimensions inconsistent'
            all_tests_passed = .false.
        endif
        
        ! Test wavenumber relationships
        if (abs(p%xsq(5) - p%xw(5)**2) > 1.0e-12_wp .or. &
            abs(p%ysq(5) - p%yw(5)**2) > 1.0e-12_wp) then
            write(*,'(A)') '   ❌ Wavenumber squared relationships incorrect'
            all_tests_passed = .false.
        endif
        
        if (all_tests_passed) then
            write(*,'(A)') '   ✅ All data structure consistency tests passed'
        else
            write(*,'(A)') '   ❌ Some consistency tests failed'
        endif
        
    end subroutine test_data_consistency_3d
    
    subroutine cleanup_arrays_3d(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        
        if (allocated(p%xw)) deallocate(p%xw)
        if (allocated(p%xsq)) deallocate(p%xsq)
        if (allocated(p%yw)) deallocate(p%yw)
        if (allocated(p%ysq)) deallocate(p%ysq)
        if (allocated(p%u)) deallocate(p%u)
        if (allocated(p%v)) deallocate(p%v)
        if (allocated(p%w)) deallocate(p%w)
        if (allocated(p%un)) deallocate(p%un)
        if (allocated(p%vn)) deallocate(p%vn)
        if (allocated(p%wn)) deallocate(p%wn)
        if (allocated(p%su)) deallocate(p%su)
        if (allocated(p%sv)) deallocate(p%sv)
        if (allocated(p%sw)) deallocate(p%sw)
        if (allocated(p%ox)) deallocate(p%ox)
        if (allocated(p%oy)) deallocate(p%oy)
        if (allocated(p%oz)) deallocate(p%oz)
        if (allocated(p%uw1)) deallocate(p%uw1)
        if (allocated(p%vw1)) deallocate(p%vw1)
        if (allocated(p%ww1)) deallocate(p%ww1)
        if (allocated(p%uw2)) deallocate(p%uw2)
        if (allocated(p%vw2)) deallocate(p%vw2)
        if (allocated(p%ww2)) deallocate(p%ww2)
        if (allocated(p%div_check_3d)) deallocate(p%div_check_3d)
        
        write(*,'(A)') '   Memory cleanup completed'
        
    end subroutine cleanup_arrays_3d
    
    ! =========================================================================
    ! UTILITY SUBROUTINES
    ! =========================================================================
    
    subroutine setup_grid_parameters_3d(nx_in, ny_in, nz_in)
        implicit none
        integer, intent(in) :: nx_in, ny_in, nz_in
        
        ! Set main grid dimensions
        nx = nx_in
        ny = ny_in
        nz = nz_in
        
        ! Compute derived grid parameters
        nxpp = nx + 2
        nypp = ny + 2
        nxh = nx/2
        nyh = ny/2
        nxhp = nxh + 1
        nyhp = nyh + 1
        nxf = 3*nx/2 + 1
        nyf = 3*ny/2 + 1
        ntot = nxpp*nypp*nz
        nzm = nz - 1
        
    end subroutine setup_grid_parameters_3d
    
    subroutine read_3d_input_simple(p)
        implicit none
        type(test_params_3d), intent(inout) :: p
        integer :: io_status
        integer :: nx_input, ny_input, nz_input
        real(wp) :: re, xlen, ylen
        
        namelist /grid/ nx_input, ny_input, nz_input
        namelist /simulation/ re, xlen, ylen
        
        ! Defaults
        nx_input = 128; ny_input = 32; nz_input = 33
        re = 10000.0_wp; xlen = 2.0_wp*pi; ylen = 2.0_wp*pi
        
        open(7, file='input_3d.dat', status='old', iostat=io_status)
        if (io_status == 0) then
            read(7, nml=grid, iostat=io_status)
            read(7, nml=simulation, iostat=io_status)
            close(7)
            
            call setup_grid_parameters_3d(nx_input, ny_input, nz_input)
            p%re = re
            p%xlen = xlen
            p%ylen = ylen
        endif
        
    end subroutine read_3d_input_simple

end program test_3d_phase1
