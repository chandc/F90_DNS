! Test wall shear stress calculation
program test_wall_shear
    implicit none
    
    ! Add a simple debug statement right at the beginning of final_output
    write(*,'(A)') ' DEBUG: Entering final_output subroutine'
    
    ! Test basic write operations
    write(*,'(A)') ' DEBUG: Testing write operations'
    write(*,'(A)') ' ============================================'
    write(*,'(A)') ' WALL SHEAR STRESS ANALYSIS (TEST)'
    write(*,'(A)') ' ============================================'
    
    write(*,'(A)') ' DEBUG: Write operations successful'
    
end program test_wall_shear
