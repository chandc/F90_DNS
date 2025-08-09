program test_lgl_only
    use lgl_module
    implicit none
    
    real(8) :: z(5), weights(5), D1(5,5)
    
    write(*,*) 'Testing LGL module...'
    call lgl_nodes_weights(5, z, weights)
    call differentiation_matrix(5, z, D1)
    write(*,*) 'LGL test completed!'
    
end program test_lgl_only
