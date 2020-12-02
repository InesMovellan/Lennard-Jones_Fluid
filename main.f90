! Main code that uses the module LJfluid to compute several properties of a fluid of Lennard-Jones
! particles using Monte Carlo techniques and periodic boundary conditions

! Author: Inés Sánchez de Movellán Sáiz

program main

    use LJfluid

    ! Declaration statements

    implicit none
    integer :: n, maxcycle
    real(kind=8), allocatable :: coord(:,:)
    real(kind=8) :: L, rc, V, V_rc, T

    ! --------------------------------------------------------------------------------------------

    ! Execution zone
    n = 100
    L = 6.d0
    call initial_geom(n, L, coord)

    rc = 3.d0
    T = 1.268
    maxcycle = 100000
    call montecarlo(n, coord, L, rc, V, V_rc, T, maxcycle)
    
end program main
