! Main code that uses the module LJfluid to compute several properties of a fluid of Lennard-Jones
! particles using Monte Carlo techniques and periodic boundary conditions

! Author: Inés Sánchez de Movellán Sáiz

program main

    use LJfluid

    ! Declaration statements

    implicit none
    integer :: n, maxcycle
    real(kind=8), allocatable :: coord(:,:)
    real(kind=8) :: threshold, L, rc, V, V_rc

    ! --------------------------------------------------------------------------------------------

    ! Execution zone
    n = 5
    L = 6.d0
    call initial_geom(n,L,coord)

    rc = 3.d0
    !call energy(n,coord,L,rc,V,V_rc)

    maxcycle = 10000
    threshold = 10.d0**(-6)
    call montecarlo(n,coord,L,rc,V,V_rc, maxcycle)
    
end program main
