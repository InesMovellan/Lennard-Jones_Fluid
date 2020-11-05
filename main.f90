! Main code that uses the module LJfluid to compute several properties of a fluid of Lennard-Jones
! particles using Monte Carlo techniques and periodic boundary conditions

! Author: Inés Sánchez de Movellán Sáiz

program main

    use LJfluid

    ! Declaration statements

    implicit none
    integer :: n, L, maxcycle
    real(kind=8), allocatable :: coord(:,:)
    real(kind=8) :: threshold

    ! --------------------------------------------------------------------------------------------

    ! Execution zone
    n = 3
    L = 5
    call initial_geom(n,L,coord)

    maxcycle = 3
    threshold = 10.d0**(-6)
    call energy(n,coord,maxcycle,threshold)
end program main
