! Main code that uses the module LJfluid to compute several properties of a fluid of Lennard-Jones
! particles using Monte Carlo techniques and periodic boundary conditions

! Author: Inés Sánchez de Movellán Sáiz

program main

    use LJfluid

    ! Declaration statements

    implicit none
    integer :: n, maxcycle
    real(kind=8), allocatable :: coord(:,:)
    real(kind=8) :: threshold, L

    ! --------------------------------------------------------------------------------------------

    ! Execution zone
    n = 3
    L = 6
    call initial_geom(n,L,coord)

    maxcycle = 1
    threshold = 10.d0**(-6)
    call energy(n,coord,L,maxcycle,threshold)
end program main
