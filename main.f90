! Main code that uses the module LJfluid to compute several properties of a fluid of Lennard-Jones
! particles using Monte Carlo techniques and periodic boundary conditions

! Author: Inés Sánchez de Movellán Sáiz

program main

    use LJfluid

    ! Declaration statements

    implicit none
    integer :: n, L
    real(kind=8), allocatable :: coord(:,:)
    real(kind=8), allocatable :: dist(:)
    real(kind=8) :: r

    ! --------------------------------------------------------------------------------------------

    ! Execution zone
    n = 3
    L = 1
    call initial_geom(n,L,coord)

    call energy(n,coord,r,dist)

end program main
