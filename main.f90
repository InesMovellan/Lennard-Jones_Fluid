! Main code that uses the module LJfluid.f90 to carry out a simulation of the Lennard-Jones fluid
! of N particles at volume V = L^3 and temperature T, i.e., in the canonical ensemble, using
! Monte Carlo techniques.

! Author: Inés Sánchez de Movellán Sáiz

program main

    ! To make use of the module LJfluid use statement is employed. 
    use LJfluid

    ! Declaration statements
    implicit none
    ! n = number of particles in the simulation box
    ! cycles = number of Monte Carlo steps
    ! therm = number of Monte Carlo steps that will be skipped for the calculation of g(r)
    ! coord(:,:) = 3 columns N rows matrix which stores (X, Y, Z) coordinates of each particle
    ! L = simulation box side length
    ! rc = cutoff radius from which the interaction between two atoms is neglected
    ! V = potential energy of the system
    ! Vrc = potential energy of the system at r = rc
    ! T = temperature
    integer :: n, cycles, therm
    real(kind=8), allocatable :: coord(:,:)
    real(kind=8) :: L, rc, V, Vrc, T


    ! Execution zone
    ! The user set the values of n, L rc, T and cycles. Magnitudes (L, rc and T) are in reduce
    ! units
    n = 100
    L = 6.d0
    rc = 3.d0
    T = 1.268
    cycles = 10**7
    therm = 10**6

    ! Call initial_geom subroutine of LJfluid to set the initial coord(:,:) matrix
    call initial_geom(n, L, coord)

    ! Call montecarlo subroutine of LJfluid to perform the simulation of the fluid
    call montecarlo(n, coord, L, rc, V, Vrc, T, cycles,therm)

    
end program main
