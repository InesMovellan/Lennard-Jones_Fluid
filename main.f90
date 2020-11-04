! Main code that uses the module LJfluid to compute several properties of a fluid of Lennard-Jones
! particles using Monte Carlo techniques and periodic boundary conditions

! Author: Inés Sánchez de Movellán Sáiz

program main

    use LJfluid

    ! Declaration statements

    implicit none
    integer :: n, L
    real(kind=8), allocatable :: x(:,:)
    ! real(kind=8), allocatable :: A(:,:)               
    ! real(kind=8) :: R2                                
    ! character(50) :: Aformat, cformat                 

    ! ----------------------------------------------------------------------------------------------

    ! Execution zone
    n = 10
    L = 10
    call initial_geom(n,L,x)
  
    ! We open the out.dat file, in which we will write the results
    !open(23, file="out.dat", action="write")
    ! The title of the output is written
    !write(23,*) " "
    
    ! The subroutine read_data is called, which reads the values of points (x,y) from input.dat file
    !call read_data(n,x,y)
    !do i = 1, n
    !    write(23,'(2f24.16)') x(i), y(i)
    !enddo

end program main
