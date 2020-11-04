! Module which carried out the Lennard-Jones simulation in a periodic box using Monte Carlo 
! algorithm. It contains several subroutines:

! SUBROUTINE 1: 

! To use this module main.f90 is needed, also contained in this directory

! Author: Inés Sánchez de Movellán Sáiz

module LJfluid

    implicit none
    public :: initial_geom
    contains

    ! ***********************************************************************************************
    !                                SUBROUTINE 1: Set the initial geometry
    
    subroutine initial_geom(n,L,pos)
    
        ! Declaration statements
        integer, intent(in) :: n, L   
        real(kind=8), allocatable, intent(inout) :: pos(:,:)
        character(20) :: pos_format

        ! Execution zone
        allocate(pos(n,n))
        write(pos_format, '(a, i0, a)') '(', n, 'f10.6)'
        call random_number(pos)
        pos = L*pos
        write(*,pos_format) pos

        return

    end subroutine initial_geom

    ! ***********************************************************************************************

    ! ***********************************************************************************************
    !                                SUBROUTINE 2: Compute the energy
    

    ! ***********************************************************************************************
end module LJfluid
