! Module which carried out the Lennard-Jones simulation in a periodic box using Monte Carlo 
! algorithm. It contains several subroutines:

! SUBROUTINE 1: 

! To use this module main.f90 is needed, also contained in this directory

! Author: Inés Sánchez de Movellán Sáiz

module LJfluid

    implicit none
    public :: initial_geom, energy
    contains

    ! *********************************************************************************************
    !                                SUBROUTINE 1: Set the initial geometry
    
    subroutine initial_geom(n,L,geom0)
    
        ! Declaration statements
        integer :: i, j
        integer, intent(in) :: n, L  ! Number of particles (n) and lenght of the cubic box (L)
        real(kind=8), allocatable, intent(inout) :: geom0(:,:)
        real(kind=8), allocatable :: xyz(:)

        ! Execution zone
        ! The initial coordinates of the particles are stored in a matrix called pos0, which
        ! has the following structure:
        ! x1  y1  z1
        ! x2  y2  z2
        ! x3  y3  z3
        !     :
        ! xn  yn  zn
        allocate(geom0(3,n))
        call random_number(geom0)
        geom0(1,1) = 1.
        geom0(2,1) = 0.
        geom0(3,1) = 0.
        geom0(1,2) = 1.
        geom0(2,2) = 1.
        geom0(3,2) = 1.
        geom0(1,3) = 1.
        geom0(2,3) = 3.
        geom0(3,3) = 3.

        geom0 = L*geom0
        ! The initial geometry is written in a .xyz file, which can be read with 3D visualization
        ! programs (vesta, vmd, avogadro). In order to do that we use a auxiliar vector xyz
        open(23, file="initial_geom.xyz", action="write")
        write(23,'(I4)') n
        write(23,*) " "
        allocate(xyz(3))
        do i = 1, n
            xyz = geom0(:,i)
            write(23,'( "H", f10.6, f10.6, f10.6)' ) xyz(1), xyz(2), xyz(3)
        enddo
        return

    end subroutine initial_geom

    ! *********************************************************************************************

    ! *********************************************************************************************
    !                                SUBROUTINE 2: Compute the energy
    
    subroutine energy(n, geom0, r12, r)

        ! Declaration statements
        integer :: i,j
        integer, intent(in) :: n
        real(kind=8), dimension(3,n), intent(in) :: geom0
        real(kind=8), allocatable, intent(inout) :: r(:,:)
        real(kind=8) :: r12

        ! Execution zone
        allocate(r(n-1,n-1))
        r = 0.d0
        do i = 2, n
            do j = 1, i-1
                r(i-1,j) = (geom0(1,i)-geom0(1,j))**2+(geom0(2,i)-geom0(2,j))**2+(geom0(3,i)-geom0(3,j))**2
                write(*,*) r(i-1,j)

            enddo
        enddo
        write(*,*) " "
        write(*,'(2f10.6)') r
        r12 = (geom0(1,2)-geom0(1,1))**2+(geom0(2,2)- geom0(2,1))**2+(geom0(3,2)- geom0(3,1))**2
        !write(*,*) r12
        !write(*,*) r

    end subroutine energy

    ! *********************************************************************************************
end module LJfluid
