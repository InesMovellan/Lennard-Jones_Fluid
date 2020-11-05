! Module which carried out the Lennard-Jones simulation in a periodic box using Monte Carlo 
! algorithm. It contains several subroutines:

! SUBROUTINE 1: 
! SUBROUTINE 2: 

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
        !call random_number(geom0)
        !geom0 = L*geom0
        geom0(1,1) = 4.567657
        geom0(2,1) = 3.299239
        geom0(3,1) = 1.499908
        geom0(1,2) = 2.784433
        geom0(2,2) = 4.387369
        geom0(3,2) = 2.163670
        geom0(1,3) = 1.344036
        geom0(2,3) = 3.873313
        geom0(3,3) = 2.126799
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
    !                 SUBROUTINE 2: Calculate the potential energy
    
    subroutine energy(n,geom0,V)

        ! Declaration statements
        integer :: i,j, maxcycle
        integer, intent(in) :: n
        real(kind=8), dimension(3,n), intent(in) :: geom0
        real(kind=8), allocatable :: r12(:,:), r6(:,:)
        real(kind=8), intent(inout) :: V
        real(kind=8) :: threshold
        character(20) :: r_format

        ! Execution zone
        allocate(r12(n-1,n-1))
        allocate(r6(n-1,n-1))
        ! The intial values of both potential energy and powers of the distance are defined as 
        ! zero
        r12 = 0.d0
        r6 = 0.d0
        V = 0.d0
        
        ! With this loop powers of distance and Lennard-Jones energy are calculated iteratively
        do i = 2, n
            do j = 1, i-1
                r12(i-1,j) = ((geom0(1,i)-geom0(1,j))**2+(geom0(2,i)-geom0(2,j))**2 &
                    +(geom0(3,i)-geom0(3,j))**2)**6
                r6(i-1,j) = ((geom0(1,i)-geom0(1,j))**2+(geom0(2,i)-geom0(2,j))**2 &
                    +(geom0(3,i)-geom0(3,j))**2)**3
                V = V + 1.d0/r12(i-1,j) - 1.d0/r6(i-1,j)
            enddo
        enddo
        V = 4.d0*V

        write(*,*) V
        !write(r_format, '(a, i0, a)') '(', n-1, 'f17.5)'
        !write(*,*) " "
        !write(*,r_format) r2
        !write(*,*) " "
        !write(*,r_format) r12
        !write(*,*) " "
        !write(*,r_format) r6
        return

    end subroutine energy

    ! *********************************************************************************************

end module LJfluid
