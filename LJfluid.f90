! Module which carried out the Lennard-Jones simulation in a periodic box using Monte Carlo 
! algorithm. It contains several subroutines:

! SUBROUTINE 1: Sets the initial geometry of the particles in a cubic box, which size is chosen
!               by the user in main.f90
! SUBROUTINE 2: Computes the potential energy of a determinated geometry configuration with
!               periodic boundary conditions in the three directions
! SUBROUTINE 3: 

! To use this module main.f90 is needed, also contained in this directory

! Author: Inés Sánchez de Movellán Sáiz

module LJfluid

    implicit none
    public :: initial_geom, energy, montecarlo
    contains

    ! *********************************************************************************************
    !                                SUBROUTINE 1: Set the initial geometry
    
    subroutine initial_geom(n,L,geom0)
    
        ! Declaration statements
        integer :: i, j
        integer, intent(in) :: n  ! Number of particles 
        real(kind=8), intent(in) :: L  ! Lenght of the cubic box
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
        geom0 = L*geom0
        !geom0(1,1) = 4.567657
        !geom0(2,1) = 3.299239
        !geom0(3,1) = 1.499908
        !geom0(1,2) = 2.784433
        !geom0(2,2) = 4.387369
        !geom0(3,2) = 2.163670
        !geom0(1,3) = 1.344036
        !geom0(2,3) = 3.873313
        !geom0(3,3) = 2.126799
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
    !                 SUBROUTINE 2: Calculate the potential energy with PBC

    subroutine energy(n, geom, L, rc, V, V_rc)
        
        ! Declaration statements
        integer :: i, j
        integer, intent(in) :: n
        real(kind=8), intent(in) :: L, rc
        real(kind=8), dimension(3,n), intent(in) :: geom
        real(kind=8), intent(out) :: V, V_rc
        real(kind=8), allocatable :: r2(:,:)

        ! Execution zone
        V = 0.d0
        V_rc = 4*(1.d0/rc**12-1.d0/rc**6)

        allocate(r2(n-1,n-1))
        do i = 2, n
            do j = 1, i-1
                if (geom(1,i)-geom(1,j).le.L/2.d0 .and. geom(2,i)-geom(2,j).le.L/2.d0 .and. &
                geom(3,i)-geom(3,j).le.L/2.d0) then
                    !write(*,*) "Todos dentro"
                    r2(i-1,j) = (geom(1,i)-geom(1,j))**2+(geom(2,i)-geom(2,j))**2 &
                    +(geom(3,i)-geom(3,j))**2
                endif
                if (geom(1,i)-geom(1,j).le.L/2.d0 .and. geom(2,i)-geom(2,j).le.L/2.d0 .and. &
                geom(3,i)-geom(3,j).gt.L/2.d0) then
                    !write(*,*) "Z fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j))**2+(geom(2,i)-geom(2,j))**2 &
                    +(geom(3,i)-geom(3,j)-L)**2
                endif
                if (geom(1,i)-geom(1,j).le.L/2.d0 .and. geom(2,i)-geom(2,j).gt.L/2.d0 .and. &
                geom(3,i)-geom(3,j).le.L/2.d0) then
                    !write(*,*) "Y fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j))**2+(geom(2,i)-geom(2,j)-L)**2 &
                    +(geom(3,i)-geom(3,j))**2
                endif
                if (geom(1,i)-geom(1,j).gt.L/2.d0 .and. geom(2,i)-geom(2,j).le.L/2.d0 .and. &
                geom(3,i)-geom(3,j).le.L/2.d0) then
                    !write(*,*) "X fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j)-L)**2+(geom(2,i)-geom(2,j))**2 &
                    +(geom(3,i)-geom(3,j))**2
                endif
                if (geom(1,i)-geom(1,j).le.L/2.d0 .and. geom(2,i)-geom(2,j).gt.L/2.d0 .and. &
                geom(3,i)-geom(3,j).gt.L/2.d0) then
                    !write(*,*) "Y+Z fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j))**2+(geom(2,i)-geom(2,j)-L)**2 &
                    +(geom(3,i)-geom(3,j)-L)**2
                endif
                if (geom(1,i)-geom(1,j).gt.L/2.d0 .and. geom(2,i)-geom(2,j).le.L/2.d0 .and. &
                geom(3,i)-geom(3,j).gt.L/2.d0) then
                    !write(*,*) "X+Z fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j)-L)**2+(geom(2,i)-geom(2,j))**2 &
                    +(geom(3,i)-geom(3,j)-L)**2
                endif
                if (geom(1,i)-geom(1,j).gt.L/2.d0 .and. geom(2,i)-geom(2,j).gt.L/2.d0 .and. &
                geom(3,i)-geom(3,j).le.L/2.d0) then
                    !write(*,*) "X+Y fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j)-L)**2+(geom(2,i)-geom(2,j)-L)**2 &
                    +(geom(3,i)-geom(3,j))**2
                endif
                if (geom(1,i)-geom(1,j).gt.L/2.d0 .and. geom(2,i)-geom(2,j).gt.L/2.d0 .and. &
                geom(3,i)-geom(3,j).gt.L/2.d0) then
                    !write(*,*) "Todos fuera"
                    r2(i-1,j) = (geom(1,i)-geom(1,j)-L)**2+(geom(2,i)-geom(2,j)-L)**2 &
                    +(geom(3,i)-geom(3,j)-L)**2
                endif

                ! Cutoff
                if (r2(i-1,j)<rc**2) then
                    V = V + (1.d0/r2(i-1,j)**6-1.d0/r2(i-1,j)**3) - V_rc
                endif
            enddo
        enddo
        V = 4.d0*V

        return

    end subroutine energy

    ! *********************************************************************************************

    ! *********************************************************************************************
    !                 SUBROUTINE 3: Simulation of the LJ fluid with MC tecniques
    
    subroutine montecarlo(n, geom0, L, rc, V, V_rc, maxcycle)

        ! Declaration statements
        integer :: i,j, counter, atom
        integer, intent(in) :: n, maxcycle
        real(kind=8), dimension(3,n) :: geom0, geomi
        real(kind=8) :: L, rc, V, V_rc, V_new, aux, delta_r, delta_V, a

        ! Execution zone
        ! The potential energy for the initial geometry is calculated
        call energy(n, geom0, L, rc, V, V_rc)
        open(24, file="out", action="write")
        open(25, file="v_out", action="write")

        ! Initialization of the counter of MC cycles
        counter = 0
        do while (counter .lt. maxcycle) 
            counter = counter + 1
            write(24,'( "CYCLE", I5)') counter
            write(24,'( "Potential V", f10.6)') V 
            write(25,*) counter, V 
            ! A random atom is chosen and displaced a random quantity between -0.5 and 0.5 units
            call random_number(aux)
            atom = int(1+aux*n)
            call random_number(delta_r)
            delta_r = delta_r-0.5
            ! The geometry in cycle i is geomi, i.e., the initial geometry with the random 
            ! displacement of atom i
            geomi = geom0
            geomi(:,atom) = geomi(:,atom)+delta_r
            write(24,'(3f10.5)') geom0
            write(24,*) " "
            write(24,'(3f10.5)') geomi
            ! Condition of restricted coordinates to a simulation box
            do i = 1, 3
                if (geomi(i,atom) .lt. 0.d0) then
                    geomi(i,atom) = geomi(i,atom) + L
                elseif (geomi(i,atom) .gt. L) then 
                    geomi(i,atom) = geomi(i,atom) - L
                endif
            enddo
            call energy(n, geomi, L, rc, V_new, V_rc) 
            write(24,'( "Potential V_new", f10.6)') V_new
            delta_V = V_new - V
            write(24,'( "Delta V", f10.6)') delta_V
            if (delta_V .lt. 0.d0) then
                write(24,*) "Accepted geometry"
                geom0 = geomi
                V = V_new
            else 
                call random_number(a)
                write(24, '( "Random a", f15.9 )') a
                write(24, '( "Boltzmann factor", f15.9 )') exp(-delta_V)
                if (a .lt. exp(-delta_V)) then
                    geom0 = geomi
                    V = V_new
                    write(24,*) "Accepted geometry"
                else
                    write(24,*) "Rejected geometry"
                endif
            endif
            write(24,*) " "
        enddo
        return

    end subroutine montecarlo

    ! *********************************************************************************************


end module LJfluid
