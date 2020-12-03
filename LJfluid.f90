! Module which carried out the Lennard-Jones simulation in a periodic box using Monte Carlo 
! algorithm. It contains several subroutines:
! Serial version of the module.

! SUBROUTINE 1: initial_geom
!               Sets the initial geometry for n Lennard-Jones particles in a cubic box, which 
!               size (and also n) is chosen by the user in main.f90

! SUBROUTINE 2: initial_V
!               Computes the initial potential energy associated with the initial geometry
!               configuration of the particles obtained with initial_geom subroutine. Since the
!               model of the system is a cubic box surrounded by an infinite number of replicas,
!               the calculation of the energy is done with periodic boundary conditions (PBC).

! SUBROUTINE 3: delta_V
!               Calculates the energy change between one trial move performed in montecarlo
!               subroutine and the previous geometry configuration (also with PBC).

! SUBROUTINE 4: montecarlo
!               Monte Carlo subroutine which performs the sampling of structures over the number
!               of cycles selected by the user in the main.f90. The pair correlation function 
!               g(r) is also computed in this subroutine.

! To use this module main.f90 is needed, also contained in this directory

! Author: Inés Sánchez de Movellán Sáiz

module LJfluid

    implicit none
    public :: initial_geom, initial_V, delta_pot, montecarlo
    contains


    ! *********************************************************************************************
    !                          SUBROUTINE 1: Set the initial geometry
    ! *********************************************************************************************
    
    subroutine initial_geom(n,L,geom0)
    
        ! Declaration statements
        implicit none
        ! integers i and j = counters for the do loops
        ! n = number of particles in the simulation box
        ! L = simulation box side length
        ! geom0 = 3 columns n rows matrix which contains the initial set of (X, Y, Z) coordinates 
        integer :: i, j
        integer, intent(in) :: n  
        real(kind=8), intent(in) :: L  
        real(kind=8), allocatable, intent(inout) :: geom0(:,:)


        ! Execution zone
        ! The initial coordinates of the particles are stored in matrix geom0, which has the
        ! following structure:
        ! x1  y1  z1
        ! x2  y2  z2
        ! x3  y3  z3
        !     :
        ! xn  yn  zn
        ! The matrix geom0 is allocated
        allocate(geom0(3,n))
        ! The matrix geom0 is filled with random coordinates of the particles, all of them within
        ! (0, L) range
        call random_number(geom0)
        geom0 = L*geom0

        ! The initial geometry is written in intial_geom.xyz file (following the format of .xyz
        ! files), which can be read with 3D visualization programs.
        open(23, file="initial_geom.xyz", action="write")
        write(23,'(I4)') n
        write(23,*) " "
        do i = 1, n
            write(23,'( "H", f10.6, f10.6, f10.6)' ) geom0(1,i), geom0(2,i), geom0(3,i) 
        enddo
        return

    end subroutine initial_geom

    ! *********************************************************************************************



    ! *********************************************************************************************
    !              SUBROUTINE 2: Calculate the initial potential energy with PBC
    ! *********************************************************************************************

    subroutine initial_V(n, geom, L, r2, rc, V, V_rc)

        ! Declaration statements
        integer :: i, j
        integer, intent(in) :: n
        real(kind=8), dimension(3,n), intent(in) :: geom
        real(kind=8), intent(in) :: L, rc
        real(kind=8), dimension(n-1,n-1), intent(out) :: r2
        real(kind=8), intent(out) :: V, V_rc

        ! Execution zone
        V = 0.d0
        V_rc = 4*(1.d0/rc**12-1.d0/rc**6)

        do i = 2, n
            do j = 1, i-1
                ! Periodic boundary conditions using the nearest image convention

                ! Cartesian coordinate X
                if (abs(geom(1,i)-geom(1,j)).gt.L/2.d0) then
                    r2(i-1,j) = (geom(1,i)-geom(1,j)-L)**2
                else
                    r2(i-1,j) = (geom(1,i)-geom(1,j))**2
                endif

                ! Cartesian coordinate Y
                if (abs(geom(2,i)-geom(2,j)).gt.L/2.d0) then
                    r2(i-1,j) = r2(i-1,j) + (geom(2,i)-geom(2,j)-L)**2
                else
                    r2(i-1,j) = r2(i-1,j) + (geom(2,i)-geom(2,j))**2
                endif

                ! Cartesian coordinate Z
                if (abs(geom(3,i)-geom(3,j)).gt.L/2.d0) then
                    r2(i-1,j) = r2(i-1,j) + (geom(3,i)-geom(3,j)-L)**2
                else
                    r2(i-1,j) = r2(i-1,j) + (geom(3,i)-geom(3,j))**2
                endif

                ! Calculation of the potential using the cutoff rc
                !if (r2(i-1,j)<rc**2) then
                    V = V + (1.d0/r2(i-1,j)**6-1.d0/r2(i-1,j)**3) !- V_rc
                !endif
            enddo
        enddo

        V = 4.d0*V
        !write(*,*) "Matrix r2 initial_V"
        !write(*,'(3f15.10)') r2
        !write(*,*) " "

        return

    end subroutine initial_V

    ! *********************************************************************************************



    ! *********************************************************************************************
    !                 SUBROUTINE 3: Simulation of the LJ fluid with MC tecniques
    ! *********************************************************************************************

    subroutine delta_pot(n, geom, L, r2, rc, atom, modify_r2, deltaV)

        ! Declaration statements
        integer :: i, j
        integer, intent(in) :: n, atom
        real(kind=8), intent(in) :: L, rc
        real(kind=8), dimension(3,n), intent(in) :: geom
        real(kind=8), intent(inout) :: deltaV
        real(kind=8), dimension(n-1,n-1), intent(out) :: r2
        real(kind=8), dimension(n-1,n-1), intent(inout) :: modify_r2

        ! Execution zone
        deltaV = 0.d0
        modify_r2 = 0.d0
        !write(*,'( "Muevo atom", 1I5)') atom
        !write(*,*) " "

        ! The evaluation of deltaV is computed using two do loops, one for the modified terms 
        ! r2(atom-1,i), i.e. r(atom-1,1), r(atom-1,2)... until i = atom-1, and the other for the
        ! terms r2(i, atom) for i > atom, i.e. r2(atom+1,atom), r2(atom+2,atom)... until i = n.
        ! This can be done also with a do loop from 1 to n and an if-else inside the loop but it 
        ! is more computationally expensive.
        do i = 1, atom-1
            ! Periodic boundary conditions using the nearest image convention apply over the 
            ! modified coordinates

            ! Cartesian coordinate X
            if (abs(geom(1,atom)-geom(1,i)).gt.L/2.d0) then
                modify_r2(atom-1,i) = (geom(1,atom)-geom(1,i)-L)**2
            else
                modify_r2(atom-1,i) = (geom(1,atom)-geom(1,i))**2
            endif

            ! Cartesian coordinate Y
            if (abs(geom(2,atom)-geom(2,i)).gt.L/2.d0) then
                modify_r2(atom-1,i) = modify_r2(atom-1,i) + (geom(2,atom)-geom(2,i)-L)**2
            else
                modify_r2(atom-1,i) = modify_r2(atom-1,i) + (geom(2,atom)-geom(2,i))**2
            endif

            ! Cartesian coordinate Z
            if (abs(geom(3,atom)-geom(3,i)).gt.L/2.d0) then
                modify_r2(atom-1,i) = modify_r2(atom-1,i) + (geom(3,atom)-geom(3,i)-L)**2
            else
                modify_r2(atom-1,i) = modify_r2(atom-1,i) + (geom(3,atom)-geom(3,i))**2
            endif

            !write(*,*) "Modify r2, corresponding r2"
            !write(*,'(2f15.10)') modify_r2(atom-1,i), r2(atom-1,i)

            !if (modify_r2(i) .le. rc**2 .and. r2(atom-1,i) .le. rc**2) then
            !if (r2(atom-1,i) .le. rc**2) then
                !write(*,*) "********************"
                !write(*,*) "dentro cutoff"
                !write(*,*) "modify_r2, old r2"
                !write(*,'(2f10.4)') modify_r2(i), r2(atom-1,i)
                !write(*,*) " "
                !write(*,*) "DELTA V"
                !write(*,*) deltaV
            deltaV = deltaV + 1.d0/modify_r2(atom-1,i)**6 - 1.d0/modify_r2(atom-1,i)**3 &
            - 1.d0/r2(atom-1,i)**6 + 1.d0/r2(atom-1,i)**3
                !write(*,*) "DELTA V"
                !write(*,*) deltaV
                !write(*,*) " "
                !write(*,*) "deltaV"
                !write(*,*) deltaV
                !write(*,*) "********************"
            !else
                !write(*,*) "fuera cutoff"
            !endif
            !write(*,*) i, deltaV
        enddo

        do i = atom+1, n
            ! Periodic boundary conditions using the nearest image convention apply over the 
            ! modified coordinates

            ! Cartesian coordinate X
            if (abs(geom(1,i)-geom(1,atom)).gt.L/2.d0) then
               modify_r2(i-1,atom) = (geom(1,i)-geom(1,atom)-L)**2
            else
               modify_r2(i-1,atom) = (geom(1,i)-geom(1,atom))**2
            endif
            
            ! Cartesian coordinate Y
            if (abs(geom(2,i)-geom(2,atom)).gt.L/2.d0) then
               modify_r2(i-1,atom) = modify_r2(i-1,atom) + (geom(2,i)-geom(2,atom)-L)**2
            else
               modify_r2(i-1,atom) = modify_r2(i-1,atom) + (geom(2,i)-geom(2,atom))**2
            endif
            
            ! Cartesian coordinate Y
            if (abs(geom(3,i)-geom(3,atom)).gt.L/2.d0) then
               modify_r2(i-1,atom) = modify_r2(i-1,atom) + (geom(3,i)-geom(3,atom)-L)**2
            else
               modify_r2(i-1,atom) = modify_r2(i-1,atom) + (geom(3,i)-geom(3,atom))**2
            endif
            !modify_r2(i-1) = (geom(1,i)-geom(1,atom))**2 + (geom(2,i)-geom(2,atom))**2 &
            !    + (geom(3,i)-geom(3,atom))**2
            !write(*,*) "Modify r2, corresponding r2"
            !write(*,'(2f15.10)') modify_r2(i-1,atom), r2(i-1,atom)
            !write(*,*) " "
            !    write(*,*) "DELTA V"
            !    write(*,*) deltaV
            deltaV = deltaV + (1.d0/modify_r2(i-1,atom)**6-1.d0/modify_r2(i-1,atom)**3- & 
            1.d0/r2(i-1,atom)**6+1.d0/r2(i-1,atom)**3)
            !    write(*,*) "DELTA V"
            !    write(*,*) deltaV
            !    write(*,*) " "
        enddo

        !write(*,*) "Modified R2 matrix"
        !write(*,'(3f15.10)') modify_r2

        deltaV = 4.d0*deltaV

        return

    end subroutine delta_pot

    ! *********************************************************************************************



    ! *********************************************************************************************
    !                 SUBROUTINE 4: Simulation of the LJ fluid with MC tecniques
    ! *********************************************************************************************
    
    subroutine montecarlo(n, geom0, L, rc, V, V_rc, temp, maxcycle)

        ! Declaration statements
        integer :: i,j,k, counter, atom, max_hist
        integer, intent(in) :: n, maxcycle
        real(kind=8), dimension(3,n) :: geom0, geomi
        real(kind=8), dimension(n-1,n-1) :: r2
        real(kind=8), dimension(n-1,n-1) :: modify_r2
        real(kind=8), dimension(3) :: xyz
        real(kind=8), allocatable :: hist(:,:)
        real(kind=8) :: L, rc, V, V_rc, V_new, aux, delta_r, delta_V, a, temp, increment, rho, pi
        real(kind=8) :: deltaV, V1

        ! Execution zone

        ! The potential energy for the initial geometry is calculated
        call initial_V(n, geom0, L, r2, rc, V, V_rc)

        ! The mean density is computed
        rho = n/(L**3)
        ! Calculate pi
        pi = 4*atan(1.0)

        open(25, file="v_out", action="write")
        open(27, file="histogram", action="write")

        increment = 0.05
        max_hist = int(L/(2*increment))+1
        ! hist is the matrix which stored the values of delta_r, 2*delta_r... in its first 
        ! column and the times in which each interval is present over the n different structures
        allocate(hist(5,max_hist))

        ! Initialization of the counter of MC cycles
        counter = 0
        do while (counter .lt. maxcycle) 
            counter = counter + 1
            write(25,*) counter, V 
            ! A random atom is chosen and displaced a random quantity between -0.5 and 0.5 units
            call random_number(aux)
            atom = int(1+aux*n)
            call random_number(delta_r)
            delta_r = (delta_r-0.5)*(L/3.d0)
            !write(*,'("delta_r, 1f10.6")') delta_r
            ! The geometry in cycle i is geomi, i.e., the initial geometry with the random 
            ! displacement of atom i
            geomi = geom0
            geomi(:,atom) = geomi(:,atom)+delta_r
            ! Condition of restricted coordinates to a simulation box when atom is moved delta_r
            do i = 1, 3
                if (geomi(i,atom) .lt. 0.d0) then
                    geomi(i,atom) = geomi(i,atom) + L
                elseif (geomi(i,atom) .gt. L) then 
                    geomi(i,atom) = geomi(i,atom) - L
                endif
            enddo

            ! The energy change associated with the trial move is computed using the subroutine
            ! delta_pot
            call delta_pot(n, geomi, L, r2, rc, atom, modify_r2, deltaV)
            
            ! The program checks if the trial move is accepted or not. 
            if (deltaV .lt. 0.d0) then
                ! If the energy change associated with the trial move is negative it means that
                ! the total energy will be lower so the new configuration is more stable. Thus,
                ! the trial move is accepted and the values of geometry, potential and r^2 array
                ! are update
                geom0 = geomi
                V = V + deltaV
                do i = 1, atom-1
                    r2(atom-1,i) = modify_r2(atom-1,i)
                enddo
                do i = atom+1, n
                    r2(i-1,atom) = modify_r2(i-1,atom)
                enddo
            else
                ! Is the energy change associated with the trial move is positive, the procedure
                ! to decide is the trial move is accepted or not consists in computes the 
                ! Boltzmann factor (which will be close to 1 if deltaV is midly positive compared
                ! to kT and 0 if it is larger) and compared it with a random number (a) between 0
                ! and 1. If a is smaller than the Boltzmann factor the trial move is accepted,
                ! otherwise it is rejected
                call random_number(a)
                if (a .lt. exp(-deltaV)/temp) then
                    ! Trial move accepted
                    geom0 = geomi
                    V = V + deltaV
                    do i = 1, atom-1
                        r2(atom-1,i) = modify_r2(atom-1,i)
                    enddo
                    do i = atom+1, n
                        r2(i-1,atom) = modify_r2(i-1,atom)
                    enddo
                endif
                    ! Trial move rejected
            endif

            ! Pair correlation function g(r)
            ! Calculation of the number of structures present on each delta_r interval
            if (mod(counter,100) .eq. 0) then
                do k = 1, max_hist
                    do i = 2, n
                        do j = 1, i-1
                            if (sqrt(r2(i-1,j)) .gt. (k-1)*increment .and. sqrt(r2(i-1,j)) .lt.&
                            k*increment) then
                                ! Compute delta_r, 2*delta_r...
                                hist(1,k) = k*increment
                                ! Compute the number of structures on each interval
                                hist(2,k) = (hist(2,k) + 1)
                                ! Compute the volume of spherical shell
                                hist(3,k) = 4*pi*r2(i-1,j)*increment
                            endif
                        enddo
                    enddo
                    ! Compute the density at r
                    hist(4,k) = hist(2,k)/hist(3,k)
                    ! Compute g(r)
                    hist(5,k) = hist(4,k)/(rho*counter)
                enddo
            endif

        enddo
        write(27, '(5f15.3)') hist
        deallocate(hist)

        ! Final geometry configuration is written on final_geom.xyz file
        open(26, file="final_geom.xyz", action="write")
        write(26,'(I4)') n
        write(26,*) " "
        do i = 1, n
            xyz = geom0(:,i)
            write(26,'( "H", f10.6, f10.6, f10.6)' ) xyz(1), xyz(2), xyz(3)
        enddo
        write(*,*) "MC stop OK"

        return

    end subroutine montecarlo

    ! *********************************************************************************************

end module LJfluid
