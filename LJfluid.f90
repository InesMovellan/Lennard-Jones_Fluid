! Module which carried out the Lennard-Jones simulation in a periodic box using Monte Carlo 
! algorithm. It contains several subroutines:
! Serial version of the module.

! SUBROUTINE 1: initial_geom
!               Sets the initial geometry for n Lennard-Jones particles in a cubic box, which 
!               size (and also n) is chosen by the user in main.f90 

! SUBROUTINE 2: pbc
!               Apply periodic boundary conditions over a pair of atoms i and j received as
!               arguments, following the nearest image convention

! SUBROUTINE 3: initial_V
!               Computes the initial potential energy associated with the initial geometry
!               configuration of the particles obtained with initial_geom subroutine. Since the
!               model of the system is a cubic box surrounded by an infinite number of replicas, 
!               the calculation of the energy is done with periodic boundary conditions (PBC).

! SUBROUTINE 4: delta_V
!               Calculates the energy change between one trial move performed in montecarlo
!               subroutine and the previous geometry configuration (also with PBC).

! SUBROUTINE 5: montecarlo
!               Monte Carlo subroutine which performs the sampling of structures over the number
!               of cycles selected by the user in the main.f90. The pair correlation function 
!               g(r) is also computed in this subroutine.

! To use this module main.f90 is needed, also contained in this directory

! Author: Inés Sánchez de Movellán Sáiz

module LJfluid

    implicit none
    public :: initial_geom, pbc, initial_V, delta_V, montecarlo
    contains
    ! List of common variables. Variables which are particular of each subroutine are described
    ! on the subroutine

    ! integers i, j and k in subroutines are counters for the do loops
    ! n = number of particles in the simulation box
    ! L = box side length
    ! dx, dy, dz = difference between x, y and z coordinates of two particles when periodic 
    !              boundary conditions are applied
    ! V = potential energy of the system
    ! dV = energy change when a trial move is performed
    ! rc = cutoff radius from which the interaction between two atoms is neglected
    ! Vrc = potential energy of the system at r = rc
    ! geom0 = 3 columns, n rows matrix which contains the initial set of (X, Y, Z) coordinates 
    ! geomi = 3 columns, n rows matrix which store the cartesian coordinates (X, Y, Z) of each
    !         particle in each row, it changes for each trial move.
    ! r2 = triangular matrix which store the value of the square of the distance between each
    !      pair of particles rij^2
    ! atom = particle over which a trial move is performed on each Monte Carlo cycle
    ! r2mod = array of length n-1 which stores the modified value of the square of the interatomic
    !         distance between a atoms over which a trial move is performed and the rest i.e. 
    !         pair atom,i


    ! *********************************************************************************************
    !                          SUBROUTINE 1: Set the initial geometry
    ! *********************************************************************************************
    
    subroutine initial_geom(n,L,geom0)
    
        ! Declaration statements
        ! r = interatomic distance
        implicit none
        integer :: i, j, cont
        integer, intent(in) :: n  
        real(kind=8), intent(in) :: L  
        real(kind=8), allocatable, intent(inout) :: geom0(:,:)
        real(kind=8) :: r, dx, dy, dz

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
        ! (0, L) range, checking that the particles are not too much closer, i.e., interatomic
        ! distances are larger than 1 (in reduced units). The geometry is set taking into 
        ! account periodic boundary conditions
        geom0 = 0.d0
        call random_number(geom0(:,1))
        geom0 = L*geom0
        do i = 2, n
            r = 0.d0
            do j = 1, i-1
                cont = 0
                do while (r .lt. 1.d0)
                    call random_number(geom0(:,i))
                    geom0(:,i) = L*geom0(:,i)
                    call pbc(n, geom0, L, i, j, dx, dy, dz)
                    r = sqrt(dx**2 + dy**2 + dz**2)
                enddo
            enddo
        enddo

        ! The initial geometry is written in intial_geom.xyz file (following the format of .xyz
        ! files), which can be read with 3D visualization programs.
        open(23, file="initial_geom.xyz", action="write")
        write(23,'(I4)') n
        write(23,*) " "
        do i = 1, n
            write(23,'( "H", f10.6, f10.6, f10.6)' ) geom0(1,i), geom0(2,i), geom0(3,i) 
        enddo
        write(*,*) "Initial geometry defined"
        return

    end subroutine initial_geom

    ! *********************************************************************************************



    ! *********************************************************************************************
    !           SUBROUTINE 2: Apply periodic boundary conditions over determinated geometry
    ! *********************************************************************************************

    subroutine pbc(n, geom, L, i, j, dx, dy, dz)

        ! Declaration statements
        ! geom = matrix which stores the cartesian coordinates of each particle
        implicit none
        integer, intent(in) :: n, i, j
        real(kind=8), dimension(3,n), intent(in) :: geom
        real(kind=8), intent(in) :: L
        real(kind=8), intent(inout) :: dx, dy, dz


        ! Execution zone
        ! Periodic boundary conditions using the nearest image convention in which
        ! we compute the square of the distance of a pair of atoms (i and j) which interatomic
        ! distance is lower than L/2.

        ! Cartesian coordinate X
        if (geom(1,i)-geom(1,j).gt.L/2.d0) then
            dx = geom(1,i)-geom(1,j)-L
        else if (geom(1,i)-geom(1,j).lt.(-L/2.d0)) then
            dx = geom(1,i)-geom(1,j)+L
        else
            dx = geom(1,i)-geom(1,j)
        endif

        ! Cartesian coordinate Y
        if (geom(2,i)-geom(2,j).gt.L/2.d0) then
            dy = geom(2,i)-geom(2,j)-L
        else if (geom(2,i)-geom(2,j).lt.(-L/2.d0)) then
            dy = geom(2,i)-geom(2,j)+L
        else
            dy = geom(2,i)-geom(2,j)
        endif

        ! Cartesian coordinate Z
        if (geom(3,i)-geom(3,j).gt.L/2.d0) then
            dz = geom(3,i)-geom(3,j)-L
        else if (geom(3,i)-geom(3,j).lt.(-L/2.d0)) then
            dz = geom(3,i)-geom(3,j)+L
        else
            dz = geom(3,i)-geom(3,j)
        endif

        return

    end subroutine pbc

    ! *********************************************************************************************


    ! *********************************************************************************************
    !              SUBROUTINE 3: Calculate the initial potential energy with PBC
    ! *********************************************************************************************

    subroutine initial_V(n, geom0, L, r2, rc, V, Vrc)

        ! Declaration statements
        implicit none
        integer :: i, j
        integer, intent(in) :: n
        real(kind=8), dimension(3,n), intent(in) :: geom0
        real(kind=8), intent(in) :: L, rc
        real(kind=8), dimension(n-1,n-1), intent(out) :: r2
        real(kind=8), intent(out) :: V, Vrc
        real(kind=8) :: dx, dy, dz


        ! Execution zone
        ! Initially set V to 0
        V = 0.d0
        ! Compute the value of V at rc
        Vrc = 4*(1.d0/rc**12-1.d0/rc**6)

        do i = 2, n
            do j = 1, i-1
                ! Call subroutine pbc to apply periodic boundary conditions
                call pbc(n, geom0, L, i, j, dx, dy, dz)

                ! Calculate the square of interatomic distance between the pair i and j
                r2(i-1,j) = dx**2 + dy**2 + dz**2

                ! Calculation of the potential using the cutoff rc
                if (r2(i-1,j)<rc**2) then
                    V = V + (1.d0/r2(i-1,j)**6-1.d0/r2(i-1,j)**3) - Vrc
                endif
            enddo
         enddo

        V = 4.d0*V
        write(*,*) "Initial potential energy calculated"

        return

    end subroutine initial_V

    ! *********************************************************************************************



    ! *********************************************************************************************
    !    SUBROUTINE 4: Compute the change in the potential energy for a trial move of one atom
    ! *********************************************************************************************

    subroutine delta_V(n, geomi, L, r2, atom, r2mod, dV)

        ! Declaration statements
        integer :: i
        integer, intent(in) :: n, atom
        real(kind=8), intent(in) :: L
        real(kind=8), dimension(3,n), intent(in) :: geomi
        real(kind=8), intent(inout) :: dV
        real(kind=8), dimension(n-1,n-1), intent(out) :: r2
        real(kind=8), dimension(n-1), intent(inout) :: r2mod
        real(kind=8) :: dx, dy, dz

        ! Execution zone
        dV = 0.d0
        r2mod = 0.d0

        ! The evaluation of dV is computed using two do loops, one for the modified terms 
        ! r2(atom-1,i), i.e. r(atom-1,1), r(atom-1,2)... until i = atom-1, and the other for the
        ! terms r2(i, atom) for i > atom, i.e. r2(atom+1,atom), r2(atom+2,atom)... until i = n.
        ! This can be done also with a do loop from 1 to n and an if-else inside the loop but it 
        ! is more computationally expensive.
        do i = 1, atom-1
            call pbc(n, geomi, L, atom, i, dx, dy, dz)
            r2mod(i) = dx**2 + dy**2 + dz**2
            ! The value of dV is computed for the implied distances
            dV = dV + 1.d0/r2mod(i)**6 - 1.d0/r2mod(i)**3 &
            - 1.d0/r2(atom-1,i)**6 + 1.d0/r2(atom-1,i)**3
        enddo

        do i = atom+1, n
            call pbc(n, geomi, L, i, atom, dx, dy, dz)
            r2mod(i-1) = dx**2 + dy**2 + dz**2
            ! The value of dV is computed for the implied distances
            dV = dV + (1.d0/r2mod(i-1)**6-1.d0/r2mod(i-1)**3- & 
            1.d0/r2(i-1,atom)**6+1.d0/r2(i-1,atom)**3)
        enddo

        dV = 4.d0*dV

        return

    end subroutine delta_V


    ! *********************************************************************************************



    ! *********************************************************************************************
    !                 SUBROUTINE 5: Simulation of the LJ fluid with MC tecniques
    ! *********************************************************************************************
    
    subroutine montecarlo(n, geom0, L, rc, V, Vrc, T, maxcycle, therm)

        ! Declaration statements
        implicit none
        ! counter = Monte Carlo step
        ! kmax = total number of r+dr intervals 
        ! samplMC = number of Monte Carlo sweeps used to calculate g(r)
        ! maxcycle = number of MC cycles, set by the user in the main.
        ! therm = MC cycles used for thermalization
        ! dr = random displacement of one particle (tagged with atom) in the trial move
        ! dat = 4 column, kmax rows matrix which stores each value of r+dr in the first column,
        !       the number of particles in each spherical shell in the second, the volume of 
        !       each spherical shell in the third and g(r) in the fourth.
        ! aux = auxiliar variable used to calculate the random atom
        ! a = random number used to accept or not a trial move when dV > 0
        ! T = temperature
        ! increment = dr used in the computation of g(r) 
        ! rho = mean density of the fluid of particles N/V 
        ! pi = 3.1416...
        ! t0 and tf = initial and final cpu_time
        integer :: i, j, k, counter, atom, kmax, samplMC
        integer, intent(in) :: n, maxcycle, therm
        real(kind=8), dimension(3,n) :: geom0, geomi
        real(kind=8), dimension(n-1,n-1) :: r2
        real(kind=8), dimension(n-1) :: r2mod
        real(kind=8), dimension(3) :: dr
        real(kind=8), allocatable :: dat(:,:)
        real(kind=8) :: aux, L, rc, V, Vrc, dV, a, T, increment, rho, pi, t0, tf

        ! Execution zone
        call cpu_time(t0)

        ! The potential energy for the initial geometry is calculated
        call initial_V(n, geom0, L, r2, rc, V, Vrc)

        ! Calculate pi
        pi = 4*atan(1.0)
        ! The mean density is computed
        rho = n/(L**3)

        open(25, file="V.out", action="write")
        write(25,*) "Potential energy of the fluid of Lennard Jones particles"
        write(25,*) " "
        write(25,*) "Column 1: Monte Carlo step"
        write(25,*) "Column 2: Potential energy in reduced units"
        write(25,*) " "
        open(26, file="g.out", action="write")
        write(26,*) "Pair correlation function of the fluid of Lennard Jones particles"
        write(26,*) " "
        write(26,*) "Column 1: interatomic distance r in reduced units (r/sigma)"
        write(26,*) "Column 2: number of particles between each interval r+dr (dN)"
        write(26,*) "Column 3: volume of each spherical shell dV in reduced units"
        write(26,*) "Column 4: pair correlation function g(r) in reduced units"
        write(26,*) " "

        ! Definition of increment and maximun r for the calculation of the pair correlation
        ! function. The maximun value of r is the half of simulation box length (L/2) divided
        ! by the increment, i.e., the number of intervals of size = increment (increment, 
        ! 2*increment, 3*increment,..., L/2)
        increment = 0.04
        kmax = int(L/(2*increment))+1

        ! The matrix dat is allocated. 
        allocate(dat(4,kmax))
        ! The column which will store the volume is initialized with a number different from 
        ! zero to avoid zero in the denominators 
        dat(3,:) = 1.d0
        ! Compute the array of interatomic distances that we will take into account
        ! increment, 2*increment...
        do k = 1, kmax
            dat(1,k) = k*increment
        enddo

        ! Initialization of the counter of MC cycles
        counter = 0
        ! Initialization of the number of MC sweeps that will be use to compute g(r)
        samplMC = 0

        write(*,*) "Monte Carlo simulation starts"
        do while (counter .lt. maxcycle) 
            ! Increment the counter of MC cycles in a unit
            counter = counter + 1

            ! A random atom is chosen and displaced a random quantity between -0.02 and 0.02 units
            call random_number(aux)
            atom = int(1+aux*n)
            call random_number(dr)
            dr = (dr-0.5)*0.04
            ! The geometry in cycle i is geomi, i.e., the initial geometry with the random 
            ! displacement of atom i, dr
            geomi = geom0
            geomi(:,atom) = geomi(:,atom)+dr
            ! Apply the condition of restricted coordinates to the simulation box 
            do i = 1, 3
                if (geomi(i,atom) .lt. 0.d0) then
                    geomi(i,atom) = geomi(i,atom) + L
                elseif (geomi(i,atom) .gt. L) then 
                    geomi(i,atom) = geomi(i,atom) - L
                endif
            enddo

            ! The energy change associated with the trial move is computed using the subroutine
            ! delta_V
            call delta_V(n, geomi, L, r2, atom, r2mod, dV)
            
            ! The program checks if the trial move is accepted or not
            if (dV .lt. 0.d0) then
                ! If the energy change associated with the trial move is negative it means that
                ! the total energy will be lower so the new configuration is more stable. Thus,
                ! the trial move is accepted and the values of geometry, potential and r^2 array
                ! are update
                geom0 = geomi
                V = V + dV
                do i = 1, atom-1
                    r2(atom-1,i) = r2mod(i)
                enddo
                do i = atom+1, n
                    r2(i-1,atom) = r2mod(i-1)
                enddo
            else
                ! Is the energy change associated with the trial move is positive, the procedure
                ! to decide is the trial move is accepted or not consists in computes the 
                ! Boltzmann factor (which will be close to 1 if dV is midly positive compared
                ! to kT and 0 if it is larger) and compared it with a random number (a) between 0
                ! and 1. If a is smaller than the Boltzmann factor the trial move is accepted,
                ! otherwise it is rejected
                call random_number(a)
                if (a .lt. exp(-dV)/T) then
                    ! Trial move accepted
                    geom0 = geomi
                    V = V + dV
                    do i = 1, atom-1
                        r2(atom-1,i) = r2mod(i)
                    enddo
                    do i = atom+1, n
                        r2(i-1,atom) = r2mod(i-1)
                    enddo
                endif
                    ! Trial move rejected
            endif

            if (mod(counter,100*n) .eq. 0) then
                write(25,*) counter, V 
            endif

            ! Pair correlation function g(r)
            ! The first MC steps are not considered (thermalization) set by the user in the main
            if (counter > therm .and. mod(counter,10*n) .eq. 0) then
                samplMC = samplMC + 1
                do i = 2, n
                    do j = 1, i-1
                        k = int((sqrt(r2(i-1,j))/increment)+1)
                        if (k .le. kmax) then
                            ! Compute the number of particles on each interval dN 
                            ! N(r+increment)-N(r)
                            dat(2,k) = dat(2,k) + 1
                            ! Compute the volume of spherical shell dV 
                            ! V(r+increment)-V(r)
                            dat(3,k) = 4.d0*pi*r2(i-1,j)*increment
                        endif
                    enddo
                enddo
            endif

        enddo ! Monte Carlo end
        write(*,*) "Monte Carlo simulation finishes"

        ! Calculate the mean number of structures between each (k-1)*increment, k*increment
        dat(2,:) = dat(2,:)/samplMC
        ! Calculation of the pair correlation function g(r) 
        dat(4,:) = dat(2,:)/(dat(3,:)*rho*n/2.d0)
        ! Write the pair correlation function and associated quantities in output file g.out
        write(26, '(4f17.10)') dat
        deallocate(dat)

        ! Final geometry configuration is written on final_geom.xyz file
        open(27, file="final_geom.xyz", action="write")
        write(27,'(I4)') n
        write(27,*) " "
        do i = 1, n
            write(27,'( "H", f10.6, f10.6, f10.6)' ) geom0(1,i), geom0(2,i), geom0(3,i)
        enddo
        call cpu_time(tf)
        print '("Time = ",f10.5," seconds")', tf-t0 

        return

    end subroutine montecarlo

    ! *********************************************************************************************

end module LJfluid
