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

    ! Include openMP library
    use omp_lib
    implicit none
    public :: initial_geom, initial_V, delta_V, montecarlo
    contains
    ! List of common variables. Variables which are particular of each subroutine are described
    ! on the subroutine

    ! integers i, j and k in subroutines are counters for the do loops
    ! n = number of particles in the simulation box
    ! V = potential energy of the system
    ! rc = cutoff radius from which the interaction between two atoms is neglected
    ! Vrc = potential energy of the system at r = rc
    ! geom0 = 3 columns, n rows matrix which contains the initial set of (X, Y, Z) coordinates 
    ! geomi = 3 columns, n rows matrix which store the cartesian coordinates (X, Y, Z) of each
    ! particle in each row, it changes for each trial move.
    ! r2 = triangular matrix which store the value of the square of the distance between each
    ! pair of particles rij^2
    ! atom = particle over which a trial move is performed on each Monte Carlo cycle
    ! r2mod = 


    ! *********************************************************************************************
    !                          SUBROUTINE 1: Set the initial geometry
    ! *********************************************************************************************
    
    subroutine initial_geom(n,L,geom0)
    
        ! Declaration statements
        implicit none
        integer :: i
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

    subroutine initial_V(n, geom0, L, r2, rc, V, Vrc)

        ! Declaration statements
        implicit none
        integer :: i, j
        integer, intent(in) :: n
        real(kind=8), dimension(3,n), intent(in) :: geom0
        real(kind=8), intent(in) :: L, rc
        real(kind=8), dimension(n-1,n-1), intent(out) :: r2
        real(kind=8), intent(out) :: V, Vrc


        ! Execution zone
        ! Initially set V to 0
        V = 0.d0
        ! Compute the value of V at rc
        Vrc = 4*(1.d0/rc**12-1.d0/rc**6)

        do i = 2, n
            do j = 1, i-1
                ! Periodic boundary conditions using the nearest image convention in which
                ! we compute the interaction between a pair of atom which interatomic distance
                ! is less than |L/2|. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REVISAR ESTO

                ! Cartesian coordinate X
                if (abs(geom0(1,i)-geom0(1,j)).gt.L/2.d0) then
                    r2(i-1,j) = (geom0(1,i)-geom0(1,j)-L)**2
                else
                    r2(i-1,j) = (geom0(1,i)-geom0(1,j))**2
                endif

                ! Cartesian coordinate Y
                if (abs(geom0(2,i)-geom0(2,j)).gt.L/2.d0) then
                    r2(i-1,j) = r2(i-1,j) + (geom0(2,i)-geom0(2,j)-L)**2
                else
                    r2(i-1,j) = r2(i-1,j) + (geom0(2,i)-geom0(2,j))**2
                endif

                ! Cartesian coordinate Z
                if (abs(geom0(3,i)-geom0(3,j)).gt.L/2.d0) then
                    r2(i-1,j) = r2(i-1,j) + (geom0(3,i)-geom0(3,j)-L)**2
                else
                    r2(i-1,j) = r2(i-1,j) + (geom0(3,i)-geom0(3,j))**2
                endif

                ! Calculation of the potential using the cutoff rc
                if (r2(i-1,j)<rc**2) then
                    V = V + (1.d0/r2(i-1,j)**6-1.d0/r2(i-1,j)**3) - Vrc
                endif
            enddo
        enddo

        V = 4.d0*V

        return

    end subroutine initial_V

    ! *********************************************************************************************



    ! *********************************************************************************************
    !                 SUBROUTINE 3: Simulation of the LJ fluid with MC tecniques
    ! *********************************************************************************************

    subroutine delta_V(n, geomi, L, r2, atom, r2mod, dV, nthr)

        ! Declaration statements
        integer :: i
        integer, intent(in) :: n, atom, nthr
        real(kind=8), intent(in) :: L
        real(kind=8), dimension(3,n), intent(in) :: geomi
        real(kind=8), intent(inout) :: dV
        real(kind=8), dimension(n-1,n-1), intent(out) :: r2
        real(kind=8), dimension(n-1), intent(inout) :: r2mod

        ! Execution zone
        dV = 0.d0
        r2mod = 0.d0
        call omp_set_num_threads(nthr)

        ! The evaluation of dV is computed using two do loops, one for the modified terms 
        ! r2(atom-1,i), i.e. r(atom-1,1), r(atom-1,2)... until i = atom-1, and the other for the
        ! terms r2(i, atom) for i > atom, i.e. r2(atom+1,atom), r2(atom+2,atom)... until i = n.
        ! This can be done also with a do loop from 1 to n and an if-else inside the loop but it 
        ! is more computationally expensive.
        !$omp parallel do shared(r2mod, r2) private(i), reduction(+:dV)
        do i = 1, atom-1
            ! Periodic boundary conditions using the nearest image convention apply over the 
            ! modified coordinates

            ! Cartesian coordinate X
            if (abs(geomi(1,atom)-geomi(1,i)).gt.L/2.d0) then
                r2mod(i) = (geomi(1,atom)-geomi(1,i)-L)**2
            else
                r2mod(i) = (geomi(1,atom)-geomi(1,i))**2
            endif

            ! Cartesian coordinate Y
            if (abs(geomi(2,atom)-geomi(2,i)).gt.L/2.d0) then
                r2mod(i) = r2mod(i) + (geomi(2,atom)-geomi(2,i)-L)**2
            else
                r2mod(i) = r2mod(i) + (geomi(2,atom)-geomi(2,i))**2
            endif

            ! Cartesian coordinate Z
            if (abs(geomi(3,atom)-geomi(3,i)).gt.L/2.d0) then
                r2mod(i) = r2mod(i) + (geomi(3,atom)-geomi(3,i)-L)**2
            else
                r2mod(i) = r2mod(i) + (geomi(3,atom)-geomi(3,i))**2
            endif

            ! The value of dV is computed for the implied distances
            dV = dV + 1.d0/r2mod(i)**6 - 1.d0/r2mod(i)**3 &
            - 1.d0/r2(atom-1,i)**6 + 1.d0/r2(atom-1,i)**3
        enddo
        !$omp end parallel do

        !$omp parallel do shared(r2mod, r2) private(i), reduction(+:dV)
        do i = atom+1, n
            ! Periodic boundary conditions using the nearest image convention apply over the 
            ! modified coordinates

            ! Cartesian coordinate X
            if (abs(geomi(1,i)-geomi(1,atom)).gt.L/2.d0) then
               r2mod(i-1) = (geomi(1,i)-geomi(1,atom)-L)**2
            else
               r2mod(i-1) = (geomi(1,i)-geomi(1,atom))**2
            endif
            
            ! Cartesian coordinate Y
            if (abs(geomi(2,i)-geomi(2,atom)).gt.L/2.d0) then
               r2mod(i-1) = r2mod(i-1) + (geomi(2,i)-geomi(2,atom)-L)**2
            else
               r2mod(i-1) = r2mod(i-1) + (geomi(2,i)-geomi(2,atom))**2
            endif
            
            ! Cartesian coordinate Y
            if (abs(geomi(3,i)-geomi(3,atom)).gt.L/2.d0) then
               r2mod(i-1) = r2mod(i-1) + (geomi(3,i)-geomi(3,atom)-L)**2
            else
               r2mod(i-1) = r2mod(i-1) + (geomi(3,i)-geomi(3,atom))**2
            endif

            dV = dV + (1.d0/r2mod(i-1)**6-1.d0/r2mod(i-1)**3- & 
            1.d0/r2(i-1,atom)**6+1.d0/r2(i-1,atom)**3)
        enddo
        !$omp end parallel do

        dV = 4.d0*dV

        return

    end subroutine delta_V


    ! *********************************************************************************************



    ! *********************************************************************************************
    !                 SUBROUTINE 4: Simulation of the LJ fluid with MC tecniques
    ! *********************************************************************************************
    
    subroutine montecarlo(n, geom0, L, rc, V, Vrc, T, maxcycle, therm, nthr)

        ! Declaration statements
        implicit none
        ! rho = mean density of the fluid of particles N/V 
        ! pi = 3.1416...
        ! maxcycle = number of MC cycles, set by the user in the main.
        ! dr = random displacement of one particle (tagged with atom) in the trial move
        ! dV = energy change for the corresponding trial move dr
        ! T = temperature
        integer :: i, j, k, counter, atom, rmax, numMC, nthr, thrn
        integer, intent(in) :: n, maxcycle, therm
        real(kind=8), dimension(3,n) :: geom0, geomi
        real(kind=8), dimension(n-1,n-1) :: r2
        real(kind=8), dimension(n-1) :: r2mod
        real(kind=8), allocatable :: dat(:,:)
        real(kind=8) :: aux, L, rc, V, Vrc, dr, dV, a, T, increment, rho, pi, t0, tf, t0cpu, tfcpu


        ! Execution zone
        ! Number of threads and thread number 
        t0 = omp_get_wtime()
        call cpu_time(t0cpu)
        thrn = omp_get_thread_num()

        call omp_set_num_threads(nthr)
        !write(*,"( 'Number of procesors is ',I5)") omp_get_num_threads()

        ! The potential energy for the initial geometry is calculated
        call initial_V(n, geom0, L, r2, rc, V, Vrc)

        ! Calculate pi
        pi = 4*atan(1.0)
        ! The mean density is computed
        rho = n/(L**3)
        !rho = n/((4.d0/3.d0)*pi*(L/2.d0)**3)

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
        increment = 0.05
        rmax = int(L/(2*increment))+1

        ! The matrix dat is allocated. The required data will be stored in this matrix, i.e., the
        ! first column will store k*increment = increment, 2*increment..., the second column
        ! the number of structures on each spherical shell, the third column the volume of each 
        ! spherical shell and the last one the value of g(r).
        allocate(dat(4,rmax))
        ! The column which will store the volume is initialized with a number different from 
        ! zero to avoid zero in the denominators 
        dat(3,:) = 1.d0

        ! Initialization of the counter of MC cycles
        counter = 0

        ! Initialization of the number of snapshots that will be use to compute g(r)
        numMC = 0
        do while (counter .lt. maxcycle) 
            ! Increment the counter of MC cycles in a unit
            counter = counter + 1

            ! A random atom is chosen and displaced a random quantity between -0.5 and 0.5 units
            call random_number(aux)
            atom = int(1+aux*n)
            call random_number(dr)
            dr = (dr-0.5)*0.1
            ! The geometry in cycle i is geomi, i.e., the initial geometry with the random 
            ! displacement of atom i, dr
            geomi = geom0
            geomi(:,atom) = geomi(:,atom)+dr
            ! Apply the condition of restricted coordinates to the simulation box when the 
            ! atom is moved dr
            do i = 1, 3
                if (geomi(i,atom) .lt. 0.d0) then
                    geomi(i,atom) = geomi(i,atom) + L
                elseif (geomi(i,atom) .gt. L) then 
                    geomi(i,atom) = geomi(i,atom) - L
                endif
            enddo

            ! The energy change associated with the trial move is computed using the subroutine
            ! delta_V
            call delta_V(n, geomi, L, r2, atom, r2mod, dV, nthr)
            
            ! The program checks if the trial move is accepted or not
            if (dV .lt. 0.d0) then
                ! If the energy change associated with the trial move is negative it means that
                ! the total energy will be lower so the new configuration is more stable. Thus,
                ! the trial move is accepted and the values of geometry, potential and r^2 array
                ! are update
                geom0 = geomi
                V = V + dV
                !$omp parallel do shared(r2, r2mod) private(i)
                do i = 1, atom-1
                    r2(atom-1,i) = r2mod(i)
                    !write(*,"( 'Number of procesors is ',I5)") omp_get_num_threads()
                enddo
                !$omp end parallel do
                !$omp parallel do shared(r2, r2mod) private(i)
                do i = atom+1, n
                    r2(i-1,atom) = r2mod(i-1)
                enddo
                !$omp end parallel do
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
                    !$omp parallel do shared(r2, r2mod) private(i)
                    do i = 1, atom-1
                        r2(atom-1,i) = r2mod(i)
                    enddo
                    !$omp end parallel do
                    !$omp parallel do shared(r2, r2mod) private(i)
                    do i = atom+1, n
                        r2(i-1,atom) = r2mod(i-1)
                    enddo
                    !$omp end parallel do
                endif
                    ! Trial move rejected
            endif

            ! Write the potential energy in output file V.out
            if (mod(counter,100*n) .eq. 0) then
                write(25,*) counter, V 
            endif


            ! Pair correlation function g(r)
            ! The firsts MC steps are not considered (thermalization) set by the user in the main
            if (counter > therm .and. mod(counter,10*n) .eq. 0) then
                numMC = numMC + 1
                !$omp parallel do shared(r2, dat) private(i,j,k)
                do k = 1, rmax
                    ! Compute the array of interatomic distances that we will take into account
                    ! increment, 2*increment...
                    dat(1,k) = k*increment
                    do i = 2, n
                        do j = 1, i-1
                            if (sqrt(r2(i-1,j)) .gt. (k-1)*increment .and. & 
                                sqrt(r2(i-1,j)) .lt. k*increment) then
                                ! Compute the number of particles on each interval dN
                                ! N(r+increment)-N(r)
                                dat(2,k) = (dat(2,k) + 1)
                                ! Compute the volume of spherical shell dV 
                                ! V(r+increment)-V(r)
                                dat(3,k) = 4.d0*pi*r2(i-1,j)*increment
                            endif
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            endif

        enddo ! Monte Carlo end

        ! Calculate the mean number of structures between each (k-1)*increment, k*increment
        dat(2,:) = dat(2,:)/numMC
        ! Calculation of the pair correlation function g(r) 
        dat(4,:) = dat(2,:)/(dat(3,:)*rho)
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
        tf = omp_get_wtime()
        call cpu_time(tfcpu)
        print '("Time OMP = ", f10.5," seconds. ")', tf-t0
        print '("Time CPU = ", f10.5," seconds. ")', tfcpu-t0cpu

        return

    end subroutine montecarlo

    ! *********************************************************************************************

end module LJfluid
