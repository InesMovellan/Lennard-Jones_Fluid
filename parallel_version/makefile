# This makefile compiles (typing make) and executes (typing make exe) the main program with the 
# module which contains the Lennard-Jones fluid simulation

main: LJfluid.o main.o
	gfortran -fopenmp LJfluid.o main.o -o main.exe -fcheck=all -Wall

main.o: main.f90
	gfortran -fopenmp -c main.f90 -fcheck=all -Wall

LJfluid.o: LJfluid.f90
	gfortran -fopenmp -c LJfluid.f90 -fcheck=all -Wall

procs:
	export omp_num_procs=6

exe:
	./main.exe

clean:
	rm *.o *.mod
