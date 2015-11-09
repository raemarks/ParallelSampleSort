all:
	mpicc -Wall -lm -o a.out main.c matrix.c random.c sample_sort.c

run:
	mpiexec -machinefile ./machinefile.txt -n 4 ./a.out

send:
	rsync machinefile1.txt machinefile2.txt machinefile4.txt run_tests.sh mpi_constants.h main.c matrix.c matrix.h random.h random.c sample_sort.h sample_sort.c Makefile rmarks@ssh1.eecs.wsu.edu:/net/u/rmarks/pvt/
