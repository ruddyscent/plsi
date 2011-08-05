/*
 * Calculate the value of pi.
 * compile command: mpicc -Wall -std=c99 -o cpi cpi.c
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline double f(double);

inline double f(double a)
{
    return 4.0 / (1.0 + a*a);
}

int main(int argc, char *argv[])
{
    int n; /* # of rectangles */
    int myid, numprocs;
    const double PI = 4.0 * atan(1.0);
    double mypi, pi, h, sum, x;
    double startwtime = 0.0, endwtime;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
  
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);
  
    fprintf(stdout, "Process %d of %d is on %s\n",
	    myid, numprocs, processor_name);
    fflush(stdout);
  
    if (myid == 0) {
	if (argc == 2) 
	    sscanf(argv[1], "%d", &n);
	else {
	    fprintf(stderr, 
		    "Usage: cpi [# of rectangles]\n");
	    return EXIT_FAILURE;
	}
	startwtime = MPI_Wtime();
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
    h = 1.0 / (double)n;
    sum = 0.0;
    /* A slightly better approach starts from large i and works back. */
    for (int i = myid; i < n; i += numprocs) {
	x = h * ((double)i + 0.5);
	sum += f(x);
    }
    mypi = h * sum;
  
    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
    if (myid == 0) {
	endwtime = MPI_Wtime();
	printf("pi is approximately %.16f, Error is %.16f\n",
	       pi, fabs(pi - PI));
	printf("wall clock time = %f\n", endwtime - startwtime);
	fflush(stdout);
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
