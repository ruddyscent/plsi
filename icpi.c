/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * This is an interactive version of cpi
 * 
 * Kyungwon Chun (kwchun@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline double f(double);
double error_bound(double);

inline double f(double a)
{
    return 4.0 / (1.0 + a * a);
}

// Error bound estimation of a Riemann middle sum.
// Reference: http://en.wikipedia.org/wiki/Riemann_sum
double error_bound(double n)
{
  return 2.0 * (1.0 - 0.0) / (24.0 * n * n);
}

int main(int argc, char *argv[])
{
    int n = 10000; /* default # of rectangles */
    int done = 0, myid, numprocs;
    const double PI = 4.0 * atan(1.0);
    double mypi, pi, h, sum, x;
    double startwtime = 0.0, endwtime;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    while (!done) {
	if (myid == 0) {
	    fprintf(stdout, "Enter the number of intervals(0 quits): ");
	    fflush(stdout);
	    if (scanf("%d", &n) == 1) {
		printf("Error will be bounded in %e\n", error_bound(n));
		startwtime = MPI_Wtime();
	    }
	    else {
		fprintf(stdout, "No number entered; quitting\n");
		n = 0;
	    }
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (n == 0)
	    done = 1;
	else {
	    h = 1.0 / (double)n;
	    sum = 0.0;
	    for (int i = myid; i < n; i += numprocs) {
		x = h * ((double)i + 0.5);
		sum += f(x);
	    }
	    mypi = h * sum;
	    
	    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    
	    if (myid == 0) {
		endwtime = MPI_Wtime();
		printf("pi is approximately %.16f, Error is %e\n",
		       pi, fabs(pi - PI));
		printf("wall clock time = %f\n", endwtime - startwtime);
		fflush(stdout);
	    }
	}
    }
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}
