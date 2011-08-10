/*
 * Calculate the pi.
 *
 * Kyungwon Chun (kwchun@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline double f(double a)
{
  return 4.0 / (1.0 + a * a);
}

// Error bound estimation of a Riemann middle sum.
// Reference: http://en.wikipedia.org/wiki/Riemann_sum
double error_bound(int n)
{
  return 1.0 / (12.0 * n * n);
}

int main(int argc, char *argv[])
{
  const double PI = 4.0 * atan(1.0);
  int n; /* # of rectangles */
  int myid, numprocs, namelen;
  double mypi, pi, h, sum = 0.0, startwtime = 0.0, endwtime;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &namelen);
  
  fprintf(stdout, "Process %d of %d is on %s\n",
	  myid, numprocs, processor_name);
  fflush(stdout);
  
  if (myid == 0) {
    if (argc == 2) {
	sscanf(argv[1], "%d", &n);
	printf("Error will be bounded in %e\n", error_bound(n));
    }
    else {
      fprintf(stderr, "Usage: cpi [# of rectangles]\n");
      return EXIT_FAILURE;
    }
    startwtime = MPI_Wtime();
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  h = 1.0 / n;
  /* A slightly better approach starts from large i and works back. */
  for (int i = myid; i < n; i += numprocs) {
    double x = h * (i + 0.5);
    sum += f(x);
  }
  mypi = h * sum;
  
  MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (myid == 0) {
    endwtime = MPI_Wtime();
    printf("The pi is approximately %.16f, Error is %e\n",
	   pi, fabs(pi - PI));
    printf("Wall clock time = %f [sec]\n", endwtime - startwtime);
    fflush(stdout);
  }

  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
