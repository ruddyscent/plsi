/*
 * Calculate the Cumulative Density Function(CDF).
 * compile command: mpicc -Wall -std=c99 -o cdf cdf.c
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double PI = 4.0 * atan(1.0);

// standard normal distribution function
inline double snd(double);

inline double snd(double x)
{
  double mu = 0; // mean
  double s2 = 1; // variance
  double exponent;

  exponent = -0.5 * pow((x - mu) / s2, 2);
  return exp(exponent) / sqrt(2 * s2 * PI);
}

int main(int argc, char *argv[])
{
  int n; /* # of rectangles */
  int myid, numprocs;
  const double PI = 4.0 * atan(1.0);
  double h, sum, x;
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
  
  startwtime = MPI_Wtime();

  for (int i = myid; i < n; i += numprocs) {
    x = h * i;
    printf("%f: %f\n", x, snd(x));
  }

  if (myid == 0) {
    endwtime = MPI_Wtime();
    printf("wall clock time = %f\n", endwtime - startwtime);
    fflush(stdout);
  }
  
  MPI_Finalize();

  return EXIT_SUCCESS;
}
