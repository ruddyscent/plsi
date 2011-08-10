/* 
 * cpi2.c - C version of the Fortran90 code calculating pi in parallel, 
 * written in the USC parallel programming camp.
 *
 * Huioon Kim (pcandme@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline int min(int a, int b)
{
  return a < b ? a : b;
}

// Error bound estimation of a Riemann middle sum.
// Reference: http://en.wikipedia.org/wiki/Riemann_sum
double error_bound(int n)
{
  return 1.0 / (12.0 * n * n);
}

int main(int argc, char *argv[])
{
  const int NUM_STEP = 100000;
  int myrank, nprocs, ista, iend, remain, stride;
  double local_pi, pi, stime, etime;
  double sum = 0.0;
  double step = 1.0 / NUM_STEP;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  stride = NUM_STEP / nprocs;
  remain = NUM_STEP % nprocs;

  ista = myrank * stride + min(myrank, remain);
  iend = (myrank + 1) * stride + min(myrank + 1, remain);

  if (myrank == 0) {
    printf("----------------------------------------------\n");
    printf(" Error bound estimation: %11.5e\n", error_bound(NUM_STEP));
  }
  stime = MPI_Wtime();

  for (int i = ista; i < iend; i++) {
    double x = (i + 0.5) * step;
    sum += 4.0 / (1.0 + x * x);
  }

  local_pi = step * sum;

  MPI_Reduce(&local_pi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  etime = MPI_Wtime();

  if (myrank == 0) {
    printf(" Pi = %17.15f (Error = %11.5e)\n", pi, fabs(acos(-1.0) - pi));
    printf(" Elapsed Time = %8.4f [sec]\n", etime - stime);
    printf("----------------------------------------------\n");
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}
