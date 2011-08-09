/* 
 * cpi2.c - C version of the Fortran90 code calculating Pi in parallel written 
 * in the USC parallel programming camp.
 *
 * Huioon Kim (pcandme@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(a, b) (((a)<(b))? (a):(b))

#define NUM_STEP 100000

int main(int argc, char *argv[])
{
  int myrank, nprocs, i, ista, iend, remain, stride;
  double sum, step, local_pi, pi, x, stime, etime;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  stride = NUM_STEP / nprocs;
  remain = NUM_STEP % nprocs;

  ista = myrank * stride + min(myrank, remain) + 1;
  iend = ista + stride - 1;

  if (remain > myrank) 
    iend = iend + 1;

  step = (1.0 / NUM_STEP);
    
  sum = 0.0;

  if (myrank == 0)
    printf("----------------------------------------------\n");

  stime = MPI_Wtime();

  for (i = ista; i <= iend; i++) {
    x = (i - 0.5) * step;
    sum = sum + 4.0 / (1.0 + x * x);
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
