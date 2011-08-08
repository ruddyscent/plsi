/*
 * Calculate the pi using the Monte Carlo method.
 * compile command: mpicc -Wall -std=c99 -o pi_montecarlo pi_montecarlo.c
 *
 * Kyungwon Chun (kwchun@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[])
{
  const double PI = acos(-1.0);
  int n; /* # of trials */
  int *local_n;
  int count_inside = 0, my_count_inside = 0;
  double pi;
  int myid, numprocs;
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
      fprintf(stderr, "Usage: montecarlo [# of trials]\n");
      return EXIT_FAILURE;
    }
    startwtime = MPI_Wtime();
  }

  srand(time(NULL));

  local_n = (int *)malloc(numprocs * sizeof(int));
  if (myid == 0) {
    for (int i = 0; i < numprocs; i++) {
      local_n[i] = n / numprocs;
      if (i < n % numprocs) 
	local_n[i]++;
    }
  }
  
  MPI_Bcast(local_n, numprocs, MPI_INT, 0, MPI_COMM_WORLD);

  for (int i = 0; i < local_n[myid]; i++) {
    double x, y, d;
    x = (double)rand() / RAND_MAX;
    y = (double)rand() / RAND_MAX;
    d = sqrt(x*x + y*y);
    if (d < 1.0) 
      my_count_inside++;
  }
  
  MPI_Reduce(&my_count_inside, &count_inside, 1, MPI_INT, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  if (myid == 0) {
    pi = 4.0 * count_inside / n;
    endwtime = MPI_Wtime();
    printf("pi si approximately %.16f, Error is %.16f\n",
	   pi, fabs(pi - PI));
    printf("wall clock time = %f\n", endwtime - startwtime);
    fflush(stdout);
  }

  free(local_n);
  local_n = NULL;

  MPI_Finalize();

  return EXIT_SUCCESS;
}
