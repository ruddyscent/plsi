/*
 * Calculate the cumulative density function(CDF) from the integration of the
 * standard normal distribution(SND) function.
 *
 * Kyungwon Chun (kwchun@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI acos(-1.0)

// standard normal distribution function
inline double snd(double);

inline double snd(double x)
{
  double mu = 0; // mean
  double s2 = 1; // variance
  double exponent;

  exponent = -0.5 * pow((x - mu) / s2, 2);
  return exp(exponent) / sqrt(2 * PI * s2);
}

int main(int argc, char *argv[])
{
  double h = 1.e-5; // step size
  double left = -15.0; // left boundary of integration
  double cdf;
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
  
  startwtime = MPI_Wtime();

  for (double right = 0; right <= 1.0; right += 0.1) {
    int i = myid;
    double sum = 0.0;
    double x = (i + 0.5) * h + left;
    double mycdf;

    while (x < right) {
      x = (i + 0.5) * h + left;
      sum += snd(x);
      i += numprocs;
    }
    mycdf = h * sum;
    MPI_Reduce(&mycdf, &cdf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid == 0)
      printf("%f\t%f\n", right, cdf);
  }

  if (myid == 0) {
    endwtime = MPI_Wtime();
    printf("wall clock time = %f\n", endwtime - startwtime);
    fflush(stdout);
  }
  
  MPI_Finalize();

  return EXIT_SUCCESS;
}
