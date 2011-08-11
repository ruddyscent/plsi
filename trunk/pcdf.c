/*
 * Calculate the cumulative density function from the integration of the
 * standard normal distribution function.
 *
 * Kyungwon Chun (kwchun@gist.ac.kr)
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// probability density function
// mu: mean
// s2: variance
inline double pdf(double x, double mu, double s2)
{
  const double PI = acos(-1.0);
  const double exponent = -0.5 * pow((x - mu) / s2, 2);
  return exp(exponent) / sqrt(2 * PI * s2);
}

int main(int argc, char *argv[])
{
  const int n = 100000; // step size
  const double da = 0.1;
  double cdf, h, startwtime = 0.0, endwtime;
  int myid, numprocs, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &namelen);
  
  fprintf(stdout, "Process %d of %d is on %s\n",
	  myid, numprocs, processor_name);
  fflush(stdout);
  
  startwtime = MPI_Wtime();

  // See the 'Integrals over infinite intervals section of
  // http://en.wikipedia.org/wiki/Numerical_integration
  h = 1.0 / n;
  for (double a = 0; a <= 1; a += da) {
    double sum = 0;
    double mycdf;
    for (int j = myid; j < n; j += numprocs) {
      const double t = (j + 0.5) * h;
      sum += pdf(a - (1 - t) / t, 0, 1) / (t * t);
    }

    mycdf = h * sum;
    MPI_Reduce(&mycdf, &cdf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid == 0)
      printf("%f\t%f\n", a, cdf);
  }

  if (myid == 0) {
    endwtime = MPI_Wtime();
    printf("wall clock time = %f\n", endwtime - startwtime);
    fflush(stdout);
  }
  
  MPI_Finalize();

  return EXIT_SUCCESS;
}
