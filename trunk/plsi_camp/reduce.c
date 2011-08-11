#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int myrank, i, ista, iend;
    double a[9], sum, tmp;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    ista = myrank * 3;
    iend = ista + 2;

    for (i = ista; i < iend + 1; i++) a[i] = i;

    sum = 0.0;

    MPI_Reduce(&sum, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    sum = tmp;

    if (myrank == 0)
    {
        printf("sum = %f\n", sum);
    }

    MPI_Finalize();

    return 0;
}
