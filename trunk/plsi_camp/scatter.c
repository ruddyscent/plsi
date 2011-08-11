#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int myrank, i, nprocs;
    int isend[3], irecv;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    for (i = 0; i < nprocs; i++) isend[i] = i + 1;

    MPI_Scatter(isend, 1, MPI_INT, &irecv, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("%d: irecv = %d\n", myrank, irecv);

    MPI_Finalize();

    return 0;
}
