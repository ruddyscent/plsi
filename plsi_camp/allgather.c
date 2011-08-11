#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int myrank, i;
    int isend, irecv[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    isend = myrank + 1;

    MPI_Allgather(&isend, 1, MPI_INT, irecv, 1, MPI_INT, MPI_COMM_WORLD);

    printf("irecv = ");
    for (i = 0; i < 3; i++)
    {
        printf("%d ", irecv[i]);
    }
    printf("\n");

    MPI_Finalize();

    return 0;
}
