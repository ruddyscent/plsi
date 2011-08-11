#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int myrank, i;
    int imsg[4];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) for (i = 0; i < 4; i++) imsg[i] = i + 1;
    else for (i = 0; i < 4; i++) imsg[i] = 0;

    printf("%d: BEFORE:", myrank);
    for (i = 0; i < 4; i++) printf("%d", imsg[i]);
    printf("\n");

    MPI_Bcast(imsg, 4, MPI_INT, 0, MPI_COMM_WORLD);
    
    printf("%d: AFTER:", myrank);
    for (i = 0; i < 4; i++) printf("%d", imsg[i]);
    printf("\n");

    MPI_Finalize();

    return 0;
}
