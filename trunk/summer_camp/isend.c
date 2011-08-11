#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int rank, i, count;
    float data[100], value[200];
    MPI_Request req;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 1)
    {
        for (i = 0; i < 100; ++i) data[i] = i;
        MPI_Isend(data, 100, MPI_FLOAT, 0, 55, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, &status);
    }
    else if (rank == 0)
    {
        for (i = 0; i < 100; ++i) value[i] = 0.0;

        MPI_Irecv(value, 200, MPI_FLOAT, 1, 55, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, &status);

        printf("P:%d Got data from processor %d\n", rank, status.MPI_SOURCE);
        MPI_Get_count(&status, MPI_FLOAT, &count);

        printf("P:%d Got %d elements\n", rank, count);
        printf("P:%d value[5] = %f\n", rank, value[5]);
    }

    MPI_Finalize();

    return 0;
}
