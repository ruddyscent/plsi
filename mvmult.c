/*
 * mvmult.c - Program for parallel matrix-vector multiplication using point-to-point communication.
 *
 * Huioon Kim (pcandme@gist.ac.kr)
 */

#include <stdio.h>
#include "mpi.h"

inline int min(int a, int b)
{
    return a < b ? a : b;
}

int main(int argc, char* argv[])
{
    const int MAX_ROWS = 1000;
    const int MAX_COLS = 1000; 

    int rows, cols;

    double a[MAX_ROWS][MAX_COLS], b[MAX_COLS], c[MAX_ROWS];
    double buffer[MAX_COLS], ans;

    int myid, master, numprocs, i, j, numsent, sender, anstype, row;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    master = 0;
    rows = 500;
    cols = 500;

    if (myid == master)
    {
        // master initializes and then dispatches...
        for (j = 0; j < cols; j++)
        {
            b[j] = 1;

            for (i = 0; i < rows; i++)
            {
                a[i][j] = i + 1;
            }
        }

        numsent = 0;

        MPI_Bcast(b, cols, MPI_DOUBLE, master, MPI_COMM_WORLD);

        for (i = 0; i < min(numprocs - 1, rows); i++)
        {
            for (j = 0; j < cols; j++)
            {
                buffer[j] = a[i][j];
            }

            MPI_Send(buffer, cols, MPI_DOUBLE, i + 1, i + 1, MPI_COMM_WORLD);

            numsent++;
        }

        for (i = 0; i < rows; i++)
        {
            MPI_Recv(&ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            sender = status.MPI_SOURCE;
            anstype = status.MPI_TAG;
            
            c[anstype - 1] = ans;

            if (numsent < rows)
            {
                for (j = 0; j < cols; j++)
                {
                    buffer[j] = a[numsent][j];
                }

                MPI_Send(buffer, cols, MPI_DOUBLE, sender, numsent + 1, MPI_COMM_WORLD);

                numsent++;
            }
            else
            {
                MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        // slaves receive b, then compute dot products until done message...
        MPI_Bcast(b, cols, MPI_DOUBLE, master, MPI_COMM_WORLD);

        if (myid <= rows)
        {
            do {
                MPI_Recv(&buffer, cols, MPI_DOUBLE, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                row = status.MPI_TAG;
                ans = 0.0;

                for (i = 0; i < cols; i++)
                {
                    ans = ans + buffer[i] * b[i];
                }

                MPI_Send(&ans, 1, MPI_DOUBLE, master, row, MPI_COMM_WORLD);
            } while (row != 0);
        }
    }

    if (myid == master)
    {
        printf("C = [");
        for (i = 0; i < rows; i++)
        {
            printf("%f", c[i]);

            if (i != rows - 1)
            {
                printf(", ");
            }
            else
                printf("]\n");
        }
    }

    MPI_Finalize();

    return 0;
}
