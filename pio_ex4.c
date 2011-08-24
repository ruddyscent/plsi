/*
 * pio_ex4.c - Parallal I/O example source code that reads a file using a explicit offset and prints it to screen.
 *
 * Huioon Kim (pcandme@gist.ac.kr)
 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

void ErrorMessage(int error, int rank, char* string)
{
    fprintf(stderr, "Process %d: Error %d in %s\n", rank, error, string);
    MPI_Finalize();
    exit(-1);
}

int main(int argc, char* argv[])
{
    int comm_size, comm_rank, itest;
    int amode, oerr;

    MPI_Aint size_int;
    MPI_Status status;
    MPI_Offset offset;
    MPI_File fh;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    amode = MPI_MODE_RDONLY;

    MPI_Type_extent(MPI_INT, &size_int);

    oerr = MPI_File_open(MPI_COMM_WORLD, "data.dat", amode, MPI_INFO_NULL, &fh);

    if (oerr != MPI_SUCCESS) ErrorMessage(oerr, comm_rank, "MPI_File_open");

    offset = comm_rank * size_int;

    MPI_File_read_at(fh, offset, &itest, 1, MPI_INT, &status);

    printf("Hello from rank %d, I read: %d.\n", comm_rank, itest);
    
    MPI_File_close(&fh);

    MPI_Finalize();

    return 0;
}
