/*
 * pio_ex1.c - Parallal I/O example source code that writes a file with individual file pointers.
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
    int comm_size, comm_rank;
    int amode, etype, filetype, oerr;

    MPI_Aint size_int;
    MPI_Status status;
    MPI_Offset disp;
    MPI_File fh;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;

    MPI_Type_extent(MPI_INT, &size_int);

    oerr = MPI_File_open(MPI_COMM_WORLD, "data.dat", amode, MPI_INFO_NULL, &fh);

    if (oerr != MPI_SUCCESS) ErrorMessage(oerr, comm_rank, "MPI_File_open");

    disp = comm_rank * size_int;
    etype = MPI_INT;
    filetype = MPI_INT;

    MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);

    MPI_File_write(fh, &comm_rank, 1, MPI_INT, &status);

    printf("Hello from rank %d, I wrote: %d.\n", comm_rank, comm_rank);

    MPI_File_close(&fh);

    MPI_Finalize();

    return 0;
}
