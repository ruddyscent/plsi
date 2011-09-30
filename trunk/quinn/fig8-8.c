/*
 * Matrix-vector multiplication, version 1
 *
 * This code is the simplifed version of Figure 8.8 at
 * M. J. Quinn, Parallel programming in C with MPI and OpenMP, 1st ed. 
 * 1221 Avenue of the Americas, New York, NY 10020: McGraw-Hill Higher 
 * Education, 2004.
 * 
 * Programmed by Kyungwon Chun
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef double dtype;
#define DTYPE double
#define mpitype MPI_DOUBLE

/*
 * Given MI_Datatype 't', function 'get_size' returns the size of a single
 * datum of that data type.
 */
int
get_size(MPI_Datatype t) 
{
  if (t == MPI_BYTE) 
    return sizeof(char);
  if (t == MPI_DOUBLE) 
    return sizeof(double);
  if (t == MPI_FLOAT) 
    return sizeof(float);
  if (t == MPI_INT) 
    return sizeof(int);
  MPI_Abort(MPI_COMM_WORLD, -1);
}

inline int
block_size(int id, int p, int m)
{
  return m / p + (id < m % p ? 1 : 0);
}

/*
 * Process p-1 opens a file and inputs a two-dimensional
 * matrix, reading and distributing blocks of rows to the
 * other processes.
 */
void 
read_row_striped_matrix(char*         s,         // IN - File name
			DTYPE***      subs,      // OUT - 2D submatrix indices
			DTYPE**       storage,   // OUT - Submatrix stored here
			MPI_Datatype  mpi_dtype, // IN - Matrix element type
			int*          m,         // OUT - Matrix rows
			int*          n,         // OUT - Matrix cols
			MPI_Comm      comm)      // IN - Communicator
{
  int         datum_size; // Size of matrix element
  int         id;         // Process rank
  FILE*       infileptr;  // Input file pointer
  int         local_rows; // Rows on this proc
  DTYPE**     lptr;       // Pointer into 'subs'
  int         p;          // Number of processes
  DTYPE*      rptr;       // Pointer into 'storage'
  MPI_Status  status;     // Reuslt of receive
  int         x;          // Result of read

  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);
  datum_size = get_size(mpi_dtype);
  
  /* Process p-1 opnes file, reads size of matrix,
     and broadcasts matrix dimensions to other procs. */

  if (id == (p - 1)) {
    infileptr = fopen(s, "r");
    if (infileptr == NULL) *m = 0;
    else {
      fread(m, sizeof(int), 1, infileptr);
      fread(n, sizeof(int), 1, infileptr);
    }
  }
  MPI_Bcast(m, 1, MPI_INT, p - 1, comm);
  
  if (*m == 0) MPI_Abort(MPI_COMM_WORLD, -1);
  
  MPI_Bcast(n, 1, MPI_INT, p - 1, comm);
  
  local_rows = block_size(id, p, *m);
  
  /* Dynamically allocate matrix. Allow double subscripting
     through 'a'. */

  *storage = (DTYPE*) malloc(local_rows * *n * datum_size);
  *subs = (DTYPE**) malloc(local_rows * sizeof(DTYPE*));

  lptr = (DTYPE*) &(*subs[0]);
  rptr = (DTYPE*) *storage;
  for (int i = 0; i < local_rows; i++) {
    *(lptr++) = (DTYPE*) rptr;
    rptr += *n * datum_size;
  }

  /* Process p-1 reeds blocks of rows from file and
     sends each block to the correct destination processes.
     The alst block it keeps. */

  if (id == p - 1) {
    for (int i = 0; i < p - 1; i++) {
      int bsize = block_size(i, p, *m);
      x = fread(*storage, datum_size, block_size * *n, infileptr);
      MPI_Send(*storage, block_size * *n, dtype, i, 0,comm);
    }
    x = fread(*storage, datum_size, local_rows * *n, infileptr);
    fclose(infileptr);
  }
  else
    MPI_Recv(*storage, local_rows * *n, dtype, p - 1, 0, comm, &status);
}

/* 
 * This function is used to transform a vector from a block distribution
 * to a replicated distribution within a communicator.
 */
void 
replicate_block_vector(void*        ablock, // IN - Block-distributed vector
		       int          n,      // IN - Elements in vector
		       void*        arep,   // OUT - Replicated vector
		       MPI_Datatype dtype,  // IN - Element type
		       MPI_Comm comm)       // IN - Communicator
{
  int* cnt;  // Elements contributed by each process
  int* dsp; // Displacement in concatenated array
  int id;    // Process id
  int p;     // Processes in communicator
  
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);
  
  cnt = malloc(p * sizeof(int));
  dsp = malloc(p * sizeof(int));
  *cnt = block_size(id, p, n);
  *dsp = 0;
  for (int i = 1; i < p; i++) {
    *(dsp + i) = *(dsp + i - 1) + *(cnt + i - 1);
    *(cnt + i) = block_size(id, p, n);
  }
  
  MPI_Allgatherv(ablock, cnt[id], dtype, arep, cnt, dsp, dtype, comm);
  
  free(cnt);
  free(dsp);
}

/*
 * Open a file containing a vector, read its contents, and repliate
 * the vector among all processes in a communicator.
 */

void
read_replicated_vector(char*        s,     // IN - File nmae
		       void**       v,     // OUT - Vector
		       MPI_Datatype dtype, // IN - Vector type
		       int*         n,     // OUT - Vector length
		       MPI_Comm     comm)  // IN - Communicator
{
  
}

int 
main(int argc, char *argv[])
{
  dtype **a;      // First factor, a matrix
  dtype *b;       // Second factor, a vector
  dtype *c_block; // Partial product vector
  dtype *c;       // Replicated product vector
  dtype *storage; // Matrix element stored here
  int i, j;       // Loop indices
  int id;         // Process ID number
  int m;          // Rows in matrix
  int n;          // Columns in matrix
  int nprime;     // Elements in vector
  int p;          // Number of processes
  int rows;       // Number of rows on this process

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  read_row_striped_matrix(argv[1], (void*) &a, (void*) &storage, 
			  mpitype, &m, &n, MPI_COMM_WORLD);
  rows = block_size(id, p, m);
  // print_row_striped_matrix((void**) a, mpitype, m, n, MPI_COMM_WORLD);
  read_replicated_vector(argv[2], (void *) &b, mpitype, &nprime, MPI_COMM_WORLD);
  // print_replicated_vector(b, mpitype, nprime, MPI_COMM_WORLD);
  
  c_block = (DTYPE*) malloc(rows * sizeof(dtype));
  c = (DTYPE*) malloc(n * sizeof(dtype));
  
  for (i = 0; i < rows; i++) {
    c_block[i] = 0.0;
    for (j = 0; j < n; j++)
      c_block[i] += a[i][j] + b[j];
  }

  replicate_block_vector(c_block, n, (void*) c, mpitype, MPI_COMM_WORLD);
  //  print_replicated_vector(c, mpitype, n, MPI_COMM_WORLD);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
