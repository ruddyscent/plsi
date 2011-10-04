/*
 * Matrix-vector multiplication, version 1
 *
 * This code is the simplifed version of Figure 8.8 at
 * M. J. Quinn, Parallel programming in C with MPI and OpenMP, 
 * 1st ed. 1221 Avenue of the Americas, New York, NY 10020: 
 * McGraw-Hill Higher 
 * Education, 2004.
 * 
 * Usage: mpiexec -n <# of procs> fig8-8 <matrix data file> <vector data file>
 * 
 * Programmed by Kyungwon Chun
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef double dtype;
#define DTYPE double
#define mpitype MPI_DOUBLE

#define OPEN_FILE_ERROR -1
#define READ_FILE_ERROR -2
#define MALLOC_ERROR -3

inline int
block_size(int id, int p, int m)
{
  return m / p + (id < m % p ? 1 : 0);
}

void
read_row_striped_matrix(const char* filename, // IN - File name
			double** storage,     // OUT - Submatrix stored here
			int* m,               // OUT - Matrix rows
			int* n,               // OUT - Matrix columns
			MPI_Comm comm)        // IN - Communicator
{
  int id;                   // Process rank
  int p;                    // Number of processes
  FILE* infileptr = NULL;   // Input file pointer
  int local_rows;           // rows on this proc
  MPI_Status status;        // Result of receive

  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);
  
  /* Process p-1 opens file, reads size of matrix,
     and broadcasts matrix dimensions to other procs. */
  if (id == p - 1) {
    infileptr = fopen(filename, "r");
    if (infileptr == NULL) *m = 0;
    else {
      size_t result;
      result = fread(m, sizeof(int), 1, infileptr);
      if (result != 1) exit(READ_FILE_ERROR);
      result = fread(n, sizeof(int), 1, infileptr);
      if (result != 1) exit(READ_FILE_ERROR);
    }
  }
  MPI_Bcast(m, 1, MPI_INT, p - 1, comm);
  if (*m == 0) MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
  MPI_Bcast(n, 1, MPI_INT, p - 1, comm);
  
  local_rows = block_size(id, p, *m);
  
  /* Dynamically allocate matrix. Allow double subscripting 
     through 'a'. */
  *storage = (double *) malloc(local_rows * *n * sizeof(double));
  if (*storage == NULL) MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
  
  /* Process p-1 reads blocks of rows from file and sends
     each block to the correct destination process. The 
     last block it keeps. */
  if (id == p - 1) {
    for (int i = 0; i < p - 1; i++) {
      size_t result;
      int bsize;
      double* block;

      bsize = block_size(i, p, *m) * *n;
      block = (double*) malloc(bsize * sizeof(double));
      result = fread(block, sizeof(double), bsize, infileptr);
      if (result != bsize) exit(READ_FILE_ERROR);
      MPI_Send(block, bsize, MPI_DOUBLE, i, 0, comm);
      free(block);
    }
    size_t result;
    result = fread(*storage, sizeof(double), local_rows * *n, infileptr);
    fclose(infileptr);
  } else
    MPI_Recv(*storage, local_rows * *n, MPI_DOUBLE, p - 1, 0, comm, &status);
}

void
print_sub_matrix(const double* const storage, // IN - Matrix elements
		 int m,                       // IN - Matrix rows
		 int n)                       // IN - Matrix columns
{
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%6.3f ", *(storage + i * n + j));
    }
    putchar('\n');
  }
}

void
print_row_striped_matrix(double* storage, // IN - Matrix elements
			 int m,           // IN - Matrix rows
			 int n,           // IN - Matrix columns
			 MPI_Comm comm)   // IN - Communicator
{
  int id;            // Process rank
  int p;             // Number of processes
  int local_rows;    // This proc's lows
  MPI_Status status; // Result of receive

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  local_rows = block_size(id, p, m);
  if (id == 0) {
    print_sub_matrix(storage, local_rows, n);
    for (int i = 1; i < p; i++) {
      double* block;
      block = (double *) malloc(block_size(i, p, m) * n * sizeof(double));
      if (block == NULL) exit(MALLOC_ERROR);

      MPI_Recv(block, block_size(i, p, m) * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      print_sub_matrix(block, block_size(i, p, m), n);

      free(block);
    }
  }
  else
    MPI_Send(storage, local_rows * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

void
read_replicated_vector(const char* const s, // IN - File name
		       double** v,          // OUT - Vector
		       int* l,              // OUT - Vector length
		       MPI_Comm comm)       // IN - Communicator
{
  int id;                 // Process id
  int p;                  // Processes in communicator
  FILE* infileptr = NULL; // input file pointer

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  if (id == p - 1) {
    infileptr = fopen(s, "r");
    if (infileptr == NULL) *l = 0;
    else {
      size_t result;
      result = fread(l, sizeof(int), 1, infileptr);
      if (result != 1) exit(READ_FILE_ERROR);
    }
  }
  MPI_Bcast(l, 1, MPI_INT, p - 1, MPI_COMM_WORLD);
  if (*l == 0) MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
  
  *v = (double*) malloc(*l * sizeof(double));

  if (id == p - 1) {
    size_t result;
    result = fread(*v, sizeof(double), *l, infileptr);
    fclose(infileptr);
  }
  MPI_Bcast(*v, *l, MPI_DOUBLE, p - 1, MPI_COMM_WORLD);
}

void
print_replicated_vector(const double* const v, // IN - Vector pointer
	     int l,                 // IN - Length of vector
	     MPI_Comm comm)         // IN - Communicator
{
  int id; // Process rank

  MPI_Comm_rank(comm, &id);
  if (id == 0) {
    for (int i = 0; i < l; i++) {
      printf("%6.3f ", *(v + i));
    }
    putchar('\n');
  }
}

void
print_striped_vector(double* v,     // IN - vector pointer
		     int l,         // IN - length of vector
		     MPI_Comm comm) // IN - communicator
{
  int id; // process rank
  int p;  // number of processes
  int local_block;
  MPI_Status status; // result of receive

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  local_block = block_size(id, p, l);
  if (id == 0) {
    for (int i = 0; i < local_block; i++)
      printf("%6.3f ", *(v + i));
    for (int i = 1; i < p; i++) {
      int bsize;
      double* block;

      bsize = block_size(i, p, l);
      block = (double*) malloc(bsize * sizeof(double));
      MPI_Recv(block, bsize, MPI_DOUBLE, i, 0, comm, &status);
      for (int j = 0; j < bsize; j++)
	printf("%6.3f ", *(block + j));
      free(block);
    }
    putchar('\n');
  } else 
    MPI_Send(v, local_block, MPI_DOUBLE, 0, 0, comm);
}

int 
main(int argc, char *argv[])
{
  int id;          // Process ID number
  int p;           // Number of processes
  int l;           // Elements in vector
  int m;           // Rows in matrix
  int n;           // Columns in matrix
  double* a;       // First factor, a matrix
  double* b;       // Second factor, a vector
  double* c;       // Matrix elements stored here
  int rows;        // Number of rows on this process

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  read_row_striped_matrix(argv[1], &a, &m, &n, MPI_COMM_WORLD);
  print_row_striped_matrix(a, m, n, MPI_COMM_WORLD);
  if (id == 0) putchar('\n');
  read_replicated_vector(argv[2], &b, &l, MPI_COMM_WORLD);
  print_replicated_vector(b, l, MPI_COMM_WORLD);
  if (id == 0) putchar('\n');

  rows = block_size(id, p, m);
  c = (double*) calloc(rows, sizeof(double));
  
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < n; j++) {
      *(c + i) += *(a + i * n + j) * *(b + j);
    }
  }

  print_striped_vector(c, l, MPI_COMM_WORLD);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
