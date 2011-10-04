/*
 * Matrix-vector multiplication, version 2
 *
 * This code is the simplifed version of Figure 8.8 at
 * M. J. Quinn, Parallel programming in C with MPI and OpenMP, 
 * 1st ed. 1221 Avenue of the Americas, New York, NY 10020: 
 * McGraw-Hill Higher 
 * Education, 2004.
 * 
 * Usage: mpiexec -n <# of procs> fig8-14 <matrix data file> <vector data file>
 * 
 * Programmed by Kyungwon Chun
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define OPEN_FILE_ERROR -1
#define READ_FILE_ERROR -2
#define MALLOC_ERROR -3

inline int
block_size(int id, // IN - process rank
	   int p,  // IN - # of processes
	   int m)  // IN - # of elements
{
  return m / p + (id < m % p ? 1 : 0);
}

/* This function creates the count and displacement arrays needed 
   by scatter and gather functions, when the number of elements
   send/received to/from othe rprocesses varies. */
void
create_mixed_xfer_arrays(int n,  // IN - total # of elements
			 int** cnt, // OUT - array of counts
			 int** dsp, // OUT - array of displacements
			 MPI_Comm comm) // IN - communicator
{
  int p; // # of procs
  int id; // this proc's rank
  
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);
  
  *cnt = (int *) malloc(p * sizeof(int));
  *dsp = (int *) malloc(p * sizeof(int));
  **cnt = block_size(0, p, n);
  **dsp = 0;

  for (int i = 1; i < p; i++) {
    (*dsp)[i] = (*dsp)[i-1] + (*cnt)[i-1];
    (*cnt)[i] = block_size(i, p, n);
  }
}

/*
 * This funciton reads a matrix from a file and allocates block of
 * columns of the matrix to the MPI processes. The first two 
 * elements of the file are integers whose values are the 
 * dimensions of the matrix ('m' rows and 'n' columns). What 
 * follows are 'm'*'n' values representing the matrix elements
 * stored in row-major order. 
 */
void 
read_col_striped_matrix(const char* const s, // IN - file name
			double** mat, // OUT - matrix elements
			int* m, // OUT - # of rows
			int* n, // OUT - # of columns
			MPI_Comm comm) // IN - communicator
{
  int id;                  // process rank
  int p;                   // # of processes
  FILE* infile_ptr = NULL; // input file pointer
  int local_cols;          // # of cols here
  int *send_cnt;           // each proc's count
  int *send_dsp;           // each proc's displacement
  double* buffer = NULL;   // flie buffer

  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);

  // Process p-1 opens file, gets # of rows and cols, and
  // broadcasts this info to other procs.
  
  if (id == p - 1) {
    infile_ptr = fopen(s, "r");
    if (infile_ptr == NULL) {
      *m = 0;
      *n = 0;
    }
    else {
      size_t result;
      result = fread(m, sizeof(int), 1, infile_ptr);
      if (result != 1) *m = -1;
      result = fread(n, sizeof(int), 1, infile_ptr);
      if (result != 1) *n = -1;
    }
  }
  MPI_Bcast(m, 1, MPI_INT, p - 1, comm);
  MPI_Bcast(n, 1, MPI_INT, p - 1, comm);
  if (*m == 0 || *n == 0) MPI_Abort(comm, OPEN_FILE_ERROR);
  if (*m == -1 || *n == -1) MPI_Abort(comm, READ_FILE_ERROR);
  
  // Process p-1 reads in the matrix one row at a time and
  // distributes each row among the MPI processes.
  local_cols = block_size(id, p, *n);

  *mat = malloc(*m * local_cols * sizeof(double));

  if (id == p - 1)
    buffer = malloc(*n * sizeof(double));

  create_mixed_xfer_arrays(*n, &send_cnt, &send_dsp, comm);
  for (int i = 0; i < *m; i++) {
    if (id == p - 1) {
      size_t result;
      result = fread(buffer, sizeof(double), *n, infile_ptr);
    }
    MPI_Scatterv(buffer, send_cnt, send_dsp, MPI_DOUBLE,
		 *mat + i * local_cols, local_cols,
  		 MPI_DOUBLE, p - 1, comm);
  }
  free(send_cnt);
  free(send_dsp);
  if (id == p - 1) free(buffer);
}

void
print_col_striped_matrix(double* a,     // IN - matrix
			 int m,         // IN - # of rows
			 int n,         // IN - # of cols
			 MPI_Comm comm) // IN - communicator
{
  int id; // proc's rank
  int p;  // # of procs
  int* rec_cnt; // # of elements received per proc
  int* rec_dsp; // # offset of eac hproc's block
  double* buffer = NULL; // room to hold 1 row
  int local_cols; // # of cols in this proc

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);
  local_cols = block_size(id, p, n);
  create_mixed_xfer_arrays(n, &rec_cnt, &rec_dsp, comm);
  
  if (id == 0) buffer = malloc(n * sizeof(double));
  
  for (int i = 0; i < m; i++) {
    MPI_Gatherv(a + i * local_cols, local_cols, MPI_DOUBLE, buffer,
		rec_cnt, rec_dsp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (id == 0) {
      for (int j = 0; j < n; j++) {
	printf("%6.3f ", buffer[j]);
      }
      putchar('\n');
    }
  }
  free(rec_cnt);
  free(rec_dsp);
  if (id == 0) {
    free(buffer);
    putchar('\n');
  }
}

void
read_block_vector(const char* const s, // IN - file name,
		  double** v,          // OUT - subvector
		  int* n,              // OUT - vector length
		  MPI_Comm comm)       // IN - communicator
{
  int p;
  int id;
  FILE* infile_ptr = NULL;
  int local_els;
  int* cnt;
  int* dsp;

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  /* Process p-1 opens file, determines number of vector elements,
     and broadcasts this value to the other processes. */
  if (id == p - 1) {
    infile_ptr = fopen(s, "r");
    if (infile_ptr == NULL) *n = 0;
    else {
      size_t result;
      result = fread(n, sizeof(int), 1, infile_ptr);
      if (result != 1) *n = -1;
    }
  }
  MPI_Bcast(n, 1, MPI_INT, p - 1, comm);
  if (*n == 0) MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
  if (*n == -1) MPI_Abort(MPI_COMM_WORLD, READ_FILE_ERROR);

  /* Block mapping of vector elemens to processes. */
  local_els = block_size(id, p, *n);
  create_mixed_xfer_arrays(*n, &cnt, &dsp, comm);

  /* Dynamically allocate vector. */
  *v = malloc(local_els * sizeof(double));
  double* buffer = NULL;
  if (id == p -1) 
    buffer = (double*) malloc(*n * sizeof(double));
  for (int i = 0; i < p - 1; i++) {
    if (id == p - 1) {
      size_t result;
      result = fread(buffer, sizeof(double), *n, infile_ptr);
    }
    MPI_Scatterv(buffer, cnt, dsp, MPI_DOUBLE,
		 *v, local_els, MPI_DOUBLE, p - 1, comm);
  }
  if (id == 0) free(buffer);
}

void
print_subvector(const double* const v, // IN
		int n) // IN
{
  for (int i = 0; i < n; i++) {
    printf("%6.3f ", *(v + i));
  }
}

void
print_block_vector(double* v, // IN
		   int n, // IN
		   MPI_Comm comm) // IN
{
  int id, p;
  MPI_Status status;

  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);

  if (id == 0) {
    print_subvector(v, block_size(id, p, n));
    double* buf;
    for (int i = 1; i < p; i++) {
      buf = malloc(block_size(i, p, n) * sizeof(double));
      MPI_Recv(buf, block_size(i, p, n), MPI_DOUBLE, i, 
	       0, comm, &status);
      print_subvector(buf, block_size(i, p, n));
      free(buf);
    }
    printf("\n\n");
  } else
    MPI_Send(v, block_size(id, p, n), MPI_DOUBLE, 0, 0, comm);
}

void
create_uniform_xfer_arrays(int n, // IN
			   int** count, // OUT
			   int** disp, // OUT
			   MPI_Comm comm) // IN
{
  int p, id;
  
  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  *count = malloc(p * sizeof(int));
  *disp = malloc(p * sizeof(int));
  
  **count = n;
  **disp = 0;

  for (int i = 0; i < p; i++) {
    (*disp)[i] = (*disp)[i-1] + n;
    (*count)[i] = n;
  }
}

int
main(int argc, char* argv[])
{
  int id; // process rank
  int p;  // # of processes
  double* a; // elements of matrix
  int m, n; // # of rows and cols
  double* b; // elements of vector
  int l; // # of vector elements
  double* c_part;
  double* c;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  read_col_striped_matrix(argv[1], &a, &m, &n, MPI_COMM_WORLD);
  print_col_striped_matrix(a, m, n, MPI_COMM_WORLD);
  read_block_vector(argv[2], &b, &l, MPI_COMM_WORLD);
  print_block_vector(b, l, MPI_COMM_WORLD);
  
  /* Each process multiplies its columns of 'a' and vector 'b',
     resulting in a partial sum of product 'c'. */
  c_part = calloc(n, sizeof(double));
  int local_els = block_size(id, p, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < local_els; j++) {
      c_part[i] += *(a + i * local_els + j) * *(b + j);
    }
  }

  c = malloc(m * sizeof(double));
  MPI_Allreduce(c_part, c, m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (id == 0) {
    print_subvector(c, m);
    putchar('\n');
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
