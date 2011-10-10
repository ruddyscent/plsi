#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define N_P 2000
#define STEP 10000                                            

const int MPI_IO_NODE = 0;

static inline int
block_size(int id,   // process rank
	   int p,    // # of processes
	   int size) // # of elements
{
  return size / p + (id < size % p ? 1 : 0);
}

void
gen_trans_arrays(int** sendcnts, // OUT
		 int** disps,    // OUT
		 int size,       // IN
		 MPI_Comm comm)  // IN
{
  int p, id;
  
  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);
  
  *sendcnts = malloc(p * sizeof(int));
  *disps = malloc(p * sizeof(int));
  
  **sendcnts = block_size(0, p, size);
  **disps = 0;
  for (int i = 1; i < p; i++) {
    (*sendcnts)[i] = block_size(i, p, size);
    (*disps)[i] = (*disps)[i-1] + (*sendcnts)[i-1];
  }
}

void
read_input(float** data,  // OUT
	   int size,      // IN
	   MPI_Comm comm) // IN
{
  int id, p;
  float *read_buffer;
  int* sendcnts;
  int* displs;
  int bsize;
  
  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  bsize = block_size(id, p, size);
  read_buffer = malloc(size * sizeof(float));
  if (id == MPI_IO_NODE) {
    for (int i = 0; i < size; i++) {
      scanf("%f", &read_buffer[i]);
      read_buffer[i] *= 10000000.0;
    }
  }
  *data = malloc(bsize * sizeof(float));
  gen_trans_arrays(&sendcnts, &displs, size, comm);
  MPI_Scatterv(read_buffer, sendcnts, displs, MPI_FLOAT, 
	       *data, bsize, MPI_FLOAT, MPI_IO_NODE, comm);
  free(read_buffer);
  free(sendcnts);
  free(displs);
}

void
print_result(float* const data,
	     int size,
	     MPI_Comm comm)
{
  int id, p;
  int bsize;

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  bsize = block_size(id, p, size);

  if (id == 0) {
    float* buf;
    int i, k;
    int dest_bsize;
    MPI_Status status;

    for (i = 0; i < bsize; i++) {
      printf("%d %f \n", i + 1, data[i]);
    }

    for (int j = 1; j < p; j++) {
      dest_bsize = block_size(i, p, size);
      buf = malloc(dest_bsize * sizeof(float));
      MPI_Recv(buf, dest_bsize, MPI_FLOAT, j, j, comm,
  	       &status);

      for (k = 0; k < dest_bsize; k++) {
  	printf("%d %f \n", k + i + 1, buf[k]);
      }
      i += k;
      free(buf);
    }
  }
  else {
    MPI_Send(data, bsize, MPI_FLOAT, 0, id, comm);
  }
}

int
main(int argc, char* argv[])
{
  float* x;
  float* force;
  float f_jk;
  int id, p;
  int bsize;
  const int size = N_P;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  bsize = block_size(id, p, size);
  read_input(&x, size, MPI_COMM_WORLD);

  if (id == 0) printf("start \n");

  force = malloc(bsize * sizeof(float));

  for (int i = 0; i < STEP; i++) {  
    for (int l = 0; l < bsize; l++) force[l] = 0.0;
    for (int j = 0; j < bsize - 1; j++) {
      for (int k = j + 1; k < bsize; k++) {
        f_jk = 1.0 / (x[k] - x[j]);
        force[j] += f_jk;
        force[k] -= f_jk;
      }
    }
    for (int k = 0; k < bsize; k++) {
      x[k] += force[k];
    }
  }

  free(force);

  print_result(x, size, MPI_COMM_WORLD);
  free(x);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
