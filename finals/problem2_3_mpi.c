#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <mpi.h>

const int IM = 500;
const int JM = 500;
const int itermax = 100000;
const float tolerance = 1.0e-7;
const float bc[] = {10.0, // left boundary
		    10.0, // right boundary
		    10.0, // top boundary
		    20.0};// bottom boundary

/*
 * utility functions and macros
 */
#define u(i,j) u[(i)*jm+(j)]
#define uo(i,j) uo[(i)*jm+(j)]
#define buf(i,j) buf[(i)*jm+(j)]

static inline int
block_size(int id, int p, int m)
{
  int result;
  
  m += 2;
  result =  m / p + (id < m % p ? 1 : 0);

  if (id != 0) result++;
  if (id != p - 1) result++;

  return result;
}

const int UPAWARD = 0;
const int DWAWARD = 1;

void
init_block(float** u_ptr,    // OUT
	   float** uo_ptr,   // OUT
	   float x_size, // IN
	   float y_size, // IN
	   const float* const bc, // IN
	   MPI_Comm comm) // IN
{
  int p, id;
  float* u;

  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &id);
  
  const int im = block_size(id, p, x_size);
  const int jm = block_size(0, 1, y_size);

  *u_ptr = calloc(im * jm, sizeof(float));
  *uo_ptr = malloc(im * jm * sizeof(float));
  
  /*
   * Set boundary conditions.
   */
  u = *u_ptr;
  // top line
  if (id == 0) {
    for (int j = 0; j < jm; j++) {
      u(0, j) = bc[0];  
    }
  }
  // bottom line
  if (id == p - 1) {
    for (int j = 0; j < jm; j++) {
      u(im - 1, j) = bc[1];
    }
  }
  // left line
  for (int i = 0; i < im; i++) {
    u(i, 0) = bc[2];
  }
  // right line
  for (int i = 0; i < im; i++) {
      u(i, jm - 1) = bc[3];
  }
}

void
jaccobi(float* const u,  // INPLACE
	float* const uo, // INPLACE
	int x_size, // IN
	int y_size, // IN
	float tolerance, // IN
	int itermax, // IN
	MPI_Comm comm) // IN
{
  float error = FLT_MAX;
  int iter = 0;

  int id, p;
  MPI_Request req_top[2], req_bot[2];
  MPI_Status stat_top[2], stat_bot[2];

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  const int im = block_size(id, p, x_size);
  const int jm = block_size(0, 1, y_size);
  
  if (id != 0) {
    MPI_Send_init(&u(1, 0), jm, MPI_FLOAT, id - 1, 
		  UPAWARD, MPI_COMM_WORLD, &req_top[0]);
    MPI_Recv_init(&u(0, 0), jm, MPI_FLOAT, id - 1, 
		  DWAWARD, MPI_COMM_WORLD, &req_top[1]);
  }

  if (id != p - 1) {
    MPI_Send_init(&u(im - 2, 0), jm, MPI_FLOAT, id + 1, 
		  DWAWARD, MPI_COMM_WORLD, &req_bot[0]);
    MPI_Recv_init(&u(im - 1, 0), jm, MPI_FLOAT, id + 1, 
		  UPAWARD, MPI_COMM_WORLD, &req_bot[1]);
  }

  while (iter <= itermax && error > tolerance) {
    // Store old data.
    memcpy(uo, u, im * jm * sizeof(float));

    // Jacobi method
    for (int i = 1; i < im - 1; i++) {
      for (int j = 1; j < jm - 1; j++) {
  	u(i, j) = 0.25 * (uo(i - 1, j) + uo(i + 1, j) + uo(i, j - 1) + uo(i, j + 1));
      }
    }

    if (id != 0) 
      MPI_Startall(2, req_top);
    if (id != p - 1) 
      MPI_Startall(2, req_bot);
    
    if (id != 0) 
      MPI_Waitall(2, req_top, stat_top);
    if (id != p - 1)
      MPI_Waitall(2, req_bot, stat_bot);

    // Calculate error.
    error = 0.0;
    int error_cnt = 0;
    for (int i = 1; i < im - 1; i++) {
      for (int j = 1; j < jm - 1; j++) {
  	error += pow(u(i, j) - uo(i, j), 2);
	error_cnt++;
      }
    }
    
    float error_buf;
    MPI_Allreduce(&error, &error_buf, 1, MPI_FLOAT, MPI_SUM,
		  MPI_COMM_WORLD);
    error = error_buf;
    iter++;
  }
}

void
print_row(const float* const row, // IN
	  int row_ord, // IN
	  int col_len) // IN
{
  for (int j = 0; j < col_len; j++) {
    printf("%d %d %f \n", row_ord, j, *(row + j));
  }
}

void
print_block(float* const u, // IN
	    int x_size,           // IN
	    int y_size,           // IN
	    MPI_Comm comm)        // IN
{
  int id, p;

  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &p);

  const int im = block_size(id, p, x_size);
  const int jm = block_size(0, 1 ,y_size);

  int i;
  // Print top nodes, first.
  if (id == 0) {
    for (i = 0; i < im - 1; i++) {
      print_row(&u(i,0), i, jm);
    }
    if (p == 1) {
      print_row(&u(i, 0), i, jm);
    }
  }
  // Print the other nodes.
  if (id == 0) {
    for (int k = 1; k < p; k++) {
      const int im_src = block_size(k, p, x_size) - 1;
      const int jm_src = block_size(0, 1, y_size);
      MPI_Status status;
      
      float* buf = malloc(im_src * jm_src * sizeof(float));
      MPI_Recv(buf, im_src * jm_src, MPI_FLOAT, k, k, comm, &status);
      
      for (int l = 0; l < im_src - 1; l++) {
  	print_row(&buf(l, 0), i, jm_src);
	i++;
      }

      // For the bottome node, print the boundary.
      if (k == p - 1) {
      	print_row(&buf(im_src - 1, 0), i, jm_src);
      }

      free(buf);
    }
  }
  else {
    MPI_Send(&u(1,0), (im - 1) * jm, MPI_FLOAT, 0, id, comm);
  }
}

int 
main(int argc, char* argv[]) {
  const int x_size = IM;
  const int y_size = JM;

  float* u;
  float* uo;

  MPI_Init(&argc, &argv);

  init_block(&u, &uo, x_size, y_size, bc, MPI_COMM_WORLD);

  jaccobi(u, uo, x_size, y_size, tolerance, itermax,
  	  MPI_COMM_WORLD);
  
  print_block(u, x_size, y_size, MPI_COMM_WORLD);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
