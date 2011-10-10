#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N_P 1000000
#define N_W 1000                                                   

int 
main(int argc, char* argv[])
{
  int id, p;
  int ista, iend;
  int part_x[N_P], part_y[N_P], tmp_x[N_P], tmp_y[N_P];
  int **pos;
  int temp;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  for (int i = 0; i < N_P; i++) {
    part_x[i] = 0;
    part_y[i] = 0;
  }

  ista = id * (N_P / p);
  iend = ista + (N_P / p);

  srand(0.5 + (double)id);

  for (int i = ista; i < iend; i++) {
    for (int j = 0; j < N_W; j++) {
      temp = rand() % 4;
      
      switch (temp) {
      case 0:
	part_x[i] += 1;
	break;
      case 1:
	part_y[i] += 1;
	break;
      case 2:
	part_x[i] -= 1;
	break;
      default:
	part_y[i] -= 1;
	break;
      }
    }
  }
  
  MPI_Reduce(part_x, tmp_x, N_P, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(part_y, tmp_y, N_P, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  for (int i = 0; i < N_P; i++) {
    part_x[i] = tmp_x[i];
    part_y[i] = tmp_y[i];
  }
  
  if (id == 0) {
    pos = malloc(sizeof(int*) * (2 * N_W + 1));
    
    for (int i = 0; i < 2 * N_W + 1; i++) {
      pos[i] = calloc(2 * N_W + 1, sizeof(int));
    }

    for (int i = 0; i < N_P; i++) {
      pos[part_x[i] + N_W][part_y[i] + N_W] += 1;
    }

    for (int i = 0; i < 2 * N_W + 1; i++) {
      for (int j = 0; j < 2 * N_W + 1; j++) {
	if (pos[i][j] != 0) {
	  printf("%d %d %d \n", i - N_W, j - N_W, pos[i][j]);
	}
      }
    }
  }

  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
