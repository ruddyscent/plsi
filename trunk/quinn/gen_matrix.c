/*
 * Generate a random matrix with m rows and n columns where m and n are 
 * given as arguments and write it to a file with the user given name
 * in a binary format.
 *
 * Usage: gen_matrix [filename] [# of rows] [# of columns]
 *
 * Programmed by Kyungwon Chun
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
  if (argc != 4) {
    fprintf(stderr, "Usage: gen_matrix [filename] [# of rows] [# of columns]\n");
    exit(EXIT_FAILURE);
  }

  char* filename = argv[1];
  int m = atoi(argv[2]), n = atoi(argv[3]);

  FILE* file_ptr;
  double* data_ptr;

  file_ptr = fopen(filename, "wb");
  if (file_ptr == NULL) {
    fprintf(stderr, "Error: opening file.\n");   
    exit(EXIT_FAILURE);
  }

  srand((int)file_ptr);
  data_ptr = (double*) malloc(m * n * sizeof(double));
  for (int i = 0; i < m * n; i++) {
    data_ptr[i] = (double)rand() / rand();
  }

  printf("(%d, %d)\n", m, n);
  for (int i = 0; i < m * n; i++) {
    printf("%f ", data_ptr[i]);
    if ((i + 1) % n == 0) printf("\n");
  }  
  
  fwrite(data_ptr, sizeof(double), m * n, file_ptr);
  fclose(file_ptr);

  return EXIT_SUCCESS;
}
