/*
 * Generate a random vector with m rows and n columns where m and n are 
 * given as arguments and write it to a file with the user given name
 * in a binary format.
 *
 * Usage: gen_vector [filename] [# of elements]
 *
 * Programmed by Kyungwon Chun
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv)
{
  if (argc != 3) {
    fprintf(stderr, "Usage: gen_vector [filename] [# of elements]\n");
    exit(EXIT_FAILURE);
  }

  const char* const filename = argv[1];
  const int m = atoi(argv[2]);

  FILE* file_ptr;
  double* data_ptr;

  file_ptr = fopen(filename, "wb");
  if (file_ptr == NULL) {
    fprintf(stderr, "Error: opening file.\n");   
    exit(EXIT_FAILURE);
  }

  srand(time(NULL));
  data_ptr = (double*) malloc(m * sizeof(double));
  for (int i = 0; i < m; i++) {
    data_ptr[i] = (double)rand() / rand();
  }

  printf("(%d)\n", m);
  for (int i = 0; i < m; i++) {
    printf("%6.3f ", data_ptr[i]);
  }
  putchar('\n');

  fwrite(&m, sizeof(int), 1, file_ptr);
  fwrite(data_ptr, sizeof(double), m, file_ptr);
  fclose(file_ptr);
  free(data_ptr);

  return EXIT_SUCCESS;
}
