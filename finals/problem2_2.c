#include <stdio.h>
#include <math.h>

#define N_P  2000
#define STEP 10000

int main(int argc, char* argv[]) {
  float x[N_P], force[N_P], f_jk;

  for (int i = 0; i < N_P; i++) {
    scanf("%f", &x[i]);
    x[i] = x[i] * 10000000.0;
  }
  printf("start \n");
  for (int i = 0; i < STEP; i++) {
    for(int j = 0; j < N_P; j++) force[j] = 0.0;
    for(int j = 0; j < N_P - 1; j++) {
      for(int k = j + 1; k < N_P; k++) {
        f_jk = 1.0 / (x[k] - x[j]);
        force[j] += f_jk;
        force[k] -= f_jk;
      }
    }
    for(int j = 0; j < N_P; j++) {
      x[j] += force[j];
    }
  }
  for(int i = 0; i < N_P; i++) printf("%d %f \n", i + 1, x[i]);
}
