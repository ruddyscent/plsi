#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N_P 10000000
#define N_W 1000

int main(int argc, char* argv[]) {
  int i,j;
  int part_x[N_P],part_y[N_P];
  int **pos;
  int temp;

  for(i=0;i<N_P;i++) {
    part_x[i]=0;
    part_y[i]=0;
  }
  srand(0.5);
  for(i=0;i<N_P;i++) {
    for(j=0;j<N_W;j++) {
      temp=rand()%4;
      if(temp == 0) {
	part_x[i]+=1;
      }
      else if(temp == 1) {
	part_y[i]+=1;
      }
      else if(temp == 2) {
	part_x[i]-=1;
      }
      else {
	part_y[i]-=1;
      }
    }
  }
  pos=malloc(sizeof(int*)*(2*N_W+1));
  for(i=0;i<2*N_W+1;i++) {
    pos[i]=malloc(sizeof(int)*(2*N_W+1));
  }
  for(i=0;i<2*N_W+1;i++) {
    for(j=0;j<2*N_W+1;j++) {
      pos[i][j]=0;
    }
  }

  for(i=0;i<N_P;i++) {
    pos[part_x[i]+N_W][part_y[i]+N_W]+=1;
  }

  for(i=0;i<2*N_W+1;i++) {
    for(j=0;j<2*N_W+1;j++) {
      if(pos[i][j] != 0) printf("%d %d %d \n",i-N_W,j-N_W,pos[i][j]);
    }
  }
}
