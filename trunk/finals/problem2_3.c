#include <stdio.h>
#include <math.h>
#define IM 500
#define JM 500

int main(int argc, char* argv[]) {
  int im1,jm1,i,j;
  int is,ie,js,je;
  int iter,itermax;
  float tolerance, error;
  float bc[4],u[IM+2][JM+2];
  float uo[IM+2][JM+2];
  itermax=100000;
  tolerance=1.0e-7;
  bc[0]=10.0;
  bc[1]=10.0;
  bc[2]=10.0;
  bc[3]=20.0;

  //initialize
  im1=IM+1;
  jm1=JM+1;
  for(i=0;i<IM+1;i++) {
    for(j=0;j<JM+1;j++) {
      u[i][j]=0.0;
    }
  }
  //boundary conditions
  for(j=0;j<jm1+1;j++) {
    u[0][j]=bc[0]; //left line
    u[im1][j]=bc[1]; //right line
  }
  for(i=0;i<im1+1;i++) {
    u[i][0]=bc[2]; //bottom line
    u[i][jm1]=bc[3]; //top line
  }
  //set computation range
  is=1;
  ie=IM;
  js=1;
  je=JM;
  //main routine
  iter=0;
  error=1000.0;
  while(iter<=itermax && error>tolerance) {
    //store old data
    for(i=0;i<IM+2;i++){
      for(j=0;j<JM+2;j++){
  	uo[i][j]=u[i][j];
      }
    }
    //jacobi
    for(i=is;i<ie+1;i++){
      for(j=js;j<je+1;j++){
  	u[i][j]=(uo[i-1][j]+uo[i+1][j]+uo[i][j-1]+uo[i][j+1])/4.0;
      }
    }
    //error
    error=0.0;
    int error_cnt = 0;
    for(i=is;i<ie+1;i++){
      for(j=js;j<je+1;j++){
  	error=error+(u[i][j]-uo[i][j])*(u[i][j]-uo[i][j]);
	error_cnt++;
      }
    }
    iter+=1;
  }
  for(i=0;i<im1+1;i++) {
    for(j=0;j<jm1+1;j++) {
      printf("%d %d %f \n",i,j,u[i][j]);
    }
  }
  //store old data
  for(i=0;i<IM+2;i++){
    for(j=0;j<JM+2;j++){
      uo[i][j]=u[i][j];
    }
  }
  //jacobi
  for(i=is;i<ie+1;i++){
    for(j=js;j<je+1;j++){
      u[i][j]=(uo[i-1][j]+uo[i+1][j]+uo[i][j-1]+uo[i][j+1])/4.0;
    }
  }
  //error
  error=0.0;
  for(i=is;i<ie+1;i++){
    for(j=js;j<je+1;j++){
      error=error+(u[i][j]-uo[i][j])*(u[i][j]-uo[i][j]);
    }
  }
}
