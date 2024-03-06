/*
Ilya Cable
CS 330
Project 2: LU Factorization Engine
10/16/2023
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"
#include <windows.h>
//creates identity matrix NxN
double **createMatrix(int N) {
  double **M = (double **) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
    M[i] = (double*) malloc(N*sizeof(double));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      M[i][j] = (i == j) ? 1.0 : 0.0;
  return M;
}

//de-allocates NxN matrix
void destroyMatrix(int N, double **M) {
  for (int i = 0; i < N; i++)
    free(M[i]);
  free(M);
}

//prints NxM matrix (for de-bugging)
void printMatrix(double **matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n"); // Move to the next row
    }
    printf("\n");
}

//performs the matrix multiplication Pxb, where P is the permutation matrix for partial pivoting.
void transform_b(LUfact *fact,const double *b, double *transformed_b){
  int N = fact->N;
  double resultant;
  //Multiply each corresponding element by each other, add, and place in transformed_b
  for(int row = 0; row<N; row++){
    resultant = 0;
    for(int col = 0; col<N; col++)
      resultant += fact->mutate[row][col] * b[col];
    transformed_b[row] = resultant;
  }
}


LUfact *LUfactor(int N, const double **A) {
  LUfact *LU = (LUfact*) malloc(sizeof(LUfact));

  //Initialize LU data structure:
  LU->N = N;
  LU->U = createMatrix(N);
  LU->L = createMatrix(N);
  LU->mutate = createMatrix(N);

  // Clone A into U
  double **A_ = LU->U;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++){
      A_[i][j] = A[i][j];
    }

  double *temp;//Used for swapping rows
  //Factorize w/ partial pivoting:
  for(int n = 0; n<N; n++){
    int index = n;
    for(int row = n; row<N; row++){
      if(LU->U[row][n]>LU->U[n][n])
        index = row;//get index of largest element in current column
    }
    if(index != n){//swap L, U, and P if necessary:
      //for L, we're going to swap rows (instead of specific indeces), then clear the upper right and re-format after we're done:
      temp = LU->L[n];
      LU->L[n] = LU->L[index];
      LU->L[index] = temp;

      temp = LU->U[n];
      LU->U[n] = LU->U[index];
      LU->U[index] = temp;

      temp = LU->mutate[n];
      LU->mutate[n] = LU->mutate[index];
      LU->mutate[index] = temp;
    }

    //Now do Gaussian Elimination:
    double coeff;
    for(int row = n+1; row<N; row++){
      coeff = LU->U[row][n]/LU->U[n][n];
      LU->L[row][n] = coeff;//place coefficient into 'L'
      for(int col = n; col<N; col++){
        LU->U[row][col] -= LU->U[n][col]*coeff;//'U' becomes upper-triangular
      } 
    }
  }

  //before we return, we need to re-format 'L' (make sure 1's on diagonal + 0's in upper right triangle):
  for(int row = 0; row<N; row++){
    for(int col = 0; col<N; col++){
      if(row==col)
        LU->L[row][col] = 1;
      else if(col>row)
        LU->L[row][col] = 0;
    }
  }
  return LU;
}

//De-allocates all memory used for LU data structure
void LUdestroy(LUfact *fact) {
  destroyMatrix(fact->N,fact->L);
  destroyMatrix(fact->N,fact->U);
  destroyMatrix(fact->N,fact->mutate);
  free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) {
  
  //start by solving Lz=b
  int N = fact->N;

  //Since we used permutation matrix 'P', we first need to multiply Pxb to get proper row locations in 'b':
  double trans_b[fact->N];
  transform_b(fact,b,trans_b);

  //we'll use 'trans_b' instead of 'b' now:
  double z[N];
  z[0] = trans_b[0];//we can make this assumption since 'L' will always have a '1' in the upper left hand corner.
  //Perform back-substitution:
  for(int row = 1; row<N; row++){
    double constant = 0;
    for(int col = 0; col<row; col++){
      constant += fact->L[row][col]*z[col];//add up all constant elements
    }
    z[row] = trans_b[row]-constant;//subtract all constants from the right side of equation.  No need to divide since diagonal will always be 1's
  }

  //now solve Ux = z

  //start from lower right corner this time:
  x[N-1] = z[N-1]/fact->U[N-1][N-1];//no guarantee that 1's will be on diagonal.
  double constant;
  double temp;//using a temp variable to hold subtraction result massively improves accuracy.  Sensitive subtraction?
  //Perform back-sub:
  for(int row = N-2; row>=0; row--){
    constant = 0;
    for(int col = row+1; col<N; col++){
      constant+=fact->U[row][col]*x[col];//add all constants
    }
    temp = z[row]-constant;
    x[row] = temp/fact->U[row][row];//unlike the previous back-sub, we don't have all 1's on the diagonal here.
  }
}