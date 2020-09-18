#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"

double **createMatrix(int N) { //this was givens
  double **M = (double **) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
    M[i] = (double*) malloc(N*sizeof(double));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      M[i][j] = (i == j) ? 1.0 : 0.0;
  return M;
}

void destroyMatrix(int N, double **M) { //this was given
  for (int i = 0; i < N; i++)
    free(M[i]);
  free(M);
}

LUfact *LUfactor(int N, const double **A) { //this will do the factorization
  LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
  LU->N = N;
  LU->LU = createMatrix(N);
  LU->mutate = (short *) malloc(N*sizeof(short));

  // Clone A into LU
  double **A_ = LU->LU;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      A_[i][j] = A[i][j];

  for (int i = 0; i < N; i++)
    LU->mutate[i] = (short) i;

  // actual factorizing goes here
  // most part of this code was based on the one available on wikepedia
  
	int i, j, k, imax;
	double maxa, *ptr, absa;
	for(i = 0; i <= N; i++)
	{
		LU->mutate[i] = i;
	}
	for(i = 0; i< N; i++)
	{
		maxa = 0;
		imax = i;
		
		for(j = i; j < N; j++)
		{
			if((absa = fabs(LU->LU[j][i])) > maxa)
			{
				maxa = absa;
				imax = j;
			}
		}
		
		if(imax != i) //it will stop the factorization when it reaches the input number

		{

			k = LU->mutate[i];
			LU->mutate[i] = LU->mutate[imax];
			LU->mutate[imax] = k;
			
			ptr = A_[i];
			A_[i] = A_[imax];
			A_[imax] = ptr;
			LU->mutate[N]++;
		}
		
		for(k = i + 1; k < N; k++) //puts the numbers on the right place
		{
			A_[k][i] /= A_[i][i];
			
			for(j = i + 1; j < N; j++)
			{
				A_[k][j] -= A_[k][i] * A_[i][j];
			}
		}
	}
  
  return LU;
}

void LUdestroy(LUfact *fact) { //release the matrix
	free(fact->mutate);
	destroyMatrix(fact->N, fact->LU);
	free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) { //solves the factorization
	int N = fact->N;
	short *P = fact->mutate;
	double **A = fact->LU;
	
    for (int i = 0; i < fact->N; i++) {
        x[i] = b[fact->mutate[i]];

        for (int k = 0; k < i; k++)
            x[i] -= fact->LU[i][k] * x[k];
    }

    for (int i = fact->N - 1; i >= 0; i--) {
        for (int k = i + 1; k < fact->N; k++)
            x[i] -= fact->LU[i][k] * x[k];

        x[i] = x[i] / fact->LU[i][i];
    }
}





