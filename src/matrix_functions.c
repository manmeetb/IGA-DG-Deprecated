

#include "matrix_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <Accelerate/Accelerate.h>



void mm_CNN(int m, int k, int n, double *A, double *B, double *C){

	/*

	A = m x k
	B = k x n
	C = m x n

	Purpose:
		Compute matrix-matrix and matrix-vector products:
			0) return C = A*B
			1) Column-major storage of A
			2) Column-major storage of B and C
		
		Comments:
			The 'CNN' function name refers to:
				C : Column-major output
	 			N : No Transpose of A
	 			N : No transpose of B

	*/

	int lda = m;
	int ldb = k;
	int ldc = m;

	cblas_dgemm (102, 111, 111, m, n, k, 1.0, A, lda, B, ldb, 0.0, C, ldc);

}

void mm_transposeR_d(int m, int n, double *A){
	
	/*
	Purpose:
		Compute the transpose of the matrix A (in column major form) which is of size
		m x n. Return the result in column major form.
	*/

	double *transA;
	int i,j,k;

	transA = malloc(m*n* sizeof *transA);

	k=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			// keep row fixed, loop through columns

			transA[k++] = A[j*m + i];
		}
	}

	A = transA;

}

void mm_inv_d(int N, double *matrix) {

	/*
	Purpose:
		Compute the inverse of a matrix using LAPACK. 
		Matrix is in column major form and result is 
		returned in column major form as well.
	*/

    int error=0;
    int *pivot = malloc(N*sizeof(int));
    double *workspace = malloc(N*sizeof(double));

    /*  LU factorisation */
    dgetrf_(&N, &N, matrix, &N, pivot, &error);

    if (error != 0) {
        free(pivot);
        free(workspace);
    }

    /*  matrix inversion */
    dgetri_(&N, matrix, &N, pivot, workspace, &N, &error);

    if (error != 0) {
        free(pivot);
        free(workspace);
    }

    // Now, inversion was done in row major form (transpose was inverted)
    // Compute the transpose of the inverse to get the inverse in column
    // major form
    mm_transposeR_d(N, N, matrix);

    free(pivot);
    free(workspace);
}

void mm_inv_d_secondInPlace(int N, double *matrix, double *matrixInv) {

	/*
	Purpose:
		Compute the inverse of a matrix using LAPACK. 
		Matrix is in column major form and result is 
		returned in column major form as well.
	*/

	int i;
	for(i=0; i<N*N; i++){
		matrixInv[i] = matrix[i];
	}

    int error=0;
    int *pivot = malloc(N*sizeof(int));
    double *workspace = malloc(N*sizeof(double));

    /*  LU factorisation */
    dgetrf_(&N, &N, matrixInv, &N, pivot, &error);

    if (error != 0) {
        free(pivot);
        free(workspace);
    }

    /*  matrix inversion */
    dgetri_(&N, matrixInv, &N, pivot, workspace, &N, &error);

    if (error != 0) {
        free(pivot);
        free(workspace);
    }

    // Now, inversion was done in row major form (transpose was inverted)
    // Compute the transpose of the inverse to get the inverse in column
    // major form
    mm_transposeR_d(N, N, matrixInv);

    free(pivot);
    free(workspace);
}

double *mm_inv_d_alloc(int N, double *matrix){
	
	/*
	Purpose:
		Compute the inverse of the matrix but do not modify
		the one that is sent as an argument.
	*/

	int i;
	double *matrix_inv;
	matrix_inv = malloc(N*N * sizeof *matrix_inv);

	for(i=0; i<N*N; i++){
		matrix_inv[i] = matrix[i];
	}

	mm_inv_d(N, matrix_inv);

	return matrix_inv;
}



