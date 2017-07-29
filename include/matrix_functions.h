
#ifndef DG__matrix_functions_h__INCLUDED
#define DG__matrix_functions_h__INCLUDED


void mm_CNN(int m, int k, int n, double *A, double *B, double *C);
void mm_inv_d(int N, double *matrix);
void mm_inv_d_secondInPlace(int N, double *matrix, double *matrixInv);
double *mm_inv_d_alloc(int N, double *matrix);
void mm_transposeR_d(int m, int n, double *A);

#endif //DG__matrix_functions_h__INCLUDED




