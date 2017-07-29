
#include "euler_flux.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void euler_flux_2D(double *W, double *F, double *G){
	
	/*
	Purpose:
		Compute the Euler Flux vector using the state vector W.
		Store results for both vectors in F and G
	*/

	double ro, u, v, e_tot, P;
	double GAMMA = 1.4;

	ro = W[0];
	u = W[1]/ro;
	v = W[2]/ro;
	e_tot = W[3]/ro;

	P = (GAMMA-1)*ro*(e_tot - 0.5*(u*u+ v*v));

	// F Vector:
	F[0] = ro*u;
	F[1] = ro*u*u + P;
	F[2] = ro*u*v;
	F[3] = ro*u*e_tot + P*u;

	// G Vector:
	G[0] = ro*v;
	G[1] = ro*u*v;
	G[2] = ro*v*v + P;
	G[3] = ro*v*e_tot + P*v;

}


void euler_flux_2D_matrix(double *W, double *F, double *G, int N){

	/*
	Purpose:
		Compute the Euler Flux vector using the state vector W at
		multiple points. Will return a matrix for F and G (multiple
		rows for the different W vectors given). Note all vectors 
		are in column major form
	*/

	// For storing the values in a row of the matrix
	double W_row[4], F_row[4], G_row[4];

	int i,j;

	for(i=0; i<N; i++){

		// Fill the W_row vector:
		for(j=0; j<4; j++){
			W_row[j] = W[j*N + i];
		}

		euler_flux_2D(W_row, F_row, G_row);

		for(j=0; j<4; j++){
			F[j*N+i] = F_row[j];
			G[j*N+i] = G_row[j];
		}

	}

}





