
#include "fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void flux_LF(double *WIn, double *WOut, double *FIn, double *FOut, double *GIn,
	double *GOut, double *FComm, double *nIn, int P, int Neq){

	/*
	Purpose:
		Compute the numerical flux at all the integration nodes on the face
		using the Lax numerical flux. Result returned is a matrix for the 
		numerical flux for all equations at each integration points :
		matrix of size nSol x nEq
	*/

	int sol_i, eq_i;

	double GAMMA = 1.4;

	double FAvg, GAvg, diffW;
	double roIn, uIn, vIn, eIn, PIn;
	double roOut, uOut, vOut, eOut, POut;
	double lambdaIn, lambdaOut, lambda;
	double nx, ny;

	for(sol_i=0; sol_i<(P+1); sol_i ++){
		// Loop over the integration nodes

		roIn = WIn[sol_i];
		uIn = WIn[sol_i + 1*(P+1)]/roIn;
		vIn = WIn[sol_i + 2*(P+1)]/roIn;
		eIn = WIn[sol_i + 3*(P+1)]/roIn;
		PIn = (GAMMA-1)*roIn*(eIn - 0.5*(uIn*uIn + vIn*vIn));

		roOut = WOut[sol_i];
		uOut = WOut[sol_i + 1*(P+1)]/roOut;
		vOut = WOut[sol_i + 2*(P+1)]/roOut;
		eOut = WOut[sol_i + 3*(P+1)]/roOut;
		POut = (GAMMA-1)*roOut*(eOut - 0.5*(uOut*uOut + vOut*vOut));

		lambdaIn = sqrt(uIn*uIn + vIn*vIn) + sqrt(fabs(GAMMA*PIn/roIn));
		lambdaOut = sqrt(uOut*uOut + vOut*vOut) + sqrt(fabs(GAMMA*POut/roOut));

		if(lambdaIn > lambdaOut){
			lambda = lambdaIn;
		} else{
			lambda = lambdaOut;
		}

		nx = nIn[sol_i];
		ny = nIn[sol_i + (P+1)];		

		for(eq_i=0; eq_i<Neq; eq_i ++){
			// Loop over the equations

			FAvg = 0.5*(FIn[eq_i*(P+1) + sol_i] + FOut[eq_i*(P+1) + sol_i]);
			GAvg = 0.5*(GIn[eq_i*(P+1) + sol_i] + GOut[eq_i*(P+1) + sol_i]);
			diffW = WIn[eq_i*(P+1) + sol_i] - WOut[eq_i*(P+1) + sol_i];

			FComm[eq_i*(P+1) + sol_i] = nx*FAvg + ny*GAvg + 0.5*lambda*diffW;

		}
	}

}


