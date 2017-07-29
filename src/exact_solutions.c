

#include "exact_solutions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "S_DB.h"

void uniform_solution_InternalSubsonic(double *XYZ, double*W){
	
	double rho, u, v, w, P, eTot;

	rho = DB.rhoInf;
	u = DB.MInf*DB.cInf;
	v = 0.0;
	w = 0.0;
	P = DB.pInf;

	eTot = P/((GAMMA-1)*rho) + (u*u + v*v)/2.;

	// Compute the state vector
	W[0] = rho;
	W[1] = rho*u;
	W[2] = rho*v;
	W[3] = rho*eTot;

}

void exact_solution_IsentropicVortex(double *XYZ, double *W){

	/*
	Purpose:
		- Compute the W vector at the given location

	Parameters:
		XYZ : Vector for the point on the physical domain
		W : The returned state vector at the XYZ point
	*/

	double 	CONST_Gamma, eStrength, x0, y0, x, y,
			r, uinf, vinf, u, v, delT, ro, P, eTot;

	CONST_Gamma = 1.4;

	//-------------------------------------------------
	//				Vortex Properties

	eStrength = 5.0;

	x0 = 0.0;
	y0 = 0.0;

	//-------------------------------------------------


	x = XYZ[0];
	y = XYZ[1];

	r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

	uinf = 1.0;
	vinf = 0.0;

	u = uinf - ((eStrength)/(2.*M_PI))*exp(0.5*(1.-r*r))*(y-y0);
	v = vinf + ((eStrength)/(2.*M_PI))*exp(0.5*(1.-r*r))*(x-x0);

	delT = -(((CONST_Gamma-1.)*(eStrength*eStrength))/(8.*CONST_Gamma*M_PI*M_PI))*exp(1.-r*r);

	ro = pow((1.+delT),(1./(CONST_Gamma-1)));
	P = pow(ro,(CONST_Gamma));

	eTot = P/((CONST_Gamma-1)*ro) + (u*u + v*v)/2.;

	// Compute the state vector
	W[0] = ro;
	W[1] = ro*u;
	W[2] = ro*v;
	W[3] = ro*eTot;

}


