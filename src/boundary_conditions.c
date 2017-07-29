
#include "boundary_conditions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "euler_flux.h"
#include "S_DB.h"
#include "Parameters.h"

void boundary_SlipWall(double *WL, double *nL, double *WB, 
	double *FB, double *GB, int n){

	/*
	Implement the slip wall boundary condition at the wall. This is done
	weakly through the numerical flux by making the state in a ghost cell
	flipped relative to the wall for the velocity so that the flux in
	velocity normal to the wall is 0.
	- Note, due to how the boundary state is computed, it will not need
		to be flipped when computing the numerical flux.

	n = number of integration nodes for this face
	*/


	// Standard datatypes
	int i, j;

	// The boundary state
	double 	*rhoB, *rhouB, *rhovB, *rhoEB, rhoVL;
	// The inner volume state
	double	*rhoL, *rhouL, *rhovL, *rhoEL;

	// For computing the euler flux vector at each boundary integration node point
	double W[4], F[4], G[4];

	double *nL_x, *nL_y;

	nL_x = &nL[0];
	nL_y = &nL[n];

	// Conservative State Vector at each integration node
	rhoL  = &WL[n*0];
	rhouL = &WL[n*1];
	rhovL = &WL[n*2];
	rhoEL = &WL[n*3];

	rhoB  = &WB[n*0];
	rhouB = &WB[n*1];
	rhovB = &WB[n*2];
	rhoEB = &WB[n*3];

	// Density and total energy are equivalent in ghost cell
	for (i = 0; i < n; i++) {
		rhoB[i] = rhoL[i];
		rhoEB[i]   = rhoEL[i];
	}

	// Set velocity components of the boundary condition
	for (i = 0; i < n; i++) {
		rhoVL = nL_x[i]*rhouL[i]+nL_y[i]*rhovL[i];

		rhouB[i] = rhouL[i]-2.0*rhoVL*nL_x[i];
		rhovB[i] = rhovL[i]-2.0*rhoVL*nL_y[i];
	}
	
	// Set up the flux vectors for the boundary state
	for (i=0; i<n; i++){
		
		W[0] = rhoB[i];
		W[1] = rhouB[i];
		W[2] = rhovB[i];
		W[3] = rhoEB[i];

		euler_flux_2D(W, F, G);

		for(j=0; j<4; j++){
			FB[j*n + i] = F[j];
			GB[j*n + i] = G[j];
		}
	}

	if(DB.Testing == 1){
		// Print the state at the wall:
		printf("Boundary State (SLIPWALL): \n");

		for (i=0; i<n; i++){
			printf("	Face Node : %d \n", i);
			printf("		[nx, ny] : [%f, %f] \n", nL_x[i], nL_y[i]);
			printf("		rhoL, rhoB : %f,  %f \n", rhoL[i], rhoB[i]);
			printf("		rhouL, rhouB : %f,  %f \n", rhouL[i], rhouB[i]);
			printf("		rhovL, rhovB : %f,  %f \n", rhovL[i], rhovB[i]);
			printf("		rhoEL, rhoEB : %f,  %f \n", rhoEL[i], rhoEB[i]);
		}

	}


}

void boundary_BackPressure(double *WL, double *nL, double *WB, 
	double *FB, double *GB, int n) {

	/*
	 *	Purpose:
	 *		Impose back Pressure (outflow) boundary condition.
	 *
	 *	References:
	 *		Carlson(2011): 2.4
	 */


	// Standard datatypes
	int i, j;

	// The boundary state
	double 	*rhoB, *rhouB, *rhovB, *rhoEB;
	// The inner volume state
	double	*rhoL, *rhouL, *rhovL, *rhoEL;

	double T_i, p_i, rho_i, rho_i_inv, 
			V_i, u_i, v_i, eTot_i, p_b, rho_b;

	double c_i, c_i_2;

	// For computing the euler flux vector at each boundary integration node point
	double W[4], F[4], G[4];

	double *nL_x, *nL_y;

	nL_x = &nL[0];
	nL_y = &nL[n];

	// Conservative State Vector at each integration node
	rhoL  = &WL[n*0];
	rhouL = &WL[n*1];
	rhovL = &WL[n*2];
	rhoEL = &WL[n*3];

	rhoB  = &WB[n*0];
	rhouB = &WB[n*1];
	rhovB = &WB[n*2];
	rhoEB = &WB[n*3];

	for (i=0; i<n; i++){
		// Loop over all the n face integration nodes


		// Get inner volume's states
		rho_i = rhoL[i];
		rho_i_inv = 1./rho_i;

		u_i = rhouL[i]*rho_i_inv;
		v_i = rhovL[i]*rho_i_inv;

		V_i = u_i*u_i + v_i*v_i;
		V_i = sqrt(V_i);

		eTot_i = rhoEL[i]*rho_i_inv;

		p_i = (GAMMA-1)*rho_i*(eTot_i - 0.5*V_i*V_i);

		c_i_2 = GAMMA*p_i/rho_i;  // speed of sound squared
		c_i = sqrt(c_i_2);

		if(fabs(V_i) >= c_i){
			// Supersonic Flow
			// 		p_b = p_i (set pressure has no influence)
			p_b = p_i;

		} else{
			// Subsonic Flow
			// 		p_b = p_set

			p_b = DB.pBack;
		}	

		// Set the boundary state:
		//	- T_b = T_i (adiabatic flow with V constant)

		rho_b = GAMMA*p_b/c_i_2;

		rhoB[i] = rho_b;
		rhouB[i] = rho_b*u_i;
		rhovB[i] = rho_b*v_i;

		rhoEB[i] = rho_b*(p_b/((GAMMA-1)*rho_b) + 0.5*V_i*V_i);

	}
	
	// Set up the flux vectors for the boundary state
	for (i=0; i<n; i++){
		
		W[0] = rhoB[i];
		W[1] = rhouB[i];
		W[2] = rhovB[i];
		W[3] = rhoEB[i];

		euler_flux_2D(W, F, G);

		for(j=0; j<4; j++){
			FB[j*n + i] = F[j];
			GB[j*n + i] = G[j];
		}
	}

	if(DB.Testing == 1){
		// Print the state at the wall:
		printf("Boundary State (BACK PRESSURE): \n");

		for (i=0; i<n; i++){
			printf("	Face Node : %d \n", i);
			printf("		[nx, ny] : [%f, %f] \n", nL_x[i], nL_y[i]);
			printf("		rhoL, rhoB : %f,  %f \n", rhoL[i], rhoB[i]);
			printf("		rhouL, rhouB : %f,  %f \n", rhouL[i], rhouB[i]);
			printf("		rhovL, rhovB : %f,  %f \n", rhovL[i], rhovB[i]);
			printf("		rhoEL, rhoEB : %f,  %f \n", rhoEL[i], rhoEB[i]);
		}

	}

}

void boundary_Total_TP(double *WL, double *nL, double *WB, 
	double *FB, double *GB, int n){

	/*
	 *	Purpose:
	 *		Impose total (P)ressure/(T)emperature (inflow) boundary condition.
	 *
	 *	Comments:
	 *		eq. (38/47) in Carlson(2011) implies that the velocity should be normal to the boundary. As the direction of
	 *		the flow velocity cannot be known, this implies that this boundary condition is not physically correct...
	 *
	 *	References:
	 *		Carlson(2011): 2.7
	 *		Toro(2009): (3.9), (8.58)
	 */


	// Standard datatypes
	int i, j;

	// The boundary state
	double 	*rhoB, *rhouB, *rhovB, *rhoEB;
	// The inner volume state
	double	*rhoL, *rhouL, *rhovL, *rhoEL;

	double T_i, p_i, rho_i, rho_i_inv, 
			V_i, Vn_i, u_i, v_i, eTot_i;

	double H_i, R_i; 

	double 	p_Total = DB.p_Total,
			T_Total = DB.T_Total;

	double c_i, c_i_2;

	// For computing the euler flux vector at each boundary integration node point
	double W[4], F[4], G[4];

	double *nL_x, *nL_y;

	nL_x = &nL[0];
	nL_y = &nL[n];

	// Conservative State Vector at each integration node
	rhoL  = &WL[n*0];
	rhouL = &WL[n*1];
	rhovL = &WL[n*2];
	rhoEL = &WL[n*3];

	rhoB  = &WB[n*0];
	rhouB = &WB[n*1];
	rhovB = &WB[n*2];
	rhoEB = &WB[n*3];

	for (i=0; i<n; i++){
		
		// Get inner volume's states
		rho_i = rhoL[i];
		rho_i_inv = 1./rho_i;

		u_i = rhouL[i]*rho_i_inv;
		v_i = rhovL[i]*rho_i_inv;

		V_i = u_i*u_i + v_i*v_i;
		V_i = sqrt(V_i);

		eTot_i = rhoEL[i]*rho_i_inv;

		p_i = (GAMMA-1)*rho_i*(eTot_i - 0.5*V_i*V_i);

		c_i_2 = GAMMA*p_i/rho_i;  // speed of sound squared
		c_i = sqrt(c_i_2);

		// Total Enthalpy (constant outside domain as well)
		H_i = (p_i/rho_i)*(GAMMA/(GAMMA-1)) + 0.5*V_i*V_i;

		// Normal component of velocity 
		Vn_i = u_i*nL_x[i] + v_i*nL_y[i];

		R_i = Vn_i + 2.0*c_i/(GAMMA-1);

		// Solve for c_b

		double 	aQ, bQ, cQ, term1, term2, cM, cP, 
				c_b, V_b, M_b, T_b, p_b, rho_b, u_b, v_b, eTot_b;

		aQ = 1+2.0/(GAMMA-1);
		bQ = -2.0*R_i;
		cQ = 0.5*(GAMMA-1)*(R_i*R_i - 2.0*H_i);

		term1 = -bQ/(2.0*aQ);
		term2 = sqrt(bQ*bQ-4.0*aQ*cQ)/(2.0*aQ);

		cM = term1-term2;
		cP = term1+term2;

		// c = max(cM,cP)
		if (cM > cP){
			c_b = cM;
		} else{
			c_b = cP;
		}

		// Same reimann invariant of inner state is that of boundary state
		V_b = R_i - 2.0*c_b/(GAMMA-1);

		// Mach number of boundary flow.
		M_b = V_b/c_b;

		// Use isentropic relations to get static conditions now at the boundary
		T_b = T_Total/(1.0+0.5*(GAMMA-1)*M_b*M_b);
		p_b = p_Total*pow(T_b/T_Total, GAMMA/(GAMMA-1));

		rho_b = p_b/(DB.Rg*T_b);
		u_b = V_b*nL_x[i];
		v_b = V_b*nL_y[i];

		eTot_b = p_b/((GAMMA-1)*rho_b) + 0.5*(u_b*u_b + v_b*v_b);

		// Get boundary state (conservative variables)
		rhoB[i] = rho_b;
		rhouB[i] = rho_b*u_b;
		rhovB[i] = rho_b*v_b;
		rhoEB[i] = rho_b*eTot_b;

	}
	
	// Set up the flux vectors for the boundary state
	for (i=0; i<n; i++){
		
		W[0] = rhoB[i];
		W[1] = rhouB[i];
		W[2] = rhovB[i];
		W[3] = rhoEB[i];

		euler_flux_2D(W, F, G);

		for(j=0; j<4; j++){
			FB[j*n + i] = F[j];
			GB[j*n + i] = G[j];
		}
	}

	if(DB.Testing == 1){
		// Print the state at the wall:
		printf("Boundary State (TOTAL TP): \n");

		for (i=0; i<n; i++){
			printf("	Face Node : %d \n", i);
			printf("		[nx, ny] : [%f, %f] \n", nL_x[i], nL_y[i]);
			printf("		rhoL, rhoB : %f,  %f \n", rhoL[i], rhoB[i]);
			printf("		rhouL, rhouB : %f,  %f \n", rhouL[i], rhouB[i]);
			printf("		rhovL, rhovB : %f,  %f \n", rhovL[i], rhovB[i]);
			printf("		rhoEL, rhoEB : %f,  %f \n", rhoEL[i], rhoEB[i]);
		}

	}


}



