
#include "setup_operators.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "memory_constructors.h"
#include "cubature.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "bases.h"
#include "matrix_functions.h"


/*
 *  Purpose:
 *	Set up operators to be used throughout the code for the reference 
 *	elements.
 * 
 *	Notation:
 *
 *  Grad)Chi(Ref)(1)_(3)(4) : 
 *		(Grad)ient (optional) (Ref)erence (optional) Basis functions (Chi) of type
 *			(1) evaluated at (3) nodes of (4)
 *			(1/4): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *			(3): (v)olume, (f)ace, (e)dge
 *
 *	I_(1)(2)_(3)(4) : (I)nterpolation operator from (1) nodes of type (2) to (4) nodes of type (5)
 *		(1/3): (v)olume, (f)ace
 *		(2/4): (P)lotting, (G)eometry, (I)ntegration, (S)olution
 *	
 */


static void setup_Interpolation_operators(void){

	/*
	Purpose:
		Setup the interpolation operators for the reference element
		using the Chi operators which have been set up already. 

		NOTE: Interpolation is always from the nodal values, to their modal
		values first. Then, the modal values are multiplied to the basis evaluated
		at a set of points to get the nodal values at those points. 

			nodal_1 -> modal -> nodal_2

		ex : Therefore, Interpolation operators are of the form:
			I_vG_vS = Chi_vS * inv(Chi_vG)
	*/

	struct S_ELEMENT *ELEMENT;
	ELEMENT = DB.ELEMENT;

	// -------------------------------------------
	//			I_vG_vS and I_vS_vG
	// -------------------------------------------
	double *I_vG_vS, *I_vS_vG, *Chi_vS, *ChiInv_vG;

	// Interpolation Operator
	I_vG_vS = ELEMENT->I_vG_vS;

	// Chi Operators
	Chi_vS = ELEMENT->Chi_vS;
	ChiInv_vG = mm_inv_d_alloc(ELEMENT->NvnG, ELEMENT->Chi_vG);  // free

	// I_vG_vS:
	mm_CNN(ELEMENT->NvnS, ELEMENT->NvnG, ELEMENT->NvnG, Chi_vS, ChiInv_vG, I_vG_vS);	
	
	// I_vS_vG:
	mm_inv_d_secondInPlace(ELEMENT->NvnS, I_vG_vS, ELEMENT->I_vS_vG);

	free(ChiInv_vG);

	// -------------------------------------------
	//					I_vG_vP
	// -------------------------------------------
	double *I_vG_vP, *Chi_vP;

	// Interpolation Operator
	I_vG_vP = ELEMENT->I_vG_vP;

	// Chi Operators
	Chi_vP = ELEMENT->Chi_vP;
	ChiInv_vG = mm_inv_d_alloc(ELEMENT->NvnG, ELEMENT->Chi_vG);  // free

	// I_vG_vP:
	mm_CNN(ELEMENT->NvnP, ELEMENT->NvnG, ELEMENT->NvnG, Chi_vP, ChiInv_vG, I_vG_vP);	

	free(ChiInv_vG);

	// -------------------------------------------
	//					I_vS_fI
	// -------------------------------------------
	double *I_vS_fI, *Chi_fI, *ChiInv_vS;

	// Interpolation Operator
	I_vS_fI = ELEMENT->I_vS_fI;

	// Chi Operators
	Chi_fI = ELEMENT->Chi_fI;
	ChiInv_vS = ELEMENT->ChiInv_vS;

	// I_vS_fI:
	mm_CNN(ELEMENT->NfnI, ELEMENT->NvnG, ELEMENT->NvnG, Chi_fI, ChiInv_vS, I_vS_fI);	

}


static void setup_Chi_operators(void){
	/*
	Purpose:
		Set up the chi operators (basis function evaluated 
		at different points on the computational domain). These operators
		will then be used to make the interpolation operators and solve the
		flow
	*/

	struct S_ELEMENT *ELEMENT;
	ELEMENT = DB.ELEMENT;

	int basis_i, basis_j, basis_index, node_index,i,j;

	// -------------------------------------------
	//			Chi_vS and ChiInv_vS
	// -------------------------------------------
	double *Chi_vS, *ChiInv_vS, *XiEtaZeta_S_xi, *XiEtaZeta_S_eta;

	// Create the Chi_vS (Vandermonde matrix): 
	// 	V = Interpolation from geometry node points to volume nodes of type solution

	// There are NvnG basis functions (one for each geometry node point due to the mapping).
	// There are NvnS solution points. Therefore, Chi_vS is of size NvnS x NvnG. All matrices
	// are stored in column major form.

	Chi_vS = ELEMENT->Chi_vS;
	
	// The xi,eta values of the solution nodes on comp domain
	XiEtaZeta_S_xi = &(ELEMENT->XiEtaZeta_S[0]);
	XiEtaZeta_S_eta = &(ELEMENT->XiEtaZeta_S[ELEMENT->NvnS]);

	basis_index = 0;
	for(basis_j=0; basis_j<(ELEMENT->P+1); basis_j++){
		for(basis_i=0; basis_i<(ELEMENT->P+1); basis_i++){
			// Loop over all the basis functions (order is same as ordering
			// of geometry node points (all j=0, j=1, ...))

			for(node_index=0 ;node_index<ELEMENT->NvnS; node_index++){

				//Loop over the solution points on reference element. Order is all
				// j=0, j=1, ...

				// Store the value in V. node_index is the row and basis_index is 
				// the column of the entry of the matrix.
				if(strstr(DB.BasisType,"Polynomial")){
					Chi_vS[basis_index*ELEMENT->NvnS + node_index] = basis_TP_Lagrange_2D(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->XiEtaZeta_G, XiEtaZeta_S_xi[node_index], XiEtaZeta_S_eta[node_index]);
				}else if(strstr(DB.BasisType,"NURBS")){
					Chi_vS[basis_index*ELEMENT->NvnS + node_index] = basis_TP_NURBS_2D(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->xiVector, ELEMENT->etaVector, XiEtaZeta_S_xi[node_index], XiEtaZeta_S_eta[node_index]);
				} else{
					printf("UNRECOGNIZED BASIS \n");
					exit(1);
				}

			}

			basis_index++;
		}
	}

	ChiInv_vS = ELEMENT->ChiInv_vS;

	for(i=0; i<ELEMENT->NvnS*ELEMENT->NvnG; i++){
		ChiInv_vS[i] = ELEMENT->Chi_vS[i];
	}

	mm_inv_d(ELEMENT->NvnG, ChiInv_vS);

	// -------------------------------------------
	//		GradChi_vS_xi and GradChi_vS_eta
	// -------------------------------------------

	double *GradChi_vS_xi, *GradChi_vS_eta;
	double *grad;

	GradChi_vS_xi = ELEMENT->GradChi_vS_xi;
	GradChi_vS_eta = ELEMENT->GradChi_vS_eta;

	// The xi,eta values of the solution nodes on comp domain. Already
	// declared when setting the previous operators (Chi_vS)
	XiEtaZeta_S_xi = &(ELEMENT->XiEtaZeta_S[0]);
	XiEtaZeta_S_eta = &(ELEMENT->XiEtaZeta_S[ELEMENT->NvnS]);

	basis_index = 0;
	for(basis_j=0; basis_j<(ELEMENT->P+1); basis_j++){
		for(basis_i=0; basis_i<(ELEMENT->P+1); basis_i++){
			// Loop over all the basis functions (order is same as ordering
			// of geometry node points (all j=0, j=1, ...))

			for(node_index=0; node_index<ELEMENT->NvnS; node_index++){

				//Loop over the solution points on reference element. Order is all
				// j=0, j=1, ...

				// grad is allocated in function
				if(strstr(DB.BasisType,"Polynomial")){
					grad = basis_TP_Lagrange_2D_Grad(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->XiEtaZeta_G, XiEtaZeta_S_xi[node_index], XiEtaZeta_S_eta[node_index]);
				}else if(strstr(DB.BasisType,"NURBS")){
					grad = basis_TP_NURBS_2D_Grad(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->xiVector, ELEMENT->etaVector, XiEtaZeta_S_xi[node_index], XiEtaZeta_S_eta[node_index]);
				} else{
					printf("UNRECOGNIZED BASIS \n");
					exit(1);
				}

				GradChi_vS_xi[basis_index*ELEMENT->NvnS + node_index] = grad[0];
				GradChi_vS_eta[basis_index*ELEMENT->NvnS + node_index] = grad[1];

				free(grad);

			}

			basis_index++;
		}
	}

	// -------------------------------------------
	//		GradChi_fI_xi and GradChi_fI_eta
	// -------------------------------------------

	// Compute gradient at the face integration nodes (for finding
	// the metric terms)
	double *GradChi_fI_xi, *GradChi_fI_eta, *XiEtaZeta_F_xi, *XiEtaZeta_F_eta;

	GradChi_fI_xi = ELEMENT->GradChi_fI_xi;
	GradChi_fI_eta = ELEMENT->GradChi_fI_eta;

	// Xi Eta values at face integration nodes on reference domain
	XiEtaZeta_F_xi = &(ELEMENT->XiEtaZeta_F[0]);
	XiEtaZeta_F_eta = &(ELEMENT->XiEtaZeta_F[ELEMENT->NfnI]);

	basis_index = 0;
	for(basis_j=0; basis_j<(ELEMENT->P+1); basis_j++){
		for(basis_i=0; basis_i<(ELEMENT->P+1); basis_i++){
			// Loop over all the basis functions (order is same as ordering
			// of geometry node points (all j=0, j=1, ...))

			for(node_index=0; node_index<ELEMENT->NfnI; node_index++){

				//Loop over the face integration points on reference element.

				// grad is allocated in function
				if(strstr(DB.BasisType,"Polynomial")){
					grad = basis_TP_Lagrange_2D_Grad(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->XiEtaZeta_G, XiEtaZeta_F_xi[node_index], XiEtaZeta_F_eta[node_index]);
				}else if(strstr(DB.BasisType,"NURBS")){
					grad = basis_TP_NURBS_2D_Grad(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->xiVector, ELEMENT->etaVector, XiEtaZeta_F_xi[node_index], XiEtaZeta_F_eta[node_index]);
				} else{
					printf("UNRECOGNIZED BASIS \n");
					exit(1);
				}

				GradChi_fI_xi[basis_index*ELEMENT->NfnI + node_index] = grad[0];
				GradChi_fI_eta[basis_index*ELEMENT->NfnI + node_index] = grad[1];

				free(grad);

			}

			basis_index++;
		}
	}



	// -------------------------------------------
	//			Chi_vG (not used in NURBS)
	// -------------------------------------------

	double *Chi_vG, *XiEtaZeta_G_xi, *XiEtaZeta_G_eta;

	Chi_vG = ELEMENT->Chi_vG;

	// The xi,eta values of the geometry nodes on comp domain
	XiEtaZeta_G_xi = &(ELEMENT->XiEtaZeta_G[0]);
	XiEtaZeta_G_eta = &(ELEMENT->XiEtaZeta_G[ELEMENT->NvnG]);

	basis_index = 0;
	for(basis_j=0; basis_j<(ELEMENT->P+1); basis_j++){
		for(basis_i=0; basis_i<(ELEMENT->P+1); basis_i++){
			// Loop over all the basis functions (order is same as ordering
			// of geometry node points (all j=0, j=1, ...))

			for(node_index=0; node_index<ELEMENT->NvnG; node_index++){

				//Loop over the solution points on reference element. Order is all
				// j=0, j=1, ...

				// Store the value in V. node_index is the row and basis_index is 
				// the column of the entry of the matrix.
				Chi_vG[basis_index*ELEMENT->NvnG + node_index] = basis_TP_Lagrange_2D(ELEMENT->P, basis_i, basis_j, 
					ELEMENT->XiEtaZeta_G, XiEtaZeta_G_xi[node_index], XiEtaZeta_G_eta[node_index]);

			}

			basis_index++;
		}
	}

	// -------------------------------------------
	//					Chi_vP
	// -------------------------------------------

	double *Chi_vP, *XiEtaZeta_P_xi, *XiEtaZeta_P_eta;

	Chi_vP = ELEMENT->Chi_vP;

	// The xi,eta values of the geometry nodes on comp domain
	XiEtaZeta_P_xi = &(ELEMENT->XiEtaZeta_P[0]);
	XiEtaZeta_P_eta = &(ELEMENT->XiEtaZeta_P[ELEMENT->NvnP]);

	basis_index = 0;
	for(basis_j=0; basis_j<(ELEMENT->P+1); basis_j++){
		for(basis_i=0; basis_i<(ELEMENT->P+1); basis_i++){
			// Loop over all the basis functions (order is same as ordering
			// of geometry node points (all j=0, j=1, ...))

			for(node_index=0; node_index<ELEMENT->NvnP; node_index++){

				//Loop over the solution points on reference element. Order is all
				// j=0, j=1, ...

				// Store the value in V. node_index is the row and basis_index is 
				// the column of the entry of the matrix.
				if(strstr(DB.BasisType,"Polynomial")){
					Chi_vP[basis_index*ELEMENT->NvnP + node_index] = basis_TP_Lagrange_2D(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->XiEtaZeta_G, XiEtaZeta_P_xi[node_index], XiEtaZeta_P_eta[node_index]);
				}else if(strstr(DB.BasisType,"NURBS")){
					Chi_vP[basis_index*ELEMENT->NvnP + node_index] = basis_TP_NURBS_2D(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->xiVector, ELEMENT->etaVector, XiEtaZeta_P_xi[node_index], XiEtaZeta_P_eta[node_index]);
				} else{
					printf("UNRECOGNIZED BASIS \n");
					exit(1);
				}

			}

			basis_index++;
		}
	}

	// -------------------------------------------
	//					Chi_fI
	// -------------------------------------------
	double *Chi_fI;

	Chi_fI = ELEMENT->Chi_fI;

	// The xi,eta values of the geometry nodes on comp domain
	XiEtaZeta_F_xi = &(ELEMENT->XiEtaZeta_F[0]);
	XiEtaZeta_F_eta = &(ELEMENT->XiEtaZeta_F[ELEMENT->NfnI]);

	basis_index = 0;
	for(basis_j=0; basis_j<(ELEMENT->P+1); basis_j++){
		for(basis_i=0; basis_i<(ELEMENT->P+1); basis_i++){
			// Loop over all the basis functions (order is same as ordering
			// of geometry node points (all j=0, j=1, ...))

			for(node_index=0; node_index<ELEMENT->NfnI; node_index++){

				//Loop over the solution points on reference element. Order is all
				// j=0, j=1, ...

				// Store the value in V. node_index is the row and basis_index is 
				// the column of the entry of the matrix.
				if(strstr(DB.BasisType,"Polynomial")){
					Chi_fI[basis_index*ELEMENT->NfnI + node_index] = basis_TP_Lagrange_2D(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->XiEtaZeta_G, XiEtaZeta_F_xi[node_index], XiEtaZeta_F_eta[node_index]);
				} else if(strstr(DB.BasisType,"NURBS")){
					Chi_fI[basis_index*ELEMENT->NfnI + node_index] = basis_TP_NURBS_2D(ELEMENT->P, basis_i, basis_j, 
						ELEMENT->xiVector, ELEMENT->etaVector, XiEtaZeta_F_xi[node_index], XiEtaZeta_F_eta[node_index]);
				} else{
					printf("UNRECOGNIZED BASIS \n");
					exit(1);
				}

			}

			basis_index++;
		}
	}

}

static void setup_reference_element(){
	
	/*
	Purpose:
		Set up the reference element struct (S_ELEMENT)
	*/

	int i,j,k;
	struct S_ELEMENT *ELEMENT;
	ELEMENT = DB.ELEMENT;
	double 	*XiEtaZeta_S, // Reference Element Nodes at Solution points
			*XiEtaZeta_S_x, *XiEtaZeta_S_y,
			*XiEtaZeta_G, // Reference Element Nodes at Geometry points
			*XiEtaZeta_G_x, *XiEtaZeta_G_y,
			*XiEtaZeta_P, // Reference Element Nodes at Plotting points
			*XiEtaZeta_P_x, *XiEtaZeta_P_y,
			*XiEtaZeta_F, // Reference Element Nodes at Face Integration points
			*XiEtaZeta_F_x, *XiEtaZeta_F_y;


	double xi, eta, deltaXi_G, deltaEta_G, deltaXi_P, deltaEta_P;

	// -----------------------------------
	//			Properties
	// -----------------------------------

	ELEMENT->P = DB.P;
	ELEMENT->d = DB.d;

	// Get number of geometry and solution nodes
	ELEMENT->NvnG = (ELEMENT->P+1);
	for(i=0; i<ELEMENT->d-1; i++){
		ELEMENT->NvnG = (ELEMENT->NvnG)*(ELEMENT->P+1);
	}
	ELEMENT->NvnS = ELEMENT->NvnG;
	ELEMENT->NvnP = ELEMENT->NvnG;
	ELEMENT->NfnI = 4*(ELEMENT->P+1);

	// -----------------------------------
	//			Solution Points
	// -----------------------------------

	double *Cubature_xi, *Cubature_wi;
	Cubature_xi = malloc((ELEMENT->P+1)* sizeof *Cubature_xi);
	Cubature_wi = malloc((ELEMENT->P+1)* sizeof *Cubature_wi);

	cubature_literature(ELEMENT->P, DB.NodeType, Cubature_xi, Cubature_wi);

	// Compute Solution Point Locations on Ref. Elem. This is a matrix of 
	// size NvnS x d stored in column major form.
	XiEtaZeta_S = malloc(ELEMENT->NvnS*ELEMENT->d* sizeof *XiEtaZeta_S);

	XiEtaZeta_S_x = &XiEtaZeta_S[0];
	XiEtaZeta_S_y = &XiEtaZeta_S[ELEMENT->NvnS];

	// Ordering of reference solution nodes is all eta = -1, then increase eta
	k = 0;
	for(j=0; j<(ELEMENT->P+1); j++){
		for(i=0; i<(ELEMENT->P+1); i++){

			xi = Cubature_xi[i];
			eta = Cubature_xi[j];

			XiEtaZeta_S_x[k] = xi;
			XiEtaZeta_S_y[k] = eta;

			k++;
		}
	}

	// -----------------------------------
	//	Geometry Points (Used in polynomial basis)
	// -----------------------------------

	// Compute Geometry Node locations on reference element. Geometry
	// node points are equally spaced out on the element. The ordering of
	// the points is all eta = -1,...

	XiEtaZeta_G = malloc(ELEMENT->NvnG*ELEMENT->d* sizeof *XiEtaZeta_G);
	XiEtaZeta_G_x = &XiEtaZeta_G[0];
	XiEtaZeta_G_y = &XiEtaZeta_G[ELEMENT->NvnG];

	deltaXi_G = (2.)/((double)(ELEMENT->P));
	deltaEta_G = (2.)/((double)(ELEMENT->P));

	k = 0;
	for(j=0; j<(ELEMENT->P+1); j++){
		for(i=0; i<(ELEMENT->P+1); i++){

			xi = -1. + i*deltaXi_G;
			eta = -1. + j*deltaEta_G;

			XiEtaZeta_G_x[k] = xi;
			XiEtaZeta_G_y[k] = eta;

			k++;
		}
	}

	// -----------------------------------
	//			Knot Vectors
	// -----------------------------------

	// Even if the mesh is not using a NURBS basis, set up the knot vector if it 
	// were being used based on the mesh. Now, to make the NURBS work on the
	// reference element, simply set the limits to be at -1 to 1 with the correct
	// multiplicity. Note that widen the knot vector range by a small amount
	// in the case we are evaluating points exactly at the limits of the element

	double *xiVector, *etaVector;
	xiVector = malloc((ELEMENT->P+1)*2* sizeof *xiVector);
	etaVector = malloc((ELEMENT->P+1)*2* sizeof *etaVector);

	for(i=0; i<ELEMENT->P+1; i++){
		xiVector[i] = -1 - 1E-14;
		etaVector[i] = -1 - 1E-14;
	}

	for(i=0; i<ELEMENT->P+1; i++){
		xiVector[i+ELEMENT->P+1] = 1 + 1E-14;
		etaVector[i+ELEMENT->P+1] = 1 + 1E-14;
	}



	// -----------------------------------
	//			Plotting Points
	// -----------------------------------

	// Compute the plotting point locations on the reference element. There
	// will be P+1 plotting points along each coordinate direction spaced equally.

	XiEtaZeta_P = malloc(ELEMENT->NvnP*ELEMENT->d* sizeof *XiEtaZeta_P);
	XiEtaZeta_P_x = &XiEtaZeta_P[0];
	XiEtaZeta_P_y = &XiEtaZeta_P[ELEMENT->NvnP];

	deltaXi_P = (2.)/((double)(ELEMENT->P));
	deltaEta_P = (2.)/((double)(ELEMENT->P));

	k = 0;
	for(j=0; j<(ELEMENT->P+1); j++){
		for(i=0; i<(ELEMENT->P+1); i++){

			xi = -1. + i*deltaXi_P;
			eta = -1. + j*deltaEta_P;

			XiEtaZeta_P_x[k] = xi;
			XiEtaZeta_P_y[k] = eta;

			k++;
		}
	}

	// -----------------------------------
	//		Face Integration Points
	// -----------------------------------

	// Compute the face integration point locations. Order of the nodes
	// in the array is counterclockwise starting from (xi,eta) = (-1, -1).
	// Note, cubature array is in increasing order.

	XiEtaZeta_F = malloc(ELEMENT->NfnI*ELEMENT->d* sizeof *XiEtaZeta_F);
	XiEtaZeta_F_x = &XiEtaZeta_F[0];
	XiEtaZeta_F_y = &XiEtaZeta_F[ELEMENT->NfnI];
	 
	int faceNodeCumul, face, faceNode;

	faceNodeCumul = 0;
	for(face = 0; face<4; face++){
		for(faceNode = 0; faceNode < ELEMENT->P+1; faceNode++){

			if(face == 0){
				xi = Cubature_xi[faceNode];
				eta = -1;

			} else if(face == 1){
				xi = 1;
				eta = Cubature_xi[faceNode];

			} else if(face == 2){
				xi = -Cubature_xi[faceNode];
				eta = 1;

			} else{
				xi = -1;
				eta = -Cubature_xi[faceNode];

			}

			XiEtaZeta_F_x[faceNodeCumul] = xi;
			XiEtaZeta_F_y[faceNodeCumul] = eta;
			faceNodeCumul++;
		}
	}

	// -----------------------------------
	//				Testing
	// -----------------------------------

	if (DB.Testing == 1 || DB.Testing == 4){
		// Print the coordinates of the solution nodes on the reference element
		printf("Geometry: \n");
		for(k=0; k<ELEMENT->NvnG; k++){
			printf("(xi,eta) = (%f, %f) \n", 	XiEtaZeta_G[k], 
												XiEtaZeta_G[ELEMENT->NvnG+k]);
		}

		printf("Solution: \n");
		for(k=0; k<ELEMENT->NvnS; k++){
			printf("(xi,eta) = (%f, %f) \n", 	XiEtaZeta_S[k], 
												XiEtaZeta_S[ELEMENT->NvnS+k]);
		}

		printf("Plotting: \n");
		for(k=0; k<ELEMENT->NvnP; k++){
			printf("(xi,eta) = (%f, %f) \n", 	XiEtaZeta_P[k], 
												XiEtaZeta_P[ELEMENT->NvnP+k]);
		}

	}

	// Allocate arrays for the different operators to be used

	// Chi Operators
	ELEMENT->Chi_vS = malloc(ELEMENT->NvnS*ELEMENT->NvnG* sizeof *ELEMENT->Chi_vS);
	ELEMENT->ChiInv_vS = malloc(ELEMENT->NvnS*ELEMENT->NvnG* sizeof *ELEMENT->ChiInv_vS);
	ELEMENT->Chi_vG = malloc(ELEMENT->NvnG*ELEMENT->NvnG* sizeof *ELEMENT->Chi_vG);
	ELEMENT->Chi_vP = malloc(ELEMENT->NvnP*ELEMENT->NvnG* sizeof *ELEMENT->Chi_vP);
	ELEMENT->Chi_fI = malloc(ELEMENT->NvnG*ELEMENT->NfnI* sizeof *ELEMENT->Chi_fI);

	// Derivative Operators
	ELEMENT->GradChi_vS_xi = malloc(ELEMENT->NvnS*ELEMENT->NvnG* sizeof *ELEMENT->GradChi_vS_xi);
	ELEMENT->GradChi_vS_eta = malloc(ELEMENT->NvnS*ELEMENT->NvnG* sizeof *ELEMENT->GradChi_vS_eta);	
	ELEMENT->GradChi_fI_xi = malloc(ELEMENT->NfnI*ELEMENT->NvnG* sizeof *ELEMENT->GradChi_fI_xi);
	ELEMENT->GradChi_fI_eta = malloc(ELEMENT->NfnI*ELEMENT->NvnG* sizeof *ELEMENT->GradChi_fI_eta);

	// Interpolation Operators
	ELEMENT->I_vG_vP = malloc(ELEMENT->NvnP*ELEMENT->NvnG* sizeof *ELEMENT->I_vG_vP);
	ELEMENT->I_vG_vS = malloc(ELEMENT->NvnG*ELEMENT->NvnS* sizeof *ELEMENT->I_vG_vS);
	ELEMENT->I_vS_vG = malloc(ELEMENT->NvnS*ELEMENT->NvnG* sizeof *ELEMENT->I_vS_vG);
	ELEMENT->I_vS_fI = malloc(ELEMENT->NvnS*ELEMENT->NfnI* sizeof *ELEMENT->I_vS_fI);

	// Save reference element node locations and other cubature data
	ELEMENT->nodes_xi = Cubature_xi;
	ELEMENT->nodes_wi = Cubature_wi;

	ELEMENT->XiEtaZeta_S = XiEtaZeta_S;
	ELEMENT->XiEtaZeta_G = XiEtaZeta_G;
	ELEMENT->XiEtaZeta_P = XiEtaZeta_P;
	ELEMENT->XiEtaZeta_F = XiEtaZeta_F;

	ELEMENT->xiVector = xiVector;
	ELEMENT->etaVector = etaVector;

}


void setup_operators(void){

	printf("Setup Operators \n");

	struct S_ELEMENT *ELEMENT;
	ELEMENT = New_ELEMENT();
	DB.ELEMENT = ELEMENT;

	setup_reference_element();

	setup_Chi_operators();

	setup_Interpolation_operators();

	if(DB.Testing == 1 || DB.Testing == 3){
		//Print V (Vandermonde matrix)

		int r, c;

		printf("Chi_vS : Vandermonde \n");
		for(r=0; r<ELEMENT->NvnS; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %.15f ", ELEMENT->Chi_vS[c*ELEMENT->NvnS + r]);

			}
			printf("\n");
		}

		printf("ChiInv_vS : Inverse Vandermonde \n");

		for(r=0; r<ELEMENT->NvnS; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %.15f ", ELEMENT->ChiInv_vS[c*ELEMENT->NvnS + r]);

			}
			printf("\n");
		}

		printf("GradChi_vS_xi : Gradient Computational Domain \n");
		for(r=0; r<ELEMENT->NvnS; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->GradChi_vS_xi[c*ELEMENT->NvnS + r]);

			}
			printf("\n");
		}

		printf("GradChi_vS_eta : Gradient Computational Domain \n");
		for(r=0; r<ELEMENT->NvnS; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->GradChi_vS_eta[c*ELEMENT->NvnS + r]);

			}
			printf("\n");
		}

		printf("Chi_vG :  \n");
		for(r=0; r<ELEMENT->NvnG; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->Chi_vG[c*ELEMENT->NvnG + r]);

			}
			printf("\n");
		}

		printf("Chi_vP :  \n");
		for(r=0; r<ELEMENT->NvnP; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->Chi_vP[c*ELEMENT->NvnP + r]);

			}
			printf("\n");
		}

		printf("Chi_fI :  \n");
		for(r=0; r<ELEMENT->NfnI; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->Chi_fI[c*ELEMENT->NfnI + r]);
			}
			printf("\n");
		}

		printf("I_vG_vS :  \n");
		for(r=0; r<ELEMENT->NvnS; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->I_vG_vS[c*ELEMENT->NvnS + r]);
			}
			printf("\n");
		}

		printf("I_vS_vG :  \n");
		for(r=0; r<ELEMENT->NvnG; r++){
			for(c=0; c<ELEMENT->NvnS; c++){
				printf(" %f ", ELEMENT->I_vS_vG[c*ELEMENT->NvnG + r]);
			}
			printf("\n");
		}

		printf("I_vG_vP :  \n");
		for(r=0; r<ELEMENT->NvnP; r++){
			for(c=0; c<ELEMENT->NvnG; c++){
				printf(" %f ", ELEMENT->I_vG_vP[c*ELEMENT->NvnP + r]);
			}
			printf("\n");
		}

		printf("I_vS_fI :  \n");
		for(r=0; r<ELEMENT->NfnI; r++){
			for(c=0; c<ELEMENT->NvnS; c++){
				printf(" %f ", ELEMENT->I_vS_fI[c*ELEMENT->NvnS + r]);
			}
			printf("\n");
		}

	}

}




