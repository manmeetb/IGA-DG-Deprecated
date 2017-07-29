
#include "setup_normals.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "S_ELEMENT.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "matrix_functions.h"

void setup_normals(struct S_FACE *FACE){

	/*
	Purpose:
		Set up the physical normal vectors at each face integration point.
		Outward normals will be set up for the VIn volume.
		To do this the following will be done
			1) Compute the metric terms at the face integration nodes. 
				Store these values at in a structure.
			2) Use the metric terms at the face integration nodes to compute
				the normal vector to each point. The convention used is that 
				all vertices will be ordered in counterclockwise order in the
				mapping so the rotation matrix to be used for each tangent vector
				is simply a function of which mapped face is being worked with.
	*/

	int fin, i;
	double *n_fI, *n_fI_x, *n_fI_y, norm_n;
	double y_eta, min_y_xi, min_x_eta, x_xi; 

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VIn;

	// Metric terms at the face integration nodes for this face
	double *C_fI, *xy_xi, *xy_eta;

	fin = FACE->fin;
	VIn = FACE->VIn;

	C_fI = FACE->C_fI;
	xy_xi = malloc(VIn->NfnI*2* sizeof *xy_xi);  // free
	xy_eta = malloc(VIn->NfnI*2* sizeof *xy_eta);  // free

	mm_CNN(DB.ELEMENT->NfnI, DB.ELEMENT->NvnG, 2, DB.ELEMENT->GradChi_fI_xi, VIn->XYZ, xy_xi);
	mm_CNN(DB.ELEMENT->NfnI, DB.ELEMENT->NvnG, 2, DB.ELEMENT->GradChi_fI_eta, VIn->XYZ, xy_eta);

	// Store the values of the metric for the current face
	for(i=0; i<FACE->P+1; i++){
		// There are P+1 face integration nodes to loop over. Use fIn to 
		// find the correct sequence of nodes

		// Store order to be C matrix in row major form
		C_fI[i] = xy_eta[VIn->NfnI + fin*(VIn->P+1) + i];  // C11 = y_eta
		C_fI[i+1*(FACE->P+1)] = -1*xy_xi[VIn->NfnI + fin*(VIn->P+1) + i];  // C12 = - y_xi
		C_fI[i+2*(FACE->P+1)] = -1*xy_eta[fin*(VIn->P+1) + i]; // C21 = - x_eta
		C_fI[i+3*(FACE->P+1)] = xy_xi[fin*(VIn->P+1) + i]; // C22 = x_xi

	}

	// Compute outward normal vector (to VIn) using metric terms
	n_fI = FACE->n_fI;
	n_fI_x = &n_fI[0];
	n_fI_y = &n_fI[(FACE->P+1)];

	for(i=0; i<(FACE->P+1); i++){

		y_eta = C_fI[i];
		min_y_xi = C_fI[i+1*(FACE->P+1)];
		min_x_eta = C_fI[i+2*(FACE->P+1)];
		x_xi = C_fI[i+3*(FACE->P+1)];

		if(fin == 0){
			// f = 0 face : Constant eta line (eta = -1)
			n_fI_x[i] = -1*min_y_xi;
			n_fI_y[i] = -1*x_xi;

		}else if(fin == 1){
			// f = 1 face : Constant xi line (xi = 1)
			n_fI_x[i] = y_eta;
			n_fI_y[i] = min_x_eta;


		} else if(fin == 2){
			// f = 2 face : Constant eta line (eta = 1)
			n_fI_x[i] = min_y_xi;
			n_fI_y[i] = x_xi;


		} else{
			// f = 3 face : Constant xi line (xi = -1)
			n_fI_x[i] = -1*y_eta;
			n_fI_y[i] = -1*min_x_eta;
		}

		// Normalize Normal Vector
		norm_n = n_fI_x[i]*n_fI_x[i] + n_fI_y[i]*n_fI_y[i];
		norm_n = sqrt(norm_n);

		n_fI_x[i] = n_fI_x[i]/norm_n;
		n_fI_y[i] = n_fI_y[i]/norm_n;

	}

	free(xy_xi);
	free(xy_eta);

}


