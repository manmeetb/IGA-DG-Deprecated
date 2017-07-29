
#include "setup_geom_factors.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "S_DB.h"
#include "S_ELEMENT.h"
#include "matrix_functions.h"


void setup_geom_factors(struct S_VOLUME *VOLUME){

	/*
	Purpose:
		Setup the metric terms and the Jacobian. For now, everything is
		collocated so the metric terms will be found at the integration
		points which are at the solution nodes.

	*/

	int i;
	double *GradChi_vS_xi, *GradChi_vS_eta;
	double *x_xi, *x_eta, *y_xi, *y_eta;
	double *XYZ_x, *XYZ_y;
	double *C_vS, *C_vS_11, *C_vS_12, *C_vS_21, *C_vS_22, *detJV_vS;

	XYZ_x = &VOLUME->XYZ[0];
	XYZ_y = &VOLUME->XYZ[VOLUME->NvnG];

	GradChi_vS_xi = DB.ELEMENT->GradChi_vS_xi;
	GradChi_vS_eta = DB.ELEMENT->GradChi_vS_eta;

	x_xi = malloc(VOLUME->NvnS*1* sizeof *x_xi);  // free
	x_eta = malloc(VOLUME->NvnS*1* sizeof *x_xi);  // free
	y_xi = malloc(VOLUME->NvnS*1* sizeof *x_xi);  // free
	y_eta = malloc(VOLUME->NvnS*1* sizeof *x_xi);  // free

	mm_CNN(VOLUME->NvnS, VOLUME->NvnG, 1, GradChi_vS_xi, XYZ_x, x_xi);
	mm_CNN(VOLUME->NvnS, VOLUME->NvnG, 1, GradChi_vS_xi, XYZ_y, y_xi);
	mm_CNN(VOLUME->NvnS, VOLUME->NvnG, 1, GradChi_vS_eta, XYZ_x, x_eta);
	mm_CNN(VOLUME->NvnS, VOLUME->NvnG, 1, GradChi_vS_eta, XYZ_y, y_eta);

	// Store the results in the Cofactor matrix and Jacobian vector

	C_vS_11 = &VOLUME->C_vS[0];
	C_vS_12 = &VOLUME->C_vS[1*VOLUME->NvnS];
	C_vS_21 = &VOLUME->C_vS[2*VOLUME->NvnS];
	C_vS_22 = &VOLUME->C_vS[3*VOLUME->NvnS];

	for(i=0; i<VOLUME->NvnS; i++){
		C_vS_11[i] = y_eta[i];
		C_vS_12[i] = -y_xi[i];
		C_vS_21[i] = -x_eta[i];
		C_vS_22[i] = x_xi[i];

		VOLUME->detJV_vS[i] = x_xi[i]*y_eta[i] - y_xi[i]*x_eta[i];
	}

	

	free(x_xi);
	free(x_eta);
	free(y_xi);
	free(y_eta);

}

