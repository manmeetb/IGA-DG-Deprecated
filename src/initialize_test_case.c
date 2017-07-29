
#include "initialize_test_case.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "exact_solutions.h"
#include "matrix_functions.h"
#include "S_VOLUME.h"
#include "S_DB.h"
#include "S_ELEMENT.h"


void initialize_test_case(void){

	/*
	Purpose:
		Setup the What vector for all volumes with the initial condition
	*/

	struct S_VOLUME *VOLUME;
	int i;
	double *XYZ_Sx, *XYZ_Sy;
	double WTemp[4], XYZTemp[2];
	double *W_Sol, *WHat, *ro, *ro_u, *ro_v, *ro_e;

	for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
		// Loop over all the volumes
		XYZ_Sx = &VOLUME->XYZ_S[0];
		XYZ_Sy = &VOLUME->XYZ_S[VOLUME->NvnS];

		W_Sol = malloc(VOLUME->NVar*VOLUME->NvnS* sizeof *W_Sol); //free

		ro = &W_Sol[0*VOLUME->NvnS];
		ro_u = &W_Sol[1*VOLUME->NvnS];
		ro_v = &W_Sol[2*VOLUME->NvnS];
		ro_e = &W_Sol[3*VOLUME->NvnS];

		for(i=0; i<VOLUME->NvnS; i++){
			XYZTemp[0] = XYZ_Sx[i];
			XYZTemp[1] = XYZ_Sy[i];

			// Set the initial state at the solution node:
			if (strstr(DB.TestType,"InviscidChannel")){
				uniform_solution_InternalSubsonic(XYZTemp, WTemp);
			} else if (strstr(DB.TestType,"PeriodicVortex")){
				exact_solution_IsentropicVortex(XYZTemp, WTemp);
			} else{
				printf("Case Not Implemented \n");
				exit(1);
			}

			ro[i] = WTemp[0];
			ro_u[i] = WTemp[1];
			ro_v[i] = WTemp[2];
			ro_e[i] = WTemp[3];

		}

		// Use Interpolation operator to find WHat (modal coefficients)
		mm_CNN(VOLUME->NvnG, VOLUME->NvnS, VOLUME->NVar, DB.ELEMENT->ChiInv_vS, W_Sol, VOLUME->What);

		free(W_Sol);

	}
}

