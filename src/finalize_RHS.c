#include "finalize_RHS.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "matrix_functions.h"

double finalize_RHS(void){

	/*
	Purpose:
		Compute the RHS for the volume using the RHS_FACE and RHS_VOL
		for each volume. Then, multiply by MInv to the RHS to get 
		the rate of change of each modal coefficient.
	*/

	int i;
	struct S_VOLUME *VOLUME;

	double maxRHS = 0;

	double *temp_RHS; // temporary RHS matrix (before MInv multiplied)

	double resVec[4];

	for(i=0; i<4; i++){
		resVec[i] = 0;
	}

	for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){

		temp_RHS = malloc(VOLUME->NvnG*VOLUME->NVar* sizeof *temp_RHS);  // free

		// Add both RHS contributions
		for(i=0; i<VOLUME->NvnG*VOLUME->NVar; i++){
			temp_RHS[i] = VOLUME->RHS_VOL[i] - VOLUME->RHS_FACE[i];	
			if(fabs(temp_RHS[i]) > maxRHS && i<VOLUME->NvnG){
				// Compute the LInf Norm of RHS for the density
				maxRHS = fabs(temp_RHS[i]);
			}
		}

		// Multiply MInv to get RHS
		mm_CNN(VOLUME->NvnG, VOLUME->NvnG, VOLUME->NVar, VOLUME->MInv, temp_RHS, VOLUME->RHS);

		free(temp_RHS);

	}

	return maxRHS;

}



