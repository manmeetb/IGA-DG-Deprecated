
#include "compute_errors.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_ELEMENT.h"

#include "exact_solutions.h"
#include "matrix_functions.h"

void compute_errors_global(){

	/*
	Purpose:
		Compute the Error of the flow with the exact solution
	*/

	struct S_VOLUME *VOLUME;
	int i, j, r, c;
	double *XVals, *YVals, *roVec, *ro_uVec, *ro_vVec, *eVec, x, y, 
			ro, ro_u, ro_v, e;
	double *W_Sol;
	double total_dof;

	double error_L2[4];
	for(i=0; i<4; i++){
		error_L2[i] = 0;
	}

	double W_exact[4], XYZ_vec[2];

	total_dof = 0;

	for (VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next) {

		total_dof = total_dof + VOLUME->NvnS;

		// Solution node locations. Recall that values are stored in 
		// column major form.	
		XVals = &(VOLUME->XYZ_S[0]);
		YVals = &(VOLUME->XYZ_S[1*VOLUME->NvnS]);

		W_Sol = malloc(VOLUME->NvnS*VOLUME->NVar* sizeof *W_Sol);

		mm_CNN(VOLUME->NvnS, VOLUME->NvnG, VOLUME->NVar, DB.ELEMENT->Chi_vS, VOLUME->What, W_Sol);

		// Loop over all the solution nodes on this volume and get the difference
		i = 0;
		for(r = 0; r<VOLUME->P+1; r++){
			for(c=0; c<VOLUME->P+1; c++){
				x = XVals[i];
				y = YVals[i];

				XYZ_vec[0] = x;
				XYZ_vec[1] = y;

				exact_solution_IsentropicVortex(XYZ_vec, W_exact);

				for(j=0; j<4; j++){
					error_L2[j] = error_L2[j] + 
						(W_Sol[j*VOLUME->NvnS + i] - W_exact[j])*(W_Sol[j*VOLUME->NvnS + i] - W_exact[j])*
						DB.ELEMENT->nodes_wi[r]*DB.ELEMENT->nodes_wi[c]*VOLUME->detJV_vS[i];
				}
				i++;
			}
		}

		free(W_Sol);
	}

	printf("L2 ERROR : \n");
	printf("	numDOF : %f \n", total_dof);
	// Normalize the Error and Print the results:
	for(i=0; i<4; i++){
		error_L2[i] = sqrt(error_L2[i]/(16.*16.));

		printf("	w : %d  -> %.14e \n", i, error_L2[i]);

	}


}





