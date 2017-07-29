
#include "update_VOLUMEs.h"

#include <stdlib.h>
#include <stdio.h>

#include "S_VOLUME.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "matrix_functions.h"

void compute_inverse_mass(struct S_VOLUME *VOLUME){

	/*
	Purpose:
		Compute the inverse Mass matrix of the volume.
	*/

	int i, j, k, iBasisSolPt, jBasisSolPt, iSolPt, jSolPt, nSolPt;
	double *nodes_wi, *wBasis_i, *wBasis_j;
	double intVal; 

	// Mass and mass inverse matrix stored in column major form.
	// ordering of modal coefficients that it multiplies is all
	// eta = -1, ...
	double *M;
	
	nodes_wi = DB.ELEMENT->nodes_wi;

	wBasis_i = malloc(VOLUME->NvnS* sizeof *wBasis_i);  // free
	wBasis_j = malloc(VOLUME->NvnS* sizeof *wBasis_j);  // free

	M = VOLUME->MInv;  // keep (after inverse)

	k = 0;
	// Loop order = Fix column and go down all rows
	for(j=0; j<VOLUME->NvnG; j++){
		for(i=0; i<VOLUME->NvnG; i++){

			// i,j entry involves w_i and w_j

			// Fill the w_i and w_j values at the solution nodes. Need
			// To fill the ith and jth column of V matrix
			// w_i: 
			for(iBasisSolPt=0; iBasisSolPt < VOLUME->NvnS; iBasisSolPt++){
				wBasis_i[iBasisSolPt] = DB.ELEMENT->Chi_vS[i*DB.ELEMENT->NvnS + iBasisSolPt];
			}
			// w_j: 
			for(jBasisSolPt=0; jBasisSolPt < VOLUME->NvnS; jBasisSolPt++){
				wBasis_j[jBasisSolPt] = DB.ELEMENT->Chi_vS[j*DB.ELEMENT->NvnS + jBasisSolPt];
			}


			// Compute Integral using Gaussian Quadrature
			intVal = 0;
			nSolPt = 0;
			for(iSolPt=0; iSolPt<VOLUME->P+1; iSolPt++){
				for(jSolPt=0; jSolPt<VOLUME->P+1; jSolPt++){
					intVal = intVal + wBasis_i[nSolPt]*wBasis_j[nSolPt]*VOLUME->detJV_vS[nSolPt]
						*nodes_wi[iSolPt]*nodes_wi[jSolPt];
					nSolPt++;

				}
			}

			M[j*VOLUME->NvnG + i] = intVal;
			k++;
		}
	}

	mm_inv_d(VOLUME->NvnG, M);

	free(wBasis_i);
	free(wBasis_j);

}


void update_VOLUME_Ops(void){
	/*
	Purpose:
		Update the operators needed for each volume
	*/

	struct S_VOLUME    *VOLUME;

	int i,j;

	for (VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->update) {
			VOLUME->update = 0;
			compute_inverse_mass(VOLUME);

			if (DB.Testing == 1){

				printf("Inverse Mass Matrix : \n");

				for(i=0; i<VOLUME->NvnG; i++){
					for(j=0; j<VOLUME->NvnG; j++){
						printf(" %f ", VOLUME->MInv[j*VOLUME->NvnG + i]);
					}
					printf("\n");
				}
			}

		}
	}

}

