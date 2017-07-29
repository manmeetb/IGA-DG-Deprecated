
#include "explicit_VOLUME_info.h"

#include <stdlib.h>
#include <stdio.h>

#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_ELEMENT.h"
#include "matrix_functions.h"
#include "euler_flux.h"


void explicit_VOLUME_info(void){

	struct S_VOLUME *VOLUME;

	int iSolPt, jSolPt, nSolPt;
	int iBasis, jBasis, nBasis;

	double intVal;
	double *W_Sol;  // Solution at each solution node as a matrix
	double *nodes_wi;

	// The flux vectors (each component). Will be a matrix in column 
	// major form with F or G for each of the NVar equations at a solution 
	// point (row in matrix) in the different columns.
	double *Flux_F, *Flux_G, Flux_F_Eq, Flux_G_Eq;  

	double *RHS_VOL, *C_vS_y_eta, *min_C_vS_y_xi, *min_C_vS_x_eta, *C_vS_x_xi;
	double *gradChi_xi, *gradChi_eta;

	int iVar;

	W_Sol = malloc(DB.VOLUME_HEAD->NvnS*DB.VOLUME_HEAD->NVar* sizeof *W_Sol);  // free
	Flux_F = malloc(DB.VOLUME_HEAD->NvnS*DB.VOLUME_HEAD->NVar* sizeof *Flux_F);  // free
	Flux_G = malloc(DB.VOLUME_HEAD->NvnS*DB.VOLUME_HEAD->NVar* sizeof *Flux_G);  // free

	nodes_wi = DB.ELEMENT->nodes_wi;

	for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){

		RHS_VOL = VOLUME->RHS_VOL;

		// Get W_Sol (approximate solution at solution nodes)
		// Ordering is always all bottom row (j=1) of nodes in the 
		// volume and in positive j direction. This is a matrix (NVar cols)

		mm_CNN(VOLUME->NvnS, VOLUME->NvnG, VOLUME->NVar, 
			DB.ELEMENT->Chi_vS, VOLUME->What, W_Sol);

		// Compute Flux vector at each solution node using solution
		euler_flux_2D_matrix(W_Sol, Flux_F, Flux_G, VOLUME->NvnS);

		// metric terms
		C_vS_y_eta = &VOLUME->C_vS[0];  
		min_C_vS_y_xi = &VOLUME->C_vS[VOLUME->NvnS];  
		min_C_vS_x_eta = &VOLUME->C_vS[2*VOLUME->NvnS];  
		C_vS_x_xi = &VOLUME->C_vS[3*VOLUME->NvnS];  

		// Build the volume contribution vector using
		// Gaussian Quadrature
		for(iVar=0; iVar < VOLUME->NVar; iVar++){
			nBasis = 0;
			for(jBasis=0; jBasis<VOLUME->P+1; jBasis++){
				for(iBasis=0; iBasis<VOLUME->P+1; iBasis++){

					// Gradient of the i,j th basis function at all the solution 
					// nodes
					gradChi_xi = &DB.ELEMENT->GradChi_vS_xi[nBasis*VOLUME->NvnS];
					gradChi_eta = &DB.ELEMENT->GradChi_vS_eta[nBasis*VOLUME->NvnS];

					intVal = 0;
					nSolPt = 0;
					for(jSolPt=0; jSolPt<VOLUME->P+1; jSolPt++){
						for(iSolPt=0; iSolPt<VOLUME->P+1; iSolPt++){
							// Loop over solution points all j=1 (first row) ...

							// Find flux at this solution node
							Flux_F_Eq = Flux_F[iVar*VOLUME->NvnS + nSolPt];
							Flux_G_Eq = Flux_G[iVar*VOLUME->NvnS + nSolPt];

							intVal = intVal	+ 
								(
								(C_vS_y_eta[nSolPt]*Flux_F_Eq + min_C_vS_x_eta[nSolPt]*Flux_G_Eq)*gradChi_xi[nSolPt]
							+ 	(min_C_vS_y_xi[nSolPt]*Flux_F_Eq + C_vS_x_xi[nSolPt]*Flux_G_Eq)*gradChi_eta[nSolPt]
								)*nodes_wi[iSolPt]*nodes_wi[jSolPt];
							
							nSolPt++;
						}
					}

					RHS_VOL[iVar*VOLUME->NvnG + nBasis] = intVal;

					nBasis++;
				}
			}
		}

	}

	free(W_Sol);
	free(Flux_F);
	free(Flux_G);

}

