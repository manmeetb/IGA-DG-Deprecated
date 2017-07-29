
#include "explicit_FACE_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "S_DB.h"
#include "S_FACE.h"
#include "S_VOLUME.h"
#include "S_ELEMENT.h"
#include "matrix_functions.h"
#include "euler_flux.h"
#include "fluxes_inviscid.h"
#include "boundary_conditions.h"
#include "Parameters.h"

static void flipMatrixOrder(double *A, int rowNum, int colNum){
	/*
	Purpose:
		Flip the row order of a matrix A
	*/

	int r, c, i;

	double *TempMatrix;

	TempMatrix = malloc(rowNum*colNum* sizeof *TempMatrix);  // free

	for(r=0; r<rowNum; r++){
		// Loop over the rows

		for(c=0; c<colNum; c++){
			// Loop over the variables

			TempMatrix[c*rowNum + r] = A[c*rowNum + rowNum-r-1];
		}
	}

	for(i=0; i<rowNum*colNum; i++){
		A[i] = TempMatrix[i];
	}

	free(TempMatrix);

}


void explicit_FACE_info(void){

	/*
	Purpose:
		Compute RHS_FACE for all volumes
	*/

	// General Variables:
	int ifI, i, j, iVar, i_face;
	int iBasis, jBasis, nBasis, nodeF_index, numfI; 

	// Temporary vectors for finding flux at integration nodes
	double WNode[4], FNode[4], GNode[4];

	double 	*C_fI_all, // Matrix of interpolated metric terms at all faces of Volume
			*Face_integral_met, // integral metric term multiplier from change of variable
			met1, met2, // metric terms arising in line integral 
			*Chi_f, *F_comm_1Var, cumulFaceIntVal;  

	double *fI_wi; // The weights at the integration nodes (one face)

	double *F_comm; // Common (numerical) flux

	struct S_FACE *FACE;
	struct S_VOLUME *VOLUME;
	struct S_ELEMENT *ELEMENT;

	// VIn Variables:
	struct S_VOLUME *VIn;
	int fin;
	// - Matrix of variables at face of inner volume
	double *What_in_fAll, *W_fIn, *F_in, *G_in; 

	// VOut Variables:
	struct S_VOLUME *VOut;
	int fout;
	// - Matrix of variables at face of outer volume
	double *What_out_fAll, *W_fOut, *F_out, *G_out; 


	// Set up arrays and structures:
	ELEMENT = DB.ELEMENT;

	// Loop over all faces
	for(FACE = DB.FACE_HEAD; FACE; FACE = FACE->next){
		if (FACE->Boundary != 1){
			// Pointers to the inner and outer volumes of the face
			VIn = FACE->VIn;
			VOut = FACE->VOut;

			// What face this is for each volume
			fin = FACE->fin;
			fout = FACE->fout;

			// Number of face integration nodes for this face
			numfI = (VIn->P+1);

			// VIn:
			// - Get solution matrix W_fIn at the face integration nodes
			// - Get the flux vector at each integration node using this vector

			What_in_fAll = malloc(VIn->NfnI*VIn->NVar* sizeof *What_in_fAll);  // free
			W_fIn = malloc(numfI*VIn->NVar* sizeof *W_fIn);  // free
			F_in = malloc(numfI*VIn->NVar* sizeof *F_in);  // free
			G_in = malloc(numfI*VIn->NVar* sizeof *G_in);  // free

			mm_CNN(VIn->NfnI, VIn->NvnG, VIn->NVar, ELEMENT->Chi_fI, VIn->What, What_in_fAll);

			// Fill the matrices of the solution and flux at the integration nodes
			// of the face for the inner volume
			for(ifI=0; ifI<numfI; ifI++){
				// Loop over the face integration nodes for ONE face

				for(iVar=0; iVar < VIn->NVar; iVar++){
					// Loop over the variables at the integration node

					W_fIn[iVar*numfI + ifI] = 
						What_in_fAll[iVar*(VIn->NfnI) // column
									+ fin*numfI + ifI]; // row

					WNode[iVar] = W_fIn[iVar*numfI + ifI];

				}

				// Compute the flux vector at each integration node
				euler_flux_2D(WNode, FNode, GNode);

				// Fill the F and G flux vectors for all equations 
				// for this integration node.
				for(iVar=0; iVar<VIn->NVar; iVar++){
					F_in[iVar*numfI + ifI] = FNode[iVar];
					G_in[iVar*numfI + ifI] = GNode[iVar];
				}
			}

			// VOut:
			// - Get solution matrix W_fOut at the face integration nodes
			//		- Flip the ordering of the nodes since we need the integration
			//			nodes to line up between the faces.
			// - Get the flux vector at each integration node using this vector

			What_out_fAll = malloc(VOut->NfnI*VOut->NVar* sizeof *What_out_fAll);  // free
			W_fOut = malloc(numfI*VOut->NVar* sizeof *W_fOut);  // free
			F_out = malloc(numfI*VOut->NVar* sizeof *F_out);  // free
			G_out = malloc(numfI*VOut->NVar* sizeof *G_out);  // free

			mm_CNN(VOut->NfnI, VOut->NvnG, VOut->NVar, ELEMENT->Chi_fI, VOut->What, What_out_fAll);

			// Fill the matrices of the solution and flux at the integration nodes
			// of the face for the inner volume
			for(ifI=0; ifI<numfI; ifI++){
				// Loop over the face integration nodes for ONE face

				for(iVar=0; iVar < VIn->NVar; iVar++){
					// Loop over the variables at the integration node

					W_fOut[iVar*numfI + ifI] = 
						What_out_fAll[iVar*(VOut->NfnI) // column
									+ fout*numfI + ifI]; // row

					WNode[iVar] = W_fOut[iVar*numfI + ifI];

				}

				// Compute the flux vector at each integration node
				euler_flux_2D(WNode, FNode, GNode);

				// Fill the F and G flux vectors for all equations 
				// for this integration node.
				for(iVar=0; iVar<VIn->NVar; iVar++){
					F_out[iVar*numfI + ifI] = FNode[iVar];
					G_out[iVar*numfI + ifI] = GNode[iVar];
				}
			}

			// Flip order of all VOut matrices
			flipMatrixOrder(W_fOut, numfI, VOut->NVar);
			flipMatrixOrder(F_out, numfI, VOut->NVar);
			flipMatrixOrder(G_out, numfI, VOut->NVar);

			// Compute the Numerical Flux
			F_comm = malloc(numfI*VIn->NVar* sizeof *F_comm);  // free
			flux_LF(W_fIn, W_fOut, F_in, F_out, G_in, G_out, F_comm, FACE->n_fI, VIn->P, VIn->NVar);

			if (DB.Testing == 1){
				// Print the common flux at the face
				
				printf("NEW FACE F_Comm: \n");
				for(i=0; i<numfI; i++){
					// Row loop
					for(j=0; j<(VIn->NVar); j++){
						// Column loop
						printf(" %.14f ", F_comm[j*(VIn->P+1) + i]);
					}
					printf("\n");

				}
			}

			// Store the F_Comm values for this face for both volumes
			for(i=0; i<numfI; i++){
				// Row loop
				for(j=0; j<(VIn->NVar); j++){
					// Column loop

					VIn->F_Comm[j*VIn->NfnI + fin*numfI + i] = F_comm[j*numfI + i];

					// order must be switched for outer volume
					VOut->F_Comm[j*VOut->NfnI + fout*numfI + i] = -1*F_comm[j*numfI + numfI-i-1];
				}
			}

			// Free temporary vectors
			free(F_comm);

			free(What_in_fAll);
			free(W_fIn);
			free(F_in);
			free(G_in);

			free(What_out_fAll);
			free(W_fOut);
			free(F_out);
			free(G_out);
		} else{
			// External Boundary Face:

			// Pointer to the inner volume
			VIn = FACE->VIn;
			fin = FACE->fin;
			
			// Number of face integration nodes for this face
			numfI = (VIn->P+1);

			// VIn:
			// - Get solution matrix W_fIn at the face integration nodes
			// - Get the flux vector at each integration node using this vector

			What_in_fAll = malloc(VIn->NfnI*VIn->NVar* sizeof *What_in_fAll);  // free
			W_fIn = malloc(numfI*VIn->NVar* sizeof *W_fIn);  // free
			F_in = malloc(numfI*VIn->NVar* sizeof *F_in);  // free
			G_in = malloc(numfI*VIn->NVar* sizeof *G_in);  // free

			mm_CNN(VIn->NfnI, VIn->NvnG, VIn->NVar, ELEMENT->Chi_fI, VIn->What, What_in_fAll);

			// Fill the matrices of the solution and flux at the integration nodes
			// of the face for the inner volume
			for(ifI=0; ifI<numfI; ifI++){
				// Loop over the face integration nodes for ONE face

				for(iVar=0; iVar < VIn->NVar; iVar++){
					// Loop over the variables at the integration node

					W_fIn[iVar*numfI + ifI] = 
						What_in_fAll[iVar*(VIn->NfnI) // column
									+ fin*numfI + ifI]; // row

					WNode[iVar] = W_fIn[iVar*numfI + ifI];

				}

				// Compute the flux vector at each integration node
				euler_flux_2D(WNode, FNode, GNode);

				// Fill the F and G flux vectors for all equations 
				// for this integration node.
				for(iVar=0; iVar<VIn->NVar; iVar++){
					F_in[iVar*numfI + ifI] = FNode[iVar];
					G_in[iVar*numfI + ifI] = GNode[iVar];
				}
			}

			// VOut:
			// - Get solution matrix W_fOut at the face integration nodes
			// 	using the proper boundary condition.
			//		- No flipping needed for matrix since boundary condition
			//			comes in correct order
			// - Get the flux vector at each integration node using this vector

			W_fOut = malloc(numfI*VOut->NVar* sizeof *W_fOut);  // free
			F_out = malloc(numfI*VOut->NVar* sizeof *F_out);  // free
			G_out = malloc(numfI*VOut->NVar* sizeof *G_out);  // free

			if(FACE->BCType == BC_SLIPWALL){
				boundary_SlipWall(W_fIn, FACE->n_fI, W_fOut, F_out, G_out, numfI);
			} else if (FACE->BCType == BC_BACKPRESSURE){
				boundary_BackPressure(W_fIn, FACE->n_fI, W_fOut, F_out, G_out, numfI);
			} else if (FACE->BCType == BC_TOTAL_TP){
				boundary_Total_TP(W_fIn, FACE->n_fI, W_fOut, F_out, G_out, numfI);
			} else{
				printf("Unsupported BC Type \n");
				exit(1);
			}

			// Compute the Numerical Flux
			F_comm = malloc(numfI*VIn->NVar* sizeof *F_comm);  // free
			flux_LF(W_fIn, W_fOut, F_in, F_out, G_in, G_out, F_comm, FACE->n_fI, VIn->P, VIn->NVar);

			if (DB.Testing == 1){
				// Print the common flux at the face
				
				printf("NEW FACE F_Comm: \n");
				for(i=0; i<numfI; i++){
					// Row loop
					for(j=0; j<(VIn->NVar); j++){
						// Column loop
						printf(" %.14f ", F_comm[j*(VIn->P+1) + i]);
					}
					printf("\n");
				}
			}

			// Store the F_Comm values for this face for both volumes
			for(i=0; i<numfI; i++){
				// Row loop
				for(j=0; j<(VIn->NVar); j++){
					// Column loop

					VIn->F_Comm[j*VIn->NfnI + fin*numfI + i] = F_comm[j*numfI + i];
				}
			}

			// Free temporary vectors
			free(F_comm);

			free(What_in_fAll);
			free(W_fIn);
			free(F_in);
			free(G_in);

			free(W_fOut);
			free(F_out);
			free(G_out);

		}
	}

	C_fI_all = malloc(4*DB.VOLUME_HEAD->NfnI* sizeof *C_fI_all);  // free
	Face_integral_met = malloc((DB.VOLUME_HEAD->P+1) * sizeof *Face_integral_met);  // free
	fI_wi = DB.ELEMENT->nodes_wi;

	// Compute the surface integral for each volume using the F_Comm
	// terms:
	int r,c;
	for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){


		// Interpolate the metric terms to the face integration nodes of all the
		// faces and get the metric jacobian term that comes from the transformation
		// of the line integral
		mm_CNN(VOLUME->NfnI, VOLUME->NvnS, 4, ELEMENT->I_vS_fI, VOLUME->C_vS, C_fI_all);

		if (DB.Testing == 1){
			printf("\n C_fI : \n");
			for(i=0; i<VIn->NfnI; i++){
				for(j=0; j<4; j++){
					printf("%.14f ", C_fI_all[j*VIn->NfnI + i]);
				}
				printf("\n");
			}
			printf("\n");
		}

		// Compute the RHS_FACE components for the volume
		for(iVar=0; iVar<VIn->NVar; iVar++){
			// Loop through the variables

			nBasis = 0;
			for(jBasis=0; jBasis<VIn->P+1; jBasis++){
				for(iBasis=0; iBasis<VIn->P+1; iBasis++){
					// Loop through the basis functions (nbasis = row of RHS_FACE)
					
					// Value of the i,j th basis function at the face integration
					// nodes
					Chi_f = &DB.ELEMENT->Chi_fI[nBasis*VOLUME->NfnI];
					
					cumulFaceIntVal = 0;
					for(i_face = 0; i_face < 4; i_face++){
						// Loop over all the faces
						numfI = VOLUME->P+1;
						
						F_comm_1Var = malloc(numfI*1* sizeof *F_comm_1Var);  // free
						
						for(i=0; i<numfI; i++){
							// Loop over the integration nodes

							// Get the face integration metrics for the nodes for the current face
							
							// Choose the correct metric terms depending on the face. 
							// Care needs to be taken in choosing the correct row in the matrix
							if (i_face == 0 || i_face == 2){
								// Constant eta faces. Take the xi metric terms
								met1 = C_fI_all[1*VOLUME->NfnI + i_face*numfI + i];
								met2 = C_fI_all[3*VOLUME->NfnI + i_face*numfI + i];
							} else {
								// Constant xi faces. Take the eta metric terms
								met1 = C_fI_all[0*VOLUME->NfnI + i_face*numfI + i];
								met2 = C_fI_all[2*VOLUME->NfnI + i_face*numfI + i];
							}

							Face_integral_met[i] = sqrt(met1*met1 + met2*met2);

							// Get the F_Comm values for the integration nodes for this
							// equation
							F_comm_1Var[i] = VOLUME->F_Comm[iVar*VOLUME->NfnI + // col
															i_face*numfI + i]; // row

						}

						// Add integral of this face to the integral value
						for(i=0; i < numfI; i++){
							cumulFaceIntVal = cumulFaceIntVal + 
							Chi_f[i_face*numfI + i]*F_comm_1Var[i]*Face_integral_met[i]*fI_wi[i];
						}

						free(F_comm_1Var);

					}

					// Assign value to RHS_FACE
					VOLUME->RHS_FACE[iVar*VOLUME->NvnG + nBasis] = cumulFaceIntVal;
					
					nBasis++;
				}
			}
		}

	}

	free(C_fI_all);
	free(Face_integral_met);

}






