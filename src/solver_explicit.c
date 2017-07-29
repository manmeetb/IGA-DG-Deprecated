
#include "solver_explicit.h"

#include <stdlib.h>
#include <stdio.h>

#include "update_VOLUMEs.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "explicit_VOLUME_info.h"
#include "explicit_FACE_info.h"
#include "finalize_RHS.h"
#include "output_solution.h"


void solver_explicit(void){

	/*
	Purpose:
		Solve the flow using the explicit solver.
	*/

	printf("Start Solver \n");

	// First, setup the mass matrix for each volume
	update_VOLUME_Ops();

	int time_step_type = 1; // 0 = euler, 1 = rk

	int t, i, rk;
	double *RHS_VOL, *RHS_FACE, maxRHS;
	int nBasis, iVar, iEq;
	struct S_VOLUME *VOLUME;

	double *RES, *RHS, *What;
	int NvnS, iMax;
	int Neq = 4;

	double dt = DB.dt;

	int output_solutionT = 2;
	// Temporary for testing purposes
	if (DB.numTimeSteps != 1){
		output_solutionT = (double)(DB.numTimeSteps)/3.;
	}


	// Loop through all time steps (for now do only one step)
	for(t=0; t<DB.numTimeSteps; t++){

		if (time_step_type == 0){

			// Compute the RHS:
			//	- Volume contribution
			printf("T = %d ", t);
			printf("V");  explicit_VOLUME_info();
			//	- Face contribution
			printf("F");  explicit_FACE_info();
			// 	- Combining both contributions
			printf("F "); maxRHS = finalize_RHS();

			// Perform the time stepping
			for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
				// Euler Explicit:
				for(i = 0; i < VOLUME->NvnG; i++){
					for(iEq = 0; iEq < 4; iEq ++){
						VOLUME->What[iEq*(VOLUME->NvnG) + i] = 
							VOLUME->What[iEq*(VOLUME->NvnG) + i] + dt*VOLUME->RHS[iEq*(VOLUME->NvnG) + i];
					}
				}

				if (DB.Testing == 1){
					if (VOLUME->index == 1){
					RHS_VOL = VOLUME->RHS_VOL;

					printf("VOLUME : %d \n ", VOLUME->index);
					printf("	RHS VOLUME: \n");
					for(nBasis=0; nBasis<VOLUME->NvnG; nBasis++){
						for(iVar = 0; iVar < VOLUME->NVar; iVar++){
							printf(" %.15f ", RHS_VOL[iVar*VOLUME->NvnG + nBasis]);
						}
						printf("\n");
					}

					RHS_FACE = VOLUME->RHS_FACE;
					printf("	RHS FACE: \n");
					for(nBasis=0; nBasis<VOLUME->NvnG; nBasis++){
						for(iVar = 0; iVar < VOLUME->NVar; iVar++){
							printf(" %.15f ", RHS_FACE[iVar*VOLUME->NvnG + nBasis]);
						}
						printf("\n");
					}
					}
				}

			}			

			printf(" maxRHS (no MInv): %.3e ", maxRHS);

			printf("\n");

			// Output 2 solutions in the middle of the run
			if(t%output_solutionT == 0){
				output_tecplot(t);
			}

			if(maxRHS < DB.exit_tol){
				printf("maxRHS Below %.4e. Exiting\n", DB.exit_tol);
				break;
			}

		}



		if (time_step_type == 1){

			printf("T = %d ", t);
			
			for (rk = 0; rk < 3; rk++) {
				// Build the RHS (== -Residual)
				printf("V");  explicit_VOLUME_info();
				printf("F");  explicit_FACE_info();
				printf("F "); maxRHS = finalize_RHS();

				// Update What
				for (VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next) {
					NvnS = VOLUME->NvnS;

					RES  = VOLUME->RES;
					RHS  = VOLUME->RHS;
					What = VOLUME->What;

					if (rk == 0) {
						for (iMax = Neq*NvnS; iMax--; ) {
							*RES++   = *What;
							*What++ += dt*(*RHS++);
						}
					} else if (rk == 1) {
						for (iMax = Neq*NvnS; iMax--; ) {
							*What = 0.25*(3.0*(*RES++) + *What + dt*(*RHS++));
							What++;
						}
					} else if (rk == 2) {
						for (iMax = Neq*NvnS; iMax--; ) {
							*What = (1.0/3.0)*(*RES++ + 2.0*(*What) + 2.0*dt*(*RHS++));
							What++;
						}
					}
				}
			}

			printf(" maxRHS (no MInv): %.3e ", maxRHS);

			printf("\n");

			if(maxRHS < DB.exit_tol){
				printf("maxRHS Below %.4e. Exiting\n", DB.exit_tol);
				break;
			}


		}

		if (DB.Testing){
			printf("\n ");
			printf("====================================================================== \n");
			printf("					Start Test \n");


			struct S_VOLUME *VOLUME_TEST;

			int numC, iTest, jTest;
			double *C11, *C12, *C21, *C22;

			for(VOLUME_TEST = DB.VOLUME_HEAD; VOLUME_TEST; VOLUME_TEST = VOLUME_TEST->next){
				
				printf("VOLUME : \n");
				// Print the volume grid point locations:
				for(iTest = 0; iTest < VOLUME_TEST->NvnG; iTest++){
					printf("	(x,y)_g : (%.15f, %.15f) \n", VOLUME_TEST->XYZ[iTest], VOLUME_TEST->XYZ[VOLUME_TEST->NvnG +iTest]);
				}

				// Print the metric terms
				numC = (VOLUME_TEST->P+1)*(VOLUME_TEST->P+1);

				C11 = &VOLUME_TEST->C_vS[0];
				C21 = &VOLUME_TEST->C_vS[2*numC];
				C12 = &VOLUME_TEST->C_vS[1*numC];
				C22 = &VOLUME_TEST->C_vS[3*numC];

				if(1){
					printf("- Metric Terms \n");
					for(iTest=0; iTest<numC; iTest++){
						printf("	C : (%.15f, %.15f, %.15f, %.15f) \n", 
							C11[iTest], C21[iTest], C12[iTest], C22[iTest]);
						printf("	det J : %.15f \n", VOLUME_TEST->detJV_vS[iTest]);
					}
				}

				if(1){
					printf("	RHS_VOL: \n");
					for(iTest = 0; iTest<VOLUME_TEST->NvnG; iTest++){
						// Loop through the rows
						for(jTest=0; jTest<VOLUME_TEST->NVar; jTest++){
							// Loop through the columns
							printf("%.15f ", VOLUME_TEST->RHS_VOL[jTest*VOLUME_TEST->NvnG + iTest]);
						}
						printf("\n");

					}

					printf("	RHS_FACE: \n");
					for(iTest = 0; iTest<VOLUME_TEST->NvnG; iTest++){
						// Loop through the rows
						for(jTest=0; jTest<VOLUME_TEST->NVar; jTest++){
							// Loop through the columns
							printf("%.15f ", VOLUME_TEST->RHS_FACE[jTest*VOLUME_TEST->NvnG + iTest]);
						}
						printf("\n");

					}

					printf("	RHS_VOL - RHS_FACE: \n");
					for(iTest = 0; iTest<VOLUME_TEST->NvnG; iTest++){
						// Loop through the rows
						for(jTest=0; jTest<VOLUME_TEST->NVar; jTest++){
							// Loop through the columns
							printf("%.15f ", VOLUME_TEST->RHS_VOL[jTest*VOLUME_TEST->NvnG + iTest] - VOLUME_TEST->RHS_FACE[jTest*VOLUME_TEST->NvnG + iTest]);
						}
						printf("\n");

					}

					printf("	RHS: \n");
					for(iTest = 0; iTest<VOLUME_TEST->NvnG; iTest++){
						// Loop through the rows
						for(jTest=0; jTest<VOLUME_TEST->NVar; jTest++){
							// Loop through the columns
							printf("%.15f ", VOLUME_TEST->RHS[jTest*VOLUME_TEST->NvnG + iTest]);
						}
						printf("\n");

					}

				}

			}

			printf("					End Test \n");
			printf("====================================================================== \n");
			printf("\n");
		}

	}


}


