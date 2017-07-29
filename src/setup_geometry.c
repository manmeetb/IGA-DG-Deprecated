
#include "setup_geometry.h"

#include <stdlib.h>
#include <stdio.h>

#include "setup_geom_factors.h"
#include "matrix_functions.h"
#include "setup_normals.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_ELEMENT.h"
#include "S_FACE.h"



void setup_geometry(void){
	printf("Setup Geometry \n");

	struct S_VOLUME *VOLUME;
	struct S_FACE *FACE;

	double *XYZ_Sx, *XYZ_Sy;
	double *C_vS_11, *C_vS_12, *C_vS_21, *C_vS_22;
	int i;

	// ----------------------------------------
	//		Volume/Plotting Solution Nodes
	// ----------------------------------------
	// Set up XYZ_S and XYZ_P

	for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){

		mm_CNN(VOLUME->NvnS, VOLUME->NvnG, VOLUME->d, DB.ELEMENT->Chi_vS, VOLUME->XYZ, VOLUME->XYZ_S);
		mm_CNN(VOLUME->NvnP, VOLUME->NvnG, VOLUME->d, DB.ELEMENT->Chi_vP, VOLUME->XYZ, VOLUME->XYZ_P);

		// Print solution node location of each volume
		if(DB.Testing == 1 || DB.Testing == 4){
			XYZ_Sx = &VOLUME->XYZ_S[0];
			XYZ_Sy = &VOLUME->XYZ_S[VOLUME->NvnS];

			printf("Volume: %d \n" , VOLUME->index);
			for(i=0; i<VOLUME->NvnS; i++){
				printf("Sol (x,y): (%.14f, %.14f) \n", XYZ_Sx[i], XYZ_Sy[i]);
			}
		}

		// Setup Geometric Factors
		setup_geom_factors(VOLUME);

		if(DB.Testing == 1){

			C_vS_11 = &VOLUME->C_vS[0];
			C_vS_12 = &VOLUME->C_vS[1*VOLUME->NvnS];
			C_vS_21 = &VOLUME->C_vS[2*VOLUME->NvnS];
			C_vS_22 = &VOLUME->C_vS[3*VOLUME->NvnS];

			printf("Volume: %d \n" , VOLUME->index);
			for(i=0; i<VOLUME->NvnS; i++){
				printf("C : (%f, %f, %f, %f) \n", 	C_vS_11[i], C_vS_12[i],
													C_vS_21[i], C_vS_22[i]);

				printf("det J : %f \n", VOLUME->detJV_vS[i]);
			}

		}

	}
	
		// Setup Normals
	for(FACE=DB.FACE_HEAD; FACE; FACE = FACE->next){
		setup_normals(FACE);
	}

	if(DB.Testing == 1){
		printf("Volume Face Normals \n");

		int i, j;
		double *nx, *ny;

		struct S_FACE  *FACE;
		for(VOLUME=DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
			for(i=0; i<4; i++){
				printf("Face Processed : %d \n", VOLUME->FIndFound[i]);

				FACE = VOLUME->FACE[i];

				nx = &FACE->n_fI[0];
				ny = &FACE->n_fI[FACE->P+1];

				printf("face in = %d \n", FACE->fin);

				for(j=0; j<FACE->P+1; j++){
					printf("	nx = %f ny = %f \n", nx[j], ny[j]);
				}

			}
		}

	}






}


