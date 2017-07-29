
#include "setup_volumes.h"

#include <stdlib.h>
#include <stdio.h>

#include "S_DB.h"
#include "S_VOLUME.h"
#include "memory_constructors.h"


/*
 *	Purpose:
 *		Holds the methods for setting up the structs for the Volumes and 
 *		Faces. 
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void setup_volumes(void){
	printf("Setup Volumes \n");

	// Variables:
	int i, VIndex;
	double *XVals, *YVals;
	struct S_VOLUME *VOLUME, *PrevVOLUME;

	// Create first the NV volume structs linked list
	VOLUME = New_VOLUME();
	DB.VOLUME_HEAD = VOLUME;

	PrevVOLUME = DB.VOLUME_HEAD;

	// Create and link the next NV-1 volumes
	for(i=0; i<DB.NV-1; i++){
		VOLUME = New_VOLUME();
		VOLUME->parent = PrevVOLUME;
		PrevVOLUME->next = VOLUME;
		PrevVOLUME = VOLUME;
	}

	// Fill the mesh information for each volume and initialize any 
	// arrays
	VIndex=0;
	for (VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
		VOLUME->index = VIndex;

		// Vertices Connectivity Information (ccw)
		// Fill the indeces of the vertices for each volume
		for(i=0; i<4; i++){
			VOLUME->VeInd[i] = DB.VeCon[VIndex*4+i];
		}

		// Properties of each Volume
		VOLUME->P = DB.P;
		VOLUME->d = DB.d;
		VOLUME->NVar = 4;
		VOLUME->update = 1; // Volume just created so it's ops must be updated

		VOLUME->NvnG = (DB.P+1)*(DB.P+1);
		VOLUME->NvnS = VOLUME->NvnG;
		VOLUME->NvnP = VOLUME->NvnG;
		VOLUME->NfnI = 4*(VOLUME->P+1);

		// Geometry Structures:
		VOLUME->XYZ = malloc(DB.d*VOLUME->NvnG * sizeof *VOLUME->XYZ);
		VOLUME->XYZ_S = malloc(DB.d*VOLUME->NvnS * sizeof *VOLUME->XYZ_S);
		VOLUME->XYZ_P = malloc(DB.d*VOLUME->NvnP * sizeof *VOLUME->XYZ_P);

		// - Setup the geometry node points
		XVals = &VOLUME->XYZ[0];
		YVals = &VOLUME->XYZ[VOLUME->NvnG];
		for(i=0; i<DB.NGConPerV; i++){
			XVals[i] = DB.XYZ_G[DB.GeoCon[VIndex*DB.NGConPerV + i]];
			YVals[i] = DB.XYZ_G[DB.GeoCon[VIndex*DB.NGConPerV + i] + DB.NGP];
		}

		// Set initial processed face to be 0 for all faces
		for(i=0; i<4; i++){
			VOLUME->FIndFound[i] = 0;
		}

		VOLUME->FACE = malloc(4*sizeof *VOLUME->FACE);

		// Metric Structures:
		VOLUME->C_vS = malloc(VOLUME->NvnS*4* sizeof *VOLUME->C_vS);
		VOLUME->detJV_vS = malloc(VOLUME->NvnS * sizeof *VOLUME->detJV_vS);

		// Solution Structures:
		VOLUME->What = malloc(VOLUME->NVar*VOLUME->NvnS* sizeof *VOLUME->What);
		VOLUME->RHS = malloc(VOLUME->NvnG*VOLUME->NVar* sizeof *VOLUME->RHS);
		VOLUME->RES = malloc(VOLUME->NvnG*VOLUME->NVar* sizeof *VOLUME->RES);
		VOLUME->RHS_VOL = malloc(VOLUME->NvnG*VOLUME->NVar* sizeof *VOLUME->RHS_VOL);
		VOLUME->RHS_FACE = malloc(VOLUME->NvnG*VOLUME->NVar* sizeof *VOLUME->RHS_FACE);
		VOLUME->MInv = malloc(VOLUME->NvnG*VOLUME->NvnG* sizeof *VOLUME->MInv);
		VOLUME->F_Comm = malloc(VOLUME->NfnI*VOLUME->NVar* sizeof *VOLUME->F_Comm);

		VIndex++;
	}

	if (DB.Testing == 3){
		// Print the Geometry node points for each volume
		for (VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
			printf("\n New Volume \n");
			for(i=0; i<4; i++){
				printf("%d ", VOLUME->VeInd[i]);
			}
			printf("\n");

			for(i=0; i<(DB.P+1)*(DB.P+1); i++){
				printf("(x,y): (%f, %f) \n", VOLUME->XYZ[i], VOLUME->XYZ[(DB.P+1)*(DB.P+1) + i]);
			}
		}
	}

}








