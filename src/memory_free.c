

#include "memory_free.h"

#include <stdlib.h>
#include <stdio.h>

#include "memory_destructors.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_ELEMENT.h"

void memory_free(void){
	int i;

	// DB Parameters:
	free(DB.NodeType);
	free(DB.MeshFileName);
	free(DB.BasisType);
	free(DB.TestType);

	// Structures Freed in setup_mesh:
	//free(DB.XYZ_G);
	//free(DB.XYZ_Ve);
	//free(DB.GeoCon);
	//free(DB.VeCon);
	//free(DB.PCon);

	// Structs:
	// - VOLUME Linked List
	struct S_VOLUME *VOLUME, *VOLUME_NEXT;
	VOLUME = DB.VOLUME_HEAD; 
	while(1){
		VOLUME_NEXT = VOLUME->next;

		// Free allocated arrays

		// Geometry: 
		free(VOLUME->XYZ);
		free(VOLUME->XYZ_S);
		free(VOLUME->XYZ_P);

		// Metric Terms
		free(VOLUME->C_vS);
		free(VOLUME->detJV_vS);

		// Solution: 
		free(VOLUME->What);
		free(VOLUME->MInv);
		free(VOLUME->RHS_FACE);
		free(VOLUME->RHS_VOL);
		free(VOLUME->RHS);
		free(VOLUME->RES);
		free(VOLUME->F_Comm);

		// Other Structs
		free(VOLUME->FACE);

		// Last VOLUME case
		if(VOLUME_NEXT == NULL){
			memory_destructor_V(VOLUME);
			break;
		}
		
		memory_destructor_V(VOLUME);
		VOLUME = VOLUME_NEXT;
	}

	// - FACE Linked List
	struct S_FACE *FACE, *FACE_NEXT;
	FACE = DB.FACE_HEAD;
	while(1){
		FACE_NEXT = FACE->next;

		// Free allocated arrays
		free(FACE->n_fI);
		free(FACE->C_fI);

		// Last FACE case
		if(FACE_NEXT == NULL){
			memory_destructor_F(FACE);
			break;
		}
		
		memory_destructor_F(FACE);
		FACE = FACE_NEXT;
	}

	// - ELEMENT
	// Allocated Arrays
	free(DB.ELEMENT->nodes_xi);
	free(DB.ELEMENT->nodes_wi);

	// Operators

	// - Chi Operators
	free(DB.ELEMENT->Chi_vS);
	free(DB.ELEMENT->ChiInv_vS);
	free(DB.ELEMENT->Chi_vG);
	free(DB.ELEMENT->Chi_vP);
	free(DB.ELEMENT->Chi_fI);

	// - Derivative Operators
	free(DB.ELEMENT->GradChi_vS_xi);
	free(DB.ELEMENT->GradChi_vS_eta);
	free(DB.ELEMENT->GradChi_fI_xi);
	free(DB.ELEMENT->GradChi_fI_eta);

	// - Interpolation Operators
	free(DB.ELEMENT->I_vG_vS);
	free(DB.ELEMENT->I_vS_vG);
	free(DB.ELEMENT->I_vG_vP);
	free(DB.ELEMENT->I_vS_fI);

	free(DB.ELEMENT->XiEtaZeta_S);
	free(DB.ELEMENT->XiEtaZeta_G);
	free(DB.ELEMENT->XiEtaZeta_P);
	free(DB.ELEMENT->XiEtaZeta_F);

	free(DB.ELEMENT->xiVector);
	free(DB.ELEMENT->etaVector);

	memory_destructor_E(DB.ELEMENT);
	
}



