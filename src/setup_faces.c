
#include "setup_faces.h"

#include <stdlib.h>
#include <stdio.h>

#include "memory_constructors.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_BC.h"
#include "Parameters.h"

/*
 *	Purpose:
 *		Set up the faces between the volumes on the mesh. This 
 *		will therefore create the S_FACE structures which 
 *		will hold references to the volumes that compose each given 
 *		face.
 *
 *	Comments:
 *
 *	Notation:
 *	
 *
 *	References:
 */

void setup_faces(void){

	printf("Setup Faces \n");

	struct S_VOLUME *VOLUME, *VIn, *VOut;
	struct S_FACE *FACE, *PrevFACE; 
	struct S_BC *BC, *BCLoop;

	int ni0, ni1, ni2, ni3, fi_index, // Nodes for each face information
		no0, no1, no2, no3, fo_index,
		fi[4][2], fo[4][2],	
		i, j, check, 
		perLIndex, perRIndex,  // (Per)iodic (L)eft/(R)ight Vol Index
		perLF, perRF, //(Per)iodic (L)eft/(R)ight (F)ace
		nPer1, nPer2; // node vertices for periodic faces

	int ExternalBCIndex, ExternalBCF;

	/*
	Volumme - Volume inner face connections:

	Loop through all the inner volumes. For each volume,
	loop through all other volumes. If another volume is found with 
	two nodes in common, then we have found a face. There will be a 
	convention for the face ordering on each element based on the node
	ordering:

		  f2
		n3	 n2
		*----*
	 f3	|    |  f1
		*----*
		n0	 n1
		  f0
	
	That is, given the node ordering, if the two nodes found to form a face
	for instance are the first two nodes in the list of node vertices for the
	element, this will correspond to face 0 for that element. The face index
	is important because the mapping of the node vertices to the mapped element
	on the computational domain will always be the same (n0 => xi = -1, eta = -1
	and n2 => xi = 1, eta = 1). 

	Hold all S_FACE structs in a linked list.

	NOTE: All volumes will have a counter clockwise node ordering in the mesh.
		So, integration node ordering of inner and outer volumes will need to flipped (one will)
		during computations.
	*/

	// Search through the BCs structs to find if there are any periodic BCs
	BC = NULL;
	for(BCLoop = DB.BC_HEAD; BCLoop; BCLoop = BCLoop->next){
		if(BCLoop->BCType == BC_PERIODIC){
			BC = BCLoop;
			break;
		}
	}

	DB.FACE_HEAD = NULL;

	for(VIn = DB.VOLUME_HEAD; VIn; VIn = VIn->next){
		
		// Indeces of the vertex nodes for VIn
		ni0 = VIn->VeInd[0];
		ni1 = VIn->VeInd[1];
		ni2 = VIn->VeInd[2];
		ni3 = VIn->VeInd[3];

		// Faces for VIn (using index convention provided above)
		fi[0][0] = ni0; fi[0][1] = ni1;
      	fi[1][0] = ni1; fi[1][1] = ni2;
      	fi[2][0] = ni2; fi[2][1] = ni3;
      	fi[3][0] = ni3; fi[3][1] = ni0;

		for(VOut = DB.VOLUME_HEAD; VOut; VOut = VOut->next){
			
			// Indeces of the vertex node for VOut
			no0 = VOut->VeInd[0];
			no1 = VOut->VeInd[1];
			no2 = VOut->VeInd[2];
			no3 = VOut->VeInd[3];	

			// Faces for VOut
			fo[0][0] = no0; fo[0][1] = no1;
	      	fo[1][0] = no1; fo[1][1] = no2;
	      	fo[2][0] = no2; fo[2][1] = no3;
	      	fo[3][0] = no3; fo[3][1] = no0;

			if (VIn != VOut){
				// Have two non common volumes

				if(BC != NULL){
					// 			Periodic Faces
					// Consider the periodic faces here. If a periodic face
					// does exist between these elements, then set the 
					// VOut volume to have the same nodes as the VIn volume for
					// the corresponding face in order for a match to be found.
					// Recall: Periodic data given as
					//  Vol1_index Vol2_index Face1_index Face2_index

					for(i=0; i<BC->nBC_Con_row; i++){
						perLIndex = BC->BC_Con[BC->nBC_Con_col*i];
						perRIndex = BC->BC_Con[BC->nBC_Con_col*i+1];
						perLF = BC->BC_Con[BC->nBC_Con_col*i+2];
						perRF = BC->BC_Con[BC->nBC_Con_col*i+3];

						if(VIn->index == perLIndex && VOut->index == perRIndex){
							// Have two periodic elements
							
							// VIn nodes for periodic face
							nPer1 = fi[perLF][0];
							nPer2 = fi[perLF][1];

							// Set VOut face nodes to be same as VIn due to periodicity
							fo[perRF][0] = nPer1;
							fo[perRF][1] = nPer2;
						}

						if(VIn->index == perRIndex && VOut->index == perLIndex){
							// Have two periodic elements
							
							// VIn nodes for periodic face
							nPer1 = fi[perRF][0];
							nPer2 = fi[perRF][1];

							// Set VOut face nodes to be same as VIn due to periodicity
							fo[perLF][0] = nPer1;
							fo[perLF][1] = nPer2;
						}
					}
				}

				// Loop through the faces to find the common inner faces
				for(fi_index=0; fi_index<4; fi_index++){
					for(fo_index=0; fo_index<4; fo_index++){
						
						// Flag that will become 1 if a face has been found
						check = 0;

						if(((fi[fi_index][0] == fo[fo_index][0]) && (fi[fi_index][1] == fo[fo_index][1])) ||
						   ((fi[fi_index][1] == fo[fo_index][0]) && (fi[fi_index][0] == fo[fo_index][1]))){
							// Found a matching face
							check = 1;
						}

						if(check==1){
							// Have a matching face

							if(VIn->FIndFound[fi_index]==0 && VOut->FIndFound[fo_index]==0){
								// If this face on the volume has not been processed yet

								VIn->FIndFound[fi_index] = 1;
								VOut->FIndFound[fo_index] = 1;

								FACE = New_FACE();

								if(DB.FACE_HEAD == NULL){

									// This is the first face being processed
									FACE->parent = NULL;
									DB.FACE_HEAD = FACE;
									PrevFACE = DB.FACE_HEAD;

								} else{

									FACE->parent = PrevFACE;
									PrevFACE->next = FACE;
									PrevFACE = FACE;
								}

								FACE->fin = fi_index;
								FACE->fout = fo_index;

								FACE->VIn = VIn;
								FACE->VOut = VOut;

								//Store face references
								VIn->FACE[fi_index] = FACE;
								VOut->FACE[fo_index] = FACE;

								// Face Properties:
								FACE->P = VIn->P;
								FACE->Boundary = 0;
								FACE->BCType = BC_INTERNAL;

								// Create any array structures needed by the face
								FACE->n_fI = malloc((FACE->P+1)*DB.d* sizeof *FACE->n_fI);
								FACE->C_fI = malloc((FACE->P+1)*4* sizeof *FACE->C_fI); 
							}

						} // End If Check

					} // End for fo_index
				} // End for fi_index

			} // End If VIn != VOut
		} // End for VOut
	} // End for VIn


	if(PrevFACE == NULL || FACE == NULL){
		printf("Need at least one internal face \n");
		exit(1);
	}

	// Set external boundary faces
	for(BC = DB.BC_HEAD; BC; BC = BC->next){
		if(BC->BCType != BC_PERIODIC){
			// Build faces for all other non periodic faces
			
			for(i=0; i<BC->nBC_Con_row; i++){
				ExternalBCIndex = BC->BC_Con[BC->nBC_Con_col*i];
				ExternalBCF = BC->BC_Con[BC->nBC_Con_col*i + 1];

				// Find the volume and create this face
				for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
					if(VOLUME->index == ExternalBCIndex){

						// Create the Face
						FACE = New_FACE();
						FACE->parent = PrevFACE;
						PrevFACE->next = FACE;
						PrevFACE = FACE;

						FACE->fin = ExternalBCF;
						FACE->fout = -1;

						FACE->VIn = VOLUME;
						FACE->VOut = NULL;

						//Store face references
						VOLUME->FACE[ExternalBCF] = FACE;
						VOLUME->FIndFound[ExternalBCF] = 1;

						// Face Properties:
						FACE->P = VOLUME->P;
						FACE->Boundary = 1;
						FACE->BCType = BC->BCType;

						// Create any array structures needed by the face
						FACE->n_fI = malloc((FACE->P+1)*DB.d* sizeof *FACE->n_fI);
						FACE->C_fI = malloc((FACE->P+1)*4* sizeof *FACE->C_fI); 

						break;

					}

				}

			}

		}
	}


	if(DB.Testing == 1 || DB.Testing == 2){
		i=0;
		for(FACE = DB.FACE_HEAD; FACE; FACE = FACE->next){
			printf("FACE: %d \n", i);

			printf("	fi: %d fo: %d \n", FACE->fin, FACE->fout);

			printf("	Boundary : %d \n", FACE->Boundary);
			printf("	BCType : %d \n", FACE->BCType);

			if(!FACE->Boundary){
				printf("	VIn: ");
				for(j=0; j<4; j++){
					printf(" %d ", FACE->VIn->VeInd[j]);
				}
				printf("\n");
				printf("	VOut: ");
				for(j=0; j<4; j++){
					printf(" %d ", FACE->VOut->VeInd[j]);
				}
				printf("\n");
			}

			i++;
		}
	}


}

