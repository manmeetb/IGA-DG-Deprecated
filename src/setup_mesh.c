
#include "setup_mesh.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "S_DB.h"
#include "Parameters.h"
#include "setup_volumes.h"
#include "setup_faces.h"
#include "S_BC.h"
#include "memory_constructors.h"
#include "memory_destructors.h"

/*
 *	Purpose:
 *		Set up mesh related parameters. This module will read the 
 *		mesh file with the geometry node and connectivity data.
 *
 *	Comments:
 *
 *	Notation:
 *		NV = (N)umber of (V)olumes
 *
 *		XYZ_Ve = XYZ coordinates of (Ve)rtices (for whole mesh)
 *		VeCon = (Ve)rtices (Con)nectivity arrays for each volume
 *		NVe = (N)umber of (Ve)rtices
 *	
 *		XYZ_G = XYZ coordinates of (Geometry) nodes (for whole mesh)
 *		GeoCon = (Geo)metry Nodes (Con)nectivity arrays for each volume
 *		NGConPerV = (N)umber of (G)eometry (Con)nectivity entries (Per)
 *			(V)olume
 *		NGP = (N)umber of (G)eometry node (P)oints
 *
 *		NPF = (N)umber of (P)eriodic (F)aces
 *		PCon = (P)eriodic (Con)nectivity arrays
 *	
 *
 *	References:
 */

static void setup_readMeshFile(){

	/*
	Mesh File:
	- Node coordinates are all listed one after another (geometry node points)
	- Connectivity array gives the value of each geometry node point in the element. The 
		points are listed with all j = 0, then j = 1, ... where i=0, j=0 is at the 
		bottom left corner of the reference element which the physical element maps to. 
	- Periodic boundary conditions list what faces (between which vertices) each element is 
		connected. Using this, we can find what xi or eta line is touching between elements
		to know which is the right and which is the left element.
	*/

	printf("Process Mesh File \n");

	char StringRead[200];
	int i, j, NGConPerV;
	FILE *file_id;

	if ((file_id = fopen(DB.MeshFileName,"r")) == NULL){
		printf("Mesh file: %s not present.\n",DB.MeshFileName);
		exit(1);
	}

	// ----------------------------------------------------------
	//					Mesh Properties
	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	fscanf(file_id, "%[^\n]\n", StringRead);
	sscanf(StringRead, "%d", &DB.NGP); // Num Grid Points

	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	fscanf(file_id, "%[^\n]\n", StringRead);
	sscanf(StringRead, "%d", &DB.NVe); // Num Vertices

	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	fscanf(file_id, "%[^\n]\n", StringRead);
	sscanf(StringRead, "%d", &DB.NV); // Num Quads

	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	fscanf(file_id, "%[^\n]\n", StringRead);
	sscanf(StringRead, "%d", &DB.P); // Order P
	// ----------------------------------------------------------

	// ----------------------------------------------------------
	//						Vertices
	// Read the vertices:
	double *XYZ_Ve, *X_Ve, *Y_Ve;

	XYZ_Ve = malloc(DB.NVe*DB.d* sizeof *XYZ_Ve);
	X_Ve = &XYZ_Ve[0];
	Y_Ve = &XYZ_Ve[DB.NVe];

	fscanf(file_id, "%[^\n]\n", StringRead); // Header
	for(i=0; i<DB.NVe; i++){
		fscanf(file_id, "%[^\n]\n", StringRead);
		sscanf(StringRead, "%lf %lf", &X_Ve[i], &Y_Ve[i]);
	}

	// Read Vertices Connectivity:
	// - Matrix stored in column major form (each column for each element)
	int *VeCon;
	VeCon = malloc(DB.NV*4* sizeof *VeCon);
	fscanf(file_id, "%[^\n]\n", StringRead); // Header
	for(i=0; i<DB.NV; i++){
		fscanf(file_id, "%[^\n]\n", StringRead);

		// read in the information into this quad
		for(j=0; j<4; j++){
			sscanf(StringRead, "%d %[^\t\n]", &VeCon[i*4+j], StringRead);
		}
	}
	// ----------------------------------------------------------

	// ----------------------------------------------------------
	//					Geometry Nodes
	// - For storing all the xyz geometry node points in column major form.
	double * XYZ_G, *X_G, *Y_G; 
	XYZ_G = malloc(DB.NGP*DB.d* sizeof *XYZ_G);
	X_G = &XYZ_G[0];
	Y_G = &XYZ_G[DB.NGP];
	
	// Read all geometry node points
	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	for(i=0; i<DB.NGP; i++){
		fscanf(file_id, "%[^\n]\n", StringRead);
		sscanf(StringRead, "%lf %lf", &X_G[i], &Y_G[i]);
	}


	// Setup temporary array to hold quad connectivity information
	int *GeoCon;
	NGConPerV = DB.P+1;
	for(i=0; i<DB.d-1; i++){
		NGConPerV = NGConPerV*(DB.P+1);
	}

	// Matrix of size (NGConPerV x NV).
	// Stored in column major form.
	GeoCon = 
		malloc(DB.NV*NGConPerV* sizeof *GeoCon);

	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	for(i=0; i<DB.NV; i++){
		fscanf(file_id, "%[^\n]\n", StringRead);

		// read in the information into this quad
		for(j=0; j<NGConPerV; j++){
			sscanf(StringRead, "%d %[^\t\n]", &GeoCon[i*NGConPerV+j], StringRead);
		}
	}
	// ----------------------------------------------------------

	// ----------------------------------------------------------
	// 					Boundary Conditions

	int numBCTypes, index_BCType, numBCRows, numBCCols;
	struct S_BC *BC, *BC_Next;

	// Find the number of BCs present in the mesh
	fscanf(file_id, "%[^\n]\n", StringRead); // Header Line
	fscanf(file_id, "%[^\n]\n", StringRead); // Num BC Types
	sscanf(StringRead, "%d", &numBCTypes);

	BC = New_BC();
	DB.BC_HEAD = BC;

	for(i=0; i<numBCTypes-1; i++){
		BC_Next = New_BC();
		BC->next = BC_Next;
		BC = BC_Next;
	}

	// Load all the BC information now
	BC = DB.BC_HEAD;
	for(index_BCType=0; index_BCType < numBCTypes; index_BCType++){

		fscanf(file_id, "%[^\n]\n", StringRead); // BC Type Line
		
		// Find the type of BC this is and set that parameter in the BC struct
		if(strstr(StringRead, "Periodic")){
			BC->BCType = BC_PERIODIC;
		} else if(strstr(StringRead, "SlipWall")){
			BC->BCType = BC_SLIPWALL;
		} else if(strstr(StringRead, "TotalTemperaturePressure")){
			BC->BCType = BC_TOTAL_TP;
		} else if(strstr(StringRead, "BackPressure")){
			BC->BCType = BC_BACKPRESSURE;
		} else{
			printf("Unsupported BC Type \n");
			exit(1);
		}

		fscanf(file_id, "%[^\n]\n", StringRead); // Number of BC Rows
		sscanf(StringRead, "%d", &numBCRows);

		// Number of columns depends on the type of BC. Periodic will have 4 
		// columns for the BC matrix (n bcs x 4 (elem_1_index, elem_1_face, ...))
		// All other BCs will have only 2
		if(BC->BCType == BC_PERIODIC){
			numBCCols = 4;
		} else{
			numBCCols = 2;
		}

		BC->nBC_Con_row = numBCRows;
		BC->nBC_Con_col = numBCCols;

		BC->BC_Con = malloc(BC->nBC_Con_row*BC->nBC_Con_col* sizeof *BC->BC_Con);

		// Read all the BC lines into the matrix in row major form
		for(i=0; i<numBCRows; i++){
			fscanf(file_id, "%[^\n]\n", StringRead);
			for(j=0; j<numBCCols; j++){
				sscanf(StringRead, "%d %[^\t\n]", 
								&(BC->BC_Con[i*BC->nBC_Con_col+j]), 
								StringRead);
			}
		}

		// Load the next data into the next BC struct
		BC = BC->next;
	}


	// ----------------------------------------------------------

	DB.XYZ_G = XYZ_G;
	DB.XYZ_Ve = XYZ_Ve;

	DB.NGConPerV = NGConPerV;

	DB.GeoCon = GeoCon;
	DB.VeCon = VeCon;

}

void setup_mesh(){

	// Read the mesh file
	setup_readMeshFile();

	// Setup the volumes
	setup_volumes();

	// Setup the faces
	setup_faces();


	int i,j;
	struct S_BC *BC, *BC_NEXT;
	
	if (DB.Testing == 1 || DB.Testing == 2){
		// Print all the information
		printf("Vertices \n");
		for(i=0; i<DB.NVe; i++){
			printf("(x,y): %f %f \n", DB.XYZ_Ve[i], DB.XYZ_Ve[DB.NVe + i]);
		}

		printf("Vertices Connectivity \n");
		for(i=0; i<DB.NV; i++){
			printf("Connect: ");
			for(j=0; j<4; j++){
				printf("%d ", DB.VeCon[i*4+j]);
			}
			printf("\n");
		}

		printf("Geometry Nodes \n");
		for(i=0; i<DB.NGP; i++){
			printf("(x,y): %f %f \n", DB.XYZ_G[i], DB.XYZ_G[DB.NGP + i]);
		}	

		printf("Geometry Nodes Connectivity \n");
		for(i=0; i<DB.NV; i++){
			printf("Connect: ");
			for(j=0; j<DB.NGConPerV; j++){
				printf("%d ", DB.GeoCon[i*DB.NGConPerV+j]);
			}
			printf("\n");
		}

		// BCs:
		for(BC=DB.BC_HEAD; BC; BC = BC->next){
			printf("BC Type : %d \n", BC->BCType);

			for(i=0; i<BC->nBC_Con_row; i++){
				for(j=0; j<BC->nBC_Con_col; j++){
					printf("%d ", BC->BC_Con[i*BC->nBC_Con_col+j]);
				}
				printf("\n");
			}
		}
	}


	// Free temporary memory
	free(DB.XYZ_G);
	free(DB.XYZ_Ve);
	free(DB.GeoCon);
	free(DB.VeCon);

	// - BC Linked List
	BC = DB.BC_HEAD;
	while(1){
		BC_NEXT = BC->next;

		// Free allocated arrays
		free(BC->BC_Con);

		// Last BC case
		if(BC_NEXT == NULL){
			memory_destructor_BC(BC);
			break;
		}
		
		memory_destructor_BC(BC);
		BC = BC_NEXT;
	}

}





