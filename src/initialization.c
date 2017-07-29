
#include "initialization.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "S_DB.h"
#include "Parameters.h"

/*
* Purpose:
*	Module to be used to read in the control file when we have one.
* 	For now, we will simply set everything using constants.

* Notation:
 *
 *		Code Parameters
 *			Dimension  : d = 2  (Will add d=3 after)
 *			MeshType   : Optional
 *							- Will have a flag to curve the mesh potentially if needed. Most 
 *								likely will do this in the mesh generation file itself.
 *
 *			Form       : Form of the equations (i.e. how many times they are itnegrated by parts)
 *			             Options: Weak (Only implement the weak form for now)
 *
 *			NodeType   : Type of VOLUME nodes to use for different element types
 *			             Options: (T)ensor(P)roduct : (G)auss(L)egendre -> Only option implemented currently
 *			                                          (G)auss(L)obatto(L)egendre
 *			BasisType  : Type of basis functions
 *			             Options: Modal
 *			P          : Polynomial order to be used (not used if p-adaptation is enabled)
 *
 *			Testing    : Run tests for standard checks.
 *			             Options: 0 (No testing)
 *			                      1 (Testing)
*/

void initialization(int nargc, char **argv){

	/*
	Read and process the control file. The control file
	should be given as a command line argument to the solver
	
	Parameters:
		int nargc: Number of command line arguments
		char **argv: Array of strings for each command line argument
	
	Return:

	*/

	// Validation:
	// - There should be two arguments

	if (nargc != 2){
		printf("NEED 2 COMMAND LINE ARGUMENTS \n");
		exit(1);
	}

	char *CtrlFile; 

	// Setup Control File Parameters
	CtrlFile = malloc(200 * sizeof *CtrlFile);  // free
	strcpy(CtrlFile, argv[1]);
	printf("Control File : %s \n", CtrlFile);

	// -----------------------------------------------------------------
	//						Read Control File
	// -----------------------------------------------------------------

	// Setup the parameters for the test
	DB.NodeType = malloc(200 * sizeof *DB.NodeType);  // keep
	DB.MeshFileName = malloc(200 * sizeof *DB.MeshFileName);  // keep
	DB.BasisType = malloc(200 * sizeof *DB.BasisType);  // keep
	DB.TestType = malloc(200* sizeof *DB.TestType);  // keep

	strcpy(DB.NodeType, "GL");
	strcpy(DB.MeshFileName, "/Users/jm-034232/Documents/McGill/Thesis/IGA-DG/meshes/6x1_P1_Deform_InviscidChannel_V4.2.msh");
	strcpy(DB.BasisType, "Polynomial");  // NURBS or Polynomial
	strcpy(DB.TestType, "InviscidChannel"); // PeriodicVortex, InviscidChannel

	// Start reading the initialization file
	printf("COMPLETED READING \n");
	exit(0);

	DB.d = 2;
	DB.Testing = 0;
	DB.numTimeSteps = 100000;
	DB.printSolFreq = 50;
	DB.shapeFuncType = 0;

	DB.dt = 1e-2;
	DB.exit_tol = 1e-5;


	//TODO: Place this in method initialize_test_case_parameters
	// 	in this same file

	// -----------------------------------------------------------------
	//						Setup Flow Properties
	// -----------------------------------------------------------------

	// Flow Properties: (used for the subsonic channel in Philip's Code)

	DB.p_Total = 1.0;
	DB.T_Total = 1.0;
	DB.Rg      = 1.0;
	DB.pBack   = 0.99*DB.p_Total;

	DB.rhoInf = DB.p_Total/(DB.Rg*DB.T_Total);
	DB.pInf   = DB.p_Total;

	DB.MInf   = sqrt(2.0/GM1*(pow((DB.pBack/DB.p_Total),-GM1/GAMMA)-1.0));
	DB.cInf   = sqrt(GAMMA*DB.pInf/DB.rhoInf);


	// Free Memory:
	free(CtrlFile); 

}





