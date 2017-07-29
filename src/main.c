
#include <stdlib.h>
#include <stdio.h>

#include "setup_mesh.h"
#include "initialization.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "memory_free.h"
#include "setup_operators.h"
#include "setup_geometry.h"
#include "matrix_functions.h"
#include "initialize_test_case.h"
#include "output_solution.h"
#include "solver_explicit.h"
#include "compute_errors.h"

struct S_DB   DB;

/*
Purpose:
	Solve the 2D and 3D Euler Equations using the Discontinuous 
	Galerkin method. 
*/


/*
---------------------------------------------------------------------

Philip Code:
- Get test case working on new branch of Philip's code

- test_integration_Euler
	(Done)- Only do the convergence order testing.
		- Note: When value is one in ctrl file, means perform that test. Otherwise
			it will be 0.
	- Compare certain aspects between the polynomial and NURBS basis to ensure
		everything has been implemented correctly.
		- Figure out what should be compared.

- Implement NURBS mesh to be read by Philip's code. 
	- Create a gmsh file with the same skeleton. Then, if IGA is 
		an option when being run, overwrite the geometry node points
		for each element with the geometry node points in a seperate
		mesh.
		- This way, connectivity does not have to be set since the gmsh file
			does this anyways.
		- The gmsh reader will still work the same, it will just have the 
			wring XYZ vector for each volume which will then be overwritten.

- Implement NURBS basis functions into Philip's code
	- First, study how the the polynomial basis functions are set up.
	- Using this information, set up the chi operators for the code.
	- Set the relevant interpolation operators. 
		- For instance, interpolating from geometry to any node
			doesn't make sense since there are no geometry nodes.
		- Will need to check which interpolation operators are 
			being used at each step.

- Check how Bezier implemented in Philip Code
	- Modify setup_ELEMENT_operators for the line (tensor product) 
		and these should transfer over to the quad elements.
	- Modify the setup_geom_factors to use the gradient of the 
		basis at the points and not the interpolation of metric terms
	- Modify setup_normals 

---------------------------------------------------------------------

---------------------------------------------------------------------

Code: 

(Done) - Add proper Makefile to the code
	- Follow Philip's format

- Add a ctrl file to be read by the code.
	- Should have a reference to the mesh file as well as time step
		and what type of time stepping to use.
	- Have a data file generated with as its perfix the name of the ctrl file.
	- Code will now take as command line argument the ctrl file name 
	- Output all error the same way as CPR code now and in the same format
		so that orders test can be used on it.

- Perform a polynomial and B Splines convergence study for coarsest 
	mesh in inviscid channel case.

- Adjust orders to work with finding error for square undeformed case
	without having to convect complete distance.
	- Perform orders study for both B-Spline and Polynomial basis

- Create mesh generator for B Splines
	- Code in python
	- Create the coarsest mesh and then use Bezier extraction
		to generate the bezier elements that have the exact same
		geometry.
		- Use numpy for the Bezier extraction operators
		- Use matplotlib for plotting final mesh

- Adjust orders test to now work by comparing vortex solution at any 
	time based on how much the vortex has convected.
	- Check the orders on the undeformed mesh with this case.
	- All that is needed is to modify how errors are computed for
		vortex case.

- Compute metrics at the face integration and geometry nodes.
	- That is, no longer interpolate the face metric terms

- Add negative pressure check

- Adjust the mesh generator for the polynomial to take command line
	arguments for generating the mesh in the case of making
	it automatically.


- Complete explicit time stepping (RK4) 
	- Make sure results with DPG code are identical

- Make NURBS Mesh Generator
	- Make it possible to visalize the mesh.
	- Make the Gaussian bump case.

---------------------------------------------------------------------

- Check:
	(Done) - Slip Wall Boundary Condition Algorithm
	(Done) - Total TP Boundary Condition Algorithm
	(Done) - Back Pressure Boundary Condition Algorithm

(Done) - Load the mesh manually for the coarsest case from Philip's 
	DPG code. 
	(Done) - Create a python script or modify the current one to be able
		to output the exact same mesh.
	(Done)- Load the new mesh into the code and check its convergence.

- Compare properties of solution between the DPG and DG code. Compare:
	(Done) - Volume metric terms at integration points
	- Inverse mass matrix
	- RHS_VOLUME contribution
	- RHS_FACE contribution


- Study Reimann Invariants and BC setting
- Derive Euler Equations properly
- Compute characteristics of the Euler Equations in 2D
	- Compute Reimann invariant formulas
	
*/

int main(int nargc, char **argv){

	printf("Initialization \n");
	initialization(nargc, argv);

	printf("Setup mesh \n");
	setup_mesh();

	printf("Setup operators \n");
	setup_operators();

	printf("Setup geometry \n");
	setup_geometry();

	printf("Initialize test operators \n");
	initialize_test_case();

	// Output the initial solution (function should take in
	// maybe a string)
	output_tecplot(0);

	printf("Start Solver");
	solver_explicit();

	printf("Compute Errors");
	// compute_errors_global();

	output_tecplot(DB.numTimeSteps);

	memory_free();

	return 0;

}



