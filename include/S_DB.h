

#ifndef DG__S_DB_h__INCLUDED
#define DG__S_DB_h__INCLUDED

struct S_DB {
	// Time
	double time_total;

	char *NodeType, // GL or GLL
		 *BasisType, // Modal or Nodal
		 *MeshFileName, // Path to the mesh file
		 *TestType; // Type of test being completed


	int d, // Dimension of the mesh
		P, // Order of the solution polynomials
		// 0 = No, 1 = All, 2 = Structures, 3 = Operators, 4 = SetupGeometry
		// 5 = Solver Explicit
		Testing, 
		numTimeSteps, // Total number of time steps
		printSolFreq, // Every how many time steps to print solution
		shapeFuncType; // 0 = Polynomial, 1 = NURBS

	double exit_tol;

	// Mesh Information:
	int NV; // Number of Volumes
	
	double dt; // Size of one time step

	// Connectivity Temporary data structures:
	int NGP, // Total number of geometry node points
		NVe; // Total number of vertices

	double	*XYZ_G,	// The matrix, stored in column major form,
					// of all geometry node points.
			*XYZ_Ve; // Vertices in column major form

	int	NGConPerV;	// Number of geometry connectivity ints per
					// element

	int 	*GeoCon, // Connectivity information for geometry nodes
			*VeCon; // Connectivity information for vertices

	// Flow Properties:
	double 	p_Total, T_Total, Rg, pBack, 
			rhoInf, pInf, MInf, cInf;


	// Structs
	struct S_VOLUME  *VOLUME_HEAD; // Pointer to the first volume.
	struct S_FACE *FACE_HEAD; // Pointer to the first face in linked list
	struct S_ELEMENT *ELEMENT; // Pointer to the reference element
	struct S_BC *BC_HEAD;  // Pointer to the structs with the BC information


	// NURBS Basis
	// The xi and eta vector for the reference element using the order
	// of the mesh
	double *xiVector_ref, *etaVector_ref;

};

// Make this database globally accessible.
extern struct S_DB DB;

#endif // DG__S_DB_h__INCLUDED

