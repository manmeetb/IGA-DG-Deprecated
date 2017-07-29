

#ifndef DG__S_VOLUME_h__INCLUDED
#define DG__S_VOLUME_h__INCLUDED

struct S_VOLUME {

	// Properties:
	int index, P, d, update; // The index of this volume in the mesh (from 
				// connectivity file) starting from 0

	// Geometry Information:
	int VeInd[4], // (Ve)rtices (Ind)eces in mesh vertices node list
		FIndFound[4]; // (F)aces (I)eces Found. Will be 1 if that face has been found
					// before and 0 otherwise

	int NvnG, // (N)umber of (v)olume (n)odes for (g)eometry
		NvnP, // (N)umber of (v)olume (n)odes for (p)lotting
		NvnS, // (N)umber of (v)olume (n)odes for (s)olution
		NfnI; 

	double	*XYZ, // Geometry Node Points matrix in column major form 
			*XYZ_S, // Solution Node Points matrix in column major form
			*XYZ_P; // Plotting Node Points matrix in column major form

	// structs:
	//	- Linked List of volumes for the mesh
	struct S_VOLUME *next, // Pointer to next volume in the L List
					*parent; // Pointer to parent in L List

	struct S_FACE **FACE;


	//	Metric Terms:
	/*
	C = y_eta  -y_xi
		-x_eta	x_xi

	J = x_xi * y_eta - y_xi * x_eta

	*/
	double 	*detJV_vS, // Jacobian at solution nodes as a vector
			*C_vS; // metric terms at solution nodes at a matrix. At
					// ith solution node (row of matrix) ordering 
					// of column elements is [C11, C12, C21, C22]

	// Solution Data Structures:
	int 	NVar;
	double 	*What; 	// Modal coefficients. Ordering is in same ordering in which
					// solution points are arranged in XYZ_vS.

	// Solving: 
	// 	- Data structures used during the explicit solving.
	double 	*RHS, *RES, *RHS_VOL, *RHS_FACE, *MInv,
			*F_Comm; // matrix of F_comm values at all the integration nodes
	
};


#endif // DG__S_VOLUME_h__INCLUDED

