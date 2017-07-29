
#ifndef DG_S_FACE_h__INCLUDED
#define DG_S_FACE_h__INCLUDED

struct S_FACE {
	
	// Properties:
	int P, // The order for the face (P+1 integration points on face)
		Boundary, // 1 = External Boundary, 0 = Internal Face
		VeInd[2], // The two vertices nodes that make up this volume
		fin,
		fout,
		BCType;

	// Geometry:
	double 	*n_fI; // Normal vectors in column major form at integration points.
					// Is the normal outward vector from the in volume.
					//  nx_1 ny_1
					//  nx_2 ny_2

	// Metric Terms:
	double 	*C_fI; // Metric terms at the face integration nodes for this face


	// Structures:
	struct S_VOLUME *VIn, // The pointer to the left volume
					*VOut; // The pointer to the right volume

	struct S_FACE 	*next, // Pointer to the next face in linked list
					*parent;  // Pointer to the previous face in the linked list
	
	
};

#endif // DG_S_FACE_h__INCLUDED

