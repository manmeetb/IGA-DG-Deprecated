#ifndef DG_S_BC_h__INCLUDED
#define DG_S_BC_h__INCLUDED

/*
Struct: S_BC

Purpose:
	Holding the boundary condition information for the mesh
*/

struct S_BC {

	int BCType;

	int *BC_Con, // Matrix with the BC connectivity information (row major form)
		nBC_Con_row, nBC_Con_col;  // number of rows and cols for the matrix

	struct S_BC *next;

};


#endif // DG_S_BC_h__INCLUDED