
#ifndef DG_S_ELEMENT_h__INCLUDED
#define DG_S_ELEMENT_h__INCLUDED

struct S_ELEMENT {
	// Structure for Reference Element and holding all
	// operators of the reference element

	/*
	Notation:
	I_(1)(2)_(3)(4) : (I)nterpolation operator from (1) nodes of type (2) to (4) nodes of type (5)
 		(1/3): (v)olume, (f)ace
		(2/4): (P)lotting, (G)eometry, (I)ntegration, (S)olution
	*/

	// Cubature:
	double	*nodes_xi, // Location of 1D cubature points
			*nodes_wi; // Weights at each node

	// Properties:
	int NvnG, NvnS, NvnP, NfnI, //Note: (N)umber of (f)ace (n)odes for (I)ntegration for one face
		d, P;

	// Operators:
	// - Chi Operators
	double 	*Chi_vS,  // Vandermonde
			*ChiInv_vS,  // Inverse Vandermonde
			*Chi_vG,
			*Chi_vP,
			*Chi_fI;

	// - Derivative Matrices: (on computational domain)
	double 	*GradChi_vS_xi, // Derivative of basis function evalutaed at solution nodes.
			*GradChi_vS_eta,
			*GradChi_fI_xi,  // Derivative of basis function evaluated at face integration
			*GradChi_fI_eta; // nodes

	// - Interpolation Operators:
	double	*I_vG_vP,  // Not implemented yet
			*I_vG_vS,  
			*I_vS_vG,
			*I_vS_fI;

	// Geometry:
	double	*XiEtaZeta_S, // coordinates Comp. Elem. Solution Pts.
			*XiEtaZeta_G, // coordinates Comp. Elem Geometry Pts.
			*XiEtaZeta_P, // coordinates Comp. Elem Plotting Pts.
			*XiEtaZeta_F;

	// NURBS Basis
	// The xi and eta vector for the reference element using the order
	// of the mesh
	double *xiVector, *etaVector;

};

#endif // DG_S_ELEMENT_h__INCLUDED


