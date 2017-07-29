#ifndef DG_bases_h__INCLUDED
#define DG_bases_h__INCLUDED

double basis_TP_Lagrange_2D(int P, int basis_i, int basis_j, double *nodeLocations, 
	double xi, double eta);

double *basis_TP_Lagrange_2D_Grad(int P, int basis_i, int basis_j, double *nodeLocations, 
	double xi, double eta);

double basis_TP_NURBS_2D(int P, int basis_i, int basis_j, double *xiVector, double *etaVector, 
	double xi, double eta);

double *basis_TP_NURBS_2D_Grad(int P, int basis_i, int basis_j, double *xiVector, double *etaVector, 
	double xi, double eta);

#endif