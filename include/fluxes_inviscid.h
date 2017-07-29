#ifndef DG__fluxes_inviscid_h__INCLUDED
#define DG__fluxes_inviscid_h__INCLUDED


void flux_LF(double *WIn, double *WOut, double *FIn, double *FOut, double *GIn,
	double *GOut, double *FComm, double *nL, int P, int Neq);

#endif //DG__fluxes_inviscid_h__INCLUDED