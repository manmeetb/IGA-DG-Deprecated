#ifndef DG__boundary_conditions_h__INCLUDED
#define DG__boundary_conditions_h__INCLUDED


void boundary_SlipWall(double *WL, double *nL, double *WB, 
	double *FB, double *GB, int n);

void boundary_BackPressure(double *WL, double *nL, double *WB, 
	double *FB, double *GB, int n);

void boundary_Total_TP(double *WL, double *nL, double *WB, 
	double *FB, double *GB, int n);

#endif // DG__boundary_conditions_h__INCLUDED