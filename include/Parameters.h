
#ifndef DG__Parameters_h__INCLUDED
#define DG__Parameters_h__INCLUDED


// Common variables
#define PI    3.1415926535897932
#define GAMMA 1.4
#define GM1   0.4
#define GM3   -1.6

// Boundary conditions
#define BC_RIEMANN        1
#define BC_SLIPWALL       2
#define BC_BACKPRESSURE   3
#define BC_TOTAL_TP       4
#define BC_SUPERSONIC_IN  5
#define BC_SUPERSONIC_OUT 6

#define BC_NOSLIP_T         7
#define BC_NOSLIP_ADIABATIC 8

#define BC_DIRICHLET    11
#define BC_NEUMANN      12

#define BC_INFLOW       13
#define BC_OUTFLOW      14

#define BC_PERIODIC		15
#define BC_INTERNAL		16

// Solution
#define EPS        1.0e-15

#endif // DG__Parameters_h__INCLUDED