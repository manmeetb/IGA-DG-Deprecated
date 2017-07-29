

#include "bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *	Purpose:
 *		Compute the value of the chosen basis function at the
 *		given points on the computational domain
 *
 *	Comments:
 *		- First, only the polynomial Lagrange Basis will be implemented.
 *		- Once this is complete, the NURBS basis will be added.
 *
 *	Notation:
 *	
 *
 *	References:
 */


static double Lagrange_iPrime(int iBasis, double x, double *nodeList, int numNodes);
static double Lagrange_i(int iBasis, double x, double *nodeList, int numNodes);
double basis_TP_Lagrange_2D(int P, int basis_i, int basis_j, double *nodeLocations, 
	double xi, double eta);
double *basis_TP_Lagrange_2D_Grad(int P, int basis_i, int basis_j, double *nodeLocations, 
	double xi, double eta);

static double N_ip(int iBasis, int P, double xi, double *xiVector);
static double N_ipPrime(int iBasis, int P, double xi, double *xiVector);
double basis_TP_NURBS_2D(int P, int basis_i, int basis_j, double *xiVector, double *etaVector, 
	double xi, double eta);
double *basis_TP_NURBS_2D_Grad(int P, int basis_i, int basis_j, double *xiVector, double *etaVector, 
	double xi, double eta);



static double Lagrange_iPrime(int iBasis, double x, double *nodeList, int numNodes){
	/*
	Compute the value of the ith Lagrange basis function (1D) evaluated at
	point x. 3 cases for computing the derivative:
		1) If x is an arbitrary point on the domain.
		2) If x is one of the node points (but not the ith node where basis is 1)
		3) If the x is at the ith node point

	:param iBasis : Which basis function is being used
	:param x : The value along the real number line to evaluate the basis function at
	:param *nodeList : The array of the location of the node points along the 1D number line.
	:param numNodes : The number of nodes

	:return value of the derivative of the basis function at x
	*/

	int i, j, k, commonIndex;
	double factor, w_i, w_j, val;

	commonIndex = -1;
	// Check what case this is
	for(i=0; i<numNodes; i++){
		if (fabs(x - nodeList[i]) < 1E-7){
			commonIndex = i;
			break;
		} 
	}

	if(commonIndex == -1){

		// Case 1:

		factor = 0.0;

		for(i=0; i<numNodes; i++){
			if(i!=iBasis){
				factor = factor + 1./(x-nodeList[i]);
			}
		}

		return factor*Lagrange_i(iBasis, x, nodeList, numNodes);

	} else{

		if (commonIndex != iBasis){

			// Case 2:

			j = iBasis;
			i = commonIndex;

			w_i = 1.0;
			for (k=0; k<numNodes; k++){
				if (i != k){
					w_i = w_i*(nodeList[i] - nodeList[k]);
				}
			}
			w_i = 1./w_i;

			w_j = 1.;
			for(k=0; k<numNodes; k++){
				if(j != k){
					w_j = w_j*(nodeList[j] - nodeList[k]);
				}
			}
			w_j = 1./w_j;

			return (w_j/w_i)/(nodeList[i] - nodeList[j]);

		} else{

			// Case 3:

			val = 0;

			for(j=0; j<numNodes; j++){
				if(j != commonIndex){
					val = val + Lagrange_iPrime(j, x, nodeList, numNodes);
				}
			}

			return val*-1.;

		}

	}

}

static double Lagrange_i(int iBasis, double x, double *nodeList, int numNodes){
	/*
	Compute the value of the ith 1D Lagrange basis function evaluated
	at point x.

	:param iBasis : Which basis function is being used
	:param x : The value along the real number line to evaluate the basis function at
	:param *nodeList : The array of the location of the node points along the 1D number line.
	:param numNodes: The number of nodes

	:return value of the basis function at x
	*/

	int i;
	double numer, denom;

	numer = 1.0;
	denom = 1.0;

	for(i=0; i<numNodes; i++){
		if(i!=iBasis){
			numer = numer*(x-nodeList[i]);
			denom = denom*(nodeList[iBasis] - nodeList[i]);
		}
	}

	return (numer/denom);

}

double basis_TP_Lagrange_2D(int P, int basis_i, int basis_j, double *nodeLocations, 
	double xi, double eta){

	/*
	Compute the lagrange shape function (2D) at the xi, eta value on the computational 
	domain. Xi and eta are the coordinate axes on the computational domain.

	:param P : The order of the shape functions
	:param basis_i : The i index of the basis 
	:param basis_j : The j index of the basis 
	:param *nodeLocations : The matrix of the value of the xi,eta points at the nodes.
	:param xi : The xi value at which to evaluate chosen basis
	:param eta : The eta value at which to evaluate chosen basis function

	:return value of the shape function at the chosen point
	*/

	double *nodeVectorXi, *nodeVectorEta, *nodeLocations_xi, *nodeLocations_eta;
	int i;

	// This is a tensor product element. Find the node location 
	// along each coordinate axis (xi and eta)
	nodeVectorXi = malloc((P+1)* sizeof *nodeVectorXi); //free
	nodeVectorEta = malloc((P+1)* sizeof *nodeVectorEta); //free

	// Recall: The nodeLocations matrix is stored with all eta = -1, ... It is 
	// also in column major form.
	nodeLocations_xi = &nodeLocations[0];
	nodeLocations_eta = &nodeLocations[(P+1)*(P+1)];

	for(i=0; i<(P+1); i++){
		nodeVectorXi[i] = nodeLocations_xi[i];
		nodeVectorEta[i] = nodeLocations_eta[i*(P+1)];
	}

	return Lagrange_i(basis_i, xi, nodeVectorXi, P+1)*Lagrange_i(basis_j, eta, nodeVectorEta, P+1);

}

double *basis_TP_Lagrange_2D_Grad(int P, int basis_i, int basis_j, double *nodeLocations, 
	double xi, double eta){

	/*
	Compute the gradient of the lagrange shape function (2D) at the xi, eta value on the \
	computational domain. Xi and eta are the coordinate axes on the computational domain.

	:param P : The order of the shape functions
	:param basis_i : The i index of the basis 
	:param basis_j : The j index of the basis 
	:param *nodeLocations : The matrix of the value of the xi,eta points at the nodes.
	:param xi : The xi value at which to evaluate chosen basis
	:param eta : The eta value at which to evaluate chosen basis function

	:return value of the shape function at the chosen point
	*/

	double del_by_del_xi, del_by_del_eta, *grad;
	double *nodeVectorXi, *nodeVectorEta, *nodeLocations_xi, *nodeLocations_eta;
	int i;

	grad = malloc(2* sizeof *grad);

	// This is a tensor product element. Find the node location 
	// along each coordinate axis (xi and eta)
	nodeVectorXi = malloc((P+1)* sizeof *nodeVectorXi); //free
	nodeVectorEta = malloc((P+1)* sizeof *nodeVectorEta); //free

	// Recall: The nodeLocations matrix is stored with all eta = -1, ... It is 
	// also in column major form.
	nodeLocations_xi = &nodeLocations[0];
	nodeLocations_eta = &nodeLocations[(P+1)*(P+1)];

	for(i=0; i<(P+1); i++){
		nodeVectorXi[i] = nodeLocations_xi[i];
		nodeVectorEta[i] = nodeLocations_eta[i*(P+1)];
	}

	del_by_del_xi = Lagrange_i(basis_j, eta, nodeVectorEta, P+1)*Lagrange_iPrime(basis_i,
		 xi, nodeVectorXi, P+1);
	del_by_del_eta = Lagrange_i(basis_i, xi, nodeVectorXi, P+1)*Lagrange_iPrime(basis_j,
		 eta, nodeVectorEta, P+1);

	grad[0] = del_by_del_xi;
	grad[1] = del_by_del_eta;

	return grad; 

}







// B Spline Basis Functions
static double N_ipPrime(int i, int p, double xi, double *xiVector){
	/*
	Compute the value of the ith B Spline basis function derivative (1D) evaluated at
	point xi on the mapped domain.

	:param iBasis : Which basis function is being used
	:param x : The value along the real number line to evaluate the basis function at
	:param *nodeList : The array of the location of the node points along the 1D number line.
	:param numNodes : The number of nodes

	:return value of the derivative of the basis function at x
	*/


	double xi_i, xi_iPlusP, xi_iPlusPPlus1, xi_iPlus1;
	double numer1, denom1, term1, numer2, denom2, term2;

	xi_i = xiVector[i];
	xi_iPlusP = xiVector[i+p];
	xi_iPlusPPlus1 = xiVector[i+p+1];
	xi_iPlus1 = xiVector[i+1];

	numer1 = p*N_ip(i, p-1, xi, xiVector);
	denom1 = (xi_iPlusP - xi_i);

	if(fabs(denom1) < 1E-14){
	    term1 = 0;
	} else{
	    term1 = numer1/denom1;
	}

	numer2 = p*N_ip(i+1, p-1, xi, xiVector);
	denom2 = (xi_iPlusPPlus1 - xi_iPlus1);
	if(fabs(denom2) < 1E-14){
	    term2 = 0;
	} else{
	    term2 = numer2/denom2;
	}

	return term1 - term2;

}

static double N_ip(int i, int p, double xi, double *xiVector){
	/*
	Compute the value of the ith 1D NURBS basis function evaluated
	at point xi. The NURSB domain is from -1 to 1.

	:return value of the basis function at x
	*/

	double xi_i, xi_iPlus1, xi_iPlusP, xi_iPlusPPlus1;

	double num1, denom1, num2, denom2, 
		term1, term2;  // Numerator and denominator of B-Spline basis


	if (p == 0){
		// Base Case
		
		if(xi < xiVector[i+1] && xi >= xiVector[i]){
			return 1.;
		} else{
			return 0.;
		}

	} else{
		// Recursive Case

		xi_i = xiVector[i];
    	xi_iPlus1 = xiVector[i+1];
    	xi_iPlusP = xiVector[i+p];
    	xi_iPlusPPlus1 = xiVector[i+p+1];

    	// - First Term
    	num1 = (xi - xi_i)*(N_ip(i,p-1, xi, xiVector));
    	denom1 = xi_iPlusP - xi_i;

    	if (fabs(num1) < 1E-14 && fabs(denom1) < 1E-14){
    		term1 = 0;
    	} else{
    		term1 = num1/denom1;
    	}

    	// - Second Term
    	num2 = (xi_iPlusPPlus1-xi)*N_ip(i+1,p-1,xi,xiVector);
    	denom2 = xi_iPlusPPlus1 - xi_iPlus1;
    	if(fabs(num2) < 1E-14 && fabs(denom2) < 1E-14){
       		term2 = 0;
    	}
    	else{
        	term2 = num2/denom2;
    	}

    	return term1 + term2;
	}

}

double basis_TP_NURBS_2D(int P, int basis_i, int basis_j, double *xiVector, double *etaVector, 
	double xi, double eta){

	/*
	Compute the NURBS shape function (2D) at the xi, eta value on the computational 
	domain. Xi and eta are the coordinate axes on the computational domain.

	:param P : The order of the shape functions
	:param basis_i : The i index of the basis 
	:param basis_j : The j index of the basis 
	:param *xiVector : The knot vector along the xi direction
	:param *etaVector : The knot vector along the eta direction	
	:param xi : The xi value at which to evaluate chosen basis
	:param eta : The eta value at which to evaluate chosen basis function

	:return value of the shape function at the chosen point
	*/

	return N_ip(basis_i, P, xi, xiVector)*N_ip(basis_j, P, eta, etaVector);

}

double *basis_TP_NURBS_2D_Grad(int P, int basis_i, int basis_j, double *xiVector, double *etaVector, 
	double xi, double eta){

	/*
	Compute the gradient of the NURBS shape function (2D) at the xi, eta value on the \
	computational domain. Xi and eta are the coordinate axes on the computational domain.

	:param P : The order of the shape functions
	:param basis_i : The i index of the basis 
	:param basis_j : The j index of the basis 
	:param *xiVector : The knot vector along the xi direction
	:param *etaVector : The knot vector along the eta direction	
	:param xi : The xi value at which to evaluate chosen basis
	:param eta : The eta value at which to evaluate chosen basis function

	:return value of the gradient of the shape function at the chosen point
	*/

	double *grad;
	double del_by_del_xi, del_by_del_eta;

	grad = malloc(2* sizeof *grad);

	del_by_del_xi = N_ip(basis_j, P, eta, etaVector)*N_ipPrime(basis_i, P, xi, xiVector);
	del_by_del_eta = N_ip(basis_i, P, xi, xiVector)*N_ipPrime(basis_j, P, eta, etaVector);

	grad[0] = del_by_del_xi;
	grad[1] = del_by_del_eta;

	return grad;

}




