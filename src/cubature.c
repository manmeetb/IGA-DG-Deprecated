
#include "cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void cubature_literature(int P, char *NodeType, double *Nodes, double *Weights){
	/*
	 *	Purpose:
	 *		Return the 1D nodes and weights for the chosen node type
	 */

	int i;

	if(strstr(NodeType, "GL")){
		
		if (P+1 == 2){
			double xi[2] = {-0.5773502691896257, 0.5773502691896257};
			double wi[2] = { 1.0000000000000000, 1.0000000000000000};

			for(i=0; i<(P+1); i++){
				Nodes[i] = xi[i];
				Weights[i] = wi[i];
			}
		} else if(P+1 == 3){

			double xi[3] = {-0.7745966692414834, 0.0000000000000000, 0.7745966692414834};
			double wi[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

			for(i=0; i<(P+1); i++){
				Nodes[i] = xi[i];
				Weights[i] = wi[i];
			}
		} else if(P+1 == 4){

			double xi[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
			double wi[4] = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};

			for(i=0; i<(P+1); i++){
				Nodes[i] = xi[i];
				Weights[i] = wi[i];
			}
		} else if(P+1 == 5){

			double xi[5] = {-0.9061798459386640, -0.5384693101056831, 0.0000000000000000, 0.5384693101056831, 0.9061798459386640};
			double wi[5] = { 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

			for(i=0; i<(P+1); i++){
				Nodes[i] = xi[i];
				Weights[i] = wi[i];
			}
		} else {
			printf("UNSUPPORTED CUBATURE N : %d \n", P+1);
			exit(1);
		}

	} else{
		printf("UNSUPPORTED CUBATURE TYPE \n");
		exit(1);
	}

}

