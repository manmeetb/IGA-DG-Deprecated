
#include "output_solution.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "S_VOLUME.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "matrix_functions.h"


void output_tecplot(int t){
	/*
	Purpose:
		Output the solution in tecplot format in order to 
		visualize it
	*/

	struct S_VOLUME *VOLUME;
	int i;
	double *XVals, *YVals, *W_Sol_P, *roVec, *ro_uVec, *ro_vVec, *eVec, x, y, 
		ro, u, v, e_tot;
   	FILE *fp;

   	char *filename;
   	filename = malloc(200 * sizeof *filename);
   	char buffer[100];
	sprintf(buffer, "%d", t);

   	strcpy(filename, "../output_tecplot/output");
	strcat(filename,buffer);
	strcat(filename,".dat");

   	fp = fopen(filename, "w+");

   	for(VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next){
   		fprintf(fp,"VARIABLES = \"X\", \"Y\", \"ro\", \"u\", \"v\", \"e_tot\"");
   		fprintf(fp, " \n");
   		fprintf(fp,"ZONE I=%d, J=%d, DATAPACKING=POINT\n",VOLUME->P+1,VOLUME->P+1);

   		XVals = &VOLUME->XYZ_P[0];
   		YVals = &VOLUME->XYZ_P[VOLUME->NvnP];

		W_Sol_P = malloc(VOLUME->NvnP*VOLUME->NVar* sizeof *W_Sol_P);  // free
		mm_CNN(VOLUME->NvnP, VOLUME->NvnG, VOLUME->NVar, DB.ELEMENT->Chi_vP, VOLUME->What, W_Sol_P);

		// Get the solution wHat at the solution points
		roVec = &(W_Sol_P[0]);
		ro_uVec = &(W_Sol_P[1*VOLUME->NvnP]);
		ro_vVec = &(W_Sol_P[2*VOLUME->NvnP]);
		eVec = &(W_Sol_P[3*VOLUME->NvnP]);

		// Loop over all the solution nodes on this volume and print the 
		// the w vector at each solution node.
		for (i=0; i<VOLUME->NvnP; i++){
			x = XVals[i];
			y = YVals[i];

			ro = roVec[i];
			u = ro_uVec[i]/ro;
			v = ro_vVec[i]/ro;
			e_tot = eVec[i]/ro;

			fprintf(fp, "%e %e %e %e %e %e\n", x,y,ro,u,v,e_tot);
		}

		free(W_Sol_P);
   	}
   	fclose(fp);

   	free(filename);

}


void outputWSol(void){

	/*
	Purpose:
		Output the solution at each solution node point
	*/

	struct S_VOLUME *VOLUME;
	int i;
	double *XVals, *YVals, *roVec, *ro_uVec, *ro_vVec, *eVec, x, y, 
			ro, ro_u, ro_v, e;
	double *W_Sol;
	
	for (VOLUME = DB.VOLUME_HEAD; VOLUME; VOLUME = VOLUME->next) {

		// Solution node locations. Recall that values are stored in 
		// column major form.	
		XVals = &(VOLUME->XYZ_S[0]);
		YVals = &(VOLUME->XYZ_S[1*VOLUME->NvnS]);

		W_Sol = malloc(VOLUME->NvnS*VOLUME->NVar* sizeof *W_Sol);

		mm_CNN(VOLUME->NvnS, VOLUME->NvnG, VOLUME->NVar, DB.ELEMENT->Chi_vS, VOLUME->What, W_Sol);

		// Get the solution wHat at the solution points
		roVec = &(W_Sol[0]);
		ro_uVec = &(W_Sol[1*VOLUME->NvnS]);
		ro_vVec = &(W_Sol[2*VOLUME->NvnS]);
		eVec = &(W_Sol[3*VOLUME->NvnS]);

		// Loop over all the solution nodes on this volume and print the 
		// the w vector at each solution node.
		for (i=0; i<VOLUME->NvnS; i++){
			x = XVals[i];
			y = YVals[i];

			ro = roVec[i];
			ro_u = ro_uVec[i];
			ro_v = ro_vVec[i];
			e = eVec[i];

			printf("x: %.14f, y: %.14f \n", x, y);
			printf("    [ro, ro*u, ro*v, e] = [%.14f, %.14f, %.14f, %.14f] \n", ro, ro_u, ro_v, e);

		}

		free(W_Sol);
	}
}



