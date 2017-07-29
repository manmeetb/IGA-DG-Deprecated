

#include "memory_constructors.h"

#include <stdlib.h>
#include <stdio.h>

#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_ELEMENT.h"
#include "S_BC.h"

struct S_ELEMENT *New_ELEMENT(void){

	struct S_ELEMENT *ELEMENT;
	ELEMENT = malloc(sizeof *ELEMENT);

	return ELEMENT;

}

struct S_VOLUME *New_VOLUME(void){

	struct S_VOLUME *VOLUME;
	VOLUME = malloc(sizeof *VOLUME);

	//Structs
	VOLUME->next = NULL;
	VOLUME->parent = NULL;

	return VOLUME;
	
}

struct S_FACE *New_FACE(void){

	struct S_FACE *FACE;
	FACE = malloc(sizeof *FACE);

	//Structs
	FACE->next = NULL;
	FACE->parent = NULL;

	return FACE;
	
}

struct S_BC *New_BC(void){
	struct S_BC *BC;
	BC = malloc(sizeof *BC);

	//Structs
	BC->next = NULL;

	return BC;
}



