

#include "memory_destructors.h"

#include <stdlib.h>
#include <stdio.h>

#include "S_DB.h"


void memory_destructor_E (struct S_ELEMENT *ELEMENT){
	free(ELEMENT);
}

void memory_destructor_V (struct S_VOLUME *VOLUME){
	free(VOLUME);
}

void memory_destructor_F (struct S_FACE *FACE){
	free(FACE);
}

void memory_destructor_BC (struct S_BC *BC){
	free(BC);
}
