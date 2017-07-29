#ifndef DG__memory_destructors_h__INCLUDED
#define DG__memory_destructors_h__INCLUDED

#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_ELEMENT.h"
#include "S_BC.h"

void memory_destructor_V (struct S_VOLUME *VOLUME);
void memory_destructor_F (struct S_FACE *FACE);
void memory_destructor_E (struct S_ELEMENT *ELEMENT);
void memory_destructor_BC (struct S_BC *BC);

#endif // DG__memory_destructors_h__INCLUDED