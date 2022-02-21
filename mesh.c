#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"readTables.h"
#include"mesh.h"

MESH* meshInit(char* fileName)
{

    TABLELIST* tl = fReadTables(fileName);
    MESH* mesh = malloc(sizeof(MESH));
    
    mesh->elem = tl->tables[tl->N-2]->values;
    mesh->p = tl->tables[tl->N-1]->values;

    mesh->Nelem = tl->tables[tl->N-2]->Nrow;
    mesh->Np = tl->tables[tl->N-1]->Nrow;   
    
    return mesh;

}


