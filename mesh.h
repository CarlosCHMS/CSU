#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct
{

    int Nelem;
    int** elem;
    char name[50];

} MESHBC;

typedef struct
{

    int Ndim;
    int Nelem;
    int Np;
    int Nmark;
    
    int** elem;
    
    double** p;
    
    MESHBC** bc;

} MESH;

char meshGetWord(FILE* ff, char* s);
 
MESHBC* meshBCread(FILE* ff);

MESH* meshInit(char* fileName);

void meshPrintBC(MESHBC* bc);

void meshPrint(MESH* mesh);

void meshBCFree(MESHBC* bc);

void meshFree(MESH* mesh);


