#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"utils.h"
#include"mesh.h"


char meshGetWord(FILE* ff, char* s)
{

    int stop;
    int ii;
    char c;
    
    stop = 0;
    ii = 0;
    while(!stop)
    {
        c = fgetc(ff);
                
        if(c == ' ')
        {
            stop = 1;
        }
        else if(c == '\n')
        {
            stop = 1;
        }
        else
        {
            s[ii] = c;
            ii++;
        }
        
    }
    s[ii] = '\0';

    return c;

} 

MESHBC* meshBCread(FILE* ff) 
{
    int jj;
    char s[100];
    MESHBC* bc = malloc(sizeof(MESHBC));
    
    meshGetWord(ff, s);
    meshGetWord(ff, bc->name);
    meshGetWord(ff, s);
    if(strcmp(s, "MARKER_ELEMS=")==0)
    {
        meshGetWord(ff, s);
        bc->Nelem = atoi(s);
    }    

    bc->elem = tableMallocInt(bc->Nelem, 2);
    jj = 0;
    while(jj<bc->Nelem)
    {
        meshGetWord(ff, s);
        meshGetWord(ff, s);
        bc->elem[jj][0] = atoi(s);
        meshGetWord(ff, s);        
        bc->elem[jj][1] = atoi(s);
        meshGetWord(ff, s);
        jj++;
    }

    return bc;

}

MESH* meshInit(char* fileName)
{

    MESH* mesh = malloc(sizeof(MESH));
    FILE* ff = fopen(fileName, "r");
    char s[100];
    int ii, jj;
    
    meshGetWord(ff, s);
    if(strcmp(s, "NDIME=")==0)
    {
        meshGetWord(ff, s);
        mesh->Ndim = atoi(s);
    }

    meshGetWord(ff, s);
    if(strcmp(s, "NELEM=")==0)
    {
        meshGetWord(ff, s);
        mesh->Nelem = atoi(s);
    }

    mesh->elem = tableMallocInt(mesh->Nelem, 3);
    
    ii = 0;
    while(ii<mesh->Nelem)
    {
        meshGetWord(ff, s);
        meshGetWord(ff, s);
        mesh->elem[ii][0] = atoi(s);
        meshGetWord(ff, s);        
        mesh->elem[ii][1] = atoi(s);
        meshGetWord(ff, s);        
        mesh->elem[ii][2] = atoi(s); 
        meshGetWord(ff, s);
        ii++;
    }

    meshGetWord(ff, s);
    if(strcmp(s, "NPOIN=")==0)
    {
        meshGetWord(ff, s);
        mesh->Np = atoi(s);
    }
    
    mesh->p = tableMallocDouble(mesh->Np, 2);
    
    ii = 0;
    while(ii<mesh->Np)
    {
        meshGetWord(ff, s);
        mesh->p[ii][0] = strtod(s, NULL);
        meshGetWord(ff, s);        
        mesh->p[ii][1] = strtod(s, NULL);
        meshGetWord(ff, s);        
        ii++;
    }    
    
    meshGetWord(ff, s);
    if(strcmp(s, "NMARK=")==0)
    {
        meshGetWord(ff, s);
        mesh->Nmark = atoi(s);
    }    
    
    mesh->bc = (MESHBC**)malloc(mesh->Nmark*sizeof(MESHBC*));
    
    ii = 0;
    while(ii<mesh->Nmark)
    {
        mesh->bc[ii] = meshBCread(ff);
        ii++;
    }
    
    fclose(ff);

    return mesh;

}

void meshPrintBC(MESHBC* bc)
{

    printf("%s\n", bc->name);
    printf("%i\n", bc->Nelem);

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        printf("%i, %i\n", bc->elem[ii][0], bc->elem[ii][1]);
    }

}

void meshPrint(MESH* mesh)
{

    int ii;

    printf("%i\n", mesh->Nelem);

    for(ii=0; ii<mesh->Nelem; ii++)
    {
        printf("%i, %i, %i\n", mesh->elem[ii][0], mesh->elem[ii][1], mesh->elem[ii][2]);
    }

    printf("%i\n", mesh->Np);
    
    for(ii=0; ii<mesh->Np; ii++)
    {
        printf("%.10e, %.10e\n", mesh->p[ii][0], mesh->p[ii][1]);
    }

    printf("%i\n", mesh->Nmark);

    for(int ii=0; ii<mesh->Nmark; ii++)
    {
        meshPrintBC(mesh->bc[ii]);
    }

}

void meshBCFree(MESHBC* bc)
{

    tableFreeInit(bc->elem, bc->Nelem);
    free(bc->name);

}

void meshFree(MESH* mesh)
{

    tableFreeInit(mesh->elem, mesh->Nelem);
    tableFreeDouble(mesh->p, mesh->Np);
    
    for(int ii; ii<mesh->Nmark; ii++)
    {

        //meshBCFree(mesh->bc[ii]);
        //free(mesh->bc[ii]);

    }

    //free(mesh);

}
    


