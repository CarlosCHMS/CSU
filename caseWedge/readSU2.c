#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct
{

    int Nelem;
    int** elem;

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

char getWord(FILE* ff, char* s)
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
    
double** tableMallocDouble(int Nrow, int Ncol)
{

    double** M = malloc(Nrow*sizeof(double*));

    for(int ii=0; ii<Nrow; ii++)
    {    
        M[ii] = malloc(Ncol*sizeof(double));
    }

    return M;

}    
 
int** tableMallocInt(int Nrow, int Ncol)
{

    int** M = malloc(Nrow*sizeof(int*));

    for(int ii=0; ii<Nrow; ii++)
    {    
        M[ii] = malloc(Ncol*sizeof(int));
    }

    return M;

} 
 
MESHBC* readBC(FILE* ff)
{
    int jj;
    int s[100];
    MESHBC* bc = malloc(sizeof(MESHBC));
    
    getWord(ff, s);
    getWord(ff, s);
    getWord(ff, s);
    if(strcmp(s, "MARKER_ELEMS=")==0)
    {
        getWord(ff, s);
        bc->Nelem = atoi(s);
    }    

    bc->elem = tableMallocInt(bc->Nelem, 2);
    jj = 0;
    while(jj<bc->Nelem)
    {
        getWord(ff, s);
        getWord(ff, s);
        bc->elem[jj][0] = atoi(s);
        getWord(ff, s);        
        bc->elem[jj][1] = atoi(s);
        jj++;
    }

    return bc;

}
    
int main()
{

    MESH* mesh = malloc(sizeof(MESH));
    FILE* ff = fopen("wedge.su2", "r");
    char s[100];
    int ii, jj;
    
    getWord(ff, s);
    if(strcmp(s, "NDIME=")==0)
    {
        getWord(ff, s);
        mesh->Ndim = atoi(s);
    }

    getWord(ff, s);
    if(strcmp(s, "NELEM=")==0)
    {
        getWord(ff, s);
        mesh->Nelem = atoi(s);
    }

    mesh->elem = tableMallocInt(mesh->Nelem, 3);
    
    ii = 0;
    while(ii<mesh->Nelem)
    {
        getWord(ff, s);
        getWord(ff, s);
        mesh->elem[ii][0] = atoi(s);
        getWord(ff, s);        
        mesh->elem[ii][1] = atoi(s);
        getWord(ff, s);        
        mesh->elem[ii][2] = atoi(s); 
        getWord(ff, s);
        ii++;
    }

    getWord(ff, s);
    if(strcmp(s, "NPOIN=")==0)
    {
        getWord(ff, s);
        mesh->Np = atoi(s);
    }
    
    mesh->p = tableMallocDouble(mesh->Np, 2);
    
    ii = 0;
    while(ii<mesh->Np)
    {
        getWord(ff, s);
        mesh->p[ii][0] = strtod(s, NULL);
        getWord(ff, s);        
        mesh->p[ii][1] = strtod(s, NULL);
        getWord(ff, s);        
        ii++;
    }    
    
    getWord(ff, s);
    if(strcmp(s, "NMARK=")==0)
    {
        getWord(ff, s);
        mesh->Nmark = atoi(s);
    }    
    
    mesh->bc = (MESHBC**)malloc(mesh->Nmark*sizeof(MESHBC*));
    
    ii = 0;
    while(ii<mesh->Nmark)
    {
        mesh->bc[ii] = readBC(ff);
        ii++;
    }     
    
    printf("%i, %i, %i, %i\n", mesh->Ndim, mesh->Nelem, mesh->Np, mesh->Nmark);
    
    /*
    
    for(ii=0; ii<mesh->Nelem; ii++)
    {
        printf("%i, %i, %i\n", mesh->elem[ii][0], mesh->elem[ii][1], mesh->elem[ii][2]);
    }

    */
    
    for(ii=0; ii<mesh->Np; ii++)
    {
        printf("%.10e, %.10e\n", mesh->p[ii][0], mesh->p[ii][1]);
    }    
    
    fclose(ff);

    return 0;

}

