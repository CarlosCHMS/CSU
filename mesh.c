#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
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
    int ii;
    
    printf("mesh: reading elements.\n");
    
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

    printf("mesh: reading points.\n");

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
    
    printf("mesh: reading marks.\n");
    
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


    printf("mesh: calculating conections.\n");
    meshCalcConnection(mesh);
    
    printf("mesh: calculating domain marks.\n");
    for(ii=0; ii<mesh->Nmark; ii++)
    {
        meshBCDomain(mesh->bc[ii], mesh);
    }

    return mesh;

}

void meshPrintBC(MESHBC* bc)
{

    printf("%s\n", bc->name);
    printf("%i\n", bc->Nelem);

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        printf("%i, %i, %i\n", bc->elem[ii][0], bc->elem[ii][1], bc->domain[ii]);
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

    /*
    printf("%i\n", mesh->Ncon);
    for(ii=0; ii<mesh->Ncon; ii++)
    {
        printf("%i, %i, %i, %i\n", mesh->con[ii][0], mesh->con[ii][1], mesh->con[ii][2], mesh->con[ii][3]);
    }
    */


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
    tableFreeInit(mesh->con, mesh->Ncon);
    
    if(mesh->order == 2)
    {
        tableFreeInit(mesh->nei, mesh->Nelem);
        tableFreeInit(mesh->neip0, mesh->Nelem);
        tableFreeInit(mesh->neip1, mesh->Nelem);        
        free(mesh->neiN);
    }
    
    for(int ii; ii<mesh->Nmark; ii++)
    {

        meshBCFree(mesh->bc[ii]);
        //free(mesh->bc[ii]);

    }

    //free(mesh);

}

void meshElemCenter(MESH* mesh, int ii, double* x, double* y)    
{

    *x = (mesh->p[mesh->elem[ii][0]][0] + mesh->p[mesh->elem[ii][1]][0] + mesh->p[mesh->elem[ii][2]][0])/3;
    *y = (mesh->p[mesh->elem[ii][0]][1] + mesh->p[mesh->elem[ii][1]][1] + mesh->p[mesh->elem[ii][2]][1])/3;

}

double meshCalcDSlateral(MESH* mesh, int ii)
{
 
    double x1 = mesh->p[mesh->elem[ii][0]][0];
    double x2 = mesh->p[mesh->elem[ii][1]][0];
    double x3 = mesh->p[mesh->elem[ii][2]][0];

    double y1 = mesh->p[mesh->elem[ii][0]][1];
    double y2 = mesh->p[mesh->elem[ii][1]][1];
    double y3 = mesh->p[mesh->elem[ii][2]][1];

    return 0.5*((x1 - x2)*(y1 + y2) + (x2 - x3)*(y2 + y3) + (x3 - x1)*(y3 + y1));

}

double meshCalcOmega(MESH* mesh, int ii)
{
 
    double x1 = mesh->p[mesh->elem[ii][0]][0];
    double x2 = mesh->p[mesh->elem[ii][1]][0];
    double x3 = mesh->p[mesh->elem[ii][2]][0];

    double y1 = mesh->p[mesh->elem[ii][0]][1];
    double y2 = mesh->p[mesh->elem[ii][1]][1];
    double y3 = mesh->p[mesh->elem[ii][2]][1];

    double ans = 0.5*((x1 - x2)*(y1 + y2) + (x2 - x3)*(y2 + y3) + (x3 - x1)*(y3 + y1));

    if(mesh->axi == 1)
    {
        ans = ans*(y1 + y2 + y3)/3;
    }

    return ans;

}

double meshIsConnected(MESH* mesh, int ii, int jj, int* p0, int* p1)
{

    int kk, mm;
    int link = 0;
    int ans = 0;
    int aux;
    
    for(kk=0; kk<3; kk++)
    {
        for(mm=0; mm<3; mm++)
        {
            if(mesh->elem[ii][kk] == mesh->elem[jj][mm])
            {
                link += 1;

                if(link==1)
                {
                    *p0 = mesh->elem[ii][kk];
                }
                else if(link==2)
                {
                    *p1 = mesh->elem[ii][kk];
                }
            }
            if(link==2)
            {
                break;
            }
        }
        if(link==2)
        {
            break;
        }
    }

    //This if ensurres correct orientation relatively to the first element
    if((*p0 == mesh->elem[ii][0]) & (*p1 == mesh->elem[ii][2]))
    {
        aux = *p0;
        *p0 = *p1;
        *p1 = aux;
    }

    if(link == 2)
    {
        ans = 1;
    }

    return ans;

}

void meshCalcConnection(MESH* mesh)
{

    int p0, p1, kk, mm;

    mesh->Ncon = 0;

    for(int ii=0; ii<mesh->Nelem-1; ii++)
    {
        mm = 0;
        for(int jj=ii+1; jj<mesh->Nelem; jj++)
        {
            if(meshIsConnected(mesh, ii, jj, &p0, &p1))
            {
                mesh->Ncon += 1;
                mm++;
            }

            if(mm==3)
            {
                break;
            }
        }
    }

    mesh->con = tableMallocInt(mesh->Ncon, 4);

    kk = 0;
    for(int ii=0; ii<mesh->Nelem-1; ii++)
    {
        mm = 0;
        for(int jj=ii+1; jj<mesh->Nelem; jj++)
        {
            if(meshIsConnected(mesh, ii, jj, &p0, &p1))
            {
                mesh->con[kk][0] = ii;
                mesh->con[kk][1] = jj;
                mesh->con[kk][2] = p0;
                mesh->con[kk][3] = p1;
                kk++;
                mm++;
            }

            if(mm==3)
            {
                break;
            }
        }
    }
}

void meshCalcDS(MESH* mesh, int p0, int p1, double* dSx, double* dSy)
{

    double x0 = mesh->p[p0][0];
    double y0 = mesh->p[p0][1];
    
    double x1 = mesh->p[p1][0];
    double y1 = mesh->p[p1][1];

    *dSx = y1 - y0;
    *dSy = -(x1 - x0);
    
    double y = (y0 + y1)/2;
    
    if(mesh->axi == 1)
    {
        *dSx = (*dSx)*y;
        *dSy = (*dSy)*y;        
    }

}

int meshBCIsConnect(int* BCp, int* p)
{

    int link = 0;
    int ans = 0;
    
    for(int ii=0; ii<2; ii++)
    {
        for(int jj=0; jj<3; jj++)
        {
            if(BCp[ii]==p[jj])
            {
                link += 1;
            }
        }
    }
    
    if(link==2)
    {
        ans = 1;
    }
    
    
    return ans;
}

void meshBCDomain(MESHBC* bc, MESH* mesh)
{
    bc->domain = malloc(bc->Nelem*sizeof(int));
    
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        for(int jj=0; jj<mesh->Nelem; jj++)
        {
            if(meshBCIsConnect(bc->elem[ii], mesh->elem[jj]))
            {
                bc->domain[ii] = jj;
                break;
            }
        }
    }   
}

double meshEdgeLength(MESH* mesh, int p0, int p1)
{

    double x0 = mesh->p[p0][0];
    double y0 = mesh->p[p0][1];
    
    double x1 = mesh->p[p1][0];
    double y1 = mesh->p[p1][1];

    return sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));

}

double meshMinEdge(MESH* mesh)
{
    int p0, p1;
    double ans;

    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        p0 = mesh->con[ii][2];
        p1 = mesh->con[ii][3];

        if(ii==0)
        {
            ans = meshEdgeLength(mesh, p0, p1);
        }
        else
        {
            if(ans > meshEdgeLength(mesh, p0, p1))
            {
                ans = meshEdgeLength(mesh, p0, p1);
            }
        }        
    }    
}

void meshPrintDStotal(MESH* mesh)
{

    double dSx, dSy;
    double dSxt, dSyt;

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        dSxt = 0.0;
        dSyt = 0.0;

        meshCalcDS(mesh, mesh->elem[ii][0], mesh->elem[ii][1], &dSx, &dSy);
        dSxt += dSx;
        dSyt += dSy;

        meshCalcDS(mesh, mesh->elem[ii][1], mesh->elem[ii][2], &dSx, &dSy);
        dSxt += dSx;
        dSyt += dSy;

        meshCalcDS(mesh, mesh->elem[ii][2], mesh->elem[ii][0], &dSx, &dSy);
        dSxt += dSx;
        dSyt += dSy;
       
        if(mesh->axi==1)
        {
        
            dSyt -= meshCalcDSlateral(mesh, ii);
        
        }
       
        printf(" %.4e, %.4e,\n", dSxt, dSyt);
        
    }
}

void meshCheckUse(MESH* mesh)
{
    
    int e0, e1;
    int* use = malloc(mesh->Nelem*sizeof(int));

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        use[ii] = 0;
    }

    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        e0 = mesh->con[ii][0];
        e1 = mesh->con[ii][1];
    
        use[e0] += 1;
        use[e1] += 1;
    }

    for(int ii=0; ii<mesh->Nmark; ii++)
    {
        for(int jj=0; jj<mesh->bc[ii]->Nelem; jj++)
        {
            e0 = mesh->bc[ii]->domain[jj];
            use[e0] += 1;
        }
    }

    printf("Use check: cases different of 3,\n");
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        if(use[ii] != 3)
        {
            printf("%i,\n", use[ii]);
        }
    }
    printf("End of use check\n");
    free(use);

}

int meshPOri(MESH* mesh, int e, int p0, int p1)
{

    int ans = -1;

    if(mesh->elem[e][0] == p0 & mesh->elem[e][1] == p1)
    {
        ans = 1;
    }

    if(mesh->elem[e][1] == p0 & mesh->elem[e][2] == p1)
    {
        ans = 1;
    }
    
    if(mesh->elem[e][2] == p0 & mesh->elem[e][0] == p1)
    {
        ans = 1;
    }

    return ans;

}

void meshCheckBorderOrientation(MESHBC* bc, MESH* mesh)
{
    int e0, p0, p1;
    for(int ii=0; ii<bc->Nelem; ii++)
    {

        e0 = bc->domain[ii];
        p0 = bc->elem[ii][0];
        p1 = bc->elem[ii][1];
        
        printf("%i,\n", meshPOri(mesh, e0, p0, p1));

    }
}

void meshCalcNeighbors(MESH* mesh)
{

    int e0, e1, p0, p1;

    mesh->neiN = malloc(mesh->Nelem*sizeof(int));
    mesh->nei = tableMallocInt(mesh->Nelem, 3);
    mesh->neip0 = tableMallocInt(mesh->Nelem, 3);
    mesh->neip1 = tableMallocInt(mesh->Nelem, 3);    
    
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        mesh->neiN[ii] = 0;
    }
    
    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        e0 = mesh->con[ii][0];
        e1 = mesh->con[ii][1];
        p0 = mesh->con[ii][2];
        p1 = mesh->con[ii][3];
        
        mesh->nei[e0][mesh->neiN[e0]] = e1;
        mesh->nei[e1][mesh->neiN[e1]] = e0;

        mesh->neip0[e0][mesh->neiN[e0]] = p0;
        mesh->neip1[e0][mesh->neiN[e0]] = p1;
        mesh->neip0[e1][mesh->neiN[e1]] = p1;
        mesh->neip1[e1][mesh->neiN[e1]] = p0;
          
        mesh->neiN[e0] += 1;
        mesh->neiN[e1] += 1;
    }
}    

void meshCheckNei(MESH* mesh)
{

    double dSx, dSy;
    double dSxt, dSyt;

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        dSxt = 0.0;
        dSyt = 0.0;

        meshCalcDS(mesh, mesh->neip0[ii][0], mesh->neip1[ii][0], &dSx, &dSy);
        dSxt += dSx;
        dSyt += dSy;

        meshCalcDS(mesh, mesh->neip0[ii][1], mesh->neip1[ii][1], &dSx, &dSy);
        dSxt += dSx;
        dSyt += dSy;

        meshCalcDS(mesh, mesh->neip0[ii][2], mesh->neip1[ii][2], &dSx, &dSy);
        dSxt += dSx;
        dSyt += dSy;
       
        if(mesh->axi==1)
        {
        
            dSyt -= meshCalcDSlateral(mesh, ii);
        
        }
       
        if(mesh->neiN[ii]<3)
        {
            printf("%i, %.4e, %.4e,\n", mesh->neiN[ii], dSxt, dSyt);
        }
    }
}
