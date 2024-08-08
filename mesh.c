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
                
        if(c == ' ' || c == '\r')
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

ELEMENT* meshElementMalloc(int type, int Nvar)
{

    ELEMENT* e = malloc(sizeof(ELEMENT));
    
    if(type==5)
    {
        e->Np = 3;
        e->p = malloc(e->Np*sizeof(int));
        e->neiL = malloc(e->Np*sizeof(ELEMENT));
        e->f = malloc(e->Np*sizeof(int));
        e->neiN = 0;        
    }
    else if(type==9)
    {
        e->Np = 4;
        e->p = malloc(e->Np*sizeof(int));
        e->neiL = malloc(e->Np*sizeof(ELEMENT));        
        e->f = malloc(e->Np*sizeof(int));
        e->neiN = 0;        
    }
    else if(type==0)
    {
        e->Np = 2;
        e->p = malloc(e->Np*sizeof(int));
        e->neiL = malloc(sizeof(ELEMENT));                
        e->neiN = 0;        
    }
  
    e->P = malloc((Nvar+1)*sizeof(double));
  
    return e;

}

MESHBC* meshBCread(FILE* ff, int Nvar) 
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

    bc->elemL = malloc(bc->Nelem*sizeof(ELEMENT));
    
    jj = 0;
    while(jj<bc->Nelem)
    {
        bc->elemL[jj] = meshElementMalloc(0, Nvar);   
        bc->elemL[jj]->ii = jj;
    
        meshGetWord(ff, s);
        meshGetWord(ff, s);
        bc->elemL[jj]->p[0] = atoi(s);
        meshGetWord(ff, s);        
        bc->elemL[jj]->p[1] = atoi(s);        
        meshGetWord(ff, s);
        jj++;
    }

    return bc;

}

MESH* meshInit(char* fileName, int Nvar, int axi)
{

    MESH* mesh = malloc(sizeof(MESH));
    FILE* ff = fopen(fileName, "r");
    char s[100];
    int ii;
    int type;
    mesh->axi = axi;
    
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

    mesh->elemL = malloc(mesh->Nelem*sizeof(ELEMENT*));
    
    ii = 0;
    while(ii<mesh->Nelem)
    {
        meshGetWord(ff, s);
        type = atoi(s);
        
        if(type==5)
        {       
            mesh->elemL[ii] = meshElementMalloc(5, Nvar);
         
            meshGetWord(ff, s);
            mesh->elemL[ii]->p[0] = atoi(s);
            meshGetWord(ff, s);        
            mesh->elemL[ii]->p[1] = atoi(s);            
            meshGetWord(ff, s);        
            mesh->elemL[ii]->p[2] = atoi(s);            
        }
        else if(type==9)
        {       
            mesh->elemL[ii] = meshElementMalloc(9, Nvar);
         
            meshGetWord(ff, s);
            mesh->elemL[ii]->p[0] = atoi(s);
            meshGetWord(ff, s);        
            mesh->elemL[ii]->p[1] = atoi(s);            
            meshGetWord(ff, s);                
            mesh->elemL[ii]->p[2] = atoi(s);            
            meshGetWord(ff, s);                
            mesh->elemL[ii]->p[3] = atoi(s);                        
        }
        
        mesh->elemL[ii]->ii = ii;
        
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
        mesh->bc[ii] = meshBCread(ff, Nvar);
        ii++;
    }
    
    fclose(ff);

    printf("mesh: calculating conections.\n");   
    meshCalcConnection2(mesh);
    
    printf("mesh: calculating neighbors.\n");
    meshCalcNeighbors(mesh);
    for(ii=0; ii<mesh->Nmark; ii++)
    {
        meshBCneighbors(mesh->bc[ii], mesh);
    }

    printf("mesh: fix border orientation.\n");
    meshFixBorderOrientation(mesh);

    meshCalcFaces(mesh);   
    
    meshUpdateOmega(mesh);

    return mesh;

}

void meshPrintBC(MESHBC* bc)
{

    printf("%s\n", bc->name);
    printf("%i\n", bc->Nelem);

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        printf("%i, %i\n", bc->elemL[ii]->p[0], bc->elemL[ii]->p[1]);
    }

}

void meshPrint(MESH* mesh)
{

    int ii, jj;

    printf("%i\n", mesh->Nelem);
    for(ii=0; ii<mesh->Nelem; ii++)
    {
        for(jj=0; jj<mesh->elemL[ii]->Np; jj++)
        {
            printf(" %i,", mesh->elemL[ii]->p[jj]);
        }
        printf("\n");
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

void meshElementFree(ELEMENT* e)
{
    if(e->Np==2)
    {
        free(e->p);
        free(e->neiL);        
    }        
    else
    {
        free(e->p);
        free(e->neiL);
        free(e->f);    
    }
    
    free(e->P);
    free(e);
}

void meshBCFree(MESHBC* bc)
{   
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        meshElementFree(bc->elemL[ii]);
    }
    free(bc->elemL);
}

void meshFree(MESH* mesh)
{

    tableFreeDouble(mesh->p, mesh->Np);

    tableFreeInit(mesh->con, mesh->Ncon);
        
    for(int ii=0; ii<mesh->Nelem; ii++)    
    {    
        meshElementFree(mesh->elemL[ii]);
    }
    
    free(mesh->elemL);
    
        
    for(int ii=0; ii<mesh->Nmark; ii++)
    {
        meshBCFree(mesh->bc[ii]);
    }
    
    free(mesh->bc);
    free(mesh);

}

void elementCenter(ELEMENT* E, MESH* mesh, double* x, double* y)    
{

    *x = 0.0;
    *y = 0.0;
    
    for(int jj=0; jj<E->Np; jj++)
    {
        *x += mesh->p[E->p[jj]][0];
        *y += mesh->p[E->p[jj]][1];
    }

    *x /= E->Np;
    *y /= E->Np;

}

double meshCalcOmegaTri(MESH* mesh, int p0, int p1, int p2)
{
 
    double x1 = mesh->p[p0][0];
    double x2 = mesh->p[p1][0];
    double x3 = mesh->p[p2][0];

    double y1 = mesh->p[p0][1];
    double y2 = mesh->p[p1][1];
    double y3 = mesh->p[p2][1];

    return 0.5*((x1 - x2)*(y1 + y2) + (x2 - x3)*(y2 + y3) + (x3 - x1)*(y3 + y1));

}


double meshCalcDSlateral(MESH* mesh, int ii)
{
 
    double ans;
 
    ans = meshCalcOmegaTri(mesh, mesh->elemL[ii]->p[0], mesh->elemL[ii]->p[1], mesh->elemL[ii]->p[2]);
 
    if(mesh->elemL[ii]->Np==4)
    {
        ans += meshCalcOmegaTri(mesh, mesh->elemL[ii]->p[2], mesh->elemL[ii]->p[3], mesh->elemL[ii]->p[0]);
    }

    return ans;

}

double meshCalcOmega(MESH* mesh, int ii)
{
 
    double x, y;
    ELEMENT* E = mesh->elemL[ii];
    
    double ans = meshCalcOmegaTri(mesh, E->p[0], E->p[1], E->p[2]);
 
    if(mesh->elemL[ii]->Np==4)
    {
        ans += meshCalcOmegaTri(mesh, E->p[2], E->p[3], E->p[0]);
    }

    if(mesh->axi == 1)
    {
        elementCenter(E, mesh, &x, &y);
        ans *= y;
    }

    return ans;

}

void meshUpdateOmega(MESH* mesh)
{
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        mesh->elemL[ii]->omega = meshCalcOmega(mesh, ii);
    }
}

double elementIsConnected(ELEMENT* e0, ELEMENT* e1, int* p0, int* p1)
{

    int kk, mm;
    int link = 0;
    int ans = 0;
    int aux;
    
    for(kk=0; kk<e0->Np; kk++)
    {
        for(mm=0; mm<e1->Np; mm++)
        {
            if(e0->p[kk] == e1->p[mm])
            {
                link += 1;

                if(link==1)
                {
                    *p0 = e0->p[kk];
                }
                else if(link==2)
                {
                    *p1 = e0->p[kk];
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
    if((*p0 == e0->p[0]) & (*p1 == e0->p[e0->Np-1]))
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

int meshSameFace(FACETYPE* f0, FACETYPE* f1)
{
    int ans = 0;
    if(!f0->full)
    {    
        if(f0->p0 == f1->p1)
        {
            if(f0->p1 == f1->p0)
            {
                ans = 1;
            }    
        }
        /*
        else if(f0->p0 == f1->p0)
        {
            if(f0->p1 == f1->p1)
            {
                ans = 1;
            }    
        }
        */
    }

    return ans;    
}



void meshCalcConnection1(MESH* mesh)
{

    int p0, p1;

    mesh->Ncon = 0;
    
    CONNECTION* initCon = malloc(sizeof(CONNECTION));
    CONNECTION* con = initCon;

    for(int ii=0; ii<mesh->Nelem-1; ii++)
    {
        for(int jj=ii+1; jj<mesh->Nelem; jj++)
        {
            if(elementIsConnected(mesh->elemL[ii], mesh->elemL[jj], &p0, &p1))
            {
                mesh->Ncon += 1;
                con->data[0] = ii;
                con->data[1] = jj;
                con->data[2] = p0;
                con->data[3] = p1;
                con->next = malloc(sizeof(CONNECTION));
                con = con->next;
            }
        }
    }

    mesh->con = tableMallocInt(mesh->Ncon, 4);
    
    CONNECTION* next;
    
    con = initCon;
    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        mesh->con[ii][0] = con->data[0];
        mesh->con[ii][1] = con->data[1];
        mesh->con[ii][2] = con->data[2];
        mesh->con[ii][3] = con->data[3];
        next = con->next;
        free(con);
        con = next;        
    } 
}

void meshCalcConnection2(MESH* mesh)
{
    mesh->Ncon = 0;
    int Nface = 0;  
    int flagCommon = 0;  
    FACETYPE* initFace = malloc(sizeof(FACETYPE));
    FACETYPE* face = initFace;
    FACETYPE* face0;

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        for(int jj=0; jj<mesh->elemL[ii]->Np; jj++)
        {
            if(ii == 0 && jj == 0)
            {
                face->p0 = mesh->elemL[ii]->p[jj];
                face->p1 = mesh->elemL[ii]->p[(jj+1)%mesh->elemL[ii]->Np];
                face->e0 = ii;
                face->full = 0;
                
                Nface += 1;                
                face->next = malloc(sizeof(FACETYPE));
                face = face->next;    
            }
            else
            {
                face->p0 = mesh->elemL[ii]->p[jj];
                face->p1 = mesh->elemL[ii]->p[(jj+1)%mesh->elemL[ii]->Np];
                face->e0 = ii;
                face->full = 0;
                
                face0 = initFace;
                flagCommon = 0;
                for(int kk=0; kk<Nface; kk++)
                {
                    if(meshSameFace(face0, face))
                    {
                        face0->e1 = face->e0;
                        face0->full = 1;
                        flagCommon = 1;
                        mesh->Ncon += 1;
                        break;
                    }
                    face0 = face0->next;   
                }
            
                if(!flagCommon)
                {
                    Nface += 1;
                    face->next = malloc(sizeof(FACETYPE));
                    face = face->next;                    
                }                            
            }
        }
    }

    mesh->con = tableMallocInt(mesh->Ncon, 4);
    
    face = initFace;
    int jj = 0;
    for(int ii=0; ii<Nface; ii++)
    {
        if(face->full)
        {
            mesh->con[jj][0] = face->e0;
            mesh->con[jj][1] = face->e1;
            mesh->con[jj][2] = face->p0;
            mesh->con[jj][3] = face->p1;
            jj += 1;        
        }
        face0 = face->next;
        free(face);
        face = face0;
    }
    free(face);
}

void meshPrintConnection(MESH* mesh, int N)
{
    for(int ii=0; ii<N; ii++)
    {
        printf("%i, %i, %i, %i,\n", mesh->con[ii][0], mesh->con[ii][1], mesh->con[ii][2], mesh->con[ii][3]);
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

void meshCalcDS2(MESH* mesh, int p0, int p1, double* nx, double* ny, double* dS)
{

    double x0 = mesh->p[p0][0];
    double y0 = mesh->p[p0][1];
    
    double x1 = mesh->p[p1][0];
    double y1 = mesh->p[p1][1];

    double dSx = y1 - y0;
    double dSy = -(x1 - x0);
    
    *dS = sqrt(dSx*dSx + dSy*dSy);
    
    *nx = dSx/(*dS);
    *ny = dSy/(*dS);
    
    if(mesh->axi == 1)
    {
        *dS = (*dS)*(y0 + y1)*0.5;
    }

}

int meshBCIsConnect(ELEMENT* BCe, ELEMENT* e)
{

    int link = 0;
    int ans = 0;
    
    for(int ii=0; ii<2; ii++)
    {
        for(int jj=0; jj<e->Np; jj++)
        {
            if(BCe->p[ii]==e->p[jj])
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

void meshCalcNeighbors(MESH* mesh)
{

    int e0, e1;
    
    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        e0 = mesh->con[ii][0];
        e1 = mesh->con[ii][1];

        mesh->elemL[e0]->neiL[mesh->elemL[e0]->neiN] = mesh->elemL[e1];
        mesh->elemL[e1]->neiL[mesh->elemL[e1]->neiN] = mesh->elemL[e0];
               
        mesh->elemL[e0]->neiN += 1;
        mesh->elemL[e1]->neiN += 1;
        
    }
}    


void meshBCneighbors(MESHBC* bc, MESH* mesh)
{
    
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        for(int jj=0; jj<mesh->Nelem; jj++)
        {
            if(meshBCIsConnect(bc->elemL[ii], mesh->elemL[jj]))
            {   
                bc->elemL[ii]->neiL[0] = mesh->elemL[jj];
                mesh->elemL[jj]->neiL[mesh->elemL[jj]->neiN] = bc->elemL[ii];
                
                bc->elemL[ii]->neiN += 1;
                mesh->elemL[jj]->neiN += 1;
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
    double ans = 0.0;

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
    
    return ans;   
}

void meshPrintDStotal(MESH* mesh)
{

    double dSx, dSy;
    double dSxt, dSyt;
    ELEMENT* e;

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        dSxt = 0.0;
        dSyt = 0.0;

        e = mesh->elemL[ii];

        for(int jj=0; jj<e->Np; jj++)
        {

            meshCalcDS(mesh, e->p[jj], e->p[(jj+1)%e->Np], &dSx, &dSy);
            dSxt += dSx;
            dSyt += dSy;
            
        }
       
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
            e0 = mesh->bc[ii]->elemL[jj]->neiL[0]->ii;
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

int meshPOri(MESH* mesh, ELEMENT* e, int p0, int p1)
{

    int ans = -1;

    for(int ii=0; ii<e->Np; ii++)
    {
        if((e->p[ii] == p0) & (e->p[(ii+1)%e->Np] == p1))
        {
            ans = 1;
        }
    }

    return ans;

}

void meshCheckBorderOrientation(MESH* mesh)
{
    int p0, p1;
    ELEMENT* e0;
    MESHBC* bc;
    
    for(int jj=0; jj<mesh->Nmark; jj++)
    {
        bc = mesh->bc[jj];
        for(int ii=0; ii<bc->Nelem; ii++)
        {

            e0 = bc->elemL[ii]->neiL[0];
            p0 = bc->elemL[ii]->p[0];
            p1 = bc->elemL[ii]->p[1];
            
            printf("%i,\n", meshPOri(mesh, e0, p0, p1));

        }
    }
}

void meshFixBorderOrientation(MESH* mesh)
{
    int p0, p1;
    ELEMENT* e0;
    MESHBC* bc;
    
    for(int jj=0; jj<mesh->Nmark; jj++)
    {
        bc = mesh->bc[jj];
        for(int ii=0; ii<bc->Nelem; ii++)
        {
            e0 = bc->elemL[ii]->neiL[0];
            p0 = bc->elemL[ii]->p[0];
            p1 = bc->elemL[ii]->p[1];
            
            if(meshPOri(mesh, e0, p0, p1)<0)
            {
                bc->elemL[ii]->p[0] = p1;
                bc->elemL[ii]->p[1] = p0;            
            }
        }
    }
}

void meshCalcFaces(MESH* mesh)
{

    int e0, e1;
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        mesh->elemL[ii]->neiN = 0;
        for(int kk=0; kk<mesh->elemL[ii]->Np; kk++)
        {
            mesh->elemL[ii]->f[kk] = 0;
        }
    }
    
    
    for(int ii=0; ii<mesh->Ncon; ii++)
    {
        e0 = mesh->con[ii][0];
        e1 = mesh->con[ii][1];
       
        mesh->elemL[e0]->f[mesh->elemL[e0]->neiN] = ii+1;
        mesh->elemL[e1]->f[mesh->elemL[e1]->neiN] = -(ii+1);
        
        mesh->elemL[e0]->neiN += 1;
        mesh->elemL[e1]->neiN += 1;
        
    }
}

