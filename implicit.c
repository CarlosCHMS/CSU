#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include <stdbool.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"boundary.h"
#include"sa.h"
#include"sst.h"
#include"implicit.h"
#include"gasprop.h"


IMPLICIT* implicitInit(INPUT* input, SOLVER* solver)
{
    IMPLICIT* implicit = malloc(sizeof(IMPLICIT));
    
    implicit->timeScheme = solver->timeScheme;
    
    if(inputNameIsInput(input, "wImp"))
    {
        implicit->wImp = strtod(inputGetValue(input, "wImp"), NULL);
    }
    else
    {
        implicit->wImp = 1.0;
    }
    
    return implicit;
}


void implicitMalloc(IMPLICIT* implicit, SOLVER* solver)
{
    if(implicit->timeScheme == 1 || implicit->timeScheme == 2)
    {
        implicit->dW0 = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
        implicit->dW1 = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
        implicit->D = malloc(solver->mesh->Nelem*sizeof(double));
        implicit->dtL = malloc(solver->mesh->Nelem*sizeof(double));
    }
    
    if(implicit->timeScheme == 2)
    {    
        implicit->BB = malloc(solver->mesh->Nelem*sizeof(BLOCK*));
        implicitInitDPLUR(solver);
    }
}


void implicitFree(IMPLICIT* implicit, SOLVER* solver)
{
    if(implicit->timeScheme == 1 || implicit->timeScheme == 2)
    {
        free(implicit->dW0);
        free(implicit->dW1);
        free(implicit->D);
        free(implicit->dtL);
    }
    
    if(implicit->timeScheme == 2)
    {    
        implicitFreeDPLUR(implicit, solver);
        free(implicit->BB);
    }
    
    free(implicit);
}


void implicitCalcD(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        double Lc = 0;
        double Lv = 0;
        ELEMENT* E = mesh->elemL[ii];

        double mi;
        
        for(int jj=0; jj<E->neiN; jj++)
        {
            int face = E->f[jj];
            int face1 = 0;
            if(face > 0)
            {
                face1 = face-1;
            }
            else if(face < 0)
            {
                face1 = -face-1;            
            }
            
            double dSx, dSy, dS;
            double nx, ny;

            int e0 = mesh->con[face1][0];
            int e1 = mesh->con[face1][1];
            int p0 = mesh->con[face1][2];
            int p1 = mesh->con[face1][3];

            ELEMENT* E0 = mesh->elemL[e0];
            ELEMENT* E1 = mesh->elemL[e1];

            meshCalcDS(mesh, p0, p1, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
                        
            double r = (E0->P[0] + E1->P[0])*0.5;
            double u = (E0->P[1] + E1->P[1])*0.5;
            double v = (E0->P[2] + E1->P[2])*0.5;
            double T = (E0->P[4] + E1->P[4])*0.5;

            double c = gasprop_T2c(solver->gas, T);

            nx = dSx/dS;
            ny = dSy/dS;

            Lc += (fabs(nx*u + ny*v) + c)*dS;
            
            if(solver->laminar)
            {            
                mi = gaspropSutherland(T);
                Lv += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr)*dS*dS;
            }
            
            if(solver->sa1->active || solver->sst->active)
            {
                mi = gaspropSutherland(T);
                Lv += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr + solver->miT[face1]/solver->gas->Pr_t)*dS*dS;
            }
            
        }
        
        implicit->dtL[ii] = Lc;
        implicit->D[ii] = 0.5*implicit->wImp*Lc;        
        if(solver->laminar || solver->sa1->active || solver->sst->active)
        {
            implicit->D[ii] += Lv/mesh->omega[ii];
        }
    }
    
    //Complement from the boundaries
    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {

        double Lc;
        double Lv = 0;
        double mi;
        double dSx, dSy, dS;
        double nx, ny;
        int p0, p1, e1;
        
        MESHBC* bc = solver->mesh->bc[jj];
        for(int ii=0; ii<bc->Nelem; ii++)
        {
 
            e1 = bc->elemL[ii]->neiL[0]->ii;
            p0 = bc->elemL[ii]->p[0];
            p1 = bc->elemL[ii]->p[1];
 
            meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);

            ELEMENT* E0 = bc->elemL[ii];
            //ELEMENT* E1 = mesh->elemL[e1];

            double r = E0->P[0];
            double u = E0->P[1];
            double v = E0->P[2];
            double T = E0->P[4];

            double c = gasprop_T2c(solver->gas, T);
            
            if(dS > 0)
            {
                nx = dSx/dS;
                ny = dSy/dS;
            }
            else
            {
                nx = 0.0;
                ny = 1.0;
            }

            Lc = (fabs(nx*u + ny*v) + c)*dS;
            
            if(solver->laminar)
            {            
                mi = gaspropSutherland(T);
                Lv = fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr)*dS*dS;
            }
            
            if(solver->sa1->active || solver->sst->active)
            {
                mi = gaspropSutherland(T);
                Lv = fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr)*dS*dS;
            }

            implicit->dtL[e1] += Lc;
            implicit->D[e1] += 0.5*implicit->wImp*Lc;        
            if(solver->laminar || solver->sa1->active || solver->sst->active)
            {
                implicit->D[e1] += Lv/mesh->omega[e1];
            }
        }
    }
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        implicit->dtL[ii] = mesh->omega[ii]/implicit->dtL[ii];
    }

    // Time step calculation
    double dt = implicit->dtL[0];
    for(int ii=1; ii<mesh->Nelem; ii++)
    {
        dt = fmin(dt, implicit->dtL[ii]);
    }

    dt *= solver->CFL;

    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        implicit->D[ii] += mesh->omega[ii]/dt;
    }  
}


void implicitLUSGS_L(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;
    
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            implicit->dW0[kk][ii] = -solver->R[kk][ii];
        }
        
        ELEMENT* E = mesh->elemL[ii];
                
        for(int jj=0; jj<E->neiN; jj++)
        {
            int face = E->f[jj];
            int face1 = 0;
            int e0 = 0, e1 = 0, p0 = 0, p1 = 0;
            if(face > 0)
            {
                face1 = face-1;
                e0 = mesh->con[face1][0];
                e1 = mesh->con[face1][1];
                p0 = mesh->con[face1][2];
                p1 = mesh->con[face1][3];
            }
            else if(face < 0)
            {
                face1 = -face-1;
                e1 = mesh->con[face1][0];
                e0 = mesh->con[face1][1];
                p1 = mesh->con[face1][2];
                p0 = mesh->con[face1][3];
            }
            
            if(e1 < e0)
            {   
                implicitFunc(solver, e0, e1, p0, p1, face1, implicit->dW0);
            }
            
        }
    
        if(solver->sst->active)
        {
            for(int kk=0; kk<4; kk++)
            {
                implicit->dW0[kk][ii] /= implicit->D[ii];            
            }
            
            double A = implicit->dW0[4][ii] + solver->sst->dQkdr[ii]*implicit->dW0[0][ii];
            double B = implicit->dW0[5][ii] + solver->sst->dQodr[ii]*implicit->dW0[0][ii];
            double a = implicit->D[ii] - solver->sst->dQkdrk[ii];
            double b = -solver->sst->dQkdro[ii];
            double c = -solver->sst->dQodrk[ii];
            double d = implicit->D[ii] - solver->sst->dQodro[ii];
            
            long double det = a*d - b*c;
            
            if(fabs(det) < 1e-14)
            {
                implicit->dW0[4][ii] /= implicit->D[ii];
                implicit->dW0[5][ii] /= implicit->D[ii];
            }
            else
            {
                implicit->dW0[4][ii] = (A*d - B*b)/det;
                implicit->dW0[5][ii] = (-A*c + B*a)/det;   
            }
            
        }
        else
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                implicit->dW0[kk][ii] /= implicit->D[ii];
            } 
        }        
    }
    
}

void implicitLUSGS_U(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        if(solver->sst->active)
        {
            for(int kk=0; kk<4; kk++)
            {
                implicit->dW1[kk][ii] = implicit->D[ii]*implicit->dW0[kk][ii];
            }
            double a = implicit->D[ii] - solver->sst->dQkdrk[ii];
            double b = -solver->sst->dQkdro[ii];
            double c = -solver->sst->dQodrk[ii];
            double d = implicit->D[ii] - solver->sst->dQodro[ii];            
            
            long double det = a*d - b*c;
            if(fabs(det) < 1e-14)
            {
                implicit->dW1[4][ii] = implicit->D[ii]*implicit->dW0[4][ii];
                implicit->dW1[5][ii] = implicit->D[ii]*implicit->dW0[5][ii];
            }
            else
            {
                implicit->dW1[4][ii] = a*implicit->dW0[4][ii] + b*implicit->dW0[5][ii] - solver->sst->dQkdr[ii]*implicit->dW0[0][ii];
                implicit->dW1[5][ii] = c*implicit->dW0[4][ii] + d*implicit->dW0[5][ii] - solver->sst->dQodr[ii]*implicit->dW0[0][ii];
            }
        }
        else
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                implicit->dW1[kk][ii] = implicit->D[ii]*implicit->dW0[kk][ii];
            }
        }
    }
    
    for(int ii=mesh->Nelem-1; ii>=0; ii--)
    {   
        ELEMENT* E = mesh->elemL[ii];
        
        for(int jj=0; jj<E->neiN; jj++)
        {
            int face = E->f[jj];
            int face1 = 0;
            int e0 = 0, e1 = 0, p0 = 0, p1 = 0;
            if(face > 0)
            {
                face1 = face-1;
                e0 = mesh->con[face1][0];
                e1 = mesh->con[face1][1];
                p0 = mesh->con[face1][2];
                p1 = mesh->con[face1][3];
            }
            else if(face < 0)
            {
                face1 = -face-1;
                e1 = mesh->con[face1][0];
                e0 = mesh->con[face1][1];
                p1 = mesh->con[face1][2];
                p0 = mesh->con[face1][3];
            }            

            if(e1 > e0)
            {
                implicitFunc(solver, e0, e1, p0, p1, face1, implicit->dW1);               
            }
        }
        
        if(solver->sst->active)
        {
            for(int kk=0; kk<4; kk++)
            {
                implicit->dW1[kk][ii] /= implicit->D[ii];            
            }
            
            double A = implicit->dW1[4][ii] + solver->sst->dQkdr[ii]*implicit->dW1[0][ii];
            double B = implicit->dW1[5][ii] + solver->sst->dQodr[ii]*implicit->dW1[0][ii];
            double a = implicit->D[ii] - solver->sst->dQkdrk[ii];
            double b = -solver->sst->dQkdro[ii];
            double c = -solver->sst->dQodrk[ii];
            double d = implicit->D[ii] - solver->sst->dQodro[ii];
            
            long double det = a*d - b*c;
            if(fabs(det) < 1e-14)
            {
                implicit->dW1[4][ii] /= implicit->D[ii];
                implicit->dW1[5][ii] /= implicit->D[ii];
            }
            else
            {
                implicit->dW1[4][ii] = (A*d - B*b)/det;
                implicit->dW1[5][ii] = (-A*c + B*a)/det;   
            }

        }
        else
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                implicit->dW1[kk][ii] /= implicit->D[ii];
            }        
        }
    }
}


void implicitAuxCalcFlux(SOLVER* solver, double* U, double p, double nx, double ny, double* F)
{
    double V = (nx*U[1] + ny*U[2])/U[0];      
    
    F[0] = U[0]*V;
    F[1] = U[1]*V + nx*p;    
    F[2] = U[2]*V + ny*p;
    F[3] = (U[3] + p)*V;
    
    for(int kk=4; kk<solver->Nvar; kk++)
    {
        F[kk] = U[kk]*V;
    }
    
}

void implicitCalcDeltaFlux(SOLVER* solver, double* P, double* dW, double nx, double ny, double* dF)
{

    double rho = P[0];
    double u = P[1];  
    double v = P[2];  
    double p = P[3];
    
    double k = 0;
    if(solver->sst->active)
    {
        k = P[4];
    }
    
    double F0[6];
    double F1[6];

    double T = p/(solver->gas->R*rho);
    double e0 = gasprop_T2e(solver->gas, T);
	double E = e0 + (u*u + v*v)*0.5 + k;

    double U[6];

    U[0] = rho;
    U[1] = rho*u;
    U[2] = rho*v;
    U[3] = rho*E;
    
    for(int kk=4; kk<solver->Nvar; kk++)
    {
        U[kk] = rho*P[kk];
    }

    implicitAuxCalcFlux(solver, U, p, nx, ny, F0);

    for(int kk=0; kk<solver->Nvar; kk++)
    {
        U[kk] += dW[kk];
    }
    
    rho = U[0];
    rho = fmax(rho, solver->rLim);
    U[0] = rho;
    
    u = U[1]/rho;
    v = U[2]/rho;

    if(solver->sa1->active)
    {
        if(U[4] < 0)
        {
            U[4] = 0;
        }
    }

    k = 0;
    if(solver->sst->active)
    {
        k = U[4]/rho;
        if(k < 1e-14)
        {
            k = 1e-14;
            U[4] = rho*k;
        }

        double om = U[5]/rho;
        if(om < solver->omLim)
        {
            om = solver->omLim;
            U[5] = rho*om;
        }
    }

    E = U[3]/rho;
    double e1 = E - (u*u + v*v)*0.5 - k;
    
    e1 = fmax(e1, 7.176325e+02);

    T = gasprop_e2Taprox(solver->gas, e0, T, e1);
    if(T < 1)
    {
        T = 1;
        e1 = gasprop_T2e(solver->gas, T);
	    E = e1 + (u*u + v*v)*0.5 + k;
        U[3] = rho*E;
    }
    
    p = solver->gas->R*T*rho;
    
    implicitAuxCalcFlux(solver, U, p, nx, ny, F1);

    for(int kk=0; kk<solver->Nvar; kk++)
    {
        dF[kk] = F1[kk] - F0[kk];
    }
}


void implicitFunc(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW)
{
    IMPLICIT* implicit = solver->implicit;

    double dSx, dSy, dS;
    double nx, ny;
    double x0, y0, x1, y1;
    
    ELEMENT* E0 = solver->mesh->elemL[e0];
    ELEMENT* E1 = solver->mesh->elemL[e1];        

    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    dS = sqrt(dSx*dSx + dSy*dSy);
        
    nx = dSx/dS;
    ny = dSy/dS;
    
    double dF[6];
    
    double P1[6];
    double dW1[6];
    
    P1[0] = E1->P[0];
    P1[1] = E1->P[1];    
    P1[2] = E1->P[2];
    P1[3] = E1->P[3];    
	for(int kk=4; kk<solver->Nvar; kk++)
    {
	    P1[kk] = E1->P[kk+1];
    }    
    
    for(int kk=0; kk<solver->Nvar; kk++)
    {
	    dW1[kk] = dW[kk][e1];
    } 
    
    implicitCalcDeltaFlux(solver, P1, dW1, nx, ny, dF);

    double T = E1->P[4];
    
    double c = gasprop_T2c(solver->gas, T);
    double ra = implicit->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;
    
    if(solver->laminar)
    {
        double r = E1->P[0];                

        double mi = gaspropSutherland(T);

        elementCenter(E0, solver->mesh, &x0, &y0);
        elementCenter(E1, solver->mesh, &x1, &y1);
        
        double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr)*dS/d;
    }    
    else if(solver->sa1->active || solver->sst->active)
    {
        double r = E1->P[0];                

        double mi = gaspropSutherland(T);

        elementCenter(E0, solver->mesh, &x0, &y0);
        elementCenter(E1, solver->mesh, &x1, &y1);
        
        double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr + solver->miT[face1]/solver->gas->Pr_t)*dS/d;
    }
    
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        dW[kk][e0] -= 0.5*(dF[kk]*dS - ra*dW[kk][e1]);
    }
}


void implicitTest(SOLVER* solver)
{
    IMPLICIT* implicit = solver->implicit;
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            implicit->dW1[kk][ii] = -solver->R[kk][ii]/implicit->D[ii];
        }
    }
}


void implicitInitDPLUR(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;

    #pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        ELEMENT* E = solver->mesh->elemL[ii];
        BLOCK* B;
                
        for(int jj=0; jj<E->neiN; jj++)
        {
            int face = E->f[jj];
            int e1 = 0, face1 = 0;
            if(face > 0)
            {
                face1 = face-1;
                e1 = mesh->con[face1][1];
            }
            else if(face < 0)
            {
                face1 = -face-1;
                e1 = mesh->con[face1][0];
            }

            if(jj == 0)
            {
                B = malloc(sizeof(BLOCK));
                B->A = tableMallocDouble(solver->Nvar, solver->Nvar);
                B->ii = e1;            
                implicit->BB[ii] = B;
            }
            else
            {
                B->next = malloc(sizeof(BLOCK));
                B = B->next;
                B->A = tableMallocDouble(solver->Nvar, solver->Nvar);
                B->ii = e1;            
            }
        }
    }
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            implicit->dW1[kk][ii] = 0.0;
        }
    }    
}


void implicitFreeDPLUR(IMPLICIT* implicit, SOLVER* solver)
{
    #pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        ELEMENT* E = solver->mesh->elemL[ii];
        BLOCK* B = NULL;
        BLOCK* Bold;
        
        for(int jj=0; jj<E->neiN; jj++)
        {
            if(jj == 0)
            {
                B = implicit->BB[ii];
                tableFreeDouble(B->A, solver->Nvar);
            }
            else
            {
                Bold = B;
                B = Bold->next;
                free(Bold);
                tableFreeDouble(B->A, solver->Nvar);
            }
        }
        free(B);
    }
}

void implicitUpdateA(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;
    
    double h = 1.0e-6;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        BLOCK* B;
        ELEMENT* E = mesh->elemL[ii];
                
        for(int jj=0; jj<E->neiN; jj++)
        {        
            if(jj==0)
            {
                B = implicit->BB[ii];
            }
            else
            {
                B = B->next;
            }
            
            int face = E->f[jj];
            int face1 = 0;
            int e0 = 0, e1 = 0, p0 = 0, p1 = 0;
            if(face > 0)
            {
                face1 = face-1;
                e0 = mesh->con[face1][0];
                e1 = mesh->con[face1][1];
                p0 = mesh->con[face1][2];
                p1 = mesh->con[face1][3];
            }
            else if(face < 0)
            {
                face1 = -face-1;
                e1 = mesh->con[face1][0];
                e0 = mesh->con[face1][1];
                p1 = mesh->con[face1][2];
                p0 = mesh->con[face1][3];
            }            
            
            double dSx, dSy, dS;
            double nx, ny;
            double x0, y0, x1, y1;
            
            ELEMENT* E0 = solver->mesh->elemL[e0];
            ELEMENT* E1 = solver->mesh->elemL[e1];        

            meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
                
            nx = dSx/dS;
            ny = dSy/dS;
            
            double F0[6];
            double F[6];
            double U[6];

            double rho, u, v, E, p, T;           
            
            rho = solver->U[0][e1];
            u = solver->U[1][e1]/rho;
            v = solver->U[2][e1]/rho;
            E = solver->U[3][e1]/rho;

            for(int kk=0; kk<solver->Nvar; kk++)
            {
                U[kk] = solver->U[kk][e1];
            } 

            implicitAuxCalcFlux(solver, U, E1->P[3], nx, ny, F0);

            double ee0 = E - 0.5*(u*u + v*v);

            for(int mm=0; mm<4; mm++)
            {            
                for(int nn=0; nn<4; nn++)
                {   
                    if(mm == nn)
                    {
                        U[nn] = solver->U[nn][e1] + h;
                    }
                    else
                    {
                        U[nn] = solver->U[nn][e1];
                    }
                } 
                
                rho = U[0];
                u = U[1]/rho;
                v = U[2]/rho;
                E = U[3]/rho;

                double ee1 = E - 0.5*(u*u + v*v);

                T = gasprop_e2Taprox(solver->gas, ee0, E1->P[4], ee1);
                p = solver->gas->R*rho*T;            

                implicitAuxCalcFlux(solver, U, p, nx, ny, F);
                
                for(int nn=0; nn<4; nn++)
                {
                    B->A[nn][mm] = 0.5*(F[nn] - F0[nn])*dS/h;
                }
            }

            double c = gasprop_T2c(solver->gas, E1->P[4]);
            double ra = implicit->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;

            if(solver->laminar)
            {
                double r = (E0->P[0] + E1->P[0])*0.5;                
                T = (E0->P[4] + E1->P[4])*0.5;
                double mi = gaspropSutherland(T);

                elementCenter(E0, solver->mesh, &x0, &y0);
                elementCenter(E1, solver->mesh, &x1, &y1);
                
                double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr)*dS/d;                    
            }

            for(int nn=0; nn<4; nn++)
            {
                B->A[nn][nn] -= 0.5*ra;
            }              
        }
    }
}

void implicitUpdateA_sa(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;
    double h = 1.0e-6;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        BLOCK* B;
        ELEMENT* E = mesh->elemL[ii];
                
        for(int jj=0; jj<E->neiN; jj++)
        {        
            if(jj==0)
            {
                B = implicit->BB[ii];
            }
            else
            {
                B = B->next;
            }
            
            int face = E->f[jj];
            int face1 = 0;
            int e0 = 0, e1 = 0, p0 = 0, p1 = 0;
            if(face > 0)
            {
                face1 = face-1;
                e0 = mesh->con[face1][0];
                e1 = mesh->con[face1][1];
                p0 = mesh->con[face1][2];
                p1 = mesh->con[face1][3];
            }
            else if(face < 0)
            {
                face1 = -face-1;
                e1 = mesh->con[face1][0];
                e0 = mesh->con[face1][1];
                p1 = mesh->con[face1][2];
                p0 = mesh->con[face1][3];
            }            
            
            double dSx, dSy, dS;
            double nx, ny;
            double x0, y0, x1, y1;
            
            ELEMENT* E0 = solver->mesh->elemL[e0];
            ELEMENT* E1 = solver->mesh->elemL[e1];        

            meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
                
            nx = dSx/dS;
            ny = dSy/dS;
            
            double F0[6];
            double F[6];
            double U[6];

            double rho, u, v, E, p, T;           
            
            rho = solver->U[0][e1];
            u = solver->U[1][e1]/rho;
            v = solver->U[2][e1]/rho;
            E = solver->U[3][e1]/rho;
            
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                U[kk] = solver->U[kk][e1];
            } 

            implicitAuxCalcFlux(solver, U, E1->P[3], nx, ny, F0);

            double ee0 = E - 0.5*(u*u + v*v);

            for(int mm=0; mm<5; mm++)
            {            
                for(int nn=0; nn<5; nn++)
                {   
                    if(mm == nn)
                    {
                        U[nn] = solver->U[nn][e1] + h;
                    }
                    else
                    {
                        U[nn] = solver->U[nn][e1];
                    }
                } 
                
                rho = U[0];
                u = U[1]/rho;
                v = U[2]/rho;
                E = U[3]/rho;

                double ee1 = E - 0.5*(u*u + v*v);

                T = gasprop_e2Taprox(solver->gas, ee0, E1->P[4], ee1);
                p = solver->gas->R*rho*T;            

                implicitAuxCalcFlux(solver, U, p, nx, ny, F);
                
                for(int nn=0; nn<5; nn++)
                {
                    B->A[nn][mm] = 0.5*(F[nn] - F0[nn])*dS/h;
                }
            }

            double c = gasprop_T2c(solver->gas, E1->P[4]);
            double ra = implicit->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;
            
            double r = (E0->P[0] + E1->P[0])*0.5;                
            T = (E0->P[4] + E1->P[4])*0.5;
            double mi = gaspropSutherland(T);

            elementCenter(E0, solver->mesh, &x0, &y0);
            elementCenter(E1, solver->mesh, &x1, &y1);
            
            double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->gas->Pr + solver->miT[face1]/solver->gas->Pr_t)*dS/d;
    
            for(int nn=0; nn<5; nn++)
            {
                B->A[nn][nn] -= 0.5*ra;
            }              
        }
    }
}

void implicitMultA(SOLVER* solver, double** x, double** y)
{
    MESH* mesh = solver->mesh;
    IMPLICIT* implicit = solver->implicit;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        BLOCK* B = NULL;

        for(int kk=0; kk<solver->Nvar; kk++)
        {
            y[kk][ii] = implicit->D[ii]*x[kk][ii];
        }

        for(int jj=0; jj<mesh->elemL[ii]->neiN; jj++)
        {
            if(jj==0)
            {
                B = implicit->BB[ii];
            }
            else
            {
                B = B->next;
            }

            for(int kk=0; kk<solver->Nvar; kk++)
            {
                for(int nn=0; nn<solver->Nvar; nn++)
                {
                    y[kk][ii] += B->A[kk][nn]*x[nn][B->ii];
                }
            }
        }
    }
}


double implicitProdInter(SOLVER* solver, double** x, double** y)
{
    double ans = 0;
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            ans += x[kk][ii]*y[kk][ii];
        }
    }

    return ans;    
}


void implicitCalcDPLUR(SOLVER* solver)
{
    IMPLICIT* implicit = solver->implicit;
    
    implicitCalcD(solver);

    if(solver->sa1->active)
    {
        implicitUpdateA_sa(solver);
    }
    else
    {
        implicitUpdateA(solver);
    }

    double w = 2./3.;
    double** aux;
    int mm;

    for(mm=0; mm<solver->Nlinear; mm++)
    {
        aux = implicit->dW1;
        implicit->dW1 = implicit->dW0;
        implicit->dW0 = aux;    
        
        implicitMultA(solver, implicit->dW0, implicit->dW1);
            
        #pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                implicit->dW1[kk][ii] = w*(-solver->R[kk][ii] - implicit->dW1[kk][ii])/implicit->D[ii] + implicit->dW0[kk][ii];
            }
        }
        
        double ans = fabs(implicit->dW1[0][0] - implicit->dW0[0][0]);
        for(int ii=1; ii<solver->mesh->Nelem; ii++)
        {
            ans = fmax(ans, fabs(implicit->dW1[0][ii] - implicit->dW0[0][ii]));
        }
        
        if(ans<1e-14)
        {
            break;
        }
    }
}

