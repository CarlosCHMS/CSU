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
#include"flux.h"
#include"boundary.h"
#include"sa.h"
#include"implicit.h"
#include"gasprop.h"


void implicitCalcD(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    
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
                mi = sutherland(T);
                Lv += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr)*dS*dS;
            }
            
            if(solver->sa || solver->sstFlag)
            {
                mi = sutherland(T);
                Lv += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr + solver->miT[face1]/solver->Pr_t)*dS*dS;
            }
            
        }
        
        solver->dtL[ii] = Lc;
        solver->D[ii] = 0.5*solver->wImp*Lc;        
        if(solver->laminar || solver->sa || solver->sstFlag)
        {
            solver->D[ii] += Lv/E->omega;
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

            //ELEMENT* E0 = bc->elemL[ii];
            ELEMENT* E1 = mesh->elemL[e1];

            double r = E1->P[0];
            double u = E1->P[1];
            double v = E1->P[2];
            double T = E1->P[4];

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
                mi = sutherland(T);
                Lv = fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr)*dS*dS;
            }
            
            if(solver->sa || solver->sstFlag)
            {
                mi = sutherland(T);
                Lv = fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr)*dS*dS;
            }

            solver->dtL[e1] += Lc;
            solver->D[e1] += 0.5*solver->wImp*Lc;        
            if(solver->laminar || solver->sa || solver->sstFlag)
            {
                solver->D[e1] += Lv/E1->omega;
            }
        }
    }
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        solver->dtL[ii] = mesh->elemL[ii]->omega/solver->dtL[ii];
    }

    // Time step calculation
    double dt = solver->dtL[0];
    for(int ii=1; ii<mesh->Nelem; ii++)
    {
        dt = fmin(dt, solver->dtL[ii]);
    }

    dt *= solver->CFL;

    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        solver->D[ii] += mesh->elemL[ii]->omega/dt;
    }  
}

void implicitAuxCalcFlux2(SOLVER* solver, double U0, double U1, double U2, double U3, double p, double nx, double ny, double* F)
{
    double V = (nx*U1 + ny*U2)/U0;  
    
    F[0] = U0*V;
    F[1] = U1*V + nx*p;    
    F[2] = U2*V + ny*p;
    F[3] = (U3 + p)*V;
}


void implicitCalcDeltaFlux(SOLVER* solver, double rho, double u, double v, double p, double d0, double d1, double d2, double d3, double nx, double ny, double* dF)
{
    double F0[4];
    double F1[4];

    double T = p/(solver->gas->R*rho);
    double e0 = gasprop_T2e(solver->gas, T);
	double E = e0 + (u*u + v*v)/2;

    double U0 = rho;
    double U1 = rho*u;
    double U2 = rho*v;
    double U3 = rho*E;  

    implicitAuxCalcFlux2(solver, U0, U1, U2, U3, p, nx, ny, F0);

    U0 = U0 + d0;
    U1 = U1 + d1;
    U2 = U2 + d2;
    U3 = U3 + d3;
    
    rho = U0;
    u = U1/rho;
    v = U2/rho;
    E = U3/rho;

    double e1 = E - 0.5*(u*u + v*v);

    T = gasprop_e2Taprox(solver->gas, e0, T, e1);
    p = solver->gas->R*T*rho;

    implicitAuxCalcFlux2(solver, U0, U1, U2, U3, p, nx, ny, F1);

    for(int kk=0; kk<4; kk++)
    {
        dF[kk] = F1[kk] - F0[kk];
    }
}

void implicitFunc(SOLVER* solver, int e0, int e1, int p0, int p1, double** dW)
{
    double dSx, dSy, dS;
    double nx, ny;
    double x0, y0, x1, y1;
    
    ELEMENT* E0 = solver->mesh->elemL[e0];
    ELEMENT* E1 = solver->mesh->elemL[e1];        

    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    dS = sqrt(dSx*dSx + dSy*dSy);
        
    nx = dSx/dS;
    ny = dSy/dS;
    
    double dF[4];
    
    implicitCalcDeltaFlux(solver, E1->P[0], E1->P[1], E1->P[2], E1->P[3], dW[0][e1], dW[1][e1], dW[2][e1], dW[3][e1], nx, ny, dF);

	double c = gasprop_T2c(solver->gas, E1->P[4]);
    double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;

    if(solver->laminar)
    {
        double r = E1->P[0];                
        double T = E1->P[4];
        double mi = sutherland(T);

        elementCenter(E0, solver->mesh, &x0, &y0);
        elementCenter(E1, solver->mesh, &x1, &y1);
        
        double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr)*dS/d;                    
    }
    
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        dW[kk][e0] -= 0.5*(dF[kk]*dS - ra*dW[kk][e1]);
    }   
}

void implicitLUSGS_L(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->dW0[kk][ii] = -solver->R[kk][ii];
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
                if(solver->sstFlag)
                {
                    implicitFunc_sst(solver, e0, e1, p0, p1, face1, solver->dW0);                
                } 
                else if(solver->sa)
                {
                    implicitFunc_sa(solver, e0, e1, p0, p1, face1, solver->dW0);                
                }
                else
                {
                    implicitFunc(solver, e0, e1, p0, p1, solver->dW0);
                }
            }
            
        }
    
        if(solver->sstFlag)
        {
            for(int kk=0; kk<4; kk++)
            {
                solver->dW0[kk][ii] /= solver->D[ii];            
            }
            
            double A = solver->dW0[4][ii] + solver->dQkdr[ii]*solver->dW0[0][ii];
            double B = solver->dW0[5][ii] + solver->dQodr[ii]*solver->dW0[0][ii];
            double a = solver->D[ii] - solver->dQkdrk[ii];
            double b = -solver->dQkdro[ii];
            double c = -solver->dQodrk[ii];
            double d = solver->D[ii] - solver->dQodro[ii];
            
            long double det = a*d - b*c;
            
            if(fabs(det) < 1e-14)
            {
                solver->dW0[4][ii] /= solver->D[ii];
                solver->dW0[5][ii] /= solver->D[ii];
            }
            else
            {
                solver->dW0[4][ii] = (A*d - B*b)/det;
                solver->dW0[5][ii] = (-A*c + B*a)/det;   
            }
            
        }
        else
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                solver->dW0[kk][ii] /= solver->D[ii];
            } 
        }        
    }
    
}

void implicitLUSGS_U(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        if(solver->sstFlag)
        {
            for(int kk=0; kk<4; kk++)
            {
                solver->dW1[kk][ii] = solver->D[ii]*solver->dW0[kk][ii];
            }
            double a = solver->D[ii] - solver->dQkdrk[ii];
            double b = -solver->dQkdro[ii];
            double c = -solver->dQodrk[ii];
            double d = solver->D[ii] - solver->dQodro[ii];            
            
            long double det = a*d - b*c;
            if(fabs(det) < 1e-14)
            {
                solver->dW1[4][ii] = solver->D[ii]*solver->dW0[4][ii];
                solver->dW1[5][ii] = solver->D[ii]*solver->dW0[5][ii];
            }
            else
            {
                solver->dW1[4][ii] = a*solver->dW0[4][ii] + b*solver->dW0[5][ii] - solver->dQkdr[ii]*solver->dW0[0][ii];
                solver->dW1[5][ii] = c*solver->dW0[4][ii] + d*solver->dW0[5][ii] - solver->dQodr[ii]*solver->dW0[0][ii];
            }
        }
        else
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                solver->dW1[kk][ii] = solver->D[ii]*solver->dW0[kk][ii];
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
                if(solver->sstFlag)
                {
                    implicitFunc_sst(solver, e0, e1, p0, p1, face1, solver->dW1);                
                } 
                else if(solver->sa)
                {
                    implicitFunc_sa(solver, e0, e1, p0, p1, face1, solver->dW1);                
                }
                else
                {
                    implicitFunc(solver, e0, e1, p0, p1, solver->dW1);
                }                
            }
        }
        
        if(solver->sstFlag)
        {
            for(int kk=0; kk<4; kk++)
            {
                solver->dW1[kk][ii] /= solver->D[ii];            
            }
            
            double A = solver->dW1[4][ii] + solver->dQkdr[ii]*solver->dW1[0][ii];
            double B = solver->dW1[5][ii] + solver->dQodr[ii]*solver->dW1[0][ii];
            double a = solver->D[ii] - solver->dQkdrk[ii];
            double b = -solver->dQkdro[ii];
            double c = -solver->dQodrk[ii];
            double d = solver->D[ii] - solver->dQodro[ii];
            
            long double det = a*d - b*c;
            if(fabs(det) < 1e-14)
            {
                solver->dW1[4][ii] /= solver->D[ii];
                solver->dW1[5][ii] /= solver->D[ii];
            }
            else
            {
                solver->dW1[4][ii] = (A*d - B*b)/det;
                solver->dW1[5][ii] = (-A*c + B*a)/det;   
            }

        }
        else
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                solver->dW1[kk][ii] /= solver->D[ii];
            }        
        }
    }
}

void implicitAuxCalcFlux_sa2(SOLVER* solver, double U0, double U1, double U2, double U3, double U4, double p, double nx, double ny, double* F)
{
    double V = (nx*U1 + ny*U2)/U0;      
    
    F[0] = U0*V;
    F[1] = U1*V + nx*p;    
    F[2] = U2*V + ny*p;
    F[3] = (U3 + p)*V;
    F[4] = U4*V;
}

void implicitCalcDeltaFlux_sa(SOLVER* solver, double rho, double u, double v, double p, double n, double d0, double d1, double d2, double d3, double d4, double nx, double ny, double* dF)
{

    double F0[5];
    double F1[5];

    double T = p/(solver->gas->R*rho);
    double e0 = gasprop_T2e(solver->gas, T);
	double E = e0 + (u*u + v*v)/2;

    double U0 = rho;
    double U1 = rho*u;
    double U2 = rho*v;
    double U3 = rho*E;  
    double U4 = rho*n;

    implicitAuxCalcFlux_sa2(solver, U0, U1, U2, U3, U4, p, nx, ny, F0);

    U0 = U0 + d0;
    U1 = U1 + d1;
    U2 = U2 + d2;
    U3 = U3 + d3;
    U4 = U4 + d4;    

    rho = U0;
    rho = fmax(rho, solver->rLim);
    U0 = rho;
    
    u = U1/rho;
    v = U2/rho;
    
    E = U3/rho;
    double e1 = E - (u*u + v*v)*0.5;
    e1 = fmax(e1, 7.176325e+02);

    T = gasprop_e2Taprox(solver->gas, e0, T, e1);
    if(T < 1)
    {
        T = 1;
        e1 = gasprop_T2e(solver->gas, T);
	    E = e1 + (u*u + v*v)*0.5;
        U3 = rho*E;
    }
    
    p = solver->gas->R*T*rho;

    n = U4/rho;
    n = fmax(n, 0);
    U4 = rho*n;

    implicitAuxCalcFlux_sa2(solver, U0, U1, U2, U3, U4, p, nx, ny, F1);

    for(int kk=0; kk<5; kk++)
    {
        dF[kk] = F1[kk] - F0[kk];
    }

}


void implicitFunc_sa(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW)
{
    double dSx, dSy, dS;
    double nx, ny;
    double x0, y0, x1, y1;
    
    ELEMENT* E0 = solver->mesh->elemL[e0];
    ELEMENT* E1 = solver->mesh->elemL[e1];        

    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    dS = sqrt(dSx*dSx + dSy*dSy);
        
    nx = dSx/dS;
    ny = dSy/dS;
    
    double dF[5];
    
    implicitCalcDeltaFlux_sa(solver, E1->P[0], E1->P[1], E1->P[2], E1->P[3], E1->P[5], dW[0][e1], dW[1][e1], dW[2][e1], dW[3][e1], dW[4][e1], nx, ny, dF);

    double T = E1->P[4];
    
    double c = gasprop_T2c(solver->gas, T);
    double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;
    
    double r = E1->P[0];                

    double mi = sutherland(T);

    elementCenter(E0, solver->mesh, &x0, &y0);
    elementCenter(E1, solver->mesh, &x1, &y1);
    
    double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr + solver->miT[face1]/solver->Pr_t)*dS/d;

    for(int kk=0; kk<solver->Nvar; kk++)
    {
        dW[kk][e0] -= 0.5*(dF[kk]*dS - ra*dW[kk][e1]);
    }   
}

void implicitTest(SOLVER* solver)
{
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->dW1[kk][ii] = -solver->R[kk][ii]/solver->D[ii];
        }
    }
}


void implicitInitDPLUR(SOLVER* solver)
{
    MESH* mesh = solver->mesh;

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
                solver->BB[ii] = B;
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
            solver->dW1[kk][ii] = 0.0;
        }
    }    
}


void implicitFreeDPLUR(SOLVER* solver)
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
                B = solver->BB[ii];
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
                B = solver->BB[ii];
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
            
            double F0[4];
            double F[4];
            double U[4];

            double rho, u, v, E, p, T;           
            
            rho = solver->U[0][e1];
            u = solver->U[1][e1]/rho;
            v = solver->U[2][e1]/rho;
            E = solver->U[3][e1]/rho;

            implicitAuxCalcFlux2(solver, solver->U[0][e1], solver->U[1][e1], solver->U[2][e1], solver->U[3][e1], E1->P[3], nx, ny, F0);

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

                implicitAuxCalcFlux2(solver, U[0], U[1], U[2], U[3], p, nx, ny, F);
                
                for(int nn=0; nn<4; nn++)
                {
                    B->A[nn][mm] = 0.5*(F[nn] - F0[nn])*dS/h;
                }
            }

            double c = gasprop_T2c(solver->gas, E1->P[4]);
            double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;

            if(solver->laminar)
            {
                double r = (E0->P[0] + E1->P[0])*0.5;                
                T = (E0->P[4] + E1->P[4])*0.5;
                double mi = sutherland(T);

                elementCenter(E0, solver->mesh, &x0, &y0);
                elementCenter(E1, solver->mesh, &x1, &y1);
                
                double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr)*dS/d;                    
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
                B = solver->BB[ii];
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
            
            double F0[5];
            double F[5];
            double U[5];

            double rho, u, v, E, p, T;           
            
            rho = solver->U[0][e1];
            u = solver->U[1][e1]/rho;
            v = solver->U[2][e1]/rho;
            E = solver->U[3][e1]/rho;

            implicitAuxCalcFlux_sa2(solver, solver->U[0][e1], solver->U[1][e1], solver->U[2][e1], solver->U[3][e1], solver->U[4][e1], E1->P[3], nx, ny, F0);

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

                implicitAuxCalcFlux_sa2(solver, U[0], U[1], U[2], U[3], U[4], p, nx, ny, F);
                
                for(int nn=0; nn<5; nn++)
                {
                    B->A[nn][mm] = 0.5*(F[nn] - F0[nn])*dS/h;
                }
            }

            double c = gasprop_T2c(solver->gas, E1->P[4]);
            double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;
            
            double r = (E0->P[0] + E1->P[0])*0.5;                
            T = (E0->P[4] + E1->P[4])*0.5;
            double mi = sutherland(T);

            elementCenter(E0, solver->mesh, &x0, &y0);
            elementCenter(E1, solver->mesh, &x1, &y1);
            
            double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr + solver->miT[face1]/solver->Pr_t)*dS/d;
    
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
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        BLOCK* B = NULL;

        for(int kk=0; kk<solver->Nvar; kk++)
        {
            y[kk][ii] = solver->D[ii]*x[kk][ii];
        }

        for(int jj=0; jj<mesh->elemL[ii]->neiN; jj++)
        {
            if(jj==0)
            {
                B = solver->BB[ii];
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
    implicitCalcD(solver);
    if(solver->sa)
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
        aux = solver->dW1;
        solver->dW1 = solver->dW0;
        solver->dW0 = aux;    
        
        implicitMultA(solver, solver->dW0, solver->dW1);
            
        #pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {
                solver->dW1[kk][ii] = w*(-solver->R[kk][ii] - solver->dW1[kk][ii])/solver->D[ii] + solver->dW0[kk][ii];
            }
        }
        
        double ans = fabs(solver->dW1[0][0] - solver->dW0[0][0]);
        for(int ii=1; ii<solver->mesh->Nelem; ii++)
        {
            ans = fmax(ans, fabs(solver->dW1[0][ii] - solver->dW0[0][ii]));
        }
        
        if(ans<1e-14)
        {
            break;
        }
    }
}

void implicitAuxCalcFlux_sst2(SOLVER* solver, double U0, double U1, double U2, double U3, double U4, double U5, double p, double nx, double ny, double* F)
{
    double V = (nx*U1 + ny*U2)/U0;      
    
    F[0] = U0*V;
    F[1] = U1*V + nx*p;    
    F[2] = U2*V + ny*p;
    F[3] = (U3 + p)*V;
    F[4] = U4*V;
    F[5] = U5*V;    
}

void implicitCalcDeltaFlux_sst(SOLVER* solver, double rho, double u, double v, double p, double k, double om, double d0, double d1, double d2, double d3, double d4, double d5, double nx, double ny, double* dF)
{

    double F0[6];
    double F1[6];

    double T = p/(solver->gas->R*rho);
    double e0 = gasprop_T2e(solver->gas, T);
	double E = e0 + (u*u + v*v)*0.5 + k;

    double U0 = rho;
    double U1 = rho*u;
    double U2 = rho*v;
    double U3 = rho*E;  
    double U4 = rho*k;
    double U5 = rho*om;    

    implicitAuxCalcFlux_sst2(solver, U0, U1, U2, U3, U4, U5, p, nx, ny, F0);

    U0 = U0 + d0;
    U1 = U1 + d1;
    U2 = U2 + d2;
    U3 = U3 + d3;
    U4 = U4 + d4;    
    U5 = U5 + d5;    

    rho = U0;
    rho = fmax(rho, solver->rLim);
    U0 = rho;
    
    u = U1/rho;
    v = U2/rho;

    k = U4/rho;
    if(k < 1e-14)
    {
        k = 1e-14;
        U4 = rho*k;
    }

    om = U5/rho;
    if(om < solver->omLim)
    {
        om = solver->omLim;
        U5 = rho*om;
    }

    E = U3/rho;
    double e1 = E - (u*u + v*v)*0.5 - k;
    
    e1 = fmax(e1, 7.176325e+02);

    T = gasprop_e2Taprox(solver->gas, e0, T, e1);
    if(T < 1)
    {
        T = 1;
        e1 = gasprop_T2e(solver->gas, T);
	    E = e1 + (u*u + v*v)*0.5 + k;
        U3 = rho*E;
    }
    
    p = solver->gas->R*T*rho;
    
    implicitAuxCalcFlux_sst2(solver, U0, U1, U2, U3, U4, U5, p, nx, ny, F1);

    for(int kk=0; kk<solver->Nvar; kk++)
    {
        dF[kk] = F1[kk] - F0[kk];
    }
}


void implicitFunc_sst(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW)
{
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
    
    implicitCalcDeltaFlux_sst(solver, E1->P[0], E1->P[1], E1->P[2], E1->P[3], E1->P[5], E1->P[6], dW[0][e1], dW[1][e1], dW[2][e1], dW[3][e1], dW[4][e1], dW[5][e1], nx, ny, dF);

    double T = E1->P[4];
    
    double c = gasprop_T2c(solver->gas, T);
    double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;
    
    double r = E1->P[0];                

    double mi = sutherland(T);

    elementCenter(E0, solver->mesh, &x0, &y0);
    elementCenter(E1, solver->mesh, &x1, &y1);
    
    double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    ra += fmax(4/(3*r), gasprop_T2gamma(solver->gas, T)/r)*(mi/solver->Pr + solver->miT[face1]/solver->Pr_t)*dS/d;
    
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        dW[kk][e0] -= 0.5*(dF[kk]*dS - ra*dW[kk][e1]);
    }
}



