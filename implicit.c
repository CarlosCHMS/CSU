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


void implicitCalcD(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    
    #pragma omp parallel for
    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        double Lc = 0;
        double Lv = 0;
        ELEMENT* E = mesh->elemL[ii];
        solver->D[ii] = 0.0;
        
        for(int jj=0; jj<E->neiN; jj++)
        {
            int face = E->f[jj];
            int face1;
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
                
            nx = dSx/dS;
            ny = dSy/dS;
        
            double r = (E0->P[0] + E1->P[0])*0.5;
            double u = (E0->P[1] + E1->P[1])*0.5;
            double v = (E0->P[2] + E1->P[2])*0.5;
            double p = (E0->P[3] + E1->P[3])*0.5;
            
            double c = sqrt(solver->gamma*p/r);
            
            Lc += (fabs(nx*u + ny*v) + c)*dS;
            
            double T, mi;
            
            if(solver->laminar)
            {
                T = (E0->P[4] + E1->P[4])*0.5;
                mi = sutherland(T);
                Lv += fmax(4/(3*r), solver->gamma/r)*(mi/solver->Pr)*dS*dS;
            }
            
            if(solver->sa)
            {
                T = (E0->P[4] + E1->P[4])*0.5;
                mi = sutherland(T);
                Lv += fmax(4/(3*r), solver->gamma/r)*(mi/solver->Pr + solver->miT[face1]/solver->Pr_t)*dS*dS;
            }
            
        }
        
        solver->dtL[ii] = E->omega/Lc;
        
        solver->D[ii] = 0.5*solver->wImp*Lc;        
        if(solver->laminar)
        {
            Lv /= E->omega;
            solver->D[ii] += Lv;
        }
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

void implicitAuxCalcFlux(SOLVER* solver, double rho, double u, double v, double p, double nx, double ny, double* F)
{
    double V = nx*u + ny*v;
    double E = p/(rho*(solver->gamma-1)) + 0.5*(u*u + v*v);
    double H = E + p/rho;
    
    F[0] = rho*V;
    F[1] = rho*u*V + nx*p;    
    F[2] = rho*v*V + ny*p;
    F[3] = rho*H*V;
}

void implicitCalcDeltaFlux(SOLVER* solver, double rho, double u, double v, double p, double d0, double d1, double d2, double d3, double nx, double ny, double* dF)
{
    double F0[4];
    double F1[4];    

    implicitAuxCalcFlux(solver, rho, u, v, p, nx, ny, F0);

    double E = p/(rho*(solver->gamma-1)) + 0.5*(u*u + v*v);

    double U0 = rho;
    double U1 = rho*u;
    double U2 = rho*v;
    double U3 = rho*E;  

    U0 = U0 + d0;
    U1 = U1 + d1;
    U2 = U2 + d2;
    U3 = U3 + d3;
    
    rho = U0;
    u = U1/rho;
    v = U2/rho;
    E = U3/rho;

    p = (solver->gamma-1)*rho*(E - 0.5*(u*u + v*v));

    implicitAuxCalcFlux(solver, rho, u, v, p, nx, ny, F1);    

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
    
    double c = sqrt(solver->gamma*E1->P[3]/E1->P[0]);
    double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;

    if(solver->laminar)
    {
        double r = (E0->P[0] + E1->P[0])*0.5;                
        double T = (E0->P[4] + E1->P[4])*0.5;
        double mi = sutherland(T);

        elementCenter(E0, solver->mesh, &x0, &y0);
        elementCenter(E1, solver->mesh, &x1, &y1);
        
        double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        ra += fmax(4/(3*r), solver->gamma/r)*(mi/solver->Pr)*dS/d;                    
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
            int face1;
            int e0, e1, p0, p1;
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
                if(solver->sa)
                {
                    implicitFunc_sa(solver, e0, e1, p0, p1, face1, solver->dW0);                
                }
                else
                {
                    implicitFunc(solver, e0, e1, p0, p1, solver->dW0);
                }
            }
        }
        
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->dW0[kk][ii] /= solver->D[ii];
        } 
    }
}

void implicitLUSGS_U(SOLVER* solver)
{
    MESH* mesh = solver->mesh;
    
    for(int ii=mesh->Nelem-1; ii>=0; ii--)
    {
        
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->dW1[kk][ii] = solver->D[ii]*solver->dW0[kk][ii];
        }
        
        ELEMENT* E = mesh->elemL[ii];
        double x0, y0, x1, y1;    
        
        for(int jj=0; jj<E->neiN; jj++)
        {
            int face = E->f[jj];
            int face1;
            int e0, e1, p0, p1;
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
                if(solver->sa)
                {
                    implicitFunc_sa(solver, e0, e1, p0, p1, face1, solver->dW1);                
                }
                else
                {
                    implicitFunc(solver, e0, e1, p0, p1, solver->dW1);
                }                
            }
        }
        
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->dW1[kk][ii] /= solver->D[ii];
        }
    }
}

void implicitAuxCalcFlux_sa(SOLVER* solver, double rho, double u, double v, double p, double n, double nx, double ny, double* F)
{
    double V = nx*u + ny*v;
    double E = p/(rho*(solver->gamma-1)) + 0.5*(u*u + v*v);
    double H = E + p/rho;
    
    F[0] = rho*V;
    F[1] = rho*u*V + nx*p;    
    F[2] = rho*v*V + ny*p;
    F[3] = rho*H*V;
    F[4] = rho*n*V;    
}

void implicitCalcDeltaFlux_sa(SOLVER* solver, double rho, double u, double v, double p, double n, double d0, double d1, double d2, double d3, double d4, double nx, double ny, double* dF)
{
    double F0[5];
    double F1[5];    

    implicitAuxCalcFlux_sa(solver, rho, u, v, p, n, nx, ny, F0);

    double E = p/(rho*(solver->gamma-1)) + 0.5*(u*u + v*v);

    double U0 = rho;
    double U1 = rho*u;
    double U2 = rho*v;
    double U3 = rho*E;  
    double U4 = rho*n;      

    U0 = U0 + d0;
    U1 = U1 + d1;
    U2 = U2 + d2;
    U3 = U3 + d3;
    U4 = U4 + d4;    
    
    rho = U0;
    u = U1/rho;
    v = U2/rho;
    E = U3/rho;
    n = U4/rho;

    p = (solver->gamma-1)*rho*(E - 0.5*(u*u + v*v));

    implicitAuxCalcFlux_sa(solver, rho, u, v, p, n, nx, ny, F1);    

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
    
    double c = sqrt(solver->gamma*E1->P[3]/E1->P[0]);
    double ra = solver->wImp*(fabs(nx*E1->P[1] + ny*E1->P[2]) + c)*dS;
    
    double r = (E0->P[0] + E1->P[0])*0.5;                
    double T = (E0->P[4] + E1->P[4])*0.5;
    double mi = sutherland(T);

    elementCenter(E0, solver->mesh, &x0, &y0);
    elementCenter(E1, solver->mesh, &x1, &y1);
    
    double d = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    ra += fmax(4/(3*r), solver->gamma/r)*(mi/solver->Pr + solver->miT[face1]/solver->Pr_t)*dS/d;                    

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


