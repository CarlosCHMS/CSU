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
#include"gasprop.h"

void saInitU(SOLVER* solver, CONDITION* inside)
{
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        solver->U[4][ii] = inside->Uin[4];
    }
}


void saInterFaceB(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
        double dSx, dSy;
        double x0, y0, x1, y1;

        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];

        ELEMENT* E0 = solver->mesh->elemL[e0];
        ELEMENT* E1 = solver->mesh->elemL[e1];

        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        double duxm = (solver->dPx[1][e0] + solver->dPx[1][e1])*0.5;
        double dvxm = (solver->dPx[2][e0] + solver->dPx[2][e1])*0.5;
        double dTxm = (solver->dPx[3][e0] + solver->dPx[3][e1])*0.5;
        double dnxm = (solver->dPx[4][e0] + solver->dPx[4][e1])*0.5;

        double duym = (solver->dPy[1][e0] + solver->dPy[1][e1])*0.5;
        double dvym = (solver->dPy[2][e0] + solver->dPy[2][e1])*0.5;
        double dTym = (solver->dPy[3][e0] + solver->dPy[3][e1])*0.5;
        double dnym = (solver->dPy[4][e0] + solver->dPy[4][e1])*0.5;

	    elementCenter(E0, solver->mesh, &x0, &y0);
	    elementCenter(E1, solver->mesh, &x1, &y1);

        double dx = x1 - x0;
        double dy = y1 - y0;
        double L = sqrt(dx*dx + dy*dy);

        double dul = (E1->P[1] - E0->P[1])/L;
        double dvl = (E1->P[2] - E0->P[2])/L;
        double dTl = (E1->P[4] - E0->P[4])/L;
        double dnl = (E1->P[5] - E0->P[5])/L;

        double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        double dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        double dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

        double dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
        double dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

        // Flow variables in the face
        double r = (E1->P[0] + E0->P[0])*0.5;
        double u = (E1->P[1] + E0->P[1])*0.5;
        double v = (E1->P[2] + E0->P[2])*0.5;
        double T = (E1->P[4] + E0->P[4])*0.5;
        double n = (E1->P[5] + E0->P[5])*0.5;

        double mi_L = sutherland(T);
        double n_L = mi_L/r;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, r, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*r*n;
        double mi = mi_L + mi_t;
        double k = gasprop_T2Cp(solver->gas, T)*(mi_L/solver->Pr + mi_t/solver->Pr_t);

        solver->miT[ii] = mi_t;

	    double txx = 2*mi*(dux - (dux + dvy)/3);
	    double tyy = 2*mi*(dvy - (dux + dvy)/3);
	    double txy = mi*(duy + dvx);

	    solver->faceFlux[1][ii] = txx*dSx + txy*dSy;
	    solver->faceFlux[2][ii] = txy*dSx + tyy*dSy;
	    solver->faceFlux[3][ii] = (txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy);
	    solver->faceFlux[4][ii] = tx*dSx + ty*dSy;

    }

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->Np; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=1; kk<5; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][face-1];
                }
            }
            else if(face < 0)
            {
                for(int kk=1; kk<5; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][-face-1];
                }
            }
        }
    }
}



void saInterFace(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
        double dSx, dSy;

        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];

        ELEMENT* E0 = solver->mesh->elemL[e0];
        ELEMENT* E1 = solver->mesh->elemL[e1];

        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        double dux = (solver->dPx[1][e0] + solver->dPx[1][e1])*0.5;
        double dvx = (solver->dPx[2][e0] + solver->dPx[2][e1])*0.5;
        double dTx = (solver->dPx[3][e0] + solver->dPx[3][e1])*0.5;
        double dnx = (solver->dPx[4][e0] + solver->dPx[4][e1])*0.5;

        double duy = (solver->dPy[1][e0] + solver->dPy[1][e1])*0.5;
        double dvy = (solver->dPy[2][e0] + solver->dPy[2][e1])*0.5;
        double dTy = (solver->dPy[3][e0] + solver->dPy[3][e1])*0.5;
        double dny = (solver->dPy[4][e0] + solver->dPy[4][e1])*0.5;

        // Flow variables in the face
        double r = (E1->P[0] + E0->P[0])*0.5;
        double u = (E1->P[1] + E0->P[1])*0.5;
        double v = (E1->P[2] + E0->P[2])*0.5;
        double T = (E1->P[4] + E0->P[4])*0.5;
        double n = (E1->P[5] + E0->P[5])*0.5;

        double mi_L = sutherland(T);
        double n_L = mi_L/r;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, r, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*r*n;
        double mi = mi_L + mi_t;
        double k = gasprop_T2Cp(solver->gas, T)*(mi_L/solver->Pr + mi_t/solver->Pr_t);

        solver->miT[ii] = mi_t;

	    double txx = 2*mi*(dux - (dux + dvy)/3);
	    double tyy = 2*mi*(dvy - (dux + dvy)/3);
	    double txy = mi*(duy + dvx);

	    solver->faceFlux[1][ii] = txx*dSx + txy*dSy;
	    solver->faceFlux[2][ii] = txy*dSx + tyy*dSy;
	    solver->faceFlux[3][ii] = (txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy);
	    solver->faceFlux[4][ii] = tx*dSx + ty*dSy;

    }

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->Np; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=1; kk<5; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][face-1];
                }
            }
            else if(face < 0)
            {
                for(int kk=1; kk<5; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][-face-1];
                }
            }
        }
    }
}


void saInterSource(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {

        ELEMENT* E0 = solver->mesh->elemL[ii];

        double drx = solver->dPx[0][ii];
        double dvx = solver->dPx[2][ii];
        double dnx = solver->dPx[4][ii];

        double dry = solver->dPy[0][ii];
        double duy = solver->dPy[1][ii];
        double dny = solver->dPy[4][ii];

        // Flow variables in the face
        double rho = E0->P[0];
        double T = E0->P[4];
        double n = E0->P[5];

        double d = E0->d;
        double mi_L = sutherland(T);
        double n_L = mi_L/rho;
        double S = fabs(duy - dvx);
        double Qt;

        saCalcSource(n, n_L, S, d, rho, drx, dry, dnx, dny, &Qt);
        solver->R[4][ii] -= Qt*E0->omega;

    }
}

void saInter(SOLVER* solver)
{
    if(solver->viscBlazek)
    {
        saInterFaceB(solver);
    }
    else
    {
        saInterFace(solver);
    }
    
    saInterSource(solver);
}


void saCalcFace(double ni, double ni_L, double r, double dnix, double dniy, double* fv1, double* tx, double* ty)
{

    double aux;

    double Cv1 = 7.1;
    double X = ni/ni_L;
    aux = X*X*X;
    *fv1 = aux/(aux + Cv1*Cv1*Cv1);

    double sig = 2./3.;
    *tx = r*(ni_L + ni)*dnix/sig;
    *ty = r*(ni_L + ni)*dniy/sig;

    // Transitional terms not included


}

void saCalcSource(double ni, double ni_L, double S, double d, double rho, double drx, double dry, double dnix, double dniy, double* Qt)
{

/*
    Based on: Allmaras S. R. et all, Modifications and Clarifications for the Implementation
of the Spalart-Allmaras Turbulence Model, 2012
*/

    double aux;

    double Cv1 = 7.1;
    double X = ni/ni_L;
    aux = X*X*X;
    double fv1 = aux/(aux + Cv1*Cv1*Cv1);

    double k = 0.41;
    double fv2 = 1. - X/(1. + X*fv1);
    double Sbar = ni*fv2/(k*k*d*d);
    double Stil;

    double cv2 = 0.7;
    double cv3 = 0.9;
    if(Sbar >= - cv2*S)
    {
        Stil = S + Sbar;
    }
    else
    {
        Stil = S + S*(cv2*cv2*S + cv3*Sbar)/((cv3 - 2*cv2)*S - Sbar);
    }

    if(Stil < 0)
    {
        printf("Stil: %e", Stil);
    }

    double sig = 2./3.;
    double Cb1 = 0.1355;
    double Cb2 = 0.622;
    double Cw1 = Cb1/(k*k) + (1. + Cb2)/sig;

    double Cw2 = 0.3;
    double r = ni/(Stil*k*k*d*d);
    r = fmin(r, 10);
    aux = r*r;
    aux = aux*aux*aux;
    double g = r + Cw2*(aux - r);

    double Cw3 = 2.;
    aux = Cw3*Cw3;
    aux = aux*aux*aux;
    double fw = g*pow((1. + aux)/(g*g*g*g*g*g + aux), 1./6.);

    // Transitional terms not included

    *Qt = rho*(Cb1*Stil*ni + Cb2*(dnix*dnix + dniy*dniy)/sig - Cw1*fw*(ni*ni/(d*d))) - (ni_L + ni)*(dnix*drx + dniy*dry)/sig;

}

void saBoundaryFace(SOLVER* solver, MESHBC* bc)
{
    double f[5];
    double miEddy;
    int e0;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        e0 = bc->elemL[ii]->neiL[0]->ii;

        saBoundaryFaceViscFlux(solver, bc, ii, f, &miEddy);

        solver->R[1][e0] -= f[1];
        solver->R[2][e0] -= f[2];
        solver->R[3][e0] -= f[3];
        solver->R[4][e0] -= f[4];
	}
}

void saBoundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f, double* miEddy)
{

    int e0, p0, p1;
    double x0, x1, y0, y1;
    double dux, duy, dvx, dvy, dTx, dTy, dnx, dny;

    e0 = bc->elemL[ii]->neiL[0]->ii;
    p0 = bc->elemL[ii]->p[0];
    p1 = bc->elemL[ii]->p[1];

    ELEMENT* E0 = bc->elemL[ii]->neiL[0];

    if(bc->flagBC == 0)
    {
        //symmetry
       	double nx, ny, dS;
        meshCalcDS2(solver->mesh, p0, p1, &nx, &ny, &dS);
        
        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];
        
        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];

        dux = duxm - (duxm*nx + duym*ny)*nx;
        duy = duym - (duxm*nx + duym*ny)*ny;        

        dvx = dvxm - (dvxm*nx + dvym*ny)*nx;
        dvy = dvym - (dvxm*nx + dvym*ny)*ny;        
        
        dTx = dTxm - (dTxm*nx + dTym*ny)*nx;
        dTy = dTym - (dTxm*nx + dTym*ny)*ny;        
        
        // Flow variables in the face
        double rho = E0->P[0];
        double u = E0->P[1];
        double v = E0->P[2];
        double T = E0->P[4];
        double n = E0->P[5];

        double mi_L = sutherland(T);
        double n_L = mi_L/rho;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*rho*n;
        double mi = mi_L + mi_t;
        double k = gasprop_T2Cp(solver->gas, T)*(mi_L/solver->Pr + mi_t/solver->Pr_t);          
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx);  
        
        f[1] = (txx*nx + txy*ny)*dS;
        f[2] = (txy*nx + tyy*ny)*dS;
        f[3] = (u*(txx*nx + txy*ny) + v*(txy*nx + tyy*ny) + k*(dTx*nx + dTy*ny))*dS;
        f[4] = 0.0;
        *miEddy = mi_t;
        
    }
    else if(bc->flagBC == 1)
    {
        //inlet
        double dSx, dSy;
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        elementCenter(E0, solver->mesh, &x0, &y0);

       	x1 = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
        y1 = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

        double dx = x1 - x0;
        double dy = y1 - y0;
        double L = sqrt(dx*dx + dy*dy);

        double dul = (solver->inlet->Pin[1] - E0->P[1])/L;
        double dvl = (solver->inlet->Pin[2] - E0->P[2])/L;
        double dTl = (solver->inlet->Pin[4] - E0->P[4])/L;
        double dnl = (solver->inlet->Pin[5] - E0->P[5])/L;

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];
        double dnxm = solver->dPx[4][e0];

        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];
        double dnym = solver->dPy[4][e0];

        dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

        dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
        dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

        // Flow variables in the face
        double rho = solver->inlet->Pin[0];
        double u = solver->inlet->Pin[1];
        double v = solver->inlet->Pin[2];
        double T = solver->inlet->Pin[4];
        double n = solver->inlet->Pin[5];

        double mi_L = sutherland(T);
        double n_L = mi_L/rho;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*rho*n;
        double mi = mi_L + mi_t;
        double k = gasprop_T2Cp(solver->gas, T)*(mi_L/solver->Pr + mi_t/solver->Pr_t);

        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);
        double txy = mi*(duy + dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy);
        f[4] = tx*dSx + ty*dSy;
        *miEddy = mi_t;
        
    }
    else if(bc->flagBC == 3)
    {
        //wall
        double dSx, dSy;
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        elementCenter(E0, solver->mesh, &x0, &y0);

       	x1 = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
        y1 = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

        double dx = x1 - x0;
        double dy = y1 - y0;
        double L = sqrt(dx*dx + dy*dy);

        double dul = (0 - E0->P[1])/L;
        double dvl = (0 - E0->P[2])/L;
        double dnl = (0 - E0->P[5])/L;

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dnxm = solver->dPx[4][e0];

        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dnym = solver->dPy[4][e0];

        dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
        dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

        // Flow variables in the face
        double rho = E0->P[0];
        double T = E0->P[4];
        double n = 0.0;

        double mi_L = sutherland(T);
        double n_L = mi_L/rho;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

        double mi = mi_L;

        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);
        double txy = mi*(duy + dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = 0.0;
        f[4] = tx*dSx + ty*dSy;
        *miEddy = 0.0;
        
    }
    else if(bc->flagBC == 4)
    {
        //wallT
        double dSx, dSy;
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        elementCenter(E0, solver->mesh, &x0, &y0);

       	x1 = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
        y1 = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

        double dx = x1 - x0;
        double dy = y1 - y0;
        double L = sqrt(dx*dx + dy*dy);

        double dul = (0 - E0->P[1])/L;
        double dvl = (0 - E0->P[2])/L;
        double dTl = (solver->Twall - E0->P[4])/L;            
        double dnl = (0 - E0->P[5])/L;

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];
        double dnxm = solver->dPx[4][e0];

        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];
        double dnym = solver->dPy[4][e0];

        dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

        dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
        dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

        // Flow variables in the face
        double rho = E0->P[0];
        double T = solver->Twall;
        double n = 0.0;        

        double mi_L = sutherland(T);
        double n_L = mi_L/rho;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*rho*n;
        double mi = mi_L + mi_t;
        double k = gasprop_T2Cp(solver->gas, T)*(mi_L/solver->Pr + mi_t/solver->Pr_t);

        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);
        double txy = mi*(duy + dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = k*(dTx*dSx + dTy*dSy);	        
        f[4] = tx*dSx + ty*dSy;        
        *miEddy = mi_t;
        
    }
    else
    {

        //outlet
        double dSx, dSy;
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        dux = solver->dPx[1][e0];
        dvx = solver->dPx[2][e0];
        dTx = solver->dPx[3][e0];
        dnx = solver->dPx[4][e0];

        duy = solver->dPy[1][e0];
        dvy = solver->dPy[2][e0];
        dTy = solver->dPy[3][e0];
        dny = solver->dPy[4][e0];

        // Flow variables in the face
        double rho = E0->P[0];
        double u = E0->P[1];
        double v = E0->P[2];
        double T = E0->P[4];
        double n = E0->P[5];

        double mi_L = sutherland(T);
        double n_L = mi_L/rho;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*rho*n;
        double mi = mi_L + mi_t;
        double k = gasprop_T2Cp(solver->gas, T)*(mi_L/solver->Pr + mi_t/solver->Pr_t);

        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);
        double txy = mi*(duy + dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy);
        f[4] = tx*dSx + ty*dSy;
        *miEddy = mi_t;
        
    }
}


void saBoundary(SOLVER* solver)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        saBoundaryFace(solver, solver->mesh->bc[ii]);
    }
}

void saCalcD(MESH* mesh)
{

    double xVol, yVol;

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        ELEMENT* Evol = mesh->elemL[ii];
        Evol->d = 1e6;
    }

    for(int jj=0; jj<mesh->Nmark; jj++)
    {

        //printf("%i\n", mesh->bc[jj]->flagBC);
        if(mesh->bc[jj]->flagBC==3 || mesh->bc[jj]->flagBC==4)
        {

            for(int ii=0; ii<mesh->Nelem; ii++)
            {
                ELEMENT* Evol = mesh->elemL[ii];
                elementCenter(Evol, mesh, &xVol, &yVol);

                for(int kk=0; kk<mesh->bc[jj]->Nelem; kk++)
                {
                    ELEMENT* Esurf = mesh->bc[jj]->elemL[kk];
                    double p0x = mesh->p[Esurf->p[0]][0];
                    double p0y = mesh->p[Esurf->p[0]][1];

                    double p1x = mesh->p[Esurf->p[1]][0];
                    double p1y = mesh->p[Esurf->p[1]][1];

                    double num = (p1x - p0x)*(xVol - p0x) + (p1y - p0y)*(yVol - p0y);
                    double den = (p1x - p0x)*(p1x - p0x) + (p1y - p0y)*(p1y - p0y);

                    double t = num/den;

                    if(t>1.)
                    {
                        t = 1.;
                    }
                    else if(t<0.)
                    {
                        t = 0.;
                    }

                    double dx = (p1x - p0x)*t + p0x - xVol;
                    double dy = (p1y - p0y)*t + p0y - yVol;

                    double d = sqrt(dx*dx + dy*dy);

                    if(Evol->d > d)
                    {
                        Evol->d = d;
                    }

                }

                //printf("\ntf: %f\n", tf);

            }
        }
    }

    /*
    for(int jj=0; jj<mesh->Nmark; jj++)
    {

        //printf("%i\n", mesh->bc[jj]->flagBC);
        if(mesh->bc[jj]->flagBC==3)
        {

            for(int ii=0; ii<mesh->Nelem; ii++)
            {

                for(int kk=0; kk<mesh->bc[jj]->Nelem; kk++)
                {
                    ELEMENT* Esurf = mesh->bc[jj]->elemL[kk];
                    ELEMENT* Evol = Esurf->neiL[0];
                    printf("%e\n", Evol->d);

                }

            }

        }

    }
    */

}

void saCalcTensorWall(SOLVER* solver, ELEMENT* E, double* Txx, double* Txy, double* Tyy, double* x, double* yp)
{

    double x0, y0, x1, y1, dSx, dSy;

    int e0 = E->neiL[0]->ii;
    int p0 = E->p[0];
    int p1 = E->p[1];

    ELEMENT* E0 = E->neiL[0];

    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    double dS = sqrt(dSx*dSx + dSy*dSy);


    elementCenter(E0, solver->mesh, &x0, &y0);
    elementCenter(E, solver->mesh, &x1, &y1);

    double dx = x1 - x0;
    double dy = y1 - y0;
    double L = sqrt(dx*dx + dy*dy);

    double dul = (0 - E0->P[1])/L;
    double dvl = (0 - E0->P[2])/L;

    double duxm = solver->dPx[1][e0];
    double dvxm = solver->dPx[2][e0];

    double duym = solver->dPy[1][e0];
    double dvym = solver->dPy[2][e0];

    double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
    double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

    double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
    double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

    double T = E0->P[4];
    double mi = sutherland(T);

    //printf("%f\n", E0->P[1]);

    double txx = 2*mi*(dux - (dux + dvy)/3);
    double tyy = 2*mi*(dvy - (dux + dvy)/3);
    double txy = mi*(duy + dvx);

    double fx = (txx*dSx + txy*dSy)/dS;
    double fy = (txy*dSx + tyy*dSy)/dS;
    double ft = (fx*dSy - fy*dSx)/dS;

    double ut = sqrt(fabs(ft)/E0->P[0]);

    *Txx = 0.0;
	*Txy = ft;
	*Tyy = 0.0;
	*x = x1;
	*yp = E0->P[0]*E0->d*ut/mi;
}

void saSolverWriteSurf(SOLVER* solver)
{

    MESHBC* bc;
    char s[50];
        
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "surfData.csv");        
    FILE* ff = fopen(s, "w");   

    double rin = solver->inlet->Pin[0];
    double uin = solver->inlet->Pin[1];
    double vin = solver->inlet->Pin[2];
    double pin = solver->inlet->Pin[3];
    double qdin = 0.5*rin*(uin*uin + vin*vin);


    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {
        bc = solver->mesh->bc[jj];
        if(strcmp(bc->name, solver->writeSurf) == 0)
        {
            fprintf(ff, "x,y,rho,u,v,p,T,n,Cp,Cfx,Cfy,q,yplus,mach,miEddy,\n");
        
            int Nvar = 13;
            int Nelem = bc->Nelem;
            double** D = tableMallocDouble(Nvar, Nelem);
            double** Dp = tableMallocDouble(Nvar, Nelem+1);
            int* d = malloc((Nelem+1)*sizeof(int));
            int* pn = malloc((Nelem+1)*sizeof(int));
            int* pn2 = malloc((solver->mesh->Np)*sizeof(int));
            double mi;

            for(int ii=0; ii<bc->Nelem+1; ii++)
            {                
                for(int kk=0; kk<Nvar; kk++)
                {
                    Dp[kk][ii] = 0.0;
                } 
                d[ii] = 0;                 
            }
            
            for(int ii=0; ii<bc->Nelem; ii++)
            {            
                pn[ii] = bc->elemL[ii]->p[0];                
            }
            pn[Nelem] = bc->elemL[Nelem-1]->p[1];

            for(int ii=0; ii<bc->Nelem+1; ii++)
            {            
                pn2[pn[ii]] = ii;                
            }
        
            double f[5];
        
            for(int ii=0; ii<bc->Nelem; ii++)
            {                                
                for(int kk=0; kk<Nvar; kk++)
                {
                    int p0 = bc->elemL[ii]->p[0];
                    int p1 = bc->elemL[ii]->p[1];
                    double dSx, dSy, dS;
                    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
                
                    dS = sqrt(dSx*dSx + dSy*dSy);
                
                    double r = bc->elemL[ii]->P[0];
                    double u = bc->elemL[ii]->P[1];
                    double v = bc->elemL[ii]->P[2];                    
                    double p = bc->elemL[ii]->P[3];
                    double T = bc->elemL[ii]->P[4];                    
                                
                    saBoundaryFaceViscFlux(solver, bc, ii, f, &mi);

                    double Cp = (p - pin)/qdin;                
                    double Cfx = f[1]/(dS*qdin);
                    double Cfy = f[2]/(dS*qdin);
                    double q = (f[3] - u*f[1] - v*f[2])/dS;
                    double tau = sqrt(f[1]*f[1] + f[2]*f[2])/dS;
                    double uplus = sqrt(tau/r);
                    double yplus = fabs(r*uplus*bc->elemL[ii]->neiL[0]->d)/sutherland(T);
                    double c = gasprop_T2c(solver->gas, T);
                    double mach = sqrt(u*u + v*v)/c;
                    
                    D[0][ii] = r;
                    D[1][ii] = u;
                    D[2][ii] = v;
                    D[3][ii] = p;
                    D[4][ii] = T;
                    D[5][ii] = bc->elemL[ii]->P[5];
                    D[6][ii] = Cp;
                    D[7][ii] = Cfx;
                    D[8][ii] = Cfy;
                    D[9][ii] = q;
                    D[10][ii] = yplus;
                    D[11][ii] = mach;
                    D[12][ii] = mi;
                }
            }
            
            for(int ii=0; ii<bc->Nelem; ii++)
            {                
                int p0 = bc->elemL[ii]->p[0];
                int p1 = bc->elemL[ii]->p[1];                
                for(int kk=0; kk<Nvar; kk++)
                {
                    Dp[kk][pn2[p0]] += D[kk][ii];
                    Dp[kk][pn2[p1]] += D[kk][ii];                    
                }
                d[pn2[p0]] += 1;
                d[pn2[p1]] += 1;                                               
            }
            
            for(int ii=0; ii<bc->Nelem+1; ii++)
            {                
                for(int kk=0; kk<Nvar; kk++)
                {
                    Dp[kk][ii] /= d[ii];
                }
            }
            
            for(int ii=0; ii<bc->Nelem+1; ii++)
            {
                fprintf(ff, "% .10e,", solver->mesh->p[pn[ii]][0]);
                fprintf(ff, "% .10e,", solver->mesh->p[pn[ii]][1]);
                for(int kk=0; kk<Nvar; kk++)
                {
                    fprintf(ff, "% .10e,", Dp[kk][ii]);
                }
                fprintf(ff, "\n");
            }
            
            tableFreeDouble(D, Nvar);
            tableFreeDouble(Dp, Nvar);            
            free(d);
            free(pn);
            free(pn2);            
        }
    }

    fclose(ff);
}



