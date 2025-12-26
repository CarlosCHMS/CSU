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
#include"sst.h"
#include"gasprop.h"
#include"sstTrans.h"

SST* sstInit()
{
    SST* sst = malloc(sizeof(SST));
    sst->sk1 = 0.85;
    sst->so1 = 0.5;
    sst->b1 = 0.075;
    sst->sk2 = 1.0;
    sst->so2 = 0.856;
    sst->b2 = 0.0828;
    sst->bs = 0.09;
    sst->a1 = 0.31;

    double k = 0.41;
    
    sst->g1 = sst->b1/sst->bs - sst->so1*k*k/sqrt(sst->bs);
    
    sst->g2 = sst->b2/sst->bs - sst->so2*k*k/sqrt(sst->bs);
        
    sst->trans = sstTransInit();
    
    return sst;
}

void sstFree(SST* sst)
{
    sstTransFree(sst->trans);
    free(sst);
}

void sstInitU(SOLVER* solver, CONDITION* inside)
{
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        solver->U[4][ii] = inside->Uin[4];
        solver->U[5][ii] = inside->Uin[5];
    }
}

void sstInterFace(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
        double dSx, dSy;

        SSTVAR var;

        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];

        ELEMENT* E0 = solver->mesh->elemL[e0];
        ELEMENT* E1 = solver->mesh->elemL[e1];

        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        /*
        var.dux = (solver->dPx[1][e0] + solver->dPx[1][e1])*0.5;
        var.dvx = (solver->dPx[2][e0] + solver->dPx[2][e1])*0.5;
        double dTx = (solver->dPx[3][e0] + solver->dPx[3][e1])*0.5;
        var.dkx = (solver->dPx[4][e0] + solver->dPx[4][e1])*0.5;
        var.dox = (solver->dPx[5][e0] + solver->dPx[5][e1])*0.5;        

        var.duy = (solver->dPy[1][e0] + solver->dPy[1][e1])*0.5;
        var.dvy = (solver->dPy[2][e0] + solver->dPy[2][e1])*0.5;
        double dTy = (solver->dPy[3][e0] + solver->dPy[3][e1])*0.5;
        var.dky = (solver->dPy[4][e0] + solver->dPy[4][e1])*0.5;
        var.doy = (solver->dPy[5][e0] + solver->dPy[5][e1])*0.5;        
        */

        double x0, x1, y0, y1;

        elementCenter(E0, solver->mesh, &x0, &y0);
        elementCenter(E1, solver->mesh, &x1, &y1);        

        double dx = x1 - x0;
        double dy = y1 - y0;
        double L = sqrt(dx*dx + dy*dy);
        double nx = dx/L;
        double ny = dy/L;        

        double dul = (E1->P[1] - E0->P[1])/L;
        double dvl = (E1->P[2] - E0->P[2])/L;
        double dTl = (E1->P[4] - E0->P[4])/L;
        double dkl = (E1->P[5] - E0->P[5])/L;
        double dol = (E1->P[6] - E0->P[6])/L;        

        double duxm = 0.5*(solver->dPx[1][e0] + solver->dPx[1][e1]);
        double dvxm = 0.5*(solver->dPx[2][e0] + solver->dPx[2][e1]);
        double dTxm = 0.5*(solver->dPx[3][e0] + solver->dPx[3][e1]);
        double dkxm = 0.5*(solver->dPx[4][e0] + solver->dPx[4][e1]);
        double doxm = 0.5*(solver->dPx[5][e0] + solver->dPx[5][e1]);

        double duym = 0.5*(solver->dPy[1][e0] + solver->dPy[1][e1]);
        double dvym = 0.5*(solver->dPy[2][e0] + solver->dPy[2][e1]);
        double dTym = 0.5*(solver->dPy[3][e0] + solver->dPy[3][e1]);
        double dkym = 0.5*(solver->dPy[4][e0] + solver->dPy[4][e1]);
        double doym = 0.5*(solver->dPy[5][e0] + solver->dPy[5][e1]);

        double aux;

        aux = duxm*nx + duym*ny;
        var.dux = duxm + (dul - aux)*nx;
        var.duy = duym + (dul - aux)*ny;

        aux = dvxm*nx + dvym*ny;
        var.dvx = dvxm + (dvl - aux)*nx;
        var.dvy = dvym + (dvl - aux)*ny;

        aux = dTxm*nx + dTym*ny;
        double dTx = dTxm + (dTl - aux)*nx;
        double dTy = dTym + (dTl - aux)*ny;

        aux = dkxm*nx + dkym*ny;
        var.dkx = dkxm + (dkl - aux)*nx;
        var.dky = dkym + (dkl - aux)*ny;

        aux = doxm*nx + doym*ny;
        var.dox = doxm + (dol - aux)*nx;
        var.doy = doym + (dol - aux)*ny;

        // Flow variables in the face
        var.r = (E1->P[0] + E0->P[0])*0.5;
        double u = (E1->P[1] + E0->P[1])*0.5;
        double v = (E1->P[2] + E0->P[2])*0.5;
        var.T = (E1->P[4] + E0->P[4])*0.5;
        var.k = (E1->P[5] + E0->P[5])*0.5;
        var.om = (E1->P[6] + E0->P[6])*0.5;        
        var.d = (E1->d + E0->d)*0.5;

        var.mi_L = sutherland(var.T);

        sstFlux(solver->sst, &var);
        
        double mi = var.mi_L + var.mi_t;
        double kk = gasprop_T2Cp(solver->gas, var.T)*(var.mi_L/solver->Pr + var.mi_t/solver->Pr_t);

        solver->miT[ii] = var.mi_t;

	    double txx = 2*mi*(var.dux - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double tyy = 2*mi*(var.dvy - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double txy = mi*(var.duy + var.dvx);
        
	    solver->faceFlux[1][ii] = txx*dSx + txy*dSy;
	    solver->faceFlux[2][ii] = txy*dSx + tyy*dSy;
	    solver->faceFlux[3][ii] = (txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + kk*(dTx*dSx + dTy*dSy);
	    solver->faceFlux[4][ii] = var.tkx*dSx + var.tky*dSy;
	    solver->faceFlux[5][ii] = var.tox*dSx + var.toy*dSy;
     
    }

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->neiN; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=1; kk<solver->Nvar; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][face-1];
                }
            }
            else if(face < 0)
            {
                for(int kk=1; kk<solver->Nvar; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][-face-1];
                }
            }
        }
    }
}


void sstInter(SOLVER* solver)
{
    sstInterFace(solver);
    sstInterSource(solver);
}

void sstBoundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f, double* miEddy)
{

    int e0, p0, p1;
    double x0, x1, y0, y1;
    double dTx, dTy;

    SSTVAR var;

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
        double dkxm = solver->dPx[4][e0];
        double doxm = solver->dPx[5][e0];
        
        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];
        double dkym = solver->dPy[4][e0];
        double doym = solver->dPy[5][e0];

        var.dux = duxm - (duxm*nx + duym*ny)*nx;
        var.duy = duym - (duxm*nx + duym*ny)*ny;        

        var.dvx = dvxm - (dvxm*nx + dvym*ny)*nx;
        var.dvy = dvym - (dvxm*nx + dvym*ny)*ny;        
        
        dTx = dTxm - (dTxm*nx + dTym*ny)*nx;
        dTy = dTym - (dTxm*nx + dTym*ny)*ny;        

        var.dkx = dkxm - (dkxm*nx + dkym*ny)*nx;
        var.dky = dkym - (dkxm*nx + dkym*ny)*ny;        

        var.dox = doxm - (doxm*nx + doym*ny)*nx;
        var.doy = doym - (doxm*nx + doym*ny)*ny;        

        // Flow variables in the face
        var.r = E0->P[0];
        double u = E0->P[1];
        double v = E0->P[2];
        var.T = E0->P[4];
        var.k = E0->P[5];
        var.om = E0->P[6];
        var.d = E0->d;

        var.mi_L = sutherland(var.T);

        sstFlux(solver->sst, &var);
                
        double mi = var.mi_L + var.mi_t;

	    double txx = 2*mi*(var.dux - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double tyy = 2*mi*(var.dvy - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double txy = mi*(var.duy + var.dvx);
        
        f[1] = (txx*nx + txy*ny)*dS;
        f[2] = (txy*nx + tyy*ny)*dS;
        f[3] = (u*(txx*nx + txy*ny) + v*(txy*nx + tyy*ny))*dS;
        f[4] = 0.0;
        f[5] = 0.0;
        *miEddy = var.mi_t;
        
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
        double dkl = (solver->inlet->Pin[5] - E0->P[5])/L;
        double dol = (solver->inlet->Pin[6] - E0->P[6])/L;        

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];
        double dkxm = solver->dPx[4][e0];
        double doxm = solver->dPx[5][e0];

        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];
        double dkym = solver->dPy[4][e0];
        double doym = solver->dPy[5][e0];        

        var.dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        var.duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        var.dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        var.dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

        var.dkx = dkxm + (dkl - (dkxm*dx + dkym*dy)/L)*dx/L;
        var.dky = dkym + (dkl - (dkxm*dx + dkym*dy)/L)*dy/L;

        var.dox = doxm + (dol - (doxm*dx + doym*dy)/L)*dx/L;
        var.doy = doym + (dol - (doxm*dx + doym*dy)/L)*dy/L;


        // Flow variables in the face
        var.r = solver->inlet->Pin[0];
        double u = solver->inlet->Pin[1];
        double v = solver->inlet->Pin[2];
        var.T = solver->inlet->Pin[4];
        var.k = solver->inlet->Pin[5];
        var.om = solver->inlet->Pin[6];

        var.d = E0->d;
        var.mi_L = sutherland(var.T);

        sstFlux(solver->sst, &var);
        
        double mi = var.mi_L + var.mi_t;
        double kk = gasprop_T2Cp(solver->gas, var.T)*(var.mi_L/solver->Pr + var.mi_t/solver->Pr_t);

	    double txx = 2*mi*(var.dux - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double tyy = 2*mi*(var.dvy - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double txy = mi*(var.duy + var.dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + kk*(dTx*dSx + dTy*dSy);
        f[4] = var.tkx*dSx + var.tky*dSy;
        f[5] = var.tox*dSx + var.toy*dSy;        
        *miEddy = var.mi_t;
        
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

        // Flow variables in the face
        var.r = E0->P[0];
        var.T = E0->P[4];        
        var.k = E0->P[5];
        var.om = E0->P[6];        
        var.d = E0->d;

        var.mi_L = sutherland(var.T);
        double n_L = var.mi_L/var.r;

        double owall = solver->sst->oWallFactor*6*n_L/(solver->sst->b1*var.d*var.d);
        //printf("\n2, %e, %e, %e, %e", owall, E0->P[6], n_L, var.d);

        double dul = (0 - E0->P[1])/L;
        double dvl = (0 - E0->P[2])/L;
        double dkl = (0 - E0->P[5])/L;
        double dol = (owall - E0->P[6])/L;        

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dkxm = solver->dPx[4][e0];
        double doxm = solver->dPx[5][e0];

        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dkym = solver->dPy[4][e0];
        double doym = solver->dPy[5][e0];

        var.dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        var.duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        var.dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        var.dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        var.dkx = dkxm + (dkl - (dkxm*dx + dkym*dy)/L)*dx/L;
        var.dky = dkym + (dkl - (dkxm*dx + dkym*dy)/L)*dy/L;

        var.dox = doxm + (dol - (doxm*dx + doym*dy)/L)*dx/L;
        var.doy = doym + (dol - (doxm*dx + doym*dy)/L)*dy/L;        

        sstFlux(solver->sst, &var);
        
        double mi = var.mi_L + var.mi_t;

	    double txx = 2*mi*(var.dux - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double tyy = 2*mi*(var.dvy - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double txy = mi*(var.duy + var.dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = 0.0;
        f[4] = var.tkx*dSx + var.tky*dSy;
        f[5] = var.tox*dSx + var.toy*dSy;        
        *miEddy = var.mi_t;
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

        // Flow variables in the face
        var.r = E0->P[0];
        var.T = E0->P[4];        
        var.k = E0->P[5];
        var.om = E0->P[6];        
        var.d = E0->d;

        var.mi_L = sutherland(var.T);
        double n_L = var.mi_L/var.r;

        double owall = solver->sst->oWallFactor*6*n_L/(solver->sst->b1*var.d*var.d);
        
        double dul = (0 - E0->P[1])/L;
        double dvl = (0 - E0->P[2])/L;
        double dTl = (solver->Twall - E0->P[4])/L;            
        double dkl = (0 - E0->P[5])/L;
        double dol = (owall - E0->P[6])/L;        

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];        
        double dkxm = solver->dPx[4][e0];
        double doxm = solver->dPx[5][e0];

        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];                
        double dkym = solver->dPy[4][e0];
        double doym = solver->dPy[5][e0];

        var.dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        var.duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        var.dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        var.dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

        var.dkx = dkxm + (dkl - (dkxm*dx + dkym*dy)/L)*dx/L;
        var.dky = dkym + (dkl - (dkxm*dx + dkym*dy)/L)*dy/L;

        var.dox = doxm + (dol - (doxm*dx + doym*dy)/L)*dx/L;
        var.doy = doym + (dol - (doxm*dx + doym*dy)/L)*dy/L;        

        sstFlux(solver->sst, &var);
        
        double mi = var.mi_L + var.mi_t;
        double kk = gasprop_T2Cp(solver->gas, var.T)*(var.mi_L/solver->Pr + var.mi_t/solver->Pr_t);

	    double txx = 2*mi*(var.dux - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double tyy = 2*mi*(var.dvy - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double txy = mi*(var.duy + var.dvx);

        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = kk*(dTx*dSx + dTy*dSy);
        f[4] = var.tkx*dSx + var.tky*dSy;
        f[5] = var.tox*dSx + var.toy*dSy;        
        *miEddy = var.mi_t;        
    }
    else
    {

        //outlet
        double dSx, dSy;
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        var.dux = solver->dPx[1][e0];
        var.dvx = solver->dPx[2][e0];
        dTx = solver->dPx[3][e0];
        var.dkx = solver->dPx[4][e0];
        var.dox = solver->dPx[5][e0];        

        var.duy = solver->dPy[1][e0];
        var.dvy = solver->dPy[2][e0];
        dTy = solver->dPy[3][e0];
        var.dky = solver->dPy[4][e0];
        var.doy = solver->dPy[5][e0];

        // Flow variables in the face
        var.r = E0->P[0];
        double u = E0->P[1];
        double v = E0->P[2];
        var.T = E0->P[4];
        var.k = E0->P[5];
        var.om = E0->P[6];        
        var.d = E0->d;

        var.mi_L = sutherland(var.T);

        sstFlux(solver->sst, &var);
        
        double mi = var.mi_L + var.mi_t;
        double kk = gasprop_T2Cp(solver->gas, var.T)*(var.mi_L/solver->Pr + var.mi_t/solver->Pr_t);

	    double txx = 2*mi*(var.dux - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double tyy = 2*mi*(var.dvy - (var.dux + var.dvy)/3) - 2.*var.r*var.k/3.;
	    double txy = mi*(var.duy + var.dvx);
	    
        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + kk*(dTx*dSx + dTy*dSy);
        f[4] = var.tkx*dSx + var.tky*dSy;
        f[5] = var.tox*dSx + var.toy*dSy;
        *miEddy = var.mi_t;
        
    }
}

void sstBoundaryFace(SOLVER* solver, MESHBC* bc)
{
    double f[6];
    double miEddy;
    int e0;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        e0 = bc->elemL[ii]->neiL[0]->ii;

        sstBoundaryFaceViscFlux(solver, bc, ii, f, &miEddy);

        solver->R[1][e0] -= f[1];
        solver->R[2][e0] -= f[2];
        solver->R[3][e0] -= f[3];
        solver->R[4][e0] -= f[4];
        solver->R[5][e0] -= f[5];        
	}
}

void sstBoundary(SOLVER* solver)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        sstBoundaryFace(solver, solver->mesh->bc[ii]);
    }
}

double flim(double x)
{
    double aux1 = x;
    double L = 1e-10;
    if(fabs(x) < L)
    {
        if(x >= 0)
        {
            aux1 = L;
        }
        else
        {
            aux1 = -L;
        }
    }
    return aux1;
}

void sstInterSource(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        ELEMENT* E0 = solver->mesh->elemL[ii];

        SSTVAR var;

        var.dux = solver->dPx[1][ii];
        var.dvx = solver->dPx[2][ii];
        var.dkx = solver->dPx[4][ii];
        var.dox = solver->dPx[5][ii];

        var.duy = solver->dPy[1][ii];
        var.dvy = solver->dPy[2][ii];        
        var.dky = solver->dPy[4][ii];
        var.doy = solver->dPy[5][ii];        

        double x, y;
        elementCenter(E0, solver->mesh, &x, &y);
        var.x = x;
        var.y = y;
        var.ii = E0->ii;

        // Flow variables in the face
        var.r = E0->P[0];
        var.T = E0->P[4];
        var.k = E0->P[5];
        var.om = E0->P[6];
        var.d = E0->d;
        var.l = sqrt(E0->omega);

        var.mi_L = sutherland(var.T);

        if(var.om == solver->omLim)
        {
            solver->dQkdr[ii] = 0;
            solver->dQkdrk[ii] = 0;
            solver->dQkdro[ii] = 0;

            solver->dQodr[ii] = 0;
            solver->dQodrk[ii] = 0;
            solver->dQodro[ii] = 0;
        }
        else
        {
            sstSources(solver->sst, &var);

            solver->R[4][ii] -= var.Qtk*E0->omega;
            solver->R[5][ii] -= var.Qto*E0->omega;

            double dQkdr = var.Qtk;
            double dQkdo = var.Qtk;
            double dQkdk = var.Qtk;        
            double dQodr = var.Qto;        
            double dQodo = var.Qto;
            double dQodk = var.Qto;
            
            double h = 1e-6;
            double dh = 1e-6;
            double hh;
            
            hh = fabs(var.om*h);
            hh = fmax(hh, dh);
            var.om += hh;
            sstSources(solver->sst, &var);
            dQkdo = (var.Qtk - dQkdo)/hh;    
            dQodo = (var.Qto - dQodo)/hh;
            var.om = E0->P[6];
            
            hh = fabs(var.k*h);
            hh = fmax(hh, dh);
            var.k += hh;
            sstSources(solver->sst, &var);
            dQkdk = (var.Qtk - dQkdk)/hh; 
            dQodk = (var.Qto - dQodk)/hh;
            var.k = E0->P[5];
            
            hh = fabs(var.r*h);
            hh = fmax(hh, dh);        
            var.r += hh;
            sstSources(solver->sst, &var);
            dQkdr = (var.Qtk - dQkdr)/hh;        
            dQodr = (var.Qto - dQodr)/hh;
            var.r = E0->P[0];

            solver->dQkdr[ii] = (dQkdr - dQkdk*var.k/var.r - dQkdo*var.om/var.r)*E0->omega;//var.dQkdr*E0->omega;
            solver->dQkdrk[ii] = dQkdk*E0->omega/var.r;//dQkdk*E0->omega/var.r;//var.dQkdrk*E0->omega;        
            solver->dQkdro[ii] = dQkdo*E0->omega/var.r;//var.dQkdro*E0->omega;

            solver->dQodr[ii] = (dQodr - dQodk*var.k/var.r - dQodo*var.om/var.r)*E0->omega;//var.dQodr*E0->omega;
            solver->dQodrk[ii] = dQodk*E0->omega/var.r;//var.dQodrk*E0->omega;        
            solver->dQodro[ii] = dQodo*E0->omega/var.r;//var.dQodro*E0->omega;
        }
    }
}

double sstBlend(double x1, double x2, double F1)
{
    return x1*F1 + x2*(1-F1);
}

void sstSources(SST* sst, SSTVAR* var)
{
    if(sst->trans->flag)
    {
        sstTransSources(sst->trans, var);
    }
    else
    {
        
        double n_L = var->mi_L/var->r;
        double omega = fabs(var->duy - var->dvx);

        double sqrtk_term = sqrt(var->k)/(sst->bs*var->om*var->d);
        
        double n_L_term = 500*n_L/(var->d*var->d*var->om);

        double F1 = sstF1(sst, var, n_L_term, sqrtk_term);

        double F2 = sstF2(sst, var, n_L_term, sqrtk_term);

        var->mi_t = var->r*sst->a1*var->k/fmax(sst->a1*var->om, omega*F2);
        //var->mi_t = fmin(var->mi_t, 1e6*var->mi_L);
        
        double n_t = var->mi_t/var->r;

        double g = sstBlend(sst->g1, sst->g2, F1);
        double b = sstBlend(sst->b1, sst->b2, F1);

        double txx = 2*var->mi_t*(var->dux - (var->dux + var->dvy)/3) - 2.*var->r*var->k/3;
	    double tyy = 2*var->mi_t*(var->dvy - (var->dux + var->dvy)/3) - 2.*var->r*var->k/3;
	    double txy = var->mi_t*(var->duy + var->dvx);	    

        double P = var->dux*txx + var->dvy*tyy + (var->duy + var->dvx)*txy;        
        
        P = fmin(P, 10*sst->bs*var->r*var->om*var->k);
        //P = fmin(P, var->mi_t*omega*omega);        
                
        var->F1 = F1;
        var->F2 = F2;
                
        double cd = 2*(1-F1)*var->r*sst->so2*(var->dkx*var->dox + var->dky*var->doy)/var->om;
                
        var->Qtk = P - sst->bs*var->r*var->om*var->k;
        var->Qto = g*P/n_t - b*var->r*var->om*var->om + cd;
        
        var->dQodro = (- 2*b*var->r*var->om - cd/var->om)/var->r;
        
    }
    
}

void sstFlux(SST* sst, SSTVAR* var)
{

    if(sst->trans->flag)
    {
        sstTransFlux(sst->trans, var);
    }
    else
    {

        double n_L = var->mi_L/var->r;

        double omega = fabs(var->duy - var->dvx);

        double sqrtk_term = sqrt(var->k)/(sst->bs*var->om*var->d);
        
        double n_L_term = 500*n_L/(var->d*var->d*var->om);

        double F1 = sstF1(sst, var, n_L_term, sqrtk_term);

        double F2 = sstF2(sst, var, n_L_term, sqrtk_term);

        var->mi_t = var->r*sst->a1*var->k/fmax(sst->a1*var->om, omega*F2);
        //var->mi_t = fmin(var->mi_t, 1e6*var->mi_L);

        double sk = sstBlend(sst->sk1, sst->sk2, F1);
        double so = sstBlend(sst->so1, sst->so2, F1);

        double aux;

        aux = (var->mi_L + sk*var->mi_t);
        var->tkx = aux*var->dkx;
        var->tky = aux*var->dky;
        aux = (var->mi_L + so*var->mi_t);
        var->tox = aux*var->dox;
        var->toy = aux*var->doy;
    }
}

double sstF1(SST* sst, SSTVAR* var, double n_L_term, double sqrtk_term)
{
    long double aux1, aux2;

    aux1 = 2*var->r*sst->so2*(var->dkx*var->dox + var->dky*var->doy)/var->om;
    if(var->om == 0)
    {
        aux1 = 0;
    }
    double CD = fmax(aux1, 1.0e-20);

    aux1 = sqrtk_term;
    aux2 = n_L_term;
    
    aux1 = fmax(aux1, aux2);
    aux2 = 4*var->r*sst->so2*var->k/(CD*var->d*var->d);
    
    double arg1 = fmin(aux1, aux2);
    double F1 = tanh(arg1*arg1*arg1*arg1);

    return F1;
}

double sstF2(SST* sst, SSTVAR* var, double n_L_term, double sqrtk_term)
{
    double arg2 = fmax(2*sqrtk_term, n_L_term);

    double F2 = tanh(arg2*arg2);

    return F2;
}

void sstSolverWriteSurf(SOLVER* solver)
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
                                
                    sstBoundaryFaceViscFlux(solver, bc, ii, f, &mi);

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

void sstInterMiT(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {

        ELEMENT* E0 = solver->mesh->elemL[ii];

        SSTVAR var;

        var.dux = solver->dPx[1][ii];
        var.dvx = solver->dPx[2][ii];
        var.dkx = solver->dPx[4][ii];
        var.dox = solver->dPx[5][ii];

        var.duy = solver->dPy[1][ii];
        var.dvy = solver->dPy[2][ii];        
        var.dky = solver->dPy[4][ii];
        var.doy = solver->dPy[5][ii];        

        // Flow variables in the face
        var.r = E0->P[0];
        var.T = E0->P[4];
        var.k = E0->P[5];
        var.om = E0->P[6];
        var.d = E0->d;
        
        var.mi_L = sutherland(var.T);
        
        sstSources(solver->sst, &var);

        double x0, y0;

        elementCenter(E0, solver->mesh, &x0, &y0);

        solver->miTe[ii] = var.mi_t;
        solver->F1[ii] = var.F1;
        solver->F2[ii] = var.F2;
        solver->dd[ii] = var.d;
        solver->om2[ii] = var.om;
               
    }
}


