#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<stdbool.h>
#include<omp.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"gasprop.h"
#include"laminar.h"

void laminarInter(SOLVER* solver)
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
        
        double duym = (solver->dPy[1][e0] + solver->dPy[1][e1])*0.5;
        double dvym = (solver->dPy[2][e0] + solver->dPy[2][e1])*0.5;
        double dTym = (solver->dPy[3][e0] + solver->dPy[3][e1])*0.5;
        
	    elementCenter(E0, solver->mesh, &x0, &y0);
	    elementCenter(E1, solver->mesh, &x1, &y1);
	        
        double dx = x1 - x0;		    
        double dy = y1 - y0;
        double L = sqrt(dx*dx + dy*dy);
        
        double dul = (E1->P[1] - E0->P[1])/L;
        double dvl = (E1->P[2] - E0->P[2])/L;
        double dTl = (E1->P[4] - E0->P[4])/L;
        
        double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

        double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;        

        double dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        double dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;        
	    
        double T = (E1->P[4] + E0->P[4])*0.5;
        double mi = sutherland(T);
        double k = gasprop_T2Cp(solver->gas, T)*mi/solver->Pr;
	    		        
	    double txx = 2*mi*(dux - (dux + dvy)/3);
	    double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	    double txy = mi*(duy + dvx);  
	    
	    double u = (E1->P[1] + E0->P[1])*0.5;
        double v = (E1->P[2] + E0->P[2])*0.5;
    
	    solver->faceFlux[1][ii] = txx*dSx + txy*dSy;
	    solver->faceFlux[2][ii] = txy*dSx + tyy*dSy;
	    solver->faceFlux[3][ii] = (txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy);
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->Np; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=1; kk<4; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][face-1];
                }
            }
            else if(face < 0)
            {
                for(int kk=1; kk<4; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][-face-1];
                }            
            }
        }
    } 
}

void laminarBoundaryVisc(SOLVER* solver, MESHBC* bc)
{

    int e0;
    double f[4];

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        e0 = bc->elemL[ii]->neiL[0]->ii;
 
        boundaryFaceViscFlux(solver, bc, ii, f);
 
        solver->R[1][e0] -= f[1];
        solver->R[2][e0] -= f[2];
        solver->R[3][e0] -= f[3];
    }        
}

void laminarBoundary(SOLVER* solver)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        laminarBoundaryVisc(solver, solver->mesh->bc[ii]);
    }
}


void boundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f)
{
    int e0, p0, p1;
    double x0, x1, y0, y1;
    double dux, duy, dvx, dvy, dTx, dTy;

    p0 = bc->elemL[ii]->p[0];
    p1 = bc->elemL[ii]->p[1];

    ELEMENT* E0 = bc->elemL[ii]->neiL[0];
    e0 = bc->elemL[ii]->neiL[0]->ii;

    if(bc->flagBC == 0)
    {

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
        
        double T = E0->P[4];
        double mi = sutherland(T);
        double k = gasprop_T2Cp(solver->gas, T)*mi/solver->Pr;            
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx);  
        
        double u = E0->P[1];
        double v = E0->P[2];		        
        
        f[1] = (txx*nx + txy*ny)*dS;
        f[2] = (txy*nx + tyy*ny)*dS;
        f[3] = (u*(txx*nx + txy*ny) + v*(txy*nx + tyy*ny) + k*(dTx*nx + dTy*ny))*dS;
        
    }
    else if(bc->flagBC == 1)
    {
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

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];
        
        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];

        dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

        dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;  
        
        dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;
        
        double T = E0->P[4];
        double mi = sutherland(T);
        double k = gasprop_T2Cp(solver->gas, T)*mi/solver->Pr;            
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx);  
        
        double u = E0->P[1];
        double v = E0->P[2];		        
        
        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy);
        
    }
    else if(bc->flagBC == 3)
    {
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

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        
        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];

        dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

        dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;
                        
        double T = E0->P[4];
        double mi = sutherland(T);
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx);  
        
        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        
    }
    else if(bc->flagBC == 4)
    {
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

        double duxm = solver->dPx[1][e0];
        double dvxm = solver->dPx[2][e0];
        double dTxm = solver->dPx[3][e0];
        
        double duym = solver->dPy[1][e0];
        double dvym = solver->dPy[2][e0];
        double dTym = solver->dPy[3][e0];

        dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

        dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;  
        
        dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;
                        
        double T = solver->Twall;
        double mi = sutherland(T);
        double k = gasprop_T2Cp(solver->gas, T)*mi/solver->Pr;            
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx);  
        
        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = k*(dTx*dSx + dTy*dSy);
        
    }        
    else
    {
        double dSx, dSy;
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);        
    
        dux = solver->dPx[1][e0];
        dvx = solver->dPx[2][e0];
        dTx = solver->dPx[3][e0];
        
        duy = solver->dPy[1][e0];
        dvy = solver->dPy[2][e0];
        dTy = solver->dPy[3][e0];                

        double T = E0->P[4];
        double mi = sutherland(T);
        double k = gasprop_T2Cp(solver->gas, T)*mi/solver->Pr;            
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx);  
        
        double u = E0->P[1];
        double v = E0->P[2];		        
        
        f[1] = txx*dSx + txy*dSy;
        f[2] = txy*dSx + tyy*dSy;
        f[3] = u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy);            
    }
}



void laminarWriteSurf(SOLVER* solver)
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
            fprintf(ff, "x,y,rho,u,v,p,T,n,Cp,Cfx,Cfy,q,yplus,mach,\n");
                    
            int Nvar = 12;
            int Nelem = bc->Nelem;
            double** D = tableMallocDouble(Nvar, Nelem);
            double** Dp = tableMallocDouble(Nvar, Nelem+1);
            int* d = malloc((Nelem+1)*sizeof(int));
            int* pn = malloc((Nelem+1)*sizeof(int));
            int* pn2 = malloc((solver->mesh->Np)*sizeof(int));
            
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
        
            double f[4];
        
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
                    
                    
                    boundaryFaceViscFlux(solver, bc, ii, f);

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


