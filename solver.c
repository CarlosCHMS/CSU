#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"
#include"boundary.h"

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny)
{
  
    CONDITION* cond = malloc(sizeof(CONDITION));

    cond->p = p;
    cond->T = T;
    cond->mach = mach;
    cond->nx = nx;
    cond->ny = ny;
    
    return cond;

}

void conditionState(CONDITION* cond, SOLVER* solver)
{

    double r, u, v, E, c;
    
    r = cond->p/(solver->Rgas*cond->T);
    c = sqrt(solver->gamma*solver->Rgas*cond->T);
    u = cond->nx*cond->mach*c;
    v = cond->ny*cond->mach*c;
    E = (solver->Rgas*cond->T)/(solver->gamma-1) + (u*u + v*v)/2;
    
    cond->Uin[0] = r;
    cond->Uin[1] = r*u;
    cond->Uin[2] = r*v;
    cond->Uin[3] = r*E;   

    cond->Pin[0] = r;
    cond->Pin[1] = u;
    cond->Pin[2] = v;
    cond->Pin[3] = solver->Rgas*cond->T*r;
    cond->Pin[4] = cond->T;


}

double conditionVref(CONDITION* cond, SOLVER* solver)
{
      
    double c = sqrt(solver->gamma*solver->Rgas*cond->T);
    double Vref = (cond->mach+1)*c;
        
    return Vref;

}

void solverMalloc(SOLVER* solver)
{
    solver->U = tableMallocDouble(4, solver->mesh->Nelem);
    solver->Uaux = tableMallocDouble(4, solver->mesh->Nelem);    
    solver->R = tableMallocDouble(4, solver->mesh->Nelem);
    solver->faceFlux = tableMallocDouble(4, solver->mesh->Ncon);
    solver->dPx = tableMallocDouble(4, solver->mesh->Nelem);
    solver->dPy = tableMallocDouble(4, solver->mesh->Nelem);    
    solverMallocP(solver);
}

void solverFree(SOLVER* solver)
{

    tableFreeDouble(solver->U, 4);
    tableFreeDouble(solver->Uaux, 4);
    tableFreeDouble(solver->R, 4);        
    tableFreeDouble(solver->faceFlux, 4);
    tableFreeDouble(solver->dPx, 4);
    tableFreeDouble(solver->dPy, 4);
    solverFreeP(solver);
    meshFree(solver->mesh);
    
}

void solverWrite(SOLVER* solver, char* fileName)
{

    /* 
        Based on: Adek Tasri, Accuracy of Cell Centres to Vertices 
        Interpolation for Unstructured Mesh Finite Volume Solver, 2021
    */

    FILE* ff = fopen(fileName, "w");
    double** Up = tableMallocDouble(4, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    int ii, jj, kk, p;
    double xc, yc, xp, yp, L;
    ELEMENT* e;

    for(ii=0; ii<solver->mesh->Np; ii++)
    {       
        for(kk=0; kk<4; kk++)
        {
            Up[kk][ii] = 0.;
        }   
        den[ii] = 0.;
    }

    for(ii=0; ii<solver->mesh->Nelem; ii++)
    {
        meshElemCenter(solver->mesh, ii, &xc, &yc);
        e = solver->mesh->elemL[ii];
        
        for(jj=0; jj<e->Np; jj++)
        {
            p = e->p[jj];
            xp = solver->mesh->p[p][0];
            yp = solver->mesh->p[p][1];
            L = sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
           
            for(kk=0; kk<4; kk++)
            {        
                Up[kk][p] += solver->U[kk][ii]/L;
            }
            
            den[p] += 1/L;
        }        
        
    }


    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        if(den[ii] != 0)
        {
            for(kk=0; kk<4; kk++)
            {
                Up[kk][ii] /= den[ii];    
            }
        }
    }

    fprintf(ff, "0, %i, 4,\n", solver->mesh->Np);

    for(ii=0; ii<solver->mesh->Np; ii++)
    {        
        for(jj=0; jj<4; jj++)
        {
            fprintf(ff, "%.10e, ", Up[jj][ii]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);

    tableFreeDouble(Up, 4);
    free(den);

}

void solverInitU(SOLVER* solver, CONDITION* inside)
{

    conditionState(inside, solver);

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            solver->U[kk][ii] = inside->Uin[kk];
        }
    }
}

void solverResetR(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            solver->R[kk][ii] = 0.0; 
        }
    }
}

double solverCalcP(SOLVER* solver, double** U, int ii)
{

	double u = U[1][ii]/U[0][ii];
	double v = U[2][ii]/U[0][ii];
    
	return (solver->gamma - 1)*(U[3][ii] - 0.5*(u*u + v*v)*U[0][ii]);
    
}

void solverCalcVel(SOLVER* solver, double** U, int ii, double* u, double* v, double* c)
{
    
    double E = U[3][ii]/U[0][ii];
    double aux;

    *u = U[1][ii]/U[0][ii];
    *v = U[2][ii]/U[0][ii];
    
    aux = (E - ((*u)*(*u) + (*v)*(*v))/2);
    aux *= solver->gamma - 1;
    *c = sqrt(aux*solver->gamma);
    
}

void rotation(double* U, double dSx, double dSy, double dS)
{

    double aux = (U[1]*dSx + U[2]*dSy)/dS;
    U[2] = (-U[1]*dSy + U[2]*dSx)/dS;
    U[1] = aux;
	
}

void inter(SOLVER* solver)
{

    if(solver->order == 2)
	{
        # pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
	        int jj, kk, p0, p1;
            double dPx, dPy;
            double x0, y0, d2, phi0, phi, xm, ym;
            double Pmin, Pmax; 
            double blend = 1.0;        

	        if(solver->mesh->elemL[ii]->neiN > 1)
	        {  
          	    meshElemCenter(solver->mesh, ii, &x0, &y0);
          	    
                for(kk=0; kk<4; kk++)
                {
                    solverCalcGrad2(solver, solver->P[kk], ii, &dPx, &dPy, &Pmin, &Pmax);
                    
                    for(jj=0; jj<solver->mesh->elemL[ii]->Np; jj++)
                    {
                        p0 = solver->mesh->elemL[ii]->p[jj]; 
                        p1 = solver->mesh->elemL[ii]->p[(jj+1)%solver->mesh->elemL[ii]->Np];
                        xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
                        ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;        
                        d2 = (dPx*(xm - x0) + dPy*(ym - y0));
                        //phi0 = limiterBJ(solver->P[kk][ii], Pmin, Pmax, d2);                        
                        phi0 = limiterV(solver->P[kk][ii], Pmin, Pmax, d2, solver->e);
                        
                        if(jj==0)
                        {
                            phi = phi0;
                        }
                        else
                        {
                            phi = fmin(phi, phi0);
                        }			            
                    }

                    solver->dPx[kk][ii] = phi*dPx*blend;
	                solver->dPy[kk][ii] = phi*dPy*blend;
	            }   
	        }
        }
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
	    int kk;
        double dSx, dSy, dS;
        double x0, y0, x1, y1, xm, ym;
        double PL[4];
	    double PR[4];
        double f[4];
 
        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        if(solver->order == 2)
        {
           	xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
            ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

		    meshElemCenter(solver->mesh, e0, &x0, &y0);
		    meshElemCenter(solver->mesh, e1, &x1, &y1);
		    
            for(kk=0; kk<4; kk++)
		    {                
		        if(solver->mesh->elemL[e0]->neiN > 1)
		        {		    
		            PL[kk] = solver->P[kk][e0] + solver->dPx[kk][e0]*(xm - x0) + solver->dPy[kk][e0]*(ym - y0);
		        }
		        else
		        {
		            PL[kk] = solver->P[kk][e0];
		        }

		        if(solver->mesh->elemL[e1]->neiN > 1)
		        {      	        
		            PR[kk] = solver->P[kk][e1] + solver->dPx[kk][e1]*(xm - x1) + solver->dPy[kk][e1]*(ym - y1);
		        }
		        else
		        {
		            PR[kk] = solver->P[kk][e1];
		        }
            }   
        }
        else
        {       
            for(kk=0; kk<4; kk++)
		    {
			    PL[kk] = solver->P[kk][e0];
			    PR[kk] = solver->P[kk][e1];			    
		    }
		}
		
        // Rotation of the velocity vectors
		rotation(PL, dSx, dSy, dS);
                    
        // Rotation of the velocity vectors
		rotation(PR, dSx, dSy, dS);
        
        // Flux calculation
        flux(solver, PL[0], PL[1], PL[2], PL[3], PR[0], PR[1], PR[2], PR[3], f);

        // Rotation of the flux
		rotation(f, dSx, -dSy, dS);
        
        for(kk=0; kk<4; kk++)
        {
            solver->faceFlux[kk][ii] = f[kk]*dS;
        }
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->neiN; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=0; kk<4; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][face-1];
                }
            }
            else
            {
                for(int kk=0; kk<4; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][-face-1];
                }            
            }
        }
    } 
}

void interVisc(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        int kk;
        double dPx, dPy, Pmin, Pmax;

        if(solver->mesh->elemL[ii]->neiN > 1)
        {   
      	    kk = 1;
            solverCalcGrad2(solver, solver->P[kk], ii, &dPx, &dPy, &Pmin, &Pmax);
            solver->dPx[kk][ii] = dPx;
            solver->dPy[kk][ii] = dPy;
            
            kk = 2;
            solverCalcGrad2(solver, solver->P[kk], ii, &dPx, &dPy, &Pmin, &Pmax);
            solver->dPx[kk][ii] = dPx;
            solver->dPy[kk][ii] = dPy;
            
            // grad p storage is being reused for T grad
            kk = 3;
            solverCalcGrad2(solver, solver->P[4], ii, &dPx, &dPy, &Pmin, &Pmax);
            solver->dPx[kk][ii] = dPx;
            solver->dPy[kk][ii] = dPy;
            
        }
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
        double dSx, dSy;
        double x0, y0, x1, y1;

        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];

        if((solver->mesh->elemL[e0]->neiN > 1) & (solver->mesh->elemL[e1]->neiN > 1))
		{		    
            meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
                        
            double duxm = (solver->dPx[1][e0] + solver->dPx[1][e1])*0.5;
            double dvxm = (solver->dPx[2][e0] + solver->dPx[2][e1])*0.5;        
            double dTxm = (solver->dPx[3][e0] + solver->dPx[3][e1])*0.5;
            
            double duym = (solver->dPy[1][e0] + solver->dPy[1][e1])*0.5;
            double dvym = (solver->dPy[2][e0] + solver->dPy[2][e1])*0.5;
            double dTym = (solver->dPy[3][e0] + solver->dPy[3][e1])*0.5;
            
	        meshElemCenter(solver->mesh, e0, &x0, &y0);
	        meshElemCenter(solver->mesh, e1, &x1, &y1);
		        
            double dx = x1 - x0;		    
            double dy = y1 - y0;
            double L = sqrt(dx*dx + dy*dy);
            
            double dul = (solver->P[1][e1] - solver->P[1][e0])/L;
            double dvl = (solver->P[2][e1] - solver->P[2][e0])/L;        
            double dTl = (solver->P[4][e1] - solver->P[4][e0])/L;
            
            double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

            double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;        

            double dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
            double dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;        
		    
            double T = (solver->P[4][e0] + solver->P[4][e1])*0.5;            
            double mi = sutherland(T);
            double k = solver->Cp*mi/solver->Pr;
		    		        
		    double txx = 2*mi*(dux - (dux + dvy)/3);
		    double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
		    double txy = mi*(duy + dvx);  
		    
		    double u = (solver->P[1][e0] + solver->P[1][e1])*0.5;
            double v = (solver->P[2][e0] + solver->P[2][e1])*0.5;
	    
		    solver->faceFlux[1][ii] = txx*dSx + txy*dSy;
		    solver->faceFlux[2][ii] = txy*dSx + tyy*dSy;
		    solver->faceFlux[3][ii] = (txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy);
		}
		else
		{
		    solver->faceFlux[1][ii] = 0.0;
		    solver->faceFlux[2][ii] = 0.0;
		    solver->faceFlux[3][ii] = 0.0;
		}
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->neiN; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=1; kk<4; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][face-1];
                }
            }
            else
            {
                for(int kk=1; kk<4; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][-face-1];
                }            
            }
        }
    } 
}

void interAxisPressure(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double dS;        
        dS = meshCalcDSlateral(solver->mesh, ii);
        solver->R[2][ii] -= solver->P[3][ii]*dS;
    }
}

void solverCalcR(SOLVER* solver, double** U)
{

    solverResetR(solver);

    solverCalcPrimitive(solver, U);

    inter(solver);
        
    boundary(solver); 
    
    if(solver->mesh->axi==1)
    {
        interAxisPressure(solver);
    }    
    
    if(solver->laminar==1)
    {
        interVisc(solver);
        boundaryVisc(solver);
    }
}

void solverRK(SOLVER* solver, double a)
{   
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double omega = meshCalcOmega(solver->mesh, ii);
        for(int kk=0; kk<4; kk++)
        {
            
            solver->Uaux[kk][ii] = solver->U[kk][ii] - solver->dt*a*solver->R[kk][ii]/omega;
            
        }
    }
}

void solverUpdateU(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            solver->U[kk][ii] = solver->Uaux[kk][ii];
        }
    }
}

void solverStepRK(SOLVER* solver)
{   

    if(solver->stages==3)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.1481);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==4)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.0833);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.2069);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4265);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==5)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.0533);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.1263);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.2375);

        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4414);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    
    solverUpdateU(solver);

}

void solverCalcRes(SOLVER* solver)
{
    
    
    for(int kk=0; kk<4; kk++)
    {
        solver->res[kk] = fabs(solver->R[kk][0]);
    }
    
    for(int ii=1; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            if(solver->res[kk] < fabs(solver->R[kk][ii]))
            {
                solver->res[kk] = fabs(solver->R[kk][ii]);
            }
        }
    }
    
    for(int kk=0; kk<4; kk++)
    {
        printf(" %+.4e,", solver->res[kk]);        
    }
    printf("\n");    
   
}

double solverLocalTimeStep(SOLVER* solver, int ii)
{
    int p0, p1;
    double dSx, dSy, dSxm, dSym, Lx, Ly, u, v, c;
    ELEMENT* e = solver->mesh->elemL[ii];
    
    dSxm = 0.;
    dSym = 0.;
    
    for(int jj=0; jj<e->Np; jj++)
    {
    
        p0 = e->p[jj];
        p1 = e->p[(jj+1)%e->Np];
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dSxm += fabs(dSx);
        dSym += fabs(dSy);
        
    }
        
    dSxm *= 0.5;
    dSym *= 0.5;
    
    solverCalcVel(solver, solver->U, ii, &u, &v, &c);
    
    Lx = (fabs(u) + c)*dSxm;
    Ly = (fabs(v) + c)*dSym;
    
    return meshCalcOmega(solver->mesh, ii)/(Lx + Ly);
    
}

double solverCalcDt(SOLVER* solver)
{

    double dtLocal;
    double dt = solverLocalTimeStep(solver, 0);

    for(int ii=1; ii<solver->mesh->Nelem; ii++)
    {
        dtLocal = solverLocalTimeStep(solver, ii);
        if(dtLocal < dt)
        {
            dt = dtLocal;
        }
    }

    dt *= solver->CFL*0.5*solver->stages;

    return dt;
}

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm)
{

    double x, y;

    conditionState(inside1, solver);
    conditionState(inside2, solver);

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            meshElemCenter(solver->mesh, ii, &x, &y);
            if(x<xm)
            {
                solver->U[kk][ii] = inside1->Uin[kk];
            }
            else
            {
                solver->U[kk][ii] = inside2->Uin[kk];
            }
        }
    }
}

void solverCalcGrad2(SOLVER* solver, double* U, int ii, double* dUx, double* dUy, double* Umin, double* Umax)
{

    /*
    Least squares calculation
    */

    int e;
    double xx, yy;
    double x[5];
    double y[5];    
    double u[5];    
    double a, b, c, A, B;
    int nN = solver->mesh->elemL[ii]->neiN;
        
    meshElemCenter(solver->mesh, ii, &xx, &yy);
    x[0] = xx;
    y[0] = yy;
    u[0] = U[ii];

    for(int jj=0; jj<nN; jj++)
    {
        e = solver->mesh->elemL[ii]->neiL[jj]->ii;
        meshElemCenter(solver->mesh, e, &xx, &yy);
        x[jj+1] = xx;
        y[jj+1] = yy;
        u[jj+1] = U[e];
    }
    
    a = 0;    
    b = 0;
    c = 0;
    A = 0;
    B = 0;
    for(int jj=1; jj<nN+1; jj++)
    {
        a += (x[jj] - x[0])*(x[jj] - x[0]);
        b += (x[jj] - x[0])*(y[jj] - y[0]);
        c += (y[jj] - y[0])*(y[jj] - y[0]);
        A += (x[jj] - x[0])*(u[jj] - u[0]);
        B += (y[jj] - y[0])*(u[jj] - u[0]);
    }
        
    *dUx = (c*A - b*B)/(a*c - b*b);
    *dUy = (-b*A + a*B)/(a*c - b*b);
    
    *Umin = u[0];
    *Umax = u[0];
    
    for(int ii=1; ii<nN+1; ii++)
    {
        *Umax = fmax(*Umax, u[ii]);
        *Umin = fmin(*Umin, u[ii]);        
    }
}

void solverCheckGrad(SOLVER* solver)
{

    double *U = malloc(solver->mesh->Nelem*sizeof(double));
    double *xx = malloc(solver->mesh->Nelem*sizeof(double));
    double x, y;
    double dUx, dUy;
    double Umin, Umax;
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        meshElemCenter(solver->mesh, ii, &x, &y);
        xx[ii] = x;
        U[ii] = x*x;
    }
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        if(solver->mesh->elemL[ii]->neiN == 2)
        {
            solverCalcGrad2(solver, U, ii, &dUx, &dUy, &Umin, &Umax);
            printf("%+e, %+e, %+e\n", dUx, 2*xx[ii], dUx - 2*xx[ii]);
            //printf("%+e, %+e\n", Umin, Umax);            
        }
    }

    free(U);
    free(xx);

}

double limiterBJ(double Ui, double Umin, double Umax, double d2)
{

    double ans;
    if(d2 == 0)
    {
        ans = 1;
    }
    else if(d2 > 0)
    {
        ans = fmin(1, (Umax - Ui)/d2);
    }
    else
    {
        ans = fmin(1, (Umin - Ui)/d2);
    }
    
    return ans;

}

double limiterV(double Ui, double Umin, double Umax, double d2, double e)
{

    double ans;
    double d1max = Umax - Ui;
    double d1min = Umin - Ui;
    
    if(d2 == 0)
    {
        ans = 1;
    }
    else if(d2 > 0)
    {
        ans = (d1max*d1max + e*e)*d2 + 2*d2*d2*d1max;
        ans /= d1max*d1max + 2*d2*d2 + d1max*d2 + e*e;
        ans /= d2;
    }
    else
    {
        ans = (d1min*d1min + e*e)*d2 + 2*d2*d2*d1min;
        ans /= d1min*d1min + 2*d2*d2 + d1min*d2 + e*e;
        ans /= d2;
    }
    
    return ans;

}

void solverCalcPrimitive(SOLVER* solver, double** U)
{   
    solver->P[0] = U[0]; // Reduce memory usage

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {        
        solver->P[1][ii] = U[1][ii]/U[0][ii];
        solver->P[2][ii] = U[2][ii]/U[0][ii];
        solver->P[3][ii] = solverCalcP(solver, U, ii);
        solver->P[4][ii] = solver->P[3][ii]/(solver->P[0][ii]*solver->Rgas);
    }
}

void solverCalcUfromP(SOLVER* solver, double r, double u, double v, double p, double* U0, double* U1, double* U2, double* U3)
{

    *U0 = r;
    *U1 = r*u;
    *U2 = r*v;
    *U3 = p/(solver->gamma-1) + 0.5*(u*u + v*v)*r;

}

void solverMallocP(SOLVER* solver)
{

    solver->P = malloc(5*sizeof(double*));

    for(int ii=1; ii<5; ii++)
    {    
        solver->P[ii] = malloc(solver->mesh->Nelem*sizeof(double));
    }

}    
 
void solverFreeP(SOLVER* solver)
{

    for(int ii=1; ii<5; ii++)
    {    
        free(solver->P[ii]);
    }    

    free(solver->P);

}

double sutherland(double T)
{
    return 1.45e-6*T*sqrt(T)/(T + 110.0);
}
