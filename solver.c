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
#include"readTables.h"
#include"sa.h"

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

    double r, u, v, E, c, n;
    
    r = cond->p/(solver->Rgas*cond->T);
    c = sqrt(solver->gamma*solver->Rgas*cond->T);
    u = cond->nx*cond->mach*c;
    v = cond->ny*cond->mach*c;
    E = (solver->Rgas*cond->T)/(solver->gamma-1) + (u*u + v*v)/2;
    n = solver->turbRatio*sutherland(cond->T)/r;

    cond->Uin[0] = r;
    cond->Uin[1] = r*u;
    cond->Uin[2] = r*v;
    cond->Uin[3] = r*E;   

    cond->Pin[0] = r;
    cond->Pin[1] = u;
    cond->Pin[2] = v;
    cond->Pin[3] = solver->Rgas*cond->T*r;
    cond->Pin[4] = cond->T;

    if(solver->sa == 1)
    {
        cond->Uin[4] = r*n;
        cond->Pin[5] = n;         
    }

}

double conditionVref(CONDITION* cond, SOLVER* solver)
{
      
    double c = sqrt(solver->gamma*solver->Rgas*cond->T);
    double Vref = (cond->mach+1)*c;
        
    return Vref;

}

void solverMalloc(SOLVER* solver)
{
    solver->U = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->Uaux = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);    
    solver->R = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->faceFlux = tableMallocDouble(solver->Nvar, solver->mesh->Ncon);
    solver->dPx = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->dPy = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);    
    
    if(solver->dtLocal == 1)
    {
        solver->dtL = malloc(solver->mesh->Nelem*sizeof(double));
    }
    
}

void solverFree(SOLVER* solver)
{

    tableFreeDouble(solver->U, solver->Nvar);
    tableFreeDouble(solver->Uaux, solver->Nvar);
    tableFreeDouble(solver->R, solver->Nvar);        
    tableFreeDouble(solver->faceFlux, solver->Nvar);
    tableFreeDouble(solver->dPx, solver->Nvar);
    tableFreeDouble(solver->dPy, solver->Nvar);
    meshFree(solver->mesh);
    
    if(solver->dtLocal == 1)
    {
        free(solver->dtL);
    }
    
    inputFree(solver->input);
    
}

void solverWriteSolution(SOLVER* solver)
{

    /* 
        Based on: Adek Tasri, Accuracy of Cell Centres to Vertices 
        Interpolation for Unstructured Mesh Finite Volume Solver, 2021
    */

    // Save solution
    printf("main: saving the solution.\n");    

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "solution.csv");
    
    FILE* ff = fopen(fileName, "w");
    double** Up = tableMallocDouble(solver->Nvar, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    int ii, jj, kk, p;
    double xc, yc, xp, yp, L;
    ELEMENT* E;

    for(ii=0; ii<solver->mesh->Np; ii++)
    {       
        for(kk=0; kk<solver->Nvar; kk++)
        {
            Up[kk][ii] = 0.;
        }   
        den[ii] = 0.;
    }

    for(ii=0; ii<solver->mesh->Nelem; ii++)
    {
        E = solver->mesh->elemL[ii];
        elementCenter(E, solver->mesh, &xc, &yc);
        
        for(jj=0; jj<E->Np; jj++)
        {
            p = E->p[jj];
            xp = solver->mesh->p[p][0];
            yp = solver->mesh->p[p][1];
            L = sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
           
            for(kk=0; kk<solver->Nvar; kk++)
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
            for(kk=0; kk<solver->Nvar; kk++)
            {
                Up[kk][ii] /= den[ii];    
            }
        }
    }

    fprintf(ff, "0, %i, %i,\n", solver->mesh->Np, solver->Nvar);

    for(ii=0; ii<solver->mesh->Np; ii++)
    {        
        for(jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", Up[jj][ii]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);

    tableFreeDouble(Up, solver->Nvar);
    free(den);

}

void solverWriteReestart(SOLVER* solver)
{

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "restart.csv");

    FILE* ff = fopen(fileName, "w");

    // Save reestart
    printf("main: saving the restart file.\n");    

    fprintf(ff, "1,\n");

    fprintf(ff, "0, %i, %i,\n", solver->mesh->Nelem, solver->Nvar);

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", solver->U[jj][ii]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);

}


void solverLoadRestart(SOLVER* solver, char* fileName)
{

    TABLELIST* tl = fReadTables(fileName);
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            solver->U[jj][ii] = tl->tables[0]->values[ii][jj];
        }        
    }    
    
    readTablesFree(tl);

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
        for(int kk=0; kk<solver->Nvar; kk++)
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
            ELEMENT* E = solver->mesh->elemL[ii];

      	    elementCenter(E, solver->mesh, &x0, &y0);
      	    
            for(kk=0; kk<solver->Nvar; kk++)
            {
                if(kk == 4)
                {
                    solverCalcGrad2(solver, E, 5, &dPx, &dPy, &Pmin, &Pmax);
                }
                else
                {
                    solverCalcGrad2(solver, E, kk, &dPx, &dPy, &Pmin, &Pmax);                
                }
                
                for(jj=0; jj<E->Np; jj++)
                {
                    p0 = E->p[jj]; 
                    p1 = E->p[(jj+1)%E->Np];
                    xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
                    ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;        
                    d2 = (dPx*(xm - x0) + dPy*(ym - y0));
                    //phi0 = limiterBJ(solver->mesh->elemL[ii]->P[kk], Pmin, Pmax, d2);                        
                    phi0 = limiterV(E->P[kk], Pmin, Pmax, d2, solver->e);
                    
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


    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
	    int kk;
        double dSx, dSy, dS;
        double x0, y0, x1, y1, xm, ym;
        double PL[5];
	    double PR[5];
        double f[5];
 
        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];
        
        ELEMENT* E0 = solver->mesh->elemL[e0];
        ELEMENT* E1 = solver->mesh->elemL[e1];        
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        if(solver->order == 2)
        {
           	xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
            ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

		    elementCenter(E0, solver->mesh, &x0, &y0);
		    elementCenter(E1, solver->mesh, &x1, &y1);
		    
            for(kk=0; kk<4; kk++)
		    {                
	            PL[kk] = E0->P[kk] + solver->dPx[kk][e0]*(xm - x0) + solver->dPy[kk][e0]*(ym - y0);
	            PR[kk] = E1->P[kk] + solver->dPx[kk][e1]*(xm - x1) + solver->dPy[kk][e1]*(ym - y1);
            }
            
            if(solver->Nvar == 5)
            {
                kk = 4;
                PL[kk] = E0->P[5] + solver->dPx[kk][e0]*(xm - x0) + solver->dPy[kk][e0]*(ym - y0);
	            PR[kk] = E1->P[5] + solver->dPx[kk][e1]*(xm - x1) + solver->dPy[kk][e1]*(ym - y1);
            }  
        }
        else
        {       
            for(kk=0; kk<4; kk++)
		    {
			    PL[kk] = E0->P[kk];
			    PR[kk] = E1->P[kk];
		    }
		    
		    if(solver->Nvar == 5)
            {
	            PL[4] = E0->P[5];
	            PR[4] = E1->P[5];
            }
		}

        // Rotation of the velocity vectors
		rotation(PL, dSx, dSy, dS);
                    
        // Rotation of the velocity vectors
		rotation(PR, dSx, dSy, dS);
        
        // Flux calculation
        if(solver->sa==1)
        {
            fluxAUSMDV_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4], PR[0], PR[1], PR[2], PR[3], PR[4], f);
        }
        else
        {
            flux(solver, PL[0], PL[1], PL[2], PL[3], PR[0], PR[1], PR[2], PR[3], f);
        }
        // Rotation of the flux
		rotation(f, dSx, -dSy, dS);
        
        for(kk=0; kk<solver->Nvar; kk++)
        {
            solver->faceFlux[kk][ii] = f[kk]*dS;
        }
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->elemL[ii]->Np; jj++)
        {
            int face = solver->mesh->elemL[ii]->f[jj];
            if(face > 0)
            {
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][face-1];
                }
            }
            else if(face < 0)
            {
                for(int kk=0; kk<solver->Nvar; kk++)
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
        ELEMENT* E = solver->mesh->elemL[ii];

  	    kk = 1;
        solverCalcGrad2(solver, E, kk, &dPx, &dPy, &Pmin, &Pmax);
        solver->dPx[kk][ii] = dPx;
        solver->dPy[kk][ii] = dPy;
        
        kk = 2;
        solverCalcGrad2(solver, E, kk, &dPx, &dPy, &Pmin, &Pmax);
        solver->dPx[kk][ii] = dPx;
        solver->dPy[kk][ii] = dPy;
        
        // grad p storage is being reused for T grad
        kk = 3;
        solverCalcGrad2(solver, E, 4, &dPx, &dPy, &Pmin, &Pmax);
        solver->dPx[kk][ii] = dPx;
        solver->dPy[kk][ii] = dPy;            
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
        double k = solver->Cp*mi/solver->Pr;
	    		        
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

void interAxisPressure(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double dS;        
        dS = meshCalcDSlateral(solver->mesh, ii);
        solver->R[2][ii] -= solver->mesh->elemL[ii]->P[3]*dS;
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
    
    if(solver->sa==1)
    {
        saInter(solver);
        saBoundary(solver);
    }
    
}

void solverRK(SOLVER* solver, double a)
{   
    if(solver->dtLocal == 1)
    {
        # pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {            
                solver->Uaux[kk][ii] = solver->U[kk][ii] - solver->dtL[ii]*a*solver->R[kk][ii]/solver->mesh->elemL[ii]->omega;
            }
        }    
    }
    else
    {
        # pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {            
                solver->Uaux[kk][ii] = solver->U[kk][ii] - solver->dt*a*solver->R[kk][ii]/solver->mesh->elemL[ii]->omega;
            }
        }
    }
}

void solverUpdateU(SOLVER* solver)
{
    double** aux;
    
    aux = solver->Uaux;
    solver->Uaux = solver->U;
    solver->U = aux;

    /*
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->U[kk][ii] = solver->Uaux[kk][ii];
        }
    }
    */
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
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        solver->res[kk] = fabs(solver->R[kk][0]);
    }
    
    for(int ii=1; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            if(solver->res[kk] < fabs(solver->R[kk][ii]))
            {
                solver->res[kk] = fabs(solver->R[kk][ii]);
            }
        }
    }
    
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        printf(" %+.4e,", solver->res[kk]);        
    }
    printf("\n");      
}

double solverLocalTimeStep(SOLVER* solver, int ii)
{
    int p0, p1;
    double dSx, dSy, dSxm, dSym, Lx, Ly, u, v, c;
    ELEMENT* E = solver->mesh->elemL[ii];
    
    dSxm = 0.;
    dSym = 0.;
    
    for(int jj=0; jj<E->Np; jj++)
    {
        p0 = E->p[jj];
        p1 = E->p[(jj+1)%E->Np];
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dSxm += fabs(dSx);
        dSym += fabs(dSy);   
    }
        
    dSxm *= 0.5;
    dSym *= 0.5;
    
    solverCalcVel(solver, solver->U, ii, &u, &v, &c);
    
    Lx = (fabs(u) + c)*dSxm;
    Ly = (fabs(v) + c)*dSym;
    
    return solver->mesh->elemL[ii]->omega/(Lx + Ly);
    
}

void solverCalcDt(SOLVER* solver)
{
    if(solver->dtLocal == 1)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            solver->dtL[ii] = (solver->CFL*0.5*solver->stages)*solverLocalTimeStep(solver, ii);
        }
    }
    else
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
        solver->dt = dt;
    }
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
    	    elementCenter(solver->mesh->elemL[ii], solver->mesh, &x, &y);
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

void solverCalcGrad2(SOLVER* solver, ELEMENT* E, int kk, double* dUx, double* dUy, double* Umin, double* Umax)
{
    /*
    Least squares calculation
    */

    double xx, yy;
    double x[5];
    double y[5];    
    double u[5];    
    double a, b, c, A, B;

    int nN = E->Np;
        
    elementCenter(E, solver->mesh, &xx, &yy);
    x[0] = xx;
    y[0] = yy;
    u[0] = E->P[kk];

    for(int jj=0; jj<nN; jj++)
    {
        elementCenter(E->neiL[jj], solver->mesh, &xx, &yy);
        x[jj+1] = xx;
        y[jj+1] = yy;
        u[jj+1] = E->neiL[jj]->P[kk];
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
    double *xx = malloc(solver->mesh->Nelem*sizeof(double));
    double x, y;
    double dUx, dUy;
    double Umin, Umax;
    ELEMENT* E;
    double aux;
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        E = solver->mesh->elemL[ii];
        elementCenter(E, solver->mesh, &x, &y);
        xx[ii] = x;
        E->P[0] = x*x;
    }

    for(int kk=0; kk<solver->mesh->Nmark; kk++)
    {
        for(int ii=0; ii<solver->mesh->bc[kk]->Nelem; ii++)
        {
            E = solver->mesh->bc[kk]->elemL[ii];
            elementCenter(E, solver->mesh, &x, &y);
            E->P[0] = x*x;
        }    
    }
  
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        E = solver->mesh->elemL[ii];    
        solverCalcGrad2(solver, E, 0, &dUx, &dUy, &Umin, &Umax);
        if(E->neiN==1)
        {
            //printf("%+e, %+e, %+e\n", dUx, 2*xx[ii], dUx/(2*xx[ii])-1);
            aux = xx[ii]*xx[ii];
            printf("%+e, %+e\n", aux - Umin, Umax - aux);  
        }
        
    }
    
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
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {         
        ELEMENT* E = solver->mesh->elemL[ii];
        
        E->P[0] = U[0][ii];
        E->P[1] = U[1][ii]/U[0][ii];
        E->P[2] = U[2][ii]/U[0][ii];
        E->P[3] = solverCalcP(solver, U, ii);
        E->P[4] = E->P[3]/(E->P[0]*solver->Rgas);

        if(solver->sa == 1)
        {
            E->P[5] = U[4][ii]/U[0][ii];
        }                
    }
    
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        boundaryCalcPrimitive(solver, solver->mesh->bc[ii]);
    }
}

void solverCalcUfromP(SOLVER* solver, double r, double u, double v, double p, double* U0, double* U1, double* U2, double* U3)
{

    *U0 = r;
    *U1 = r*u;
    *U2 = r*v;
    *U3 = p/(solver->gamma-1) + 0.5*(u*u + v*v)*r;

}

double sutherland(double T)
{
    return 1.45e-6*T*sqrt(T)/(T + 110.0);
}

void solverPrintP(SOLVER* solver)
{

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar+1; kk++)
        {
            printf("%e,", solver->mesh->elemL[ii]->P[kk]);
        }
        printf("\n");
    }

    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {
        for(int ii=0; ii<solver->mesh->bc[jj]->Nelem; ii++)
        {
            for(int kk=0; kk<solver->Nvar+1; kk++)
            {
                printf("%e,", solver->mesh->bc[jj]->elemL[ii]->P[kk]);
            }
            printf("\n");
        }
    }

}

void solverCalcCoeff(SOLVER* solver, double *Cx, double *Cy)
{
    int p0, p1;
    MESHBC* bc;
    double cp, dSx, dSy, fx, fy;

    *Cx = 0.0;
    *Cy = 0.0;
    
    double r = solver->inlet->Pin[0];
    double u = solver->inlet->Pin[1];
    double v = solver->inlet->Pin[2];
    double P = solver->inlet->Pin[3];
    double q = 0.5*r*(u*u + v*v);   
        
    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {
        bc = solver->mesh->bc[jj];
        if(bc->flagBC == 3)
        {
            for(int ii=0; ii<bc->Nelem; ii++)
            {
                cp = (bc->elemL[ii]->P[3] - P)/q;
                p0 = bc->elemL[ii]->p[0];
                p1 = bc->elemL[ii]->p[1];
                
                meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);                
             
                *Cx += cp*dSx;   
                *Cy += cp*dSy;                
                
                if(solver->laminar==1 || solver->sa==1)
                {
                    boundaryCalcFrictionWall(solver, bc->elemL[ii], &fx, &fy);
                    *Cx -= fx/q;
                    *Cy -= fy/q;
                }                               
            }
        }
    }
        
    *Cx /= solver->Sref;
    *Cy /= solver->Sref;
    
    if(solver->mesh->axi)
    {
        *Cx *= 2*M_PI;
    }    
}

void solverCalcCoeff2(SOLVER* solver, char* path)
{
    MESHBC* bc;
    double x, txx, txy, tyy, yp;
    char s[50];
        
    s[0] = '\0';
    strcat(s, path);
    strcat(s, "friction.csv");        
    FILE* ff = fopen(s, "w");   
        
    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {
        bc = solver->mesh->bc[jj];
        if(bc->flagBC == 3)
        {
            for(int ii=0; ii<bc->Nelem; ii++)
            {                
                if(solver->laminar==1 || solver->sa==1)
                {
                    
                    if(solver->sa==1)
                    {
                        saCalcTensorWall(solver, bc->elemL[ii], &txx, &txy, &tyy, &x, &yp);
                    }
                    else
                    {
                        boundaryCalcTensorWall(solver, bc->elemL[ii], &txx, &txy, &tyy, &x, &yp);
                    }
                    
                    fprintf(ff, "%e, %e, %e, %e, %e,\n", x, txx, txy, tyy, yp);
                }                               
            }
        }
    }
        
    fclose(ff);
}

void solverCalcCoeff3(SOLVER* solver, FILE* convFile, int Nint)
{
    int p0, p1;
    MESHBC* bc;
    double cp, dSx, dSy, fx, fy;    
    
    double r = solver->inlet->Pin[0];
    double u = solver->inlet->Pin[1];
    double v = solver->inlet->Pin[2];
    double P = solver->inlet->Pin[3];
    double q = 0.5*r*(u*u + v*v);
 
    double Cx_p = 0;
    double Cx_v = 0;       
        
    double Cy_p = 0;
    double Cy_v = 0;               
        
    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {
        bc = solver->mesh->bc[jj];
        if(bc->flagBC == 3)
        {
            for(int ii=0; ii<bc->Nelem; ii++)
            {
                cp = (bc->elemL[ii]->P[3] - P)/q;
                p0 = bc->elemL[ii]->p[0];
                p1 = bc->elemL[ii]->p[1];
                
                meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);                
             
                Cx_p += cp*dSx;
                Cy_p += cp*dSy;                                
                
                if(solver->laminar==1 || solver->sa==1)
                {
                    boundaryCalcFrictionWall(solver, bc->elemL[ii], &fx, &fy);
                    Cx_v -= fx/q;
                    Cy_v -= fy/q;
                } 
                                              
            }
        }
    }
        
    Cx_p /= solver->Sref;
    Cx_v /= solver->Sref;    

    Cy_p /= solver->Sref;
    Cy_v /= solver->Sref;    

    if(solver->mesh->axi)
    {
        Cx_p *= 2*M_PI;
        Cx_v *= 2*M_PI;    
    }

    fprintf(convFile, " %e, %e, %e, %e,", Cx_p, Cx_v, Cy_p, Cy_v); 

}

void solverSetData(SOLVER* solver, INPUT* input)
{

    solver->order = atoi(inputGetValue(input, "order"));
    solver->mesh->order = solver->order;

    //Get boundary conditions
    printf("main: get boundary conditions.\n");
    boundaryGetBC(solver->mesh, input);

    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->Pr = 0.72;
    solver->Pr_t = 0.9;
    solver->Cp = solver->gamma*solver->Rgas/(solver->gamma-1);
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
    
    if(inputNameIsInput(input, "laminar"))
    {
        solver->laminar = atoi(inputGetValue(input, "laminar"));     
    }
    else
    {
        solver->laminar = 0;
    }

    if(inputNameIsInput(input, "restart"))
    {
        solver->restart = atoi(inputGetValue(input, "restart"));     
    }
    else
    {
        solver->restart = 0;
    }

    if(inputNameIsInput(input, "Sref"))
    {
        solver->Sref = strtod(inputGetValue(input, "Sref"), NULL);     
    }
    else
    {
        solver->Sref = 1.0;
    }

    if(inputNameIsInput(input, "turbRatio"))
    {
        solver->turbRatio = strtod(inputGetValue(input, "turbRatio"), NULL);     
    }
    else
    {
        solver->turbRatio = 10.0;
    }

    if(inputNameIsInput(input, "dtLocal"))
    {
        solver->dtLocal = atoi(inputGetValue(input, "dtLocal"));
    }
    else
    {
        solver->dtLocal = 0;
    }

    if(inputNameIsInput(input, "dtLocalN"))
    {
        solver->dtLocalN = atoi(inputGetValue(input, "dtLocalN"));
    }
    else
    {
        solver->dtLocalN = -1;
    }

    // Selection of several variables
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));
    solver->CFL = strtod(inputGetValue(input, "CFL"), NULL);

    // Environmental condition
    if(inputNameIsInput(input, "pout"))
    {
        solver->pout = strtod(inputGetValue(input, "pout"), NULL);     
    }

}

void solverInitDomain(SOLVER* solver)
{
    char s[50];

    printf("main: initialize U.\n");       
    if(atoi(inputGetValue(solver->input, "tube")) == 0)
    {


        solver->inlet = conditionInit(strtod(inputGetValue(solver->input, "pressure"), NULL), 
                                      strtod(inputGetValue(solver->input, "temperature"), NULL), 
                                      strtod(inputGetValue(solver->input, "mach"), NULL), 
                                      strtod(inputGetValue(solver->input, "nx"), NULL),
                                      strtod(inputGetValue(solver->input, "ny"), NULL));

        conditionState(solver->inlet, solver);

        if(solver->restart)                                              
        {
            s[0] = '\0';
            strcat(s, solver->wd);
            strcat(s, "restart.csv");        
            solverLoadRestart(solver, s);
        }
        else
        {
            // Initialization of U
            solverInitU(solver, solver->inlet);
            if(solver->sa == 1)
            {
                saInitU(solver, solver->inlet);
            }

        } 
                
    }
    else
    {
        if(solver->restart)                                              
        {
            s[0] = '\0';
            strcat(s, solver->wd);
            strcat(s, "restart.csv");        
            solverLoadRestart(solver, s);
        }
        else
        {        
            CONDITION* inside1 = conditionInit(strtod(inputGetValue(solver->input, "pressure1"), NULL), 
                                               strtod(inputGetValue(solver->input, "temperature1"), NULL), 
                                               strtod(inputGetValue(solver->input, "mach1"), NULL), 
                                               strtod(inputGetValue(solver->input, "nx1"), NULL),
                                               strtod(inputGetValue(solver->input, "ny1"), NULL));

            CONDITION* inside2 = conditionInit(strtod(inputGetValue(solver->input, "pressure2"), NULL), 
                                               strtod(inputGetValue(solver->input, "temperature2"), NULL), 
                                               strtod(inputGetValue(solver->input, "mach2"), NULL), 
                                               strtod(inputGetValue(solver->input, "nx2"), NULL),
                                               strtod(inputGetValue(solver->input, "ny2"), NULL));      
        
            solverInitUTube(solver, inside1, inside2, strtod(inputGetValue(solver->input, "xm"), NULL));
            free(inside1);
            free(inside2);
        }   
    }
}

SOLVER* solverInit(char* wd)
{
    SOLVER* solver = malloc(sizeof(SOLVER));
    char s[50];

    // Work directory
    solver->wd = wd;

    // Load input   
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "input.dat");
    solver->input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(solver->input);

    // Set number of threads
    omp_set_num_threads(atoi(inputGetValue(solver->input, "threads")));

    // Set turbulence model
    if(inputNameIsInput(solver->input, "sa"))
    {
        solver->sa = atoi(inputGetValue(solver->input, "sa"));
    }
    else
    {
        solver->sa = 0;
    }

    // Set number of flow variables
    solver->Nvar = 4;
    if(solver->sa == 1)
    {
        solver->Nvar = 5;
    }    

    // Load mesh    
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s, solver->Nvar, atoi(inputGetValue(solver->input, "axisymmetric")));
    

    // Setting the solver   
    solverSetData(solver, solver->input);
      
    if(solver->sa == 1)
    {
        printf("main: calculating distance.\n");
        saCalcD(solver->mesh);

    }    
  
    //meshCheckNei(solver->mesh);
    //solverCheckGrad(solver);
    //meshPrint(solver->mesh);
    //meshPrintDStotal(solver->mesh);
    //meshCheckBorderOrientation(solver->mesh);

    // Memory allocation
    solverMalloc(solver);
    
    // Domain initialization
    solverInitDomain(solver);
    
    return solver;
}

void solverSolve(SOLVER* solver)
{

    char s[50];
    double Cx, Cy;
       
    if(atoi(inputGetValue(solver->input, "tube")) == 0)
    {                
        // Calculate time step        
        int Nmax = atoi(inputGetValue(solver->input, "Nmax"));

        // Convergence history file
        s[0] = '\0';
        strcat(s, solver->wd);
        strcat(s, "convergence.csv"); 
        FILE* convFile = fopen(s, "a");

        // Run the solver
        printf("\nmain: running solution:\n");
        for(int ii=0; ii<Nmax; ii++)
        {
            solverCalcDt(solver);
            solverStepRK(solver);
            
            if(ii == solver->dtLocalN)
            {
                solver->dtLocal = 0;
            }
            
            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
                solverCalcCoeff(solver, &Cx, &Cy);
                printf("%f, %f\n", Cx, Cy);
                
                // Write convergence file
                fprintf(convFile, "%i,", ii);
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    fprintf(convFile, " %+.4e,", solver->res[kk]);        
                }                
                solverCalcCoeff3(solver, convFile, ii);
                fprintf(convFile, "\n");
            }            
        }
        solverCalcCoeff2(solver, solver->wd); 
        
        fclose(convFile);
    }
    else
    {   
        double tmax = strtod(inputGetValue(solver->input, "tmax"), NULL);                

        // Run the solver
        double t = 0.0;
        printf("\nmain: running solution:\n");
        int stopLoop = 0;
        int ii = 0;
        while(stopLoop == 0)
        {
            solverCalcDt(solver);
            
            if(t + solver->dt>tmax)
            {
                solver->dt = (tmax-t);
                stopLoop = 1;
            }

            solverStepRK(solver);
            t += solver->dt;
            ii++;

            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
            }
        }
        printf("time %f s\n", t);        
    }
}

