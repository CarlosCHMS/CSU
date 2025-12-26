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
#include"flux.h"
#include"boundary.h"
#include"readTables.h"
#include"sa.h"
#include"saCC.h"
#include"implicit.h"
#include"gasprop.h"
#include"limiter.h"
#include"laminar.h"
#include"sst.h"
#include"sstTrans.h"


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

    double r, u, v, E, c, k;
    
    k = 0.0;
    
    r = cond->p/(solver->gas->R*cond->T);
    c = gasprop_T2c(solver->gas, cond->T);
    u = cond->nx*cond->mach*c;
    v = cond->ny*cond->mach*c;
    
    if(solver->sstFlag == 1)
    {
        double L = solver->sst->L;
        double U = c*cond->mach;
        double ReL = (r*U*L)/sutherland(cond->T);
        
        k = solver->sst->kFactor*U*U/ReL;
        double om = solver->sst->oFactor*U/L;
        
        cond->Uin[4] = r*k;
        cond->Uin[5] = r*om;        
        
        cond->Pin[5] = k;
        cond->Pin[6] = om;                          
        
    }
      
    E = gasprop_T2e(solver->gas, cond->T) + (u*u + v*v)/2 + k;
    
    cond->Uin[0] = r;
    cond->Uin[1] = r*u;
    cond->Uin[2] = r*v;
    cond->Uin[3] = r*E;   

    cond->Pin[0] = r;
    cond->Pin[1] = u;
    cond->Pin[2] = v;
    cond->Pin[3] = solver->gas->R*cond->T*r;
    cond->Pin[4] = cond->T;

    if(solver->sa == 1)
    {
        double n = solver->turbRatio*sutherland(cond->T)/r;
        cond->Uin[4] = r*n;
        cond->Pin[5] = n;         
    }

}

double conditionVref(CONDITION* cond, SOLVER* solver)
{
      
    double c = gasprop_T2c(solver->gas, cond->T);
    double Vref = (cond->mach+1)*c;
        
    return Vref;

}

void solverMalloc(SOLVER* solver)
{
    solver->U = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->R = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->faceFlux = tableMallocDouble(solver->Nvar, solver->mesh->Ncon);
    solver->dPx = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->dPy = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    solver->phi = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);    

    if(solver->dtLocal == 1)
    {
        solver->dtL = malloc(solver->mesh->Nelem*sizeof(double));
    }

    if(solver->timeScheme == 0)
    {
        solver->Uaux = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
    }
    
    if(solver->timeScheme == 1 || solver->timeScheme == 2)
    {
        solver->dW0 = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
        solver->dW1 = tableMallocDouble(solver->Nvar, solver->mesh->Nelem);
        solver->D = malloc(solver->mesh->Nelem*sizeof(double));
        solver->dtL = malloc(solver->mesh->Nelem*sizeof(double));
    }
    
    if(solver->timeScheme == 2)
    {    
        solver->BB = malloc(solver->mesh->Nelem*sizeof(BLOCK*));
    }

    if(solver->sa)
    {
        solver->miT = malloc(solver->mesh->Ncon*sizeof(double));
    }
    
    if(solver->sstFlag)
    {
        solver->miT = malloc(solver->mesh->Ncon*sizeof(double));
        solver->miTe = malloc(solver->mesh->Nelem*sizeof(double));
        solver->F1 = malloc(solver->mesh->Nelem*sizeof(double));
        solver->F2 = malloc(solver->mesh->Nelem*sizeof(double));                
        solver->dd = malloc(solver->mesh->Nelem*sizeof(double));
        solver->om2 = malloc(solver->mesh->Nelem*sizeof(double));
        solver->dQodro = malloc(solver->mesh->Nelem*sizeof(double));
        solver->dQodrk = malloc(solver->mesh->Nelem*sizeof(double));
        solver->dQodr = malloc(solver->mesh->Nelem*sizeof(double));                        
        solver->dQkdro = malloc(solver->mesh->Nelem*sizeof(double));                
        solver->dQkdrk = malloc(solver->mesh->Nelem*sizeof(double));
        solver->dQkdr = malloc(solver->mesh->Nelem*sizeof(double));
    }
}

void solverFree(SOLVER* solver)
{

    if(solver->timeScheme == 0)
    {
        tableFreeDouble(solver->Uaux, solver->Nvar);
    }

    if(solver->timeScheme == 1 || solver->timeScheme == 2)
    {
        tableFreeDouble(solver->dW0, solver->Nvar);
        tableFreeDouble(solver->dW1, solver->Nvar);
        free(solver->D);
        free(solver->dtL);        
    }
    
    if(solver->timeScheme == 2)
    {    
        implicitFreeDPLUR(solver);
        free(solver->BB);
    }

    tableFreeDouble(solver->U, solver->Nvar);
    tableFreeDouble(solver->R, solver->Nvar);        
    tableFreeDouble(solver->faceFlux, solver->Nvar);
    tableFreeDouble(solver->dPx, solver->Nvar);
    tableFreeDouble(solver->dPy, solver->Nvar);
    tableFreeDouble(solver->phi, solver->Nvar);
    meshFree(solver->mesh);
    
    if(solver->dtLocal == 1)
    {
        free(solver->dtL);
    }    

    if(solver->sa)
    {
        free(solver->miT);
    }
    
    if(solver->sstFlag)
    {
        free(solver->miT);
        free(solver->miTe);
        free(solver->F1);
        free(solver->F2);
        free(solver->dd);
        free(solver->om2);
        sstFree(solver->sst);
        free(solver->dQodr);
        free(solver->dQodrk);        
        free(solver->dQodro);
        free(solver->dQkdr);
        free(solver->dQkdrk);        
        free(solver->dQkdro);
    }    
    
    inputFree(solver->input);
    
    gaspropFree(solver->gas);
    
    limiterFree(solver->limiter);
    
    free(solver);

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

    int p0, p1;

    if(solver->laminar || solver->sa || solver->sstFlag)
    {
        for(int ii=0; ii<solver->mesh->Nmark; ii++)
        {
            if(solver->mesh->bc[ii]->flagBC == 3 || solver->mesh->bc[ii]->flagBC == 4)
            {
                for(int jj=0; jj<solver->mesh->bc[ii]->Nelem; jj++)
                {
                    p0 = solver->mesh->bc[ii]->elemL[jj]->p[0];
                    p1 = solver->mesh->bc[ii]->elemL[jj]->p[1];

                    Up[1][p0] = 0.0;
                    Up[2][p0] = 0.0;
                    Up[1][p1] = 0.0;
                    Up[2][p1] = 0.0;
                }
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

    double k = 0.0;
    if(solver->sstFlag)
    {
        k = U[4][ii]/U[0][ii];
    }

	double u = U[1][ii]/U[0][ii];
	double v = U[2][ii]/U[0][ii];
	
    double e = U[3][ii]/U[0][ii] - 0.5*(u*u + v*v) - k;
    double T = gasprop_e2T(solver->gas, e);
    
	return solver->gas->R*T*U[0][ii];
    
}

void solverCalcVel(SOLVER* solver, double** U, int ii, double* u, double* v, double* c)
{
    
    double k = 0;
    if(solver->sstFlag)
    {
        k = U[4][ii]/U[0][ii];
    }
    
    double E = U[3][ii]/U[0][ii];
    double aux;

    *u = U[1][ii]/U[0][ii];
    *v = U[2][ii]/U[0][ii];
    
    aux = (E - ((*u)*(*u) + (*v)*(*v))/2 - k);
    double T = gasprop_e2T(solver->gas, aux);
    *c = gasprop_T2c(solver->gas, T);
    
}

void rotation(double* U, double dSx, double dSy, double dS)
{

    double aux = (U[1]*dSx + U[2]*dSy)/dS;
    U[2] = (-U[1]*dSy + U[2]*dSx)/dS;
    U[1] = aux;
	
}


void solverUpdateGrad(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double dPx, dPy;
        ELEMENT* E = solver->mesh->elemL[ii];
  	    
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            if(kk > 3)
            {
                solverCalcGrad3(solver, E, kk+1, &dPx, &dPy);
            }
            else
            {
                solverCalcGrad3(solver, E, kk, &dPx, &dPy);
            }
            
            solver->dPx[kk][ii] = dPx;
            solver->dPy[kk][ii] = dPy;
        }   
    }
}


void solverGrad_T(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double dPx, dPy;
        ELEMENT* E = solver->mesh->elemL[ii];
  	    
        solverCalcGrad3(solver, E, 4, &dPx, &dPy);
        solver->dPx[3][ii] = dPx;
        solver->dPy[3][ii] = dPy;

    }
}


void inter(SOLVER* solver)
{

    if(solver->order == 2)
	{
	    limiterUpdate(solver->limiter, solver);

        # pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
	        int jj, kk, mm, p0, p1;
            double x0, y0, d2, phi0, phi = 0, xm, ym;
            double Pmin, Pmax; 
            ELEMENT* E = solver->mesh->elemL[ii];
            double Pref2[6];

            limiterCalc(solver->limiter, solver, ii, Pref2);

      	    elementCenter(E, solver->mesh, &x0, &y0);
      	    
            for(kk=0; kk<solver->Nvar; kk++)
            {
                if(kk > 3)
                {
                    mm = kk + 1;
                }
                else
                {
                    mm = kk;
                }
                
                solverCalcMinMax(solver, E, mm, &Pmin, &Pmax);
                
                for(jj=0; jj<E->Np; jj++)
                {
                    p0 = E->p[jj]; 
                    p1 = E->p[(jj+1)%E->Np];
                    xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
                    ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;        
                    d2 = (solver->dPx[kk][ii]*(xm - x0) + solver->dPy[kk][ii]*(ym - y0));
                    phi0 = limiterV2(E->P[mm], Pmin, Pmax, d2, Pref2[kk]);
                    
                    if(jj==0)
                    {
                        phi = phi0;
                    }
                    else
                    {
                        phi = fmin(phi, phi0);
                    }			            
                }

                solver->phi[kk][ii] = phi;
            }   
        }
    }

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
	    int kk;
        double dSx, dSy, dS;
        double x0, y0, x1, y1, xm, ym;
        double PL[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	    double PR[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double f[6];
 
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
		    
		    for(kk=0; kk<solver->Nvar; kk++)
		    {		    
		        if(kk > 3)
		        {
		            if(solver->turb1order)
		            {
		                PL[kk] = E0->P[kk+1];
        			    PR[kk] = E1->P[kk+1];
		            }
		            else
		            {
		                PL[kk] = E0->P[kk+1] + (solver->dPx[kk][e0]*(xm - x0) + solver->dPy[kk][e0]*(ym - y0))*solver->phi[kk][e0];
    			        PR[kk] = E1->P[kk+1] + (solver->dPx[kk][e1]*(xm - x1) + solver->dPy[kk][e1]*(ym - y1))*solver->phi[kk][e1];
    			    }
			    }
			    else
			    {
			        PL[kk] = E0->P[kk] + (solver->dPx[kk][e0]*(xm - x0) + solver->dPy[kk][e0]*(ym - y0))*solver->phi[kk][e0];
    			    PR[kk] = E1->P[kk] + (solver->dPx[kk][e1]*(xm - x1) + solver->dPy[kk][e1]*(ym - y1))*solver->phi[kk][e1];
			    }
		    }
        }
        else
        {       
            for(kk=0; kk<solver->Nvar; kk++)
		    {
		        if(kk > 3)
		        {
		            PL[kk] = E0->P[kk+1];
    			    PR[kk] = E1->P[kk+1];
			    }
			    else
			    {
			        PL[kk] = E0->P[kk];
    			    PR[kk] = E1->P[kk];
			    }
		    }
		}

        // Rotation of the velocity vectors
		rotation(PL, dSx, dSy, dS);
                    
        // Rotation of the velocity vectors
		rotation(PR, dSx, dSy, dS);
        
        // Flux calculation
        if(solver->sstFlag==1)
        {
            flux_sst(solver, PL[0], PL[1], PL[2], PL[3], PL[4], PL[5], PR[0], PR[1], PR[2], PR[3], PR[4], PR[5], f);
        }
        else if(solver->sa==1)
        {
            flux_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4], PR[0], PR[1], PR[2], PR[3], PR[4], f);
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
        for(int jj=0; jj<solver->mesh->elemL[ii]->neiN; jj++)
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

void interAxisPressure(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double dS;        
        dS = meshCalcDSlateral(solver->mesh, ii);
        solver->R[2][ii] -= solver->mesh->elemL[ii]->P[3]*dS;
        if(solver->sstFlag)
        {
            solver->R[2][ii] -= (2./3.)*solver->mesh->elemL[ii]->P[0]*solver->mesh->elemL[ii]->P[5]*dS;
        }
    }
}

void solverCalcR(SOLVER* solver, double** U)
{

    solverResetR(solver);
    solverCalcPrimitive(solver, U);
    solverUpdateGrad(solver);

    inter(solver); 
    
    /*
    for(int ii=0; ii<10; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            printf(" %e,", solver->R[kk][ii]);
        }
        printf("\n");
    }
    */ 
    boundary(solver); 
    
    if(solver->mesh->axi==1)
    {
        interAxisPressure(solver);
    }    
    
    if(solver->laminar==1)
    {
        solverGrad_T(solver);
        laminarInter(solver);
        laminarBoundary(solver);
    }
    
    if(solver->sa==1)
    {
        solverGrad_T(solver);
        if(solver->saCC)
        {
            saCC_Inter(solver);
            saCC_Boundary(solver);
        }
        else
        {
            saInter(solver);
            saBoundary(solver);
        }
    }
    
    if(solver->sstFlag==1)
    {
        solverGrad_T(solver);
        sstInter(solver);
        sstBoundary(solver);
    }
}

void solverRK(SOLVER* solver, double a)
{   
    double sigma = 0.5*solver->stages;
    if(solver->dtLocal == 1)
    {
        # pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {            
                solver->Uaux[kk][ii] = solver->U[kk][ii] - solver->dtL[ii]*sigma*a*solver->R[kk][ii]/solver->mesh->elemL[ii]->omega;
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
                solver->Uaux[kk][ii] = solver->U[kk][ii] - solver->dt*sigma*a*solver->R[kk][ii]/solver->mesh->elemL[ii]->omega;
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
    if(solver->stages==1)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==3)
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

void solverCalcResOld(SOLVER* solver)
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

void solverCalcRes(SOLVER* solver)
{
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        solver->res[kk] = 0.0;
    }
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {         
            solver->res[kk] += solver->R[kk][ii]*solver->R[kk][ii];         
        }
    }

    for(int kk=0; kk<solver->Nvar; kk++)
    {
        solver->res[kk] = sqrt(solver->res[kk]/solver->mesh->Nelem);
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
            solver->dtL[ii] = solver->CFL*solverLocalTimeStep(solver, ii);
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

        dt *= solver->CFL;
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

void solverCalcGrad3(SOLVER* solver, ELEMENT* E, int kk, double* dUx, double* dUy)
{
    /*
    Least squares calculation
    */

    double xx, yy;
    double x[5];
    double y[5];    
    double u[5];    
    long double a, b, c, A, B;

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
}

void solverCalcMinMax(SOLVER* solver, ELEMENT* E, int kk, double* Umin, double* Umax)
{    
    *Umin = E->P[kk];
    *Umax = E->P[kk];
    
    for(int ii=0; ii<E->Np; ii++)
    {
        *Umax = fmax(*Umax, E->neiL[ii]->P[kk]);
        *Umin = fmin(*Umin, E->neiL[ii]->P[kk]);        
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

void solverCalcPrimitive(SOLVER* solver, double** U)
{   
    solver->omLim = 0;
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {   
        ELEMENT* E = solver->mesh->elemL[ii];
          
        E->P[0] = U[0][ii];

        E->P[0] = fmax(E->P[0], solver->rLim);
        U[0][ii] = E->P[0];

        double k = 0;
        if(solver->sstFlag == 1)
        {
            k = U[4][ii]/E->P[0];
            if(k < 1e-14)
            {
                k = 1e-14;
                U[4][ii] = E->P[0]*k;
            }
            E->P[5] = k;
            
            double om = U[5][ii]/E->P[0];
            if(om <= solver->omLim)
            {
                om = solver->omLim;
                U[5][ii] = E->P[0]*om;
            }
            E->P[6] = om;
        }                        

        double u = U[1][ii]/E->P[0];
        double v = U[2][ii]/E->P[0];
        
        E->P[1] = u;
        E->P[2] = v;

        double e = U[3][ii]/E->P[0] - 0.5*(u*u + v*v) - k;
        //printf("\n%e\n", gasprop_T2e(solver->gas, 1));
        e = fmax(e, 7.176325e+02);
        double T = gasprop_e2T(solver->gas, e);

        if(T < 1)
        {
            T = 1;
            e = gasprop_T2e(solver->gas, T);        
            U[3][ii] = (e + 0.5*(u*u + v*v) + k)*E->P[0];
        }    
        E->P[4] = T;

        E->P[3] = solver->gas->R*E->P[4]*E->P[0];

        if(solver->sa == 1)
        {
            E->P[5] = U[4][ii]/E->P[0];
            E->P[5] = fmax(E->P[5], 0);
            U[4][ii] = E->P[5]*E->P[0];
        }                        
    }
    
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        boundaryCalcPrimitive(solver, solver->mesh->bc[ii]);
    }
    
}

double sutherland(double T)
{
    return 1.458e-6*T*sqrt(T)/(T + 110.4);
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
        if(bc->flagBC == 3 || bc->flagBC == 4)
        {
            for(int ii=0; ii<bc->Nelem; ii++)
            {
                cp = (bc->elemL[ii]->P[3] - P)/q;
                p0 = bc->elemL[ii]->p[0];
                p1 = bc->elemL[ii]->p[1];
                
                meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);                
             
                *Cx += cp*dSx;   
                *Cy += cp*dSy;                
                
                if(solver->laminar==1 || solver->sa==1 || solver->sstFlag==1)
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
        if(bc->flagBC == 3 || bc->flagBC == 4)
        {
            for(int ii=0; ii<bc->Nelem; ii++)
            {                
                if(solver->laminar==1 || solver->sa==1 || solver->sstFlag==1)
                {
                    
                    if(solver->sa==1 || solver->sstFlag==1)
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
        if(bc->flagBC == 3 || bc->flagBC == 4)
        {
            for(int ii=0; ii<bc->Nelem; ii++)
            {
                cp = (bc->elemL[ii]->P[3] - P)/q;
                p0 = bc->elemL[ii]->p[0];
                p1 = bc->elemL[ii]->p[1];
                
                meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);                
             
                Cx_p += cp*dSx;
                Cy_p += cp*dSy;                                
                
                if(solver->laminar==1 || solver->sa==1 || solver->sstFlag==1)
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

    int TP;

    if(inputNameIsInput(input, "TP"))
    {
        TP = atoi(inputGetValue(input, "TP"));     
    }
    else
    {
        TP = 0;
    }

    // Constants    
    solver->gas = gaspropInit(1.4, 287.0530, TP);

    //printf("\noi%f\n", gasprop_e2T(solver->gas, gasprop_T2e(solver->gas, 3000)));
    //exit(0);

    solver->Pr = 0.72;
    solver->Pr_t = 0.9;
    solver->eFix = 0.1;
    
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
    
    if(inputNameIsInput(input, "K3"))
    {
        solver->K3 = strtod(inputGetValue(input, "K3"), NULL);     
    }
    else
    {
        solver->K3 = 0.1;
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
    
    if(inputNameIsInput(input, "wImp"))
    {
        solver->wImp = strtod(inputGetValue(input, "wImp"), NULL);
    }
    else
    {
        solver->wImp = 1.0;
    }    

    if(inputNameIsInput(input, "turb1order"))
    {
        solver->turb1order = atoi(inputGetValue(input, "turb1order"));
    }
    else
    {
        solver->turb1order = 0;
    }

    if(inputNameIsInput(input, "timeScheme"))
    {
        solver->timeScheme = solverTimeSchemeChoice(inputGetValue(input, "timeScheme"));
    }
    else
    {
        solver->timeScheme = 0;
    }

    if(inputNameIsInput(input, "Nlinear"))
    {
        solver->Nlinear = atoi(inputGetValue(input, "Nlinear"));
    }
    else
    {
        solver->Nlinear = 20;
    }

    if(inputNameIsInput(input, "rLim"))
    {
        solver->rLim = strtod(inputGetValue(input, "rLim"), NULL);
    }
    else
    {
        solver->rLim = 1e-6;
    }

    if(inputNameIsInput(input, "pLim"))
    {
        solver->pLim = strtod(inputGetValue(input, "pLim"), NULL);
    }
    else
    {
        solver->pLim = 1.0;
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
    else
    {
        solver->pout = 1.0e5;
    }
    
    if(inputNameIsInput(input, "viscBlazek"))
    {
        solver->viscBlazek = atoi(inputGetValue(input, "viscBlazek"));
    }
    else
    {
        solver->viscBlazek = 0;
    }
    
    if(inputNameIsInput(input, "Twall"))
    {
        solver->Twall = strtod(inputGetValue(input, "Twall"), NULL);     
    }
    else
    {
        solver->Twall = 300;
    }
    
    solver->writeSurf[0] = '\0';
    
    if(inputNameIsInput(input, "writeSurf"))
    {
        strcat(solver->writeSurf, inputGetValue(input, "writeSurf"));
    }
    else
    {
        strcat(solver->writeSurf, "wall");
    }
    
    int limType;
    if(inputNameIsInput(solver->input, "limiter"))
    {
        limType = atoi(inputGetValue(solver->input, "limiter"));
    }
    else
    {
        limType = 0;
    }
    
    double limK;
    if(inputNameIsInput(solver->input, "limK"))
    {
        limK = strtod(inputGetValue(solver->input, "limK"), NULL);     
    }
    else
    {
        limK = 1.0;
    }

    solver->limiter = limiterInit(limType, limK, solver);

    if(solver->sstFlag)
    {
        if(inputNameIsInput(solver->input, "kFactor"))
        {
            solver->sst->kFactor = strtod(inputGetValue(solver->input, "kFactor"), NULL);     
        }
        else
        {
            solver->sst->kFactor = 1.125;
        }

        if(inputNameIsInput(solver->input, "oFactor"))
        {
            solver->sst->oFactor = strtod(inputGetValue(solver->input, "oFactor"), NULL);     
        }
        else
        {
            solver->sst->oFactor = 125.0;
        }

        if(inputNameIsInput(solver->input, "oWallFactor"))
        {
            solver->sst->oWallFactor = strtod(inputGetValue(solver->input, "oWallFactor"), NULL);     
        }
        else
        {
            solver->sst->oWallFactor = 10;
        }

        if(inputNameIsInput(solver->input, "Lsst"))
        {
            solver->sst->L = strtod(inputGetValue(solver->input, "Lsst"), NULL);     
        }
        else
        {
            solver->sst->L = 1.0;
        }
    
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
            else if(solver->sstFlag == 1)
            {
                sstInitU(solver, solver->inlet);
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
            solver->inlet = inside1;
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

    if(inputNameIsInput(solver->input, "saCC"))
    {
        solver->saCC = atoi(inputGetValue(solver->input, "saCC"));
    }
    else
    {
        solver->saCC = 0;
    }

    if(inputNameIsInput(solver->input, "sst"))
    {
        solver->sstFlag = atoi(inputGetValue(solver->input, "sst"));
    }
    else
    {
        solver->sstFlag = 0;
    }

    if(solver->sstFlag)
    {
        solver->sst = sstInit();
        if(inputNameIsInput(solver->input, "sstTrans"))
        {
            solver->sst->trans->flag = atoi(inputGetValue(solver->input, "sstTrans"));
        }
        else
        {
            solver->sst->trans->flag = 0;
        }
    }   

    // Set number of flow variables
    solver->Nvar = 4;
    if(solver->sa == 1)
    {
        solver->Nvar = 5;
    }    

    if(solver->sstFlag == 1)
    {
        solver->Nvar = 6;
    }

    // Load mesh    
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s, solver->Nvar, atoi(inputGetValue(solver->input, "axisymmetric")));
    
    // Setting the solver   
    solverSetData(solver, solver->input);
      
    if((solver->sa == 1) || (solver->sstFlag == 1))
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
    
    if(solver->timeScheme == 2)
    {
        implicitInitDPLUR(solver);
    }
        
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
        FILE* convFile;
        
        if(solver->restart)
        {
            convFile = fopen(s, "a");
        }
        else
        {
            convFile = fopen(s, "w");
        }

        // Run the solver
        printf("\nmain: running solution:\n");
        for(int ii=0; ii<Nmax; ii++)
        {
            if(solver->timeScheme == 0)
            {
                solverCalcDt(solver);
                solverStepRK(solver);            
            }
            else if(solver->timeScheme == 1)
            {
                solverCalcR(solver, solver->U);
                implicitCalcD(solver);  
                implicitLUSGS_L(solver);              
                implicitLUSGS_U(solver);    
                solverUpdateUImplicit(solver);        
            }
            else if(solver->timeScheme == 2)
            {
                solverCalcR(solver, solver->U);
                implicitCalcDPLUR(solver);
                solverUpdateUImplicit(solver);
            }
                    
            if(ii == solver->dtLocalN)
            {
                solver->dtLocal = 0;
            }
            
            if(ii%1 == 0)
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
            t += solver->dt*solver->stages/2.0;
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

void solverUpdateUImplicit(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {            
            solver->U[kk][ii] += solver->dW1[kk][ii];
        }
    }
}

int solverTimeSchemeChoice(char* s)
{
    int ans;

    if(strcmp(s, "RK") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "LUSGS") == 0)
    {
        ans = 1;
    }
    else if(strcmp(s, "DPLUR") == 0)
    {
        ans = 2;
    }
    else
    {
        printf("Error in time scheme choice: %s.\n", s);
        exit(0);
    }
    
    return ans;

}


void solverWriteSurf(SOLVER* solver)
{
    if(solver->laminar)
    {
        laminarWriteSurf(solver);
    }

    if(solver->sa)
    {
        if(solver->saCC)
        {
            saCC_SolverWriteSurf(solver);
        }
        else
        {
            saSolverWriteSurf(solver);
        }
    }
    
    if(solver->sstFlag)
    {
        sstSolverWriteSurf(solver);
    }
}

void solverWriteSolution2(SOLVER* solver)
{

    /* 
        Based on: Adek Tasri, Accuracy of Cell Centres to Vertices 
        Interpolation for Unstructured Mesh Finite Volume Solver, 2021
    */

    // Save solution
    printf("main: saving the solution 2.\n");    

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "solution2.csv");
    
    FILE* ff = fopen(fileName, "w");
    int Naux = 3;
    if(solver->sstFlag)
    {
        Naux = 8;
    }
    double** P = tableMallocDouble(solver->Nvar+1, solver->mesh->Np);
    double** Q = tableMallocDouble(Naux, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    int ii, jj, kk, p;
    double xc, yc, xp, yp, L;
    ELEMENT* E;

    solverCalcPrimitive(solver, solver->U);

    for(ii=0; ii<solver->mesh->Np; ii++)
    {       
        for(kk=0; kk<solver->Nvar+1; kk++)
        {
            P[kk][ii] = 0.;
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
           
            for(kk=0; kk<solver->Nvar+1; kk++)
            {        
                P[kk][p] += E->P[kk]/L;
            }                        
            
            den[p] += 1/L;
        }        
        
    }

    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        if(den[ii] != 0)
        {
            for(kk=0; kk<solver->Nvar+1; kk++)
            {
                P[kk][ii] /= den[ii];    
            }
        }
    }

    int p0, p1;

    if(solver->laminar || solver->sa || solver->sstFlag)
    {
        for(int ii=0; ii<solver->mesh->Nmark; ii++)
        {
            if(solver->mesh->bc[ii]->flagBC == 3)
            {
                for(int jj=0; jj<solver->mesh->bc[ii]->Nelem; jj++)
                {
                    p0 = solver->mesh->bc[ii]->elemL[jj]->p[0];
                    p1 = solver->mesh->bc[ii]->elemL[jj]->p[1];

                    P[1][p0] = 0.0;
                    P[2][p0] = 0.0;
                    P[1][p1] = 0.0;
                    P[2][p1] = 0.0;
                    
                    if(solver->sa)
                    {
                        P[5][p0] = 0.0;
                        P[5][p1] = 0.0;
                        if(solver->Twall >= 0)
                        {
                            P[4][p0] = solver->Twall;
                            P[4][p1] = solver->Twall;
                        }
                    }
                    
                    if(solver->sstFlag)
                    {
                        P[5][p0] = 0.0;
                        P[5][p1] = 0.0;
                        //P[6][p0] = 0.0;
                        //P[6][p1] = 0.0;

                        if(solver->Twall >= 0)
                        {
                            P[4][p0] = solver->Twall;
                            P[4][p1] = solver->Twall;
                        }
                    }
                }
            }
        }
    }

    if(solver->sstFlag)
    {
        fprintf(ff, "r,u,v,p,T,k,o,mach,H,s,miEddy,F1,F2,d,om2,\n"); 
        sstInterMiT(solver);
        
        for(ii=0; ii<solver->mesh->Np; ii++)
        {       
            
            Q[3][ii] = 0.;
            Q[4][ii] = 0.;
            Q[5][ii] = 0.;
            Q[6][ii] = 0.;
            Q[7][ii] = 0.;
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
               
                Q[3][p] += solver->miTe[ii]/L;
                Q[4][p] += solver->F1[ii]/L;
                Q[5][p] += solver->F2[ii]/L;
                Q[6][p] += solver->dd[ii]/L;
                Q[7][p] += solver->om2[ii]/L;
                den[p] += 1/L;
            }        
            
        }

        for(ii=0; ii<solver->mesh->Np; ii++)
        {
            if(den[ii] != 0)
            {
                 Q[3][ii] /= den[ii];    
                 Q[4][ii] /= den[ii];                 
                 Q[5][ii] /= den[ii];
                 Q[6][ii] /= den[ii];
                 Q[7][ii] /= den[ii];
            }
        }
    }
    else if(solver->sa)
    {
        fprintf(ff, "r,u,v,p,T,n,mach,H,s,\n");
    }
    else
    {
        fprintf(ff, "r,u,v,p,T,mach,H,s,\n");    
    }

    //fprintf(ff, "0, %i, %i,\n", solver->mesh->Np, solver->Nvar+1+Naux);

    double s0 = gasprop_Tp2entropy(solver->gas, solver->inlet->T, solver->inlet->p);
    for(ii=0; ii<solver->mesh->Np; ii++)
    {   
        //Mach calc     
        double V2 = P[1][ii]*P[1][ii] + P[2][ii]*P[2][ii];
        Q[0][ii] = sqrt(V2)/gasprop_T2c(solver->gas, P[4][ii]);
        Q[1][ii] = gasprop_T2e(solver->gas, P[4][ii]) + P[3][ii]/P[0][ii] + 0.5*V2;
        Q[2][ii] = gasprop_Tp2entropy(solver->gas, P[4][ii], P[3][ii]) - s0;
    }
    
    for(ii=0; ii<solver->mesh->Np; ii++)
    {        
        for(jj=0; jj<solver->Nvar+1; jj++)
        {
            fprintf(ff, "%.10e, ", P[jj][ii]);
        }
                
        for(jj=0; jj<Naux; jj++)
        {
            fprintf(ff, "%.10e, ", Q[jj][ii]);
        }        
        fprintf(ff, "\n");        
    }

    fclose(ff);

    tableFreeDouble(P, solver->Nvar+1);
    tableFreeDouble(Q, Naux);
    free(den);
    
}

