#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"utils.h"
#include"mesh.h"
#include"solver.h"

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
    

}

double conditionVref(CONDITION* cond, SOLVER* solver)
{
      
    double c = sqrt(solver->gamma*solver->Rgas*cond->T);
    double Vref = (cond->mach+1)*c;
        
    return Vref;

}

void solverFree(SOLVER* solver)
{

    meshFree(solver->mesh);
    tableFreeDouble(solver->U, 4);
    tableFreeDouble(solver->Uaux, 4);
    tableFreeDouble(solver->R, 4);        

}

void solverWrite(SOLVER* solver, char* fileName)
{

    FILE* ff = fopen(fileName, "w");
    double** Up = tableMallocDouble(4, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    double aux;
    int ii, jj;

    for(ii=0; ii<solver->mesh->Np; ii++)
    {       
        for(jj=0; jj<4; jj++)
        {
            Up[jj][ii] = 0.;
    
        }   
        den[ii] = 0.;
    }

    for(ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(jj=0; jj<4; jj++)
        {
            Up[jj][solver->mesh->elem[ii][0]] += solver->U[jj][ii];
            Up[jj][solver->mesh->elem[ii][1]] += solver->U[jj][ii];
            Up[jj][solver->mesh->elem[ii][2]] += solver->U[jj][ii];
    
        }

        den[solver->mesh->elem[ii][0]] += 1;
        den[solver->mesh->elem[ii][1]] += 1;
        den[solver->mesh->elem[ii][2]] += 1;
    }

    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        for(jj=0; jj<4; jj++)
        {
            Up[jj][ii] /= den[ii] ;
    
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
    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            solver->R[kk][ii] = 0.0; 
        }
    }
}

double solverCalcP(SOLVER* solver, double*** U, int ii, int jj)
{

	double u = U[1][ii][jj]/U[0][ii][jj];
	double v = U[2][ii][jj]/U[0][ii][jj];
    
	return (solver->gamma - 1)*(U[3][ii][jj] - 0.5*(u*u + v*v)*U[0][ii][jj]);
    
}

void solverCalcVel(SOLVER* solver, double*** U, int ii, int jj, double* u, double* v, double* c)
{
    
    double E = U[3][ii][jj]/U[0][ii][jj];
    double aux;

    *u = U[1][ii][jj]/U[0][ii][jj];
    *v = U[2][ii][jj]/U[0][ii][jj];
    
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
