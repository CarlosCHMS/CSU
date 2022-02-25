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
        if(den[ii] != 0)
        {
            for(jj=0; jj<4; jj++)
            {
                Up[jj][ii] /= den[ii];    
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
    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
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

void inter(SOLVER* solver, double **U)
{
    
	int ii, kk;
    double dSx, dSy, dS;
    double aux;
    double UL[4];
	double UR[4];
    double f[4];
    double delta;
    int e0, e1, p0, p1;

    for(ii=0; ii<solver->mesh->Ncon; ii++)
    {
 
        e0 = solver->mesh->con[ii][0];
        e1 = solver->mesh->con[ii][1];
        p0 = solver->mesh->con[ii][2];
        p1 = solver->mesh->con[ii][3];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			UL[kk] = U[kk][e0];
		}
        
        for(kk=0; kk<4; kk++)
		{
			UR[kk] = U[kk][e1];
		}
		
        // Rotation of the velocity vectors
		rotation(UL, dSx, dSy, dS);
                    
        // Rotation of the velocity vectors
		rotation(UR, dSx, dSy, dS);
        
        // Flux calculation
        flux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);

        // Rotation of the flux
		rotation(f, dSx, -dSy, dS);
                         
        for(kk=0; kk<4; kk++)
        {
            aux = f[kk]*dS;
            solver->R[kk][e0] += aux;
            solver->R[kk][e1] -= aux;
        } 
    }
}

void boundaryCalc(SOLVER* solver, double **U, MESHBC* bc)
{
    
	int kk;
    double dSx, dSy, dS;
    double aux;
    double UL[4];
	double UR[4];
    double f[4];
    double delta;
    int e0, p0, p1;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
 
        e0 = bc->domain[ii];
        p0 = bc->elem[ii][0];
        p1 = bc->elem[ii][1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			UL[kk] = U[kk][e0];
		}      		
        
        if(bc->flagBC == 0)
        {

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
        
            // Reflexive
            flux(solver, UL[0], UL[1], UL[2], UL[3], UL[0], -UL[1], UL[2], UL[3], f);

        }
        else if(bc->flagBC == 1)
        {

            for(kk=0; kk<4; kk++)
		    {
			    UR[kk] = solver->inlet->Uin[kk];
		    }

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
            rotation(UR, dSx, dSy, dS);

            flux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);
        
        }
        else if(bc->flagBC == 2)
        {           

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
            fluxFree(solver, UL[0], UL[1], UL[2], UL[3], f);
		
        }
        else if(bc->flagBC == 3)
        {
        
            /*
            Based on: Blazek J., Computacional Fluid Dynamics, Principles and Applications (2001)
            */        

		    double p2 = solverCalcP(solver, U, e0);
	
            // Outlet
            f[0] = .0;
            f[2] = .0;
            f[3] = .0;

            f[1] =  p2;
        }       

        // Rotation of the flux
		rotation(f, dSx, -dSy, dS);
                
        for(kk=0; kk<4; kk++)
        {
            aux = f[kk]*dS;
            solver->R[kk][e0] += aux;
        } 
    }
}

void boundary(SOLVER* solver, double **U)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        boundaryCalc(solver, U, solver->mesh->bc[ii]);
    }
}

void solverCalcR(SOLVER* solver, double** U)
{

    solverResetR(solver);

    inter(solver, U);
    
    boundary(solver, U); 

}

void solverRK(SOLVER* solver, double a)
{   
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

    solverCalcR(solver, solver->U);
    solverRK(solver, 1.0);

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

int boundaryChoice(char* s)
{
    int ans;

    if(strcmp(s, "symmetry") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "inlet") == 0)
    {
        ans = 1;
    }
    else if(strcmp(s, "outlet") == 0)
    {
        ans = 2;
    }
    else if(strcmp(s, "wall") == 0)
    {
        ans = 3;
    }
 
    return ans;

}

void boundaryGetBC(MESH* mesh, INPUT* input)
{
    char s[50];
    for(int ii=0; ii<mesh->Nmark; ii++)
    {
        s[0] = '\0';
        strcat(s, "BC:");
        strcat(s, mesh->bc[ii]->name);
        mesh->bc[ii]->flagBC = boundaryChoice(inputGetValue(input, s));

    }
}
