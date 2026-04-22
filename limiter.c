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
#include"limiter.h"

LIMITER* limiterInit(INPUT* input, SOLVER* solver)
{

    LIMITER* limiter = malloc(sizeof(LIMITER));
    
    //Determination of the limiter function
    limiter->name = malloc(50*sizeof(char));
    limiter->name[0] = '\0';
    
    if(inputNameIsInput(input, "limiterName"))    
    {
        strcat(limiter->name, inputGetValue(input, "limiterName"));
    }
    else
    {
        strcat(limiter->name, "VK");
    }
    
    if(strcmp(limiter->name, "VK") == 0)
    {        
        limiter->func = limiterFuncV2;
    }
    else if(strcmp(limiter->name, "BJ") == 0)
    {
        limiter->func = limiterFuncBJ;
    }
    else
    {
        printf("\nError: Limiter name incorrect.");
        exit(0);
    }    
        
    // Determination of the limiter factor for Venkatakrishnan limiter    
    if(inputNameIsInput(input, "limiter"))
    {
        limiter->type = atoi(inputGetValue(input, "limiter"));
    }
    else
    {
        limiter->type = 0;
    }
    
    // Determination of the K parameter for Venkatakrishnan limiter
    if(inputNameIsInput(input, "limK"))
    {
        limiter->K = strtod(inputGetValue(input, "limK"), NULL);     
    }
    else
    {
        limiter->K = 1.0;
    }

    return limiter;
}

void limiterMalloc(LIMITER* limiter, SOLVER* solver)
{
    limiter->Pref20 = malloc(solver->Nvar*sizeof(double));    
    limiter->phi = tableMallocDouble(solver->Nvar, solver->mesh->Nelem); 
}

void limiterFree(LIMITER* limiter, SOLVER* solver)
{
    tableFreeDouble(limiter->phi, solver->Nvar);
    free(limiter->Pref20);    
    free(limiter);
}

void limiterUpdate(LIMITER* limiter, SOLVER* solver)
{

    double Pmin[5];
    double Pmax[5];
    int mm;
    double aux;

    if(limiter->type == 0)
    {        
        limiter->Pref20[0] = solver->inlet->Pin[0]*solver->inlet->Pin[0];
	    limiter->Pref20[1] = solver->inlet->Pin[1]*solver->inlet->Pin[1] + solver->inlet->Pin[2]*solver->inlet->Pin[2];
	    limiter->Pref20[2] = limiter->Pref20[1];
	    limiter->Pref20[3] = solver->inlet->Pin[3]*solver->inlet->Pin[3];
	    if(solver->Nvar > 4)
	    {
	        for(int kk = 4; kk<solver->Nvar; kk++)
	        {
	            limiter->Pref20[kk] = solver->inlet->Pin[kk+1]*solver->inlet->Pin[kk+1];
	        }
	    }
    }
    else if(limiter->type == 1)
    {
   
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            mm = jj;
            if(jj>3)
            {
                mm = jj + 1;
            }
            Pmin[jj] = solver->mesh->elemL[0]->P[mm];
            Pmax[jj] = solver->mesh->elemL[0]->P[mm]; 
        }
        
	    for(int ii=1; ii<solver->mesh->Nelem; ii++)
        {
            for(int jj=0; jj<solver->Nvar; jj++)
            {
                mm = jj;
                if(jj>3)
                {
                    mm = jj+1;
                }

                Pmin[jj] = fmin(Pmin[jj], solver->mesh->elemL[ii]->P[mm]);
                Pmax[jj] = fmax(Pmax[jj], solver->mesh->elemL[ii]->P[mm]);                            
            }
	    }        

        for(int jj=0; jj<solver->Nvar; jj++)
        {
            aux = (Pmax[jj] - Pmin[jj])*limiter->K;
            limiter->Pref20[jj] = aux*aux;
        }
    }
}

void limiterCalc(LIMITER* limiter, SOLVER* solver, int ii, double* Pref2)
{
    double aux;
    if(limiter->type == 0)
    {
        aux = limiter->K*sqrt(solver->mesh->omega[ii]/solver->mesh->volMax);
        aux = aux*aux*aux;
        
	    for(int jj=0; jj<solver->Nvar; jj++)
        {
            Pref2[jj] = limiter->Pref20[jj]*aux;
        }    
    }
    else if(limiter->type == 1)
    {
	    for(int jj=0; jj<solver->Nvar; jj++)
        {
            Pref2[jj] = limiter->Pref20[jj];
        }    
    }
    
    limiter->Pref = Pref2;

}

double limiterFuncV2(LIMITER* limiter, double Ui, double Umin, double Umax, double d2, double ee)
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
        ans = (d1max*d1max + ee)*d2 + 2*d2*d2*d1max;
        ans /= d1max*d1max + 2*d2*d2 + d1max*d2 + ee;
        ans /= d2;
    }
    else
    {
        ans = (d1min*d1min + ee)*d2 + 2*d2*d2*d1min;
        ans /= d1min*d1min + 2*d2*d2 + d1min*d2 + ee;
        ans /= d2;
    }
    
    return ans;

}

double limiterFuncBJ(LIMITER* limiter, double Ui, double Umin, double Umax, double d2, double ee)
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

