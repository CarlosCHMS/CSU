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

LIMITER* limiterInit(int type, double K, SOLVER* solver)
{

    LIMITER* limiter = malloc(sizeof(LIMITER));
    limiter->type = type;
    limiter->K = K;
    limiter->Pref20 = malloc(solver->Nvar*sizeof(double));
    
    return limiter;
}

void limiterFree(LIMITER* limiter)
{

    free(limiter->Pref20);
    free(limiter->Pref2);
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
        limiter->volMax = 0.0;
	    for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
	        limiter->volMax = fmax(solver->mesh->elemL[ii]->omega, limiter->volMax);
	    }        
        limiter->Pref20[0] = solver->inlet->Pin[0]*solver->inlet->Pin[0];
	    limiter->Pref20[1] = solver->inlet->Pin[1]*solver->inlet->Pin[1] + solver->inlet->Pin[2]*solver->inlet->Pin[2];
	    limiter->Pref20[2] = limiter->Pref20[1];
	    limiter->Pref20[3] = solver->inlet->Pin[3]*solver->inlet->Pin[3];
	    if(solver->Nvar == 5)
	    {
	        limiter->Pref20[4] = solver->inlet->Pin[5]*solver->inlet->Pin[5];
	    }
    }
    else if(limiter->type == 1)
    {
   
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            mm = jj;
            if(jj==4)
            {
                mm = 5;
            }
            Pmin[jj] = solver->mesh->elemL[0]->P[mm];
            Pmax[jj] = solver->mesh->elemL[0]->P[mm]; 
        }
        
	    for(int ii=1; ii<solver->mesh->Nelem; ii++)
        {
            for(int jj=0; jj<solver->Nvar; jj++)
            {
                mm = jj;
                if(jj==4)
                {
                    mm = 5;
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
        aux = limiter->K*sqrt(solver->mesh->elemL[ii]->omega/limiter->volMax);
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

}

double limiterV2(double Ui, double Umin, double Umax, double d2, double ee)
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

