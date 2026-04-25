#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<stdbool.h>
#include<omp.h>
#include"shockTube.h"
#include"solver.h"
#include"input.h"
#include"mesh.h"

SHOCKTUBE* shockTubeInit(INPUT* input)
{

    SHOCKTUBE* shockTube = malloc(sizeof(SHOCKTUBE));
    
    if(inputNameIsInput(input, "tube"))
    {
        shockTube->active = atoi(inputGetValue(input, "tube"));     
    }
    else
    {
        shockTube->active = false;
    }

    if(inputNameIsInput(input, "xm"))
    {
        shockTube->xm = strtod(inputGetValue(input, "xm"), NULL);     
    }
    else
    {
        shockTube->xm = 1;
    }

    if(inputNameIsInput(input, "tmax"))
    {
        shockTube->tmax = strtod(inputGetValue(input, "tmax"), NULL);     
    }
    else
    {
        shockTube->tmax = 1;
    }

    return shockTube;

}

void shockTubeInitU(SOLVER* solver)
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

    conditionState(inside1, solver);
    conditionState(inside2, solver);
    
    double x, y;

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
    	    elementCenter(solver->mesh->elemL[ii], solver->mesh, &x, &y);
            if(x < solver->shockTube->xm)
            {
                solver->U[kk][ii] = inside1->Uin[kk];
            }
            else
            {
                solver->U[kk][ii] = inside2->Uin[kk];
            }
        }
    }        
        
    solver->inlet = inside1;
    free(inside2);
}

void shockTubeSolve(SOLVER* solver)
{
    char s[50];
    double tmax = solver->shockTube->tmax;

    // Convergence history file
    strcpy(s, solver->wd);
    strcat(s, "convergence.csv"); 
    FILE* convFile;
    
    if(solver->restart)
    {
        convFile = fopen(s, "a");
    }
    else
    {
        convFile = fopen(s, "w");
        solverPrintConvReader(solver, convFile);
    }

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

        if(ii%1 == 0)
        {
            printf("%i, ", ii);
            solverCalcRes(solver);
            
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
    printf("time %f s\n", t);        
    
    fclose(convFile);
}

