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
#include"sa.h"

void mainSetData(SOLVER* solver, INPUT* input)
{

    solver->order = atoi(inputGetValue(input, "order"));
    solver->mesh->order = solver->order;
    solver->mesh->axi = atoi(inputGetValue(input, "axisymmetric"));

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
        solver->Sref = atoi(inputGetValue(input, "Sref"));     
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

int main(int argc, char **argv)
{

    char s[50];
    double Cx, Cy;
    SOLVER* solver = malloc(sizeof(SOLVER));

    // Load input   
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "input.dat");
    INPUT* input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(input);

    // Set number of threads
    omp_set_num_threads(atoi(inputGetValue(input, "threads")));

    // Set turbulence model
    if(inputNameIsInput(input, "sa"))
    {
        solver->sa = atoi(inputGetValue(input, "sa"));
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
    strcat(s, argv[1]);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s, solver->Nvar);
    

    // Setting the solver   
    mainSetData(solver, input);
      
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
       
    if(atoi(inputGetValue(input, "tube")) == 0)
    {

        printf("main: initialize U.\n");
        solver->inlet = conditionInit(strtod(inputGetValue(input, "pressure"), NULL), 
                                      strtod(inputGetValue(input, "temperature"), NULL), 
                                      strtod(inputGetValue(input, "mach"), NULL), 
                                      strtod(inputGetValue(input, "nx"), NULL),
                                      strtod(inputGetValue(input, "ny"), NULL));

        conditionState(solver->inlet, solver);

        if(solver->restart)                                              
        {
            s[0] = '\0';
            strcat(s, argv[1]);
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
                
        // Calculate time step        
        int Nmax = atoi(inputGetValue(input, "Nmax"));

        // Run the solver
        printf("\nmain: running solution:\n");
        for(int ii=0; ii<Nmax; ii++)
        {
            solver->dt = solverCalcDt(solver);
            solverStepRK(solver);
            
            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
                solverCalcCoeff(solver, &Cx, &Cy);
                printf("%f, %f\n", Cx, Cy);
            }            
        }
        solverCalcCoeff2(solver, argv[1]); 
    }
    else
    {
        printf("main: initialize U.\n");
        if(solver->restart)                                              
        {
            s[0] = '\0';
            strcat(s, argv[1]);
            strcat(s, "restart.csv");        
            solverLoadRestart(solver, s);
        }
        else
        {        
            CONDITION* inside1 = conditionInit(strtod(inputGetValue(input, "pressure1"), NULL), 
                                               strtod(inputGetValue(input, "temperature1"), NULL), 
                                               strtod(inputGetValue(input, "mach1"), NULL), 
                                               strtod(inputGetValue(input, "nx1"), NULL),
                                               strtod(inputGetValue(input, "ny1"), NULL));

            CONDITION* inside2 = conditionInit(strtod(inputGetValue(input, "pressure2"), NULL), 
                                               strtod(inputGetValue(input, "temperature2"), NULL), 
                                               strtod(inputGetValue(input, "mach2"), NULL), 
                                               strtod(inputGetValue(input, "nx2"), NULL),
                                               strtod(inputGetValue(input, "ny2"), NULL));      
        
            solverInitUTube(solver, inside1, inside2, strtod(inputGetValue(input, "xm"), NULL));
            free(inside1);
            free(inside2);
        }
        
        double tmax = strtod(inputGetValue(input, "tmax"), NULL);                

        // Run the solver
        double t = 0.0;
        printf("\nmain: running solution:\n");
        int stopLoop = 0;
        int ii = 0;
        while(stopLoop == 0)
        {
            solver->dt = solverCalcDt(solver);
            
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

    // Save solution
    printf("main: saving the solution.\n");    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWrite(solver, s);

    // Save reestart
    printf("main: saving the restart file.\n");    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "restart.csv");
    solverWriteReestart(solver, s);
      
    solverFree(solver);
    inputFree(input);    

    return 0;

}
