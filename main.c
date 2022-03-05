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


int main(int argc, char **argv)
{

    char s[50];
    SOLVER* solver = malloc(sizeof(SOLVER));

    // Load input   
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "input.dat");
    INPUT* input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(input);

    solver->order = atoi(inputGetValue(input, "order"));

    // Set number of threads
    omp_set_num_threads(atoi(inputGetValue(input, "threads")));

    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s); 
    
    solver->mesh->order = solver->order;
    if(solver->mesh->order == 2)
    {
        meshCalcNeighbors(solver->mesh);
    }
    
    //meshCheckNei(solver->mesh);
    
    //solverCheckGrad(solver);
    
    // axisymmetric
    solver->mesh->axi = atoi(inputGetValue(input, "axisymmetric"));
    
    //meshPrint(solver->mesh);
    //meshPrintDStotal(solver->mesh);

    //Get boundary conditions
    printf("main: get boundary conditions.\n");
    boundaryGetBC(solver->mesh, input);

    // Memory allocation
    solver->U = tableMallocDouble(4, solver->mesh->Nelem);
    solver->Uaux = tableMallocDouble(4, solver->mesh->Nelem);    
    solver->R = tableMallocDouble(4, solver->mesh->Nelem);
    solver->faceFlux = tableMallocDouble(4, solver->mesh->Ncon);
    solver->dPx = tableMallocDouble(4, solver->mesh->Nelem);
    solver->dPy = tableMallocDouble(4, solver->mesh->Nelem);    
    solverMallocP(solver);      

    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
        
    // Selection of several variables
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));
    solver->CFL = strtod(inputGetValue(input, "CFL"), NULL);

    // Environmental condition
    if(inputNameIsInput(input, "pout"))
    {
        solver->pout = strtod(inputGetValue(input, "pout"), NULL);     
    }
       
    if(atoi(inputGetValue(input, "tube")) == 0)
    {

        printf("main: initialize U.\n");
        solver->inlet = conditionInit(strtod(inputGetValue(input, "pressure"), NULL), 
                                      strtod(inputGetValue(input, "temperature"), NULL), 
                                      strtod(inputGetValue(input, "mach"), NULL), 
                                      strtod(inputGetValue(input, "nx"), NULL),
                                      strtod(inputGetValue(input, "ny"), NULL));
        
        // Initialization of U
        solverInitU(solver, solver->inlet);
        
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
            }            
        } 
    }
    else
    {
        printf("main: initialize U.\n");
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
      
    solverFree(solver);
    inputFree(input);    

    return 0;

}
