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

    // Set number of threads
    //omp_set_num_threads(atoi(inputGetValue(input, "threads")));

    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s); 
    
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

    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
        
    // Selection of several variables
    solver->MUSCL = atoi(inputGetValue(input, "order")) - 1;
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));
    solver->CFL = strtod(inputGetValue(input, "CFL"), NULL);

    // Environmental condition
    printf("main: Initialize U.\n");
    solver->inlet = conditionInit(strtod(inputGetValue(input, "pressure"), NULL), 
                                  strtod(inputGetValue(input, "temperature"), NULL), 
                                  strtod(inputGetValue(input, "mach"), NULL), 
                                  strtod(inputGetValue(input, "nx"), NULL),
                                  strtod(inputGetValue(input, "ny"), NULL));
    
    conditionState(solver->inlet, solver);
    solverInitU(solver, solver->inlet);        
    
    if(inputNameIsInput(input, "pout"))
    {
        solver->pout = strtod(inputGetValue(input, "pout"), NULL);     
    }
    
    //Integration       
    printf("main: start solution calculation.\n");
    int Nmax = atoi(inputGetValue(input, "Nmax"));
    for(int ii=0; ii<Nmax; ii++)
    {
        solver->dt = solver->CFL*0.5*solver->stages*solverCalcDt(solver);        
        solverStepRK(solver);

        if(ii%100==0)
        {
            printf("%i, ", ii);
            solverCalcRes(solver);
        }
    }

    // Save solution
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWrite(solver, s);
      
    solverFree(solver);
    inputFree(input);    

    return 0;

}
