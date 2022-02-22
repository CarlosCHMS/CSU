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

void test(SOLVER* solver)
{

    double x, y, aux;

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {   
        meshElemCenter(solver->mesh, ii, &x, &y);
        solver->U[0][ii] = cos(1.5*x)*cos(1.5*y);
    }    

}


int main(int argc, char **argv)
{

    char s[50];
    char* ss;
    SOLVER* solver = malloc(sizeof(SOLVER));

    // Load input   
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "input.dat");
    INPUT* input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(input);

    //ss = inputGetValue2(input, "interpE");

    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s);    

    // Memory allocation
    solver->U = tableMallocDouble(4, solver->mesh->Nelem);
    solver->Uaux = tableMallocDouble(4, solver->mesh->Nelem);    
    solver->R = tableMallocDouble(4, solver->mesh->Nelem);

    // Inlet condition
    solver->inlet = conditionInit(strtod(inputGetValue(input, "pressure"), NULL), 
                                      strtod(inputGetValue(input, "temperature"), NULL), 
                                      strtod(inputGetValue(input, "mach"), NULL), 
                                      strtod(inputGetValue(input, "nx"), NULL),
                                      strtod(inputGetValue(input, "ny"), NULL));
            
    solverInitU(solver, solver->inlet);        
                                      
    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
    
    if(inputNameIsInput(input, "CFL"))
    {
        solver->CFL = strtod(inputGetValue(input, "CFL"), NULL);     
    }
    else
    {
        solver->CFL = 1.0;
    }
    
    // Seletion of MUSCL and flux
    solver->MUSCL = atoi(inputGetValue(input, "order")) - 1;
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<4; jj++)
        {
            solver->U[jj][ii] = 0.0;
    
        }

    } 

    test(solver);

    //meshPrint(solver->mesh);
    
    // Save solution
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWrite(solver, s);
        
    solverFree(solver);
    inputFree(input);    

    return 0;

}
