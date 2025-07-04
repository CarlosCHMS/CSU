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
#include"flux.h"
#include"boundary.h"
#include"sa.h"
#include"implicit.h"


int main(int argc, char **argv)
{

    struct timeval start, stop;
    
    gettimeofday(&start,NULL);
   
    SOLVER* solver = solverInit(argv[1]);
    
    solverSolve(solver);
        
    solverWriteSolution(solver);
    
    solverWriteSolution2(solver);

    solverWriteReestart(solver);
    
    solverWriteSurf(solver);
      
    solverFree(solver); 

    gettimeofday(&stop,NULL);
    
    printf("\nDuration %f s\n", duration(start, stop));

    return 0;

}
