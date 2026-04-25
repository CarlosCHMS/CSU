#ifndef SHOCKTUBE_H
#define SHOCKTUBE_H

typedef struct SOLVER SOLVER;

typedef struct INPUT INPUT;

typedef struct SHOCKTUBE
{

    bool active;
    
    double xm;
    
    double tmax;    
    

} SHOCKTUBE;

SHOCKTUBE* shockTubeInit(INPUT* input);

void shockTubeInitU(SOLVER* solver);

void shockTubeSolve(SOLVER* solver);

#endif

