#ifndef BOUNDARY_H
#define BOUNDARY_H

typedef struct BOUNDARY
{

    MESHBC* bc;
    char type[50];
    
    void (*primitive) (BOUNDARY*, SOLVER*);
    void (*convective) (BOUNDARY*, SOLVER*);
    void (*viscousFlux) (BOUNDARY*, SOLVER*, int, double*, double*);

} BOUNDARY;

BOUNDARY* boundaryInit(INPUT* input, MESHBC* bc, SOLVER* solver);

void boundaryInlet(SOLVER* solver, double* Pa, double* Pd, double* Pb, double nx, double ny);

void boundaryOutlet(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);

void boundaryWall(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);

void boundary1(SOLVER* solver);

void boundaryGetBC(MESH* mesh, INPUT* input);

int boundaryChoice(char* s);

void boundaryCalcPrimitive(SOLVER* solver, MESHBC* bc);

void boundaryCalcFrictionWall(SOLVER* solver, ELEMENT* E, double* fx, double* fy);

void boundaryPrimitiveSymmetry(BOUNDARY* boundary, SOLVER* solver);

void boundaryPrimitiveInlet(BOUNDARY* boundary, SOLVER* solver);

void boundaryPrimitiveOutlet(BOUNDARY* boundary, SOLVER* solver);

void boundaryPrimitiveWall(BOUNDARY* boundary, SOLVER* solver);

void boundaryPrimitiveWallT(BOUNDARY* boundary, SOLVER* solver);

void boundaryConvectiveSymmetry(BOUNDARY* boundary, SOLVER* solver);

void boundaryConvectiveGeneral(BOUNDARY* boundary, SOLVER* solver);

void boundaryViscous(BOUNDARY* boundary, SOLVER* solver);

#endif
