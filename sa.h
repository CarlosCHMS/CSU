#ifndef SA_H
#define SA_H

typedef struct SOLVER SOLVER;

typedef struct CONDITION CONDITION;

typedef struct MESHBC MESHBC;

typedef struct SA
{
    bool active;

    double Cv1;
    double Cv1_3;    
    double sig;
    double k;
    double cv2;
    double cv3;
    double Cb1;
    double Cb2;
    double Cw1;
    double Cw2;
    double Cw3;
    double Cw3_6;

} SA;

SA* saInit();

void saFree(SA* sa);

void saInitU(SOLVER* solver, CONDITION* inside);

void saInterFaceB(SOLVER* solver);

void saInterFace(SOLVER* solver);

void saInterSource(SOLVER* solver);

void saInter(SOLVER* solver);

void saCalcFace(SA* sa, double ni, double ni_L, double r, double dnix, double dniy, double* fv1, double* tx, double* ty);

void saCalcSource(SA* sa, double ni, double ni_L, double S, double d, double rho, double drx, double dry, double dnix, double dniy, double* Qt);

void saBoundaryFace(SOLVER* solver, MESHBC* bc);

void saBoundary(SOLVER* solver);

void saBoundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f, double* miEddy);

void saSolverWriteSurf(SOLVER* solver);

#endif
