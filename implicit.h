#ifndef IMPLICIT_H
#define IMPLICIT_H


typedef struct BLOCK{

    int ii;
    double** A;
    struct BLOCK* next;

} BLOCK;


typedef struct IMPLICIT
{
    int timeScheme;

    double wImp;

    double *D;
    double *dtL;    
    
    double **dW0;
    double **dW1;
    
    BLOCK** BB;

} IMPLICIT;

IMPLICIT* implicitInit(INPUT* input, SOLVER* solver);

void implicitFree(IMPLICIT* implicit, SOLVER* solver);

void implicitCalcD(SOLVER* solver);

void implicitLUSGS_L(SOLVER* solver);

void implicitLUSGS_U(SOLVER* solver);

void implicitAuxCalcFlux(SOLVER* solver, double* U, double p, double nx, double ny, double* F);

void implicitCalcDeltaFlux(SOLVER* solver, double* P, double* dW, double nx, double ny, double* dF);

void implicitFunc(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW);

void implicitTest(SOLVER* solver);

void implicitInitDPLUR(SOLVER* solver);

void implicitFreeDPLUR(IMPLICIT* implicit, SOLVER* solver);

void implicitUpdateA(SOLVER* solver);

void implicitUpdateA_sa(SOLVER* solver);

void implicitMultA(SOLVER* solver, double** x, double** y);

void implicitCalcDPLUR(SOLVER* solver);

#endif
