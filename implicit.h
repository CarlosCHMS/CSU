#ifndef IMPLICIT_H
#define IMPLICIT_H


typedef struct BLOCK{

    int ii;
    double** A;
    struct BLOCK* next;

} BLOCK;


typedef struct IMPLICIT
{
    bool isMatrix;
    
    int timeScheme;

    double wImp;

    double *D;
    double *dtL;
    
    double **dW0;
    double **dW1;
    double **dW10;
    double **U0;
    double **R0;
    double **r;
    double **w;
    
    double ***v;
    
    BLOCK** BB;
    
    void (*LUSGSinv)(SOLVER*, double**, double**, double);

} IMPLICIT;

IMPLICIT* implicitInit(INPUT* input, SOLVER* solver);

void implicitMalloc(IMPLICIT* implicit, SOLVER* solver);

void implicitFree(IMPLICIT* implicit, SOLVER* solver);

void implicitCalcD(SOLVER* solver);

void implicitInitMatrix(SOLVER* solver);

void implicitFreeMatrix(IMPLICIT* implicit, SOLVER* solver);

void implicitUpdateA(SOLVER* solver);

void implicitMultA2(SOLVER* solver, double** U0, double** R0, double** x, double** y);

void implicitCalcLUSGS_matrix(SOLVER* solver, double** b, double** dW1, double sig);

void implicitDPLUR(SOLVER* solver);

void implicitCopy(SOLVER* solver, double** x0, double** x1);

void implicitGMRES(SOLVER* solver);

void implicitGMRES_solveMinimization(int dim, double beta, double** H, double* y);

void implicitCalcJacobi(SOLVER* solver, double* U, double** A, double nx, double ny, double T);

void implicitCalcLUSGS_matrixFree(SOLVER* solver, double** b, double** dW1, double sig);

void implicitLUSGS(SOLVER* solver);

#endif
