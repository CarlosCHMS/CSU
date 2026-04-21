#ifndef SOLVER_H
#define SOLVER_H

typedef struct INPUT INPUT;

typedef struct MESH MESH;

typedef struct ELEMENT ELEMENT;

typedef struct GASPROP GASPROP;

typedef struct FLUX FLUX;

typedef struct LIMITER LIMITER;

typedef struct SA SA;

typedef struct SST SST;

typedef struct IMPLICIT IMPLICIT;

typedef struct CONDITION{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[6];  
    double Pin[7];     

} CONDITION;

typedef struct SOLVER{

    int Nvar;
    int Nrow;
    int Ncol;
    int pOutFlag;
    int order;
    int stages;
    int laminar;
    int restart;   
    int dtLocal;
    int turb1order;
    int timeScheme;
    int Nlinear;
    int viscBlazek;
    int tube;

    char* wd;
    char writeSurf[50];

    double k4;
    double dt;
    double pout; 
    double k; 
    double res[6];
    double CFL;
    double Sref;
    double dtLocalN;
    double tol;
    double rLim;
    double pLim; 
    double Twall;
    double omLim;       

    double *dtL;
    double *miT;
            
    double **U;
    double **R;
    double **Uaux;     
    double **faceFlux;       
    double **dPx;
    double **dPy;
    double **uD;
    double **rD;
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    INPUT* input;
    
    FLUX* flux1;

    GASPROP* gas;
    
    LIMITER* limiter;

    SA* sa1;
    
    SST* sst;
    
    IMPLICIT* implicit;

} SOLVER;

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny);

void conditionState(CONDITION* cond, SOLVER* solver);

double conditionVref(CONDITION* cond, SOLVER* solver);

void solverMalloc(SOLVER* solver);

void solverFree(SOLVER* solver);

void solverWriteSolution(SOLVER* solver);

void solverWriteReestart(SOLVER* solver);

void solverLoadRestart(SOLVER* solver, char* fileName);

void solverInitU(SOLVER* solver, CONDITION* inside);

void solverResetR(SOLVER* solver);

double solverCalcP(SOLVER* solver, double** U, int ii);

void solverCalcVel(SOLVER* solver, double** U, int ii, double* u, double* v, double* c);

void rotation(double* U, double dSx, double dSy, double dS);

void solverUpdateGrad(SOLVER* solver);

void solverGrad_T(SOLVER* solver);

void inter(SOLVER* solver);

void interAxisPressure(SOLVER* solver);

void solverCalcR(SOLVER* solver, double** U);

void solverRK(SOLVER* solver, double a);

void solverUpdateU(SOLVER* solver);

void solverStepRK(SOLVER* solver);

void solverCalcRes(SOLVER* solver);

double solverLocalTimeStep(SOLVER* solver, int ii);

void solverCalcDt(SOLVER* solver);

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm);

void solverCalcGrad2(SOLVER* solver, ELEMENT* E, int kk, double* dUx, double* dUy, double* Umin, double* Umax);

void solverCalcGrad3(SOLVER* solver, ELEMENT* E, int kk, double* dUx, double* dUy);

void solverCalcMinMax(SOLVER* solver, ELEMENT* E, int kk, double* Umin, double* Umax);

void solverCheckGrad(SOLVER* solver);

void solverCalcPrimitive(SOLVER* solver, double** U);

int solverTimeSchemeChoice(char* s);

void solverPrintP(SOLVER* solver);

void solverCalcCoeff(SOLVER* solver, double *Fx, double *Fy);

void solverCalcCoeff3(SOLVER* solver, FILE* convFile, int Nint);

void solverSetData(SOLVER* solver, INPUT* input);

SOLVER* solverInit(char* wd);

void solverInitDomain(SOLVER* solver);

void solverSolve(SOLVER* solver);

void solverUpdateUImplicit(SOLVER* solver);

int solverTimeSchemeChoice(char* s);

void solverWriteSurf(SOLVER* solver);

void solverWriteSolution2(SOLVER* solver);

void solverPrintConvReader(SOLVER* solver, FILE* convFile);

void inviscidWriteSurf(SOLVER* solver);

void solverFaceRes(SOLVER* solver);

#endif
