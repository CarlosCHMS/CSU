
typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[5];  
    double Pin[6];     

} CONDITION;

typedef struct {

    int Nvar;
    int Nrow;
    int Ncol;
    int pOutFlag;
    int order;
    int flux;
    int stages;
    int laminar;
    int restart;
    int sa;
    int dtLocal;

    char* wd;

    double Rgas;
    double gamma;
    double k4;
    double dt;
    double pout;
    double turbRatio; 
    double eFix;       
    double e; 
    double k; 
    double res[5];
    double CFL;
    double Cp;
    double Pr;
    double Pr_t;
    double Sref;
    double dtLocalN;

    double *dtL;
            
    double **U;
    double **R;
    double **Uaux;     
    double **faceFlux;       
    double **dPx;
    double **dPy;  
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    INPUT* input;

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

void inter(SOLVER* solver);

void interVisc(SOLVER* solver);

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

void solverCheckGrad(SOLVER* solver);

double limiterBJ(double Ui, double Umin, double Umax, double d2);

double limiterV(double Ui, double Umin, double Umax, double d2, double e);

void solverCalcPrimitive(SOLVER* solver, double** U);

void solverCalcUfromP(SOLVER* solver, double r, double u, double v, double p, double* U0, double* U1, double* U2, double* U3);

double sutherland(double T);

void solverPrintP(SOLVER* solver);

void solverCalcCoeff(SOLVER* solver, double *Fx, double *Fy);

void solverCalcCoeff2(SOLVER* solver, char* path);

void solverCalcCoeff3(SOLVER* solver, FILE* convFile, int Nint);

void solverSetData(SOLVER* solver, INPUT* input);

SOLVER* solverInit(char* wd);

void solverInitDomain(SOLVER* solver);

void solverSolve(SOLVER* solver);
