
typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[4];  
    double Pin[5];     

} CONDITION;

typedef struct {

    int Nrow;
    int Ncol;
    int pOutFlag;
    int order;
    int flux;
    int stages;
    int laminar;

    double Rgas;
    double gamma;
    double k4;
    double dt;
    double pout; 
    double eFix;       
    double e; 
    double k; 
    double res0[4];
    double res[4];
    double CFL;
    double Cp;
    double Pr;
            
    double **U;
    double **R;
    double **Uaux;     
    double **faceFlux;       
    double **dPx;
    double **dPy;    
    double **P;    
    
    CONDITION* inlet;
        
    MESH* mesh;

} SOLVER;

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny);

void conditionState(CONDITION* cond, SOLVER* solver);

double conditionVref(CONDITION* cond, SOLVER* solver);

void solverMalloc(SOLVER* solver);

void solverFree(SOLVER* solver);

void solverWrite(SOLVER* solver, char* fileName);

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

double solverCalcDt(SOLVER* solver);

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm);

void solverCalcGrad2(SOLVER* solver, double* U, int ii, double* dUx, double* dUy, double* Umin, double* Umax);

void solverCheckGrad(SOLVER* solver);

double limiterBJ(double Ui, double Umin, double Umax, double d2);

double limiterV(double Ui, double Umin, double Umax, double d2, double e);

void solverCalcPrimitive(SOLVER* solver, double** U);

void solverCalcUfromP(SOLVER* solver, double r, double u, double v, double p, double* U0, double* U1, double* U2, double* U3);

void solverMallocP(SOLVER* solver);

void solverFreeP(SOLVER* solver);

double sutherland(double T);
