
typedef struct BAUX{

    int ii;
    double** A;
    struct BAUX* next;

} BLOCK;

typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[6];  
    double Pin[7];     

} CONDITION;

typedef struct{

    double gamma;
    double R;
    double Cp;
    double Cv;

    double* cc;
    int N;
    int TP;

} GASPROP;

typedef struct{

    int type;
    double K;
    double volMax;
    
    double* Pref20;
    double* Pref2;

} LIMITER;


typedef struct{
    int flag;

    double sk1;
    double so1;
    double b1;
    
    double sk2;
    double so2;
    double b2;   
    double bs2;
    double alphas2;
    double alpha2;
    
    double a1;
} SST_TRANS;


typedef struct{   
    double L;
    double kFactor;
    double oFactor;
    double oWallFactor;
    double g1;
    double g2;
    double sk1;
    double so1;
    double b1;   
    double sk2;
    double so2;
    double b2;   
    double bs;
    double a1;
    
    SST_TRANS* trans;
} SST;


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
    int saCC;    
    int dtLocal;
    int turb1order;
    int timeScheme;
    int Nlinear;
    int viscBlazek;
    int sstFlag;

    char* wd;
    char writeSurf[50];

    double k4;
    double dt;
    double pout;
    double turbRatio; 
    double eFix; 
    double k; 
    double res[6];
    double CFL;
    double Pr;
    double Pr_t;
    double Sref;
    double dtLocalN;
    double K3;
    double wImp;
    double tol;
    double rLim;
    double pLim; 
    double Twall;
    double omLim;       

    double *dtL;
    double *D;
    double *miT;
    double *miTe;
    double *F1;
    double *F2;  
    double *dd;
    double *om2; 
    double *dQodr;     
    double *dQodrk;     
    double *dQodro; 
    double *dQkdr;     
    double *dQkdrk;     
    double *dQkdro;
            
    double **U;
    double **R;
    double **Uaux;     
    double **faceFlux;       
    double **dPx;
    double **dPy;  
    double **phi;
    double **dW0;
    double **dW1;
    double **uD;
    double **rD;    
    
    BLOCK** BB;
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    INPUT* input;

    GASPROP* gas;
    
    LIMITER* limiter;
    
    SST* sst;

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

double limiterBJ(double Ui, double Umin, double Umax, double d2);

void solverCalcPrimitive(SOLVER* solver, double** U);

double sutherland(double T);

int solverTimeSchemeChoice(char* s);

void solverPrintP(SOLVER* solver);

void solverCalcCoeff(SOLVER* solver, double *Fx, double *Fy);

void solverCalcCoeff2(SOLVER* solver, char* path);

void solverCalcCoeff3(SOLVER* solver, FILE* convFile, int Nint);

void solverSetData(SOLVER* solver, INPUT* input);

SOLVER* solverInit(char* wd);

void solverInitDomain(SOLVER* solver);

void solverSolve(SOLVER* solver);

void solverUpdateUImplicit(SOLVER* solver);

int solverTimeSchemeChoice(char* s);

void solverWriteSurf(SOLVER* solver);

void solverWriteSolution2(SOLVER* solver);

