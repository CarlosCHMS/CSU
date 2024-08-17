
typedef struct E1 ELEMENT1;

typedef struct F1{
    
    bool b;
    double* flux;
    double* dS;
    double* c;

    ELEMENT1** eL;

} FACE1;

typedef struct E1{

    int Nf;

    double omega;
    
    bool* nei;

    double* P;    
    double* U;
    double* Uaux;
    double* R;
    double* c;

    double** grad;    
    
    FACE1** fL;

} ELEMENT1;

typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[5];  
    double Pin[6];     

} CONDITION;

typedef struct{

    int Ne;
    int Nf;
    
    char name[50];

    int flagBC;

    ELEMENT1** eL;
    
    FACE1** fL;


} SOLVERBC;

typedef struct {

    int Ndim;
    int Nvar;
    int Ne;
    int Nf;
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
    int Nbc;
    int axi;
    int mallocType;

    FILE* fb;

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
    double *dSLateral;
    double *dist;
            
    double **U;
    double **R;
    double **Uaux;     
    double **faceFlux;       
    double **dPx;
    double **dPy;  
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    INPUT* input;
    
    ELEMENT1** eL;
    
    FACE1** fL;
    
    SOLVERBC** bcL;

} SOLVER;

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny);

void conditionState(CONDITION* cond, SOLVER* solver);

void solverMalloc(SOLVER* solver);

void solverWriteSolution(SOLVER* solver);

void solverResetR(SOLVER* solver);

void rotation(double* U, double dSx, double dSy, double dS);

void solverCalcRes(SOLVER* solver);

double limiterBJ(double Ui, double Umin, double Umax, double d2);

double limiterV(double Ui, double Umin, double Umax, double d2, double e);

double sutherland(double T);

void solverSetData(SOLVER* solver, INPUT* input);

SOLVER* solverInit(char* wd);

void solverSolve(SOLVER* solver);

void newSolverMalloc(SOLVER* solver);

void newSolverMalloc2(SOLVER* solver);

void newSolverFree(SOLVER* solver);

void newSolverGetMeshData(SOLVER* solver);

FACE1* newSolverInitFace(SOLVER* solver);

void newSolverInitDomain(SOLVER* solver);

void newSolverInitU(SOLVER* solver, CONDITION* inside);

double newSolverLocalTimeStep(SOLVER* solver, ELEMENT1* e);

void newSolverCalcVel(SOLVER* solver, ELEMENT1* e, double* u, double* v, double* c);

void newSolverCheckOrientation(SOLVER* solver);

void newBoundaryCalcPrimitive(SOLVER* solver, SOLVERBC* bc);

void newSolverUpdatePrimitives(SOLVER* solver);

double newSolverCalcP(SOLVER* solver, ELEMENT1* e);

void newSolverInter(SOLVER* solver);

double solverModule(double* v);

void newSolverInterFace(SOLVER* solver);

void newSolverInterFace_sa(SOLVER* solver);

void newSolverCalcR(SOLVER* solver);

void newSolverStep(SOLVER* solver);

void newSolverSolve(SOLVER* solver);

void newSolverEuler(SOLVER* solver);

void newSolverUaux(SOLVER* solver);

void checkNei(SOLVER* solver);

void newSolverCalcDt(SOLVER* solver);

void newBoundaryCalc(SOLVER* solver, SOLVERBC* bc);

void solverWriteR(SOLVER* solver);

void solverWriteFlux(SOLVER* solver);

void solverWriteFluxBC(SOLVER* solver, double* f);

void solverClassifyFace(SOLVER* solver);

void newSolverRK(SOLVER* solver, double a);

void newSolverStepRK(SOLVER* solver);

void newSolverCalcCoeff3(SOLVER* solver, FILE* convFile, int Nint);

void newSolverCalcGrad2(ELEMENT1* E, int kk, int grad_ii, double* Umin, double* Umax);

void newSolverReconstruct(SOLVER* solver);

void solverWriteGrad(SOLVER* solver);

void newSolverLoadRestart(SOLVER* solver, char* fileName);

void newSolverWriteRestart(SOLVER* solver);

void newSolverWriteElemBoundary(SOLVER* solver);

void newSolverWriteBCP(SOLVER* solver);

void solverWriteU(SOLVER* solver);

void solverInterAxisPressure(SOLVER* solver);

void newSolverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm);

void newSolverInterVisc(SOLVER* solver);

void newSolverDistributeFlux(SOLVER* solver);

void newSolverBoundaryCalcVisc2(SOLVER* solver, SOLVERBC* bc);

void newSolverBoundaryCalcVisc(SOLVER* solver, SOLVERBC* bc);

void newBoundaryCalcFrictionWall(SOLVER* solver, FACE1* f, double* fx, double* fy);

void newSolverElementFree(SOLVER* solver, ELEMENT1* e, bool internal);

void newSolverFaceFree(SOLVER* solver, FACE1* f);
