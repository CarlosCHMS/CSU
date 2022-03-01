
typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[4];  

} CONDITION;

typedef struct {

    int Nrow;
    int Ncol;
    int pOutFlag;
    int order;
    int flux;
    int stages;

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
            
    double **U;
    double **R;
    double **Uaux;     
    double **faceFlux;       
    double **dUx;
    double **dUy;    
    
    CONDITION* inlet;
        
    MESH* mesh;

} SOLVER;

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny);

void conditionState(CONDITION* cond, SOLVER* solver);

double conditionVref(CONDITION* cond, SOLVER* solver);

void solverFree(SOLVER* solver);

void solverWrite(SOLVER* solver, char* fileName);

void solverInitU(SOLVER* solver, CONDITION* inside);

void solverResetR(SOLVER* solver);

double solverCalcP(SOLVER* solver, double** U, int ii);

void solverCalcVel(SOLVER* solver, double** U, int ii, double* u, double* v, double* c);

void rotation(double* U, double dSx, double dSy, double dS);

void inter(SOLVER* solver, double **U);

void boundaryInlet(SOLVER* solver, double* Ua, double* Ud, double* Ub, double nx, double ny);

void boundaryOutlet(SOLVER* solver, double* Ud, double* Ub, double nx, double ny);

void boundaryWall(SOLVER* solver, double* Ud, double* Ub, double nx, double ny);

void boundaryCalc(SOLVER* solver, double **U, MESHBC* bc);

void boundary(SOLVER* solver, double **U);

void interAxisPressure(SOLVER* solver, double **U);

void solverCalcR(SOLVER* solver, double** U);

void solverRK(SOLVER* solver, double a);

void solverUpdateU(SOLVER* solver);

void solverStepRK(SOLVER* solver);

void solverCalcRes(SOLVER* solver);

int boundaryChoice(char* s);

void boundaryGetBC(MESH* mesh, INPUT* input);

double solverLocalTimeStep(SOLVER* solver, int ii);

double solverCalcDt(SOLVER* solver);

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm);

void solverCalcGrad(SOLVER* solver, double* U, int ii, double* dUx, double* dUy);

void solverCalcGrad2(SOLVER* solver, double* U, int ii, double* dUx, double* dUy, double* Umin, double* Umax);

void solverCheckGrad(SOLVER* solver);

double limiterBJ(double Ui, double Umin, double Umax, double d2);

double limiterV(double Ui, double Umin, double Umax, double d2, double e);
