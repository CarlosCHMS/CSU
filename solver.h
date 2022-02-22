
typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[4];  

} CONDITION;

typedef struct{

    char* down;
    char* up;
    char* left;
    char* right;

    int Ndown;
    int Nup;
    int Nleft;
    int Nright;

} BOUNDARY;


typedef struct {

    int Nrow;
    int Ncol;
    int pOutFlag;
    int MUSCL;
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
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    BOUNDARY* bc;

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
