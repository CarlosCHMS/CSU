
typedef struct{

    int ii;

    double dux;
    double duy;
    double dvx;
    double dvy;    
    double dkx;
    double dky;
    double dox;
    double doy;
    double r;
    double T;
    double k;
    double om;
    double d;
    double mi_L;
    double aux;
    double l;
    double omMin;
    double x;
    double y;
    
    double F1;
    double F2;
    
    double Qtk;
    double Qto;

    double dQkdr;
    double dQkdrk;    
    double dQkdro;
    
    double dQodr;
    double dQodrk;    
    double dQodro;
    
    double mi_t;
    double tkx;
    double tky;
    double tox;
    double toy;

} SSTVAR;


SST* sstInit();

void sstFree(SST* sst);

void sstInitU(SOLVER* solver, CONDITION* inside);

void sstInterFace(SOLVER* solver);

void sstInter(SOLVER* solver);

void sstBoundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f, double* miEddy);

void sstBoundaryFace(SOLVER* solver, MESHBC* bc);

void sstBoundary(SOLVER* solver);

double sstBlend(double x1, double x2, double F1);

void sstInterSource(SOLVER* solver);

void sstSources(SST* sst, SSTVAR* var);

void sstFlux(SST* sst, SSTVAR* var);

double sstF1(SST* sst, SSTVAR* var, double n_L_term, double sqrtk_term);

double sstF2(SST* sst, SSTVAR* var, double n_L_term, double sqrtk_term);

void sstSolverWriteSurf(SOLVER* solver);

void sstInterMiT(SOLVER* solver);

