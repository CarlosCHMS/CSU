#ifndef SST_H
#define SST_H

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


typedef struct SST{ 

    bool active;
  
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
    
    SST_TRANS* trans;
} SST;

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


SST* sstInit(INPUT* input);

void sstMalloc(SST* sst, int Nelem);

void sstFree(SST* sst);

void sstInitU(SOLVER* solver, CONDITION* inside);

void sstInterFace(SOLVER* solver);

void sstInter(SOLVER* solver);

double sstBlend(double x1, double x2, double F1);

void sstInterSource(SOLVER* solver);

void sstSources(SST* sst, SSTVAR* var);

void sstFlux(SST* sst, SSTVAR* var);

double sstF1(SST* sst, SSTVAR* var, double n_L_term, double sqrtk_term);

double sstF2(SST* sst, SSTVAR* var, double n_L_term, double sqrtk_term);

void sstSolverWriteSurf(SOLVER* solver);

void sstInterMiT(SOLVER* solver);

void sstBoundaryViscousFluxSymmetry(BOUNDARY* boundary, SOLVER* solver, int ii, double* f, double* miEddy);

void sstBoundaryViscousFluxGeneral(BOUNDARY* boundary, SOLVER* solver, int ii, double* f, double* miEddy);

void sstBoundaryViscousFluxWall(BOUNDARY* boundary, SOLVER* solver, int ii, double* f, double* miEddy);

#endif
