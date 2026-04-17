#ifndef FLUX_H
#define FLUX_H

typedef struct GASPROP GASPROP;

typedef struct SOLVER SOLVER;

typedef struct FLUX
{
    char* type;
    
    int Nvar;
    
    double eFix;
    double Minf;
    
    bool extraVar;
    
    void (*func)(struct FLUX*, GASPROP*, double*, double*, double*);

} FLUX;

FLUX* fluxInit(INPUT* input, SOLVER* solver);

void fluxFree1(FLUX* flux);

void fluxFuncRoe(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f);
	           
void fluxEntropyFix(FLUX* flux, double *l);	           
	           
void fluxFuncAUSM(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f);

void fluxFuncAUSMDV(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f);

void fluxFuncAUSMpup(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f);

void fluxFuncAUSMpup2(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f);

#endif
