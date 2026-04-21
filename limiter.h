#ifndef LIMITER_H
#define LIMITER_H

typedef struct INPUT INPUT;

typedef struct LIMITER{

    char* name;

    int type;
    double K;
    double volMax;
    
    double* Pref20;
    double* Pref;

    double** phi;

    double (*func)(struct LIMITER*, double, double, double, double, double);

} LIMITER;

LIMITER* limiterInit(INPUT* input, SOLVER* solver);

void limiterFree(LIMITER* limiter, SOLVER* solver);

void limiterUpdate(LIMITER* limiter, SOLVER* solver);

void limiterCalc(LIMITER* limiter, SOLVER* solver, int ii, double* Pref2);

double limiterFuncV2(LIMITER* limiter, double Ui, double Umin, double Umax, double d2, double ee);

double limiterFuncBJ(LIMITER* limiter, double Ui, double Umin, double Umax, double d2, double ee);

#endif
