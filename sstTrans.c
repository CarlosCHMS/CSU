#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include <stdbool.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"
#include"boundary.h"
#include"sst.h"
#include"sstTrans.h"
#include"gasprop.h"

SST_TRANS* sstTransInit()
{
    SST_TRANS* sstTrans = malloc(sizeof(SST_TRANS));
    sstTrans->sk1 = 0.85;
    sstTrans->so1 = 0.5;
    sstTrans->b1 = 0.075;
    
    sstTrans->sk2 = 1.0;
    sstTrans->so2 = 0.856;
    sstTrans->b2 = 0.0828;
    sstTrans->bs2 = 0.09;
    sstTrans->alphas2 = 1.0;
    sstTrans->alpha2 = 0.4403;
    
    sstTrans->a1 = 0.31;
    
    return sstTrans;
}

void sstTransFree(SST_TRANS* sstTrans)
{
    free(sstTrans);
}

void sstTransFlux(SST_TRANS* trans, SSTVAR* var)
{
    double n_L = var->mi_L/var->r;

    double omega = fabs(var->duy - var->dvx);

    double sqrtk = sqrt(var->k);

    double sqrtk_term = sqrtk/(trans->bs2*var->om*var->d);
    
    double n_L_term = 500*n_L/(var->d*var->d*var->om);

    double Fs1 = sstTransF1(trans, var, n_L_term, sqrtk_term);

    double Ry = 1e6;//var->r*var->d*sqrtk/var->mi_L;

    double aux = Ry/120;
    aux = aux*aux;
    aux = aux*aux;
    aux = aux*aux;

    double F1 = fmax(Fs1, exp(-aux));

    double F2 = sstTransF2(trans, var, n_L_term, sqrtk_term);

    double Rt = 1e6;//var->r*var->k/(var->mi_L*var->om);

    double alphas1 = sstTransAlphas1(Rt);

    double alphas = sstBlend(alphas1, trans->alphas2, F1);    
    double sk = sstBlend(trans->sk1, trans->sk2, F1);
    double so = sstBlend(trans->so1, trans->so2, F1);

    var->mi_t = fmin(alphas*var->r*var->k/var->om, trans->a1*var->r*var->k/(omega*F2));

    aux = (var->mi_L + sk*var->mi_t);
    var->tkx = aux*var->dkx;
    var->tky = aux*var->dky;
    aux = (var->mi_L + so*var->mi_t);
    var->tox = aux*var->dox;
    var->toy = aux*var->doy;
}

double sstTransAlphas1(double Rt)
{
    return (.025 + Rt/6)/(1 + Rt/6);
}

double sstTransBs1(double Rt)
{
    double aux = Rt/8;
    aux = aux*aux;
    aux = aux*aux;

    return 0.09*(5./18. + aux)/(1 + aux);
}

double sstTransF1(SST_TRANS* sst, SSTVAR* var, double n_L_term, double sqrtk_term)
{
    long double aux1, aux2;

    aux1 = 2*var->r*sst->so2*(var->dkx*var->dox + var->dky*var->doy)/var->om;
    double CD = fmax(aux1, 1.0e-20);

    aux1 = sqrtk_term;
    aux2 = n_L_term;
    
    aux1 = fmax(aux1, aux2);
    aux2 = 4*var->r*sst->so2*var->k/(CD*var->d*var->d);
    
    double arg1 = fmin(aux1, aux2);
    double F1 = tanh(arg1*arg1*arg1*arg1);

    return F1;
}

double sstTransF2(SST_TRANS* sst, SSTVAR* var, double n_L_term, double sqrtk_term)
{
    double arg2 = fmax(2*sqrtk_term, n_L_term);

    double F2 = tanh(arg2*arg2);

    return F2;
}

void sstTransSources(SST_TRANS* trans, SSTVAR* var)
{
    double n_L = var->mi_L/var->r;

    double omega = fabs(var->duy - var->dvx);

    double sqrtk = sqrt(var->k);

    double sqrtk_term = sqrtk/(trans->bs2*var->om*var->d);
    
    double n_L_term = 500*n_L/(var->d*var->d*var->om);

    double Fs1 = sstTransF1(trans, var, n_L_term, sqrtk_term);

    double Ry = 1e6;//var->r*var->d*sqrtk/var->mi_L;

    double aux = Ry/120;
    aux = aux*aux;
    aux = aux*aux;
    aux = aux*aux;

    double F1 = fmax(Fs1, exp(-aux));

    double F2 = sstTransF2(trans, var, n_L_term, sqrtk_term);

    double Rt = 1e6;//var->r*var->k/(var->mi_L*var->om);

    double alpha1 = sstTransAlpha1(Rt);
    double alphas1 = sstTransAlphas1(Rt);
    double bs1 = sstTransBs1(Rt);    

    double alpha = sstBlend(alpha1, trans->alpha2, F1);    
    double alphas = sstBlend(alphas1, trans->alphas2, F1);    
    double b = sstBlend(trans->b1, trans->b2, F1);    
    double bs = sstBlend(bs1, trans->bs2, F1);

    var->mi_t = fmin(alphas*var->r*var->k/var->om, trans->a1*var->r*var->k/(omega*F2));
    
    double n_t = var->mi_t/var->r;
    
    double P = var->mi_t*omega*omega;
        
    var->F1 = F1;
    var->F2 = F2;
        
    var->Qtk = P - bs*var->r*var->om*var->k;
    var->Qto = alpha*P/n_t - b*var->r*var->om*var->om + 2*(1-F1)*var->r*trans->so2*(var->dkx*var->dox + var->dky*var->doy)/var->om;
}

double sstTransAlpha1(double Rt)
{
    return (5./9.)*(0.1 + Rt/2.7)/(1  + Rt/2.7);
}
