#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<stdbool.h>
#include<omp.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"gasprop.h"

GASPROP* gaspropInit(double gamma, double R, int TP)
{

    GASPROP* gas = malloc(sizeof(GASPROP));
    gas->gamma = gamma;
    gas->R = R;
    gas->Cp = gas->gamma*gas->R/(gas->gamma - 1);
    gas->Cv = gas->Cp - gas->R;

    gas->N = 6;
    gas->TP = TP;

    gas->cc = malloc(gas->N*sizeof(double));

    gas->cc[0] = 1.81716e3;
    gas->cc[1] = 9.890504e2 - 287.0530;
    gas->cc[2] = -9.595592e-3;
    gas->cc[3] = 1.041469e-4;
    gas->cc[4] = -4.433065e-8;
    gas->cc[5] = 5.879263e-12;
    
    return gas;

}

void gaspropFree(GASPROP* gas)
{

    free(gas->cc);
    free(gas);

}

double gasprop_T2e(GASPROP* gas, double T)
{
    double aux;
    double P;

    if(gas->TP)
    {
        P = 1;
        aux = gas->cc[0];
        for(int ii=1; ii<gas->N; ii++)
        {
            P *= T;
            aux += gas->cc[ii]*P;
        }
    }
    else
    {
        aux = gas->Cv*T;
    }   

    return aux;

}

double gasprop_T2Cv(GASPROP* gas, double T)
{

    double aux;
    double P;

    if(gas->TP)
    {
        P = 1;
        aux = gas->cc[1];
        for(int ii=2; ii<gas->N; ii++)
        {
            P *= T;
            aux += ii*gas->cc[ii]*P;
        }
    }
    else
    {
        aux = gas->Cv;
    }   

    return aux;

}

void gasprop_T2eCv(GASPROP* gas, double T, double* e, double* Cv)
{
    double P;

    P = 1;
    *e = gas->cc[0];
    *Cv = gas->cc[1];
    for(int ii=1; ii<gas->N; ii++)
    {
        P *= T;
        *e += gas->cc[ii]*P;
        if(ii + 1 < gas->N)
        {
            *Cv += (ii+1)*gas->cc[ii+1]*P;
        }
    }
}

double gasprop_e2T(GASPROP* gas, double e)
{

    double T = e/gas->Cv;

    if(gas->TP)
    {
        double error = 1;
        double e1, Cv;
        int ii = 0;
        while((fabs(error) > 1e-14) & (ii<8))
        {
            gasprop_T2eCv(gas, T, &e1, &Cv);
            error = (e - e1)/Cv;
            T += error;
            ii++;
            //printf("%e\n", error);
        }
    }

    return T;

}

double gasprop_T2c(GASPROP* gas, double T)
{
    return sqrt(gasprop_T2gamma(gas, T)*gas->R*T);
}


double gasprop_T2Cp(GASPROP* gas, double T)
{
    return gasprop_T2Cv(gas, T) + gas->R;
}


double gasprop_T2gamma(GASPROP* gas, double T)
{
    double Cv = gasprop_T2Cv(gas, T);
    return (Cv + gas->R)/Cv;
}

double gasprop_e2Taprox(GASPROP* gas, double e0, double T0, double e1)
{
    return T0 + (e1 - e0)/gasprop_T2Cv(gas, T0);
}

double gasprop_T2fentro(GASPROP* gas, double T)
{
    double aux;
    double P;

    if(gas->TP)
    {
        P = 1;
        aux = (gas->cc[1] + gas->R)*log(T);
        for(int ii=2; ii<gas->N; ii++)
        {
            P *= T;
            aux += ii*gas->cc[ii]*P/(ii-1);
        }
    }
    else
    {
        aux = gas->Cp*log(T);
    }   

    return aux;
}

double gasprop_Tp2entropy(GASPROP* gas, double T, double p)
{
    return gasprop_T2fentro(gas, T) - gas->R*log(p);
}

double gasprop_critic_T2entalpy(GASPROP* gas, double T)
{
    double e, Cv;
    gasprop_T2eCv(gas, T, &e, &Cv);

    double RT = gas->R*T;
    double gamma = (Cv + gas->R)/Cv;
    return e + RT + (gamma*RT)*0.5;

}

void gasprop_critic_T2entalpyDer(GASPROP* gas, double T, double *H, double *dH)
{
    double e, Cv;
    gasprop_T2eCv(gas, T, &e, &Cv);

    double RT = gas->R*T;
    double gamma = (Cv + gas->R)/Cv;
    *H = e + RT + gamma*RT*0.5;
    *dH = Cv + gas->R + gamma*gas->R*0.5; //Aproximated 
}


double gasprop_critic_H2c(GASPROP* gas, double H)
{

    double c2 = 2*(gas->gamma-1)*H/(gas->gamma+1);
    double ans;
    double H0, T, dHdT;
    double error = 1;
    int ii = 0;

    if(gas->TP)
    {
        T = c2/(gas->gamma*gas->R);
        gasprop_critic_T2entalpyDer(gas, T, &H0, &dHdT);
        error = (H - H0)/dHdT;
        T += error;
        
        while((fabs(error) > 1.0e-12) & (ii < 3))
        {        
            gasprop_critic_T2entalpyDer(gas, T, &H0, &dHdT);
            error = (H - H0)/dHdT;
            T += error;
            ii++;
        }
      
        ans = gasprop_T2c(gas, T);
    }
    else
    {
        ans = sqrt(c2);
    }
    
    return ans;
}
