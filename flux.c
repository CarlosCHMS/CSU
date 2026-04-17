#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include <stdbool.h>
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"
#include"gasprop.h"


FLUX* fluxInit(INPUT* input, SOLVER* solver)
{

    FLUX* flux = malloc(sizeof(FLUX));
    flux->type = malloc(10*sizeof(char));
    flux->type[0] = '\0';
    strcat(flux->type, inputGetValue(input, "flux"));
    
    flux->Nvar = solver->Nvar;
    
    flux->Minf = solver->inlet->mach;
    
    flux->eFix = 0.1;
    
    flux->extraVar = false;
    if(flux->Nvar > 4)
    {
        flux->extraVar = true;
    }
    
    if(strcmp(flux->type, "ROE") == 0)
    {        
        flux->func = fluxFuncRoe;
    }
    else if(strcmp(flux->type, "AUSM") == 0)
    {
        flux->func = fluxFuncAUSM;
    }
    else if(strcmp(flux->type, "AUSMDV") == 0)
    {
        flux->func = fluxFuncAUSMDV;
    }
    else if(strcmp(flux->type, "AUSMpup") == 0)
    {
        flux->func = fluxFuncAUSMpup;        
    }
    else if(strcmp(flux->type, "AUSMpup2") == 0)
    {
        flux->func = fluxFuncAUSMpup2;        
    }    
    else
    {
        printf("Error in flux choice: %s.\n", flux->type);
        exit(0);
    }    
    
    return flux;
}


void fluxFree1(FLUX* flux)
{

    free(flux->type);
    free(flux);
}


void fluxFuncRoe(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f)
{

    /*
    Based on: P. L. ROE, Riemann Solvers, Parameter Vectors, and Difference Schemes, (1981)
    */
    
    double rL = PL[0];
    double uL = PL[1];
    double vL = PL[2];
    double pL = PL[3];
    
    double rR = PR[0];
    double uR = PR[1];
    double vR = PR[2];
    double pR = PR[3];
      
	double U0L = rL;
	double U1L = rL*uL;	
	double U2L = rL*vL;	
	double U3L = pL/(gas->gamma - 1) + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

	double U0R = rR;
	double U1R = rR*uR;	
	double U2R = rR*vR;	
	double U3R = pR/(gas->gamma - 1) + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;

    // Mean values calculation
	double rqL = sqrt(rL);
    double rqR = sqrt(rR);

	double ub = (rqL*uL + rqR*uR)/(rqL + rqR);
	double vb = (rqL*vL + rqR*vR)/(rqL + rqR);
	double Hb  = (rqL*HL + rqR*HR)/(rqL + rqR);
	double ab  = sqrt((gas->gamma-1) * (Hb - (ub*ub + vb*vb)/2));

    // Eigenvalues
	double l1 = ub - ab;
	double l2 = ub;
	double l4 = ub;
	double l5 = ub + ab;	

    // Eigenvectors
	double e1[4] = {1.0, ub-ab, vb, Hb-ub*ab};
	double e2v[4] = {0.0, 0.0, 1.0, vb};
	double e4[4] = {1.0, ub, vb, 0.5 * (ub*ub + vb*vb)};
	double e5[4] = {1.0, ub+ab, vb, Hb + ub*ab};

    // Diferences
	double d1 = U0R - U0L;
	double d2 = U1R - U1L;
	double d3 = U2R - U2L;
	double d5 = U3R - U3L;

    // Projections
    double a4 = (Hb - (ub*ub + vb*vb))*d1 + ub*d2 + vb*d3 - d5;
    a4 /= (ab*ab)/(gas->gamma-1);    
    double a2v = d3 - d1*vb;
    double a5 = ((d1 - a4) + (d2 - ub*d1)/ab)*0.5;
    double a1 = ((d1 - a4) - (d2 - ub*d1)/ab)*0.5;

    // Fluxes
	double fL[4] = {U1L, U1L*uL + pL, U1L*vL, uL*(U3L + pL)};
	double fR[4] = {U1R, U1R*uR + pR, U1R*vR, uR*(U3R + pR)};

    // Entropy fix
    fluxEntropyFix(flux, &l1);
    fluxEntropyFix(flux, &l2);
    fluxEntropyFix(flux, &l4);
    fluxEntropyFix(flux, &l5);    
    
	
	for (int ii = 0; ii < 4; ++ii) {
		f[ii] = 0.5 * (fR[ii] + fL[ii] - a1*fabs(l1)*e1[ii] - a2v*fabs(l2)*e2v[ii] - a4*fabs(l4)*e4[ii] - a5*fabs(l5)*e5[ii]);
	}
}


void fluxEntropyFix(FLUX* flux, double *l)
{

    // Harten Hyman entropy fix
    if((*l < flux->eFix) & (*l > -flux->eFix))
    {
        *l = 0.5*(*l * *l/flux->eFix + flux->eFix);
    }

}


void fluxFuncAUSM(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f)
{

    double rL = PL[0];
    double uL = PL[1];
    double vL = PL[2];
    double pL = PL[3];
    
    double rR = PR[0];
    double uR = PR[1];
    double vR = PR[2];
    double pR = PR[3];

    double TL = pL/(gas->R*rL);
	double U3L = gasprop_T2e(gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(gas->R*rR);
	double U3R = gasprop_T2e(gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(gas, TL);
    double cR = gasprop_T2c(gas, TR);

	double ML = uL/cL;
	double MR = uR/cR;

    double Mplus;
    double Mminus;

    double Pplus;
    double Pminus;

    double M2p;
    double M2m;
    
    if(fabs(ML) >= 1)
    {
        Mplus = 0.5*(ML + fabs(ML));
        Pplus = Mplus/ML;
    } 
    else
    {
        M2p = 0.25*(ML + 1)*(ML + 1);
        Mplus = M2p;
        Pplus = M2p*(2 - ML);
    }

    if(fabs(MR) >= 1)
    {
        Mminus = 0.5*(MR - fabs(MR));
        Pminus = Mminus/MR;
    } 
    else
    {
        M2m = -0.25*(MR - 1)*(MR - 1);    
        Mminus = M2m;
        Pminus = -M2m*(2 + MR);
    }

    double Mm = Mplus + Mminus;

    double pm = Pplus*pL + Pminus*pR;

    if(Mm > 0)
    {
	    f[0] = rL*cL*Mm;
	    f[1] = f[0]*uL + pm;
	    f[2] = f[0]*vL;
	    f[3] = f[0]*HL;
	    
	    if(flux->extraVar)
	    {
	        for(int kk=4; kk<flux->Nvar; kk++)
	        {
	            f[kk] = f[0]*PL[kk];
	        }
	    }
    }
    else
    {
	    f[0] = rR*cR*Mm;
	    f[1] = f[0]*uR + pm;
	    f[2] = f[0]*vR;
	    f[3] = f[0]*HR;
	    
	    if(flux->extraVar)
	    {
	        for(int kk=4; kk<flux->Nvar; kk++)
	        {
	            f[kk] = f[0]*PR[kk];
	        }
	    }
    }	
}


void fluxFuncAUSMDV(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f)
{

    /*
    Based on: YASUHIRO WADA † AND MENG-SING LIOU, A Flux Splitting Scheme 
    With High-Resolution and Robustness for Discontinuities, (1994)
    */
    
    double rL = PL[0];
    double uL = PL[1];
    double vL = PL[2];
    double pL = PL[3];
    
    double rR = PR[0];
    double uR = PR[1];
    double vR = PR[2];
    double pR = PR[3];    
    
    double aux;
	double TL = pL/(gas->R*rL);
	double U3L = gasprop_T2e(gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(gas->R*rR);
	double U3R = gasprop_T2e(gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(gas, TL);
    double cR = gasprop_T2c(gas, TR);
	double cm = fmax(cL, cR);

	double alphaL = (2.0*pL/rL)/(pL/rL + pR/rR);
	double alphaR = (2.0*pR/rR)/(pL/rL + pR/rR);

	double uLPlus, pLPlus;
	aux = 0.5*(uL + fabs(uL));
	if (fabs(uL) < cm) 
	{
		uLPlus = alphaL*(0.25*(uL + cm)*(uL + cm)/cm - aux) + aux;
		pLPlus = 0.25*pL*(uL + cm)*(uL + cm)*(2.0 - uL/cm)/(cm*cm);
	} else {
		uLPlus = aux;
		pLPlus = pL*aux/uL;
	}

	double uRMinus, pRMinus;
	aux = 0.5*(uR - fabs(uR));
	if (fabs(uR) < cm) {
		uRMinus = alphaR*(-0.25*(uR - cm)*(uR - cm)/cm - aux) + aux;
		pRMinus = 0.25*pR*(uR - cm)*(uR - cm)*(2.0 + uR/cm)/(cm*cm);
	} else {
		uRMinus = aux;
		pRMinus = pR*aux/uR;
	}

    double psiL[6];
    
    psiL[0] = 1;
    psiL[1] = uL;
    psiL[2] = vL;
    psiL[3] = HL;
    if(flux->extraVar)
    {
        for(int kk=4; kk<flux->Nvar; kk++)
        {
            psiL[kk] = PL[kk];
        }
    }    

    double psiR[6];
    
    psiR[0] = 1;
    psiR[1] = uR;
    psiR[2] = vR;
    psiR[3] = HR;
    if(flux->extraVar)
    {
        for(int kk=4; kk<flux->Nvar; kk++)
        {
            psiR[kk] = PR[kk];
        }
    }    

	double rU = uLPlus*rL + uRMinus*rR;
	
    for(int kk=0; kk<flux->Nvar; kk++)
    {
        f[kk] = 0.5*(rU * (psiR[kk] + psiL[kk]) - fabs(rU) * (psiR[kk] - psiL[kk]));
    }
	
	double f1AUSMD = f[1];
	double f1AUSMV = uLPlus*rL*uL + uRMinus*rR*uR;
	
	double s = 0.5*fmin(1, 10*fabs(pR - pL)/fmin(pL, pR));

	f[1] = (pLPlus + pRMinus);
	f[1] += (0.5 + s)*f1AUSMV + (0.5 - s)*f1AUSMD;
	
	// entropy fix 
	bool caseA = (uL - cL < 0.0) & (uR - cR > 0.0);
	bool caseB = (uL + cL < 0.0) & (uR + cR > 0.0);
	
	if (caseA & !caseB) {
	    aux = 0.125*((uR - cR) - (uL - cL));
	    for(int kk = 0; kk < flux->Nvar; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}
	else if (!caseA & caseB) {
	    aux = 0.125*((uR + cR) - (uL + cL));
    	for(int kk = 0; kk < flux->Nvar; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}	
}


void fluxFuncAUSMpup(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f)
{

    double rL = PL[0];
    double uL = PL[1];
    double vL = PL[2];
    double pL = PL[3];
    
    double rR = PR[0];
    double uR = PR[1];
    double vR = PR[2];
    double pR = PR[3];    

	double TL = pL/(gas->R*rL);
	double U3L = gasprop_T2e(gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(gas->R*rR);
	double U3R = gasprop_T2e(gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double astar;
    
    astar = gasprop_critic_H2c(gas, HL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(gas, HR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double Ku = 0.75;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = flux->Minf;
    double M0 = sqrt(fmin(1, fmax(Mbar*Mbar, Minf*Minf)));

    double fa = M0*(2 - M0);
    if(fa < 1e-3)
    {
        fa = 1e-3;
    }
    double alpha = (3.0/16.0)*(-4 + 5*fa*fa);
    
    double rhom = 0.5*(rL + rR);

    double Mplus;
    double Mminus;

    double Pplus;
    double Pminus;

    double M2p;
    double M2m;
    
    if(fabs(ML) >= 1)
    {
        Mplus = 0.5*(ML + fabs(ML));
        Pplus = Mplus/ML;
    } 
    else
    {
        M2p = 0.25*(ML + 1)*(ML + 1);
        M2m = -0.25*(ML - 1)*(ML - 1);
        Mplus = M2p*(1 - 16*beta*M2m);
        Pplus = M2p*((2 - ML) - 16*alpha*ML*M2m);
    }

    if(fabs(MR) >= 1)
    {
        Mminus = 0.5*(MR - fabs(MR));
        Pminus = Mminus/MR;
    } 
    else
    {
        M2p = 0.25*(MR + 1)*(MR + 1);
        M2m = -0.25*(MR - 1)*(MR - 1);    
        Mminus = M2m*(1 + 16*beta*M2p);
        Pminus = M2m*((-2 - MR) + 16*alpha*MR*M2p);
    }

    double Mm = Mplus + Mminus - (Kp/fa)*fmax(1 - sig*Mbar*Mbar, 0)*(pR - pL)/(rhom*am*am);

    double pm = Pplus*pL + Pminus*pR - Ku*Pplus*Pminus*(rL + rR)*(fa*am)*(uR-uL);

    double mm;
    if(Mm > 0)
    {   
        mm = am*Mm*rL;
    }
    else
    {
        mm = am*Mm*rR;        
    }

    if(mm > 0)
    {
	    f[0] = mm;
	    f[1] = mm*uL + pm;
	    f[2] = mm*vL;
	    f[3] = mm*HL;
	    
        if(flux->extraVar)
	    {
	        for(int kk=4; kk<flux->Nvar; kk++)
	        {
	            f[kk] = mm*PL[kk];
	        }
	    }
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
	    
	    if(flux->extraVar)
	    {
	        for(int kk=4; kk<flux->Nvar; kk++)
	        {
	            f[kk] = mm*PR[kk];
	        }
	    }
    }
}


void fluxFuncAUSMpup2(FLUX* flux, GASPROP* gas, double* PL, double* PR, double* f)
{

    double rL = PL[0];
    double uL = PL[1];
    double vL = PL[2];
    double pL = PL[3];
    
    double rR = PR[0];
    double uR = PR[1];
    double vR = PR[2];
    double pR = PR[3]; 
    
	double TL = pL/(gas->R*rL);
	double U3L = gasprop_T2e(gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(gas->R*rR);
	double U3R = gasprop_T2e(gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double astar;
    
    astar = gasprop_critic_H2c(gas, HL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(gas, HR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = flux->Minf;
    double M0 = sqrt(fmin(1, fmax(Mbar*Mbar, Minf*Minf)));

    double fa = M0*(2 - M0);
    if(fa < 1e-3)
    {
        fa = 1e-3;
    }
    double alpha = (3.0/16.0)*(-4 + 5*fa*fa);
    
    double rhom = 0.5*(rL + rR);

    double Mplus;
    double Mminus;

    double Pplus;
    double Pminus;

    double M2p;
    double M2m;
    
    if(fabs(ML) >= 1)
    {
        Mplus = 0.5*(ML + fabs(ML));
        Pplus = Mplus/ML;
    } 
    else
    {
        M2p = 0.25*(ML + 1)*(ML + 1);
        M2m = -0.25*(ML - 1)*(ML - 1);
        Mplus = M2p*(1 - 16*beta*M2m);
        Pplus = M2p*((2 - ML) - 16*alpha*ML*M2m);
    }

    if(fabs(MR) >= 1)
    {
        Mminus = 0.5*(MR - fabs(MR));
        Pminus = Mminus/MR;
    } 
    else
    {
        M2p = 0.25*(MR + 1)*(MR + 1);
        M2m = -0.25*(MR - 1)*(MR - 1);    
        Mminus = M2m*(1 + 16*beta*M2p);
        Pminus = M2m*((-2 - MR) + 16*alpha*MR*M2p);
    }

    double Mm = Mplus + Mminus - (Kp/fa)*fmax(1 - sig*Mbar*Mbar, 0)*(pR - pL)/(rhom*am*am);

    double pm = 0.5*(pL + pR) + 0.5*(Pplus - Pminus)*(pL - pR) + sqrt(0.5*(uL*uL + vL*vL + uR*uR + vR*vR))*(Pplus + Pminus - 1)*0.5*(pL + pR)/am;

    double mm;
    if(Mm > 0)
    {   
        mm = am*Mm*rL;
    }
    else
    {
        mm = am*Mm*rR;        
    }

    if(mm > 0)
    {
	    f[0] = mm;
	    f[1] = mm*uL + pm;
	    f[2] = mm*vL;
	    f[3] = mm*HL;
	    
        if(flux->extraVar)
	    {
	        for(int kk=4; kk<flux->Nvar; kk++)
	        {
	            f[kk] = mm*PL[kk];
	        }
	    }	    
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
	    
	    if(flux->extraVar)
	    {
	        for(int kk=4; kk<flux->Nvar; kk++)
	        {
	            f[kk] = mm*PR[kk];
	        }
	    }
    }	
}


