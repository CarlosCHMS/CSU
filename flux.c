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

int fluxChoice(char* s)
{
    int ans;

    if(strcmp(s, "ROE") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "AUSM") == 0)
    {
        ans = 1;
    }
    else if(strcmp(s, "AUSMDV") == 0)
    {
        ans = 2;
    }
    else if(strcmp(s, "AUSMpup") == 0)
    {
        ans = 3;
    }
    else if(strcmp(s, "AUSMpup2") == 0)
    {
        ans = 4;
    }    
    else
    {
        printf("Error in flux choice: %s.\n", s);
        exit(0);
    }
    
    return ans;

}

void entropyFix(SOLVER* solver, double *l)
{

    // Harten Hyman entropy fix
    if((*l < solver->eFix) & (*l > -solver->eFix))
    {
        *l = 0.5*(*l * *l/solver->eFix + solver->eFix);
    }

}

void fluxRoe(SOLVER* solver, 
               double rL, double uL, double vL, double pL, 
               double rR, double uR, double vR, double pR,
	           double* f)
{

    /*
    Based on: P. L. ROE, Riemann Solvers, Parameter Vectors, and Difference Schemes, (1981)
    */
      
	double U0L = rL;
	double U1L = rL*uL;	
	double U2L = rL*vL;	
	double U3L = pL/(solver->gas->gamma - 1) + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

	double U0R = rR;
	double U1R = rR*uR;	
	double U2R = rR*vR;	
	double U3R = pR/(solver->gas->gamma - 1) + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;

    // Mean values calculation
	double rqL = sqrt(rL);
    double rqR = sqrt(rR);

	double ub = (rqL*uL + rqR*uR)/(rqL + rqR);
	double vb = (rqL*vL + rqR*vR)/(rqL + rqR);
	double Hb  = (rqL*HL + rqR*HR)/(rqL + rqR);
	double ab  = sqrt((solver->gas->gamma-1) * (Hb - (ub*ub + vb*vb)/2));

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
    a4 /= (ab*ab)/(solver->gas->gamma-1);    
    double a2v = d3 - d1*vb;
    double a5 = ((d1 - a4) + (d2 - ub*d1)/ab)*0.5;
    double a1 = ((d1 - a4) - (d2 - ub*d1)/ab)*0.5;

    // Fluxes
	double fL[4] = {U1L, U1L*uL + pL, U1L*vL, uL*(U3L + pL)};
	double fR[4] = {U1R, U1R*uR + pR, U1R*vR, uR*(U3R + pR)};

    // Entropy fix
    entropyFix(solver, &l1);
    entropyFix(solver, &l2);
    entropyFix(solver, &l4);
    entropyFix(solver, &l5);    

	
	for (int ii = 0; ii < 4; ++ii) {
		f[ii] = 0.5 * (fR[ii] + fL[ii] - a1*fabs(l1)*e1[ii] - a2v*fabs(l2)*e2v[ii] - a4*fabs(l4)*e4[ii] - a5*fabs(l5)*e5[ii]);
	}
}


void fluxAUSM(SOLVER* solver, 
               double rL, double uL, double vL, double pL,
               double rR, double uR, double vR, double pR,
	           double* f)
{

    double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(solver->gas, TL);
    double cR = gasprop_T2c(solver->gas, TR);

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
    }
    else
    {
	    f[0] = rR*cR*Mm;
	    f[1] = f[0]*uR + pm;
	    f[2] = f[0]*vR;
	    f[3] = f[0]*HR;
    }	
}


void fluxAUSMDV(SOLVER* solver, 
               double rL, double uL, double vL, double pL, 
               double rR, double uR, double vR, double pR,
	           double* f)
{

    /*
    Based on: YASUHIRO WADA † AND MENG-SING LIOU, A Flux Splitting Scheme 
    With High-Resolution and Robustness for Discontinuities, (1994)
    */
    
    double aux;
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(solver->gas, TL);
    double cR = gasprop_T2c(solver->gas, TR);
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

	double rU = uLPlus*rL + uRMinus*rR;
	f[0] = rU;
	f[1] = (pLPlus + pRMinus);
	f[2] = 0.5*(rU * (vR + vL) - fabs(rU) * (vR - vL));
	f[3] = 0.5*(rU * (HR + HL) - fabs(rU) * (HR - HL));

	double f1AUSMD = 0.5*(rU * (uR + uL) - fabs(rU) * (uR - uL));	
	double f1AUSMV = uLPlus*rL*uL + uRMinus*rR*uR;
	
	double s = 0.5*fmin(1, 10*fabs(pR - pL)/fmin(pL, pR));
	
	f[1] += (0.5 + s)*f1AUSMV + (0.5 - s)*f1AUSMD;
	
	// entropy fix 
	int caseA = (uL - cL < 0.0) & (uR - cR > 0.0);
	int caseB = (uL + cL < 0.0) & (uR + cR > 0.0);
	double psiL[4] = {1.0, uL, vL, HL};
	double psiR[4] = {1.0, uR, vR, HR};
	if (caseA & ~caseB) {
	    aux = 0.125*((uR - cR) - (uL - cL));
	    for(int kk = 0; kk < 4; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}
	else if (~caseA & caseB) {
	    aux = 0.125*((uR + cR) - (uL + cL));
    	for(int kk = 0; kk < 4; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}	
	
}


void fluxAUSMpup(SOLVER* solver, 
               double rL, double uL, double vL, double pL,
               double rR, double uR, double vR, double pR,
	           double* f)
{
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double astar;
    
    astar = gasprop_critic_H2c(solver->gas, HL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(solver->gas, HR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double Ku = 0.75;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = solver->inlet->mach;
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
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
    }
	
}


void fluxAUSMpup2(SOLVER* solver, 
               double rL, double uL, double vL, double pL,
               double rR, double uR, double vR, double pR,
	           double* f)
{
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double astar;
    
    astar = gasprop_critic_H2c(solver->gas, HL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(solver->gas, HR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = solver->inlet->mach;
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
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
    }
	
}



void fluxAUSM_sa(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double nL,
               double rR, double uR, double vR, double pR, double nR,
	           double* f)
{
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(solver->gas, TL);
    double cR = gasprop_T2c(solver->gas, TR);

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
	    f[4] = f[0]*nL;
    }
    else
    {
	    f[0] = rR*cR*Mm;
	    f[1] = f[0]*uR + pm;
	    f[2] = f[0]*vR;
	    f[3] = f[0]*HR;
	    f[4] = f[0]*nR;
    }
	
}


void fluxAUSMDV_sa(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double nL,
               double rR, double uR, double vR, double pR, double nR,
	           double* f)
{

    /*
    Based on: YASUHIRO WADA † AND MENG-SING LIOU, A Flux Splitting Scheme 
    With High-Resolution and Robustness for Discontinuities, (1994)
    */
    
    double aux;
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(solver->gas, TL);
    double cR = gasprop_T2c(solver->gas, TR);
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

	double rU = uLPlus*rL + uRMinus*rR;
	f[0] = rU;
	f[1] = (pLPlus + pRMinus);
	f[2] = 0.5*(rU * (vR + vL) - fabs(rU) * (vR - vL));
	f[3] = 0.5*(rU * (HR + HL) - fabs(rU) * (HR - HL));
	f[4] = 0.5*(rU * (nR + nL) - fabs(rU) * (nR - nL));

	double f1AUSMD = 0.5*(rU * (uR + uL) - fabs(rU) * (uR - uL));	
	double f1AUSMV = uLPlus*rL*uL + uRMinus*rR*uR;
	
	double s = 0.5*fmin(1, 10*fabs(pR - pL)/fmin(pL, pR));
	
	f[1] += (0.5 + s)*f1AUSMV + (0.5 - s)*f1AUSMD;
	
	// entropy fix 
	int caseA = (uL - cL < 0.0) & (uR - cR > 0.0);
	int caseB = (uL + cL < 0.0) & (uR + cR > 0.0);
	double psiL[5] = {1.0, uL, vL, HL, nL};
	double psiR[5] = {1.0, uR, vR, HR, nR};
	if (caseA & ~caseB) {
	    aux = 0.125*((uR - cR) - (uL - cL));
	    for(int kk = 0; kk < 5; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}
	else if (~caseA & caseB) {
	    aux = 0.125*((uR + cR) - (uL + cL));
    	for(int kk = 0; kk < 5; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}	
	
}


void fluxAUSMpup_sa(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double nL,
               double rR, double uR, double vR, double pR, double nR,
	           double* f)
{
    
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    double astar;
    
    astar = gasprop_critic_H2c(solver->gas, HL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(solver->gas, HR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double Ku = 0.75;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = solver->inlet->mach;
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

    double pm = Pplus*pL + Pminus*pR - Ku*Pplus*Pminus*(rL + rR)*(am*fa)*(uR - uL);

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
	    f[4] = mm*nL;    
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
	    f[4] = mm*nR;        
    }
	
}


void fluxAUSMpup2_sa(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double nL,
               double rR, double uR, double vR, double pR, double nR,
	           double* f)
{
    
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;
    
    //Verify this calculation of sound velocity with thermaly perfect gas
    double astar;
    
    astar = gasprop_critic_H2c(solver->gas, HL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(solver->gas, HR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = solver->inlet->mach;
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
	    f[4] = mm*nL;    
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
	    f[4] = mm*nR;        
    }
	
}

void fluxAUSM_sst(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double kL, double omL,
               double rR, double uR, double vR, double pR, double kR, double omR,
	           double* f)
{

    double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2 + kL*rL;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2 + kR*rR;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(solver->gas, TL);
    double cR = gasprop_T2c(solver->gas, TR);

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
	    f[4] = f[0]*kL;
	    f[5] = f[0]*omL;	    
    }
    else
    {
	    f[0] = rR*cR*Mm;
	    f[1] = f[0]*uR + pm;
	    f[2] = f[0]*vR;
	    f[3] = f[0]*HR;
	    f[4] = f[0]*kR;
	    f[5] = f[0]*omR;
    }	
}


void fluxAUSMDV_sst(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double kL, double oL,
               double rR, double uR, double vR, double pR, double kR, double oR,
	           double* f)
{

    /*
    Based on: YASUHIRO WADA † AND MENG-SING LIOU, A Flux Splitting Scheme 
    With High-Resolution and Robustness for Discontinuities, (1994)
    */
    
    double aux;
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2 + kL;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2 + kR;
    double HR = (U3R + pR)/rR;
    
    double cL = gasprop_T2c(solver->gas, TL);
    double cR = gasprop_T2c(solver->gas, TR);
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

	double rU = uLPlus*rL + uRMinus*rR;
	f[0] = rU;
	f[1] = (pLPlus + pRMinus);
	f[2] = 0.5*(rU * (vR + vL) - fabs(rU) * (vR - vL));
	f[3] = 0.5*(rU * (HR + HL) - fabs(rU) * (HR - HL));
	f[4] = 0.5*(rU * (kR + kL) - fabs(rU) * (kR - kL));
	f[5] = 0.5*(rU * (oR + oL) - fabs(rU) * (oR - oL));	

	double f1AUSMD = 0.5*(rU * (uR + uL) - fabs(rU) * (uR - uL));	
	double f1AUSMV = uLPlus*rL*uL + uRMinus*rR*uR;
	
	double s = 0.5*fmin(1, 10*fabs(pR - pL)/fmin(pL, pR));
	
	f[1] += (0.5 + s)*f1AUSMV + (0.5 - s)*f1AUSMD;
	
	// entropy fix 
	int caseA = (uL - cL < 0.0) & (uR - cR > 0.0);
	int caseB = (uL + cL < 0.0) & (uR + cR > 0.0);
	double psiL[6] = {1.0, uL, vL, HL, kL, oL};
	double psiR[6] = {1.0, uR, vR, HR, kR, oR};
	if (caseA & ~caseB) {
	    aux = 0.125*((uR - cR) - (uL - cL));
	    for(int kk = 0; kk < 6; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}
	else if (~caseA & caseB) {
	    aux = 0.125*((uR + cR) - (uL + cL));
    	for(int kk = 0; kk < 6; kk++)
	    {
		    f[kk] -= aux*(rR*psiR[kk] - rL*psiL[kk]);
		}		
	}	
	
}



void fluxAUSMpup_sst(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double kL, double oL, 
               double rR, double uR, double vR, double pR, double kR, double oR,
	           double* f)
{
    
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2 + kL*rL;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2 + kR*rR;
    double HR = (U3R + pR)/rR;
    
    double astar;
    
    astar = gasprop_critic_H2c(solver->gas, HL - kL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(solver->gas, HR - kR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double Ku = 0.75;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = solver->inlet->mach;
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

    double pm = Pplus*pL + Pminus*pR - Ku*Pplus*Pminus*(rL + rR)*(am*fa)*(uR - uL);

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
	    f[4] = mm*kL;    
	    f[5] = mm*oL;
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
	    f[4] = mm*kR;
	    f[5] = mm*oR;	            
    }
	
}



void fluxAUSMpup2_sst(SOLVER* solver, 
               double rL, double uL, double vL, double pL, double kL, double oL,
               double rR, double uR, double vR, double pR, double kR, double oR,
	           double* f)
{
    
	double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2 + kL*rL;
    double HL = (U3L + pL)/rL;

    double TR = pR/(solver->gas->R*rR);
	double U3R = gasprop_T2e(solver->gas, TR)*rR + (uR*uR + vR*vR)*rR/2 + kR*rR;
    double HR = (U3R + pR)/rR;
    
    //Verify this calculation of sound velocity with thermaly perfect gas
    double astar;
    
    astar = gasprop_critic_H2c(solver->gas, HL - kL);
    double ahL = astar*astar/fmax(astar, fabs(uL));

    astar = gasprop_critic_H2c(solver->gas, HR - kR);
    double ahR = astar*astar/fmax(astar, fabs(uR));
    
	double am = fmin(ahL, ahR);

    double Kp = 0.25;
    double sig = 1.0;
    double beta = 1.0/8.0;

	double ML = uL/am;
	double MR = uR/am;

    double Mbar = sqrt((uL*uL + uR*uR)/(2*am*am));
    double Minf = solver->inlet->mach;
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
	    f[4] = mm*kL;    
	    f[5] = mm*oL;
    }
    else
    {
	    f[0] = mm;
	    f[1] = mm*uR + pm;
	    f[2] = mm*vR;
	    f[3] = mm*HR;
	    f[4] = mm*kR;        
	    f[5] = mm*oR;	    
    }
	
}


void flux(SOLVER* solver, double rL, double uL, double vL, double pL,
                              double rR, double uR, double vR, double pR, double* f)
{
	if(solver->flux == 0)
	{
        fluxRoe(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
    else if(solver->flux == 1)     
    {
        fluxAUSM(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
    else if(solver->flux == 2)     
    {
        fluxAUSMDV(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
    else if(solver->flux == 3)
    {
        fluxAUSMpup(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
    else if(solver->flux == 4)
    {
        fluxAUSMpup2(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
}	

void flux_sa(SOLVER* solver, double rL, double uL, double vL, double pL, double nL,
                              double rR, double uR, double vR, double pR, double nR, double* f)
{
    if(solver->flux == 1)     
    {
        fluxAUSM_sa(solver, rL, uL, vL, pL, nL, rR, uR, vR, pR, nR, f);
    }
    if(solver->flux == 2)     
    {
        fluxAUSMDV_sa(solver, rL, uL, vL, pL, nL, rR, uR, vR, pR, nR, f);
    }
    else if(solver->flux == 3)
    {
        fluxAUSMpup_sa(solver, rL, uL, vL, pL, nL, rR, uR, vR, pR, nR, f);
    }
    else if(solver->flux == 4)
    {
        fluxAUSMpup2_sa(solver, rL, uL, vL, pL, nL, rR, uR, vR, pR, nR, f);
    }
    
}	

void flux_sst(SOLVER* solver, double rL, double uL, double vL, double pL, double kL, double oL,
                              double rR, double uR, double vR, double pR, double kR, double oR, double* f)
{
    if(solver->flux == 1)
    {
        fluxAUSM_sst(solver, rL, uL, vL, pL, kL, oL, rR, uR, vR, pR, kR, oR, f);
    }
    else if(solver->flux == 2)
    {
        fluxAUSMDV_sst(solver, rL, uL, vL, pL, kL, oL, rR, uR, vR, pR, kR, oR, f);
    }
    else if(solver->flux == 3)
    {
        fluxAUSMpup_sst(solver, rL, uL, vL, pL, kL, oL, rR, uR, vR, pR, kR, oR, f);
    }
    else if(solver->flux == 4)
    {
        fluxAUSMpup2_sst(solver, rL, uL, vL, pL, kL, oL, rR, uR, vR, pR, kR, oR, f);   
    }
}	


void fluxFree(SOLVER* solver, double rL, double uL, double vL, double pL, double* f)
{

    double TL = pL/(solver->gas->R*rL);
	double U3L = gasprop_T2e(solver->gas, TL)*rL + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

    f[0] = rL*uL;   
    f[1] = rL*uL*uL + pL;
    f[2] = rL*uL*vL;
    f[3] = rL*uL*HL;
}
