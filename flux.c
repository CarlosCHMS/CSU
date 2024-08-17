#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"

int fluxChoice(char* s)
{
    int ans;

    if(strcmp(s, "ROE") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "AUSMD") == 0)
    {
        ans = 1;
    }
    else if(strcmp(s, "AUSMDV") == 0)
    {
        ans = 2;
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
	double U3L = pL/(solver->gamma - 1) + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

	double U0R = rR;
	double U1R = rR*uR;	
	double U2R = rR*vR;	
	double U3R = pR/(solver->gamma - 1) + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;

    // Mean values calculation
	double rqL = sqrt(rL);
    double rqR = sqrt(rR);

	double ub = (rqL*uL + rqR*uR)/(rqL + rqR);
	double vb = (rqL*vL + rqR*vR)/(rqL + rqR);
	double Hb  = (rqL*HL + rqR*HR)/(rqL + rqR);
	double ab  = sqrt((solver->gamma-1) * (Hb - (ub*ub + vb*vb)/2));

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
    a4 /= (ab*ab)/(solver->gamma-1);    
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
	double U3L = pL/(solver->gamma - 1) + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

	double U3R = pR/(solver->gamma - 1) + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;

    double cL = sqrt(solver->gamma*pL/rL);
    double cR = sqrt(solver->gamma*pR/rR);
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
	double U3L = pL/(solver->gamma - 1) + (uL*uL + vL*vL)*rL/2;
    double HL = (U3L + pL)/rL;

	double U3R = pR/(solver->gamma - 1) + (uR*uR + vR*vR)*rR/2;
    double HR = (U3R + pR)/rR;

    double cL = sqrt(solver->gamma*pL/rL);
    double cR = sqrt(solver->gamma*pR/rR);
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


void flux(SOLVER* solver, double rL, double uL, double vL, double pL,
                              double rR, double uR, double vR, double pR, double* f)
{
	if(solver->flux == 0)
	{
        fluxRoe(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
    else if(solver->flux == 2)     
    {
        fluxAUSMDV(solver, rL, uL, vL, pL, rR, uR, vR, pR, f);
    }
    
}	


void fluxFree(SOLVER* solver, double rL, double uL, double vL, double pL, double* f)
{
    f[0] = rL*uL;   
    f[1] = rL*uL*uL + pL;
    f[2] = rL*uL*vL;
    f[3] = uL*(solver->gamma*pL/(solver->gamma-1) + 0.5*(uL*uL + vL*vL)*rL);    
}
