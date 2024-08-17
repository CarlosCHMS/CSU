#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include<stdbool.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"
#include"boundary.h"
#include"sa.h"

void saInitU(SOLVER* solver, CONDITION* inside)
{
    for(int ii=0; ii<solver->Ne; ii++)
    {
        solver->eL[ii]->U[4] = inside->Uin[4];
    }
}

void saGrad(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        double Pmin, Pmax;
        ELEMENT1* E = solver->eL[ii];

        newSolverCalcGrad2(E, 0, 0, &Pmin, &Pmax);

        newSolverCalcGrad2(E, 1, 1, &Pmin, &Pmax);
        
        newSolverCalcGrad2(E, 2, 2, &Pmin, &Pmax);

        newSolverCalcGrad2(E, 4, 3, &Pmin, &Pmax);
        
        newSolverCalcGrad2(E, 5, 4, &Pmin, &Pmax);        
    }
}

void saInterFace(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->Nf; ii++)
    {
        double dSx = solver->fL[ii]->dS[0];
        double dSy = solver->fL[ii]->dS[1];
    
        ELEMENT1* E0 = solver->fL[ii]->eL[0];
        ELEMENT1* E1 = solver->fL[ii]->eL[1];

        double duxm = (E0->grad[1][0] + E1->grad[1][0])*0.5;
        double dvxm = (E0->grad[2][0] + E1->grad[2][0])*0.5;
        double dTxm = (E0->grad[3][0] + E1->grad[3][0])*0.5;
        double dnxm = (E0->grad[4][0] + E1->grad[4][0])*0.5;
        
        double duym = (E0->grad[1][1] + E1->grad[1][1])*0.5;
        double dvym = (E0->grad[2][1] + E1->grad[2][1])*0.5;
        double dTym = (E0->grad[3][1] + E1->grad[3][1])*0.5;
        double dnym = (E0->grad[4][1] + E1->grad[4][1])*0.5;

        double dx = E1->c[0] - E0->c[0];
        double dy = E1->c[1] - E0->c[1];
        double L = sqrt(dx*dx + dy*dy);
        
        double dul = (E1->P[1] - E0->P[1])/L;
        double dvl = (E1->P[2] - E0->P[2])/L;
        double dTl = (E1->P[4] - E0->P[4])/L;
        double dnl = (E1->P[5] - E0->P[5])/L;

        double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

        double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

        double dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        double dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

        double dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
        double dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

        // Flow variables in the face
        double r = (E1->P[0] + E0->P[0])*0.5;
        double u = (E1->P[1] + E0->P[1])*0.5;
        double v = (E1->P[2] + E0->P[2])*0.5;
        double T = (E1->P[4] + E0->P[4])*0.5;
        double n = (E1->P[5] + E0->P[5])*0.5;

        double mi_L = sutherland(T);
        double n_L = mi_L/r;

        double fv1;
        double tx;
        double ty;

        saCalcFace(n, n_L, r, dnx, dny, &fv1, &tx, &ty);

        double mi_t = fv1*r*n;
        double mi = mi_L + mi_t;
        double k = solver->Cp*(mi_L/solver->Pr + mi_t/solver->Pr_t);

	    double txx = 2*mi*(dux - (dux + dvy)/3);
	    double tyy = 2*mi*(dvy - (dux + dvy)/3);
	    double txy = mi*(duy + dvx);

	    solver->fL[ii]->flux[1] += -(txx*dSx + txy*dSy);
	    solver->fL[ii]->flux[2] += -(txy*dSx + tyy*dSy);
	    solver->fL[ii]->flux[3] += -((txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy));
        solver->fL[ii]->flux[4] += -(tx*dSx + ty*dSy);
    }

    
    for(int ii=0; ii<solver->Nbc; ii++)
    {
        saBoundaryFace(solver, solver->bcL[ii]);
    }
    
    
    //newSolverDistributeFlux(solver);
}

void saInterSource(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        ELEMENT1* E0 = solver->eL[ii];

        double drx = E0->grad[0][0];
        double dvx = E0->grad[2][0];
        double dnx = E0->grad[4][0];

        double dry = E0->grad[0][1];
        double duy = E0->grad[1][1];
        double dny = E0->grad[4][1];

        // Flow variables in the face
        double rho = E0->P[0];
        double T = E0->P[4];
        double n = E0->P[5];

        double mi_L = sutherland(T);
        double n_L = mi_L/rho;
        double S = fabs(duy - dvx);
        double Qt;

        saCalcSource(n, n_L, S, solver->dist[ii], rho, drx, dry, dnx, dny, &Qt);

        E0->R[4] -= Qt*E0->omega;
    }
}

void saInter(SOLVER* solver)
{
    saGrad(solver);
    saInterFace(solver);
    saInterSource(solver);
}


void saCalcFace(double ni, double ni_L, double r, double dnix, double dniy, double* fv1, double* tx, double* ty)
{

    double aux;

    double Cv1 = 7.1;
    double X = ni/ni_L;
    aux = X*X*X;
    *fv1 = aux/(aux + Cv1*Cv1*Cv1);

    double sig = 2./3.;
    *tx = r*(ni_L + ni)*dnix/sig;
    *ty = r*(ni_L + ni)*dniy/sig;

    // Transitional terms not included


}

void saCalcSource(double ni, double ni_L, double S, double d, double rho, double drx, double dry, double dnix, double dniy, double* Qt)
{

/*
    Based on: Allmaras S. R. et all, Modifications and Clarifications for the Implementation
of the Spalart-Allmaras Turbulence Model, 2012
*/

    double aux;

    double Cv1 = 7.1;
    double X = ni/ni_L;
    aux = X*X*X;
    double fv1 = aux/(aux + Cv1*Cv1*Cv1);

    double k = 0.41;
    double fv2 = 1. - X/(1. + X*fv1);
    double Sbar = ni*fv2/(k*k*d*d);
    double Stil;

    double cv2 = 0.7;
    double cv3 = 0.9;
    if(Sbar >= - cv2*S)
    {
        Stil = S + Sbar;
    }
    else
    {
        Stil = S + S*(cv2*cv2*S + cv3*Sbar)/((cv3 - 2*cv2)*S - Sbar);
    }

    if(Stil < 0)
    {
        printf("Stil: %e", Stil);
    }

    double sig = 2./3.;
    double Cb1 = 0.1355;
    double Cb2 = 0.622;
    double Cw1 = Cb1/(k*k) + (1. + Cb2)/sig;

    double Cw2 = 0.3;
    double r = ni/(Stil*k*k*d*d);
    r = fmin(r, 10);
    aux = r*r;
    aux = aux*aux*aux;
    double g = r + Cw2*(aux - r);

    double Cw3 = 2.;
    aux = Cw3*Cw3;
    aux = aux*aux*aux;
    double fw = g*pow((1. + aux)/(g*g*g*g*g*g + aux), 1./6.);

    // Transitional terms not included

    *Qt = rho*(Cb1*Stil*ni + Cb2*(dnix*dnix + dniy*dniy)/sig - Cw1*fw*(ni*ni/(d*d))) - (ni_L + ni)*(dnix*drx + dniy*dry)/sig;

}

void saBoundaryFace(SOLVER* solver, SOLVERBC* bc)
{
    # pragma omp parallel for
    for(int ii=0; ii<bc->Nf; ii++)
    {
        FACE1* f = bc->fL[ii];
        ELEMENT1* E0 = f->eL[0];
        ELEMENT1* E1 = f->eL[1];
        double dux, duy, dvx, dvy, dTx, dTy, dnx, dny;
        double F = 1;

        if(bc->flagBC == 0)
        {
            //symmetry
           	double nx, ny, dS;
           	dS = f->dS[2];
            nx = f->dS[0]/dS;
            ny = f->dS[1]/dS;            
            
            double duxm = E1->grad[1][0];
            double dvxm = E1->grad[2][0];
            double dTxm = E1->grad[3][0];
            
            double duym = E1->grad[1][1];
            double dvym = E1->grad[2][1];
            double dTym = E1->grad[3][1];

            dux = duxm - (duxm*nx + duym*ny)*nx;
            duy = duym - (duxm*nx + duym*ny)*ny;        

            dvx = dvxm - (dvxm*nx + dvym*ny)*nx;
            dvy = dvym - (dvxm*nx + dvym*ny)*ny;        
            
            dTx = dTxm - (dTxm*nx + dTym*ny)*nx;
            dTy = dTym - (dTxm*nx + dTym*ny)*ny;        
            
            // Flow variables in the face
            double rho = E1->P[0];
            double u = E1->P[1];
            double v = E1->P[2];
            double T = E1->P[4];
            double n = E1->P[5];

            double mi_L = sutherland(T);
            double n_L = mi_L/rho;

            double fv1;
            double tx;
            double ty;

            saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

            double mi_t = fv1*rho*n;
            double mi = mi_L + mi_t;
            double k = solver->Cp*(mi_L/solver->Pr + mi_t/solver->Pr_t);          
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        f->flux[1] += F*(txx*nx + txy*ny)*dS;
	        f->flux[2] += F*(txy*nx + tyy*ny)*dS;
	        f->flux[3] += F*(u*(txx*nx + txy*ny) + v*(txy*nx + tyy*ny) + k*(dTx*nx + dTy*ny))*dS;
            
        }
        else if(bc->flagBC == 1)
        {
            //inlet
            double dSx, dSy;
            dSx = f->dS[0];
            dSy = f->dS[1];            

            double dx = E1->c[0] - E0->c[0];
            double dy = E1->c[1] - E0->c[1];
            double L = sqrt(dx*dx + dy*dy);

            double dul = (E1->P[1] - solver->inlet->Pin[1])/L;
            double dvl = (E1->P[2] - solver->inlet->Pin[2])/L;
            double dTl = (E1->P[4] - solver->inlet->Pin[4])/L;
            double dnl = (E1->P[5] - solver->inlet->Pin[5])/L;

            double duxm = E1->grad[1][0];
            double dvxm = E1->grad[2][0];
            double dTxm = E1->grad[3][0];
            double dnxm = E1->grad[4][0];
            
            double duym = E1->grad[1][1];
            double dvym = E1->grad[2][1];
            double dTym = E1->grad[3][1];
            double dnym = E1->grad[4][1];

            dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

            dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

            dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
            dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;

            dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
            dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

            // Flow variables in the face
            double rho = solver->inlet->Pin[0];
            double u = solver->inlet->Pin[1];
            double v = solver->inlet->Pin[2];
            double T = solver->inlet->Pin[4];
            double n = solver->inlet->Pin[5];

            double mi_L = sutherland(T);
            double n_L = mi_L/rho;

            double fv1;
            double tx;
            double ty;

            saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

            double mi_t = fv1*rho*n;
            double mi = mi_L + mi_t;
            double k = solver->Cp*(mi_L/solver->Pr + mi_t/solver->Pr_t);

	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);
	        double txy = mi*(duy + dvx);

	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);
	        f->flux[3] += F*(u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy));
            f->flux[4] += F*(tx*dSx + ty*dSy);
            
        }
        else if(bc->flagBC == 3)
        {

            //wall
            double dSx, dSy;
            dSx = f->dS[0];
            dSy = f->dS[1];            

            double dx = E1->c[0] - E0->c[0];
            double dy = E1->c[1] - E0->c[1];
            double L = sqrt(dx*dx + dy*dy);

            double dul = (E1->P[1] - 0)/L;
            double dvl = (E1->P[2] - 0)/L;
            double dnl = (E1->P[5] - 0)/L;

            double duxm = E1->grad[1][0];
            double dvxm = E1->grad[2][0];
            double dnxm = E1->grad[4][0];
            
            double duym = E1->grad[1][1];
            double dvym = E1->grad[2][1];
            double dnym = E1->grad[4][1];

            dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

            dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

            dnx = dnxm + (dnl - (dnxm*dx + dnym*dy)/L)*dx/L;
            dny = dnym + (dnl - (dnxm*dx + dnym*dy)/L)*dy/L;

            // Flow variables in the face
            double rho = E1->P[0];
            double T = E1->P[4];
            double n = 0.0;

            double mi_L = sutherland(T);
            double n_L = mi_L/rho;

            double fv1;
            double tx;
            double ty;

            saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

            double mi = mi_L;

	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);
	        double txy = mi*(duy + dvx);

	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);
            f->flux[4] += F*(tx*dSx + ty*dSy);
            
        }
        else
        {

            //outlet
            double dSx, dSy;
            dSx = f->dS[0];
            dSy = f->dS[1];

            double dux = E1->grad[1][0];
            double dvx = E1->grad[2][0];
            double dTx = E1->grad[3][0];
            double dnx = E1->grad[4][0];
            
            double duy = E1->grad[1][1];
            double dvy = E1->grad[2][1];
            double dTy = E1->grad[3][1];
            double dny = E1->grad[4][1];
            
            // Flow variables in the face
            double rho = E1->P[0];
            double u = E1->P[1];
            double v = E1->P[2];
            double T = E1->P[4];
            double n = E1->P[5];

            double mi_L = sutherland(T);
            double n_L = mi_L/rho;

            double fv1;
            double tx;
            double ty;

            saCalcFace(n, n_L, rho, dnx, dny, &fv1, &tx, &ty);

            double mi_t = fv1*rho*n;
            double mi = mi_L + mi_t;
            double k = solver->Cp*(mi_L/solver->Pr + mi_t/solver->Pr_t);

	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);
	        double txy = mi*(duy + dvx);

	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);
	        f->flux[3] += F*(u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy));
            f->flux[4] += F*(tx*dSx + ty*dSy);
        }
	}
}

void saCalcD(MESH* mesh)
{

    double xVol, yVol;

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        ELEMENT* Evol = mesh->elemL[ii];
        Evol->d = 1e6;
    }

    for(int jj=0; jj<mesh->Nmark; jj++)
    {

        //printf("%i\n", mesh->bc[jj]->flagBC);
        if(mesh->bc[jj]->flagBC==3)
        {
            double tf = 0;

            for(int ii=0; ii<mesh->Nelem; ii++)
            {
                ELEMENT* Evol = mesh->elemL[ii];
                elementCenter(Evol, mesh, &xVol, &yVol);

                for(int kk=0; kk<mesh->bc[jj]->Nelem; kk++)
                {
                    ELEMENT* Esurf = mesh->bc[jj]->elemL[kk];
                    double p0x = mesh->p[Esurf->p[0]][0];
                    double p0y = mesh->p[Esurf->p[0]][1];

                    double p1x = mesh->p[Esurf->p[1]][0];
                    double p1y = mesh->p[Esurf->p[1]][1];

                    double num = (p1x - p0x)*(xVol - p0x) + (p1y - p0y)*(yVol - p0y);
                    double den = (p1x - p0x)*(p1x - p0x) + (p1y - p0y)*(p1y - p0y);

                    double t = num/den;

                    if(t>1.)
                    {
                        t = 1.;
                    }
                    else if(t<0.)
                    {
                        t = 0.;
                    }

                    double dx = (p1x - p0x)*t + p0x - xVol;
                    double dy = (p1y - p0y)*t + p0y - yVol;

                    double d = sqrt(dx*dx + dy*dy);

                    if(Evol->d > d)
                    {
                        Evol->d = d;
                        tf = t;
                    }

                }

                //printf("\ntf: %f\n", tf);

            }
        }
    }

    /*
    for(int jj=0; jj<mesh->Nmark; jj++)
    {

        //printf("%i\n", mesh->bc[jj]->flagBC);
        if(mesh->bc[jj]->flagBC==3)
        {

            for(int ii=0; ii<mesh->Nelem; ii++)
            {

                for(int kk=0; kk<mesh->bc[jj]->Nelem; kk++)
                {
                    ELEMENT* Esurf = mesh->bc[jj]->elemL[kk];
                    ELEMENT* Evol = Esurf->neiL[0];
                    printf("%e\n", Evol->d);

                }

            }

        }

    }
    */

}

void saCalcTensorWall(SOLVER* solver, ELEMENT* E, double* Txx, double* Txy, double* Tyy, double* x, double* yp)
{

    double x0, y0, x1, y1, dSx, dSy;

    int e0 = E->neiL[0]->ii;
    int p0 = E->p[0];
    int p1 = E->p[1];

    ELEMENT* E0 = E->neiL[0];

    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    double dS = sqrt(dSx*dSx + dSy*dSy);


    elementCenter(E0, solver->mesh, &x0, &y0);
    elementCenter(E, solver->mesh, &x1, &y1);

    double dx = x1 - x0;
    double dy = y1 - y0;
    double L = sqrt(dx*dx + dy*dy);

    double dul = (0 - E0->P[1])/L;
    double dvl = (0 - E0->P[2])/L;

    double duxm = solver->dPx[1][e0];
    double dvxm = solver->dPx[2][e0];

    double duym = solver->dPy[1][e0];
    double dvym = solver->dPy[2][e0];

    double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
    double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;

    double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
    double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;

    double T = E0->P[4];
    double mi = sutherland(T);

    //printf("%f\n", E0->P[1]);

    double txx = 2*mi*(dux - (dux + dvy)/3);
    double tyy = 2*mi*(dvy - (dux + dvy)/3);
    double txy = mi*(duy + dvx);

    double fx = (txx*dSx + txy*dSy)/dS;
    double fy = (txy*dSx + tyy*dSy)/dS;
    double ft = (fx*dSy - fy*dSx)/dS;

    double ut = sqrt(fabs(ft)/E0->P[0]);

    *Txx = 0.0;
	*Txy = ft;
	*Tyy = 0.0;
	*x = x1;
	*yp = E0->P[0]*E0->d*ut/mi;
}

void saBoundary(SOLVER* solver)
{

}

void saBoundaryCalc(SOLVER* solver, SOLVERBC* bc)
{  
    # pragma omp parallel for
    for(int ii=0; ii<bc->Nf; ii++)
    { 
        int kk;
        double dSx, dSy, dS;
        double PL[5];
	    double Pb[4];  
        FACE1* f;
        double F = 1;

        f = bc->fL[ii];
 
        dSx = f->dS[0];
        dSy = f->dS[1];
        dS = f->dS[2];
       
        if(f->dS[2] > 0)
        {
        
            for(kk=0; kk<4; kk++)
		    {
			    PL[kk] = f->eL[1]->P[kk];
		    }      		
			PL[4] = f->eL[1]->P[5];

            if(bc->flagBC == 0)
            {
                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
            
                // Reflexive
                fluxAUSMDV_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4], PL[0], -PL[1], PL[2], PL[3], PL[4], f->flux);
            }
            else if(bc->flagBC == 1)
            {
                boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
	            rotation(Pb, f->dS[0], f->dS[1], f->dS[2]);

                fluxAUSMDV_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4],  Pb[0], Pb[1], Pb[2], Pb[3], solver->inlet->Pin[5], f->flux);
            }
            else if(bc->flagBC == 2)
            {           
                boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
	            rotation(Pb, f->dS[0], f->dS[1], f->dS[2]);

                fluxAUSMDV_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4], Pb[0], Pb[1], Pb[2], Pb[3], PL[4], f->flux);
            }
            else if(bc->flagBC == 3)
            {
                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
                
                fluxAUSMDV_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4], PL[0], -PL[1], -PL[2], PL[3], PL[4], f->flux);
            }       

            // Rotation of the flux
		    rotation(f->flux, f->dS[0], -f->dS[1], f->dS[2]);
		                        
            for(kk=0; kk<solver->Nvar; kk++)
            {
                f->flux[kk] *= - F*f->dS[2];
            }
        }
        else
        {
            for(kk=0; kk<solver->Nvar; kk++)
            {
                f->flux[kk] = 0.0;
            }
        }
    }
}

