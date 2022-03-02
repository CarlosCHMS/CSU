#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny)
{
  
    CONDITION* cond = malloc(sizeof(CONDITION));

    cond->p = p;
    cond->T = T;
    cond->mach = mach;
    cond->nx = nx;
    cond->ny = ny;
    
    return cond;

}

void conditionState(CONDITION* cond, SOLVER* solver)
{

    double r, u, v, E, c;
    
    r = cond->p/(solver->Rgas*cond->T);
    c = sqrt(solver->gamma*solver->Rgas*cond->T);
    u = cond->nx*cond->mach*c;
    v = cond->ny*cond->mach*c;
    E = (solver->Rgas*cond->T)/(solver->gamma-1) + (u*u + v*v)/2;
    
    cond->Uin[0] = r;
    cond->Uin[1] = r*u;
    cond->Uin[2] = r*v;
    cond->Uin[3] = r*E;   

}

double conditionVref(CONDITION* cond, SOLVER* solver)
{
      
    double c = sqrt(solver->gamma*solver->Rgas*cond->T);
    double Vref = (cond->mach+1)*c;
        
    return Vref;

}

void solverFree(SOLVER* solver)
{

    tableFreeDouble(solver->U, 4);
    tableFreeDouble(solver->Uaux, 4);
    tableFreeDouble(solver->R, 4);        
    tableFreeDouble(solver->faceFlux, 4);
    tableFreeDouble(solver->dUx, 4);
    tableFreeDouble(solver->dUy, 4);
    meshFree(solver->mesh);
    
}

void solverWrite(SOLVER* solver, char* fileName)
{

    /* 
        Based on: Adek Tasri, Accuracy of Cell Centres to Vertices 
        Interpolation for Unstructured Mesh Finite Volume Solver, 2021
    */

    FILE* ff = fopen(fileName, "w");
    double** Up = tableMallocDouble(4, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    int ii, jj, kk, p;
    double xc, yc, xp, yp, L;

    for(ii=0; ii<solver->mesh->Np; ii++)
    {       
        for(kk=0; kk<4; kk++)
        {
            Up[kk][ii] = 0.;
        }   
        den[ii] = 0.;
    }

    for(ii=0; ii<solver->mesh->Nelem; ii++)
    {
        meshElemCenter(solver->mesh, ii, &xc, &yc);
        
        for(jj=0; jj<3; jj++)
        {
            p = solver->mesh->elem[ii][jj];
            xp = solver->mesh->p[p][0];
            yp = solver->mesh->p[p][1];
            L = sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
           
            for(kk=0; kk<4; kk++)
            {        
                Up[kk][p] += solver->U[kk][ii]/L;
            }
            
            den[p] += 1/L;
        }        
        
    }


    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        if(den[ii] != 0)
        {
            for(kk=0; kk<4; kk++)
            {
                Up[kk][ii] /= den[ii];    
            }
        }
    }

    fprintf(ff, "0, %i, 4,\n", solver->mesh->Np);

    for(ii=0; ii<solver->mesh->Np; ii++)
    {        
        for(jj=0; jj<4; jj++)
        {
            fprintf(ff, "%.10e, ", Up[jj][ii]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);

    tableFreeDouble(Up, 4);
    free(den);

}

void solverInitU(SOLVER* solver, CONDITION* inside)
{

    conditionState(inside, solver);

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            solver->U[kk][ii] = inside->Uin[kk];
        }
    }
}

void solverResetR(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            solver->R[kk][ii] = 0.0; 
        }
    }
}

double solverCalcP(SOLVER* solver, double** U, int ii)
{

	double u = U[1][ii]/U[0][ii];
	double v = U[2][ii]/U[0][ii];
    
	return (solver->gamma - 1)*(U[3][ii] - 0.5*(u*u + v*v)*U[0][ii]);
    
}

void solverCalcVel(SOLVER* solver, double** U, int ii, double* u, double* v, double* c)
{
    
    double E = U[3][ii]/U[0][ii];
    double aux;

    *u = U[1][ii]/U[0][ii];
    *v = U[2][ii]/U[0][ii];
    
    aux = (E - ((*u)*(*u) + (*v)*(*v))/2);
    aux *= solver->gamma - 1;
    *c = sqrt(aux*solver->gamma);
    
}

void rotation(double* U, double dSx, double dSy, double dS)
{

    double aux = (U[1]*dSx + U[2]*dSy)/dS;
    U[2] = (-U[1]*dSy + U[2]*dSx)/dS;
    U[1] = aux;
	
}

void inter(SOLVER* solver, double **U)
{

    if(solver->order == 2)
	{
        # pragma omp parallel for
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
	        int jj, kk, p0, p1;
            double dUx, dUy;
            double x0, y0, d2, phi0, phi, xm, ym;
            double Umin, Umax; 
            double blend = 1.0;        

	        if(solver->mesh->neiN[ii] > 1)
	        {  
          	    meshElemCenter(solver->mesh, ii, &x0, &y0);
          	    
                for(kk=0; kk<4; kk++)
                {
                    solverCalcGrad2(solver, U[kk], ii, &dUx, &dUy, &Umin, &Umax);
                    
                    for(jj=0; jj<3; jj++)
                    {
                        p0 = solver->mesh->elem[ii][jj]; 
                        p1 = solver->mesh->elem[ii][(jj+1)%3];
                        xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
                        ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;        
                        d2 = (dUx*(xm - x0) + dUy*(ym - y0));
                        //phi0 = limiterBJ(U[kk][ii], Umin, Umax, d2);                        
                        phi0 = limiterV(U[kk][ii], Umin, Umax, d2, solver->e);
                        
                        if(jj==0)
                        {
                            phi = phi0;
                        }
                        else
                        {
                            phi = fmin(phi, phi0);
                        }			            
                    }

                    solver->dUx[kk][ii] = phi*dUx*blend;
	                solver->dUy[kk][ii] = phi*dUy*blend;
	            }   
	        }
        }
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Ncon; ii++)
    {
	    int kk;
        double dSx, dSy, dS;
        double x0, y0, x1, y1, xm, ym;
        double UL[4];
	    double UR[4];
        double f[4];
 
        int e0 = solver->mesh->con[ii][0];
        int e1 = solver->mesh->con[ii][1];
        int p0 = solver->mesh->con[ii][2];
        int p1 = solver->mesh->con[ii][3];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        if(solver->order == 2)
        {
           	xm = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
            ym = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

		    meshElemCenter(solver->mesh, e0, &x0, &y0);
		    meshElemCenter(solver->mesh, e1, &x1, &y1);
		    
            for(kk=0; kk<4; kk++)
		    {                
		        if(solver->mesh->neiN[e0] > 1)
		        {		    
		            UL[kk] = U[kk][e0] + solver->dUx[kk][e0]*(xm - x0) + solver->dUy[kk][e0]*(ym - y0);
		        }
		        else
		        {
		            UL[kk] = U[kk][e0];
		        }

		        if(solver->mesh->neiN[e1] > 1)
		        {      	        
		            UR[kk] = U[kk][e1] + solver->dUx[kk][e1]*(xm - x1) + solver->dUy[kk][e1]*(ym - y1);
		        }
		        else
		        {
		            UR[kk] = U[kk][e1];
		        }
            }   
        }
        else
        {       
            for(kk=0; kk<4; kk++)
		    {
			    UL[kk] = U[kk][e0];
			    UR[kk] = U[kk][e1];			    
		    }
		}
		
        // Rotation of the velocity vectors
		rotation(UL, dSx, dSy, dS);
                    
        // Rotation of the velocity vectors
		rotation(UR, dSx, dSy, dS);
        
        // Flux calculation
        flux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);

        // Rotation of the flux
		rotation(f, dSx, -dSy, dS);
        
        for(kk=0; kk<4; kk++)
        {
            solver->faceFlux[kk][ii] = f[kk]*dS;
        }
    }
    
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<solver->mesh->faceN[ii]; jj++)
        {
            int face = solver->mesh->elemFace[ii][jj];
            if(face > 0)
            {
                for(int kk=0; kk<4; kk++)
                {
                    solver->R[kk][ii] += solver->faceFlux[kk][face-1];
                }
            }
            else
            {
                for(int kk=0; kk<4; kk++)
                {
                    solver->R[kk][ii] -= solver->faceFlux[kk][-face-1];
                }            
            }
        }
    } 
}

void boundaryInlet(SOLVER* solver, double* Ua, double* Ud, double* Ub, double nx, double ny)
{
    double rd = Ud[0];
    double ud = Ud[1]/Ud[0];
    double vd = Ud[2]/Ud[0];
    double pd = (solver->gamma - 1)*(Ud[3] - 0.5*(ud*ud + vd*vd)*rd);
    double c0 = sqrt(solver->gamma*pd/rd);    
    double m = sqrt(ud*ud + vd*vd)/c0;
    
    if(m < 1.)
    {
        double ra = Ua[0];
        double ua = Ua[1]/Ua[0];
        double va = Ua[2]/Ua[0];
        double pa = (solver->gamma - 1)*(Ua[3] - 0.5*(ua*ua + va*va)*ra);

        double pb = 0.5*(pa + pd - rd*c0*(nx*(ua-ud) + ny*(va-vd)));
        double rb = ra + (pb - pa)/(c0*c0);
        double ub = ua - nx*(pa - pb)/(rd*c0);
        double vb = va - ny*(pa - pb)/(rd*c0);

        Ub[0] = rb;
        Ub[1] = rb*ub;
        Ub[2] = rb*vb;
        Ub[3] = pb/(solver->gamma-1) + 0.5*(ub*ub + vb*vb)*rb;
    }
    else
    {
        for(int ii=0; ii<4; ii++)
        {
            Ub[ii] = Ua[ii];
        }
    }    
}

void boundaryOutlet(SOLVER* solver, double* Ud, double* Ub, double nx, double ny)
{

    double rd = Ud[0];
    double ud = Ud[1]/Ud[0];
    double vd = Ud[2]/Ud[0];
    double pd = (solver->gamma - 1)*(Ud[3] - 0.5*(ud*ud + vd*vd)*rd);
    double c0 = sqrt(solver->gamma*pd/rd);    
    double m = sqrt(ud*ud + vd*vd)/c0;
    
    if(m < 1.)
    {
        double pb = solver->pout;
        double rb = rd + (pb - pd)/(c0*c0);
        double ub = ud + nx*(pd - pb)/(rd*c0);
        double vb = vd + ny*(pd - pb)/(rd*c0);

        Ub[0] = rb;
        Ub[1] = rb*ub;
        Ub[2] = rb*vb;
        Ub[3] = pb/(solver->gamma-1) + 0.5*(ub*ub + vb*vb)*rb;
    }
    else
    {
        for(int ii=0; ii<4; ii++)
        {
            Ub[ii] = Ud[ii];
        }
    }    
}

void boundaryWall(SOLVER* solver, double* Ud, double* Ub, double nx, double ny)
{

    double rd = Ud[0];
    double ud = Ud[1]/Ud[0];
    double vd = Ud[2]/Ud[0];
    double pd = (solver->gamma - 1)*(Ud[3] - 0.5*(ud*ud + vd*vd)*rd);
    double c0 = sqrt(solver->gamma*pd/rd);
    
    double pb = pd + rd*c0*(nx*ud + ny*vd);
    double rb = rd + (pb - pd)/(c0*c0);
    double ub = ud - nx*(nx*ud + ny*vd);
    double vb = vd - ny*(nx*ud + ny*vd);

    Ub[0] = rb;
    Ub[1] = rb*ub;
    Ub[2] = rb*vb;
    Ub[3] = pb/(solver->gamma-1) + 0.5*(ub*ub + vb*vb)*rb;
    
}

void boundaryCalc(SOLVER* solver, double **U, MESHBC* bc)
{
    
	int kk;
    double dSx, dSy, dS;
    double aux;
    double UL[4];
	double Ub[4];
    double f[4];
    int e0, p0, p1;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
 
        e0 = bc->domain[ii];
        p0 = bc->elem[ii][0];
        p1 = bc->elem[ii][1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			UL[kk] = U[kk][e0];
		}      		
        
        if(bc->flagBC == 0)
        {

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
        
            // Reflexive
            flux(solver, UL[0], UL[1], UL[2], UL[3], UL[0], -UL[1], UL[2], UL[3], f);

        }
        else if(bc->flagBC == 1)
        {

            boundaryInlet(solver, solver->inlet->Uin, UL, Ub, dSx/dS, dSy/dS);

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
	        rotation(Ub, dSx, dSy, dS);

            //solverFlux(solver, UL[0], UL[1], UL[2], UL[3], Ub[0], Ub[1], Ub[2], Ub[3], f);
		    fluxFree(solver, Ub[0], Ub[1], Ub[2], Ub[3], f);
        
        }
        else if(bc->flagBC == 2)
        {           

            boundaryOutlet(solver, UL, Ub, dSx/dS, dSy/dS);

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
	        rotation(Ub, dSx, dSy, dS);
        
            flux(solver, UL[0], UL[1], UL[2], UL[3], Ub[0], Ub[1], Ub[2], Ub[3], f);
		
        }
        else if(bc->flagBC == 3)
        {

            /*
            boundaryWall(solver, UL, Ub, dSx/dS, dSy/dS);

            // Rotation of the velocity vectors
            rotation(UL, dSx, dSy, dS);
	        rotation(Ub, dSx, dSy, dS);
        
            flux(solver, UL[0], UL[1], UL[2], UL[3], Ub[0], Ub[1], Ub[2], Ub[3], f);
            //fluxFree(solver, Ub[0], Ub[1], Ub[2], Ub[3], f);
            */

		    double p2 = solverCalcP(solver, U, e0);
	
            // Outlet
            f[0] = .0;
            f[2] = .0;
            f[3] = .0;

            f[1] =  p2;
        }       

        // Rotation of the flux
		rotation(f, dSx, -dSy, dS);
        
        if(dS > 0)
        {             
            for(kk=0; kk<4; kk++)
            {
                aux = f[kk]*dS;
                solver->R[kk][e0] += aux;
            }
        } 
    }
}

void boundary(SOLVER* solver, double **U)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        boundaryCalc(solver, U, solver->mesh->bc[ii]);
    }
}

void interAxisPressure(SOLVER* solver, double **U)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double dS;        
        dS = meshCalcDSlateral(solver->mesh, ii);
        solver->R[2][ii] -= solverCalcP(solver, U, ii)*dS;
    }
}

void solverCalcR(SOLVER* solver, double** U)
{

    solverResetR(solver);

    inter(solver, U);
    
    boundary(solver, U); 
    
    if(solver->mesh->axi==1)
    {
        interAxisPressure(solver, U);
    }    

}

void solverRK(SOLVER* solver, double a)
{   
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        double omega = meshCalcOmega(solver->mesh, ii);
        for(int kk=0; kk<4; kk++)
        {
            
            solver->Uaux[kk][ii] = solver->U[kk][ii] - solver->dt*a*solver->R[kk][ii]/omega;
            
        }
    }
}

void solverUpdateU(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            solver->U[kk][ii] = solver->Uaux[kk][ii];
        }
    }
}

void solverStepRK(SOLVER* solver)
{   

    if(solver->stages==3)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.1481);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==4)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.0833);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.2069);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4265);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==5)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.0533);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.1263);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.2375);

        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4414);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    
    solverUpdateU(solver);

}

void solverCalcRes(SOLVER* solver)
{
    
    
    for(int kk=0; kk<4; kk++)
    {
        solver->res[kk] = fabs(solver->R[kk][0]);
    }
    
    for(int ii=1; ii<solver->mesh->Nelem; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            if(solver->res[kk] < fabs(solver->R[kk][ii]))
            {
                solver->res[kk] = fabs(solver->R[kk][ii]);
            }
        }
    }
    
    for(int kk=0; kk<4; kk++)
    {
        printf(" %+.4e,", solver->res[kk]);        
    }
    printf("\n");    
   
}

int boundaryChoice(char* s)
{
    int ans;

    if(strcmp(s, "symmetry") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "inlet") == 0)
    {
        ans = 1;
    }
    else if(strcmp(s, "outlet") == 0)
    {
        ans = 2;
    }
    else if(strcmp(s, "wall") == 0)
    {
        ans = 3;
    }
    else
    {
        printf("Error: incorrent input of bc: %s\n", s);
        exit(0);
    }
 
    return ans;

}

void boundaryGetBC(MESH* mesh, INPUT* input)
{
    char s[50];
    for(int ii=0; ii<mesh->Nmark; ii++)
    {
        s[0] = '\0';
        strcat(s, "BC:");
        strcat(s, mesh->bc[ii]->name);
        mesh->bc[ii]->flagBC = boundaryChoice(inputGetValue(input, s));

    }
}

double solverLocalTimeStep(SOLVER* solver, int ii)
{
    int p0, p1;
    double dSx, dSy, dSxm, dSym, Lx, Ly, u, v, c;
    
    dSxm = 0.;
    dSym = 0.;
    
    p0 = solver->mesh->elem[ii][0];
    p1 = solver->mesh->elem[ii][1];        
    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    dSxm += fabs(dSx);
    dSym += fabs(dSy);
    
    p0 = solver->mesh->elem[ii][1];
    p1 = solver->mesh->elem[ii][2];        
    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    dSxm += fabs(dSx);
    dSym += fabs(dSy);

    p0 = solver->mesh->elem[ii][2];
    p1 = solver->mesh->elem[ii][0];        
    meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
    dSxm += fabs(dSx);
    dSym += fabs(dSy);
    
    dSxm *= 0.5;
    dSym *= 0.5;
    
    solverCalcVel(solver, solver->U, ii, &u, &v, &c);
    
    Lx = (fabs(u) + c)*dSxm;
    Ly = (fabs(v) + c)*dSym;
    
    return meshCalcOmega(solver->mesh, ii)/(Lx + Ly);
    
}

double solverCalcDt(SOLVER* solver)
{

    double dtLocal;
    double dt = solverLocalTimeStep(solver, 0);

    for(int ii=1; ii<solver->mesh->Nelem; ii++)
    {
        dtLocal = solverLocalTimeStep(solver, ii);
        if(dtLocal < dt)
        {
            dt = dtLocal;
        }
    }

    dt *= solver->CFL*0.5*solver->stages;

    return dt;
}

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm)
{

    double x, y;

    conditionState(inside1, solver);
    conditionState(inside2, solver);

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->mesh->Nelem; ii++)
        {
            meshElemCenter(solver->mesh, ii, &x, &y);
            if(x<xm)
            {
                solver->U[kk][ii] = inside1->Uin[kk];
            }
            else
            {
                solver->U[kk][ii] = inside2->Uin[kk];
            }
        }
    }
}

void solverCalcGrad(SOLVER* solver, double* U, int ii, double* dUx, double* dUy)
{

    /*
    Green-Gauss calculation
    */


    double dSx, dSy, aux;
    
    *dUx = 0.0;
    *dUy = 0.0;
    
    for(int jj=0; jj<3; jj++)
    {
        meshCalcDS(solver->mesh, solver->mesh->neip0[ii][jj], solver->mesh->neip1[ii][jj], &dSx, &dSy);
        aux = (U[ii] + U[solver->mesh->nei[ii][jj]])*0.5;
        *dUx += dSx*aux;
        *dUy += dSy*aux;
    }
 
    aux = meshCalcOmega(solver->mesh, ii);
    
    *dUx /= aux;
    *dUy /= aux;

}

void solverCalcGrad2(SOLVER* solver, double* U, int ii, double* dUx, double* dUy, double* Umin, double* Umax)
{

    /*
    Least squares calculation
    */

    int e;
    double xx, yy;
    double x[4];
    double y[4];    
    double u[4];    
    double a, b, c, A, B;
    int nN = solver->mesh->neiN[ii];
        
    meshElemCenter(solver->mesh, ii, &xx, &yy);
    x[0] = xx;
    y[0] = yy;
    u[0] = U[ii];

    for(int jj=0; jj<nN; jj++)
    {
        e = solver->mesh->nei[ii][jj];
        meshElemCenter(solver->mesh, e, &xx, &yy);
        x[jj+1] = xx;
        y[jj+1] = yy;
        u[jj+1] = U[e];
    }
    
    a = 0;    
    b = 0;
    c = 0;
    A = 0;
    B = 0;
    for(int jj=1; jj<nN+1; jj++)
    {
        a += (x[jj] - x[0])*(x[jj] - x[0]);
        b += (x[jj] - x[0])*(y[jj] - y[0]);
        c += (y[jj] - y[0])*(y[jj] - y[0]);
        A += (x[jj] - x[0])*(u[jj] - u[0]);
        B += (y[jj] - y[0])*(u[jj] - u[0]);
    }
        
    *dUx = (c*A - b*B)/(a*c - b*b);
    *dUy = (-b*A + a*B)/(a*c - b*b);
    
    *Umin = u[0];
    *Umax = u[0];
    
    for(int ii=1; ii<nN+1; ii++)
    {
        *Umax = fmax(*Umax, u[ii]);
        *Umin = fmin(*Umin, u[ii]);        
    }
}

void solverCheckGrad(SOLVER* solver)
{

    double *U = malloc(solver->mesh->Nelem*sizeof(double));
    double *xx = malloc(solver->mesh->Nelem*sizeof(double));
    double x, y;
    double dUx, dUy;
    double Umin, Umax;
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        meshElemCenter(solver->mesh, ii, &x, &y);
        xx[ii] = x;
        U[ii] = x*x;
    }
    
    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        if(solver->mesh->neiN[ii] == 2)
        {
            solverCalcGrad2(solver, U, ii, &dUx, &dUy, &Umin, &Umax);
            printf("%+e, %+e, %+e\n", dUx, 2*xx[ii], dUx - 2*xx[ii]);
            //printf("%+e, %+e\n", Umin, Umax);            
        }
    }

    free(U);
    free(xx);

}

double limiterBJ(double Ui, double Umin, double Umax, double d2)
{

    double ans;
    if(d2 == 0)
    {
        ans = 1;
    }
    else if(d2 > 0)
    {
        ans = fmin(1, (Umax - Ui)/d2);
    }
    else
    {
        ans = fmin(1, (Umin - Ui)/d2);
    }
    
    return ans;

}

double limiterV(double Ui, double Umin, double Umax, double d2, double e)
{

    double ans;
    double d1max = Umax - Ui;
    double d1min = Umin - Ui;
    
    if(d2 == 0)
    {
        ans = 1;
    }
    else if(d2 > 0)
    {
        ans = (d1max*d1max + e*e)*d2 + 2*d2*d2*d1max;
        ans /= d1max*d1max + 2*d2*d2 + d1max*d2 + e*e;
        ans /= d2;
    }
    else
    {
        ans = (d1min*d1min + e*e)*d2 + 2*d2*d2*d1min;
        ans /= d1min*d1min + 2*d2*d2 + d1min*d2 + e*e;
        ans /= d2;
    }
    
    return ans;

}
