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
#include"boundary.h"

void boundaryInlet(SOLVER* solver, double* Pa, double* Pd, double* Pb, double nx, double ny)
{
    double rd = Pd[0];
    double ud = Pd[1];
    double vd = Pd[2];
    double pd = Pd[3];
    double c0 = sqrt(solver->gamma*pd/rd);    
    double m = sqrt(ud*ud + vd*vd)/c0;
    
    if(m < 1.)
    {
        double ra = Pa[0];
        double ua = Pa[1];
        double va = Pa[2];
        double pa = Pa[3];

        double pb = 0.5*(pa + pd - rd*c0*(nx*(ua-ud) + ny*(va-vd)));
        double rb = ra + (pb - pa)/(c0*c0);
        double ub = ua - nx*(pa - pb)/(rd*c0);
        double vb = va - ny*(pa - pb)/(rd*c0);

        Pb[0] = rb;
        Pb[1] = ub;
        Pb[2] = vb;
        Pb[3] = pb;
    }
    else
    {
        for(int ii=0; ii<4; ii++)
        {
            Pb[ii] = Pa[ii];
        }
    }    
}

void boundaryOutlet(SOLVER* solver, double* Pd, double* Pb, double nx, double ny)
{

    double rd = Pd[0];
    double ud = Pd[1];
    double vd = Pd[2];
    double pd = Pd[3];
    double c0 = sqrt(solver->gamma*pd/rd);    
    double m = sqrt(ud*ud + vd*vd)/c0;
    
    if(m < 1.)
    {
        double pb = solver->pout;
        double rb = rd + (pb - pd)/(c0*c0);
        double ub = ud + nx*(pd - pb)/(rd*c0);
        double vb = vd + ny*(pd - pb)/(rd*c0);

        Pb[0] = rb;
        Pb[1] = ub;
        Pb[2] = vb;
        Pb[3] = pb;
    }
    else
    {
        for(int ii=0; ii<4; ii++)
        {
            Pb[ii] = Pd[ii];
        }
    }    
}

void boundaryWall(SOLVER* solver, double* Pd, double* Pb, double nx, double ny)
{

    double rd = Pd[0];
    double ud = Pd[1];
    double vd = Pd[2];
    double pd = Pd[3];
    double c0 = sqrt(solver->gamma*pd/rd);
    
    double pb = pd + rd*c0*(nx*ud + ny*vd);
    double rb = rd + (pb - pd)/(c0*c0);
    double ub = ud - nx*(nx*ud + ny*vd);
    double vb = vd - ny*(nx*ud + ny*vd);

    Pb[0] = rb;
    Pb[1] = ub;
    Pb[2] = vb;
    Pb[3] = pb;
    
}

void boundaryCalc(SOLVER* solver, MESHBC* bc)
{
    
	int kk;
    double dSx, dSy, dS;
    double aux;
    double PL[4];
	double Pb[4];
    double f[4];
    int e0, p0, p1;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
 
        e0 = bc->elemL[ii]->neiL[0]->ii;
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = solver->mesh->elemL[e0]->P[kk];
		}      		
        
        if(bc->flagBC == 0)
        {

            // Rotation of the velocity vectors
            rotation(PL, dSx, dSy, dS);
        
            // Reflexive
            flux(solver, PL[0], PL[1], PL[2], PL[3], PL[0], -PL[1], PL[2], PL[3], f);

        }
        else if(bc->flagBC == 1)
        {

            boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

            // Rotation of the velocity vectors
            rotation(PL, dSx, dSy, dS);
	        rotation(Pb, dSx, dSy, dS);

            flux(solver, PL[0], PL[1], PL[2], PL[3], Pb[0], Pb[1], Pb[2], Pb[3], f);
        
        }
        else if(bc->flagBC == 2)
        {           

            boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

            // Rotation of the velocity vectors
            rotation(PL, dSx, dSy, dS);
	        rotation(Pb, dSx, dSy, dS);
        
            flux(solver, PL[0], PL[1], PL[2], PL[3], Pb[0], Pb[1], Pb[2], Pb[3], f);
		
        }
        else if(bc->flagBC == 3)
        {
            
            if(solver->laminar==1)
            {
                // Rotation of the velocity vectors
                rotation(PL, dSx, dSy, dS);
        
                // Reflexive
                flux(solver, PL[0], PL[1], PL[2], PL[3], PL[0], -PL[1], -PL[2], PL[3], f);
            }
            else
            {            

                
                boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);

                // Rotation of the velocity vectors
                rotation(PL, dSx, dSy, dS);
	            rotation(Pb, dSx, dSy, dS);
            
                //flux(solver, PL[0], PL[1], PL[2], PL[3], Pb[0], Pb[1], Pb[2], Pb[3], f);
                fluxFree(solver, Pb[0], Pb[1], Pb[2], Pb[3], f);
                

                /*
		        double p2 = PL[3];
	    
                // Outlet
                f[0] = .0;
                f[2] = .0;
                f[3] = .0;

                f[1] =  p2;
                */
            }
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

void boundaryCalcVisc(SOLVER* solver, MESHBC* bc)
{
    double dSx, dSy;
    int e0, p0, p1;
    double x0, x1, y0, y1;
    double dux, duy, dvx, dvy, dTx, dTy;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        e0 = bc->elemL[ii]->neiL[0]->ii;
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];

        ELEMENT* E0 = bc->elemL[ii]->neiL[0];

        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

        if(bc->flagBC == 1)
        {
            elementCenter(E0, solver->mesh, &x0, &y0);

           	x1 = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
            y1 = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

            double dx = x1 - x0;		    
            double dy = y1 - y0;
            double L = sqrt(dx*dx + dy*dy);
            
            double dul = (solver->inlet->Pin[1] - E0->P[1])/L;
            double dvl = (solver->inlet->Pin[2] - E0->P[2])/L;            
            double dTl = (solver->inlet->Pin[4] - E0->P[4])/L; 

            double duxm = solver->dPx[1][e0];
            double dvxm = solver->dPx[2][e0];
            double dTxm = solver->dPx[3][e0];
            
            double duym = solver->dPy[1][e0];
            double dvym = solver->dPy[2][e0];
            double dTym = solver->dPy[3][e0];

            dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

            dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;  
            
            dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
            dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;
            
            double T = E0->P[4];
            double mi = sutherland(T);
            double k = solver->Cp*mi/solver->Pr;            
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        double u = E0->P[1];
	        double v = E0->P[2];		        
	        
	        solver->R[1][e0] -= txx*dSx + txy*dSy;
	        solver->R[2][e0] -= txy*dSx + tyy*dSy;
	        solver->R[3][e0] -= u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy);
            
        }
        else if(bc->flagBC == 3)
        {
            elementCenter(E0, solver->mesh, &x0, &y0);

           	x1 = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
            y1 = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;

            double dx = x1 - x0;		    
            double dy = y1 - y0;
            double L = sqrt(dx*dx + dy*dy);
            
            double dul = (0 - E0->P[1])/L;
            double dvl = (0 - E0->P[2])/L;            

            double duxm = solver->dPx[1][e0];
            double dvxm = solver->dPx[2][e0];
            
            double duym = solver->dPy[1][e0];
            double dvym = solver->dPy[2][e0];

            dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

            dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;
                            
            double T = E0->P[4];
            double mi = sutherland(T);
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        solver->R[1][e0] -= txx*dSx + txy*dSy;
	        solver->R[2][e0] -= txy*dSx + tyy*dSy;
            
        }
        else
        {
            dux = solver->dPx[1][e0];
            dvx = solver->dPx[2][e0];
            dTx = solver->dPx[3][e0];
            
            duy = solver->dPy[1][e0];
            dvy = solver->dPy[2][e0];
            dTy = solver->dPy[3][e0];                

            double T = E0->P[4];
            double mi = sutherland(T);
            double k = solver->Cp*mi/solver->Pr;            
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        double u = E0->P[1];
	        double v = E0->P[2];		        
	        
	        solver->R[1][e0] -= txx*dSx + txy*dSy;
	        solver->R[2][e0] -= txy*dSx + tyy*dSy;
	        solver->R[3][e0] -= u*(txx*dSx + txy*dSy) + v*(txy*dSx + tyy*dSy) + k*(dTx*dSx + dTy*dSy);            
        }
	}        
}



void boundary(SOLVER* solver)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        boundaryCalc(solver, solver->mesh->bc[ii]);
    }
}

void boundaryVisc(SOLVER* solver)
{
    for(int ii=0; ii<solver->mesh->Nmark; ii++)
    {
        boundaryCalcVisc(solver, solver->mesh->bc[ii]);
    }
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


void boundaryCalcPrimitive(SOLVER* solver, MESHBC* bc)
{
    
	int kk;
    double dSx, dSy, dS;
    double PL[4];
	double Pb[4];
    int p0, p1;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
 
        ELEMENT* E0 = bc->elemL[ii]->neiL[0];
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E0->P[kk];
		}      		
        
        if(bc->flagBC == 0)
        {
            if(dS > 0)
            {
                rotation(PL, dSx, dSy, dS);
                PL[1] = 0.0;
                rotation(PL, dSx, -dSy, dS);
            }
            else
            {
                PL[2] = 0.0;                
            }
            
            for(kk=0; kk<4; kk++)
            {
                bc->elemL[ii]->P[kk] = PL[kk];
            }
            
            bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->Rgas);
        
        }
        else if(bc->flagBC == 1)
        {

            boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

            for(kk=0; kk<4; kk++)
            {
                bc->elemL[ii]->P[kk] = Pb[kk];
            }
            
            bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->Rgas);            
        
        }
        else if(bc->flagBC == 2)
        {           

            boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

            for(kk=0; kk<4; kk++)
            {
                bc->elemL[ii]->P[kk] = Pb[kk];
            }
            
            bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->Rgas);            
		
        }
        else if(bc->flagBC == 3)
        {
            
            if(solver->laminar==1)
            {
                bc->elemL[ii]->P[0] = PL[0];
                bc->elemL[ii]->P[1] = 0.0;
                bc->elemL[ii]->P[2] = 0.0;
                bc->elemL[ii]->P[3] = PL[3];

                bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->Rgas);
            }
            else
            {                        
                if(dS > 0)
                {
                    //rotation(PL, dSx, dSy, dS);
                    //PL[1] = 0.0;
                    //rotation(PL, dSx, -dSy, dS);
                    
                    boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
                    for(kk=0; kk<4; kk++)
                    {
                        bc->elemL[ii]->P[kk] = Pb[kk];
                    }                    
                }
                else
                {
                    PL[2] = 0.0;                
                    for(kk=0; kk<4; kk++)
                    {
                        bc->elemL[ii]->P[kk] = PL[kk];
                    }                    
                }
                
                
                bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->Rgas);
            }
        }       
    }
}

void boundaryCalcFrictionWall(SOLVER* solver, ELEMENT* E, double* fx, double* fy)
{

        double x0, y0, x1, y1, dSx, dSy;

        int e0 = E->neiL[0]->ii;
        int p0 = E->p[0];
        int p1 = E->p[1];

        ELEMENT* E0 = E->neiL[0];
         
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);

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
        
        *fx = txx*dSx + txy*dSy;
		*fy = txy*dSx + tyy*dSy;
 
}
