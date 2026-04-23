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
#include"boundary.h"
#include"gasprop.h"
#include"laminar.h"
#include"sa.h"
#include"sst.h"


BOUNDARY* boundaryInit(INPUT* input, MESHBC* bc, SOLVER* solver)
{

    BOUNDARY* boundary = malloc(sizeof(BOUNDARY));
    boundary->bc = bc;
    
    char s[50];
    strcpy(s, "BC:");
    strcat(s, boundary->bc->name);
    if(inputNameIsInput(input, s))
    {
        strcpy(boundary->type, inputGetValue(input, s));
    }
    else
    {
        printf("Error: boundary type %s not available in input file", s);
        exit(0);
    }
    
    if(strcmp(boundary->type, "symmetry") == 0)
    {
        boundary->primitive = boundaryPrimitiveSymmetry;
        boundary->convective = boundaryConvectiveSymmetry;
        if(solver->laminar)
        {
            boundary->viscousFlux = laminarBoundaryViscousFluxSymmetry;
        }
        else if(solver->sa1->active)
        {
            boundary->viscousFlux = saBoundaryViscousFluxSymmetry;
        }
    }
    else if(strcmp(boundary->type, "inlet") == 0)
    {
        boundary->primitive = boundaryPrimitiveInlet;
        boundary->convective = boundaryConvectiveGeneral;
        if(solver->laminar)
        {
            boundary->viscousFlux = laminarBoundaryViscousFluxGeneral;
        }
        else if(solver->sa1->active)
        {
            boundary->viscousFlux = saBoundaryViscousFluxGeneral;
        }

    }
    else if(strcmp(boundary->type, "outlet") == 0)
    {
        boundary->primitive = boundaryPrimitiveOutlet;
        boundary->convective = boundaryConvectiveGeneral;
        if(solver->laminar)
        {
            boundary->viscousFlux = laminarBoundaryViscousFluxGeneral;
        }
        else if(solver->sa1->active)
        {
            boundary->viscousFlux = saBoundaryViscousFluxGeneral;
        }

    }
    else if(strcmp(boundary->type, "wall") == 0)
    {
        boundary->primitive = boundaryPrimitiveWall;
        boundary->convective = boundaryConvectiveGeneral;
        if(solver->laminar)
        {
            boundary->viscousFlux = laminarBoundaryViscousFluxWall;
        }
        else if(solver->sa1->active)
        {
            boundary->viscousFlux = saBoundaryViscousFluxWall;
        }

    }
    else if(strcmp(boundary->type, "wallT") == 0)
    {
        boundary->primitive = boundaryPrimitiveWallT;
        boundary->convective = boundaryConvectiveGeneral;
        if(solver->laminar)
        {
            boundary->viscousFlux = laminarBoundaryViscousFluxGeneral;
        }
        else if(solver->sa1->active)
        {
            boundary->viscousFlux = saBoundaryViscousFluxGeneral;
        }
    }
    else
    {
        printf("Error: incorrent input of bc: %s\n", s);
        exit(0);
    }
    
    return boundary;
}



void boundaryInlet(SOLVER* solver, double* Pa, double* Pd, double* Pb, double nx, double ny)
{
    double rd = Pd[0];
    double ud = Pd[1];
    double vd = Pd[2];
    double pd = Pd[3];
    double T = pd/rd/solver->gas->R;
    double c0 = gasprop_T2c(solver->gas, T);    
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
    double T = pd/rd/solver->gas->R;
    double c0 = gasprop_T2c(solver->gas, T); 
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
    double T = pd/rd/solver->gas->R;
    double c0 = gasprop_T2c(solver->gas, T); 
    
    double pb = pd + rd*c0*(nx*ud + ny*vd);
    double rb = rd + (pb - pd)/(c0*c0);
    double ub = ud - nx*(nx*ud + ny*vd);
    double vb = vd - ny*(nx*ud + ny*vd);

    Pb[0] = rb;
    Pb[1] = ub;
    Pb[2] = vb;
    Pb[3] = pb;
    
}

void boundaryConvectiveGeneral(BOUNDARY* boundary, SOLVER* solver)
{
    
    MESHBC* bc = boundary->bc;
    
    int kk;
    double dSx, dSy, dS;
    double aux;
    double PL[6];
    double Pb[6];
    double f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int e0, p0, p1;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
 
        e0 = bc->elemL[ii]->neiL[0]->ii;
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<solver->Nvar; kk++)
	    {
	        if(kk > 3)
	        {
		        PL[kk] = solver->mesh->elemL[e0]->P[kk+1];
		        Pb[kk] = bc->elemL[ii]->P[kk+1];
		    }
		    else
		    {
		        PL[kk] = solver->mesh->elemL[e0]->P[kk];
		        Pb[kk] = bc->elemL[ii]->P[kk];
		    }
	    }      		
                
        // Rotation of the velocity vectors
        rotation(PL, dSx, dSy, dS);
        rotation(Pb, dSx, dSy, dS);
        
        // Flux calculation
        solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);

        // Rotation of the flux
	    rotation(f, dSx, -dSy, dS);
        
        if(dS > 0)
        {             
            for(kk=0; kk<solver->Nvar; kk++)
            {
                aux = f[kk]*dS;
                solver->R[kk][e0] += aux;
            }
        } 
    }
}


void boundaryConvectiveSymmetry(BOUNDARY* boundary, SOLVER* solver)
{
    
    MESHBC* bc = boundary->bc;;
    
    int kk;
    double dSx, dSy, dS;
    double aux;
    double PL[6];
    double Pb[6];
    double f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int e0, p0, p1;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
 
        e0 = bc->elemL[ii]->neiL[0]->ii;
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<solver->Nvar; kk++)
	    {
	        if(kk > 3)
	        {
		        PL[kk] = solver->mesh->elemL[e0]->P[kk+1];
		    }
		    else
		    {
		        PL[kk] = solver->mesh->elemL[e0]->P[kk];			
		    }
	    }      		
        
        // Rotation of the velocity vectors
        rotation(PL, dSx, dSy, dS);
    
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            Pb[kk] = PL[kk];
        }
        Pb[1] *= -1;
    
        solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);
        
        // Rotation of the flux
	    rotation(f, dSx, -dSy, dS);
        
        if(dS > 0)
        {             
            for(kk=0; kk<solver->Nvar; kk++)
            {
                aux = f[kk]*dS;
                solver->R[kk][e0] += aux;
            }
        } 
    }
}


void boundary1(SOLVER* solver)
{

    for(int jj=0; jj<solver->mesh->Nmark; jj++)
    {
        MESHBC* bc = solver->mesh->bc[jj];
        
	    int kk;
        double dSx, dSy, dS;
        double aux;
        double PL[6];
	    double Pb[6];
        double f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        int e0, p0, p1;
        if(jj==2)
        {
            solver->boundaryL[jj]->convective(solver->boundaryL[jj], solver);
            if(0)//jj==2)
            {
                printf("\n%s\n", solver->boundaryL[jj]->type);
            }
        }
        else
        {

        for(int ii=0; ii<bc->Nelem; ii++)
        {
     
            e0 = bc->elemL[ii]->neiL[0]->ii;
            p0 = bc->elemL[ii]->p[0];
            p1 = bc->elemL[ii]->p[1];
     
            meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
            
            for(kk=0; kk<solver->Nvar; kk++)
		    {
		        if(kk > 3)
		        {
			        PL[kk] = solver->mesh->elemL[e0]->P[kk+1];
			    }
			    else
			    {
			        PL[kk] = solver->mesh->elemL[e0]->P[kk];			
			    }
		    }      		
            
            if(bc->flagBC == 0)
            {
                //Symmetry
                
                // Rotation of the velocity vectors
                rotation(PL, dSx, dSy, dS);
            
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    Pb[kk] = PL[kk];
                }
                Pb[1] *= -1;
            
                solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);

            }
            else if(bc->flagBC == 1)
            {
                //Inlet
                boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

                for(kk=0; kk<solver->Nvar; kk++)
	            {
	                if(kk > 3)
	                {
		                Pb[kk] = bc->elemL[ii]->P[kk+1];
		            }
		            else
		            {
		                Pb[kk] = bc->elemL[ii]->P[kk];
		            }
	            }

                // Rotation of the velocity vectors
                rotation(PL, dSx, dSy, dS);
	            rotation(Pb, dSx, dSy, dS);

                if(solver->flux1->extraVar)
                {
                    for(int kk=4; kk<solver->Nvar; kk++)
                    {
                        Pb[kk] = solver->inlet->Pin[kk+1];
                    }
                }
                
                solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);
            }
            else if(bc->flagBC == 2)
            {           
                //Outlet
                //boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);
                
                for(kk=0; kk<solver->Nvar; kk++)
	            {
	                if(kk > 3)
	                {
		                Pb[kk] = bc->elemL[ii]->P[kk+1];
		            }
		            else
		            {
		                Pb[kk] = bc->elemL[ii]->P[kk];
		            }
	            }

                // Rotation of the velocity vectors
                rotation(PL, dSx, dSy, dS);
	            rotation(Pb, dSx, dSy, dS);
            
                if(solver->flux1->extraVar)
                {
                    for(int kk=4; kk<solver->Nvar; kk++)
                    {
                        Pb[kk] = PL[kk];
                    }
                }
                
                solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);
            }
            else if((bc->flagBC == 3) || (bc->flagBC == 4))
            {            
                //Wall

                /*
                boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
                
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    f[kk] = 0;
                }
                f[1] = Pb[3];            
                */
                
                /*
                // Rotation of the velocity vectors
                rotation(PL, dSx, dSy, dS);
            
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    Pb[kk] = PL[kk];
                }
                Pb[1] *= -1;
                if(solver->sa || solver->sst->active)
                {
                    Pb[2] *= -1;
                }
            
                solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);
                */
                
                // Rotation of the velocity vectors
                //boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
                
                for(kk=0; kk<solver->Nvar; kk++)
	            {
	                if(kk > 3)
	                {
		                Pb[kk] = bc->elemL[ii]->P[kk+1];
		            }
		            else
		            {
		                Pb[kk] = bc->elemL[ii]->P[kk];
		            }
	            }
	            
                rotation(PL, dSx, dSy, dS);
	            rotation(Pb, dSx, dSy, dS);
            
                if(solver->flux1->extraVar)
                {
                    for(int kk=4; kk<solver->Nvar; kk++)
                    {
                        Pb[kk] = bc->elemL[ii]->P[kk+1];                    
                        //Pb[kk] = PL[kk];
                    }
                }
                
                solver->flux1->func(solver->flux1, solver->gas, PL, Pb, f);
                
            }       

            // Rotation of the flux
		    rotation(f, dSx, -dSy, dS);
            
            if(dS > 0)
            {             
                for(kk=0; kk<solver->Nvar; kk++)
                {
                    aux = f[kk]*dS;
                    solver->R[kk][e0] += aux;
                }
            } 
        }
        }
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
    else if(strcmp(s, "wallT") == 0)
    {
        ans = 4;
    }
    else
    {
        printf("Error: incorrent input of bc: %s\n", s);
        exit(0);
    }
 
    return ans;

}

void boundaryPrimitiveSymmetry(BOUNDARY* boundary, SOLVER* solver)
{
    MESHBC* bc = boundary->bc;

    # pragma omp parallel for
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        int kk;
        double dSx, dSy, dS;
        double PL[4];

        int p0, p1;
 
        ELEMENT* E0 = bc->elemL[ii]->neiL[0];
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E0->P[kk];
		}      		
        
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
        
        for(kk=0; kk<solver->Nvar+1; kk++)
        {
            bc->elemL[ii]->P[kk] = PL[kk];
        }
    }
}

void boundaryPrimitiveInlet(BOUNDARY* boundary, SOLVER* solver)
{
    MESHBC* bc = boundary->bc;

    # pragma omp parallel for
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        int kk;
        double dSx, dSy, dS;
        double PL[4];
	    double Pb[4];
        int p0, p1;    
        
        ELEMENT* E0 = bc->elemL[ii]->neiL[0];
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E0->P[kk];
		}      		
        
        boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

        for(kk=0; kk<4; kk++)
        {
            bc->elemL[ii]->P[kk] = Pb[kk];
        }
        
        bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);
        
        for(kk=5; kk<solver->Nvar+1; kk++)
        {
            bc->elemL[ii]->P[kk] = solver->inlet->Pin[kk];
        }
    }
}


void boundaryPrimitiveOutlet(BOUNDARY* boundary, SOLVER* solver)
{
    MESHBC* bc = boundary->bc;

    # pragma omp parallel for
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        int kk;
        double dSx, dSy, dS;
        double PL[4];
	    double Pb[4];
        int p0, p1;    

        ELEMENT* E0 = bc->elemL[ii]->neiL[0];
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E0->P[kk];
		}      		
        
        boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

        for(kk=0; kk<4; kk++)
        {
            bc->elemL[ii]->P[kk] = Pb[kk];
        }

        bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);
        
        for(kk=5; kk<solver->Nvar+1; kk++)
        {
            bc->elemL[ii]->P[kk] = E0->P[kk];
        }
    }
}


void boundaryPrimitiveWall(BOUNDARY* boundary, SOLVER* solver)
{
    MESHBC* bc = boundary->bc;

    # pragma omp parallel for
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        int kk;
        double dSx, dSy, dS;
        double PL[4];
	    double Pb[4];
        int p0, p1;    

        ELEMENT* E0 = bc->elemL[ii]->neiL[0];
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E0->P[kk];
		}      		
        
        if(solver->laminar==1 || solver->sa1->active || solver->sst->active)
        {
            boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
    
            bc->elemL[ii]->P[0] = Pb[0];
            bc->elemL[ii]->P[1] = 0.0;
            bc->elemL[ii]->P[2] = 0.0;
            bc->elemL[ii]->P[3] = Pb[3];
        }
        else
        {                              
            boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
            
            for(kk=0; kk<4; kk++)
            {
                bc->elemL[ii]->P[kk] = Pb[kk];
            } 
        }

        bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);
        
        if(solver->sa1->active)
        {            
            bc->elemL[ii]->P[5] = 0.0;
        }
        
        if(solver->sst->active)
        {
            double n = gaspropSutherland(bc->elemL[ii]->P[4])/bc->elemL[ii]->P[0];
            double d = solver->mesh->d[E0->ii];
            double owall = solver->sst->oWallFactor*6*n/(solver->sst->b1*d*d);
            bc->elemL[ii]->P[5] = 1e-14;
            bc->elemL[ii]->P[6] = owall;
            //printf("\n1,%e, %e, %e, %e", owall, E0->P[6], n, d);
        } 
    } 
}

void boundaryPrimitiveWallT(BOUNDARY* boundary, SOLVER* solver)
{

    MESHBC* bc = boundary->bc;
    # pragma omp parallel for
    for(int ii=0; ii<bc->Nelem; ii++)
    {
        int kk;
        double dSx, dSy, dS;
        double PL[4];
	    double Pb[4];
        int p0, p1;    
    
        ELEMENT* E0 = bc->elemL[ii]->neiL[0];
        p0 = bc->elemL[ii]->p[0];
        p1 = bc->elemL[ii]->p[1];
 
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        dS = sqrt(dSx*dSx + dSy*dSy);
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E0->P[kk];
		}      		
        
        if(solver->laminar==1 || solver->sa1->active || solver->sst->active)
        {
            boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
    
            bc->elemL[ii]->P[0] = Pb[3]/(solver->gas->R*solver->Twall);
            bc->elemL[ii]->P[1] = 0.0;
            bc->elemL[ii]->P[2] = 0.0;
            bc->elemL[ii]->P[3] = Pb[3];
            bc->elemL[ii]->P[4] = solver->Twall;

        }
        
        if(solver->sa1->active)
        {            
            bc->elemL[ii]->P[5] = 0.0;
        }
        
        if(solver->sst->active)
        {
            double n = gaspropSutherland(bc->elemL[ii]->P[4])/bc->elemL[ii]->P[0];
            double d = solver->mesh->d[E0->ii];
            double owall = solver->sst->oWallFactor*6*n/(solver->sst->b1*d*d);
            bc->elemL[ii]->P[5] = 1e-14;
            bc->elemL[ii]->P[6] = owall;
            //printf("\n1,%e, %e, %e, %e", owall, E0->P[6], n, d);
        } 
    } 
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
            
            bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);
        
        }
        else if(bc->flagBC == 1)
        {

            boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

            for(kk=0; kk<4; kk++)
            {
                bc->elemL[ii]->P[kk] = Pb[kk];
            }
            
            bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);            
        
        }
        else if(bc->flagBC == 2)
        {           

            boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

            for(kk=0; kk<4; kk++)
            {
                bc->elemL[ii]->P[kk] = Pb[kk];
            }
            
            bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);            
		
        }
        else if((bc->flagBC == 3) || (bc->flagBC == 4))
        {            
            if(solver->laminar==1 || solver->sa1->active || solver->sst->active)
            {
                boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);
            
                if(bc->flagBC == 4)
                {
                    bc->elemL[ii]->P[0] = Pb[3]/(solver->gas->R*solver->Twall);
                    bc->elemL[ii]->P[1] = 0.0;
                    bc->elemL[ii]->P[2] = 0.0;
                    bc->elemL[ii]->P[3] = Pb[3];
                    bc->elemL[ii]->P[4] = solver->Twall;
                }
                else
                {
                    bc->elemL[ii]->P[0] = Pb[0];
                    bc->elemL[ii]->P[1] = 0.0;
                    bc->elemL[ii]->P[2] = 0.0;
                    bc->elemL[ii]->P[3] = Pb[3];                
                    bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);
                }
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
                
                
                bc->elemL[ii]->P[4] = bc->elemL[ii]->P[3]/(bc->elemL[ii]->P[0]*solver->gas->R);
            }
        } 

        if(solver->sa1->active == 1)
        {
            if(bc->flagBC == 0)
            {
                //sym
                bc->elemL[ii]->P[5] = E0->P[5];
            }
            else if(bc->flagBC == 1)
            {        
                //inlet    
                bc->elemL[ii]->P[5] = solver->inlet->Pin[5];
            }
            else if(bc->flagBC == 2)
            {
                //out
                bc->elemL[ii]->P[5] = E0->P[5];
            }
            else if((bc->flagBC == 3) || (bc->flagBC == 4))
            {
                //wall wallT
                bc->elemL[ii]->P[5] = 0.0;
            }

        }
        
        if(solver->sst->active)
        {
            if(bc->flagBC == 0)
            {
                //sym
                bc->elemL[ii]->P[5] = E0->P[5];
                bc->elemL[ii]->P[6] = E0->P[6];
            }
            else if(bc->flagBC == 1)
            {        
                //inlet    
                bc->elemL[ii]->P[5] = solver->inlet->Pin[5];
                bc->elemL[ii]->P[6] = solver->inlet->Pin[6];                
            }
            else if(bc->flagBC == 2)
            {
                //out
                bc->elemL[ii]->P[5] = E0->P[5];
                bc->elemL[ii]->P[6] = E0->P[6];                
            }
            else if((bc->flagBC == 3) || (bc->flagBC == 4))
            {
                //wall wallT
                double n = gaspropSutherland(bc->elemL[ii]->P[4])/bc->elemL[ii]->P[0];
                double d = solver->mesh->d[E0->ii];
                double owall = solver->sst->oWallFactor*6*n/(solver->sst->b1*d*d);
                bc->elemL[ii]->P[5] = 1e-14;
                bc->elemL[ii]->P[6] = owall;
                //printf("\n1,%e, %e, %e, %e", owall, E0->P[6], n, d);
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
        double mi = gaspropSutherland(T);
        
        //printf("%f\n", E0->P[1]);
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx); 
        
        *fx = txx*dSx + txy*dSy;
		*fy = txy*dSx + tyy*dSy;

}

void boundaryViscous(BOUNDARY* boundary, SOLVER* solver)
{

    int e0;
    double f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MESHBC* bc = boundary->bc;
    double miEddy;

    for(int ii=0; ii<bc->Nelem; ii++)
    {
        e0 = bc->elemL[ii]->neiL[0]->ii;
 
        boundary->viscousFlux(boundary, solver, ii, f, &miEddy);
 
        for(int kk=1; kk<solver->Nvar; kk++)
        {
            solver->R[kk][e0] -= f[kk];
        }   
    }        
}

