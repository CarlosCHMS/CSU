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
#include"readTables.h"
#include"sa.h"

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

    double r, u, v, E, c, n;
    
    r = cond->p/(solver->Rgas*cond->T);
    c = sqrt(solver->gamma*solver->Rgas*cond->T);
    u = cond->nx*cond->mach*c;
    v = cond->ny*cond->mach*c;
    E = (solver->Rgas*cond->T)/(solver->gamma-1) + (u*u + v*v)/2;
    n = solver->turbRatio*sutherland(cond->T)/r;

    cond->Uin[0] = r;
    cond->Uin[1] = r*u;
    cond->Uin[2] = r*v;
    cond->Uin[3] = r*E;   

    cond->Pin[0] = r;
    cond->Pin[1] = u;
    cond->Pin[2] = v;
    cond->Pin[3] = solver->Rgas*cond->T*r;
    cond->Pin[4] = cond->T;

    if(solver->sa == 1)
    {
        cond->Uin[4] = r*n;
        cond->Pin[5] = n;         
    }

}

void solverMalloc(SOLVER* solver)
{
    if(solver->mallocType == 2)
    {
        newSolverMalloc2(solver);
    }
    else
    {
        newSolverMalloc(solver);    
    }
}

void solverWriteSolution(SOLVER* solver)
{

    /* 
        Based on: Adek Tasri, Accuracy of Cell Centres to Vertices 
        Interpolation for Unstructured Mesh Finite Volume Solver, 2021
    */

    // Save solution
    printf("main: saving the solution.\n");    

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "solution.csv");
    
    FILE* ff = fopen(fileName, "w");
    double** Up = tableMallocDouble(solver->Nvar, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    int ii, jj, kk, p;
    double xc, yc, xp, yp, L;
    ELEMENT* E;
    ELEMENT1* E1;

    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        for(kk=0; kk<solver->Nvar; kk++)
        {
            Up[kk][ii] = 0.;
        }   
        den[ii] = 0.;
    }

    for(ii=0; ii<solver->Ne; ii++)
    {
        E = solver->mesh->elemL[ii];
        E1 = solver->eL[ii];
        xc = E1->c[0];
        yc = E1->c[1];        
        
        for(jj=0; jj<E->Np; jj++)
        {
            p = E->p[jj];
            xp = solver->mesh->p[p][0];
            yp = solver->mesh->p[p][1];
            L = sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
           
            for(kk=0; kk<solver->Nvar; kk++)
            {        
                Up[kk][p] += E1->U[kk]/L;
            }
            
            den[p] += 1/L;
        }        
        
    }


    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        if(den[ii] != 0)
        {
            for(kk=0; kk<solver->Nvar; kk++)
            {
                Up[kk][ii] /= den[ii];    
            }
        }
    }

    fprintf(ff, "0, %i, %i,\n", solver->mesh->Np, solver->Nvar);

    for(ii=0; ii<solver->mesh->Np; ii++)
    {        
        for(jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", Up[jj][ii]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);

    tableFreeDouble(Up, solver->Nvar);
    free(den);

}

void solverResetR(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->eL[ii]->R[kk] = 0.0; 
        }
    }
}

void rotation(double* U, double dSx, double dSy, double dS)
{

    double aux = (U[1]*dSx + U[2]*dSy)/dS;
    U[2] = (-U[1]*dSy + U[2]*dSx)/dS;
    U[1] = aux;
	
}

void solverCalcRes(SOLVER* solver)
{
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        solver->res[kk] = fabs(solver->eL[0]->R[kk]);
    }
    
    for(int ii=1; ii<solver->Ne; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            if(solver->res[kk] < fabs(solver->eL[ii]->R[kk]))
            {
                solver->res[kk] = fabs(solver->eL[ii]->R[kk]);
            }
        }
    }
    
    for(int kk=0; kk<solver->Nvar; kk++)
    {
        printf(" %+.4e,", solver->res[kk]);        
    }
    printf("\n");      
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

double sutherland(double T)
{
    return 1.45e-6*T*sqrt(T)/(T + 110.0);
}

void solverSetData(SOLVER* solver, INPUT* input)
{

    solver->order = atoi(inputGetValue(input, "order"));

    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->Pr = 0.72;
    solver->Pr_t = 0.9;
    solver->Cp = solver->gamma*solver->Rgas/(solver->gamma-1);
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
    
    if(inputNameIsInput(input, "laminar"))
    {
        solver->laminar = atoi(inputGetValue(input, "laminar"));     
    }
    else
    {
        solver->laminar = 0;
    }

    if(inputNameIsInput(input, "restart"))
    {
        solver->restart = atoi(inputGetValue(input, "restart"));     
    }
    else
    {
        solver->restart = 0;
    }

    if(inputNameIsInput(input, "Sref"))
    {
        solver->Sref = strtod(inputGetValue(input, "Sref"), NULL);     
    }
    else
    {
        solver->Sref = 1.0;
    }

    if(inputNameIsInput(input, "turbRatio"))
    {
        solver->turbRatio = strtod(inputGetValue(input, "turbRatio"), NULL);     
    }
    else
    {
        solver->turbRatio = 10.0;
    }

    if(inputNameIsInput(input, "dtLocal"))
    {
        solver->dtLocal = atoi(inputGetValue(input, "dtLocal"));
    }
    else
    {
        solver->dtLocal = 0;
    }

    if(inputNameIsInput(input, "dtLocalN"))
    {
        solver->dtLocalN = atoi(inputGetValue(input, "dtLocalN"));
    }
    else
    {
        solver->dtLocalN = -1;
    }

    if(inputNameIsInput(input, "mallocType"))
    {
        solver->mallocType = atoi(inputGetValue(input, "mallocType"));
    }
    else
    {
        solver->mallocType = 2;
    }


    // Selection of several variables
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));
    solver->CFL = strtod(inputGetValue(input, "CFL"), NULL);

    // Environmental condition
    if(inputNameIsInput(input, "pout"))
    {
        solver->pout = strtod(inputGetValue(input, "pout"), NULL);     
    }
    
    // Set turbulence model
    if(inputNameIsInput(solver->input, "sa"))
    {
        solver->sa = atoi(inputGetValue(solver->input, "sa"));
    }
    else
    {
        solver->sa = 0;
    }

    // Set number of flow variables
    solver->Nvar = 4;
    if(solver->sa == 1)
    {
        solver->Nvar = 5;
    }
    
    solver->axi = atoi(inputGetValue(solver->input, "axisymmetric"));

}

SOLVER* solverInit(char* wd)
{
    SOLVER* solver = malloc(sizeof(SOLVER));
    char s[50];
    
    solver->Ndim = 2;

    // Work directory
    solver->wd = wd;

    // Load input   
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "input.dat");
    solver->input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(solver->input);

    // Setting the solver   
    solverSetData(solver, solver->input);

    // Set number of threads
    omp_set_num_threads(atoi(inputGetValue(solver->input, "threads")));

    // Load mesh    
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s, solver->Nvar, solver->axi);

    solver->mesh->order = solver->order;

    //Get boundary conditions
    printf("main: get boundary conditions.\n");
    boundaryGetBC(solver->mesh, solver->input);

      
    if(solver->sa == 1)
    {
        printf("main: calculating distance.\n");
        saCalcD(solver->mesh);
    }    
    
    //meshCheckNei(solver->mesh);
    //solverCheckGrad(solver);
    //meshPrint(solver->mesh);
    //meshPrintDStotal(solver->mesh);
    //meshCheckBorderOrientation(solver->mesh);

    // Memory allocation
    solverMalloc(solver);
    
    // Domain initialization
    //solverInitDomain(solver);
    
    // Get data from the mesh
    newSolverGetMeshData(solver);
     
    newSolverInitDomain(solver);
    
    return solver;
}

void newSolverMalloc(SOLVER* solver)
{

    solver->Ne = solver->mesh->Nelem;
    solver->Nf = solver->mesh->Ncon;

    solver->eL = malloc(solver->Ne*sizeof(ELEMENT1*));
    solver->fL = malloc(solver->Nf*sizeof(FACE1*));

    if(solver->axi)
    {
        solver->dSLateral = malloc(solver->Ne*sizeof(double));
    }
    
    if(solver->dtLocal)
    {
        solver->dtL = malloc(solver->Ne*sizeof(double));
    }
    
    if(solver->sa)
    {
        solver->dist = malloc(solver->Ne*sizeof(double));
    }
    
    // Allocation of the elements
    for(int ii=0; ii<solver->Ne; ii++)
    {
        solver->eL[ii] = malloc(sizeof(ELEMENT1));

        solver->eL[ii]->Nf = 0;
        solver->eL[ii]->fL = malloc(solver->mesh->elemL[ii]->Np*sizeof(FACE1*));
        solver->eL[ii]->nei = malloc(solver->mesh->elemL[ii]->Np*sizeof(bool));
        
        solver->eL[ii]->P = malloc((solver->Nvar+1)*sizeof(double));
        solver->eL[ii]->R = malloc(solver->Nvar*sizeof(double));        
        solver->eL[ii]->U = malloc(solver->Nvar*sizeof(double));
        solver->eL[ii]->Uaux = malloc(solver->Nvar*sizeof(double));

        solver->eL[ii]->c = malloc(solver->Ndim*sizeof(double));
        
        solver->eL[ii]->grad = malloc(solver->Nvar*sizeof(double*));        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            solver->eL[ii]->grad[jj] = malloc(solver->Ndim*sizeof(double));
        }
    }

    // Allocation of the faces
    for(int ii=0; ii<solver->Nf; ii++)
    {
        solver->fL[ii] = newSolverInitFace(solver);        
    }
        
    // Alocation of the boundary conditions
    solver->Nbc = solver->mesh->Nmark;
    solver->bcL = malloc(solver->Nbc*sizeof(SOLVERBC*));
    for(int ii=0; ii<solver->Nbc; ii++)
    {
        solver->bcL[ii] = malloc(sizeof(SOLVERBC));
        SOLVERBC* bc = solver->bcL[ii];

        bc->Ne = solver->mesh->bc[ii]->Nelem;
        bc->Nf = solver->mesh->bc[ii]->Nelem;
            
        bc->eL = malloc(bc->Ne*sizeof(ELEMENT1*));
        bc->fL = malloc(bc->Nf*sizeof(FACE1*));
        
        for(int jj=0; jj<bc->Ne; jj++)
        {
            bc->eL[jj] = malloc(sizeof(ELEMENT1));
        
            bc->eL[jj]->P = malloc((solver->Nvar+1)*sizeof(double));            
            bc->eL[jj]->fL = malloc(sizeof(FACE1*));
            bc->eL[jj]->nei = malloc(sizeof(bool));
            bc->eL[jj]->c = malloc(solver->Ndim*sizeof(double));
                        
            bc->fL[jj] = newSolverInitFace(solver);
        }                
    }
}

void newSolverMalloc2(SOLVER* solver)
{
    solver->Ne = solver->mesh->Nelem;
    solver->Nf = solver->mesh->Ncon;

    solver->eL = malloc(solver->Ne*sizeof(ELEMENT1*));
    solver->fL = malloc(solver->Nf*sizeof(FACE1*));

    if(solver->axi)
    {
        solver->dSLateral = malloc(solver->Ne*sizeof(double));
    }
    
    if(solver->dtLocal)
    {
        solver->dtL = malloc(solver->Ne*sizeof(double));
    }
    
    if(solver->sa)
    {
        solver->dist = malloc(solver->Ne*sizeof(double));
    }
        
    // Allocation of the elements
    for(int ii=0; ii<solver->Ne; ii++)
    {
        solver->eL[ii] = malloc(sizeof(ELEMENT1));        
        solver->eL[ii]->Nf = 0;
        solver->eL[ii]->nei = malloc(solver->mesh->elemL[ii]->Np*sizeof(bool));
        solver->eL[ii]->fL = malloc(solver->mesh->elemL[ii]->Np*sizeof(FACE1*));
    }

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        solver->eL[ii]->P = malloc((solver->Nvar+1)*sizeof(double));        
    }
    
    for(int ii=0; ii<solver->Ne; ii++)
    {        
        solver->eL[ii]->U = malloc(solver->Nvar*sizeof(double));
    }    
    
    for(int ii=0; ii<solver->Ne; ii++)
    {        
        solver->eL[ii]->Uaux = malloc(solver->Nvar*sizeof(double));
    }    

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        solver->eL[ii]->R = malloc(solver->Nvar*sizeof(double));
    }    

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        solver->eL[ii]->c = malloc(solver->Ndim*sizeof(double));
    }    

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        solver->eL[ii]->grad = malloc(solver->Nvar*sizeof(double*));
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            solver->eL[ii]->grad[jj] = malloc(solver->Ndim*sizeof(double));
        }
    }    

    // Allocation of the faces
    for(int ii=0; ii<solver->Nf; ii++)
    {
        solver->fL[ii] = malloc(sizeof(FACE1));
        solver->fL[ii]->eL = malloc(2*sizeof(ELEMENT1*));
    }

    for(int ii=0; ii<solver->Nf; ii++)
    {
        solver->fL[ii]->flux = malloc(solver->Nvar*sizeof(double));
    }
     
    for(int ii=0; ii<solver->Nf; ii++)
    {
        solver->fL[ii]->dS = malloc((solver->Ndim+1)*sizeof(double));
    }        
    
    for(int ii=0; ii<solver->Nf; ii++)
    {
        solver->fL[ii]->c = malloc(solver->Ndim*sizeof(double));
    }        
        
    // Alocation of the boundary conditions
    solver->Nbc = solver->mesh->Nmark;
    solver->bcL = malloc(solver->Nbc*sizeof(SOLVERBC*));
    for(int ii=0; ii<solver->Nbc; ii++)
    {
        solver->bcL[ii] = malloc(sizeof(SOLVERBC));
        SOLVERBC* bc = solver->bcL[ii];

        bc->Ne = solver->mesh->bc[ii]->Nelem;
        bc->Nf = solver->mesh->bc[ii]->Nelem;
            
        bc->eL = malloc(bc->Ne*sizeof(ELEMENT1*));
        bc->fL = malloc(bc->Nf*sizeof(FACE1*));
        
        // elements
        for(int jj=0; jj<bc->Ne; jj++)
        {
            bc->eL[jj] = malloc(sizeof(ELEMENT1));                
            bc->eL[jj]->fL = malloc(sizeof(FACE1*));
            bc->eL[jj]->nei = malloc(sizeof(bool));
        }

        for(int jj=0; jj<bc->Ne; jj++)
        {
            bc->eL[jj]->P = malloc((solver->Nvar+1)*sizeof(double));            
        }
        
        for(int jj=0; jj<bc->Ne; jj++)
        {
            bc->eL[jj]->c = malloc(solver->Ndim*sizeof(double)); 
        }
        
        // faces
        for(int jj=0; jj<bc->Nf; jj++)
        {
            bc->fL[jj] = malloc(sizeof(FACE1));
            bc->fL[jj]->eL = malloc(2*sizeof(ELEMENT1*));
        }
        
        for(int jj=0; jj<bc->Nf; jj++)
        {                    
            bc->fL[jj]->flux = malloc(solver->Nvar*sizeof(double));
        }

        for(int jj=0; jj<bc->Nf; jj++)
        {
            bc->fL[jj]->dS = malloc((solver->Ndim+1)*sizeof(double));
        }
        
        for(int jj=0; jj<bc->Nf; jj++)
        {        
            bc->fL[jj]->c = malloc(solver->Ndim*sizeof(double));
        }
    }
}


void newSolverFree(SOLVER* solver)
{
    meshFree(solver->mesh);
    
    inputFree(solver->input);

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        newSolverElementFree(solver, solver->eL[ii], true);
    }
    free(solver->eL);

    for(int jj=0; jj<solver->Nbc; jj++)
    {
        for(int ii=0; ii<solver->bcL[jj]->Ne; ii++)
        {        
            newSolverElementFree(solver, solver->bcL[jj]->eL[ii], false);
        }
        free(solver->bcL[jj]->eL);
    }
    
    for(int ii=0; ii<solver->Nf; ii++)
    {
        newSolverFaceFree(solver, solver->fL[ii]);
    }
    free(solver->fL);
    
    for(int jj=0; jj<solver->Nbc; jj++)
    {
        for(int ii=0; ii<solver->bcL[jj]->Nf; ii++)
        {        
            newSolverFaceFree(solver, solver->bcL[jj]->fL[ii]);
        }
        free(solver->bcL[jj]->fL);
    }    
    
    free(solver->bcL);
    
    if(solver->axi)
    {
        free(solver->dSLateral);
    }
    
    if(solver->dtLocal)
    {
        free(solver->dtL);
    }
    
    if(solver->sa)
    {
        free(solver->dist);
    }
}

void newSolverElementFree(SOLVER* solver, ELEMENT1* e, bool internal)
{
    free(e->P);
    if(internal)
    {
        free(e->U);
        free(e->Uaux);
        free(e->R);
    
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            free(e->grad[jj]);
        }
        
        free(e->grad);
    }
    
    free(e->fL);
    free(e->c);
    free(e->nei);
    free(e);
}

void newSolverFaceFree(SOLVER* solver, FACE1* f)
{
    free(f->flux);
    free(f->dS);
    free(f->c);
    free(f->eL);
    free(f);
}

void newSolverGetMeshData(SOLVER* solver)
{        
    int p0, p1;
    double x, y, dSx, dSy;
    ELEMENT1* e;    

    for(int ii=0; ii<solver->Ne; ii++)
    {   
        solver->eL[ii]->omega = solver->mesh->elemL[ii]->omega;
        elementCenter(solver->mesh->elemL[ii], solver->mesh, &x, &y);
        solver->eL[ii]->c[0] = x;
        solver->eL[ii]->c[1] = y;
        if(solver->sa)
        {
            solver->dist[ii] = solver->mesh->elemL[ii]->d;
        }
    }
    
    for(int ii=0; ii<solver->Nf; ii++)
    {
        p0 = solver->mesh->con[ii][2];
        p1 = solver->mesh->con[ii][3];
        meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
        solver->fL[ii]->dS[0] = dSx;
        solver->fL[ii]->dS[1] = dSy;  
        solver->fL[ii]->dS[2] = sqrt(dSx*dSx + dSy*dSy);
        solver->fL[ii]->c[0] = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
        solver->fL[ii]->c[1] = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;
        solver->fL[ii]->eL[0] = solver->eL[solver->mesh->con[ii][0]];
        solver->fL[ii]->eL[1] = solver->eL[solver->mesh->con[ii][1]];        
    }
    
    for(int ii=0; ii<solver->Nf; ii++)
    {
        e = solver->fL[ii]->eL[0];
        e->fL[e->Nf] = solver->fL[ii];
        e->nei[e->Nf] = false;
        e->Nf += 1;
        
        e = solver->fL[ii]->eL[1];
        e->fL[e->Nf] = solver->fL[ii];
        e->nei[e->Nf] = true;
        e->Nf += 1;        
    }  

    for(int ii=0; ii<solver->Nbc; ii++)
    {
        SOLVERBC* bc = solver->bcL[ii];
        bc->flagBC = solver->mesh->bc[ii]->flagBC;
        
        for(int jj=0; jj<bc->Ne; jj++)
        {                                    
            bc->fL[jj]->eL[0] = bc->eL[jj];
            bc->fL[jj]->eL[0]->fL[bc->fL[jj]->eL[0]->Nf] = bc->fL[jj];
            bc->fL[jj]->eL[0]->nei[bc->fL[jj]->eL[0]->Nf] = false;
            bc->fL[jj]->eL[0]->Nf += 1;
                        
            bc->fL[jj]->eL[1] = solver->eL[solver->mesh->bc[ii]->elemL[jj]->neiL[0]->ii];
            bc->fL[jj]->eL[1]->fL[bc->fL[jj]->eL[1]->Nf] = bc->fL[jj];
            bc->fL[jj]->eL[1]->nei[bc->fL[jj]->eL[1]->Nf] = true;
            bc->fL[jj]->eL[1]->Nf += 1;
            
            p0 = solver->mesh->bc[ii]->elemL[jj]->p[0];
            p1 = solver->mesh->bc[ii]->elemL[jj]->p[1];
            meshCalcDS(solver->mesh, p0, p1, &dSx, &dSy);
            bc->fL[jj]->dS[0] = dSx;
            bc->fL[jj]->dS[1] = dSy;
            bc->fL[jj]->dS[2] = sqrt(dSx*dSx + dSy*dSy);
            bc->fL[jj]->c[0] = (solver->mesh->p[p0][0] + solver->mesh->p[p1][0])*0.5;
            bc->fL[jj]->c[1] = (solver->mesh->p[p0][1] + solver->mesh->p[p1][1])*0.5;
            
            bc->eL[jj]->c[0] = bc->fL[jj]->c[0];
            bc->eL[jj]->c[1] = bc->fL[jj]->c[1];
        }                
    }
    
    if(solver->axi)
    {
        for(int ii=0; ii<solver->Ne; ii++)
        {
            solver->dSLateral[ii] = meshCalcDSlateral(solver->mesh, ii);
        }
    }

    solverClassifyFace(solver);
}

FACE1* newSolverInitFace(SOLVER* solver)
{
    FACE1* f = malloc(sizeof(FACE1));
        
    f->flux = malloc(solver->Nvar*sizeof(double));
    f->dS = malloc((solver->Ndim+1)*sizeof(double));
    f->c = malloc(solver->Ndim*sizeof(double));
    f->eL = malloc(2*sizeof(ELEMENT1*));
    
    return f;
}

void newSolverInitDomain(SOLVER* solver)
{
    char s[50];

    printf("main: initialize U.\n");       
    if(atoi(inputGetValue(solver->input, "tube")) == 0)
    {
        solver->inlet = conditionInit(strtod(inputGetValue(solver->input, "pressure"), NULL), 
                                      strtod(inputGetValue(solver->input, "temperature"), NULL), 
                                      strtod(inputGetValue(solver->input, "mach"), NULL), 
                                      strtod(inputGetValue(solver->input, "nx"), NULL),
                                      strtod(inputGetValue(solver->input, "ny"), NULL));

        conditionState(solver->inlet, solver);

        if(solver->restart)                                              
        {
            s[0] = '\0';
            strcat(s, solver->wd);
            strcat(s, "restart.csv");        
            newSolverLoadRestart(solver, s);
        }
        else
        {    
            newSolverInitU(solver, solver->inlet);
            if(solver->sa == 1)
            {
                saInitU(solver, solver->inlet);
            }
        }
    }
    else
    {
        if(solver->restart)                                              
        {
            s[0] = '\0';
            strcat(s, solver->wd);
            strcat(s, "restart.csv");        
            newSolverLoadRestart(solver, s);
        }
        else
        {        
            CONDITION* inside1 = conditionInit(strtod(inputGetValue(solver->input, "pressure1"), NULL), 
                                               strtod(inputGetValue(solver->input, "temperature1"), NULL), 
                                               strtod(inputGetValue(solver->input, "mach1"), NULL), 
                                               strtod(inputGetValue(solver->input, "nx1"), NULL),
                                               strtod(inputGetValue(solver->input, "ny1"), NULL));

            CONDITION* inside2 = conditionInit(strtod(inputGetValue(solver->input, "pressure2"), NULL), 
                                               strtod(inputGetValue(solver->input, "temperature2"), NULL), 
                                               strtod(inputGetValue(solver->input, "mach2"), NULL), 
                                               strtod(inputGetValue(solver->input, "nx2"), NULL),
                                               strtod(inputGetValue(solver->input, "ny2"), NULL));      
        
            newSolverInitUTube(solver, inside1, inside2, strtod(inputGetValue(solver->input, "xm"), NULL));
            free(inside1);
            free(inside2);
        }
    }
}

void newSolverInitU(SOLVER* solver, CONDITION* inside)
{
    conditionState(inside, solver);
    for(int ii=0; ii<solver->Ne; ii++)
    {
        for(int kk=0; kk<4; kk++)
        {
            solver->eL[ii]->U[kk] = inside->Uin[kk];
        }
    }
}

double newSolverLocalTimeStep(SOLVER* solver, ELEMENT1* e)
{
    double dSxm, dSym, Lx, Ly, u, v, c;

    dSxm = 0.;
    dSym = 0.;
    
    for(int jj=0; jj<e->Nf; jj++)
    {
        dSxm += fabs(e->fL[jj]->dS[0]);
        dSym += fabs(e->fL[jj]->dS[1]);   
    }
        
    dSxm *= 0.5;
    dSym *= 0.5;
    
    newSolverCalcVel(solver, e, &u, &v, &c);
    
    Lx = (fabs(u) + c)*dSxm;
    Ly = (fabs(v) + c)*dSym;
    
    return e->omega/(Lx + Ly);   
}

void newSolverCalcVel(SOLVER* solver, ELEMENT1* e, double* u, double* v, double* c)
{
    double E = e->U[3]/e->U[0];
    double aux;

    *u = e->U[1]/e->U[0];
    *v = e->U[2]/e->U[0];
    
    aux = (E - ((*u)*(*u) + (*v)*(*v))/2);
    aux *= solver->gamma - 1;
    *c = sqrt(aux*solver->gamma);    
}

void newSolverCheckOrientation(SOLVER* solver)
{
    int pos = 0;
    int neg = 0;
    double dx, dy, ans;
    FACE1* f;
    for(int ii=0; ii<solver->Nf; ii++)
    {
        f = solver->fL[ii];
        dx = f->eL[1]->c[0] - f->eL[0]->c[0];
        dy = f->eL[1]->c[1] - f->eL[0]->c[1];
        ans = f->dS[0]*dx + f->dS[1]*dy;
        if(ans > 0)
        {
            pos += 1;
        }
        else
        {
            neg += 1;
        }
    }
    
    for(int jj=0; jj<solver->Nbc; jj++)
    {
        SOLVERBC* bc = solver->bcL[jj];
        for(int ii=0; ii<bc->Nf; ii++)
        {
            f = bc->fL[ii];
            dx = f->eL[1]->c[0] - f->eL[0]->c[0];
            dy = f->eL[1]->c[1] - f->eL[0]->c[1];
            ans = f->dS[0]*dx + f->dS[1]*dy;
            if(ans > 0)
            {
                pos += 1;
            }
            else
            {
                neg += 1;
            }
        }
    }
    
    printf("\npos = %i, neg = %i\n", pos, neg);
}

void newBoundaryCalcPrimitive(SOLVER* solver, SOLVERBC* bc)
{    
    # pragma omp parallel for
    for(int ii=0; ii<bc->Ne; ii++)
    {
        double dSx, dSy, dS;
        double PL[4];
	    double Pb[4];
	    int kk;


        ELEMENT1* E0 = bc->eL[ii];
        ELEMENT1* E1 = E0->fL[0]->eL[1];
        dSx = E0->fL[0]->dS[0];
        dSy = E0->fL[0]->dS[1];        
        dS = E0->fL[0]->dS[2];
        
        for(kk=0; kk<4; kk++)
		{
			PL[kk] = E1->P[kk];
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
                E0->P[kk] = PL[kk];
            }
        }
        else if(bc->flagBC == 1)
        {
            boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

            for(kk=0; kk<4; kk++)
            {
                E0->P[kk] = Pb[kk];
            }
        }
        else if(bc->flagBC == 2)
        {           
            boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

            for(kk=0; kk<4; kk++)
            {
                E0->P[kk] = Pb[kk];
            }		
        }
        else if(bc->flagBC == 3)
        {  
            if(solver->laminar==1 || solver->sa==1)
            {
                E0->P[0] = PL[0];
                E0->P[1] = 0.0;
                E0->P[2] = 0.0;
                E0->P[3] = PL[3];
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
                        E0->P[kk] = Pb[kk];
                    }                    
                }
                else
                {
                    PL[2] = 0.0;                
                    for(kk=0; kk<4; kk++)
                    {
                        E0->P[kk] = PL[kk];
                    }                    
                }
            }
        }
        
        E0->P[4] = bc->eL[ii]->P[3]/(bc->eL[ii]->P[0]*solver->Rgas);         
        
        if(solver->sa == 1)
        {
            if(bc->flagBC == 0)
            {
                //sym
                E0->P[5] = E1->P[5];
            }
            else if(bc->flagBC == 1)
            {        
                //inlet    
                E0->P[5] = solver->inlet->Pin[5];
            }
            else if(bc->flagBC == 2)
            {
                //out
                E0->P[5] = E1->P[5];
            }
            else if(bc->flagBC == 3)
            {
                //wall
                E0->P[5] = 0.0;
            }

        }
    }
}


void newSolverUpdatePrimitives(SOLVER* solver)
{   
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {         
        ELEMENT1* E = solver->eL[ii];
        
        E->P[0] = E->U[0];
        E->P[1] = E->U[1]/E->U[0];
        E->P[2] = E->U[2]/E->U[0];
        E->P[3] = newSolverCalcP(solver, E);
        E->P[4] = E->P[3]/(E->P[0]*solver->Rgas);
        if(solver->sa)
        {
            E->P[5] = E->U[4]/E->U[0];
            //printf("\n%f, %f, %f", E->P[5], E->U[4], E->U[0]);            
        }
    }
    
    for(int ii=0; ii<solver->Nbc; ii++)
    {
        newBoundaryCalcPrimitive(solver, solver->bcL[ii]);
    }
}

double newSolverCalcP(SOLVER* solver, ELEMENT1* e)
{

	double u = e->U[1]/e->U[0];
	double v = e->U[2]/e->U[0];
    
	return (solver->gamma - 1)*(e->U[3] - 0.5*(u*u + v*v)*e->U[0]);
    
}

void newSolverInter(SOLVER* solver)
{

    newSolverInterFace(solver);

    for(int jj=0; jj<solver->Nbc; jj++)
    {
        if(solver->sa)
        {
            saBoundaryCalc(solver, solver->bcL[jj]);
        }
        else
        {
            newBoundaryCalc(solver, solver->bcL[jj]);
        }
    }    
 
    //newSolverDistributeFlux(solver);
}

double solverModule(double* v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1]);
}

void newSolverInterFace(SOLVER* solver)
{
    
    if(solver->order==2)
    {
        newSolverReconstruct(solver);
    }
    
    /*
    solverWriteGrad(solver);
    newSolverWriteElemBoundary(solver);
    newSolverWriteBCP(solver);
    exit(0);
    */

    # pragma omp parallel for
    for(int ii=0; ii<solver->Nf; ii++)
    {
        FACE1* f; 
        ELEMENT1* e0;
        ELEMENT1* e1;        
        int kk, mm;
        double PL[5];
        double PR[5];

        f = solver->fL[ii];
        e0 = f->eL[0];
        e1 = f->eL[1];        

        if(solver->order==2)
        {
            for(kk=0; kk<solver->Nvar; kk++)
            {
                if(kk == 4)
                {
                    mm = 5;
                }
                else
                {
                    mm = kk;
                }
	            PL[kk] = e0->P[mm] + e0->grad[kk][0]*(f->c[0] - e0->c[0]) + e0->grad[kk][1]*(f->c[1] - e0->c[1]);
	            PR[kk] = e1->P[mm] + e1->grad[kk][0]*(f->c[0] - e1->c[0]) + e1->grad[kk][1]*(f->c[1] - e1->c[1]);;
            }        
        }
        else
        {         
            for(kk=0; kk<solver->Nvar; kk++)
            {
                if(kk == 4)
                {
                    mm = 5;
                }
                else
                {
                    mm = kk;
                }
	            PL[kk] = e0->P[mm];
	            PR[kk] = e1->P[mm];
            }
        }

        // Rotation of the velocity vectors
	    rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
                    
        // Rotation of the velocity vectors
	    rotation(PR, f->dS[0], f->dS[1], f->dS[2]);
        
        // Flux calculation
        if(solver->sa)
        {
            fluxAUSMDV_sa(solver, PL[0], PL[1], PL[2], PL[3], PL[4], PR[0], PR[1], PR[2], PR[3], PR[4], f->flux);
        }
        else
        {
            flux(solver, PL[0], PL[1], PL[2], PL[3], PR[0], PR[1], PR[2], PR[3], f->flux);
        }

        // Rotation of the flux
	    rotation(f->flux, f->dS[0], -f->dS[1], f->dS[2]);
        
        for(kk=0; kk<solver->Nvar; kk++)
        {
            f->flux[kk] *= f->dS[2];
        }
    }
}

void newSolverCalcR(SOLVER* solver)
{
    solverResetR(solver);
    
    newSolverUpdatePrimitives(solver);   
    
    newSolverInter(solver);
    
    if(solver->axi)
    {
        solverInterAxisPressure(solver);
    }
    
    if(solver->laminar)
    {
        newSolverInterVisc(solver);
    }

    if(solver->sa)
    {
        saInter(solver);
    }

    newSolverDistributeFlux(solver);

    //solverWriteGrad(solver);
    //solverWriteFlux(solver);
    //solverWriteR(solver);
    //solverWriteP(solver);
    //solverWriteU(solver);
    //solverWriteFlux2(solver);
    //exit(0);
}

void newSolverStep(SOLVER* solver)
{
    newSolverUaux(solver);

    newSolverCalcR(solver); 
  
    newSolverEuler(solver);
}

void newSolverSolve(SOLVER* solver)
{       
    char s[50];

    // Calculate time step        
    int Nmax;

    // Convergence history file
    s[0] = '\0';
    strcat(s, solver->wd);
    strcat(s, "convergence.csv"); 
    FILE* convFile = fopen(s, "a");

    // Run the solver
    printf("\nmain: running solution:\n");
    if(atoi(inputGetValue(solver->input, "tube")) == 0)
    {        
        Nmax = atoi(inputGetValue(solver->input, "Nmax"));        
        for(int ii=0; ii<Nmax; ii++)
        {
            newSolverCalcDt(solver);
            newSolverStepRK(solver);
            
            if(ii == solver->dtLocalN)
            {
                solver->dtLocal = 0;
            }
                    
            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
                // Write convergence file
                fprintf(convFile, "%i,", ii);
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    fprintf(convFile, " %+.4e,", solver->res[kk]);        
                }
                newSolverCalcCoeff3(solver, convFile, ii);
                fprintf(convFile, "\n");
            }        
        }
    }
    else
    {
        double tmax = strtod(inputGetValue(solver->input, "tmax"), NULL);                

        // Run the solver
        double t = 0.0;        
        int stopLoop = 0;
        int ii = 0;
        while(stopLoop == 0)
        {
            newSolverCalcDt(solver);
            
            if(t + solver->dt>tmax)
            {
                solver->dt = (tmax-t);
                stopLoop = 1;
            }

            newSolverStepRK(solver);
            t += solver->dt;
            ii++;
            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
            }
        }
        printf("time %f s\n", t);
    }
    fclose(convFile);
}

void newSolverEuler(SOLVER* solver)
{   
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {            
            solver->eL[ii]->U[kk] = solver->eL[ii]->Uaux[kk] - solver->dt*solver->eL[ii]->R[kk]/solver->eL[ii]->omega;
        }
    }
}

void newSolverUaux(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        for(int kk=0; kk<solver->Nvar; kk++)
        {
            solver->eL[ii]->Uaux[kk] = solver->eL[ii]->U[kk];
        }
    }
}

void checkNei(SOLVER* solver)
{
    for(int ii=0; ii<solver->Ne; ii++)
    {
        printf("\n%i, %i, %i, %i", ii, solver->eL[ii]->nei[0], solver->eL[ii]->nei[1], solver->eL[ii]->nei[2]);
    }
}

void newSolverCalcDt(SOLVER* solver)
{
    double dt;
    double dtLocal;
    
    if(solver->dtLocal)
    {
        # pragma omp parallel for
        for(int ii=0; ii<solver->Ne; ii++)
        {
            solver->dtL[ii] = (solver->CFL*0.5*solver->stages)*newSolverLocalTimeStep(solver, solver->eL[ii]);
            if(solver->stages == 1)
            {
                solver->dtL[ii] *= 2;
            }
        }
    }
    else
    {
        dt = newSolverLocalTimeStep(solver, solver->eL[0]);
        for(int ii=1; ii<solver->Ne; ii++)
        {
            dtLocal = newSolverLocalTimeStep(solver, solver->eL[ii]);
            if(dtLocal < dt)
            {
                dt = dtLocal;
            }

        }

        dt *= solver->CFL*0.5*solver->stages;
        if(solver->stages == 1)
        {
            dt *= 2;
        }
        solver->dt = dt;
    }
    
}

void newBoundaryCalc(SOLVER* solver, SOLVERBC* bc)
{  
    # pragma omp parallel for
    for(int ii=0; ii<bc->Nf; ii++)
    { 
        int kk;
        double dSx, dSy, dS;
        double PL[4];
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

            if(bc->flagBC == 0)
            {

                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
            
                // Reflexive
                flux(solver, PL[0], PL[1], PL[2], PL[3], PL[0], -PL[1], PL[2], PL[3], f->flux);

            }
            else if(bc->flagBC == 1)
            {

                boundaryInlet(solver, solver->inlet->Pin, PL, Pb, dSx/dS, dSy/dS);

                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
	            rotation(Pb, f->dS[0], f->dS[1], f->dS[2]);

                flux(solver, PL[0], PL[1], PL[2], PL[3], Pb[0], Pb[1], Pb[2], Pb[3], f->flux);
            
            }
            else if(bc->flagBC == 2)
            {           

                boundaryOutlet(solver, PL, Pb, dSx/dS, dSy/dS);

                // Rotation of the velocity vectors
                rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
	            rotation(Pb, f->dS[0], f->dS[1], f->dS[2]);

                flux(solver, PL[0], PL[1], PL[2], PL[3], Pb[0], Pb[1], Pb[2], Pb[3], f->flux);
		    
            }
            else if(bc->flagBC == 3)
            {
                if(solver->laminar==1)
                {
                    // Rotation of the velocity vectors
                    rotation(PL, dSx, dSy, dS);
            
                    // Reflexive
                    flux(solver, PL[0], PL[1], PL[2], PL[3], PL[0], -PL[1], -PL[2], PL[3], f->flux);
                }
                else
                { 
                    boundaryWall(solver, PL, Pb, dSx/dS, dSy/dS);

                    // Rotation of the velocity vectors
                    rotation(PL, f->dS[0], f->dS[1], f->dS[2]);
	                rotation(Pb, f->dS[0], f->dS[1], f->dS[2]);
                
                    //flux(solver, PL[0], PL[1], PL[2], PL[3], Pb[0], Pb[1], Pb[2], Pb[3], f);
                    fluxFree(solver, Pb[0], Pb[1], Pb[2], Pb[3], f->flux);
                }
            }       

            // Rotation of the flux
		    rotation(f->flux, f->dS[0], -f->dS[1], f->dS[2]);
		                        
            for(kk=0; kk<4; kk++)
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

void solverWriteR(SOLVER* solver)
{

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../R.dat");

    FILE* ff = fopen(fileName, "w");

    // Save reestart
    printf("main: saving the restart file.\n");    

    fprintf(ff, "1,\n");

    fprintf(ff, "0, %i, %i,\n", solver->Ne, solver->Nvar);

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", solver->eL[ii]->R[jj]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);
}

void solverWriteFlux(SOLVER* solver)
{
    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../flux.dat");

    FILE* ff = fopen(fileName, "w");

    fprintf(ff, "1,\n");

    fprintf(ff, "0, %i, %i,\n", solver->Nf, solver->Nvar);

    for(int ii=0; ii<solver->Nf; ii++)
    {        
        for(int jj=1; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", -solver->fL[ii]->flux[jj]);
        }
        fprintf(ff, "\n");        
    }
    fclose(ff);
}

void solverWriteFlux2(SOLVER* solver)
{
    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../fluxB.dat");

    FILE* ff = fopen(fileName, "w");

    for(int kk=0; kk<solver->Nbc; kk++)
    {
        for(int ii=0; ii<solver->bcL[kk]->Nf; ii++)
        {        
            for(int jj=1; jj<solver->Nvar; jj++)
            {
                fprintf(ff, "%.10e, ", solver->bcL[kk]->fL[ii]->flux[jj]);
            }
            fprintf(ff, "\n");        
        }
    }
    fclose(ff);
}


void solverWriteFluxBC(SOLVER* solver, double* f)
{
    for(int jj=0; jj<solver->Nvar; jj++)
    {
        fprintf(solver->fb, "%.10e, ", f[jj]);
    }
    fprintf(solver->fb, "\n");        
}

void solverClassifyFace(SOLVER* solver)
{
    for(int ii=0; ii<solver->Nf; ii++)
    {
        solver->fL[ii]->b= false;
    }

    for(int jj=0; jj<solver->Nbc; jj++)
    {
        for(int ii=0; ii<solver->bcL[jj]->Nf; ii++)
        {
            solver->bcL[jj]->fL[ii]->b = true;
        }
    }
}

void newSolverRK(SOLVER* solver, double a)
{   
    if(solver->dtLocal)
    {
        # pragma omp parallel for
        for(int ii=0; ii<solver->Ne; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {            
                solver->eL[ii]->U[kk] = solver->eL[ii]->Uaux[kk] - solver->dtL[ii]*a*solver->eL[ii]->R[kk]/solver->eL[ii]->omega;
            }
        }
    }
    else
    {
        # pragma omp parallel for
        for(int ii=0; ii<solver->Ne; ii++)
        {
            for(int kk=0; kk<solver->Nvar; kk++)
            {            
                solver->eL[ii]->U[kk] = solver->eL[ii]->Uaux[kk] - solver->dt*a*solver->eL[ii]->R[kk]/solver->eL[ii]->omega;
            }
        }    
    }
}

void newSolverStepRK(SOLVER* solver)
{   

    newSolverUaux(solver);

    if(solver->stages==1)
    {
        newSolverCalcR(solver);
        newSolverRK(solver, 1.0);
    }
    else if(solver->stages==3)
    {
        newSolverCalcR(solver);
        newSolverRK(solver, 0.1481);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 0.4);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 1.0);
    }
    else if(solver->stages==4)
    {
        newSolverCalcR(solver);
        newSolverRK(solver, 0.0833);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 0.2069);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 0.4265);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 1.0);
    }
    else if(solver->stages==5)
    {
        newSolverCalcR(solver);
        newSolverRK(solver, 0.0533);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 0.1263);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 0.2375);

        newSolverCalcR(solver);
        newSolverRK(solver, 0.4414);
        
        newSolverCalcR(solver);
        newSolverRK(solver, 1.0);
    }

    //solverWriteU(solver);
    //exit(0);
}

void newSolverCalcCoeff3(SOLVER* solver, FILE* convFile, int Nint)
{
    SOLVERBC* bc;
    double cp, fx, fy;    
    
    double r = solver->inlet->Pin[0];
    double u = solver->inlet->Pin[1];
    double v = solver->inlet->Pin[2];
    double P = solver->inlet->Pin[3];
    double q = 0.5*r*(u*u + v*v);
 
    double Cx_p = 0;
    double Cx_v = 0;       
        
    double Cy_p = 0;
    double Cy_v = 0;               
        
    for(int jj=0; jj<solver->Nbc; jj++)
    {
        bc = solver->bcL[jj];
        if(bc->flagBC == 3)
        {
            for(int ii=0; ii<bc->Ne; ii++)
            {
                cp = (bc->eL[ii]->P[3] - P)/q;
                Cx_p += cp*bc->eL[ii]->fL[0]->dS[0];
                Cy_p += cp*bc->eL[ii]->fL[0]->dS[1];
                
                if(solver->laminar==1 || solver->sa==1)
                {
                    newBoundaryCalcFrictionWall(solver, bc->eL[ii]->fL[0], &fx, &fy);
                    Cx_v -= fx/q;
                    Cy_v -= fy/q;
                }
                                              
            }
        }
    }
        
    Cx_p /= solver->Sref;
    Cx_v /= solver->Sref;    

    Cy_p /= solver->Sref;
    Cy_v /= solver->Sref;    

    if(solver->mesh->axi)
    {
        Cx_p *= 2*M_PI;
        Cx_v *= 2*M_PI;    
    }

    fprintf(convFile, " %e, %e, %e, %e,", Cx_p, Cx_v, Cy_p, Cy_v); 

}

void newSolverCalcGrad2(ELEMENT1* E, int kk, int grad_ii, double* Umin, double* Umax)
{
    /*
    Least squares calculation
    */

    double x[5];
    double y[5];    
    double u[5];    
    double a, b, c, A, B;
    ELEMENT1* E1;

    x[0] = E->c[0];
    y[0] = E->c[1];
    u[0] = E->P[kk];

    for(int jj=0; jj<E->Nf; jj++)
    {
        if(E->nei[jj])
        {
            E1 = E->fL[jj]->eL[0];
        }
        else
        {
            E1 = E->fL[jj]->eL[1];
        }
    
        x[jj+1] = E1->c[0];
        y[jj+1] = E1->c[1];
        u[jj+1] = E1->P[kk];
    }
    
    a = 0;    
    b = 0;
    c = 0;
    A = 0;
    B = 0;
    for(int jj=1; jj<E->Nf+1; jj++)
    {
        a += (x[jj] - x[0])*(x[jj] - x[0]);
        b += (x[jj] - x[0])*(y[jj] - y[0]);
        c += (y[jj] - y[0])*(y[jj] - y[0]);
        A += (x[jj] - x[0])*(u[jj] - u[0]);
        B += (y[jj] - y[0])*(u[jj] - u[0]);
    }
        
    E->grad[grad_ii][0] = (c*A - b*B)/(a*c - b*b);
    E->grad[grad_ii][1] = (-b*A + a*B)/(a*c - b*b);
    
    *Umin = u[0];
    *Umax = u[0];
    
    for(int ii=1; ii<E->Nf+1; ii++)
    {
        *Umax = fmax(*Umax, u[ii]);
        *Umin = fmin(*Umin, u[ii]);        
    }
}

void newSolverReconstruct(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        int jj, kk, mm;
        double d2, phi0, phi;
        double Pmin, Pmax;        
        ELEMENT1* E = solver->eL[ii];
  	    
        for(kk=0; kk<solver->Nvar; kk++)
        {      
              
            if(kk == 4)
            {
                mm = 5;
            }
            else
            {
                mm = kk;
            }

            newSolverCalcGrad2(E, mm, kk, &Pmin, &Pmax);       
            for(jj=0; jj<E->Nf; jj++)
            {         
                d2 = (E->grad[kk][0]*(E->fL[jj]->c[0] - E->c[0]) + E->grad[kk][1]*(E->fL[jj]->c[1] - E->c[1]));
                phi0 = limiterV(E->P[mm], Pmin, Pmax, d2, solver->e);   
                
                if(jj==0)
                {
                    phi = phi0;
                }
                else
                {
                    phi = fmin(phi, phi0);
                }			            
            }

            E->grad[kk][0] *= phi;
            E->grad[kk][1] *= phi;
        }   
    }
    
    //newSolverWriteBCP(solver);
    //solverWriteGrad(solver);
    //exit(0);
}

void solverWriteGrad(SOLVER* solver)
{

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../grad.dat");

    FILE* ff = fopen(fileName, "w");

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", solver->eL[ii]->grad[jj][0]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);
}

void newSolverLoadRestart(SOLVER* solver, char* fileName)
{
    TABLELIST* tl = fReadTables(fileName);
    
    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            solver->eL[ii]->U[jj] = tl->tables[0]->values[ii][jj];
        }        
    }    
    readTablesFree(tl);

}

void newSolverWriteRestart(SOLVER* solver)
{
    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "restart.csv");

    FILE* ff = fopen(fileName, "w");

    // Save reestart
    printf("main: saving the restart file.\n");    

    fprintf(ff, "1,\n");

    fprintf(ff, "0, %i, %i,\n", solver->mesh->Nelem, solver->Nvar);

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", solver->eL[ii]->U[jj]);
        }
        fprintf(ff, "\n");        
    }
    fclose(ff);
}

void newSolverWriteElemBoundary(SOLVER* solver)
{
    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "eb.dat");

    FILE* ff = fopen(fileName, "w");

    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->eL[ii]->Nf; jj++)
        {
            if(solver->eL[ii]->fL[jj]->b)
            {
                fprintf(ff, "%i\n", ii);
                break;
            }
        } 
    }
    fclose(ff);
}

void newSolverWriteBCP(SOLVER* solver)
{
    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../BCP.csv");

    FILE* ff = fopen(fileName, "w");

    for(int kk=0; kk<solver->Nbc; kk++)
    {

        for(int ii=0; ii<solver->bcL[kk]->Ne; ii++)
        {        
            for(int jj=0; jj<solver->Nvar+1; jj++)
            {
                fprintf(ff, "%.10e, ", solver->bcL[kk]->eL[ii]->P[jj]);
            }
            fprintf(ff, "\n");        
        }
    }
    fclose(ff);
}

void solverWriteU(SOLVER* solver)
{

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../U.dat");

    FILE* ff = fopen(fileName, "w");
    
    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->Nvar; jj++)
        {
            fprintf(ff, "%.10e, ", solver->eL[ii]->U[jj]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);
}

void solverWriteP(SOLVER* solver)
{

    char fileName[50];
    fileName[0] = '\0';
    strcat(fileName, solver->wd);
    strcat(fileName, "../../P.dat");

    FILE* ff = fopen(fileName, "w");
    
    for(int ii=0; ii<solver->Ne; ii++)
    {        
        for(int jj=0; jj<solver->Nvar+1; jj++)
        {
            fprintf(ff, "%.10e, ", solver->eL[ii]->P[jj]);
        }
        fprintf(ff, "\n");        
    }

    fclose(ff);
}


void solverInterAxisPressure(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        solver->eL[ii]->R[2] -= solver->eL[ii]->P[3]*solver->dSLateral[ii];
    }
}

void newSolverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm)
{
    conditionState(inside1, solver);
    conditionState(inside2, solver);

    for(int ii=0; ii<solver->Ne; ii++)
    {
        if(solver->eL[ii]->c[0]<xm)
        {
            for(int kk=0; kk<4; kk++)
            {
                solver->eL[ii]->U[kk] = inside1->Uin[kk];
            }
        }
        else
        {
           for(int kk=0; kk<4; kk++)
            {
                solver->eL[ii]->U[kk] = inside2->Uin[kk];
            }
        }
    }
}

void newSolverInterVisc(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        double Pmin, Pmax;
        ELEMENT1* E = solver->eL[ii];

        newSolverCalcGrad2(E, 1, 1, &Pmin, &Pmax);
        
        newSolverCalcGrad2(E, 2, 2, &Pmin, &Pmax);
        // grad p storage is being reused for T grad
        newSolverCalcGrad2(E, 4, 3, &Pmin, &Pmax);
    }
    
    //solverWriteGrad(solver);
    //exit(0);
        
    # pragma omp parallel for
    for(int ii=0; ii<solver->Nf; ii++)
    {
        double F = -1;
                
        ELEMENT1* E0 = solver->fL[ii]->eL[0];
        ELEMENT1* E1 = solver->fL[ii]->eL[1];
                    
        double duxm = (E0->grad[1][0] + E1->grad[1][0])*0.5;
        double dvxm = (E0->grad[2][0] + E1->grad[2][0])*0.5;
        double dTxm = (E0->grad[3][0] + E1->grad[3][0])*0.5;
        
        double duym = (E0->grad[1][1] + E1->grad[1][1])*0.5;
        double dvym = (E0->grad[2][1] + E1->grad[2][1])*0.5;
        double dTym = (E0->grad[3][1] + E1->grad[3][1])*0.5;
        	        
        double dx = E1->c[0] - E0->c[0];
        double dy = E1->c[1] - E0->c[1];
        double L = sqrt(dx*dx + dy*dy);
        
        double dul = (E1->P[1] - E0->P[1])/L;
        double dvl = (E1->P[2] - E0->P[2])/L;
        double dTl = (E1->P[4] - E0->P[4])/L;
        
        double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

        double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;        

        double dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
        double dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;        
	    
        double T = (E1->P[4] + E0->P[4])*0.5;
        double mi = sutherland(T);
        double k = solver->Cp*mi/solver->Pr;
	    		        
	    double txx = 2*mi*(dux - (dux + dvy)/3);
	    double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	    double txy = mi*(duy + dvx);  
	    
	    double u = (E1->P[1] + E0->P[1])*0.5;
        double v = (E1->P[2] + E0->P[2])*0.5;
    
        double dSx = solver->fL[ii]->dS[0];
        double dSy = solver->fL[ii]->dS[1];        
        
	    solver->fL[ii]->flux[1] += F*(txx*dSx + txy*dSy);
	    solver->fL[ii]->flux[2] += F*(txy*dSx + tyy*dSy);
	    solver->fL[ii]->flux[3] += F*((txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy));
    }
    
    
    for(int ii=0; ii<solver->Nbc; ii++)
    {
        newSolverBoundaryCalcVisc(solver, solver->bcL[ii]);
    }
    
    //newSolverDistributeFlux(solver);
}

void newSolverDistributeFlux(SOLVER* solver)
{
    # pragma omp parallel for
    for(int ii=0; ii<solver->Ne; ii++)
    {
        for(int jj=0; jj<solver->eL[ii]->Nf; jj++)
        {
            if(solver->eL[ii]->nei[jj])
            {
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    solver->eL[ii]->R[kk] -= solver->eL[ii]->fL[jj]->flux[kk];
                }
            }
            else
            {
                for(int kk=0; kk<solver->Nvar; kk++)
                {
                    solver->eL[ii]->R[kk] += solver->eL[ii]->fL[jj]->flux[kk];
                }
            }
        }
    }
}

void newSolverBoundaryCalcVisc(SOLVER* solver, SOLVERBC* bc)
{
    # pragma omp parallel for
    for(int ii=0; ii<bc->Nf; ii++)
    {
        double dux, duy, dvx, dvy, dTx, dTy;
        
        FACE1* f = bc->fL[ii];
        ELEMENT1* E0 = f->eL[0];
        ELEMENT1* E1 = f->eL[1];

        double dSx = bc->fL[ii]->dS[0];
        double dSy = bc->fL[ii]->dS[1];        

        double F = 1;

        if(bc->flagBC == 0)
        {
           	double nx, ny, dS;
           	dS = bc->fL[ii]->dS[2];
            nx = dSx/dS;
            ny = dSy/dS;
            
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
            
            double T = E1->P[4];
            double mi = sutherland(T);
            double k = solver->Cp*mi/solver->Pr;            
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        double u = E1->P[1];
	        double v = E1->P[2];		        
	        
	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);
	        f->flux[3] += F*((txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy));       
        }
        else if(bc->flagBC == 1)
        {
            double dx = E0->c[0] - E1->c[0];		    
            double dy = E0->c[1] - E1->c[1];
            double L = sqrt(dx*dx + dy*dy);
            
            double dul = (solver->inlet->Pin[1] - E1->P[1])/L;
            double dvl = (solver->inlet->Pin[2] - E1->P[2])/L;            
            double dTl = (solver->inlet->Pin[4] - E1->P[4])/L; 

            double duxm = E1->grad[1][0];
            double dvxm = E1->grad[2][0];
            double dTxm = E1->grad[3][0];
            
            double duym = E1->grad[1][1];
            double dvym = E1->grad[2][1];
            double dTym = E1->grad[3][1];

            dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

            dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;  
            
            dTx = dTxm + (dTl - (dTxm*dx + dTym*dy)/L)*dx/L;
            dTy = dTym + (dTl - (dTxm*dx + dTym*dy)/L)*dy/L;
            
            double T = E1->P[4];
            double mi = sutherland(T);
            double k = solver->Cp*mi/solver->Pr;            
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        double u = E1->P[1];
	        double v = E1->P[2];		        
	        
	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);
	        f->flux[3] += F*((txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy));            
	                
        }
        else if(bc->flagBC == 3)
        {
            double dx = E1->c[0] - E0->c[0];		    
            double dy = E1->c[1] - E0->c[1];
            double L = sqrt(dx*dx + dy*dy);
            
            double dul = (E1->P[1] - 0)/L;
            double dvl = (E1->P[2] - 0)/L;            

            double duxm = E1->grad[1][0];
            double dvxm = E1->grad[2][0];
            
            double duym = E1->grad[1][1];
            double dvym = E1->grad[2][1];

            dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
            duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

            dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
            dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;
                            
            double T = E1->P[4];
            double mi = sutherland(T);
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);

        }
        else
        {        
            double dux = E1->grad[1][0];
            double dvx = E1->grad[2][0];
            double dTx = E1->grad[3][0];
            
            double duy = E1->grad[1][1];
            double dvy = E1->grad[2][1];
            double dTy = E1->grad[3][1];

            double T = E1->P[4];
            double mi = sutherland(T);
            double k = solver->Cp*mi/solver->Pr;            
	            
	        double txx = 2*mi*(dux - (dux + dvy)/3);
	        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
	        double txy = mi*(duy + dvx);  
	        
	        double u = E1->P[1];
	        double v = E1->P[2];		        
	        
	        f->flux[1] += F*(txx*dSx + txy*dSy);
	        f->flux[2] += F*(txy*dSx + tyy*dSy);
	        f->flux[3] += F*((txx*dSx + txy*dSy)*u + (txy*dSx + tyy*dSy)*v + k*(dTx*dSx + dTy*dSy));            
	        
        }
	}	
}

void newBoundaryCalcFrictionWall(SOLVER* solver, FACE1* f, double* fx, double* fy)
{
        ELEMENT1* E0 = f->eL[0];
        ELEMENT1* E1 = f->eL[1];

        double dSx = f->dS[0];
        double dSy = f->dS[1];        

        double dx = E1->c[0] - E0->c[0];		    
        double dy = E1->c[1] - E0->c[1];
        double L = sqrt(dx*dx + dy*dy);
        
        double dul = (E1->P[1] - 0.0)/L;
        double dvl = (E1->P[2] - 0.0)/L;            

        double duxm = E1->grad[1][0];
        double dvxm = E1->grad[2][0];
        
        double duym = E1->grad[1][1];
        double dvym = E1->grad[2][1];

        double dux = duxm + (dul - (duxm*dx + duym*dy)/L)*dx/L;
        double duy = duym + (dul - (duxm*dx + duym*dy)/L)*dy/L;        

        double dvx = dvxm + (dvl - (dvxm*dx + dvym*dy)/L)*dx/L;
        double dvy = dvym + (dvl - (dvxm*dx + dvym*dy)/L)*dy/L;
                        
        double T = E1->P[4];
        double mi = sutherland(T);
            
        double txx = 2*mi*(dux - (dux + dvy)/3);
        double tyy = 2*mi*(dvy - (dux + dvy)/3);		    
        double txy = mi*(duy + dvx); 
        
        *fx = txx*dSx + txy*dSy;
		*fy = txy*dSx + tyy*dSy;

}

