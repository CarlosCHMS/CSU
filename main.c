#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"utils.h"
#include"mesh.h"

typedef struct
{

    MESH* mesh;
    double** U;

} SOLVER;


void solverFree(SOLVER* solver)
{

    meshFree(solver->mesh);
    tableFreeDouble(solver->U, 4);

}

void solverWrite(SOLVER* solver, char* fileName)
{

    FILE* ff = fopen(fileName, "w");
    double** Up = tableMallocDouble(4, solver->mesh->Np);
    double* den = malloc(solver->mesh->Np*sizeof(double));
    double aux;
    int ii, jj;

    for(ii=0; ii<solver->mesh->Np; ii++)
    {       
        for(jj=0; jj<4; jj++)
        {
            Up[jj][ii] = 0.;
    
        }   
        den[ii] = 0.;
    }

    for(ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(jj=0; jj<4; jj++)
        {
            Up[jj][solver->mesh->elem[ii][0]] += solver->U[jj][ii];
            Up[jj][solver->mesh->elem[ii][1]] += solver->U[jj][ii];
            Up[jj][solver->mesh->elem[ii][2]] += solver->U[jj][ii];
    
        }

        den[solver->mesh->elem[ii][0]] += 1;
        den[solver->mesh->elem[ii][1]] += 1;
        den[solver->mesh->elem[ii][2]] += 1;
    }

    for(ii=0; ii<solver->mesh->Np; ii++)
    {
        for(jj=0; jj<4; jj++)
        {
            Up[jj][ii] /= den[ii] ;
    
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

void test(SOLVER* solver)
{

    double x, y, aux;

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {   
        meshElemCenter(solver->mesh, ii, &x, &y);
        solver->U[0][ii] = cos(1.5*x)*cos(1.5*y);
    }    

}

int main(int argc, char **argv)
{

    char s[50];

    SOLVER* solver = malloc(sizeof(solver));

    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s);    

    solver->U = tableMallocDouble(4, solver->mesh->Nelem);

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<4; jj++)
        {
            solver->U[jj][ii] = 0.0;
    
        }

    } 

    test(solver);

    meshPrint(solver->mesh);
    
    // Save solution
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWrite(solver, s);
    
    solverFree(solver);

    return 0;

}
