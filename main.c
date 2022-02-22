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
    double* num= malloc(solver->mesh->Np*sizeof(double));
    double* den= malloc(solver->mesh->Np*sizeof(double));
    double aux;

    
    fprintf(ff, "0, %i, 4,\n", solver->mesh->Np);

    for(int ii; ii<solver->mesh->Np; ii++)
    {

        aux = cos(1.5*solver->mesh->p[ii][0])*cos(1.5*solver->mesh->p[ii][1]);
        fprintf(ff, "%f, ", aux);
        fprintf(ff, "%f, ", aux);
        fprintf(ff, "%f, ", aux);
        fprintf(ff, "%f, ", aux);
        fprintf(ff, "\n");        

    }

    fclose(ff);

    free(num);
    free(den);

}

int main()
{

    SOLVER* solver = malloc(sizeof(solver));

    solver->mesh = meshInit("./caseWedge/wedge.su2");

    solver->U = tableMallocDouble(4, solver->mesh->Nelem);

    

    meshPrint(solver->mesh);
    
    solverWrite(solver, "./caseWedge/solution.csv");
    
    solverFree(solver);

    return 0;

}
