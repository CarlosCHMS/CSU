#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"


typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[4];  

} CONDITION;

typedef struct{

    char* down;
    char* up;
    char* left;
    char* right;

    int Ndown;
    int Nup;
    int Nleft;
    int Nright;

} BOUNDARY;


typedef struct {

    int Nrow;
    int Ncol;
    int pOutFlag;
    int MUSCL;
    int flux;
    int stages;

    double Rgas;
    double gamma;
    double k4;
    double dt;
    double pout; 
    double eFix;       
    double e; 
    double k; 
    double res0[4];
    double res[4];
    double CFL;
            
    double **U;
    double **R;
    double **Uaux;        
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    BOUNDARY* bc;

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

char* inputGetValue2(INPUT* input, char* name)
{
    
    
    char* s;
    int found = 0;
    
    for(int ii=0; ii<input->N; ii++)
    {                
        if(strcmp(name, input->name[ii]) == 0)
        {
            s = input->value[ii];
            found = 1;
        }   
    }
    
    if(found == 0)
    {
        printf("Error: Value not found in the input file.\n");
        exit(0);
    }

    return s;
    
}

int main(int argc, char **argv)
{

    char s[50];
    char* ss;
    SOLVER* solver = malloc(sizeof(SOLVER));

    // Load input   
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "input.dat");
    INPUT* input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(input);

    //ss = inputGetValue2(input, "interpE");

    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.su2");
    solver->mesh = meshInit(s);    

    solver->U = tableMallocDouble(4, solver->mesh->Nelem);

    //printf("%s\n", inputGetValue(input, "interpE"));
    

    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
    
    if(inputNameIsInput(input, "CFL"))
    {
        solver->CFL = strtod(inputGetValue(input, "CFL"), NULL);     
    }
    else
    {
        solver->CFL = 1.0;
    }
    
    // Seletion of MUSCL and flux
    solver->MUSCL = atoi(inputGetValue(input, "order")) - 1;
    //solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));

    for(int ii=0; ii<solver->mesh->Nelem; ii++)
    {
        for(int jj=0; jj<4; jj++)
        {
            solver->U[jj][ii] = 0.0;
    
        }

    } 

    test(solver);

    //meshPrint(solver->mesh);
    
    // Save solution
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWrite(solver, s);
        
    solverFree(solver);
    inputFree(input);    

    return 0;

}
