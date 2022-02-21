#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"input.h"
#include"mesh.h"

int main(int argc, char **argv)
{

    MESH* mesh;
    mesh = meshInit("./caseWedge/mesh.csv");

    for(int ii=0; ii<mesh->Np; ii++)
    {
        printf("%i, %f, %f\n", ii, mesh->p[ii][0], mesh->p[ii][1]);
    
    }

    for(int ii=0; ii<mesh->Nelem; ii++)
    {
        printf("%i, %i, %i, %i\n", ii, mesh->elem[ii][0], mesh->elem[ii][1], mesh->elem[ii][2]);
    
    }


    return 0;

}

