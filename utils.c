#include<stdlib.h>

double** tableMallocDouble(int Nrow, int Ncol)
{

    double** M = malloc(Nrow*sizeof(double*));

    for(int ii=0; ii<Nrow; ii++)
    {    
        M[ii] = malloc(Ncol*sizeof(double));
    }

    return M;

}    
 
int** tableMallocInt(int Nrow, int Ncol)
{

    int** M = malloc(Nrow*sizeof(int*));

    for(int ii=0; ii<Nrow; ii++)
    {    
        M[ii] = malloc(Ncol*sizeof(int));
    }

    return M;

} 

void tableFreeDouble(double** M, int Nrow)
{

    for(int ii=0; ii<Nrow; ii++)
    {    
        free(M[ii]);
    }    

    free(M);

}

void tableFreeInit(int** M, int Nrow)
{

    for(int ii=0; ii<Nrow; ii++)
    {    
        free(M[ii]);
    }    

    free(M);

}

double duration(struct timeval start, struct timeval stop){

    double tstart, tstop;
    tstart = (double)start.tv_sec + start.tv_usec/1000000.;
    tstop = (double)stop.tv_sec + stop.tv_usec/1000000.;

    return tstop - tstart;

}

