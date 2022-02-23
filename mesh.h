
typedef struct
{

    int Nelem;
    int** elem;
    char name[50];
    
    int* domain;

} MESHBC;

typedef struct
{

    int Ndim;
    int Nelem;
    int Np;
    int Nmark;
    int Ncon;
    
    int** elem;
    
    double** p;
    
    MESHBC** bc;

    int** con;

} MESH;

char meshGetWord(FILE* ff, char* s);
 
MESHBC* meshBCread(FILE* ff);

MESH* meshInit(char* fileName);

void meshPrintBC(MESHBC* bc);

void meshPrint(MESH* mesh);

void meshBCFree(MESHBC* bc);

void meshFree(MESH* mesh);

void meshElemCenter(MESH* mesh, int ii, double* x, double* y);

double meshCalcOmega(MESH* mesh, int ii);

double meshIsConnected(MESH* mesh, int ii, int jj, int* p0, int* p1);

double meshCalcConnection(MESH* mesh);

void meshCalcDS(MESH* mesh, int p0, int p1, double* dSx, double* dSy);

int meshBCIsConnect(int* BCp, int* p);

void meshBCDomain(MESHBC* bc, MESH* mesh);

double meshEdgeLength(MESH* mesh, int p0, int p1);

double meshMinEdge(MESH* mesh);
