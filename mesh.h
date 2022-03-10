
typedef struct
{

    int Nelem;
    int** elem;
    char name[50];
    
    int* domain;

    int flagBC;


} MESHBC;

typedef struct
{

    int Np;
    
    int neiN;
    
    int* p;
    int* nei;
    int* f;
    
} ELEMENT;

typedef struct
{

    int Ndim;
    int Nelem;
    int Np;
    int Nmark;
    int Ncon;
    int axi;
    int order;
    
    int** con;    

    double** p;

    ELEMENT** elemL;
    
    MESHBC** bc;

} MESH;

char meshGetWord(FILE* ff, char* s);
 
MESHBC* meshBCread(FILE* ff);

MESH* meshInit(char* fileName);

void meshPrintBC(MESHBC* bc);

void meshPrint(MESH* mesh);

void meshBCFree(MESHBC* bc);

void meshFree(MESH* mesh);

void meshElemCenter(MESH* mesh, int ii, double* x, double* y);

double meshCalcOmegaTri(MESH* mesh, int p0, int p1, int p2);

double meshCalcDSlateral(MESH* mesh, int ii);

double meshCalcOmega(MESH* mesh, int ii);

double elementIsConnected(ELEMENT* e0, ELEMENT* e1, int* p0, int* p1);

void meshCalcConnection(MESH* mesh);

void meshCalcDS(MESH* mesh, int p0, int p1, double* dSx, double* dSy);

int meshBCIsConnect(int* BCp, ELEMENT* e);

void meshBCDomain(MESHBC* bc, MESH* mesh);

double meshEdgeLength(MESH* mesh, int p0, int p1);

double meshMinEdge(MESH* mesh);

void meshPrintDStotal(MESH* mesh);

void meshCheckUse(MESH* mesh);

int meshPOri(MESH* mesh, ELEMENT* e, int p0, int p1);

void meshCheckBorderOrientation(MESHBC* bc, MESH* mesh);

void meshCalcNeighbors(MESH* mesh);
