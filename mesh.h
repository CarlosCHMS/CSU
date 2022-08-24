
typedef struct ELEM
{

    int ii;
    int Np;
    
    int neiN;
    
    int* p;

    struct ELEM** neiL;

    int* f;
    
    double* P; 

    double d;   
    
} ELEMENT;

typedef struct
{

    int Nelem;
    
    char name[50];

    int flagBC;

    ELEMENT** elemL;

} MESHBC;


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
 
MESHBC* meshBCread(FILE* ff, int Nvar);

MESH* meshInit(char* fileName, int Nvar);

void meshPrintBC(MESHBC* bc);

void meshPrint(MESH* mesh);

void meshElementFree(ELEMENT* e);

void meshBCFree(MESHBC* bc);

void meshFree(MESH* mesh);

void elementCenter(ELEMENT* E, MESH* mesh, double* x, double* y);

double meshCalcOmegaTri(MESH* mesh, int p0, int p1, int p2);

double meshCalcDSlateral(MESH* mesh, int ii);

double meshCalcOmega(MESH* mesh, int ii);

double elementIsConnected(ELEMENT* e0, ELEMENT* e1, int* p0, int* p1);

void meshCalcConnection(MESH* mesh);

void meshCalcDS(MESH* mesh, int p0, int p1, double* dSx, double* dSy);

int meshBCIsConnect(ELEMENT* BCe, ELEMENT* e);

void meshCalcNeighbors(MESH* mesh);

void meshBCneighbors(MESHBC* bc, MESH* mesh);

double meshEdgeLength(MESH* mesh, int p0, int p1);

double meshMinEdge(MESH* mesh);

void meshPrintDStotal(MESH* mesh);

void meshCheckUse(MESH* mesh);

int meshPOri(MESH* mesh, ELEMENT* e, int p0, int p1);

void meshCheckBorderOrientation(MESH* mesh);

void meshFixBorderOrientation(MESH* mesh);

void meshCalcFaces(MESH* mesh);
