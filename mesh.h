
typedef struct ELEM
{

    int ii;
    int Np;
    int neiN;
    
    int* p;
    int* f;
    
    double* P; 
    double d;       
    double omega;

    struct ELEM** neiL;    
    
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

typedef struct CON
{

    int data[4];
    struct CON* next; 

} CONNECTION;

typedef struct FACE
{
    int full;
    int p0;
    int p1;
    int e0;
    int e1;
    struct FACE* next; 

} FACETYPE;


char meshGetWord(FILE* ff, char* s);
 
MESHBC* meshBCread(FILE* ff, int Nvar);

MESH* meshInit(char* fileName, int Nvar, int axi);

void meshPrintBC(MESHBC* bc);

void meshPrint(MESH* mesh);

void meshElementFree(ELEMENT* e);

void meshBCFree(MESHBC* bc);

void meshFree(MESH* mesh);

void elementCenter(ELEMENT* E, MESH* mesh, double* x, double* y);

double meshCalcOmegaTri(MESH* mesh, int p0, int p1, int p2);

double meshCalcDSlateral(MESH* mesh, int ii);

double meshCalcOmega(MESH* mesh, int ii);

void meshUpdateOmega(MESH* mesh);

double elementIsConnected(ELEMENT* e0, ELEMENT* e1, int* p0, int* p1);

int meshSameFace(FACETYPE* f0, FACETYPE* f1);

void meshCalcConnection1(MESH* mesh);

void meshCalcConnection2(MESH* mesh);

void meshPrintConnection(MESH* mesh, int N);

void meshCalcDS(MESH* mesh, int p0, int p1, double* dSx, double* dSy);

void meshCalcDS2(MESH* mesh, int p0, int p1, double* nx, double* ny, double* dS);

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
