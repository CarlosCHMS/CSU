
typedef struct {

    int Nb;
    int Np;
    int Nelem;
    
    double** elem;
    double** p;    

} MESH;

MESH* meshInit(char* fileName);

