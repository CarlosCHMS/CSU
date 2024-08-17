

void boundaryInlet(SOLVER* solver, double* Pa, double* Pd, double* Pb, double nx, double ny);

void boundaryOutlet(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);

void boundaryWall(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);

void boundaryCalc(SOLVER* solver, MESHBC* bc);

void boundaryCalcVisc(SOLVER* solver, MESHBC* bc);

void boundary(SOLVER* solver);

void boundaryVisc(SOLVER* solver);

void boundaryGetBC(MESH* mesh, INPUT* input);

int boundaryChoice(char* s);

void boundaryCalcPrimitive(SOLVER* solver, MESHBC* bc);

void boundaryCalcFrictionWall(SOLVER* solver, ELEMENT* E, double* fx, double* fy);

void boundaryCalcTensorWall(SOLVER* solver, ELEMENT* E, double* Txx, double* Txy, double* Tyy, double* x, double* yp);

void boundaryOutlet_sa(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);
