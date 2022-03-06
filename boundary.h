

void boundaryInlet(SOLVER* solver, double* Pa, double* Pd, double* Pb, double nx, double ny);

void boundaryOutlet(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);

void boundaryWall(SOLVER* solver, double* Pd, double* Pb, double nx, double ny);

void boundaryCalc(SOLVER* solver, MESHBC* bc);

void boundary(SOLVER* solver);
