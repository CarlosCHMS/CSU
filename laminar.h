
void laminarInter(SOLVER* solver);

void boundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f);

void laminarBoundaryVisc(SOLVER* solver, MESHBC* bc);

void laminarBoundary(SOLVER* solver);

void laminarWriteSurf(SOLVER* solver);


