
void saInitU(SOLVER* solver, CONDITION* inside);

void saGrad(SOLVER* solver);

void saInterFace(SOLVER* solver);

void saInterSource(SOLVER* solver);

void saInter(SOLVER* solver);

void saCalcFace(double ni, double ni_L, double r, double dnix, double dniy, double* fv1, double* tx, double* ty);

void saCalcSource(double ni, double ni_L, double S, double d, double rho, double drx, double dry, double dnix, double dniy, double* Qt);

void saBoundaryFace(SOLVER* solver, MESHBC* bc);

void saBoundary(SOLVER* solver);

void saCalcD(MESH* mesh);

void saCalcTensorWall(SOLVER* solver, ELEMENT* E, double* Txx, double* Txy, double* Tyy, double* x, double* yp);
