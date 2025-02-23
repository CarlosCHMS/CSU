
void saCC_InterFace(SOLVER* solver);

void saCC_InterSource(SOLVER* solver);

void saCC_Inter(SOLVER* solver);

void saCC_CalcFace(double ni, double ni_L, double r, double dnix, double dniy, double drx, double dry, double* fv1, double* tx, double* ty);

void saCC_CalcSource(double ni, double ni_L, double S, double d, double rho, double drx, double dry, double dnix, double dniy, double* Qt);

void saCC_BoundaryFace(SOLVER* solver, MESHBC* bc);

void saCC_BoundaryFaceViscFlux(SOLVER* solver, MESHBC* bc, int ii, double* f, double* miEddy);

void saCC_Boundary(SOLVER* solver);

void saCC_SolverWriteSurf(SOLVER* solver);
