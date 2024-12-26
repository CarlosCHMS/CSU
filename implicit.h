

void implicitCalcD(SOLVER* solver);

void implicitAuxCalcFlux(SOLVER* solver, double rho, double u, double v, double p, double nx, double ny, double* F);

void implicitCalcDeltaFlux(SOLVER* solver, double rho, double u, double v, double p, double d0, double d1, double d2, double d3, double nx, double ny, double* dF);

void implicitFunc(SOLVER* solver, int e0, int e1, int p0, int p1, double** dW);

void implicitFunc_sa(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW);

void implicitLUSGS_L(SOLVER* solver);

void implicitLUSGS_U(SOLVER* solver);

void implicitAuxCalcFlux_sa(SOLVER* solver, double rho, double u, double v, double p, double n, double nx, double ny, double* F);

void implicitCalcDeltaFlux_sa(SOLVER* solver, double rho, double u, double v, double p, double n, double d0, double d1, double d2, double d3, double d4, double nx, double ny, double* dF);

void implicitTest(SOLVER* solver);
