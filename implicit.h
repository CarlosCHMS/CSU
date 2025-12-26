

void implicitCalcD(SOLVER* solver);

void implicitCalcDeltaFlux(SOLVER* solver, double rho, double u, double v, double p, double d0, double d1, double d2, double d3, double nx, double ny, double* dF);

void implicitFunc(SOLVER* solver, int e0, int e1, int p0, int p1, double** dW);

void implicitFunc_sa(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW);

void implicitLUSGS_L(SOLVER* solver);

void implicitLUSGS_U(SOLVER* solver);

void implicitCalcDeltaFlux_sa(SOLVER* solver, double rho, double u, double v, double p, double n, double d0, double d1, double d2, double d3, double d4, double nx, double ny, double* dF);

void implicitTest(SOLVER* solver);

void implicitInitDPLUR(SOLVER* solver);

void implicitFreeDPLUR(SOLVER* solver);

void implicitUpdateA(SOLVER* solver);

void implicitUpdateA_sa(SOLVER* solver);

void implicitMultA(SOLVER* solver, double** x, double** y);

void implicitCalcDPLUR(SOLVER* solver);

void implicitAuxCalcFlux2(SOLVER* solver, double U0, double U1, double U2, double U3, double p, double nx, double ny, double* F);

void implicitAuxCalcFlux_sa2(SOLVER* solver, double U0, double U1, double U2, double U3, double U4, double p, double nx, double ny, double* F);

void implicitAuxCalcFlux_sst2(SOLVER* solver, double U0, double U1, double U2, double U3, double U4, double U5, double p, double nx, double ny, double* F);

void implicitCalcDeltaFlux_sst(SOLVER* solver, double rho, double u, double v, double p, double k, double om, double d0, double d1, double d2, double d3, double d4, double d5, double nx, double ny, double* dF);

void implicitFunc_sst(SOLVER* solver, int e0, int e1, int p0, int p1, int face1, double** dW);
