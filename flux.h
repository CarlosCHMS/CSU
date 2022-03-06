

int fluxChoice(char* s);

void entropyFix(SOLVER* solver, double *l);

void fluxRoe(SOLVER* solver, double rL, double uL, double vL, double pL,
                             double rR, double uR, double vR, double pR, double* f);

void fluxAUSMD(SOLVER* solver, double rL, double uL, double vL, double pL,
                               double rR, double uR, double vR, double pR, double* f);

void flux(SOLVER* solver, double rL, double uL, double vL, double pL,
                          double rR, double uR, double vR, double pR, double* f);
                           
void fluxFree(SOLVER* solver, double rL, double uL, double vL, double pL, double* f);
