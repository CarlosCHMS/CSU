

int fluxChoice(char* s);

void entropyFix(SOLVER* solver, double *l);

void fluxRoe(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f);

void fluxAUSMD(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f);

void flux(SOLVER* solver, double rL, double uL, double vL, double pL,
                           double rR, double uR, double vR, double pR, double* f);
                           
void fluxFree(SOLVER* solver, double rL, double uL, double vL, double pL, double* f);
