
LIMITER* limiterInit(int type, double K, SOLVER* solver);

void limiterFree(LIMITER* limiter);

void limiterUpdate(LIMITER* limiter, SOLVER* solver);

void limiterCalc(LIMITER* limiter, SOLVER* solver, int ii, double* Pref2);

double limiterV2(double Ui, double Umin, double Umax, double d2, double ee);
