
GASPROP* gaspropInit(double gamma, double R, int TP);

void gaspropFree(GASPROP* gas);

double gasprop_T2e(GASPROP* gas, double T);

double gasprop_e2T(GASPROP* gas, double e);

double gasprop_T2c(GASPROP* gas, double T);

double gasprop_T2Cv(GASPROP* gas, double T);

double gasprop_T2Cp(GASPROP* gas, double T);

double gasprop_T2gamma(GASPROP* gas, double T);

void gasprop_T2eCv(GASPROP* gas, double T, double* e, double* Cv);

double gasprop_e2Taprox(GASPROP* gas, double e0, double T0, double e1);

double gasprop_T2fentro(GASPROP* gas, double T);

double gasprop_Tp2entropy(GASPROP* gas, double T, double p);
