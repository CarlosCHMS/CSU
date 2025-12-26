
SST_TRANS* sstTransInit();

void sstTransFree(SST_TRANS* sstTrans);

double sstTransAlphas1(double Rt);

double sstTransF1(SST_TRANS* sst, SSTVAR* var, double n_L_term, double sqrtk_term);

double sstTransF2(SST_TRANS* sst, SSTVAR* var, double n_L_term, double sqrtk_term);

void sstTransFlux(SST_TRANS* trans, SSTVAR* var);

double sstTransBs1(double Rt);

double sstTransAlpha1(double Rt);

void sstTransSources(SST_TRANS* trans, SSTVAR* var);
