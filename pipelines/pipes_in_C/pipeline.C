#include "/home/queenmab/GitHub/Math/matrixes.h"
#include "DMLimitsLib.h"

#include <iostream>

using namespace std;

void Run()
{

  Init(20,20,20,6,3);
  
  TString extended[N_ext] = {"Irf","Leptonic","Hadronic","3FHL_J0500.9-6945e","3FHL_J0530.0-6900e","3FHL_J0531.8-6639e"};
  TString point[N_ps] = {"J0537-691","J0524.5-6937","J0534.1-6732"};
  TString suf = "_KSPpointing_v2_";
  TString suf_DM = "_003-100_KSP_v2";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(0.500,"W",suf_DM);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  //Ntotal = DataSim(Obs_data);

  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  cout << calc_MaxlogL(Kpars,false) << endl;

  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  Number UpperLimit = Upper_Minimizer(Kpars);
  cout << "Upper Limit: " << UpperLimit << endl;
  cout << endl;
  Number intervals[Nbar+1] = {UpperLimit,0.00025,0.25,0.5,0.02,0.06,0.015,0.1,0.15,4.5};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  calc_CorrFactors(Kpars,intervals,Cfactors);
}

void Run2()
{
  Init(20,20,20,1,0);
  
  TString extended[N_ext] = {"Irf"};
  TString point[N_ps];
  TString suf = "_KSPpointing_v2_";
  TString suf_DM = "_003-100_KSP_v2";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(1,"W",suf_DM);
  FillContainer_Obs("Irf",true,suf);
  //Ntotal = DataSim(Obs_data);

  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  calc_MaxlogL(Kpars,false);

  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  Number UpperLimit = Upper_Minimizer(Kpars);
  cout << "Upper Limit: " << UpperLimit << endl;
  cout << endl;
  Number intervals[Nbar+1] = {UpperLimit,0.00025};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  calc_CorrFactors(Kpars,intervals,Cfactors);
}
