#include "/afs/ciemat.es/user/b/bernardos/GitHub/Math/matrixes.h"
#include "SrcAnalysisLib.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void RunT(){

  const  int Ndif = 6;
  const int Nps = 10;
  const int Nbar = Ndif+Nps;
  
  Init(20,20,20,Ndif,Nps);
  
  TString extended[Ndif] = {"Irf",
			    "Leptonic",
                            "Hadronic",
                            "3FHL_J0500.9-6945e",
                            "3FHL_J0530.0-6900e",
                            "3FHL_J0531.8-6639e"};
  
  TString point[Nps] = {"J0537-691",
			"J0524.5-6937",
                        "J0534.1-6732",
                        "J0525.2-6614",
                        "J0535.3-6559",
                        "J0454.6-6825",
                        "J0537.0-7113",
                        "J0535-691",
			"J0525-696",
			"J0509.9-6418"};
  TString suf = "_rebin_0.1x100";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  //Ntotal = DataSim(Obs_data);
    
  V Kpars; init(Kpars,Nbar);
  Number tol=HUGE_VAL;
  Number MaxlogL = 0;                                                                              
  MaxlogL = Conjugate_Gradients(Kpars);
  cout << MaxlogL << endl;
  for (int ii=0; ii<Nbar; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  MaxlogL=Expectation_Maximization(Kpars);
  cout << MaxlogL << endl;                                                                 
  for (int ii=0; ii<Nbar; ii++) cout << Kpars[ii] << "  ";   
  cout << endl;
  //Number steps[Nbar]={0.001,1,1,1,1,0.5,0.01,0.5,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  //MaxlogL = My_Minimizer(Kpars,steps,tol);
  for (int ii=0; ii<Nbar; ii++) cout << Kpars[ii] << "  ";                        
  cout << endl;
  
  Number intervals[Nbar] = {0.00025,
			    0.25,
			    0.5,
			    0.02,
			    0.06,
			    0.015,
			    0.015,
			    0.03,
			    0.1,
			    0.04,
			    0.04,
			    0.04,
			    0.3,
			    0.1,
			    0.1,
			    8};
    
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  calc_CorrFactors(Kpars,intervals,Cfactors);
}

void TS(){

  const  int Ndif = 6;
  const int Nps = 10;
  const int Nbar = Ndif+Nps;
  
  Init(20,20,20,Ndif,Nps);
  
  TString extended[Ndif] = {"Irf",
			    "Leptonic",
                            "Hadronic",
                            "3FHL_J0500.9-6945e",
                            "3FHL_J0530.0-6900e",
                            "3FHL_J0531.8-6639e"};
  
  TString point[Nps] = {"J0537-691",
                        "J0524.5-6937",
                        "J0534.1-6732",
                        "J0525.2-6614",
                        "J0535.3-6559",
                        "J0454.6-6825",
                        "J0537.0-7113",
                        "J0535-691",
			"J0525-696",
			"J0509.9-6418"};
  TString suf = "_rebin_0.1x100";
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);  
  FillContainer_Bkg(extended,point,suf);

  V Kpars; init(Kpars,Nbar);
  Number tol=HUGE_VAL;
  Number MaxlogL = 0;                                                                              
  MaxlogL = Conjugate_Gradients(Kpars);
  cout <<"MaxlogL: " << "  " <<  MaxlogL << endl;
  
  vector<TString> extended_;
  vector<TString>  point_;
  
  Init(20,20,20,Ndif-1,Nps);
    
  for(int i=1; i<Ndif; i++){
    extended_.clear();
    for (int ii=0; ii<Ndif; ii++) extended_.push_back(extended[ii]);
    extended_.erase(extended_.begin()+i);
    TString ext[Ndif];
    for (int ii=0; ii<Ndif; ii++) {ext[ii] = extended_[ii];}
    FillContainer_Bkg(ext,point,suf);
    
    //Ntotal = DataSim(Obs_data);
    
    V Kpars; init(Kpars,Nbar);
    Number tol=HUGE_VAL;
    Number NullogL = Conjugate_Gradients(Kpars);
    cout <<"Source: " << extended[i] << " TS: " <<  2*(MaxlogL-NullogL) << endl;
    //for (int ii=0; ii<Nbar; ii++) cout << Kpars[ii] << "  ";
    //cout << endl;
  }

  Init(20,20,20,Ndif,Nps-1);
    
  for(int i=0; i<Nps; i++){
    point_.clear();
    for (int ii=0; ii<Nps; ii++) point_.push_back(point[ii]);
    point_.erase(point_.begin()+i);
    TString pnt[Nps];
    for (int ii=0; ii<Nps; ii++) {pnt[ii] = point_[ii];}
    FillContainer_Bkg(extended,pnt,suf);
    
    //Ntotal = DataSim(Obs_data);
    
    V Kpars; init(Kpars,Nbar);
    Number tol=HUGE_VAL;
    Number NullogL = Conjugate_Gradients(Kpars);
    cout << "Source: " << point[i] << "  TS: " << 2*(MaxlogL-NullogL) << endl;
    //for (int ii=0; ii<Nbar; ii++) cout << Kpars[ii] << "  ";
    //cout << endl;
  }  


}
