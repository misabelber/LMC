#include "/afs/ciemat.es/user/b/bernardos/GitHub/Math/matrixes.h"
#include "DMLimitsLib.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include "TStopwatch.h"
#include "TNtuple.h"
#include "TFile.h"
using namespace std;


void RunBands(Number dm_mass)
{
  TString particle = "Tau";
  //Build wisely the name of the file
  TString masstr;
  ostringstream os;
  if (dm_mass < 1.0){
    os<<dm_mass*1000;
    masstr = os.str()+"GeV";
  }
  else{
    os<<dm_mass;
    masstr = os.str()+"TeV";
  }      
  
  TString filename = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/Bands_"+particle+masstr+"_modelHM.dat";
  const int Ndif = 6;
  const int Nps = 9;
  const int Nbar = Ndif+Nps;
  
  Init(20,20,20,Ndif,Nps);
  ofstream outfilelim;
  outfilelim.open(filename,ios::app);
  
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
                        "J0525-696"};
                        
  TString suf = "_functions_KSP_100GeV-100TeV";
  TString suf_DM = "_jfactorNFW_KSP_100GeV-100TeV";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(dm_mass,particle,suf_DM);

  
  //  VM data_model = Obs_data;
  
  int reals=300;
  
  for (int i=0; i<reals; i++){

    TString realization; 
    int digit1 = 0; int digit2 = 0; int digit3 = 0;
    if (i<10) {digit3=i;}
    if (i>=10 && i<100){digit2=i/10; digit3=i-digit2*10;}
    if (i>=100) {digit1 = i/100; digit2 = (i-digit1*100)/10; digit3 = i-digit1*100-digit2*10;}   
    realization.Form("%d%d%d",digit1,digit2,digit3);
    
    cout << "Realization number " << realization << endl;

    FillContainer_Obs("Irf+CR+DiffuseSources+PS",false,suf+realization);
    // FillContainer_Obs("Irf+CR+DiffuseSources+PS",false,suf);
    //Obs_data = data_model;
    //Ntotal = DataSim(Obs_data);
    
    //Number steps[Nbar+1]={10,0.005,1,1,1,1,0.1,0.001,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
    //Number steps[Nbar+1]={100,0.005,1,1,1,1,0.5,0.001,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
    //Number steps[Nbar+1]={100,0.005,1,1,1,1,0.5,0.001,0.001,0.001,0.001,0.001};
    Number steps[Nbar+1]={100,
			  0.001,
			  0.5,0.5,
			  0.25,0.5,0.25,
			  0.005,0.1,0.01,0.001,0.001,0.001,0.0005,0.001,0.001};
    V Kpars; init(Kpars,Nbar+1);
    
    cout << "Maximizing Likelihood..." << endl;
    cout << endl;
    cout << calc_MaxlogL(Kpars,steps,false) << endl;
    
    cout << "Maximum Likelihood parameters: " << endl;
    for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
    cout << endl;
    cout << endl;
    cout << "Calculating Upper Limit on DM normalization..." << endl;
    V nuis;
    nuis.push_back(1); nuis.push_back(2); nuis.push_back(3);
    Number UpperLimit = Upper_Finder(Kpars,0,nuis);
    cout << "Upper Limit for realization "<< i <<": "  << UpperLimit << endl;
    cout << endl;
    outfilelim << UpperLimit << endl;
  }
  
}


void Run()
{
  const int Ndif = 6;
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
                        "J0525-696"};
                        
  TString suf = "_KSP_100GeV-100TeV";
  TString suf_DM = "_jfactorNFW_KSP_100GeV-100TeV";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(0.1,"W",suf_DM);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  //Ntotal = DataSim(Obs_data);
  
  Number steps[Nbar+1]={10,0.001,1,1,1,1,0.5,0.001,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  
  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  cout << calc_MaxlogL(Kpars,steps,false) << endl;
  
  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  Number UpperLimit = Upper_Minimizer(Kpars,false);
  cout << "Upper Limit: " << UpperLimit << endl;
  cout << endl;
  Number intervals[Nbar+1] = {UpperLimit,0.00025,0.25,0.5,0.02,0.06,0.015,0.015,0.03,0.1,0.04,0.04,0.04,0.3,0.1,0.1,8};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  //  calc_CorrFactors(Kpars,intervals,Cfactors);
}

void RunTest(Number dmmass,Number range)
{
  const int Ndif = 1;
  const int Nps = 0;
  const int Nbar = Ndif+Nps;
  
  Init(20,20,20,Ndif,Nps);
  
  TString extended[Ndif] = {"Irf"};
  
  TString point[Nps] = {};
 
  //TString suf = "_KSPpointing_v2_";
  TString suf = "_rebin_0.1x100_Pointin5deg+1";
  TString suf_DM = "_jfactorNFW_rebin_0.1x100_Pointin5deg+1";
  //TString suf_DM = "_jfactorNFW";
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(dmmass,"W",suf_DM);
  FillContainer_Obs("Irf",true,suf);
  //Ntotal = DataSim(Obs_data);
  
  Number steps[Nbar+1]={10,0.001};//,0.001,0.001};
  
  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  cout << calc_MaxlogL(Kpars,steps,false) << endl;
  
  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  //Number UpperLimit = Upper_Minimizer(Kpars,false);
  //cout << "Upper Limit: " << UpperLimit << endl;
  //cout << endl;
  Number intervals[Nbar+1] = {range,0.00025};//,0.25,0.5};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  //calc_CorrFactors(Kpars,intervals,Cfactors);
}

void RunRebinned()
{
  const int Ndif = 6;
  const int Nps = 10;
  const int Nbar = Ndif+Nps;
  
  Init(20,100,100,Ndif,Nps);
  
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
  TString suf_DM = "rebin_0.1x100";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(1,"W",suf_DM);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  Ntotal = DataSim(Obs_data);
  
  Number steps[Nbar+1]={10,0.001,1,1,1,1,0.5,0.001,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  
  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  cout << calc_MaxlogL(Kpars,steps,false) << endl;
  
  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  Number UpperLimit = Upper_Minimizer(Kpars,false);
  cout << "Upper Limit: " << UpperLimit << endl;
  cout << endl;
  Number intervals[Nbar+1] = {UpperLimit,0.00025,0.25,0.5,0.02,0.06,0.015,0.1,0.15,4.5,1,4,2,1,0.1,0.1,200};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  //calc_CorrFactors(Kpars,intervals,Cfactors);
}

void CheckJFACTOR()
{
  const int Ndif = 6;
  const int Nps = 1;
  const int Nbar = Ndif+Nps;
  
  Init(20,100,100,Ndif,Nps);
    
  TString extended[Ndif] = {"Irf","Leptonic","Hadronic","3FHL_J0500.9-6945e","3FHL_J0530.0-6900e","3FHL_J0531.8-6639e"};
  TString point[Nps] = {"J0537-691"};
  TString suf = "_KSPpointing_v2_";
  TString suf_DM = "_jfactorgamma1.5";
  
  TString DMtypes[3] = {"_jfactorNFW","_jfactorgamma0.5","_jfactorgamma1.5"};
  Number DMmasses[6] = {0.100,0.500,1,10,50,100};
  TString DMparticles[4] = {"W","b","Tau","Mu"};
  
  Number steps[Nbar+1]={10,0.001,1,0.01,0.001,1,0.001};
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_Obs("Irf+CR+DiffuseSources+1PS",true,suf);
  
  for (int npart=3; npart<4; npart++)
    {
      for (int ntype=0; ntype<1; ntype++)
	{
	  TString filenamelim = "/home/queenmab/GitHub/LMC/results/Limits"+DMparticles[npart]+DMtypes[ntype];
	  TString filenamecfact = "/home/queenmab/GitHub/LMC/results/Cfactors"+DMparticles[npart]+DMtypes[ntype];
	  ofstream outfilelim;
	  ofstream outfilecfact;
	  outfilelim.open(filenamelim,ios::app);
	  outfilecfact.open(filenamecfact,ios::app);
	  
	  for (int nmass=5; nmass<6; nmass++)
	    {
	      FillContainer_DM(DMmasses[nmass],DMparticles[npart],DMtypes[ntype]);
	      V Kpars; init(Kpars,Nbar+1);
	      
	      cout << "Maximizing Likelihood..." << endl;
	      cout << endl;
	      cout << calc_MaxlogL(Kpars,steps,false) << endl;
	      
	      cout << "Maximum Likelihood parameters: " << endl;
	      for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
	      cout << endl;
	      cout << endl;
	      cout << "Calculating Upper Limit on DM normalization..." << endl;
	      Number UpperLimit = Upper_Minimizer(Kpars,false);
	      cout << "Upper Limit: " << UpperLimit << endl;
	      cout << endl;
	      outfilelim << DMmasses[nmass] << "  " << UpperLimit << endl;
	      Number intervals[Nbar+1] = {UpperLimit,0.00025,0.25,0.5,0.02,0.06,0.015,0.1};
	      V Cfactors;
	      cout << "Calculating Correlation Factors..." << endl;
	      //calc_CorrFactors(Kpars,intervals,Cfactors);
	      outfilecfact << DMmasses[nmass] << "  ";
	      for (int ii=0; ii<Nbar; ii++) outfilecfact << Cfactors[ii] << "  ";
	      outfilecfact << endl;
	    }
	}
    }
}

void Pruebas(Number dmmass){
  TStopwatch t;
  t.Start();
  const int Ndif = 6;
  const int Nps = 9;
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
                        "J0525-696"};
  
  TString suf = "_functions_KSP_100GeV-100TeV";
  TString suf_DM = "_jfactorNFW_KSP_100GeV-100TeV";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(dmmass,"W",suf_DM);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  //Ntotal = DataSim(Obs_data);
  
  Number steps[Nbar+1]={200,0.0005,0.25,0.25,0.25,0.25,0.25,0.001,0.05,0.005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005};
  
  V Kpars; init(Kpars,Nbar+1);
  
  cout << calc_MaxlogL(Kpars,steps,false) << endl;
  for (int ii=0; ii<Nbar+1; ii++){
    cout << Kpars[ii] << "  ";
  }
  cout << endl;
  V Cfactors;
  V nuis;
  nuis.push_back(1); nuis.push_back(2); nuis.push_back(3);
  Number Upperlimit = Upper_Finder(Kpars,0,nuis);
  cout << "Upper limit: " << Upperlimit << endl;
  t.Stop();
  t.Print();
  //Number intervals[Nbar+1] = {Upperlimit,0.00025,0.25,0.5,0.02,0.06,0.015,0.015,0.03,0.1,0.04,0.04,0.04,0.3,0.1,0.1,8};
  Number intervals[Nbar+1] = {fabs(Upperlimit-Kpars[0]),0.0003,0.3,0.5,0.1,0.06,0.015,0.015,0.03,0.2,0.04,0.04,0.04,0.3,0.1,0.1};
  
  TFile *PSfile = new TFile("parameterspace.root","RECREATE");
  TNtuple *ParSpace;
  calc_CorrFactors(Kpars,intervals,Cfactors,ParSpace);
  ParSpace->Write();
  
}

void Bin(Number dmmass,int first,int more){
  TStopwatch t;
  t.Start();
  const int Ndif = 6;
  const int Nps = 10;
  const int Nbar = Ndif+Nps;
  firstebin = first;
  nebins = firstebin+more;

  Init(nebins,20,20,Ndif,Nps);
  

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
  TString suf = "_rebin_0.1x100_Pointin5deg";
  TString suf_DM = "_jfactorNFW_rebin_0.1x100_Pointin5deg";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(dmmass,"W",suf_DM);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  Ntotal = DataSim(Obs_data);
  Number steps[Nbar+1]={0.01,0.005,1,1,1,1,0.1,0.001,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  
  V Kpars; init(Kpars,Nbar+1);
  

  cout << calc_MaxlogL(Kpars,steps,false) << endl;
  for (int ii=0; ii<Nbar+1; ii++){
    cout << Kpars[ii] << "  ";
  }
  cout << endl;
  cout << logL(Kpars,firstebin,nebins) << endl;
  
  V Cfactors;
  
  //  Number Upperlimit = Upper_Minimizer(Kpars);
  //cout << Upper_Function(Kpars,-0.0001) << endl;
  ofstream outfile;
  outfile.open("fluxlimits.dat",ios::app);
  V nuis;
  nuis.push_back(1); nuis.push_back(2); nuis.push_back(3);
  Number Upperlimit = Upper_Finder(Kpars,0,nuis);
  cout << "Upper limit: " << Upperlimit << endl;
  outfile << Upperlimit << endl;
  t.Stop();
  t.Print();


  
}
