#pragma once
#include "/afs/ciemat.es/user/b/bernardos/GitHub/Math/matrixes.h"
#include <vector>
#include "TString.h"
#include "TNtuple.h"
//Important Paths where data and models are stored

TString dir = "/scratch/bernardos/LMC/Obs_";

TString dir_PS = "/scratch/bernardos/LMC/Obs_PS"; //Path for pointsources models
TString dir_DM = "/scratch/bernardos/LMC/Obs_DM"; //Path for DM models 


//Global variables

//Number of bins
int firstebin=0;
int nebins=0; //Energy bins
int nxbins=0; //Spatial bins
int nybins=0;

//Number of Baryonic Background components
int Nbar=0;
int N_ext=0;
int N_ps=0;

//Data and Model Containers

VM DM_model;
vector<VM> Bkg_model;

VM Obs_data;

// Number of counts
Number N_dm = 0; //In the DM model
V N_bkg; //In the baryonic background models

Number Ntotal = 0; //In the data


void Init();
void Init(int ne,int nx,int ny,int Nex,int Nps);

void ReadFits();
Number ReadFits(VM &data, TString filename);

void FillContainer_Bkg();
void FillContainer_Bkg(TString ext[],TString ps[],TString suf);

void FillContainer_DM();
void FillContainer_DM(Number dm_mass, TString particle, TString suf);

void FillContainer_Obs();
void FillContainer_Obs(TString obsname, bool modcube, TString suf);

void DataSim();
Number DataSim(VM &data);

void logL();
Number logL(V Kpars,int firstebin,int nebins);

void calc_P(V &P, vector<VM> &p_compbin);
void Conjugate_Gradients();
Number Conjugate_Gradients(V &Kpars);

Number EM_Estimate(VM &data, VM &dm, vector<VM> &back,Number &Nest_dm, V &Nest_bkg, Number DMnorm, V P);
Number EM_Update_pars(Number Nest_dm, V Nest_bkg,V &Kpars, Number &DMnorm, V &P);
void Expectation_Maximization();
Number Expectation_Maximization(V &Kpars);

void My_Minimizer();
Number My_Minimizer(V &Kpars, Number steps[], Number tol=HUGE_VAL);

void calc_MaxlogL();
Number calc_MaxlogL(V &Kpars, Number steps[], bool BestCase=false,Number tol=HUGE_VAL);

void Upper_Minimizer();
Number Upper_Minimizer(V &Kpars,bool test=false, Number tol=HUGE_VAL);

Number Upper_Function(V Kpars,const int which_goal,const int which_nuis,Number nuis_step);
Number Upper_Finder(V Kpars,const int which_goal, V nuis);

void calc_CorrFactors();
void calc_CorrFactors(V Kpars, Number intervals[], V &Cfactors, TNtuple* &ParSpace);



