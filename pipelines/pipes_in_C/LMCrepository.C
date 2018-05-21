#include "/afs/ciemat.es/user/b/bernardos/GitHub/Math/matrixes.h"
#include "/afs/ciemat.es/user/b/bernardos/GitHub/Math/linear_fitter.h"
#include "TMatrixD.h"
#include "TFITS.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMultiGraph.h"
#include "TF2.h"
#include "TRandom.h"
#include "TFile.h"
#include "TFitter.h"
#include <ctime>


#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
#define PI 3.14159265359 


vector<VM> data_bkg;
const int N=7; //Number of components in the Baryonic Background
VM obs_data,dm_data;
double A; //DM Normalization factor
V Pars; //BKG components normalization factors

double N_dm; //Total number of counts of dark matter in the observation
V N_bkg; // Total number of counts of background components in the observation

double model_dm; //Total number of counts of dark matter in the model
V model_bkg; // Total number of counts of background components in the model
double total_data; //Total number of counts in the data

double _loglike; 
VM p_dm;
vector<VM> p_bkg;
vector<VM> p_compbin;

V P;
M Mat;

int ebinning=20;
int nebins, nxbins, nybins,firstebin; 

//Dark Matter final state and masses to try
TString particle = "W";
const int nmasses = 1;
//float masses[nmasses] = {0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1,2,3,4,5,6,7,8,9,10,50,100};
//float masses[nmasses] = {0.100,0.200,0.500,1,5,10,50,100};
//float masses[nmasses] = {50,100};
float masses[nmasses] = {1};


///////FUNCTIONS TO READ DATA AND FILL MATRIXES///////////

void ReadFits(VM &data,TString filename){
  nebins=ebinning;
  data.clear();
  //Open the model component fits file
  TFITSHDU *hdu =  new TFITSHDU(filename);
  if (hdu == 0) {
    printf("ERROR: could not access the HDU\n"); return;
  }
  int counter=0;
  for (int bin=0; bin<nebins; bin++){
    TMatrixD *mat =  hdu->ReadAsMatrix(bin); 
    if (bin==0){nxbins = mat->GetNcols(); nybins = mat->GetNrows();}

    M bin_data;
    for (int i=0; i<nxbins;i++){
      V column;
      for (int j=0; j<nybins; j++){
	column.push_back((*mat)(i,j));
      }
      bin_data.push_back(column);
    }
    data.push_back(bin_data);
    delete mat;
  }
  delete hdu;
}

void FillDataM(Number &dm_mass){
  
  Pars.clear();
  data_bkg.clear();
  obs_data.clear();
  dm_data.clear();
  N_dm=0;
  N_bkg.clear();
  _loglike=0;
  
  p_dm.clear();
  p_bkg.clear();
  p_compbin.clear();
  
  model_dm=0;
  model_bkg.clear();
  
  total_data=0;
  
  
  for (int i=0; i<N; i++){
    VM data;
    data_bkg.push_back(data);
  }
  
  
  // Read Background model components data
  // Components files
  TString dir = "/scratch/bernardos/LMC/Obs_";
  TString dir_PS = "/scratch/bernardos/LMC/Obs_PS";
  //TString components[N] = {"Irf","Leptonic","Hadronic","3FHL_J0531.8-6639e","3FHL_J0530.0-6900e","3FHL_J0500.9-6945e","J0454.6-6825","J0529.8-7242","J0525.2-6614","J0537.0-7113","J0509.9-6418","J0601.5-7237","J0524.5-6937","J0535.3-6559","J0601.2-7035","J0537-691","J0525-696","J0534.1-6732","J0535-691"};
  TString components[N] = {"Irf","Leptonic","Hadronic","3FHL_J0531.8-6639e","3FHL_J0530.0-6900e","3FHL_J0500.9-6945e","J0537-691"};
  //TString components[N] = {"Leptonic","Hadronic"};
  //TString components[N] = {"Leptonic"};
  //TString components[N] = {"Irf"};
  TString filename;
  for (int i =0; i<N; i++){
    //if (i<6) filename = dir+components[i]+"/modcube_LMC_"+components[i]+"_prod3_003-100_300h.fits";
    if (i<6) filename = dir+components[i]+"/modcube_LMC_"+components[i]+"_KSPpointing_v2_.fits";
    else {filename = dir_PS+"/modcube_LMC_"+components[i]+"_KSPpointing_v2_.fits";}
    ReadFits(data_bkg[i],filename);
  }
  
  //Read Observations data
  
  //filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+PS/cntcube_LMC_Irf+CR+DiffuseSources+PS_300h.fits";
  //filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources/cntcube_LMC_Irf+CR+DiffuseSources_003-100_300h.fits";
  //TString filename = "/scratch/bernardos/LMC/Obs_Leptonic+Hadronic/cntcube_LMC_Leptonic+Hadronic_prod3_003-100_300h.fits";
  //TString filename = "/scratch/bernardos/LMC/Obs_Leptonic/cntcube_LMC_Leptonic_prod3_003-100_300h.fits";
  //filename = "/scratch/bernardos/LMC/Obs_Irf/cntcube_LMC_Irf_prod3_003-100_300h.fits";
 
  //filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources/cntcube_LMC_Irf+CR+DiffuseSources_KSPpointing_v2_.fits";

  //filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+1PS/cntcube_LMC_Irf+CR+DiffuseSources+1PS_KSPpointing_v2_.fits";
  
  filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+PS/modcube_LMC_Irf+CR+DiffuseSources+PS_KSPpointing_v2_.fits";

  ReadFits(obs_data,filename);
  
  // Dark Matter Final State and masses
  TString masstr;
  
  //for (int mass=0; mass<nmasses; mass++){
  //Build wisely the name of the file model from DM
  if (dm_mass < 1.0){
    ostringstream os;
    os<<dm_mass*1000;
    masstr = os.str()+"GeV";      
  }
  else{
    ostringstream os;
    os<<dm_mass;
    masstr = os.str()+"TeV"; 
  }
  
  TString modelname = "dm_LMC_"+particle+masstr;
  /*if (dm_mass==1) filename = "/scratch/bernardos/LMC/Obs_DM/modcube_"+modelname+"_003-100_300h_40ebins_002binning.fits";
  else if(dm_mass==0.500) filename = "/scratch/bernardos/LMC/Obs_DM/modcube_"+modelname+"_003-100_300h_40ebins_025binning.fits";
  else*/// filename = "/scratch/bernardos/LMC/Obs_DM/modcube_"+modelname+"_003-100_300h_40ebins_01binning.fits";
  //filename = "/scratch/bernardos/LMC/Obs_DM/modcube_"+modelname+"_003-100_300h.fits";
  if (dm_mass == 0) filename = "/scratch/bernardos/LMC/Obs_DM/modcube_dm_PowerLaw_003-100_KSP.fits";
  //else filename = "/scratch/bernardos/LMC/Obs_DM/modcube_"+modelname+"_003-100_300h.fits";
  else filename = "/scratch/bernardos/LMC/Obs_DM/modcube_"+modelname+"_003-100_KSP_v2.fits";
  
  ReadFits(dm_data,filename);

  /*for (int i=0; i<nebins; i++) {
    for (int j=0; j<nxbins; j++) {
    for (int k=0; k<nybins; k++) {
    dm_data[i][j][k]=10000*dm_data[i][j][k];
    
    }
    }
    }*/
  
  

}

void DataSim(VM &data){
  gRandom->SetSeed(0); 
  //Observations Model Cube
  //TString filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources/modcube_LMC_Irf+CR+DiffuseSources_KSPpointing_v2_.fits";
  
  //TString filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+1PS/modcube_LMC_Irf+CR+DiffuseSources+1PS_KSPpointing_v2_.fits";
  
  TString filename = "/scratch/bernardos/LMC/Obs_Irf+CR+DiffuseSources+PS/modcube_LMC_Irf+CR+DiffuseSources+PS_KSPpointing_v2_.fits";

  ReadFits(data,filename);

  int ebins = data.size();
  int xbins = data[0].size();
  int ybins = data[0][0].size();

  for (int i=0; i<ebins; i++){
    for (int j=0; j<xbins; j++){
      for (int k=0; k<ybins; k++){
	Number mean = data[i][j][k];
	data[i][j][k] = gRandom->Poisson(mean);
      }
    }
  }
}


//////////////////////////////////////////////////////////////

//////////////FUNCTIONS TO CALCULATE PROBABILITIES///////////

void calc_P(){

  P.clear();
  Mat.clear();

  p_dm; p_dm = dm_data;
  p_bkg; p_bkg = data_bkg;
  
  total_data = 0;
  init(model_bkg,N);
  
  for (int ebin=firstebin; ebin<nebins; ebin++){
    for (int i=0; i<nxbins; i++){
      for (int j=0; j<nybins; j++){
	total_data+=obs_data[ebin][i][j];
	if (obs_data[ebin][i][j] < 10) continue;
	model_dm += dm_data[ebin][i][j];
	for (int ii=0; ii<N; ii++) {
	  model_bkg[ii]+=data_bkg[ii][ebin][i][j];	  
	}
      }
    }
  }  
  
  for (int ii=0; ii<N; ii++) {if (model_bkg[ii]<1e-15) model_bkg[ii] = 1e-20;}
  
  V P_n;
  init(N_bkg,N); init(P,N+1); 
  //Calculate P(i/bin) (Probability of an event belonging to model i to be in a certain bin.)
  //Loop over bins to get the total number of events belonging to each model component:

  for (int ebin=firstebin; ebin<nebins; ebin++){
      for (int i=0; i<nxbins; i++){
	for (int j=0; j<nybins; j++){
	  if (obs_data[ebin][i][j] < 10) continue;
	  p_dm[ebin][i][j]/=model_dm;
	  P[0]+=p_dm[ebin][i][j];
	  for (int ii=0; ii<N; ii++) {
	    p_bkg[ii][ebin][i][j]/=model_bkg[ii];
	    P[ii+1]+=p_bkg[ii][ebin][i][j];
	  }
	}
      }
    }  
  p_compbin.clear();
  p_compbin.push_back(p_dm); 
  for (int ii=0; ii<N; ii++) p_compbin.push_back(p_bkg[ii]);
}

////////////////////////////////////////////////////////////////

///////////////FUNCTIONS TO CALCULATE LIKELIHOODS///////////////

Number logL(V &Kpars){
  double sumlogL=0;
  for (int ebin=firstebin; ebin<nebins; ebin++){
    Number log_per_ebin=0;
    for (int i=0; i<nxbins; i++){
      for (int j=0; j<nybins; j++){
	double n = obs_data[ebin][i][j];	
	double model = Kpars[0]*dm_data[ebin][i][j];
	for (int ii=0; ii<N; ii++) model+=Kpars[ii+1]*data_bkg[ii][ebin][i][j];
	//if (n<1e-15) n=1e-10;
	if (model<1e-15) model=1e-10;
	//if (model<1e-100) continue;
	Number loglike=0;
	if (n<1e-15) loglike = -model;
	else loglike = -model+n*log(model)-n*log(n)+n;
	//loglike = -model+n*log(model)-n*log(n)+n;
	sumlogL+=loglike;			
      }
    }
  }
  
  return sumlogL;  
}

double logL_err(V &Kpars, V &Kerrors){
  init(Kerrors,N+1);
  double sumlogL=0;
  for (int ebin=firstebin; ebin<nebins; ebin++){
    Number log_per_ebin=0;
    for (int i=0; i<nxbins; i++){
      for (int j=0; j<nybins; j++){
	double n = obs_data[ebin][i][j];	
	double model = Kpars[0]*dm_data[ebin][i][j];
	//if (n<1e-15) n=1e-10;
	for (int ii=0; ii<N; ii++) model+=Kpars[ii+1]*data_bkg[ii][ebin][i][j];
	if (model<1e-15) model=1e-10;
	//if (model<1e-100) continue;
	double loglike = 0;
	if (n<1e-15) loglike = -model;
	else loglike = -model+n*log(model)-n*log(n)+n;
	//double loglike = -model+n*log(model)-n*log(n)+n;
	log_per_ebin+=loglike;
	sumlogL+=loglike;
	double n_dm = Kpars[0]*dm_data[ebin][i][j];
	//if (n_dm<1e-10) n_dm=0.000000001;
	double err_dm = (n*n_dm*n_dm)/(model*model);
	Kerrors[0] += err_dm;
	for (int ii=0; ii<N; ii++){
	  double n_model = data_bkg[ii][ebin][i][j];
	  //if (n_model<1e-10) n_model=0.00000001;
	  double error = (n*n_model*n_model)/(model*model);
	  Kerrors[ii+1] += error;
	}
			
      }
    }
  }
  for (int ii=0; ii<N+1; ii++) Kerrors[ii]=sqrt(1./Kerrors[ii]);
  return sumlogL;  
}

Number logL_gausprior(int ebin,V &Kpars, V &Kerrors, V &Kpars_perbin){
  cout << "OUTDATED" << endl;
  return -1;
}

Number logL_gausprior(V &Kpars, V &Kerrors){
  Number sumlogL=0;
  for (int ebin=firstebin; ebin<nebins; ebin++){
    for (int i=0; i<nxbins; i++){
      for (int j=0; j<nybins; j++){
	Number n = obs_data[ebin][i][j];
	Number model = Kpars[0]*dm_data[ebin][i][j];
	for (int ii=0; ii<N; ii++) model+=Kpars[ii+1]*data_bkg[ii][ebin][i][j];
	if (model<1e-15) model=1e-10;
	Number loglike=0;
	if (n<1e-15) loglike=-model;
	else loglike = -model+n*log(model)-n*log(n)+n;
	sumlogL+=loglike;
	double n_dm = Kpars[0]*dm_data[ebin][i][j];
	double err_dm = (n*n_dm*n_dm)/(model*model);
	Kerrors[0] += err_dm;
	for (int ii=0; ii<N; ii++){
	  double n_model = data_bkg[ii][ebin][i][j];
	  double error = (n*n_model*n_model)/(model*model);
	  Kerrors[ii+1] += error;
	}
      }
    }
  }
  Number gaus_prior=0;
  for (int ii=0; ii<N; ii++){
    Number sigma = 10*Kerrors[ii+1];
    gaus_prior += -pow(Kpars[ii+1]-1,2)/(2*sigma*sigma) - log(sqrt(2*PI)*sigma);
    sumlogL+=gaus_prior;
  }
  
  return sumlogL;
}

///////////////FUNCTIONS FOR MINIMIZERS/MAXIMIZERS/////////////////////////////////

Number Conjugate_Gradients(V &Kpars){ //First aproximation to the minimum using CG

init(Mat,N+1,N+1);
    
    for (int ebin=firstebin; ebin<nebins; ebin++){
      for (int i=0; i<nxbins; i++){
	for (int j=0; j<nybins; j++){
	  if (obs_data[ebin][i][j] < 10) continue;
	  for (int ii=0; ii<N+1; ii++){
	    for (int jj=0; jj<N+1; jj++){

	      Mat[ii][jj]+=p_compbin[ii][ebin][i][j]*p_compbin[jj][ebin][i][j]*total_data/obs_data[ebin][i][j];	   
	    }
	  }
	}
      }
    }

    for (int ii=0; ii<N+1; ii++) {
      for (int jj=0; jj<N+1; jj++) {
	if (Mat[ii][jj] < 1e-20) Mat[ii][jj]  = 1e-50;
      }
    }

    Number det;
    V W;
    init(W,N+1);
    for (int ii=0; ii<N+1; ii++) W[ii] = 1./(N+1);
    
    solvecholesky(Mat,P,W,det);
    Kpars[0] = W[0]*total_data/model_dm;
    for (int ii=1; ii<N+1; ii++) 
      {
	Kpars[ii] = W[ii]*total_data/model_bkg[ii-1];
      }
    V Kerrors;
    init(Kerrors,N+1);
  
    Number MaxlogL = logL(Kpars);
    //    Number MaxlogL = logL_gausprior(Kpars,Kerrors);
    return MaxlogL;
}



// Functions for Expectation-Maximization

double DM_estimate(VM &input, VM &output, vector<VM> &back){
  //This function estimates the predicted number of events for each model component.
  // input: Is the data
  // output: Is the DM predicted number of events in each bin
  // back: Is the background components number of events in each bin
  _loglike = 0;
  output = p_dm; 
  back = p_bkg;
  double sum = 0;
  double totalcounts=0;
  N_dm=0;
  init(N_bkg,N);
  //Loop over bins.
  for (int ebin=firstebin; ebin<nebins; ebin++){
    for (int i=0; i<nxbins; i++){
      for (int j=0; j<nybins; j++){
	double n = input[ebin][i][j];//Counts content in the bin
	double n_of_dm = output[ebin][i][j]*A*n;
	double denom = A*p_dm[ebin][i][j];
	
	V n_of_bkg; init(n_of_bkg,N);
	for (int ii=0; ii<N;ii++) {
	  n_of_bkg[ii] = back[ii][ebin][i][j]*Pars[ii]*n;
	  denom+=Pars[ii]*p_bkg[ii][ebin][i][j];
	}
	if (denom < 1e-10) denom=1e-10;
	n_of_dm = n_of_dm/denom;
	double pred = n_of_dm;
	N_dm+=n_of_dm;
	totalcounts+=n_of_dm;
	output[ebin][i][j]=n_of_dm;
	for (int ii=0; ii<N;ii++) {
	  n_of_bkg[ii] = n_of_bkg[ii]/denom;
	  N_bkg[ii]+=n_of_bkg[ii];
	  back[ii][ebin][i][j]=n_of_bkg[ii];
	  pred+=n_of_bkg[ii];
	  totalcounts+=n_of_bkg[ii];
	}	
	double loglike = -pred+n*log(pred)-n*log(n)+n;
	_loglike+=loglike;
	sum+=n;
      }
    }    
  }
  
  return sum;

}

double  EMupdate_pars(V &Kpars,bool fixBaryonic){
  //This function uses the estimation of the counts number for each component to recalculate the parameters
  double total_counts = N_dm;
  for (int ii=0; ii<N; ii++) {total_counts+=N_bkg[ii];}
  A = N_dm/total_counts;
  if (!fixBaryonic) for (int ii=0; ii<N; ii++){ Pars[ii] = N_bkg[ii]/total_counts;}
  
  Kpars[0] = A*total_data/model_dm; 
  for (int ii=0; ii<N; ii++) {Kpars[ii+1] = Pars[ii]*total_data/model_bkg[ii];}
  return total_counts;
}


Number Expectation_Maximization(V &Kpars, V &Kerrors, bool fixBaryonic = false){
  init(Pars,N); 
  //init(Kpars,N+1);
  double maxlogL=logL(Kpars);
  //Initialize parameters to values close to the estimated (0 for DM, 1 for bkg components);
  
  A = Kpars[0]*model_dm/total_data;
  
  for (int ii=0; ii<N; ii++){
    Pars[ii] = Kpars[ii+1]*model_bkg[ii]/total_data; 
  }   
  //Use EM to maximize likelihood:

  VM dm; //This is the predicted number of DM events
  vector<VM> back; //This is the predicted number of each background component events
  
  double tol=1e-5;
  double endiff = 100000;
 
  do{    
    double old_loglike = logL(Kpars);
    V old_Kpars = Kpars;
    DM_estimate(obs_data,dm,back);
    EMupdate_pars(Kpars,fixBaryonic);
    maxlogL = logL(Kpars);
    endiff = fabs(maxlogL-old_loglike);
    if (maxlogL <= old_loglike) {Kpars = old_Kpars; break;}
    }while(endiff>tol);
  maxlogL = logL_err(Kpars,Kerrors);
  return maxlogL;
}

Number My_Minimizer(V &Kpars, V &Kerrors, Number &maxlogL, Number tol=HUGE_VAL){
  gRandom->SetSeed(0); 
  V K_0 = Kpars;
  int niter=3000;
  V x; x = Kpars;
  Number dm_step = 1;
  Number irf_step = 0.001;
  /*Number comp_step =0.1;
    Number ps_step = 0.1;*/
  Number comp_step = 0.001;
  Number ps_step = 0.001;

  for (int i=0; i<niter; i++){
    for (int j=0; j<N+1; j++){
      Number step=0;
      if (j==0) step =-dm_step+gRandom->Uniform(2*dm_step); 
      else if (j==1) step = -irf_step+gRandom->Uniform(2*irf_step);
      else if (j==4) step = -ps_step+gRandom->Uniform(2*ps_step);
      else step = -comp_step+gRandom->Uniform(2*comp_step);
      x[j] = x[j]+step;
      if (j>0) {
	if (x[j] > K_0[j]+tol) x[j] = K_0[j]+tol; //Only for fixed bkg
	if (x[j] < K_0[j]-tol) x[j] = K_0[j]-tol; //Only for fixed bkg
      }
      Number loglike = logL(x);
      if (loglike>=maxlogL) {maxlogL = loglike; Kpars = x;}
      else {x[j] = x[j]-step;}	
    }     
  }

  Kpars = x;
  maxlogL = logL_err(Kpars,Kerrors);
  return maxlogL; 

}

Number Upper_Minimizer(V &Kpars, V &Kerrors, Number &maxlogL, Number tol=HUGE_VAL){
  gRandom->SetSeed(0);
  V K_0 = Kpars;
  int niter = 25;
  int trials = 100;
  V x; x = Kpars;
  Number dm_step = 300;
  Number irf_step = 0.001;
  Number comp_step = 0.001;
  Number ps_step = 0.001;

  Number Upperlimit = Kpars[0];
  for (int trial = 0; trial<trials; trial++) {
    x = K_0;
    Number maxdiflog = fabs(2*(maxlogL-logL(x)))-2.71;
    for (int i=0; i<niter; i++){
      for (int j=0; j<N+1; j++) {
	Number step=0;
	if (j==0){
	  if (x[j]<0) step = -gRandom->Uniform(dm_step);
	  else step = gRandom->Uniform(dm_step);
	}
	else if(j==1) step = -irf_step+gRandom->Uniform(2*irf_step);
	else if(j==4) step = -ps_step+gRandom->Uniform(2*ps_step);
	else step = -comp_step+gRandom->Uniform(2*comp_step);
	x[j] = x[j]+step;
	if (j>0) {
	  if (x[j] > K_0[j]+tol) x[j] = K_0[j]+tol; //Only for fixed bkg                             
	  if (x[j] < K_0[j]-tol) x[j] = K_0[j]-tol; //Only for fixed bkg
	}
	Number difflog = fabs(2*(maxlogL-logL(x)))-2.71;
	if (fabs(difflog)<=fabs(maxdiflog)) {
	  maxdiflog = difflog; 
	  Kpars = x;
	}
	else if (j==0) x[j] = x[j]-step;
	else  x[j] = x[j]-step; 
      }
    }
    
    if (fabs(2*(maxlogL-logL(x))) > 2.72 || fabs(2*(maxlogL-logL(x))) < 2.70) continue;
    if (fabs(x[0]) > fabs(Upperlimit)) {Upperlimit = x[0]; Kpars = x;}
  }
  //Kpars = x;
  Kpars = K_0;
  return fabs(Upperlimit); 
  
}


Number My_Minimizer_bin_by_bin(int ebin, V &Kpars, V &Kerrors, V &Kpars_perbin, Number &maxlogLbin){

  gRandom->SetSeed(0); 
  int niter=1000;
  
  V x; x =  Kpars_perbin;
  Number dm_step = 0.01;
  Number irf_step = 0.0001;
  Number comp_step = 0.001;
  
  for (int i=0; i<niter; i++){
    for (int j=0; j<N+1; j++){
      Number step=0;
      if (j==0) step =-dm_step+gRandom->Uniform(2*dm_step); 
      else if (j==1) step = -irf_step+gRandom->Uniform(2*irf_step);
      else step = -comp_step+gRandom->Uniform(2*comp_step);
      x[j] = x[j]+step;
      Number loglike = logL_gausprior(ebin,Kpars,Kerrors,x);
      if (loglike>=maxlogLbin) {maxlogLbin = loglike; Kpars_perbin = x;}
      else {x[j] = x[j]-step;}	
    }     
  }

  Kpars_perbin = x;
  maxlogLbin = logL_gausprior(ebin,Kpars,Kerrors,Kpars_perbin);
  return maxlogLbin;


}

void _minuitFunction(int& nDim, double* gout, double& result, double par[], int flg){

  V p; 
  for (int i=0; i<N+1; i++){
    p.push_back(par[i]);
  }

  result = logL(p);
  
}

Number RunMinuit(V &Kpars){

TFitter* minimizer = new TFitter(N+1);
  // MAKE IT QUIET!!
  { 
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
    }

  minimizer->SetFCN(_minuitFunction);

  for (int ii=0; ii<N+1; ii++){
    minimizer->SetParameter(ii,Form("%d",ii),Kpars[ii],100,0,0);
  }

  // minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);
  
  for (int ii=0; ii<N+1; ii++){
    Kpars[ii] = minimizer->GetParameter(ii);
  }
  
  return logL(Kpars);

}

////////////////FUNCTIONS FOR PARAMETER CALCULATIONS///////////////

Number calc_MaxlogL(V &Kpars, V &Kerrors){ //Maximum Likelihood Claculation
  Number MaxlogL = 0;

  cout << "Maximizing Likelihood............." << endl;

  MaxlogL = Conjugate_Gradients(Kpars);
  for (int ii=0; ii<N+1; ii++) {cout << Kpars[ii] << "  ";}
  cout << endl;
  for (int ii=1; ii<N+1; ii++) Kpars[ii] = 1.0; //Only for Best Case
  MaxlogL = Expectation_Maximization(Kpars,Kerrors,true); //Set true for Best Case
  for (int ii=0; ii<N+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  MaxlogL = My_Minimizer(Kpars,Kerrors,MaxlogL,0.1); //Set tol=0.1 for Best Case
  
  fprintf(stderr,"Max Likelihood:\t%lf\n",MaxlogL);

  cout << "MaxL Nuisance Parameters: " << endl;
  for (int ii=0; ii<N+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << "Nuisance Parameters errors: " << endl;
  for (int ii=0; ii<N+1; ii++) cout << Kerrors[ii] << "  ";
  cout << endl;

  return MaxlogL;
}


void calc_Correlation(V &Kpars, Number &maxlogL, M &Mcov){ //Correlation Factors Calculation
  //TNtuple *ntt=new TNtuple("nt","nt","type:logL:dm:comp");
  TNtuple *ntt=new TNtuple("nt","nt","type1:type2:logL:comp1:comp2");
  init(Mcov,N+1,N+1);
  cout << "Calculating Correlation Factors..............." << endl;
  V Kpars_min = Kpars;
  double intervals[N+1];
  /*if (firstebin==0 && nebins==ebinning) {double array[N+1] = {150,0.00025,0.25,0.5,0.02,0.06,0.015,
							      0.1,3.5,0.06,0.5,10,3,0.02,
							      0.1,80,0.01,0.1,0.1,0.12}; */

  if (firstebin==0 && nebins==ebinning) {double array[N+1] = {150,0.00025,0.25,0.5,0.02,0.06,0.015,0.1};
    
    for (int i=0; i<N+1; i++) intervals[i] = array[i];}
  
  else {double array[N+1] = {1000,0.01,5,5,1,1,0.5}; for (int i=0; i<N+1; i++) intervals[i] = array[i];}
  
  int npoints = 100; 
  
  for (int comp1=0; comp1<1/*N+1*/; comp1++){
    for (int comp2=comp1; comp2<N+1; comp2++){
      Kpars = Kpars_min; 
      Number step1 = intervals[comp1];
      Number step2 = intervals[comp2];
      if (comp1==comp2) continue;
      Number minComp1 = Kpars_min[comp1]-step1;
      Number maxComp1 = Kpars_min[comp1]+step1;
      Number minComp2 = Kpars_min[comp2]-step2;
      Number maxComp2 = Kpars_min[comp2]+step2;
      for (int i=0; i<npoints; i++){
	Number norm1 = minComp1 + i*((maxComp1-minComp1)/npoints);
	for (int j=0; j<npoints; j++){
	  Number norm2=minComp2 + j*((maxComp2-minComp2)/npoints); 
	  Kpars[comp1] = norm1;
	  Kpars[comp2] = norm2;
	  Number loglike = logL(Kpars);
	  ntt->Fill(comp1,comp2,loglike,norm1,norm2);
	}
      }
      /*ntt->Draw("logL:comp2:comp1",Form("2*(%lf-logL)<2.71  && type1==%d && type2==%d",maxlogL,maxlogL,comp1,comp2),"colz prof");
      TF2 f("func",Form("0.5*([0]*(x-%lf)*(x-%lf)+[1]*(y-%lf)*(y-%lf)+2*[2]*(x-%lf)*(y-%lf)+[3])",Kpars_min[comp1],Kpars_min[comp1],Kpars_min[comp2],Kpars_min[comp2],Kpars_min[comp1],Kpars_min[comp2]));
      TH1 *h_cor = ntt->GetHistogram();
      h_cor->SetMaximum(maxlogL);
      h_cor->Fit(&f,"Q");
      gPad->Update();
      gPad->Modified();
      //gPad->WaitPrimitive();*/ 
    }
  }
  cout << "PLEASE, CHECK THAT PLOTS ARE CORRECT!" << endl;  
  for (int ii=0; ii<N+1; ii++) cout << Kpars_min[ii] << "  ";
  cout << endl;

  V corr;

  for (int comp1=0; comp1<1/*N+1*/; comp1++){
    for (int comp2=comp1; comp2<N+1; comp2++){
      if (comp1==comp2) continue;
      ntt->Draw("logL:comp2:comp1",Form("2*(%lf-logL)<2.71  && type1==%d && type2==%d",maxlogL,maxlogL,comp1,comp2),"colz prof");
      TF2 f("func",Form("0.5*([0]*(x-%lf)*(x-%lf)+[1]*(y-%lf)*(y-%lf)+2*[2]*(x-%lf)*(y-%lf)+[3])",Kpars_min[comp1],Kpars_min[comp1],Kpars_min[comp2],Kpars_min[comp2],Kpars_min[comp1],Kpars_min[comp2]));
      TH1 *h_cor = ntt->GetHistogram();
      h_cor->SetMaximum(maxlogL);
      h_cor->Fit(&f,"Q");
      
      Number cfactor = f.GetParameter(2)/sqrt(fabs(f.GetParameter(0)*f.GetParameter(1)));
      if (comp1==0) corr.push_back(cfactor);
      
      Mcov[comp1][comp1] = f.GetParameter(0);
      Mcov[comp2][comp2] = f.GetParameter(1);
      Mcov[comp1][comp2] = f.GetParameter(2);
      Mcov[comp2][comp1] = Mcov[comp1][comp2];

      gPad->Update();
      gPad->Modified();
      gPad->WaitPrimitive();
      
    }
  }
  
  gPad->Update();
  gPad->Modified();
  TString str[N+1] ={"Irf","Leptonic","Hadronic","Diffuse 1", "Diffuse 2","Diffuse 3"}; 
  
  double binstep = (log10(100.) - log10(0.03))/ebinning;
  double binlower = pow(10,log10(0.03)+firstebin*binstep);
  double binupper = pow(10,log10(0.03)+(nebins)*binstep);
  
  /*cout << "ENERGY RANGE: " << binlower<< "-"<<binupper << " TeV" << endl;
  for (int ii=0; ii<N; ii++) cout << str[ii] << ":   " <<  corr[ii] <<endl;
  Kpars = Kpars_min;*/
  cout << masses[0] << endl;
   double energy = binlower+(binupper-binlower)/2;
  cout << energy << "  ";
  for (int ii=0; ii<N; ii++) cout << corr[ii] << "  ";
  cout << endl;

  //cfactors = corr;
}




Number f_TS(Number maxlogL,Number dmNorm, V &Kpars){
  A = dmNorm*model_dm/total_data;
  VM dm;
  vector<VM> back;
  
  double tol=1e-4;
  const int maxit = 10000;
  double diff = HUGE_VAL;
  double loglike =  10000000;

  int it=0;
  do{
    double old_loglike = loglike;
    DM_estimate(obs_data,dm,back);
    EMupdate_pars(Kpars,true);
    
    loglike = logL(Kpars);
    diff = fabs(loglike-old_loglike);
    it++;
  }while(diff > tol);
  Number TS = -2*(loglike-maxlogL)-2.71;
  return TS;  
}


void MatrixtoHisto(M input, TH2D *h){

  h->Reset();
  h->SetBins(10,0,9,10,0,9);
  for (int i=0; i<nxbins; i++){
    for (int j=0; j<nybins; j++){
      cout << input[i][j] << "  " << i << "  " << j << endl;
      h->SetBinContent(i,j,input[i][j]);
    }
  }
}

void seeProfile(V input,TH1D *h){

 h->Reset();
 h->SetBins(10,-5,5);
 for (int i=0; i<nxbins; i++) {
   h->SetBinContent(i,input[i]);
 }


}

Number get_DMUpperLimit(Number &dm_mass){
  
  FillDataM(dm_mass);

  //init(N_bkg,N);
  
  ////EM Algorithm ////
  
    // Set initial values for parametes. A is the parameter for DM and Pars[] is a vector of parameters for Background components.
    //Must satisfy: A+Sum(Pars)) = 1
    //The total number of parameters is N+1, being N the number of Backgorund components
    
    p_dm; p_dm = dm_data;
    p_bkg; p_bkg = data_bkg;
   
    V P_n;
    init(N_bkg,N); init(P,N+1); 
    //Calculate P(i/bin) (Probability of an event belonging to model i to be in a certain bin.)
    //Loop over bins to get the total number of events belonging to each model component:
        
    for (int ebin=firstebin; ebin<nebins; ebin++){
      for (int i=0; i<nxbins; i++){
	for (int j=0; j<nybins; j++){
	  if (obs_data[ebin][i][j] < 10) continue;
	  p_dm[ebin][i][j]/=model_dm;
	  P[0]+=p_dm[ebin][i][j];
	  for (int ii=0; ii<N; ii++) {
	    p_bkg[ii][ebin][i][j]/=model_bkg[ii];
	    P[ii+1]+=p_bkg[ii][ebin][i][j];
	  }
	}
      }
    }

    

    
    p_compbin.clear();
    p_compbin.push_back(p_dm); 
    for (int ii=0; ii<N; ii++) p_compbin.push_back(p_bkg[ii]);
    
    init(Mat,N+1,N+1);
    
    for (int ebin=firstebin; ebin<nebins; ebin++){
      for (int i=0; i<nxbins; i++){
	for (int j=0; j<nybins; j++){
	  if (obs_data[ebin][i][j] < 10) continue;
	  for (int ii=0; ii<N+1; ii++){
	    for (int jj=0; jj<N+1; jj++){

	      Mat[ii][jj]+=p_compbin[ii][ebin][i][j]*p_compbin[jj][ebin][i][j]*total_data/obs_data[ebin][i][j];
	      
	    }
	  }
	  
	  
	}
      }
    }

    Number det;
    V W;
    init(W,N+1);
    for (int ii=0; ii<N+1; ii++) W[ii] = 1./(N+1);   
    
    //Aproximate the Minimum using conjugate gradients

    solvecholesky(Mat,P,W,det);
    //return 0;
    V Kpars; init(Kpars, N+1);
    V Kerrors; init(Kerrors,N+1);
    //If any of the W is negative, put it to 0.
    if (W[0] < 0) W[0] = 1e-20;  

    // Compute Kpars, which depends on W

    Kpars[0] = W[0]*total_data/model_dm;
    for (int ii=1; ii<N+1; ii++) 
      {
	if (W[ii] < 0) W[ii] = 1e-20;
      Kpars[ii] = W[ii]*total_data/model_bkg[ii-1];
      }

    // Calculate MaxlogL and initial Kpars for bkg components, giving the approximate Kpars calculated before.

    double maxlogL = Expectation_Maximization(Kpars,Kerrors); //Calculates the Maximum likelihhood and updates Pars[].
    //fprintf(stderr,"%lf\n",maxlogL);
    //for (int ii=0; ii<N+1; ii++) cout << Kpars[ii] <<"  ";
    cout << "   " << maxlogL << "  " << Kpars[0] << "  " << Kpars[1] <<  endl;
    //cout << A <<"  ";
    //for (int ii=0; ii<N; ii++) cout << Pars[ii] <<"  ";
    //cout << endl;
    //cout << model_dm << "   ";
    //for (int ii=0; ii<N; ii++) cout << model_bkg[ii] <<"  ";
    //cout << total_data << endl;
    //cout << endl;

    //Now use the resulting Kpars to initialize the K pars of bkg components in the search for DM norm.
    //Secant method to calculate the zero of TS function:
    Number dmNorm = 0;
    Number x1 = Kpars[0];
    Number x2 = 100;
    if (dm_mass==1){
    x1 = -50;
    x2 = 50;
    }

    Number tol = 0.0001;

    int maxiter = 100;
    for (int i=0; i<maxiter; i++) {
      Number old_x1 = x1; 
      Number old_x2 = x2;
      Number f1 = f_TS(maxlogL,old_x1,Kpars);
      Number f2 = f_TS(maxlogL,old_x2,Kpars);
      Number diff = fabs(f1-f2);
      if (diff < 1e-10) diff = 0.000001; 

      x1 = old_x1 - f1*((old_x1-old_x2)/(f1-f2));
      x2 = old_x1;

      //cout << x1 << "  " << x2 <<"  " << f1 << "  " << f2 << endl;
      if (fabs(x1-x2) < tol) break;
    }
    dmNorm = x1;

    /*double trials = 100.;
    double start = -100.;
    double end = 100.;
    
    for (int i=0; i<trials; i++){
      
      double dmNorm = start+i*(end-start)/trials;  
      
      double TS = f_TS(maxlogL,dmNorm,Kpars);
      
      cout << "TS VALUE: " << TS << "  " << dmNorm <<  endl;
      }*/
    
    
    return dmNorm;
    //}
}


/*#include "TRandom.h"

double fakeLogL(V x){
  double sum=0;
  for(int i=0;i<x.size();i++){
    sum+=-0.5*(x[i]-3)*(x[i]-3)/3;
  }
  //  cout<<"RETURNING "<<sum<<endl;
  return sum;
}

#define logL fakeLogL


#include "TNtuple.h"
TNtuple nt("","","k0:k1:k2:lp");

void compute_all_metropolis(){

  // Obten una aproximacion al minimo
  // output es un vector de ks

  V ks;

#ifdef logL
  ks.push_back(3);
  ks.push_back(3);
  ks.push_back(3);
#endif

  for(int i=0;i<ks.size();i++) if(ks[i]<=0) ks[i]=1;


  // Comienza el metropolis

  V k_actual=ks;
  V k_maximum=ks;
  double log_prob_maximum=logL(k_actual);
  double log_prob_actual= logL(k_actual);

  // Search for the maximum
  const int maximization_steps=100;
  for(int i=0;i<maximization_steps;i++){
    V candidate=k_actual;
for(int ii=0;ii<candidate.size();ii++)  candidate[ii]=exp(log(candidate[ii]) + log(2)*gRandom->Gaus());
    double log_prob=logL(candidate);

    if(log_prob>log_prob_maximum){
      log_prob_maximum=log_prob;
      k_maximum=candidate;
    }
    
    double lu=log(gRandom->Uniform());
    if(log_prob-log_prob_actual>lu){
      k_actual=candidate;
      log_prob_actual=log_prob;
    }

nt.Fill(k_actual[0],k_actual[1],k_actual[2],log_prob_actual);    
  }
  
  
  // Search for upper limit
  const int search_steps=10000;
  const int dark_matter_component=0;
  double closer_value=HUGE_VAL;
  k_actual=k_maximum;
  log_prob_actual=log_prob_maximum;
  V k_limit=k_maximum;
  for(int i=0;i<search_steps;i++){
    V candidate=k_actual;
for(int ii=0;ii<candidate.size();ii++)  candidate[ii]=exp(log(candidate[ii]) + 0.1*gRandom->Gaus());
    double log_prob=logL(candidate);

cout<<-2*(log_prob-log_prob_maximum)<<" "<< candidate[dark_matter_component]<<" "<< k_limit[dark_matter_component]<<endl;
    if(-2*(log_prob-log_prob_maximum)<2.71 && candidate[dark_matter_component]>k_limit[dark_matter_component]){

	k_limit=candidate;
    }
    
    double lu=log(gRandom->Uniform());
    if(log_prob-log_prob_actual>lu){
      k_actual=candidate;
      log_prob_actual=log_prob;
    }
  }


cout<<"Minimum at "<<k_maximum[0]<<" "<<k_maximum[1]<<" "<<k_maximum[2]<<" "<<log_prob_maximum<<endl
<<" UPPER LIMIT "<<k_limit[0]<<" "<<k_limit[1]<<" "<<k_limit[2]<<endl;

}*/

void run(){

  double x[nmasses];
  double y[nmasses];
  double y_thermal[nmasses];

  for (int mass=0; mass < nmasses; mass++) {
    Number dmMass = masses[mass];
    Number Thermal_X = 3e-26; 
    Number dmNorm = get_DMUpperLimit(dmMass);
    
    if (dmNorm < 0) dmNorm*=-1.;

    x[mass] = dmMass;
    y[mass] = dmNorm*Thermal_X;
    y_thermal[mass] = Thermal_X;

    cout <<dmMass << "  " <<  dmNorm << "  " << dmNorm*Thermal_X << endl;
  }
  TCanvas *c =  new TCanvas();
  TMultiGraph *mg =  new TMultiGraph();
  TGraph *g = new TGraph(nmasses,x,y);
  TGraph *g_thermal = new TGraph(nmasses,x,y_thermal);
  mg->Add(g,"pl");
  mg->Add(g_thermal,"l");
  c->SetLogx();
  c->SetLogy();
  g->SetMarkerStyle(20);
  //g->Draw("apl");
  //g_thermal->Draw("same");
  mg->Draw("a");

}

void Check_Correlations(int firstbin,int howmany){
  firstebin = firstbin;
  nebins = firstebin+howmany;
  TCanvas *c = new TCanvas("multicanvas","multicanvas");
  Number dm_mass = masses[0];
  FillDataM(dm_mass); // Fill Matrixes with data
  DataSim(obs_data); 
  calc_P(); // Calculate probabilities
  
  //Minimize Kpars
  
  V Kpars; init(Kpars, N+1);
  V Kerrors; init(Kerrors,N+1);
    
  Number maxlogL = calc_MaxlogL(Kpars,Kerrors);
  M Mcov;
  calc_Correlation(Kpars,maxlogL,Mcov); 

  double binstep = (log10(100.) - log10(0.03))/20;
  double binlower = pow(10,log10(0.03)+firstebin*binstep);
  double binupper = pow(10,log10(0.03)+(nebins)*binstep);

  /*for (int ii=0; ii<N+1; ii++){
    for (int jj=0; jj<N+1; jj++){
      cout << Mcov[ii][jj] << "  ";
    }
    cout << endl;
    }*/

  /*double energy = binlower+(binupper-binlower)/2;
  cout << energy << "  ";
  for (int ii=0; ii<N; ii++) cout << cfactors[ii] << "  ";
  cout << endl;*/

}


void calc_Mcov(V &Kpars, V &Kerrors, Number &dm_mass){

  /*firstebin=0;
  nebins=firstebin+ebinning;
  Number dm_mass = masses[0];
  FillDataM(dm_mass); // Fill Matrixes with data
  calc_P(); // Calculate probabilities*/
   //Minimize Kpars
  
  //V Kpars; 
  init(Kpars, N+1);
  //V Kerrors; 
  init(Kerrors,N+1);

  Kpars[0] = 0.001;
  for (int ii=1; ii<N+1; ii++) Kpars[ii] = 1;
    
  Number maxlogL = calc_MaxlogL(Kpars,Kerrors);
  
  M Mcov;

  //Calculate Correlation Matrix
  calc_Correlation(Kpars,maxlogL,Mcov);

  TString masstr; 
  //double dm_mass = masses[0];
  if (dm_mass < 1.0){
    ostringstream os;
    os<<dm_mass*1000;
    masstr = os.str()+"GeV";      
  }
  else{
    ostringstream os;
    os<<dm_mass;
    masstr = os.str()+"TeV"; 
  }
  TString filename = "covmatrix"+masstr+".dat";

  ofstream outfile;
  outfile.open(filename);

  for (int ii=0; ii<N+1; ii++){
    for (int jj=0; jj<N+1; jj++){
      //cout << Mcov[ii][jj] << "  ";
      outfile <<  Mcov[ii][jj] << "  ";
    }
    //cout << endl;
    outfile << endl;
  }
  outfile.close();

}

Number Plot_Limits(V &Kpars,Number &dm_mass){ 
  gRandom->SetSeed(0); 
  M Mcov;init(Mcov,N+1,N+1);
  
  Number maxlogL = logL(Kpars);

  TString masstr; 
  if (dm_mass < 1.0){
    ostringstream os;
    os<<dm_mass*1000;
    masstr = os.str()+"GeV";      
  }
  else{
    ostringstream os;
    os<<dm_mass;
    masstr = os.str()+"TeV"; 
  }
  TString filename = "covmatrix"+masstr+".dat";
  ifstream openfile(filename);
  int i=0;
  while(!openfile.eof()){
    for (int j=0; j<N+1; j++){
      openfile>>Mcov[i][j];
    }
    i++;    
  }
  openfile.close();
  V K, Kmin;
  init(K,N+1);

  Kmin = Kpars;
  
  //double intervals[N+1]={500,0.0005,0.6,0.6,0.05,0.05,0.02,
  //			 0.05,0.05,0.05,0.05,0.05,0.05,0.05,
  //			 0.05,0.05,0.05,0.05,0.05,0.05};

  /*double intervals[N+1]= {100,0.0002,0.2,0.4,0.01,0.05,0.01,
			0.06,1.5,0.04,0.25,4,1,0.01,
			0.05,10,0.004,0.05,0.05,0.07};*///PARA 200, 500 GEV, 1TeV
  

  /* double intervals[N+1]= {100,0.00015,0.35,0.4,0.01,0.04,0.01,
			0.05,1.5,0.035,0.2,4,1.5,0.02,
			0.05,20,0.004,0.05,0.05,0.06}; //PARA 5TeV*/


  /* double intervals[N+1]= {250,0.0002,0.3,0.5,0.015,0.04,0.01,
			  0.06,1.5,0.04,0.25,0.5,0.5,0.015,
			  0.1,0.25,0.005,0.05,0.05,0.05}; //PARA 100 GEV*/

  /* double intervals[N+1]={55,0.0003,0.25,0.5,0.025,0.06,0.02,
			 0.1,3.5,0.1,0.5,10,3,0.03,
			 0.15,80,0.01,0.15,0.1,0.12};*/

  //double intervals[N+1]={300,0.0003,0.4,0.7,0.025,0.06,0.01,0.05};

  //  double intervals[N+1] = {300,0.0003,0.1,0.1,0.025,0.06,0.01,0.05};
  double intervals[N+1] = {150,0.0003,0.05,0.05,0.05,0.06,0.01,0.01};

  TNtuple *points1 = new TNtuple("p1","p1","logL:dm:irf:lepto:hadro:diff1:diff2:diff3:ps1");
  TNtuple *points2 = new TNtuple("p2","p2","logL:dm:ps7:ps8:ps9:ps10:ps11:ps12:ps13");
  
  
  int niters = 100000;
  
  Number upperlimit = Kpars[0];
  V Kuppers; init(Kuppers,N+1);

  int howmany=0;
  for (int iter=0; iter<niters; iter++){
    
    if (iter%10000==0)  cout<< "Iteration number: " << iter <<"," << howmany << "Points"<< endl;
    for (int i=0; i<N+1; i++){
      Number min = Kmin[i]-intervals[i];
      Number max = Kmin[i]+intervals[i];
      if (i==0 && Kmin[0] > 0) min = Kmin[i];
      if (i==0 && Kmin[0] < 0) max = Kmin[i];      
      K[i] = min+gRandom->Uniform(max-min);
      //      K[i]-=Kmin[i];
    }
    V temp;
    mult(Mcov,K,temp);
    //Number loglike = -mult(K,temp);
    
    Number loglike = 2*(maxlogL-logL(K));
    //for (int ii=0; ii<N+1; ii++)cout << K[ii] << "  ";
    //cout << endl;
    //cout << loglike <<  endl;
    //if (loglike<2.715 && loglike>2.705) {
    if (loglike<2.72 && loglike>2.70) {
      howmany++;
    //cout << loglike <<  endl;
      //for (int ii=0; ii<N+1; ii++)cout << K[ii]+Kmin[ii] << "  ";
      //cout << endl;
      //points1->Fill(loglike,K[0]+Kmin[0],K[1]+Kmin[1],K[2]+Kmin[2],K[3]+Kmin[3],K[4]+Kmin[4],K[5]+Kmin[5],K[6]+Kmin[6],K[7]+Kmin[7],K[8]+Kmin[8],K[9]+Kmin[9],K[10]+Kmin[10],K[11]+Kmin[11],K[12]+Kmin[12]);
      //points2->Fill(loglike,K[0]+Kmin[0],K[13]+Kmin[13],K[14]+Kmin[14],K[15]+Kmin[15],K[16]+Kmin[16],K[17]+Kmin[17],K[18]+Kmin[18],K[19]+Kmin[19]);
      //points1->Fill(loglike,K[0]+Kmin[0],K[1]+Kmin[1],K[2]+Kmin[2],K[3]+Kmin[3],K[4]+Kmin[4],K[5]+Kmin[5],K[6]+Kmin[6],K[7]+Kmin[7]);
      
      //points1->Fill(loglike,K[0]+Kmin[0],K[1]+Kmin[1],K[2]+Kmin[2],K[3]+Kmin[3],K[4]+Kmin[4],
      //	    K[5]+Kmin[5],K[6]+Kmin[6],K[7]+Kmin[7]);
      points1->Fill(loglike,K[0],K[1],K[2],K[3],K[4],K[5],K[6],K[7]);
      //if (fabs(K[0]+Kmin[0]) > upperlimit) {upperlimit = fabs(K[0]+Kmin[0]); Kuppers=K;}
      if (fabs(K[0]) > upperlimit) {upperlimit = fabs(K[0]); Kuppers=K;}
    }
  }
  
  cout << "UPPER LIMIT: " << upperlimit << endl;
  for (int ii=0; ii<N+1; ii++)  cout << Kuppers[ii] << "  ";
  cout << endl;
  
//for (int ii=0; ii<N+1; ii++)  cout << Kuppers[ii]+Kmin[ii] << "  ";                              
//cout << endl;


  points1->Draw("logL:dm:hadro","","prof colz");
  
  
  return upperlimit;
  
}


void Bin_by_Bin(){
  
  firstebin = 0;
  nebins = firstebin+ebinning;
  
  Number dm_mass = 0;
  FillDataM(dm_mass); // Fill Matrixes with data
  calc_P(); // Calculate probabilities
  
  V Kpars; init(Kpars,N+1);
  V Kerrors; init(Kerrors,N+1);
    
  Number maxlogL = calc_MaxlogL(Kpars,Kerrors);
  Number maxlogL_bin_by_bin=0;
  for (int ebin = 0; ebin<20; ebin++){

    firstebin=ebin;
    nebins = firstebin+1;
    calc_P(); // Calculate probabilities
    V Kpars_perbin; init(Kpars_perbin,N+1);
    V Kerrors_perbin; init(Kerrors_perbin,N+1);
    
    Kpars_perbin = Kpars;
    Number maxlogLbin = logL_gausprior(firstebin,Kpars,Kerrors,Kpars_perbin);
    cout  << maxlogLbin << endl;
    
    maxlogLbin = My_Minimizer_bin_by_bin(firstebin,Kpars,Kerrors,Kpars_perbin,maxlogLbin);
    
    cout  << maxlogLbin << endl;
    
    for (int ii=0; ii<N+1; ii++) cout << Kpars_perbin[ii] <<"  ";
    cout << endl;
    maxlogL_bin_by_bin+=maxlogLbin;
  }

  cout << maxlogL << "  " << maxlogL_bin_by_bin << endl;
  
}



void test(Number &dm_mass){
  
  firstebin=0;
  nebins=firstebin+ebinning;
  //Number dm_mass = masses[0];
  FillDataM(dm_mass); // Fill Matrixes with data    

  //Build wisely the name of the file
  TString masstr;
  if (dm_mass < 1.0){
    ostringstream os;
    os<<dm_mass*1000;
    masstr = os.str()+"GeV";      
  }
  else{
    ostringstream os;
    os<<dm_mass;
    masstr = os.str()+"TeV"; 
  }

  TString filename = "limit"+masstr+".dat";

  ofstream outfile;
  outfile.open(filename,ios::app);

  for (int i=0; i<100; i++){
  DataSim(obs_data); //Simulate a random realization of datax
  calc_P(); // Calculate probabilities
  
  V Kpars; init(Kpars,N+1);
  V Kerrors;init(Kerrors,N+1);

  //calc_Mcov(Kpars,Kerrors,dm_mass); //Calculate covariance matrix
  Number maxlogL = calc_MaxlogL(Kpars,Kerrors);
  Number UpperLimit = Plot_Limits(Kpars,dm_mass);
  outfile << UpperLimit << endl;
  }
  outfile.close();
 
}

Number calc_UpperLimit(Number &dm_mass, bool bestcase=false,double tol=0.01){
  
  firstebin = 0; 
  nebins = 20; 
  FillDataM(dm_mass); 
  DataSim(obs_data); 
  calc_P(); 
  V Kpars; init(Kpars, N+1);  
  V Kerrors; init(Kerrors,N+1);
  cout << "Maximizing Likelihood............." << endl;
  Conjugate_Gradients(Kpars); 
  Number Upperlimit=0;
  if (bestcase){
    for (int ii=1; ii<N+1; ii++) Kpars[ii] = 1.0; 
    Expectation_Maximization(Kpars,Kerrors,true);
    Number maxlogL = logL(Kpars);
    maxlogL = My_Minimizer(Kpars,Kerrors,maxlogL,tol);
    fprintf(stderr,"Max Likelihood:\t%lf\n",maxlogL);
    cout << "MaxL Nuisance Parameters: " << endl;
    for (int ii=0; ii<N+1; ii++) {cout << Kpars[ii] << "  ";}                                      
    cout << endl;
    Upperlimit = Upper_Minimizer(Kpars,Kerrors,maxlogL,tol);
    cout << "Upperlimit: " << Upperlimit << endl;
  }
  
  else{
    Expectation_Maximization(Kpars,Kerrors,false);
    Number maxlogL = logL(Kpars);
    maxlogL = My_Minimizer(Kpars,Kerrors,maxlogL);
    fprintf(stderr,"Max Likelihood:\t%lf\n",maxlogL);
    cout << "MaxL Nuisance Parameters: " << endl;
    for (int ii=0; ii<N+1; ii++) {cout << Kpars[ii] << "  ";}                                      
    cout << endl;
    Upperlimit = Upper_Minimizer(Kpars,Kerrors,maxlogL);
    cout << "Upperlimit: " << Upperlimit << endl;
  }
  return Upperlimit;
}

void Bands(Number &dm_mass, bool bestcase=true, Number tol=0.01){

  int realizations = 10;
  
  //Build wisely the name of the file
  TString masstr;
  if (dm_mass < 1.0){
    ostringstream os;
    os<<dm_mass*1000;
    masstr = os.str()+"GeV";      
  }
  else{
    ostringstream os;
    os<<dm_mass;
    masstr = os.str()+"TeV"; 
  }
  
  TString filename;
  if (!bestcase) filename = "limit"+masstr+".dat";
  else filename = "limit"+masstr+"_BestCase.dat";
  
  ofstream outfile;
  outfile.open(filename,ios::app);
  
  for (int i=0; i<realizations; i++){
    Number upperlimit = calc_UpperLimit(dm_mass,bestcase,tol);
    outfile << upperlimit << endl;
  }
  outfile.close(); 
  
  
  
}


void checksimu(){

  firstebin=0;
  nebins=firstebin+ebinning;
  Number dm_mass = masses[0];
  FillDataM(dm_mass);
  
  VM simu1;
  VM simu2;
  simu1=obs_data;
  simu2=obs_data;
  DataSim(simu1);
  //DataSim(simu2);

  TNtuple *nt = new TNtuple("nt","nt","ebin:x:y:obs:simu");

  Number chi2=0;

  for (int i=0; i<ebinning; i++) {
    for (int j=0; j<nxbins; j++){
      for (int k=0; k<nybins; k++){
	
	if (simu1[i][j][k] > 0){
	  chi2+=pow(simu2[i][j][k]-simu1[i][j][k],2)/simu1[i][j][k];}
	nt->Fill(i,j,k,simu1[i][j][k],simu2[i][j][k]);
	
	
      }
    }
  }
cout << chi2 << endl;
}


