#include "/home/queenmab/GitHub/Math/matrixes.h"
#include "DMLimitsLib.h"
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
#include <ctime>

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

void Init()
{ //Set Initial constants
  cout << "Initialize constants" << endl;
  cout << endl;
  cout << "USAGE: Init(int nebins,int nxbins,int nybins,int Nex,int Nps)" << endl;
}

void Init(int ne,int nx,int ny,int Nex,int Nps)
{
  nebins = ne;
  nxbins = nx;
  nybins = ny;

  N_ext = Nex;
  N_ps = Nps;
  
  Nbar = Nex+Nps;
  
}

void ReadFits()
{
  cout << "Read an Image contained in a fits file into a Vector of Matrixes (VM)" << endl;
  cout << endl;
  cout << "USAGE: ReadFits(VM &data, TString filename)" << endl;
  cout << endl;
  cout << "data: Empty Vector of Matrixes to store fits file binned image" << endl;
  cout << "filename: Fits file name" << endl;
  cout << endl;
}

Number ReadFits(VM &data,
              TString filename)
{
  data.clear();
  //Open the model component fits file
  TFITSHDU *hdu =  new TFITSHDU(filename);
  
  if (hdu == 0) 
    {
        printf("ERROR: could not access the HDU\n"); 
        return -1;
    }

  Number TotalCounts = 0;
  for (int bin=0; bin<nebins; bin++)
    {
      TMatrixD *mat =  hdu->ReadAsMatrix(bin); 
      if (bin==0)
        {
	  nxbins = mat->GetNcols(); 
	  nybins = mat->GetNrows();
        }
      
      M bin_data;
        for (int i=0; i<nxbins;i++)
	  {
            V column;
            for (int j=0; j<nybins; j++)
            {
	      column.push_back((*mat)(i,j));
	      TotalCounts += (*mat)(i,j);
            }
            bin_data.push_back(column);
        }
        data.push_back(bin_data);
        delete mat;
    }
  delete hdu;
  return TotalCounts;
}

void FillContainer_Bkg()
{
  cout << "Fill the container of Background models." << endl;
  cout << endl;
  cout << "USAGE: FillContainer_Bkg(TString ext[], TString> ps[], TString suf)" << endl;
  cout << endl;
  cout << "ext: Vector of strings containing the names of the baryonic extended components (Irf, Leptonoic, Hadronic, etc.)" << endl;
  cout << "ps: Vector of strings containing the names of the baryonic point sourcces." << endl;
  cout << "suf: String with suffix of fits file" << endl;
}

void FillContainer_Bkg(TString ext[], TString ps[], TString suf)
{
  N_bkg.clear();
  Bkg_model.clear();
  
  for (int i=0; i<N_ext; i++)
    {
      if (N_ext==0) break;
      TString filename = dir + ext[i] + "/modcube_LMC_" + ext[i] + suf + ".fits";

      VM data;
      N_bkg.push_back(ReadFits(data,filename));
      Bkg_model.push_back(data);
    }

  for (int i=0; i<N_ps; i++)
    {
      if (N_ps==0) break;
      TString filename =  dir_PS + "/modcube_LMC_" + ps[i] + suf + ".fits";
      VM data;
      N_bkg.push_back(ReadFits(data,filename));
      Bkg_model.push_back(data);
    }
}

void FillContainer_DM()
{
  cout << "Fill the container of DM model." << endl;
  cout << endl;
  cout << "USAGE: FillContainer_DM(Number dm_mass, TString particle, TString suf)" << endl;
  cout << endl;
  cout << "dm_mass: Dark Matter particle mass in TeV" << endl;
  cout << "particle: Final SM state of DM annihilation (W,b,Tau...)" << endl;
  cout << "suf: String with suffix of fits file" << endl;
}

void FillContainer_DM(Number dm_mass, TString particle, TString suf)
{
  N_dm = 0;
  DM_model.clear();
  TString masstr;
  
  //Build wisely the name of the file model from DM
  if (dm_mass < 1.0)
    {
      ostringstream os;
      os << dm_mass*1000;
      masstr = os.str() + "GeV";
    }
    else
      {
        ostringstream os;
        os << dm_mass;
        masstr = os.str() + "TeV"; 
      }
  
    TString modelname = "dm_LMC_" + particle + masstr;
    TString filename = dir_DM + "/modcube_" + modelname + suf + ".fits";

    N_dm = ReadFits(DM_model,filename);
    
}

void FillContainer_Obs()
{
   cout << "Fill the container of the Observed data" << endl;
  cout << endl;
  cout << "USAGE: FillContainer_Obs(TString obsname)" << endl;
  cout << endl;
}

void FillContainer_Obs(TString obsname, bool modcube, TString suf)
{
  Ntotal = 0;
  Obs_data.clear();
  TString filename;
  if (modcube)
    {
      filename = dir + obsname + "/modcube_LMC_" + obsname + suf + ".fits";
    }
  else
    {
      filename = dir + obsname + "/cntcube_LMC_" + obsname + suf + ".fits";     
    }
  Ntotal = ReadFits(Obs_data,filename);
  
}

void DataSim()
{
  cout << "Simulates a Poisson realization of data from a Modelcube" << endl;
  cout << endl;
  cout << "USAGE: DataSim(VM &data)" << endl;
}

Number DataSim(VM &data)
{
  gRandom -> SetSeed(0);
  int ebin = data.size();
  int xbin = data[0].size();
  int ybin = data[0][0].size();
  
  Number TotalCounts = 0;
  for (int i=0; i<ebin; i++)
    {
      for (int j=0; j<xbin; j++)
        {
	  for (int k=0; k<ybin; k++)
            {
	      Number mean = data[i][j][k];
	      data[i][j][k] = gRandom->Poisson(mean);
	      TotalCounts += data[i][j][k];
            }
        }
    }
  return TotalCounts;
}

void logL()
{
  cout << "Calculates Poisson logLikelihood for a given set of Normalization parameters" << endl;
  cout << endl;
  cout << "USAGE: logL(V Kpars)" << endl;
}

Number logL(V Kpars,int firstebin,int nebins)
{
    Number sumlogL=0;
    for (int ebin=firstebin; ebin<nebins; ebin++)
    {
        for (int i=0; i<nxbins; i++)
        {
            for (int j=0; j<nybins; j++)
            {
                Number n = Obs_data[ebin][i][j];
                Number model = Kpars[0]*DM_model[ebin][i][j];
                for (int ii=0; ii<Nbar; ii++) 
                {
                    model+=Kpars[ii+1]*Bkg_model[ii][ebin][i][j];
                }
                if (model<1e-15) 
                {
                    model=1e-10;
                }
                Number loglike=0;
                if (n<1e-15) 
                {
                    loglike = -model;
                }
                else 
                {
                    loglike = -model+n*log(model)-n*log(n)+n;
                }
                sumlogL+=loglike;
            }
        }
    }
    return sumlogL;  
}

/////////////////MINIMIZERS/////////////////////////////////////////

/*-------------CONJUGATE GRADIENTS--------------------*/


void calc_P(V &P, vector<VM> &p_compbin)
{
  p_compbin.clear();
  P.clear();
  VM p_dm = DM_model;
  vector<VM> p_bkg = Bkg_model;
  
  init(P,Nbar+1);
  for (int ebin=0; ebin<nebins; ebin++)
    {
      for (int i=0; i<nxbins; i++)
        {
	  for (int j=0; j<nybins; j++)
            {
	      if (Obs_data[ebin][i][j] < 10)
                {
		  continue;
                }
	      p_dm[ebin][i][j]/=N_dm;
	      P[0]+=p_dm[ebin][i][j];
	      for (int ii=0; ii<Nbar; ii++) 
		{
		  p_bkg[ii][ebin][i][j]/=N_bkg[ii];
		  P[ii+1]+=p_bkg[ii][ebin][i][j];
                }
            }
        }
    }
  
  p_compbin.push_back(p_dm); 
  for (int ii=0; ii<Nbar; ii++) 
    {
      p_compbin.push_back(p_bkg[ii]);
    }
}

void Conjugate_Gradients()
{
  cout << "Calculates an approximation of the maximum logLikelihood using Conjugate Gradients method" << endl;
  cout << endl;
  cout << "USAGE: Conjugate_Gradients(V &Kpars)" << endl;
}

Number Conjugate_Gradients(V &Kpars)
{ //First aproximation to the minimum using CG
  V P;
  vector<VM> p_compbin; 
  calc_P(P,p_compbin);
  
  M Mat;
  init(Mat,Nbar+1,Nbar+1);
  
  for (int ebin=0; ebin<nebins; ebin++)
    {
      for (int i=0; i<nxbins; i++)
        {
	  for (int j=0; j<nybins; j++)
            {
	      if (Obs_data[ebin][i][j] < 10)
                {
		  continue;
                }
	      for (int ii=0; ii<Nbar+1; ii++)
                {
		  for (int jj=0; jj<Nbar+1; jj++)
                    {
		      Mat[ii][jj]+= p_compbin[ii][ebin][i][j]*
			p_compbin[jj][ebin][i][j]*
			Ntotal/Obs_data[ebin][i][j];
                    }
                }
            }
        }
    }
  
  for (int ii=0; ii<Nbar+1; ii++) 
    {
      for (int jj=0; jj<Nbar+1; jj++) 
        {
	  if (Mat[ii][jj] < 1e-20) 
            {
	      Mat[ii][jj]  = 1e-50;
            }
        }
    }
  Number det;
  V W;
  init(W,Nbar+1);
  for (int ii=0; ii<Nbar+1; ii++) 
    {
      W[ii] = 1./(Nbar+1);
    }
  
  solvecholesky(Mat,P,W,det);
  
  Kpars[0] = W[0]*Ntotal/N_dm;
  for (int ii=1; ii<Nbar+1; ii++) 
    {
      Kpars[ii] = W[ii]*Ntotal/N_bkg[ii-1];
    }
    
  Number MaxlogL = logL(Kpars,0,nebins);
  // Return
  return MaxlogL;
}

Number EM_Estimate(VM &data, VM &dm, vector<VM> &back,Number &Nest_dm, V &Nest_bkg, Number DMnorm, V P)
{
    //This function estimates the predicted number of events for each model component.
  dm = DM_model;
  back = Bkg_model;
  Number sumlogL = 0;
  init(Nest_bkg, Nbar);
  //Loop over bins.
  for (int ebin=0; ebin<nebins; ebin++)
    {
      for (int i=0; i<nxbins; i++)
        {
	  for (int j=0; j<nybins; j++)
            {
	      Number n = data[ebin][i][j];//Counts content in the bin
	      Number n_of_dm = (DM_model[ebin][i][j]/N_dm)*DMnorm*n;
	      Number denom = DMnorm*DM_model[ebin][i][j]/N_dm;
	      V n_of_bkg; 
	      init(n_of_bkg,Nbar);
	      for (int ii=0; ii<Nbar;ii++) 
                {
		  n_of_bkg[ii] = (Bkg_model[ii][ebin][i][j]/N_bkg[ii])*P[ii]*n;
		  denom+=P[ii]*Bkg_model[ii][ebin][i][j]/N_bkg[ii];
                }
	      if (denom < 1e-10) 
                {
		  denom=1e-10;
                }
	      n_of_dm = n_of_dm/denom;
	      
	      Number pred = n_of_dm;
	      Nest_dm+=n_of_dm;
	      
	      dm[ebin][i][j]=n_of_dm;
	      for (int ii=0; ii<Nbar;ii++) 
                {
		  n_of_bkg[ii] = n_of_bkg[ii]/denom;
		  Nest_bkg[ii]+=n_of_bkg[ii];
		  back[ii][ebin][i][j]=n_of_bkg[ii];
		  pred+=n_of_bkg[ii];
		  
                }
	      Number loglike = -pred+n*log(pred)-n*log(n)+n;
	      sumlogL+=loglike;
	    }
        }
    }
  // Return
  return sumlogL;
}

Number EM_Update_pars(Number Nest_dm, V Nest_bkg,V &Kpars, Number &DMnorm, V &P)
{
    // This function uses the estimation of the counts number for each component
    // to recalculate the parameters
    Number total_counts = Nest_dm;
    for (int ii=0; ii<Nbar; ii++) 
    {
      total_counts+=Nest_bkg[ii];
    }
    
    DMnorm = Nest_dm/total_counts;
    for (int ii=0; ii<Nbar; ii++)
      {
	P[ii] = Nest_bkg[ii]/total_counts;
      }
    Kpars[0] = DMnorm*Ntotal/N_dm;
    for (int ii=0; ii<Nbar; ii++) 
    {
         Kpars[ii+1] = P[ii]*Ntotal/N_bkg[ii];
    }
    // Return
    return total_counts;
}

void Expectation_Maximization()
{
  cout << "Approximates Maximum logLikelihood using Expectation Maximization methods" << endl;
cout << endl;
 cout << "USAGE: Expectation_Maximization(V &Kpars)" << endl;
}

Number Expectation_Maximization(V &Kpars)
{
  Number maxlogL=logL(Kpars,0,nebins);
  // Initialize parameters to values close to the estimated 
  
  Number DMnorm;
  V P; init(P,Nbar); 
  DMnorm = Kpars[0]*N_dm/Ntotal;
  for (int ii=0; ii<Nbar; ii++)
    {
      P[ii] = Kpars[ii+1]*N_bkg[ii]/Ntotal;
    }
  
  //Use EM to maximize likelihood:
  
  VM dm; //This is the predicted number of DM events
  //This is the predicted number of each background component events
  vector<VM> back; 
  
  Number tol    = 1e-5;
  Number endiff = 100000;

  Number Nest_dm;
  V Nest_bkg;
  
  do
    {    
      Number old_loglike = logL(Kpars,0,nebins);
      V old_Kpars        = Kpars;
      EM_Estimate(Obs_data,dm,back,Nest_dm,Nest_bkg,DMnorm,P);
      EM_Update_pars(Nest_dm,Nest_bkg,Kpars,DMnorm,P);
      maxlogL            = logL(Kpars,0,nebins);
      endiff             = fabs(maxlogL-old_loglike);
      if (maxlogL <= old_loglike) 
        {
	  Kpars = old_Kpars; 
	  break;
	}
    }
  while(endiff>tol);
  maxlogL = logL(Kpars,0,nebins);
  // Return
  return maxlogL;
}

void My_Minimizer()
{
  cout << "Finds Maximum likelihood using random step approximations" << endl;
  cout << endl;
  cout << "USAGE: My_Minimizer(V &Kpars, Number tol)" << endl;
}

Number My_Minimizer(V &Kpars,  
                    Number tol)
{
    // Randomly scan normalization parameter space in order to find the 
    // maximum likelihood
  
  Number maxlogL = logL(Kpars,0,nebins);
  
  gRandom          -> SetSeed(0); 
  V K_0            = Kpars;
  int niter        = 5000;
  V x; 
  x                = Kpars;
  Number dm_step   = 10;
  // Each baryonic sources is optimized for varying within a given range
  // (This was manually checked)
  Number irf_step  = 0.001;
  Number comp_step = 0.001;
  Number ps_step = 0.001;
  Number steps[Nbar+1]={10,0.001,1,0.01,0.001,1,0.001,0.001};
  if (irf_step > tol) irf_step = tol;
  if (comp_step > tol) comp_step = tol;
 
   // May be interesting to implement a convergence condition and leave
    //the loop if convergence is achieve, i.e. loglike - maxlogL < convergence
   
  for (int i=0; i<niter; i++)
    {
      for (int j=0; j<Nbar+1; j++)
        {
	  Number step=0;
	  /*if (j==0) 
            {
	      step = -dm_step + gRandom->Uniform(2*dm_step); 
            }
	  else if (j==1) 
            {
	      step = -irf_step + gRandom->Uniform(2*irf_step);
            }
	  else if(j>6)
	    {
	      step = -ps_step + gRandom->Uniform(2*ps_step);
	    }
	  else 
            {
	      step = -comp_step+gRandom->Uniform(2*comp_step);
	      }*/
	  step = -steps[j]+gRandom->Uniform(2*steps[j]);
	  x[j] = x[j]+step;
	  if (j>0) 
            {
	      if (x[j] > K_0[j]+tol) 
                {
		  x[j] = K_0[j]+tol;
		}
	      if (x[j] < K_0[j]-tol) 
                {
		  x[j] = K_0[j]-tol; 
                }
            }
	  
	  Number loglike = logL(x,0,nebins);
	  if (loglike >= maxlogL) 
	    {
	      maxlogL = loglike; 
	      Kpars   = x;
	    }
	  else 
	    {
	      x[j] = x[j]-step;
	    }
        }
    }
  maxlogL = logL(Kpars,0,nebins);
  // Return
  return maxlogL; 
}

void calc_MaxlogL()
{
  cout << "Calculate Maximum logLikelihood applying three minimization methods" << endl;
  cout << endl;
  cout << "USAGE: calc_MaxlogL(V &Kpars)" << endl;
  
}
Number calc_MaxlogL(V &Kpars,bool BestCase,Number tol)
{
  
  Number MaxlogL = 0;
  Conjugate_Gradients(Kpars);
  Expectation_Maximization(Kpars);
  if (BestCase)
    {
      for (int ii=1; ii<Nbar+1; ii++) Kpars[ii] = 1.;
    }
  MaxlogL = My_Minimizer(Kpars,tol);

  return MaxlogL;
  
}

void Upper_Minimizer()
{
  cout << "Finds the maximum DM normalization parameter that is compatible with a Likelihood ratio of 2.71" << endl;
  cout << endl;
  cout << "USAGE: Upper_Minimizer(V &Kpars, Number tol)" << endl;
}

Number Upper_Minimizer(V &Kpars,   
                       Number tol)
{
  //This function searches for the maximum DM normalization value that gives 
  // a likelihood ratio = 2.71 (Say, the upper limit)
  Number maxlogL = logL(Kpars,0,nebins);
  gRandom          -> SetSeed(0);
  V K_0             = Kpars;
  int niter         = 300;
  int trials        = 30;
  V x               = Kpars;
  Number dm_step    = 50;
  Number irf_step   = 0.001;
  Number comp_step  = 0.001;

  if (irf_step > tol) irf_step = tol;
  if (comp_step > tol) comp_step = tol;
  
  Number Upperlimit = Kpars[0];
  for (int trial = 0; trial<trials; trial++) 
    {
      x = K_0;
      //for (int ii=1; ii<Nbar+1; ii++) x[ii] = K_0[ii]; 
      Number maxdiflog = fabs(2*(maxlogL-logL(x,0,nebins)))-2.71;
      for (int i=0; i<niter; i++)
        {
	  for (int j=0; j<Nbar+1; j++) 
	    {
                Number step = 0;
                if (j==0)
		  {
                    if (x[j]<0) 
		      {
                        step = -gRandom->Uniform(dm_step);
		      }
                    else 
		      {
                        step = gRandom->Uniform(dm_step);
		      }
		  }
                else if(j==1) 
		  {
                    step = -irf_step+gRandom->Uniform(2*irf_step);
		  }
		else 
		  {
                    step = -comp_step+gRandom->Uniform(2*comp_step);
		  }
                x[j] = x[j]+step;
                if (j>0) 
		  {
                    if (x[j] > K_0[j]+tol) 
		      {
			x[j] = K_0[j]+tol; //Only for fixed bkg
		      }
                    if (x[j] < K_0[j]-tol) 
		      {
                        x[j] = K_0[j]-tol; //Only for fixed bkg
		      }
		  }
                Number difflog = fabs(2*(maxlogL-logL(x,0,nebins)))-2.71;
                if (fabs(difflog)<=fabs(maxdiflog)) 
		  {
                    maxdiflog = difflog;
                    Kpars     = x; 
// ========================================================================== 
                    Number loglike = logL(x,0,nebins);
                    if (fabs(x[0]) > fabs(Upperlimit))
                    {
                        if (fabs(2*(maxlogL-loglike))>2.715 || 
                            fabs(2*(maxlogL-loglike))<2.705)
                        {
                            continue;
                        }
                        if (fabs(x[0]) > fabs(Upperlimit)) 
                        {
                            Upperlimit = x[0]; 
                            Kpars      = x;
                            cout << "trial (within iter)= " << trial << "  " 
                                 << Kpars[0] << "  " << Kpars[1] << endl;
                            continue;
                        }
                    }
// ========================================================================== 
		  }
                else if (j==0) 
		  {
                    x[j] = x[j]-step;
		  }
                else  
		  {
                    x[j] = x[j]-step; 
		  }
	      }
        }
    }
  
  Kpars = K_0;
  // Return
  return fabs(Upperlimit); 
}

void calc_CorrFactors()
{
  cout << "Calculate the Correlation Factor between Dark Matter component and Baryonic Sources components" << endl;
  cout << endl;
  cout << "USAGE: calc_CorrFactors(V Kpars, Number intervals[], V &Cfactors" << endl;
}

void calc_CorrFactors(V Kpars, Number intervals[], V &Cfactors)
{
  TNtuple *ParSpace = new TNtuple("ParSpace","ParSpace","type1:type2:logL:comp1:comp2");
  
  Number maxlogL = logL(Kpars,0,nebins);
  
  int npoints = 100; 

  for (int comp1=0; comp1<1/*Nbar+1*/; comp1++)
    {
      for (int comp2=comp1; comp2<Nbar+1; comp2++)
        {
	  V K = Kpars;
	  Number step1 = intervals[comp1];
	  Number step2 = intervals[comp2];
	  if (comp1==comp2) 
            {
	      continue;
            }
	  Number minComp1 = Kpars[comp1]-step1;
	  Number maxComp1 = Kpars[comp1]+step1;
	  Number minComp2 = Kpars[comp2]-step2;
	  Number maxComp2 = Kpars[comp2]+step2;
	  for (int i=0; i<npoints; i++)
            {
	      Number norm1 = minComp1 + i*((maxComp1-minComp1)/npoints);
	      for (int j=0; j<npoints; j++)
                {
		  Number norm2=minComp2 + j*((maxComp2-minComp2)/npoints); 
		  K[comp1] = norm1;
		  K[comp2] = norm2;
		  Number loglike = logL(K,0,nebins);
		  ParSpace->Fill(comp1,comp2,loglike,norm1,norm2);
                }
            }
        }
    }
  cout << "PLEASE, CHECK THAT PLOTS ARE CORRECT!" << endl;  
  for (int comp1=0; comp1<1/*Nbar+1*/; comp1++)
    {
      for (int comp2=comp1; comp2<Nbar+1; comp2++)
	{
	  if (comp1==comp2)
	    {
	      continue;
	    }
	  ParSpace -> Draw("logL:comp2:comp1",
			   Form("2*(%lf-logL)<2.71  && type1==%d && type2==%d",
				maxlogL,comp1,comp2),"colz prof");
	  TF2 f("func",
		Form("0.5*([0]*(x-%lf)*(x-%lf)+[1]*(y-%lf)*(y-%lf)+2*[2]*(x-%lf)*(y-%lf)+[3])",
		     Kpars[comp1], Kpars[comp1], Kpars[comp2],
		     Kpars[comp2],Kpars[comp1],Kpars[comp2]));
	  TH1 *h_cor = ParSpace->GetHistogram();
	  h_cor -> SetMaximum(maxlogL);
	  h_cor -> Fit(&f,"Q");
	  
	  Number cfactor = f.GetParameter(2)/sqrt(fabs(f.GetParameter(0)*f.GetParameter(1)));
	  if (comp1==0) 
	    {
	      Cfactors.push_back(cfactor);
	    }
	  
	  gPad->Update();
	  gPad->Modified();
	  gPad->WaitPrimitive();
	}
    }
  cout << "Correlation Factors: ";
  for (int ii=0; ii<Nbar; ii++) 
    {
      cout << Cfactors[ii] << "  ";
    }
  cout << endl;
}

  
