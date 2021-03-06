{
  TFile *infile =  new TFile("/afs/ciemat.es/user/b/bernardos/GitHub/LMC/pipelines/pipes_in_C/ParSpace_Baryonics.root");
  TNtuple *ParSpace = (TNtuple*)infile->Get("ParSpace");
  const int Nbar = 15;
  double maxlogL=-16023.2;
  vector<double> Cfactors;
  double CMatrix[Nbar+1][Nbar+1];
  for (int i=0; i<Nbar+1; i++) CMatrix[i][i] = 1;
    

  //double Kpars[Nbar+1] ={68.9802,1.00003,0.79,1.37212,1.022,0.959367,1.008,0.988131,0.971637,1.04593,0.984007,1.02561,0.943566,1.026,0.88204,1.00782,-1.81};
  double Kpars[Nbar] = {0.99961,-0.16,4.34,3.75138,1.128,1.00341,1.20123,0.948311,2.42713,0.99101,0.994571,0.992237,0.943513,0.995791,0.83824};
  
  /*double Kpars[Nbar+1]; Kpars[0] = -0.1;
  Kpars[1] = 1;
  for (int i=2; i<Nbar+1; i++) Kpars[i] = 1; 
  Kpars[Nbar] = 0.97;*/  
  
  for (int comp1=0; comp1<6; comp1++)
    {
  
      for (int comp2=comp1; comp2<6; comp2++)
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
	  double cfactor = f.GetParameter(2)/sqrt(fabs(f.GetParameter(0)*f.GetParameter(1)));
	  if (comp1==0) 
	    {
	      Cfactors.push_back(cfactor);
	    }	    
	  
	  CMatrix[comp1][comp2] = cfactor;
	  CMatrix[comp2][comp1] = cfactor;
	  
	  gPad->Update();
	  gPad->Modified();
	  gPad->WaitPrimitive();
	}
    }
    
    for (int i=0; i<Nbar+1; i++){
      for (int j=0; j<Nbar+1; j++){
	cout << CMatrix[i][j] << "  ";
      }
      cout << endl;
    }
    
    for (int i=0; i<Nbar; i++)  cout << Cfactors[i] << " ";
    cout << endl;
    
}
