{
  int nfiles = 12;
  TString PATH = "/home/queenmab/GitHub/LMC/results/";
  TString files[nfiles] = {"limit10TeV.dat","limit10TeV_DiffuseTruncated.dat","limitW10TeVBestCase_DiffuseTruncated_tol_0-01.dat","limitW10TeVBestCase_tol_0-01.dat","limitW10TeVBestCase_DiffuseTruncated.dat","limitW10TeVBestCase.dat","limit10TeVDiffuseTruncated1TeV.dat","limitW10TeVBestCase_DiffuseTruncated1TeV_tol0-01.dat","limitW10TeVBestCase_DiffuseTruncated1TeV.dat","limit10TeV_1PS_DiffuseTruncated500GeV.dat","limitW10TeVBestCase_1PSDiffuseTruncated500GeV.dat","limitW10TeVBestCase_1PSDiffuseTruncated500GeV_tol0-01.dat"};

  for (int i=0; i<nfiles; i++){
    TString filename = PATH+files[i];
    ifstream openfile(filename);
      
    double mean;
    vector<double> limits;
    int n=0;
    while(!openfile.eof()){
      double limit=0;
      openfile>>limit;
      limit=limit;//*3e-26;
      limits.push_back(limit);
      mean+=limit;
      n++;

    }
    openfile.close();
    mean/=n;
    
    double sstdev=0;
    for (int ii=0; ii<n; ii++){
      double dif = limits[ii]-mean;
      sstdev+=dif*dif;
    }
    
    sstdev/=n;
    cout << files[i] << "  " << mean << "  " << sqrt(sstdev)<< endl;
    
  }
  
}
