{

  const int  nmasses = 16;
  //double masses[6] = {0.100,0.200,0.500,1,5,10};
  //  double masses[nmasses] = {0.100,0.200,0.500,1,5,10,50,100};
  double masses[nmasses] = {0.100,0.200,0.300,0.400,0.500,0.600,0.800,1,4,5,8,10,40,50,80,100};
  //double masses[nmasses] = {0.100,0.200,0.300,0.400,0.500,0.600,0.800,1,4,5,8,10,40,50,80,100};
  double meanlimit[nmasses];
  double sstdevs[nmasses];
  double sstdevs2[nmasses];
  double thermalX[nmasses];
  TLegend *l =  new TLegend(0.6,0.28,0.86,0.46);

  TString particle = "b";
  for (int i=0; i<nmasses; i++){

    double dm_mass = masses[i];
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
    
    TString filename = "/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/limit"+particle+masstr+".dat";
    
    ifstream openfile(filename);
    
    double mean;
    vector<double> limits;
    int n=0;
    while(!openfile.eof()){
      double limit=0;
      openfile>>limit;
      limit=limit*3e-26;
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
    cout << mean << "  " << sqrt(sstdev)<< endl;
    meanlimit[i] = mean;
    sstdevs[i] = 1.96*sqrt(sstdev);
    sstdevs2[i]  = 1.4*sqrt(sstdev);
    thermalX[i] = 3e-26;
  }

  TCanvas *c =  new TCanvas();
  c->SetLogx();
  c->SetLogy();
  TMultiGraph *mg = new TMultiGraph();
  TGraph *g1 =  new TGraph(nmasses,masses,meanlimit);
  mg->Add(g1);
  TGraph *g1_ = new TGraph(nmasses,masses,thermalX);
  
  l->AddEntry(g1,"Upper limit","l");
  l->AddEntry(g1_,"Thermal <#sigma v>","l");
  g1_->SetLineStyle(9);
  g1_->SetLineWidth(2);
  g1_->SetLineColor(2);
  mg->Add(g1_);
  mg->Draw("a");

  TGraphErrors g(nmasses,masses,meanlimit,0,sstdevs);
  g.SetMarkerStyle(21);
  mg->SetTitle("Upper Limits on DM for "+particle+particle+" channel");
  mg->GetXaxis()->SetTitle("DM Mass (TeV)");
  mg->GetYaxis()->SetTitle("<#sigma v> (cm^{3} s^{-1})");
  g.SetFillColor(5);
  g.Draw("3 same");
  TGraphErrors g2(nmasses,masses,meanlimit,0,sstdevs2);
  g2.SetFillColor(3);
  g2.Draw("3 same");
  g.SetLineWidth(2);
  g.Draw("lX same");
  l->Draw();

  

}
