{
  

  TMultiGraph *mg =  new TMultiGraph();
  
  TLegend *l =  new TLegend(0.1,0.7,0.48,0.9);

  char *comp[7] = {"Irf","Leptonic","Hadronic","E1","E2","E3","PS1"};

  TString files[4] = {"corrFactors_1PS.dat","corrFactors_1PS_TruncDiff10TeV.dat","corrFactors_1PS_TruncDiff1TeV.dat","corrFactors_1PS_TruncDiff500GeV.dat"};

  
  for (int ii=0;ii<4; ii++){
    TNtuple *nt =  new TNtuple("nt","nt","e:e1:e2:e3:e4:e5:e6:e7");
    nt->ReadFile("../../results/"+files[ii]);
    for (int i=2; i<=3; i++){
      nt->Draw(Form("e%d:e",i));
      TGraph *g =  new TGraph(5,nt->GetV2(),nt->GetV1());
      const char *name = comp[i-1];
      g->SetName(name);
      g->SetMarkerColor(ii+1);
      g->SetLineColor(ii+1);
      g->SetLineWidth(2);
      if (i==3) g->SetLineStyle(7);
      g->SetMarkerStyle(7);
      l->AddEntry(g,name,"l");
      mg->Add(g,"l");
    }
  }
  mg->Draw("a");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0.100,100);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(-1,1);
  mg->GetHistogram()->GetXaxis()->SetTitle("Mass[TeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Correlation Factor");
  gPad->SetLogx();
  l->Draw();
}
