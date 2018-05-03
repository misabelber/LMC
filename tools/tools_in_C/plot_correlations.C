{
  TNtuple *nt =  new TNtuple("nt","nt","e:e1:e2:e3:e4:e5:e6");
  nt->ReadFile("corrpermassKSP_v2.dat");
  TMultiGraph *mg =  new TMultiGraph();
  
  TLegend *l =  new TLegend(0.1,0.7,0.48,0.9);

  char *comp[6] = {"Irf","Leptonic","Hadronic","E1","E2","E3"};

  cout << nt->GetEntries() << endl;
  for (int i=1; i<=6; i++){
    nt->Draw(Form("e%d:e",i));
    TGraph *g =  new TGraph(8,nt->GetV2(),nt->GetV1());
    const char *name = comp[i-1];
    g->SetName(name);
    g->SetMarkerColor(i);
    g->SetLineColor(i);
    g->SetLineWidth(2);
    g->SetMarkerStyle(7);
    l->AddEntry(g,name,"l");
    mg->Add(g,"l");
  }
  mg->Draw("a");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0.100,100);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(-1,1);
  mg->GetHistogram()->GetXaxis()->SetTitle("Mass[TeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Correlation Factor");
  gPad->SetLogx();
  l->Draw();
}
