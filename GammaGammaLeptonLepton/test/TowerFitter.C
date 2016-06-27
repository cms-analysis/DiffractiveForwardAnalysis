Double_t fitfunc(Double_t *v, Double_t *par)
{
  Double_t arg = 0;
  arg = v[0];

  Double_t fitval = par[0] + TMath::Exp(-arg * par[1]);
  //  Double_t fitval = par[0] + (arg)* par[1];
  return fitval;
}

void TowerFitter()
{
  TFile *f1 = new TFile("dymumu.610.anal.root");
  TFile *f2 = new TFile("dymumu.1040.anal.root");
  TFile *f3 = new TFile("dymumu.40.anal.root");
  TFile *f4 = new TFile("gamgammumu.lpair.anal.root");
  TFile *f5 = new TFile("gamgammumu.lpairinelastic.anal.root");
  TFile *f6 = new TFile("upsilonmumu.starlight.anal.root");
  TFile *f7 = new TFile("upsilonmumu.starlight2s.anal.root");
  TFile *f8 = new TFile("upsilonmumu.starlight3s.anal.root");

  TTree *t1 = f1->Get("ntp1");
  TTree *t2 = f2->Get("ntp1");
  TTree *t3 = f3->Get("ntp1");
  TTree *t4 = f4->Get("ntp1");
  TTree *t5 = f5->Get("ntp1");
  TTree *t6 = f6->Get("ntp1");
  TTree *t7 = f7->Get("ntp1");
  TTree *t8 = f8->Get("ntp1");


  TH1F *h1 = new TH1F("h1","h1",10,0,25);
  TH1F *h2 = new TH1F("h2","h2",10,0,25);
  TH1F *h3 = new TH1F("h3","h3",10,0,25);
  TH1F *h4 = new TH1F("h4","h4",10,0,25);
  TH1F *h5 = new TH1F("h5","h5",10,0,25);
  TH1F *h6 = new TH1F("h6","h6",10,0,25);
  TH1F *h7 = new TH1F("h7","h7",10,0,25);
  TH1F *h8 = new TH1F("h8","h8",10,0,25);
  TH1F *h9 = new TH1F("h9","h9",10,0,25);

  t1->Draw("nExtraCaloTowersE5 >> h1","nTrackCand < 5","e");
  t2->Draw("nExtraCaloTowersE5 >> h2","nTrackCand < 5","e");
  t3->Draw("nExtraCaloTowersE5 >> h3","nTrackCand < 5","e");
  t4->Draw("nExtraCaloTowersE5 >> h4","nTrackCand < 5","e");
  t5->Draw("nExtraCaloTowersE5 >> h5","nTrackCand < 5","e");
  t6->Draw("nExtraCaloTowersE5 >> h6","nTrackCand < 5","e");
  t7->Draw("nExtraCaloTowersE5 >> h7","nTrackCand < 5","e");
  t8->Draw("nExtraCaloTowersE5 >> h8","nTrackCand < 5","e");

  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h4->Sumw2();
  h5->Sumw2();
  h6->Sumw2();
  h7->Sumw2();
  h8->Sumw2();
  h9->Sumw2();

  double lumi = 1.0;
  double xsec1 = 0.5 * 37820.0 * lumi / 10000.0;
  double xsec2 = 0.5 * 5951.0 * lumi / 38800.0;
  double xsec3 = 0.5 * 1797.0 * lumi / 22000.0;
  double xsec4 = 74.65 * lumi / 100000.0;
  double xsec5 = 76.2 * lumi / 20000.0;
  double xsec6 = 39.0 * lumi / 10000.0;
  double xsec7 = 13.0 * lumi / 3359.0;
  double xsec8 = 10.0 * lumi / 2636.0;

  h1->Scale(xsec1);
  h2->Scale(xsec2);
  h3->Scale(xsec3);
  h4->Scale(xsec4);
  h5->Scale(xsec5);
  h6->Scale(xsec6);
  h7->Scale(xsec7);
  h8->Scale(xsec8);

  h1->Add(h2); h1->Add(h3); h1->Add(h4); h1->Add(h5); h1->Add(h6); h1->Add(h7); h1->Add(h8); h1->Add(h9);
  h4->Add(h5); h4->Add(h6); h4->Add(h7); h4->Add(h8);
  h4->SetFillColor(4);

  h6->SetFillColor(3);
  h7->SetFillColor(6);
  h8->SetFillColor(1);
  h9->SetFillColor(2);

  h5->SetFillColor(5);

  //  h1->SetLineColor(2);
  h1->SetFillColor(0);
  h1->SetXTitle("N(Extra towers)");
  h1->SetYTitle("Events");
  h1->Draw("e");
  
  TF1 *func = new TF1("thefit",fitfunc,-10,10,2);
  func->SetParameters(0.05,-1.3);
  func->SetParNames("a","b");  
  h1->Fit("thefit","","",5,30);
  Double_t par1 = func->GetParameter("a");
  Double_t par2 = func->GetParameter("b");
  TString strpar1 = Form("%f", par1);
  TString strpar2 = Form("%f", par2);

  cout << strpar1 << ", " << strpar2 << endl;

  TF1 *func2 = new TF1("func1",strpar1 + " + TMath::Exp(-x * " + strpar2 + ")",0,5);
  func->SetLineColor(2);
  func->SetLineWidth(3);
  func2->SetLineColor(2);
  func2->SetLineStyle(2);
  func2->SetLineWidth(3);
  func2->Draw("same");
  cout << "Integral = " << func->Integral(0,5) << endl;

  h1->Draw("e");
  h4->Draw("histsame");
  h5->Draw("histsame");
  //  h6->Draw("histsame");
  //  h7->Draw("histsame");
  //  h8->Draw("histsame");

  h2->SetMarkerStyle(20);
  h2->SetLineWidth(3);
  h1->SetMarkerStyle(20);
  h1->SetLineWidth(3);

  h1->Draw("esame");
  func->Draw("same");
  func2->Draw("same");

  TLegend *l1 = new TLegend(0.6,0.6,0.8,0.8);
  l1->AddEntry(h1,"Signal plus backgrounds");
  l1->AddEntry(h4,"Elastic signal (#gamma #gamma and #Upsilon)");
  l1->AddEntry(h5,"Singly inelastic #gamma #gamma #rightarrow #mu^{+} #mu^{-}");
  l1->AddEntry(func,"Fit");
  l1->Draw("same");

}

void ETowerFitter()
{
  TFile *f1 = new TFile("dyee.610.anal.root");
  TFile *f2 = new TFile("dyee.1040.anal.root");
  TFile *f3 = new TFile("dyee.40.anal.root");
  TFile *f4 = new TFile("gamgamee.lpair.anal.root");
  TFile *f5 = new TFile("gamgamee.lpairinelastic.anal.root");

  TTree *t1 = f1->Get("ntp1");
  TTree *t2 = f2->Get("ntp1");
  TTree *t3 = f3->Get("ntp1");
  TTree *t4 = f4->Get("ntp1");
  TTree *t5 = f5->Get("ntp1");


  TH1F *h1 = new TH1F("h1","h1",10,0,25);
  TH1F *h2 = new TH1F("h2","h2",10,0,25);
  TH1F *h3 = new TH1F("h3","h3",10,0,25);
  TH1F *h4 = new TH1F("h4","h4",10,0,25);
  TH1F *h5 = new TH1F("h5","h5",10,0,25);

  t1->Draw("nExtraCaloTowersE5 >> h1","nTrackCand < 5","e");
  t2->Draw("nExtraCaloTowersE5 >> h2","nTrackCand < 5","e");
  t3->Draw("nExtraCaloTowersE5 >> h3","nTrackCand < 5","e");
  t4->Draw("nExtraCaloTowersE5 >> h4","nTrackCand < 5","e");
  t5->Draw("nExtraCaloTowersE5 >> h5","nTrackCand < 5","e");

  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h4->Sumw2();
  h5->Sumw2();

  double lumi = 1.0;
  double xsec1 = 0.5 * 37820.0 * lumi / 10000.0;
  double xsec2 = 0.5 * 5951.0 * lumi / 38800.0;
  double xsec3 = 0.5 * 1797.0 * lumi / 22000.0;
  double xsec4 = 10.35 * lumi / 100000.0;
  double xsec5 = 13.60 * lumi / 18000.0;

  h1->Scale(xsec1);
  h2->Scale(xsec2);
  h3->Scale(xsec3);
  h4->Scale(xsec4);
  h5->Scale(xsec5);

  h1->Add(h2); h1->Add(h3); h1->Add(h4); h1->Add(h5); 
  h4->Add(h5); 
  h4->SetFillColor(4);


  h5->SetFillColor(5);

  //  h1->SetLineColor(2);
  h1->SetFillColor(0);
  h1->SetXTitle("N(Extra towers)");
  h1->SetYTitle("Events");
  h1->Draw("e");
  
  TF1 *func = new TF1("thefit",fitfunc,-10,10,2);
  func->SetParameters(0.05,-1.3);
  func->SetParNames("a","b");  
  h1->Fit("thefit","","",5,25);
  Double_t par1 = func->GetParameter("a");
  Double_t par2 = func->GetParameter("b");
  TString strpar1 = Form("%f", par1);
  TString strpar2 = Form("%f", par2);

  cout << strpar1 << ", " << strpar2 << endl;

  TF1 *func2 = new TF1("func1",strpar1 + " + TMath::Exp(-x * " + strpar2 + ")",0,5);
  func->SetLineColor(2);
  func->SetLineWidth(3);
  func2->SetLineColor(2);
  func2->SetLineStyle(2);
  func2->SetLineWidth(3);
  func2->Draw("same");
  cout << "Integral = " << func->Integral(0,5) << endl;

  h1->Draw("e");
  h4->Draw("histsame");
  h5->Draw("histsame");
  //  h6->Draw("histsame");
  //  h7->Draw("histsame");
  //  h8->Draw("histsame");

  h2->SetMarkerStyle(20);
  h2->SetLineWidth(3);
  h1->SetMarkerStyle(20);
  h1->SetLineWidth(3);

  h1->Draw("esame");
  func->Draw("same");
  func2->Draw("same");

  TLegend *l1 = new TLegend(0.6,0.6,0.8,0.8);
  l1->AddEntry(h1,"Signal plus backgrounds");
  l1->AddEntry(h4,"Elastic signal (#gamma #gamma and #Upsilon)");
  l1->AddEntry(h5,"Singly inelastic #gamma #gamma #rightarrow e^{+} e^{-}");
  l1->AddEntry(func,"Fit");
  l1->Draw("same");
}
