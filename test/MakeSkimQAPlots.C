void MakeSkimQAPlots()
{
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineWidth(3);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);

  TFile *f1 = new TFile("mumu.recolevel.root");
  TFile *f2 = new TFile("ee.recolevel.root");

  TTree *t1 = f1->Get("ntp1");
  TTree *t2 = f2->Get("ntp1");

  TH1F *hpsim = new TH1F("hpsim","hpsim",50,2.5,4.5);
  TH1F *hpsie = new TH1F("hpsie","hpsie",50,2.5,4.5);
  TH1F *hupsm = new TH1F("hupsm","hupsm",50,8,12);
  TH1F *hupse = new TH1F("hupse","hupse",50,8,12);
  TH1F *hzm = new TH1F("hzm","hzm",50,60,120);
  TH1F *hze = new TH1F("hze","hze",50,60,120);

  TCanvas *c1 = new TCanvas("#mu^{+}#mu^{-}","#mu^{+}#mu^{-}");
  c1->Divide(3,3);
  c1->cd(1);
  hpsim->SetXTitle("#psi mass");
  t1->Draw("MuMu_mass >> hpsim");
  c1->cd(2);
  hupsm->SetXTitle("#Upsilon mass");
  t1->Draw("MuMu_mass >> hupsm");
  c1->cd(3);
  hzm->SetXTitle("Z mass");
  t1->Draw("MuMu_mass >> hzm");
  c1->cd(4);
  t1->Draw("MuonCand_pt");
  c1->cd(5);
  t1->Draw("MuonCand_pt[0]-MuonCand_pt[1]");
  c1->cd(6);
  t1->Draw("MuMu_dphi");
  c1->cd(7);
  t1->Draw("CaloTower_eta");
  c1->cd(8);
  c1->SetLogy();
  t1->Draw("CaloTower_e");
  c1->SetLogy(0);

  TCanvas *c2 = new TCanvas("e^{+}e^{-}","e^{+}e^{-}");
  c2->Divide(3,3);
  c2->cd(1);
  hpsie->SetXTitle("#psi mass");
  t2->Draw("ElEl_mass >> hpsie");
  c2->cd(2);
  hupse->SetXTitle("#Upsilon mass");
  t2->Draw("ElEl_mass >> hupse");
  c2->cd(3);
  hze->SetXTitle("Z mass");
  t2->Draw("ElEl_mass >> hze");
  c2->cd(4);
  t2->Draw("EleCand_et");
  c2->cd(5);
  t2->Draw("EleCand_et[0]-EleCand_et[1]");
  c2->cd(6);
  t2->Draw("ElEl_dphi");
  c2->cd(7);
  t2->Draw("CaloTower_eta");
  c2->cd(8);
  c2->SetLogy();
  t2->Draw("CaloTower_e");
  c2->SetLogy(0);
}
