void MakeSkimQAPlots(TString dataset = "")
{
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineWidth(3);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleXOffset(0.7);
  gStyle->SetTitleYOffset(1.05);


  TString fname1("mumu.recolevel." + dataset + ".root");
  TString fname2("ee.recolevel." + dataset + ".root");

  //  TFile *f1 = new TFile(fname1);
  //  TFile *f2 = new TFile(fname2);

  //  TTree *t1 = f1->Get("ntp1");
  //  TTree *t2 = f2->Get("ntp1");

  TChain *t1 = new TChain("ntp1");
  t1->Add("crab_0_071121_090253/res/mumu.recolevel_*.root");
  TChain *t2 = new TChain("ntp1");
  t2->Add("crab_0_071121_102409/res/ee.recolevel_*.root");

  TH1F *hpsim = new TH1F("hpsim","hpsim",50,2.5,4.5);
  TH1F *hpsie = new TH1F("hpsie","hpsie",50,2.5,4.5);
  TH1F *hupsm = new TH1F("hupsm","hupsm",50,8,12);
  TH1F *hupse = new TH1F("hupse","hupse",50,8,12);
  TH1F *hzm = new TH1F("hzm","hzm",50,60,120);
  TH1F *hze = new TH1F("hze","hze",50,60,120);

  TCanvas *c1 = new TCanvas("c1","#mu^{+}#mu^{-}",1400,800);
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
  htemp->SetXTitle("p_{T} (#mu)");
  htemp->Draw();
  c1->cd(5);
  t1->Draw("MuonCand_pt[0]-MuonCand_pt[1]");
  htemp->SetXTitle("#Delta p_{T}");
  htemp->Draw();
  c1->cd(6);
  t1->Draw("MuMu_dphi");
  htemp->SetXTitle("#Delta #phi");
  htemp->Draw();
  c1->cd(7);
  t1->Draw("CaloTower_eta");
  htemp->SetXTitle("#eta (towers)");
  htemp->Draw();
  c1->cd(8);
  c1->SetLogy();
  t1->Draw("CaloTower_e");
  htemp->SetXTitle("E (towers)");
  htemp->Draw();
  c1->SetLogy(0);
  c1->cd(9);
  t1->Draw("nExtraCaloTowersE5");
  htemp->SetXTitle("Extra towers (E > 5 GeV)");
  htemp->Draw();
  TString muname = "dimuons.skimqa." + dataset + ".png";
  c1->SaveAs(muname);

  TCanvas *c2 = new TCanvas("c2","e^{+}e^{-}",1400,800);
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
  htemp->SetXTitle("E_{T} (e)");
  htemp->Draw();
  c2->cd(5);
  t2->Draw("EleCand_et[0]-EleCand_et[1]");
  htemp->SetXTitle("#Delta E_{T}");
  htemp->Draw();
  c2->cd(6);
  t2->Draw("ElEl_dphi");
  htemp->SetXTitle("#Delta #phi");
  htemp->Draw();
  c2->cd(7);
  t2->Draw("CaloTower_eta");
  htemp->SetXTitle("#eta (towers)");
  htemp->Draw();
  c2->cd(8);
  c2->SetLogy();
  t2->Draw("CaloTower_e");
  htemp->SetXTitle("E (towers)");
  htemp->Draw();
  c2->SetLogy(0);
  c2->cd(9);
  t2->Draw("nExtraCaloTowersE5");
  htemp->SetXTitle("Extra towers (E > 5 GeV)");
  htemp->Draw();
  TString elname = "dielectrons.skimqa." + dataset + ".png";
  c2->SaveAs(elname);


}
