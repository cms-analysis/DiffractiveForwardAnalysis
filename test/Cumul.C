void Cumul()
{
gROOT->Reset();

gStyle->SetPalette(1);

#define pi 3.14159265359
//definition des fichiers + Tree
  TFile *f0 = new TFile("MuPrompt-361p4-138560_140076.root"); // 
  TTree *t0 = f0->Get("ntp1");
  TFile *f1 = new TFile("../MuMu_0PU.root"); //
  TTree *t1 = f1->Get("ntp1");
  TFile *f2 = new TFile("../InelEl_MuMu.root"); //
  TTree *t2 = f2->Get("ntp1");
  TFile *f3 = new TFile("../InelInel_MuMu.root"); //
  TTree *t3 = f3->Get("ntp1");
  TFile *f4 = new TFile("../Upsilon_MuMu.root"); //
  TTree *t4 = f4->Get("ntp1");
  
// definitions : example for mass:
//                /\
//               |  |  /\    /\  <<----Upsilon resonnances (4)
//               |  | |  |  |  |  
//         |--------------------|
//         |    El-El (3)       |
//         |--------------------|
//         |    Inel-El (2)     |
//         |--------------------|
//         |    Inel-Inel (1)   |
//         |--------------------|
//					+ Data (0)

  TH1F* hEB0 = new TH1F("neb_data","",51,-1.,50.);
  TH1F* hEE0 = new TH1F("nee_data","",51,-1.,50.);
  TH1F* hHB0 = new TH1F("nhb_data","",51,-1.,50.);
  TH1F* hHE0 = new TH1F("nhe_data","",51,-1.,50.);
  TH1F* hHFp0 = new TH1F("nhfp_data","",51,-1.,50.);
  TH1F* hHFm0 = new TH1F("ncastor_data","",31,-1.,30.);
  TH1F* nTower0 = new TH1F("nTower_data","",31,-1.,30.);
  TH1F* hEB1 = new TH1F("neb_cumulInelInel","",51,-1.,50.);
  TH1F* hEE1 = new TH1F("nee_cumulInelInel","",51,-1.,50.);
  TH1F* hHB1 = new TH1F("nhb_cumulInelInel","",51,-1.,50.);
  TH1F* hHE1 = new TH1F("nhe_cumulInelInel","",51,-1.,50.);
  TH1F* hHFp1 = new TH1F("nhfp_cumulInelInel","",51,-1.,50.);
  TH1F* hHFm1 = new TH1F("ncastor_cumulInelInel","",31,-1.,30.);
  TH1F* nTower1 = new TH1F("nTower_cumulInelInel","",31,-1.,30.);
  TH1F* hEB2 = new TH1F("neb_cumulInelEl","",51,-1.,50.);
  TH1F* hEE2 = new TH1F("nee_cumulInelEl","",51,-1.,50.);
  TH1F* hHB2 = new TH1F("nhb_cumulInelEl","",51,-1.,50.);
  TH1F* hHE2 = new TH1F("nhe_cumulInelEl","",51,-1.,50.);
  TH1F* hHFp2 = new TH1F("nhfp_cumulInelEl","",51,-1.,50.);
  TH1F* hHFm2 = new TH1F("ncastor_cumulInelEl","",31,-1.,30.);
  TH1F* nTower2 = new TH1F("nTower_cumulInelEl","",31,-1.,30.);
  TH1F* hEB3 = new TH1F("neb_cumulElEl","",51,-1.,50.);
  TH1F* hEE3 = new TH1F("nee_cumulElEl","",51,-1.,50.);
  TH1F* hHB3 = new TH1F("nhb_cumulElEl","",51,-1.,50.);
  TH1F* hHE3 = new TH1F("nhe_cumulElEl","",51,-1.,50.);
  TH1F* hHFp3 = new TH1F("nhfp_cumulElEl","",51,-1.,50.);
  TH1F* hHFm3 = new TH1F("ncastor_cumulElEl","",31,-1.,30.);
  TH1F* nTower3 = new TH1F("nTower_cumulElEl","",31,-1.,30.);
  TH1F* hEB4 = new TH1F("neb_cumulUps","",51,-1.,50.);
  TH1F* hEE4 = new TH1F("nee_cumulUps","",51,-1.,50.);
  TH1F* hHB4 = new TH1F("nhb_cumulUps","",51,-1.,50.);
  TH1F* hHE4 = new TH1F("nhe_cumulUps","",51,-1.,50.);
  TH1F* hHFp4 = new TH1F("nhfp_cumulUps","",51,-1.,50.);
  TH1F* hHFm4 = new TH1F("ncastor_cumulUps","",31,-1.,30.);
  TH1F* nTower4 = new TH1F("nTower_cumulUps","",31,-1.,30.);

  TH1F* nTrack0 = new TH1F("ntrack_data","",101,-1.,100.);
  TH1F* nTrack1 = new TH1F("ntrack_cumulInelInel","",101,-1.,100.);
  TH1F* nTrack2 = new TH1F("ntrack_cumulInelEl","",101,-1.,100.);
  TH1F* nTrack3 = new TH1F("ntrack_cumulElEl","",101,-1.,100.);
  TH1F* nTrack4 = new TH1F("ntrack_cumulUps","",101,-1.,100.);

  TH1F* MuMuMass0 = new TH1F("mass_data","",120.,0.,40.);
  TH1F* MuMuMass1 = new TH1F("mass_cumulInelInel","",120.,0.,40.);
  TH1F* MuMuMass2 = new TH1F("mass_cumulInelEl","",120.,0.,40.);
  TH1F* MuMuMass3 = new TH1F("mass_cumulElEl","",120.,0.,40.);
  TH1F* MuMuMass4 = new TH1F("mass_cumulUps","",120.,0.,40.);

  TH1F* MuMuMassUps0 = new TH1F("massUps_data","",120.,8.,12.);
  TH1F* MuMuMassUps1 = new TH1F("massUps_cumulInelInel","",120.,8.,12.);
  TH1F* MuMuMassUps2 = new TH1F("massUps_cumulInelEl","",120.,8.,12.);
  TH1F* MuMuMassUps3 = new TH1F("massUps_cumulElEl","",120.,8.,12.);
  TH1F* MuMuMassUps4 = new TH1F("massUps_cumulUps","",120.,8.,12.);


  TH1F* MuMudpt0 = new TH1F("dpt_data","",110.,-1.,10.);
  TH1F* MuMudpt1 = new TH1F("dpt_cumulInelInel","",110.,-1.,10.);
  TH1F* MuMudpt2 = new TH1F("dpt_cumulInelEl","",110.,-1.,10.);
  TH1F* MuMudpt3 = new TH1F("dpt_cumulElEl","",110.,-1.,10.);
  TH1F* MuMudpt4 = new TH1F("dpt_cumulUps","",110.,-1.,10.);

  TH1F* MuMudphi0 = new TH1F("dphi_data","",60.,-0.1,1.1);
  TH1F* MuMudphi1 = new TH1F("dphi_cumulInelInel","",60.,-0.1,1.1);
  TH1F* MuMudphi2 = new TH1F("dphi_cumulInelEl","",60.,-0.1,1.1);
  TH1F* MuMudphi3 = new TH1F("dphi_cumulElEl","",60.,-0.1,1.1);
  TH1F* MuMudphi4 = new TH1F("dphi_cumulUps","",60.,-0.1,1.1);

  TH1F* MuMudeta0 = new TH1F("deta_data","",28.,-1.,6.);
  TH1F* MuMudeta1 = new TH1F("deta_cumulInelInel","",28.,-1.,6.);
  TH1F* MuMudeta2 = new TH1F("deta_cumulInelEl","",28.,-1.,6.);
  TH1F* MuMudeta3 = new TH1F("deta_cumulElEl","",28.,-1.,6.);
  TH1F* MuMudeta4 = new TH1F("deta_cumulUps","",28.,-1.,6.);

  TH1F* ZDCemplus0 = new TH1F("zdcEm+_data","",72.,-100.,3500.);
  TH1F* ZDCemplus1 = new TH1F("zdcEm+_cumulInelInel","",72.,-100.,3500.);
  TH1F* ZDCemplus2 = new TH1F("zdcEm+_cumulInelEl","",72.,-100.,3500.);
  TH1F* ZDCemplus3 = new TH1F("zdcEm+_cumulElEl","",72.,-100.,3500.);
  TH1F* ZDCemplus4 = new TH1F("zdcEm+_cumulUps","",72.,-100.,3500.);
  TH1F* ZDChadplus0 = new TH1F("zdcHad+_data","",92.,-1000.,45000.);
  TH1F* ZDChadplus1 = new TH1F("zdcHad+_cumulInelInel","",92.,-1000.,45000.);
  TH1F* ZDChadplus2 = new TH1F("zdcHad+_cumulInelEl","",92.,-1000.,45000.);
  TH1F* ZDChadplus3 = new TH1F("zdcHad+_cumulElEl","",92.,-1000.,45000.);
  TH1F* ZDChadplus4 = new TH1F("zdcHad+_cumulUps","",92.,-1000.,45000.);
  TH1F* ZDCemminus0 = new TH1F("zdcEm-_data","",72.,-100.,3500.);
  TH1F* ZDCemminus1 = new TH1F("zdcEm-_cumulInelInel","",72.,-100.,3500.);
  TH1F* ZDCemminus2 = new TH1F("zdcEm-_cumulInelEl","",72.,-100.,3500.);
  TH1F* ZDCemminus3 = new TH1F("zdcEm-_cumulElEl","",72.,-100.,3500.);
  TH1F* ZDCemminus4 = new TH1F("zdcEm-_cumulUps","",72.,-100.,3500.);
  TH1F* ZDChadminus0 = new TH1F("zdcHad-_data","",92.,-1000.,45000.);
  TH1F* ZDChadminus1 = new TH1F("zdcHad-_cumulInelInel","",92.,-1000.,45000.);
  TH1F* ZDChadminus2 = new TH1F("zdcHad-_cumulInelEl","",92.,-1000.,45000.);
  TH1F* ZDChadminus3 = new TH1F("zdcHad-_cumulElEl","",92.,-1000.,45000.);
  TH1F* ZDChadminus4 = new TH1F("zdcHad-_cumulUps","",92.,-1000.,45000.);
  TH1F* ZDCtime0 =  new TH1F("zdcTime_data","",80.,-20.,60.);
  TH1F* ZDCtime1 =  new TH1F("zdcTime_cumulInelInel","",80.,-20.,60.);
  TH1F* ZDCtime2 =  new TH1F("zdcTime_cumulInelEl","",80.,-20.,60.);
  TH1F* ZDCtime3 =  new TH1F("zdcTime_cumulElEl","",80.,-20.,60.);
  TH1F* ZDCtime4 =  new TH1F("zdcTime_cumulUps","",80.,-20.,60.);
  TH1F* ZDCenergyEM0 =  new TH1F("zdcEnEM_data","",100.,-500.,2500.);
  TH1F* ZDCenergyEM1 =  new TH1F("zdcEnEM_cumulInelInel","",100.,-500.,2500.);
  TH1F* ZDCenergyEM2 =  new TH1F("zdcEnEM_cumulInelEl","",100.,-500.,2500.);
  TH1F* ZDCenergyEM3 =  new TH1F("zdcEnEM_cumulElEl","",100.,-500.,2500.);
  TH1F* ZDCenergyEM4 =  new TH1F("zdcEnEM_cumulUps","",100.,-500.,2500.);
  TH1F* ZDCenergyHAD0 =  new TH1F("zdcEnHAD_data","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD1 =  new TH1F("zdcEnHAD_cumulInelInel","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD2 =  new TH1F("zdcEnHAD_cumulInelEl","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD3 =  new TH1F("zdcEnHAD_cumulElEl","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD4 =  new TH1F("zdcEnHAD_cumulUps","",100.,-5000.,25000.);

  TH1F* CastorSumE0 = new TH1F("castorE_data","",140.,-500.,10300.);
  TH1F* CastorSumE1 = new TH1F("castorE_cumulInelInel","",140.,-500.,10300.);
  TH1F* CastorSumE2 = new TH1F("castorE_cumulInelEl","",140.,-500.,10300.);
  TH1F* CastorSumE3 = new TH1F("castorE_cumulElEl","",140.,-500.,10300.);
  TH1F* CastorSumE4 = new TH1F("castorE_cumulUps","",140.,-500.,10300.);

// definitions des # d'entrées
  const int NUM0 = t0->GetEntries();
  const int NUM1 = t1->GetEntries();
  const int NUM2 = t2->GetEntries();
  const int NUM3 = t3->GetEntries();
  const int NUM4 = t4->GetEntries();

  const float integrated_lumi = 120.814826; //in nb-1

  const float fac_lumi1 = 1.0820e-6*integrated_lumi;
  const float fac_lumi2 = 3.05250e-6*integrated_lumi;
  const float fac_lumi3 = 4.740e-6*integrated_lumi;
  const float fac_lumi4 = 1.350e-6*integrated_lumi;
  const float fac_lumi0 = 1.0;

//definition des variables
  Int_t hlt_d1[1], hlt_d2[1], hlt_d0[1], hlt_d3[1], hlt_d4[1];
  Int_t techBit1[1][128], techBit2[1][128], techBit0[1][128], techBit3[1][128], techBit4[1][128];
// MuonID
  Int_t var_idA1[10], var_idA2[10], var_idA0[10],var_idA3[10],var_idA4[10],var_idB1[10], var_idB2[10], var_idB0[10],var_idB3[10],var_idB4[10],var_idC1[10], var_idC2[10],var_idC0[10],var_idC3[10],var_idC4[10], var_idD1[10], var_idD2[10], var_idD0[10], var_idD3[10],var_idD4[10],var_idE1[10], var_idE2[10], var_idE0[10] , var_idE3[10], var_idE4[10];
  Int_t var_nMuon[1];
// RecoTrack
  Int_t var_nTrack1[1], var_nTrack2[1], var_nTrack0[1],var_nTrack3[1],var_nTrack4[1];
  Int_t var_nTrackQual1[1], var_nTrackQual2[1], var_nTrackQual0[1], var_nTrackQual3[1], var_nTrackQual4[1];
  Double_t var_TrackPt1[2000], var_TrackPt2[2000], var_TrackPt0[2000], var_TrackPt3[2000], var_TrackPt4[2000];
  Double_t var_TrackD1[2000],var_TrackD2[2000], var_TrackD0[2000], var_TrackD3[2000], var_TrackD4[2000];
  Double_t var_TrackQuality1[2000],var_TrackQuality2[2000], var_TrackQuality0[2000], var_TrackQuality3[2000], var_TrackQuality4[2000];
// Prim Vtx
  Int_t var_nvtx1[1],var_nvtx2[1], var_nvtx0[1], var_nvtx3[1], var_nvtx4[1];
  Int_t var_vtxTrack1[10],var_vtxTrack2[10], var_vtxTrack0[10],  var_vtxTrack3[10],  var_vtxTrack4[10];
  Double_t var_vtxZ1[10],var_vtxZ2[10], var_vtxZ0[10], var_vtxZ3[10], var_vtxZ4[10];
  Double_t var_vtxX1[10],var_vtxX2[10], var_vtxX0[10], var_vtxX3[10], var_vtxX4[10];
  Double_t var_vtxY1[10],var_vtxY2[10], var_vtxY0[10], var_vtxY3[10], var_vtxY4[10];
  Double_t var_MuMuvtxX1[10],var_MuMuvtxX2[10], var_MuMuvtxX0[10], var_MuMuvtxX3[10], var_MuMuvtxX4[10];
  Double_t var_MuMuvtxY1[10],var_MuMuvtxY2[10], var_MuMuvtxY0[10], var_MuMuvtxY3[10], var_MuMuvtxY4[10];;
  Double_t var_MuMuvtxZ1[10],var_MuMuvtxZ2[10], var_MuMuvtxZ0[10], var_MuMuvtxZ3[10], var_MuMuvtxZ4[10];
  Int_t var_MuMuvtxValid1[10],var_MuMuvtxValid2[10], var_MuMuvtxValid0[10], var_MuMuvtxValid3[10], var_MuMuvtxValid4[10];
  Double_t var_vertexChi2_1[10],var_vertexChi2_2[10], var_vertexChi2_0[10], var_vertexChi2_3[10], var_vertexChi2_4[10];
  Double_t var_vertexNdf1[10], var_vertexNdf2[10], var_vertexNdf0[10],  var_vertexNdf3[10], var_vertexNdf4[10];

// CaloTowers 
  Int_t var_ncalo1[1], var_ncalo2[1], var_ncalo0[1], var_ncalo3[1], var_ncalo4[1];
  Int_t var_tower1[1], var_tower2[1], var_tower0[1], var_tower3[1], var_tower4[1];
  Int_t var_caloId1[2000],var_caloId2[2000], var_caloId0[2000], var_caloId3[2000], var_caloId4[2000];
  Double_t var_caloEn1[2000],var_caloEn2[2000], var_caloEn0[2000],  var_caloEn3[2000],  var_caloEn4[2000];
  Double_t var_caloTime1[2000], var_caloTime2[2000], var_caloTime0[2000], var_caloTime3[2000], var_caloTime4[2000];
  Double_t var_etmiss1[1], var_etmiss2[1], var_etmiss0[1], var_etmiss3[1], var_etmiss4[1];
  Double_t var_calodR1[2000], var_calodR2[2000], var_calodR0[2000], var_calodR3[2000], var_calodR4[2000];
  Double_t var_caloZ1[2000], var_caloZ2[2000], var_caloZ0[2000], var_caloZ3[2000], var_caloZ4[2000];

// Bunch crossing
  Int_t var_bx0[1], var_run0[1], var_ls0[1], var_event0[1], var_event3[1], var_event4[1];
// MuMu kinematics
  Double_t var_mass1[5], var_mass2[5] ,var_mass0[2], var_mass3[2], var_mass4[2];
  Double_t var_dpt1[5], var_dpt2[5], var_dpt0[5], var_dpt3[5], var_dpt4[5];
  Double_t var_dphi1[5], var_dphi2[5], var_dphi0[5], var_dphi3[5], var_dphi4[5];
  Int_t var_global1[10], var_global2[10], var_global0[10], var_global3[10], var_global4[10],var_tracker1[10], var_tracker2[10],var_tracker0[10],var_tracker3[10],var_tracker4[10],var_standalone1[10], var_standalone2[10], var_standalone0[10], var_standalone3[10],var_standalone4[10];
  Double_t var_pt1[10], var_pt2[10], var_pt0[10],var_pt3[10],var_pt4[10],var_pz1[10], var_pz2[10], var_pz0[10],var_pz3[10],var_pz4[10],var_phi1[10], var_phi2[10], var_phi0[10],var_phi3[10], var_phi4[10],var_eta1[10], var_eta2[10], var_eta0[10], var_eta3[10], var_eta4[10];
  Int_t var_nhitsTrack1[10], var_nhitsTrack2[10], var_nhitsTrack0[10],  var_nhitsTrack3[10], var_nhitsTrack4[10];
  Int_t var_Pair0[2], var_Pair1[2],var_Pair2[2],var_Pair3[2],var_Pair4[2];

// ZDC
  Int_t var_nZDC1[1], var_nZDC2[1], var_nZDC0[1], var_nZDC3[1], var_nZDC4[1];
  Int_t var_zdcsection1[5000],var_zdcsection2[5000],var_zdcsection0[5000],var_zdcsection3[5000],var_zdcsection4[5000];
  Double_t var_zdcE1[5000], var_zdcE2[5000], var_zdcE0[5000], var_zdcE3[5000], var_zdcE4[5000];
  Double_t var_zdcEmMinus1[1], var_zdcEmMinus2[1],var_zdcEmMinus0[1],var_zdcEmMinus3[1], var_zdcEmMinus4[1],var_zdcHadMinus1[1], var_zdcHadMinus2[1], var_zdcHadMinus0[1], var_zdcHadMinus3[1], var_zdcHadMinus4[1];
  Double_t var_zdcEmPlus1[1], var_zdcEmPlus2[1], var_zdcEmPlus0[1], var_zdcEmPlus3[1],var_zdcEmPlus4[1],var_zdcHadPlus1[1], var_zdcHadPlus2[1], var_zdcHadPlus0[1], var_zdcHadPlus3[1], var_zdcHadPlus4[1];
  Double_t var_zdcTime1[5000], var_zdcTime2[5000], var_zdcTime0[5000], var_zdcTime3[5000], var_zdcTime4[5000];

//Castor
  Int_t var_nCastor1[1],var_nCastor2[1], var_nCastor0[1], var_nCastor3[1], var_nCastor4[1];
  Double_t var_CastorE1[1000],  var_CastorE2[1000], var_CastorE0[1000], var_CastorE3[1000], var_CastorE4[1000];
  Double_t var_CastorEta1[1000],  var_CastorEta2[1000], var_CastorEta0[1000], var_CastorEta3[1000], var_CastorEta4[1000];
  Double_t var_CastorPhi1[1000],  var_CastorPhi2[1000], var_CastorPhi0[1000], var_CastorPhi3[1000], var_CastorPhi4[1000];
  Double_t var_CastorRecHit1[1], var_CastorRecHit2[1], var_CastorRecHit0[1], var_CastorRecHit3[1], var_CastorRecHit4[1];

  t1->SetBranchAddress("HLT_DoubleMu0",hlt_d1);
  t2->SetBranchAddress("HLT_DoubleMu0",hlt_d2);
  t0->SetBranchAddress("HLT_DoubleMu0",hlt_d0);
  t3->SetBranchAddress("HLT_DoubleMu0",hlt_d3);
  t4->SetBranchAddress("HLT_DoubleMu0",hlt_d4);

  t0->SetBranchAddress("L1TechnicalTriggers",techBit0);
  t1->SetBranchAddress("L1TechnicalTriggers",techBit1);
  t2->SetBranchAddress("L1TechnicalTriggers",techBit2);
  t3->SetBranchAddress("L1TechnicalTriggers",techBit3);
  t4->SetBranchAddress("L1TechnicalTriggers",techBit4);

  t1->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA1);
  t2->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA2);
  t1->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB1);
  t2->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB2);
  t1->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC1);
  t2->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC2);
  t1->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD1);
  t2->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD2);
  t1->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE1);
  t2->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE2);
  t0->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA0);
  t0->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB0);
  t0->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC0);
  t0->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD0);
  t0->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE0);
  t3->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA3);
  t3->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB3);
  t3->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC3);
  t3->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD3);
  t3->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE3);
  t4->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA4);
  t4->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB4);
  t4->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC4);
  t4->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD4);
  t4->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE4);

  t1->SetBranchAddress("nTrackCand",var_nTrack1);
  t1->SetBranchAddress("nQualityTrackCand",var_nTrackQual1);
  t1->SetBranchAddress("TrackCand_pt",var_TrackPt1);
  t1->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD1);
  t1->SetBranchAddress("TrackCand_purity",var_TrackQuality1);
  t0->SetBranchAddress("nTrackCand",var_nTrack0);
  t0->SetBranchAddress("nQualityTrackCand",var_nTrackQual0);
  t0->SetBranchAddress("TrackCand_pt",var_TrackPt0);
  t0->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD0);
  t0->SetBranchAddress("TrackCand_purity",var_TrackQuality0);
  t2->SetBranchAddress("nTrackCand",var_nTrack2);
  t2->SetBranchAddress("nQualityTrackCand",var_nTrackQual2);
  t2->SetBranchAddress("TrackCand_pt",var_TrackPt2);
  t2->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD2);
  t2->SetBranchAddress("TrackCand_purity",var_TrackQuality2);
  t3->SetBranchAddress("nTrackCand",var_nTrack3);
  t3->SetBranchAddress("nQualityTrackCand",var_nTrackQual3);
  t3->SetBranchAddress("TrackCand_pt",var_TrackPt3);
  t3->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD3);
  t3->SetBranchAddress("TrackCand_purity",var_TrackQuality3);
  t2->SetBranchAddress("nTrackCand",var_nTrack2);
  t4->SetBranchAddress("nQualityTrackCand",var_nTrackQual4);
  t4->SetBranchAddress("TrackCand_pt",var_TrackPt4);
  t4->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD4);
  t4->SetBranchAddress("TrackCand_purity",var_TrackQuality4);

  t1->SetBranchAddress("nPrimVertexCand",var_nvtx1);
  t2->SetBranchAddress("nPrimVertexCand",var_nvtx2);
  t1->SetBranchAddress("PrimVertexCand_z",var_vtxZ1);
  t2->SetBranchAddress("PrimVertexCand_z",var_vtxZ2);
  t1->SetBranchAddress("PrimVertexCand_x",var_vtxX1);
  t2->SetBranchAddress("PrimVertexCand_x",var_vtxX2);
  t1->SetBranchAddress("PrimVertexCand_y",var_vtxY1);
  t2->SetBranchAddress("PrimVertexCand_y",var_vtxY2);
  t1->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_1);
  t1->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf1);
  t2->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_2);
  t2->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf2);
  t0->SetBranchAddress("nPrimVertexCand",var_nvtx0);
  t0->SetBranchAddress("PrimVertexCand_z",var_vtxZ0);
  t0->SetBranchAddress("PrimVertexCand_x",var_vtxX0);
  t0->SetBranchAddress("PrimVertexCand_y",var_vtxY0);
  t0->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_0);
  t0->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf0);
  t0->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack0);
  t1->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack1);
  t2->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack2);
  t3->SetBranchAddress("nPrimVertexCand",var_nvtx3);
  t3->SetBranchAddress("PrimVertexCand_z",var_vtxZ3);
  t3->SetBranchAddress("PrimVertexCand_x",var_vtxX3);
  t3->SetBranchAddress("PrimVertexCand_y",var_vtxY3);
  t3->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_3);
  t3->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf3);
  t3->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack3);
  t4->SetBranchAddress("nPrimVertexCand",var_nvtx4);
  t4->SetBranchAddress("PrimVertexCand_z",var_vtxZ4);
  t4->SetBranchAddress("PrimVertexCand_x",var_vtxX4);
  t4->SetBranchAddress("PrimVertexCand_y",var_vtxY4);
  t4->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_4);
  t4->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf4);
  t4->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack4);

  t1->SetBranchAddress("MuMu_vtxx",var_MuMuvtxX1);
  t2->SetBranchAddress("MuMu_vtxx",var_MuMuvtxX2);
  t1->SetBranchAddress("MuMu_vtxy",var_MuMuvtxY1);
  t2->SetBranchAddress("MuMu_vtxy",var_MuMuvtxY2);
  t1->SetBranchAddress("MuMu_vtxz",var_MuMuvtxZ1);
  t2->SetBranchAddress("MuMu_vtxz",var_MuMuvtxZ2);
  t0->SetBranchAddress("MuMu_vtxx",var_MuMuvtxX0);
  t0->SetBranchAddress("MuMu_vtxy",var_MuMuvtxY0);
  t0->SetBranchAddress("MuMu_vtxz",var_MuMuvtxZ0);
  t0->SetBranchAddress("MuMu_vtxisvalid",var_MuMuvtxValid0);
  t1->SetBranchAddress("MuMu_vtxisvalid",var_MuMuvtxValid1);
  t2->SetBranchAddress("MuMu_vtxisvalid",var_MuMuvtxValid2);
  t3->SetBranchAddress("MuMu_vtxx",var_MuMuvtxX3);
  t3->SetBranchAddress("MuMu_vtxy",var_MuMuvtxY3);
  t3->SetBranchAddress("MuMu_vtxz",var_MuMuvtxZ3);
  t3->SetBranchAddress("MuMu_vtxisvalid",var_MuMuvtxValid3);
  t4->SetBranchAddress("MuMu_vtxx",var_MuMuvtxX4);
  t4->SetBranchAddress("MuMu_vtxy",var_MuMuvtxY4);
  t4->SetBranchAddress("MuMu_vtxz",var_MuMuvtxZ4);
  t4->SetBranchAddress("MuMu_vtxisvalid",var_MuMuvtxValid4);

  t1->SetBranchAddress("nCaloCand",var_ncalo1);
  t2->SetBranchAddress("nCaloCand",var_ncalo2);
  t1->SetBranchAddress("CaloTower_ID",var_caloId1);
  t2->SetBranchAddress("CaloTower_ID",var_caloId2);
  t1->SetBranchAddress("CaloTower_e",var_caloEn1);
  t2->SetBranchAddress("CaloTower_e",var_caloEn2);
  t1->SetBranchAddress("CaloTower_t",var_caloTime1);
  t2->SetBranchAddress("CaloTower_t",var_caloTime2);
  t1->SetBranchAddress("CaloTower_dr",var_calodR1);
  t2->SetBranchAddress("CaloTower_dr",var_calodR2);
  t1->SetBranchAddress("Etmiss",var_etmiss1);
  t2->SetBranchAddress("Etmiss",var_etmiss2);
  t1->SetBranchAddress("nExtraCaloTowersE5",var_tower1);
  t2->SetBranchAddress("nExtraCaloTowersE5",var_tower2);
  t1->SetBranchAddress("CaloTower_z",var_caloZ1);
  t2->SetBranchAddress("CaloTower_z",var_caloZ2);
  t0->SetBranchAddress("nCaloCand",var_ncalo0);
  t0->SetBranchAddress("CaloTower_ID",var_caloId0);
  t0->SetBranchAddress("CaloTower_e",var_caloEn0);
  t0->SetBranchAddress("CaloTower_t",var_caloTime0);
  t0->SetBranchAddress("CaloTower_dr",var_calodR0);
  t0->SetBranchAddress("Etmiss",var_etmiss0);
  t0->SetBranchAddress("nExtraCaloTowersE5",var_tower0);
  t0->SetBranchAddress("CaloTower_z",var_caloZ0);
  t3->SetBranchAddress("nCaloCand",var_ncalo3);
  t3->SetBranchAddress("CaloTower_ID",var_caloId3);
  t3->SetBranchAddress("CaloTower_e",var_caloEn3);
  t3->SetBranchAddress("CaloTower_t",var_caloTime3);
  t3->SetBranchAddress("CaloTower_dr",var_calodR3);
  t3->SetBranchAddress("Etmiss",var_etmiss3);
  t3->SetBranchAddress("nExtraCaloTowersE5",var_tower3);
  t3->SetBranchAddress("CaloTower_z",var_caloZ3);
  t4->SetBranchAddress("nCaloCand",var_ncalo4);
  t4->SetBranchAddress("CaloTower_ID",var_caloId4);
  t4->SetBranchAddress("CaloTower_e",var_caloEn4);
  t4->SetBranchAddress("CaloTower_t",var_caloTime4);
  t4->SetBranchAddress("CaloTower_dr",var_calodR4);
  t4->SetBranchAddress("Etmiss",var_etmiss4);
  t4->SetBranchAddress("nExtraCaloTowersE5",var_tower4);
  t4->SetBranchAddress("CaloTower_z",var_caloZ4);

  t0->SetBranchAddress("BX",var_bx0);
  t0->SetBranchAddress("Run",var_run0);
  t0->SetBranchAddress("LumiSection",var_ls0);
  t0->SetBranchAddress("EventNum",var_event0);

  t1->SetBranchAddress("MuMu_mass",var_mass1);
  t1->SetBranchAddress("MuMu_dpt",var_dpt1);
  t1->SetBranchAddress("MuMu_dphi",var_dphi1);
  t0->SetBranchAddress("MuMu_mass",var_mass0);
  t0->SetBranchAddress("MuMu_dpt",var_dpt0);
  t0->SetBranchAddress("MuMu_dphi",var_dphi0);
  t2->SetBranchAddress("MuMu_mass",var_mass2);
  t2->SetBranchAddress("MuMu_dpt",var_dpt2);
  t2->SetBranchAddress("MuMu_dphi",var_dphi2);
  t3->SetBranchAddress("MuMu_mass",var_mass3);
  t3->SetBranchAddress("MuMu_dpt",var_dpt3);
  t3->SetBranchAddress("MuMu_dphi",var_dphi3);
  t4->SetBranchAddress("MuMu_mass",var_mass4);
  t4->SetBranchAddress("MuMu_dpt",var_dpt4);
  t4->SetBranchAddress("MuMu_dphi",var_dphi4);

  t1->SetBranchAddress("MuonCand_isglobal",var_global1);
  t2->SetBranchAddress("MuonCand_isglobal",var_global2);
  t1->SetBranchAddress("MuonCand_istracker",var_tracker1);
  t2->SetBranchAddress("MuonCand_istracker",var_tracker2);
  t1->SetBranchAddress("MuonCand_isstandalone",var_standalone1);
  t2->SetBranchAddress("MuonCand_isstandalone",var_standalone2);
  t1->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack1);
  t2->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack2);
  t0->SetBranchAddress("MuonCand_isglobal",var_global0);
  t0->SetBranchAddress("MuonCand_istracker",var_tracker0);
  t0->SetBranchAddress("MuonCand_isstandalone",var_standalone0);
  t0->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack0);
  t3->SetBranchAddress("MuonCand_isglobal",var_global3);
  t3->SetBranchAddress("MuonCand_istracker",var_tracker3);
  t3->SetBranchAddress("MuonCand_isstandalone",var_standalone3);
  t3->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack3);
  t4->SetBranchAddress("MuonCand_isglobal",var_global4);
  t4->SetBranchAddress("MuonCand_istracker",var_tracker4);
  t4->SetBranchAddress("MuonCand_isstandalone",var_standalone4);
  t4->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack4);

  t1->SetBranchAddress("MuonCand_pt",var_pt1);
  t2->SetBranchAddress("MuonCand_pt",var_pt2);
  t1->SetBranchAddress("MuonCand_pz",var_pz1);
  t2->SetBranchAddress("MuonCand_pz",var_pz2);
  t1->SetBranchAddress("MuonCand_phi",var_phi1);
  t2->SetBranchAddress("MuonCand_phi",var_phi2);
  t1->SetBranchAddress("MuonCand_eta",var_eta1);
  t2->SetBranchAddress("MuonCand_eta",var_eta2);
  t0->SetBranchAddress("MuonCand_pt",var_pt0);
  t0->SetBranchAddress("MuonCand_pz",var_pz0);
  t0->SetBranchAddress("MuonCand_phi",var_phi0);
  t0->SetBranchAddress("MuonCand_eta",var_eta0);
  t3->SetBranchAddress("MuonCand_pt",var_pt3);
  t3->SetBranchAddress("MuonCand_pz",var_pz3);
  t3->SetBranchAddress("MuonCand_phi",var_phi3);
  t3->SetBranchAddress("MuonCand_eta",var_eta3);
  t4->SetBranchAddress("MuonCand_pt",var_pt4);
  t4->SetBranchAddress("MuonCand_pz",var_pz4);
  t4->SetBranchAddress("MuonCand_phi",var_phi4);
  t4->SetBranchAddress("MuonCand_eta",var_eta4);

  t0->SetBranchAddress("MuonPairCand",var_Pair0);
  t1->SetBranchAddress("MuonPairCand",var_Pair1);
  t2->SetBranchAddress("MuonPairCand",var_Pair2);
  t3->SetBranchAddress("MuonPairCand",var_Pair3);
  t4->SetBranchAddress("MuonPairCand",var_Pair4);

  t0->SetBranchAddress("nZDChitCand",var_nZDC0);
  t1->SetBranchAddress("nZDChitCand",var_nZDC1);
  t2->SetBranchAddress("nZDChitCand",var_nZDC2);
  t3->SetBranchAddress("nZDChitCand",var_nZDC3);
  t4->SetBranchAddress("nZDChitCand",var_nZDC4);
  t0->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus0);
  t1->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus1);
  t2->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus2);
  t3->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus3);
  t4->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus4);
  t0->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus0);
  t1->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus1);
  t2->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus2);
  t3->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus3);
  t4->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus4);
  t0->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus0);
  t1->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus1);
  t2->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus2);
  t3->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus3);
  t4->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus4);
  t0->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus0);
  t1->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus1);
  t2->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus2);
  t3->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus3);
  t4->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus4);
  t0->SetBranchAddress("ZDChit_time",var_zdcTime0);
  t1->SetBranchAddress("ZDChit_time",var_zdcTime1);
  t2->SetBranchAddress("ZDChit_time",var_zdcTime2);
  t3->SetBranchAddress("ZDChit_time",var_zdcTime3);
  t4->SetBranchAddress("ZDChit_time",var_zdcTime4);
  t0->SetBranchAddress("ZDChit_energy",var_zdcE0);
  t1->SetBranchAddress("ZDChit_energy",var_zdcE1);
  t2->SetBranchAddress("ZDChit_energy",var_zdcE2);
  t3->SetBranchAddress("ZDChit_energy",var_zdcE3);
  t4->SetBranchAddress("ZDChit_energy",var_zdcE4);
  t0->SetBranchAddress("ZDChit_section",var_zdcsection0);
  t1->SetBranchAddress("ZDChit_section",var_zdcsection1);
  t2->SetBranchAddress("ZDChit_section",var_zdcsection2);
  t3->SetBranchAddress("ZDChit_section",var_zdcsection3);
  t4->SetBranchAddress("ZDChit_section",var_zdcsection4);

  t0->SetBranchAddress("nCastorTowerCand",var_nCastor0);
  t1->SetBranchAddress("nCastorTowerCand",var_nCastor1);
  t2->SetBranchAddress("nCastorTowerCand",var_nCastor2);
  t3->SetBranchAddress("nCastorTowerCand",var_nCastor3);
  t4->SetBranchAddress("nCastorTowerCand",var_nCastor4);
  t0->SetBranchAddress("CastorTower_e",var_CastorE0);
  t1->SetBranchAddress("CastorTower_e",var_CastorE1);
  t2->SetBranchAddress("CastorTower_e",var_CastorE2);
  t3->SetBranchAddress("CastorTower_e",var_CastorE3);
  t4->SetBranchAddress("CastorTower_e",var_CastorE4);
  t0->SetBranchAddress("CastorTower_eta",var_CastorEta0);
  t1->SetBranchAddress("CastorTower_eta",var_CastorEta1);
  t2->SetBranchAddress("CastorTower_eta",var_CastorEta2);
  t3->SetBranchAddress("CastorTower_eta",var_CastorEta3);
  t4->SetBranchAddress("CastorTower_eta",var_CastorEta4);
  t0->SetBranchAddress("CastorTower_phi",var_CastorPhi0);
  t1->SetBranchAddress("CastorTower_phi",var_CastorPhi1);
  t2->SetBranchAddress("CastorTower_phi",var_CastorPhi2);
  t3->SetBranchAddress("CastorTower_phi",var_CastorPhi3);
  t4->SetBranchAddress("CastorTower_phi",var_CastorPhi4);
  t0->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit0);
  t1->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit1);
  t2->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit2);
  t3->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit3);
  t4->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit4);

  t0->SetBranchAddress("L1TechnicalTriggers",techBit0);
  t1->SetBranchAddress("L1TechnicalTriggers",techBit1);
  t2->SetBranchAddress("L1TechnicalTriggers",techBit2);
  t3->SetBranchAddress("L1TechnicalTriggers",techBit3);
  t4->SetBranchAddress("L1TechnicalTriggers",techBit4);

  int filter0Gen(0);
  int filter0Track(0);
  int filter0Events(0);
  for(Int_t i = 0;i < NUM0;i++){
      	t0->GetEntry(i);
	int pair1 = var_Pair0[0]; int pair2 = var_Pair0[1];
        double muID1 = var_idA0[pair1];
        double muID2 = var_idA0[pair2];
	double muAng1 = var_idB0[pair1];
        double muAng2 = var_idB0[pair2];
        int hlt_pass = hlt_d0[0];
	int nPrimVtx = var_nvtx0[0];
	int nTrack=var_nTrack0[0];
	int nTrackQual=var_nTrackQual0[0]; 
	int nCalo=var_ncalo0[0];
	int label_vertex(99);
	double min_distance_vertex(99.0);
//	cout<<"--------------------"<<var_event0[0]<<"-----------------------"<<endl;
	if(nPrimVtx>=1){
          for(Int_t j=0; j<nPrimVtx; j++){
                double distance_vertex_z=(var_MuMuvtxZ0[0]-var_vtxZ0[j] < 0) ? -(var_MuMuvtxZ0[0]-var_vtxZ0[j]) : var_MuMuvtxZ0[0]-var_vtxZ0[j];
                double distance_vertex_x=(var_MuMuvtxX0[0]-var_vtxX0[j] < 0) ? -(var_MuMuvtxX0[0]-var_vtxX0[j]) : var_MuMuvtxX0[0]-var_vtxX0[j];
                double distance_vertex_y=(var_MuMuvtxY0[0]-var_vtxY0[j] < 0) ? -(var_MuMuvtxY0[0]-var_vtxY0[j]) : var_MuMuvtxY0[0]-var_vtxY0[j];
//              cout<<"vtx: nTracks="<<var_vtxTrack0[j]<<"\t d="<<distance_vertex_z<<endl;
                if(var_vtxTrack0[j]==2 && TMath::Prob(var_vertexChi2_0[j],var_vertexNdf0[j]+0.5)>0.001
                   && distance_vertex_z < 0.1 && distance_vertex_z<min_distance_vertex_z && fabs(var_vtxZ0[j])<15.0
                   && distance_vertex_x < 0.071 && distance_vertex_y < 0.071 && var_MuMuvtxValid0[0]==1
                   && (techBit0[0][0]==1))   {label_vertex=j;
                                               min_distance_vertex_z=distance_vertex_z;}
          }

	if(label_vertex!=99
           && var_MuMuvtxValid0[0]==1 && abs(var_MuMuvtxZ0[label_vertex])<=16.0
           && muID1==1 && muID2==1 && muAng1==1 && muAng2==1 /*&& hlt_pass==1*/ 
	   && var_nhitsTrack0[pair1]>12. && var_nhitsTrack0[pair2]>12. && (var_global0[pair1]==1 || var_global0[pair2]==1)
// 	   && var_dpt0[0]<1.5 && (var_dphi0[0]/pi) > 0.9
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(var_TrackQuality0[j]==1){
                  filter0Track++;
		  if(var_TrackD0[j]<0.5) nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0), nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR0[k]>0.3){
		  if(var_caloId0[k]==1 && var_caloEn0[k]>1.15) nEB++;
                  if(var_caloId0[k]==2 && var_caloEn0[k]>2.40) nEE++;
                  if(var_caloId0[k]==4 && var_caloEn0[k]>1.25) nHB++;
                  if(var_caloId0[k]==5 && var_caloEn0[k]>1.90) nHE++;
                  if((var_caloId0[k]==3 ||var_caloId0[k]==6) && var_caloZ0[k]>0 && var_caloEn0[k]>4.2) nHFp++;
                  if((var_caloId0[k]==3 ||var_caloId0[k]==6) && var_caloZ0[k]<0 && var_caloEn0[k]>3.5) nHFm++;
	       }
	    }

	if(nTrackExclu<1 && (nEB+nEE+nHB+nHE+nHFp+nHFm) <200){ 
          filter0Events++;
	cout<<"candidate  Run "<<var_run0[0]<<"  LS "<<var_ls0[0]<<"  Evt "<<var_event0[0]<<endl;

	  hEB0->Fill(nEB,fac_lumi0);
          hEE0->Fill(nEE,fac_lumi0);
          hHB0->Fill(nHB,fac_lumi0);
          hHE0->Fill(nHE,fac_lumi0);
          hHFp0->Fill(nHFp,fac_lumi0);
          hHFm0->Fill(nHFm,fac_lumi0);
	  nTower0->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi0);
	  nTrack0->Fill(nTrackExclu,fac_lumi0);

	  MuMuMass0->Fill(var_mass0[0],fac_lumi0);
          MuMuMassUps0->Fill(var_mass0[0],fac_lumi0);
          MuMudpt0->Fill(var_dpt0[0],fac_lumi0);
          MuMudphi0->Fill(var_dphi0[0]/pi,fac_lumi0);
          MuMudeta0->Fill(fabs(var_eta0[pair1]+var_eta0[pair2]),fac_lumi0);

	  ZDCemminus0->Fill(var_zdcEmMinus0[0],fac_lumi0); ZDCemplus0->Fill(var_zdcEmPlus0[0],fac_lumi0);
	  ZDChadminus0->Fill(var_zdcHadMinus0[0],fac_lumi0); ZDChadplus0->Fill(var_zdcHadPlus0[0],fac_lumi0);
	  for(Int_t l=0; l<var_nZDC0[0]; l++){
	     if(var_zdcsection0[l]==1 && var_zdcE0[l]>7.0){ 
                ZDCtime0->Fill(var_zdcTime0[l],fac_lumi0); ZDCenergyEM0->Fill(var_zdcE0[l],fac_lumi0);
	     }
	     if(var_zdcsection0[l]==2 && var_zdcE0[l]>120.0){
		ZDCtime0->Fill(var_zdcTime0[l],fac_lumi0); ZDCenergyHAD0->Fill(var_zdcE0[l],fac_lumi0);
	     }
	  }

	  CastorSumE0->Fill(var_CastorRecHit0[0],fac_lumi0);
	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"Data :"<<endl;
cout<<"  # Dimuon events = "<<filter0Events<<endl;



  int filter1Gen(0);
  int filter1Track(0);
  int filter1Events(0);
  for(Int_t i = 0;i < NUM1;i++){
      	t1->GetEntry(i);
	int pair1 = var_Pair1[0]; int pair2 = var_Pair1[1];
        double muID1 = var_idA1[pair1];
        double muID2 = var_idA1[pair2];
	double muAng1 = var_idB1[pair1];
        double muAng2 = var_idB1[pair2];
        int hlt_pass = hlt_d1[0];
	int nPrimVtx = var_nvtx1[0];
	int nTrack=var_nTrack1[0];
	int nTrackQual=var_nTrackQual1[0]; 
	int nCalo=var_ncalo1[0];
	int label_vertex(99);
	double min_distance_vertex(99.0);
//	cout<<"--------------------"<<var_event1[0]<<"-----------------------"<<endl;
	if(nPrimVtx>=1){
	  for(Int_t j=0; j<nPrimVtx; j++){
		double distance_vertex=(var_MuMuvtxZ1[0]-var_vtxZ1[j] < 0) ? -(var_MuMuvtxZ1[0]-var_vtxZ1[j]) : var_MuMuvtxZ1[0]-var_vtxZ1[j];
		if(var_vtxTrack1[j]==2 && TMath::Prob(var_vertexChi2_1[j],var_vertexNdf1[j]+0.5)>0.001 
                   && distance_vertex <0.1 && distance_vertex<min_distance_vertex && abs(var_vtxZ1[j])<16.0
		   && !(techBit1[0][0]==1))   {label_vertex=j;
				  	       min_distance_vertex=distance_vertex;}
	  }
	}

	if(label_vertex!=99
           && var_MuMuvtxValid1[0]==1 && abs(var_MuMuvtxZ1[label_vertex])<=16.0
           && muID1==1 && muID2==1 && muAng1==1 && muAng2==1 /*&& hlt_pass==1*/ 
	   && var_nhitsTrack1[pair1]>12. && var_nhitsTrack1[pair2]>12. && (var_global1[pair1]==1 || var_global1[pair2]==1)
// 	   && var_dpt1[0]<1.5 && (var_dphi1[0]/pi) > 0.9
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(var_TrackQuality1[j]==1){
                  filter1Track++;
		  if(var_TrackD1[j]<0.5) nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR1[k]>0.3){
		  if(var_caloId1[k]==1 && var_caloEn1[k]>1.15) nEB++;
                  if(var_caloId1[k]==2 && var_caloEn1[k]>2.40) nEE++;
                  if(var_caloId1[k]==4 && var_caloEn1[k]>1.25) nHB++;
                  if(var_caloId1[k]==5 && var_caloEn1[k]>1.90) nHE++;
                  if((var_caloId1[k]==3 ||var_caloId1[k]==6) && var_caloZ0[k]>0 && var_caloEn1[k]>4.2) nHFp++;
                  if((var_caloId1[k]==3 ||var_caloId1[k]==6) && var_caloZ0[k]<0 && var_caloEn1[k]>3.5) nHFm++;
	       }
	    }

	if(nTrackExclu<1 && (nEB+nEE+nHB+nHE+nHFp+nHFm) <200){ 
          filter1Events++;

	  hEB3->Fill(nEB,fac_lumi1);          hEB4->Fill(nEB,fac_lumi1);
          hEE3->Fill(nEE,fac_lumi1);          hEE4->Fill(nEE,fac_lumi1);
          hHB3->Fill(nHB,fac_lumi1);          hHB4->Fill(nHB,fac_lumi1);
          hHE3->Fill(nHE,fac_lumi1);          hHE4->Fill(nHE,fac_lumi1);
          hHFp3->Fill(nHFp,fac_lumi1);          hHFp4->Fill(nHFp,fac_lumi1);
          hHFm3->Fill(nHFm,fac_lumi1);          hHFm4->Fill(nHFm,fac_lumi1);
	  nTower3->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi1);          nTower4->Fill(nEB+nEE+nHB+nHE+nHFp,fac_lumi1);
	  nTrack3->Fill(nTrackExclu,fac_lumi1);          nTrack4->Fill(nTrackExclu,fac_lumi1);

	  MuMuMass3->Fill(var_mass1[0],fac_lumi1);          MuMuMass4->Fill(var_mass1[0],fac_lumi1);
          MuMuMassUps3->Fill(var_mass1[0],fac_lumi1);          MuMuMassUps4->Fill(var_mass1[0],fac_lumi1);
          MuMudpt3->Fill(var_dpt1[0],fac_lumi1);          MuMudpt4->Fill(var_dpt1[0],fac_lumi1);
          MuMudphi3->Fill(var_dphi1[0]/pi,fac_lumi1);          MuMudphi4->Fill(var_dphi1[0]/pi,fac_lumi1);
          MuMudeta3->Fill(fabs(var_eta1[pair1]+var_eta1[pair2]),fac_lumi1);          MuMudeta4->Fill(fabs(var_eta1[pair1]+var_eta1[pair2]),fac_lumi1);

	  ZDCemminus3->Fill(var_zdcEmMinus1[0],fac_lumi1);          ZDCemminus4->Fill(var_zdcEmMinus1[0],fac_lumi1);
	  ZDCemplus3->Fill(var_zdcEmPlus1[0],fac_lumi1);          ZDCemplus4->Fill(var_zdcEmPlus1[0],fac_lumi1);
	  ZDChadminus3->Fill(var_zdcHadMinus1[0],fac_lumi1);          ZDChadminus4->Fill(var_zdcHadMinus1[0],fac_lumi1);
	  ZDChadplus3->Fill(var_zdcHadPlus1[0],fac_lumi1);          ZDChadplus4->Fill(var_zdcHadPlus1[0],fac_lumi1);
	  for(Int_t l=0; l<var_nZDC1[0]; l++){
	     if(var_zdcsection1[l]==1 && var_zdcE1[l]>7.0){ 
                ZDCtime3->Fill(var_zdcTime1[l],fac_lumi1);                 ZDCtime4->Fill(var_zdcTime1[l],fac_lumi1);
                ZDCenergyEM3->Fill(var_zdcE1[l],fac_lumi1);                ZDCenergyEM4->Fill(var_zdcE1[l],fac_lumi1);
	     }
	     if(var_zdcsection1[l]==2 && var_zdcE1[l]>120.0){
		ZDCtime3->Fill(var_zdcTime1[l],fac_lumi1);                ZDCtime4->Fill(var_zdcTime1[l],fac_lumi1);
                ZDCenergyHAD3->Fill(var_zdcE1[l],fac_lumi1);                ZDCenergyHAD4->Fill(var_zdcE1[l],fac_lumi1);
	     }
	  }

	  CastorSumE3->Fill(var_CastorRecHit1[0],fac_lumi1);          CastorSumE4->Fill(var_CastorRecHit1[0],fac_lumi1);
	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"ElEl :"<<endl;
cout<<"  # Dimuon events = "<<filter1Events*fac_lumi1<<endl;

  int filter2Gen(0);
  int filter2Track(0);
  int filter2Events(0);
  for(Int_t i = 0;i < NUM2;i++){
      	t2->GetEntry(i);
	int pair2 = var_Pair2[0]; int pair2 = var_Pair2[1];
        double muID1 = var_idA2[pair1];
        double muID2 = var_idA2[pair2];
	double muAng1 = var_idB2[pair1];
        double muAng2 = var_idB2[pair2];
        int hlt_pass = hlt_d2[0];
	int nPrimVtx = var_nvtx2[0];
	int nTrack=var_nTrack2[0];
	int nTrackQual=var_nTrackQual2[0]; 
	int nCalo=var_ncalo2[0];
	int label_vertex(99);
	double min_distance_vertex(99.0);
//	cout<<"--------------------"<<var_event2[0]<<"-----------------------"<<endl;
	if(nPrimVtx>=1){
	  for(Int_t j=0; j<nPrimVtx; j++){
		double distance_vertex=(var_MuMuvtxZ2[0]-var_vtxZ2[j] < 0) ? -(var_MuMuvtxZ2[0]-var_vtxZ2[j]) : var_MuMuvtxZ2[0]-var_vtxZ2[j];
		if(var_vtxTrack2[j]==2 && TMath::Prob(var_vertexChi2_2[j],var_vertexNdf2[j]+0.5)>0.001 
                   && distance_vertex <0.1 && distance_vertex<min_distance_vertex && abs(var_vtxZ2[j])<16.0
		   && !(techBit2[0][0]==1))   {label_vertex=j;
				  	       min_distance_vertex=distance_vertex;}
	  }
	}

	if(label_vertex!=99
           && var_MuMuvtxValid2[0]==1 && abs(var_MuMuvtxZ2[label_vertex])<=16.0
           && muID1==1 && muID2==1 && muAng1==1 && muAng2==1 /*&& hlt_pass==1*/ 
	   && var_nhitsTrack2[pair1]>12. && var_nhitsTrack2[pair2]>12. && (var_global2[pair1]==1 || var_global2[pair2]==1)
// 	   && var_dpt2[0]<1.5 && (var_dphi2[0]/pi) > 0.9
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(var_TrackQuality2[j]==1){
                  filter2Track++;
		  if(var_TrackD2[j]<0.5) nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR2[k]>0.3){
		  if(var_caloId2[k]==1 && var_caloEn2[k]>1.15) nEB++;
                  if(var_caloId2[k]==2 && var_caloEn2[k]>2.40) nEE++;
                  if(var_caloId2[k]==4 && var_caloEn2[k]>1.25) nHB++;
                  if(var_caloId2[k]==5 && var_caloEn2[k]>1.90) nHE++;
                  if((var_caloId2[k]==3 ||var_caloId2[k]==6) && var_caloZ0[k]>0 && var_caloEn2[k]>4.2) nHFp++;
                  if((var_caloId2[k]==3 ||var_caloId2[k]==6) && var_caloZ0[k]<0 && var_caloEn2[k]>3.5) nHFm++;
	       }
	    }

	if(nTrackExclu<1 && (nEB+nEE+nHB+nHE+nHFp+nHFm) <200){ 
          filter2Events++;

	  hEB4->Fill(nEB,fac_lumi2); hEB3->Fill(nEB,fac_lumi2); hEB2->Fill(nEB,fac_lumi2);
          hEE4->Fill(nEE,fac_lumi2); hEE3->Fill(nEE,fac_lumi2); hEE2->Fill(nEE,fac_lumi2);
          hHB4->Fill(nHB,fac_lumi2); hHB3->Fill(nHB,fac_lumi2); hHB2->Fill(nHB,fac_lumi2);
          hHE4->Fill(nHE,fac_lumi2); hHE3->Fill(nHE,fac_lumi2); hHE2->Fill(nHE,fac_lumi2);
          hHFp4->Fill(nHFp,fac_lumi2); hHFp3->Fill(nHFp,fac_lumi2); hHFp2->Fill(nHFp,fac_lumi2);
          hHFm4->Fill(nHFm,fac_lumi2); hHFm3->Fill(nHFm,fac_lumi2); hHFm2->Fill(nHFm,fac_lumi2);
	  nTower4->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi2);nTower3->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi2);nTower2->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi2);
	  nTrack4->Fill(nTrackExclu,fac_lumi2);nTrack3->Fill(nTrackExclu,fac_lumi2);nTrack2->Fill(nTrackExclu,fac_lumi2);

	  MuMuMass4->Fill(var_mass2[0],fac_lumi2); MuMuMass3->Fill(var_mass2[0],fac_lumi2); MuMuMass2->Fill(var_mass2[0],fac_lumi2);
          MuMuMassUps4->Fill(var_mass2[0],fac_lumi2); MuMuMassUps3->Fill(var_mass2[0],fac_lumi2); MuMuMassUps2->Fill(var_mass2[0],fac_lumi2);
          MuMudpt4->Fill(var_dpt2[0],fac_lumi2);  MuMudpt3->Fill(var_dpt2[0],fac_lumi2);  MuMudpt2->Fill(var_dpt2[0],fac_lumi2);
          MuMudphi4->Fill(var_dphi2[0]/pi,fac_lumi2); MuMudphi3->Fill(var_dphi2[0]/pi,fac_lumi2);  MuMudphi2->Fill(var_dphi2[0]/pi,fac_lumi2);
          MuMudeta4->Fill(fabs(var_eta2[pair1]+var_eta2[pair2]),fac_lumi2); MuMudeta3->Fill(fabs(var_eta2[pair1]+var_eta2[pair2]),fac_lumi2); MuMudeta2->Fill(fabs(var_eta2[pair1]+var_eta2[pair2]),fac_lumi2);

	  ZDCemminus4->Fill(var_zdcEmMinus2[0],fac_lumi2); ZDCemminus3->Fill(var_zdcEmMinus2[0],fac_lumi2); ZDCemminus2->Fill(var_zdcEmMinus2[0],fac_lumi2);
	  ZDCemplus4->Fill(var_zdcEmPlus2[0],fac_lumi2); ZDCemplus3->Fill(var_zdcEmPlus2[0],fac_lumi2); ZDCemplus2->Fill(var_zdcEmPlus2[0],fac_lumi2);
	  ZDChadminus4->Fill(var_zdcHadMinus2[0],fac_lumi2); ZDChadminus3->Fill(var_zdcHadMinus2[0],fac_lumi2); ZDChadminus2->Fill(var_zdcHadMinus2[0],fac_lumi2);
	  ZDChadplus4->Fill(var_zdcHadPlus2[0],fac_lumi2); ZDChadplus3->Fill(var_zdcHadPlus2[0],fac_lumi2); ZDChadplus2->Fill(var_zdcHadPlus2[0],fac_lumi2);
	  for(Int_t l=0; l<var_nZDC2[0]; l++){
	     if(var_zdcsection2[l]==1 && var_zdcE2[l]>7.0){ 
                ZDCtime4->Fill(var_zdcTime2[l],fac_lumi2); ZDCtime3->Fill(var_zdcTime2[l],fac_lumi2); ZDCtime2->Fill(var_zdcTime2[l],fac_lumi2);
                ZDCenergyEM4->Fill(var_zdcE2[l],fac_lumi2); ZDCenergyEM3->Fill(var_zdcE2[l],fac_lumi2); ZDCenergyEM2->Fill(var_zdcE2[l],fac_lumi2);
	     }
	     if(var_zdcsection2[l]==2 && var_zdcE2[l]>120.0){
		ZDCtime4->Fill(var_zdcTime2[l],fac_lumi2); ZDCtime3->Fill(var_zdcTime2[l],fac_lumi2); ZDCtime2->Fill(var_zdcTime2[l],fac_lumi2);
                ZDCenergyHAD4->Fill(var_zdcE2[l],fac_lumi2);ZDCenergyHAD3->Fill(var_zdcE2[l],fac_lumi2);ZDCenergyHAD2->Fill(var_zdcE2[l],fac_lumi2);
	     }
	  }

	  CastorSumE4->Fill(var_CastorRecHit2[0],fac_lumi2); CastorSumE3->Fill(var_CastorRecHit2[0],fac_lumi2); CastorSumE2->Fill(var_CastorRecHit2[0],fac_lumi2);
	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"InelEl :"<<endl;
cout<<"  # Dimuon events = "<<filter2Events*fac_lumi2<<endl;

  int filter3Gen(0);
  int filter3Track(0);
  int filter3Events(0);
  for(Int_t i = 0;i < NUM3;i++){
      	t3->GetEntry(i);
	int pair1 = var_Pair3[0]; int pair2 = var_Pair3[1];
        double muID1 = var_idA3[pair1];
        double muID2 = var_idA3[pair2];
	double muAng1 = var_idB3[pair1];
        double muAng2 = var_idB3[pair2];
        int hlt_pass = hlt_d3[0];
	int nPrimVtx = var_nvtx3[0];
	int nTrack=var_nTrack3[0];
	int nTrackQual=var_nTrackQual3[0]; 
	int nCalo=var_ncalo3[0];
	int label_vertex(99);
	double min_distance_vertex(99.0);
//	cout<<"--------------------"<<var_event1[0]<<"-----------------------"<<endl;
	if(nPrimVtx>=1){
	  for(Int_t j=0; j<nPrimVtx; j++){
		double distance_vertex=(var_MuMuvtxZ3[0]-var_vtxZ3[j] < 0) ? -(var_MuMuvtxZ3[0]-var_vtxZ3[j]) : var_MuMuvtxZ3[0]-var_vtxZ3[j];
		if(var_vtxTrack3[j]==2 && TMath::Prob(var_vertexChi2_3[j],var_vertexNdf3[j]+0.5)>0.001 
                   && distance_vertex <0.1 && distance_vertex<min_distance_vertex && abs(var_vtxZ3[j])<16.0
		   && !(techBit3[0][0]==1))   {label_vertex=j;
				  	       min_distance_vertex=distance_vertex;}
	  }
	}

	if(label_vertex!=99
           && var_MuMuvtxValid3[0]==1 && abs(var_MuMuvtxZ3[label_vertex])<=16.0
           && muID1==1 && muID2==1 && muAng1==1 && muAng2==1 /*&& hlt_pass==1*/ 
	   && var_nhitsTrack3[pair1]>12. && var_nhitsTrack3[pair2]>12. && (var_global3[pair1]==1 || var_global3[pair2]==1)
// 	   && var_dpt3[0]<1.5 && (var_dphi3[0]/pi) > 0.9
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(var_TrackQuality3[j]==1){
                  filter3Track++;
		  if(var_TrackD3[j]<0.5) nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR3[k]>0.3){
		  if(var_caloId3[k]==1 && var_caloEn3[k]>1.15) nEB++;
                  if(var_caloId3[k]==2 && var_caloEn3[k]>2.40) nEE++;
                  if(var_caloId3[k]==4 && var_caloEn3[k]>1.25) nHB++;
                  if(var_caloId3[k]==5 && var_caloEn3[k]>1.90) nHE++;
                  if((var_caloId3[k]==3 ||var_caloId3[k]==6) && var_caloZ0[k]>0 && var_caloEn3[k]>4.2) nHFp++;
                  if((var_caloId3[k]==3 ||var_caloId3[k]==6) && var_caloZ0[k]<0 && var_caloEn3[k]>3.5) nHFm++;
	       }
	    }

	if(nTrackExclu<1 && (nEB+nEE+nHB+nHE+nHFp+nHFm) <200){ 
          filter3Events++;

	  hEB1->Fill(nEB,fac_lumi3); hEB2->Fill(nEB,fac_lumi3); hEB3->Fill(nEB,fac_lumi3); hEB4->Fill(nEB,fac_lumi3);
          hEE1->Fill(nEE,fac_lumi3); hEE2->Fill(nEE,fac_lumi3); hEE3->Fill(nEE,fac_lumi3); hEE4->Fill(nEE,fac_lumi3);
          hHB1->Fill(nHB,fac_lumi3); hHB2->Fill(nHB,fac_lumi3); hHB3->Fill(nHB,fac_lumi3); hHB4->Fill(nHB,fac_lumi3);
          hHE1->Fill(nHE,fac_lumi3); hHE2->Fill(nHE,fac_lumi3); hHE3->Fill(nHE,fac_lumi3); hHE4->Fill(nHE,fac_lumi3);
          hHFp1->Fill(nHFp,fac_lumi3); hHFp2->Fill(nHFp,fac_lumi3); hHFp3->Fill(nHFp,fac_lumi3); hHFp4->Fill(nHFp,fac_lumi3);
          hHFm1->Fill(nHFm,fac_lumi3); hHFm2->Fill(nHFm,fac_lumi3); hHFm3->Fill(nHFm,fac_lumi3); hHFm4->Fill(nHFm,fac_lumi3);
	  nTower1->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi3);nTower2->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi3);nTower3->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi3);;nTower4->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi3);
	  nTrack1->Fill(nTrackExclu,fac_lumi3);nTrack2->Fill(nTrackExclu,fac_lumi3);nTrack3->Fill(nTrackExclu,fac_lumi3);nTrack4->Fill(nTrackExclu,fac_lumi3);

	  MuMuMass1->Fill(var_mass3[0],fac_lumi3); MuMuMass2->Fill(var_mass3[0],fac_lumi3); MuMuMass3->Fill(var_mass3[0],fac_lumi3); MuMuMass4->Fill(var_mass3[0],fac_lumi3);
          MuMuMassUps1->Fill(var_mass3[0],fac_lumi3); MuMuMassUps2->Fill(var_mass3[0],fac_lumi3); MuMuMassUps3->Fill(var_mass3[0],fac_lumi3); MuMuMassUps4->Fill(var_mass3[0],fac_lumi3);
          MuMudpt1->Fill(var_dpt3[0],fac_lumi3);  MuMudpt2->Fill(var_dpt3[0],fac_lumi3);  MuMudpt3->Fill(var_dpt3[0],fac_lumi3);  MuMudpt4->Fill(var_dpt3[0],fac_lumi3);
          MuMudphi1->Fill(var_dphi3[0]/pi,fac_lumi3);  MuMudphi2->Fill(var_dphi3[0]/pi,fac_lumi3);  MuMudphi3->Fill(var_dphi3[0]/pi,fac_lumi3); MuMudphi4->Fill(var_dphi3[0]/pi,fac_lumi3);
          MuMudeta1->Fill(fabs(var_eta3[pair1]+var_eta3[pair2]),fac_lumi3); MuMudeta2->Fill(fabs(var_eta3[pair1]+var_eta3[pair2]),fac_lumi3); MuMudeta3->Fill(fabs(var_eta3[pair1]+var_eta3[pair2]),fac_lumi3); MuMudeta4->Fill(fabs(var_eta3[pair1]+var_eta3[pair2]),fac_lumi3);

	  ZDCemminus1->Fill(var_zdcEmMinus3[0],fac_lumi3); ZDCemminus2->Fill(var_zdcEmMinus3[0],fac_lumi3);ZDCemminus3->Fill(var_zdcEmMinus3[0],fac_lumi3);ZDCemminus4->Fill(var_zdcEmMinus3[0],fac_lumi3);
	  ZDCemplus1->Fill(var_zdcEmPlus3[0],fac_lumi3); ZDCemplus2->Fill(var_zdcEmPlus3[0],fac_lumi3); ZDCemplus3->Fill(var_zdcEmPlus3[0],fac_lumi3); ZDCemplus4->Fill(var_zdcEmPlus3[0],fac_lumi3);
	  ZDChadminus1->Fill(var_zdcHadMinus3[0],fac_lumi3); ZDChadminus2->Fill(var_zdcHadMinus3[0],fac_lumi3); ZDChadminus3->Fill(var_zdcHadMinus3[0],fac_lumi3); ZDChadminus4->Fill(var_zdcHadMinus3[0],fac_lumi3);
	  ZDChadplus1->Fill(var_zdcHadPlus3[0],fac_lumi3); ZDChadplus2->Fill(var_zdcHadPlus3[0],fac_lumi3); ZDChadplus3->Fill(var_zdcHadPlus3[0],fac_lumi3); ZDChadplus4->Fill(var_zdcHadPlus3[0],fac_lumi3);
	  for(Int_t l=0; l<var_nZDC3[0]; l++){
	     if(var_zdcsection3[l]==1 && var_zdcE3[l]>7.0){ 
                ZDCtime1->Fill(var_zdcTime3[l],fac_lumi3); ZDCtime2->Fill(var_zdcTime3[l],fac_lumi3); ZDCtime3->Fill(var_zdcTime3[l],fac_lumi3);  ZDCtime4->Fill(var_zdcTime3[l],fac_lumi3);
                ZDCenergyEM1->Fill(var_zdcE3[l],fac_lumi3);  ZDCenergyEM2->Fill(var_zdcE3[l],fac_lumi3); ZDCenergyEM3->Fill(var_zdcE3[l],fac_lumi3); ZDCenergyEM4->Fill(var_zdcE3[l],fac_lumi3);
	     }
	     if(var_zdcsection3[l]==2 && var_zdcE3[l]>120.0){
		ZDCtime1->Fill(var_zdcTime3[l],fac_lumi3); ZDCtime2->Fill(var_zdcTime3[l],fac_lumi3);ZDCtime3->Fill(var_zdcTime3[l],fac_lumi3);ZDCtime4->Fill(var_zdcTime3[l],fac_lumi3);
                ZDCenergyHAD1->Fill(var_zdcE3[l],fac_lumi3);ZDCenergyHAD2->Fill(var_zdcE3[l],fac_lumi3);ZDCenergyHAD3->Fill(var_zdcE3[l],fac_lumi3);ZDCenergyHAD4->Fill(var_zdcE3[l],fac_lumi3);
	     }
	  }

	  CastorSumE1->Fill(var_CastorRecHit3[0],fac_lumi3); CastorSumE2->Fill(var_CastorRecHit3[0],fac_lumi3);CastorSumE3->Fill(var_CastorRecHit3[0],fac_lumi3);CastorSumE4->Fill(var_CastorRecHit3[0],fac_lumi3);

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"InelInel :"<<endl;
cout<<"  # Dimuon events = "<<filter3Events*fac_lumi3<<endl;


  int filter4Gen(0);
  int filter4Track(0);
  int filter4Events(0);
  for(Int_t i = 0;i < NUM4;i++){
      	t4->GetEntry(i);
	int pair1 = var_Pair4[0]; int pair2 = var_Pair4[1];
        double muID1 = var_idA4[pair1];
        double muID2 = var_idA4[pair2];
	double muAng1 = var_idB4[pair1];
        double muAng2 = var_idB4[pair2];
        int hlt_pass = hlt_d4[0];
	int nPrimVtx = var_nvtx4[0];
	int nTrack=var_nTrack4[0];
	int nTrackQual=var_nTrackQual4[0]; 
	int nCalo=var_ncalo4[0];
	int label_vertex(99);
	double min_distance_vertex(99.0);
//	cout<<"--------------------"<<var_event4[0]<<"-----------------------"<<endl;
	if(nPrimVtx>=1){
	  for(Int_t j=0; j<nPrimVtx; j++){
		double distance_vertex=(var_MuMuvtxZ4[0]-var_vtxZ4[j] < 0) ? -(var_MuMuvtxZ4[0]-var_vtxZ4[j]) : var_MuMuvtxZ4[0]-var_vtxZ4[j];
//		cout<<"vtx: nTracks="<<var_vtxTrack4[j]<<"\t d="<<distance_vertex<<endl;
		if(var_vtxTrack4[j]>=2 && TMath::Prob(var_vertexChi2_4[j],var_vertexNdf4[j]+0.5)>0.001 
                   && distance_vertex <0.1 && distance_vertex<min_distance_vertex && abs(var_vtxZ4[j])<16.0
		   && !(techBit4[0][0]==1))   {label_vertex=j;
				  	       min_distance_vertex=distance_vertex;}
	  }
	}

	if(label_vertex!=99
           && var_MuMuvtxValid4[0]==1 && abs(var_MuMuvtxZ4[label_vertex])<=16.0
           && muID1==1 && muID2==1 && muAng1==1 && muAng2==1 /*&& hlt_pass==1*/ 
	   && var_nhitsTrack4[pair1]>12. && var_nhitsTrack4[pair2]>12. && (var_global4[pair1]==1 || var_global4[pair2]==1)
// 	   && var_dpt4[0]<1.5 && (var_dphi4[0]/pi) > 0.9
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(var_TrackQuality4[j]==1){
                  filter4Track++;
		  if(var_TrackD4[j]<0.5) nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR4[k]>0.3){
		  if(var_caloId4[k]==1 && var_caloEn4[k]>1.15) nEB++;
                  if(var_caloId4[k]==2 && var_caloEn4[k]>2.40) nEE++;
                  if(var_caloId4[k]==4 && var_caloEn4[k]>1.25) nHB++;
                  if(var_caloId4[k]==5 && var_caloEn4[k]>1.90) nHE++;
                  if((var_caloId4[k]==3 ||var_caloId4[k]==6) && var_caloZ0[k]>0 && var_caloEn4[k]>4.2) nHFp++;
                  if((var_caloId4[k]==3 ||var_caloId4[k]==6) && var_caloZ0[k]<0 && var_caloEn4[k]>4.2) nHFm++;
	       }
	    }

	if(nTrackExclu<1 && (nEB+nEE+nHB+nHE+nHFp+nHFm) <200){ 
          filter4Events++;

	  hEB4->Fill(nEB,fac_lumi4);
          hEE4->Fill(nEE,fac_lumi4);
          hHB4->Fill(nHB,fac_lumi4);
          hHE4->Fill(nHE,fac_lumi4);
          hHFp4->Fill(nHFp,fac_lumi4);
          hHFm4->Fill(nHFm,fac_lumi4);
	  nTower4->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi4);
	  nTrack4->Fill(nTrackExclu,fac_lumi4);

	  MuMuMass4->Fill(var_mass4[0],fac_lumi4);
          MuMuMassUps4->Fill(var_mass4[0],fac_lumi4);
          MuMudpt4->Fill(var_dpt4[0],fac_lumi4);
          MuMudphi4->Fill(var_dphi4[0]/pi,fac_lumi4);
          MuMudeta4->Fill(fabs(var_eta4[pair1]+var_eta4[pair2]),fac_lumi4);

	  ZDCemminus4->Fill(var_zdcEmMinus4[0],fac_lumi4); ZDCemplus4->Fill(var_zdcEmPlus4[0],fac_lumi4);
	  ZDChadminus4->Fill(var_zdcHadMinus4[0],fac_lumi4); ZDChadplus4->Fill(var_zdcHadPlus4[0],fac_lumi4);
	  for(Int_t l=0; l<var_nZDC4[0]; l++){
	     if(var_zdcsection4[l]==1 && var_zdcE4[l]>7.0){ 
                ZDCtime4->Fill(var_zdcTime4[l],fac_lumi4); ZDCenergyEM4->Fill(var_zdcE4[l],fac_lumi4);
	     }
	     if(var_zdcsection4[l]==2 && var_zdcE4[l]>120.0){
		ZDCtime4->Fill(var_zdcTime4[l],fac_lumi4); ZDCenergyHAD4->Fill(var_zdcE4[l],fac_lumi4);
	     }
	  }

	  CastorSumE4->Fill(var_CastorRecHit4[0],fac_lumi4);
	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"Upsilon :"<<endl;
cout<<"  # Dimuon events = "<<filter4Events*fac_lumi4<<endl;



if(1==1){

   ci = TColor::GetColor("#ffff99");

// Save into histo
//Draw
TCanvas *PtSpectrumEl = new TCanvas("PtSpectrumEl","Tracks",800,500);
   PtSpectrumEl->SetFillColor(0);
   PtSpectrumEl->SetBorderMode(0);
   PtSpectrumEl->SetBorderSize(2);
   PtSpectrumEl->SetFrameBorderMode(0);

//   leg->Draw();

PtSpectrumEl->cd(1);
nTrack4->GetXaxis()->SetTitle("n Extra Tracks");
nTrack4->GetYaxis()->SetTitle("# events");
nTrack0->Sumw2();
nTrack0->SetLineWidth(2);
nTrack0->SetMarkerStyle(20);
nTrack1->SetFillColor(30);
nTrack2->SetFillColor(30);
nTrack2->SetFillStyle(3001);
nTrack3->SetFillColor(ci);
nTrack4->SetFillColor(38);
nTrack4->Draw("hist");
nTrack3->Draw("same");
nTrack2->Draw("same");
nTrack1->Draw("same");
nTrack0->Draw("same");


//---------------------------
TCanvas *Calo = new TCanvas("Calo","Calorimeter",800,500);
   Calo->SetFillColor(0);
   Calo->SetBorderMode(0);
   Calo->SetBorderSize(2);
   Calo->SetFrameBorderMode(0);
Calo->Divide(3,2);
Calo->cd(4);
hEB4->GetXaxis()->SetTitle("n EB towers (E>1.15GeV)");
hEB4->GetYaxis()->SetTitle("# events");
hEB0->Sumw2();
hEB0->SetLineWidth(2);
hEB0->SetMarkerStyle(20);
hEB1->SetFillColor(30);
hEB2->SetFillColor(30);
hEB2->SetFillStyle(3001);
hEB3->SetFillColor(ci);
hEB4->SetFillColor(38);
hEB4->Draw("hist");
hEB3->Draw("same");
hEB2->Draw("same");
hEB1->Draw("same");
hEB0->Draw("same");


Calo->cd(5);
hEE4->GetXaxis()->SetTitle("n EE towers (E>2.40GeV)");
hEE4->GetYaxis()->SetTitle("# events");
hEE0->Sumw2();
hEE0->SetLineWidth(2);
hEE0->SetMarkerStyle(20);
hEE1->SetFillColor(30);
hEE2->SetFillColor(30);
hEE2->SetFillStyle(3001);
hEE3->SetFillColor(ci);
hEE4->SetFillColor(38);
hEE4->Draw("hist");
hEE3->Draw("same");
hEE2->Draw("same");
hEE1->Draw("same");
hEE0->Draw("same");

Calo->cd(1);
hHB4->GetXaxis()->SetTitle("n HB towers (E>1.25GeV)");
hHB4->GetYaxis()->SetTitle("# events");
hHB0->Sumw2();
hHB0->SetLineWidth(2);
hHB0->SetMarkerStyle(20);
hHB1->SetFillColor(30);
hHB2->SetFillColor(30);
hHB2->SetFillStyle(3001);
hHB3->SetFillColor(ci);
hHB4->SetFillColor(38);
hHB4->Draw("hist");
hHB3->Draw("same");
hHB2->Draw("same");
hHB1->Draw("same");
hHB0->Draw("same");

Calo->cd(2);
hHE4->GetXaxis()->SetTitle("n HE towers (E>1.25GeV)");
hHE4->GetYaxis()->SetTitle("# events");
hHE0->Sumw2();
hHE0->SetLineWidth(2);
hHE0->SetMarkerStyle(20);
hHE1->SetFillColor(30);
hHE2->SetFillColor(30);
hHE2->SetFillStyle(3001);
hHE3->SetFillColor(ci);
hHE4->SetFillColor(38);
hHE4->Draw("hist");
hHE3->Draw("same");
hHE2->Draw("same");
hHE1->Draw("same");
hHE0->Draw("same");

Calo->cd(3);
hHFp4->GetXaxis()->SetTitle("n HF+ towers (E>4.2GeV)");
hHFp4->GetYaxis()->SetTitle("# events");
hHFp0->Sumw2();
hHFp0->SetLineWidth(2);
hHFp0->SetMarkerStyle(20);
hHFp1->SetFillColor(30);
hHFp2->SetFillColor(30);
hHFp2->SetFillStyle(3001);
hHFp3->SetFillColor(ci);
hHFp4->SetFillColor(38);
hHFp4->Draw("hist");
hHFp3->Draw("same");
hHFp2->Draw("same");
hHFp1->Draw("same");
hHFp0->Draw("same");

Calo->cd(6);
hHFm4->GetXaxis()->SetTitle("n HF- towers (E>3.5GeV)");
hHFm4->GetYaxis()->SetTitle("# events");
hHFm0->Sumw2();
hHFm0->SetLineWidth(2);
hHFm0->SetMarkerStyle(20);
hHFm1->SetFillColor(30);
hHFm2->SetFillColor(30);
hHFm2->SetFillStyle(3001);
hHFm3->SetFillColor(ci);
hHFm4->SetFillColor(38);
hHFm4->Draw("hist");
hHFm3->Draw("same");
hHFm2->Draw("same");
hHFm1->Draw("same");
hHFm0->Draw("same");

TCanvas *Calo2 = new TCanvas("Calo2","Calorimeter 2",800,500);
   Calo2->SetFillColor(0);
   Calo2->SetBorderMode(0);
   Calo2->SetBorderSize(2);
   Calo2->SetFrameBorderMode(0);
Calo2->cd(1);
nTower4->GetXaxis()->SetTitle("n towers (E>E_{Threshold})");
nTower4->GetYaxis()->SetTitle("# events");
nTower0->Sumw2();
nTower0->SetLineWidth(2);
nTower0->SetMarkerStyle(20);
nTower1->SetFillColor(30);
nTower2->SetFillColor(30);
nTower2->SetFillStyle(3001);
nTower3->SetFillColor(ci);
nTower4->SetFillColor(38);
nTower4->Draw("hist");
nTower3->Draw("same");
nTower2->Draw("same");
nTower1->Draw("same");
nTower0->Draw("same");

TCanvas *Upsilon = new TCanvas("Upsilon","Upsilon",800,500);
   Upsilon->SetFillColor(0);
   Upsilon->SetBorderMode(0);
   Upsilon->SetBorderSize(2);
   Upsilon->SetFrameBorderMode(0);

Upsilon->cd(1);
MuMuMassUps4->GetXaxis()->SetTitle("#mu#mu mass [GeV]");
MuMuMassUps4->GetYaxis()->SetTitle("# events / 0.5 GeV");
MuMuMassUps0->Sumw2();
MuMuMassUps0->SetLineWidth(2);
MuMuMassUps0->SetMarkerStyle(20);
MuMuMassUps1->SetFillColor(30);
MuMuMassUps2->SetFillColor(30);
MuMuMassUps2->SetFillStyle(3001);
MuMuMassUps3->SetFillColor(ci);
MuMuMassUps4->SetFillColor(38);
MuMuMassUps4->Draw("hist");
MuMuMassUps3->Draw("same");
MuMuMassUps2->Draw("same");
MuMuMassUps1->Draw("same");
MuMuMassUps0->Draw("same");



TCanvas *Kinematic1 = new TCanvas("Kinematic1","Kinematic MuMu",800,500);
   Kinematic1->SetFillColor(0);
   Kinematic1->SetBorderMode(0);
   Kinematic1->SetBorderSize(2);
   Kinematic1->SetFrameBorderMode(0);

Kinematic1->Divide(2,2);
Kinematic1->cd(1);
MuMuMass4->GetXaxis()->SetTitle("#mu#mu mass [GeV]");
MuMuMass4->GetYaxis()->SetTitle("# events / 0.5 GeV");
MuMuMass0->Sumw2();
MuMuMass0->SetLineWidth(2);
MuMuMass0->SetMarkerStyle(20);
MuMuMass1->SetFillColor(30);
MuMuMass2->SetFillColor(30);
MuMuMass2->SetFillStyle(3001);
MuMuMass3->SetFillColor(ci);
MuMuMass4->SetFillColor(38);
MuMuMass4->Draw("hist");
MuMuMass3->Draw("same");
MuMuMass2->Draw("same");
MuMuMass1->Draw("same");
MuMuMass0->Draw("same");

Kinematic1->cd(2);
MuMudeta4->GetXaxis()->SetTitle("#mu#mu #Delta |#eta|");
MuMudeta4->GetYaxis()->SetTitle("# events / 0.5 ");
MuMudeta0->Sumw2();
MuMudeta0->SetLineWidth(2);
MuMudeta0->SetMarkerStyle(20);
MuMudeta1->SetFillColor(30);
MuMudeta2->SetFillColor(30);
MuMudeta2->SetFillStyle(3001);
MuMudeta3->SetFillColor(ci);
MuMudeta4->SetFillColor(38);
MuMudeta4->Draw("hist");
MuMudeta3->Draw("same");
MuMudeta2->Draw("same");
MuMudeta1->Draw("same");
MuMudeta0->Draw("same");

Kinematic1->cd(3);
MuMudpt4->GetXaxis()->SetTitle("#mu#mu |#Delta p_{T}| [GeV]");
MuMudpt4->GetYaxis()->SetTitle("# events / 0.25 GeV ");
MuMudpt0->Sumw2();
MuMudpt0->SetLineWidth(2);
MuMudpt0->SetMarkerStyle(20);
MuMudpt1->SetFillColor(30);
MuMudpt2->SetFillColor(30);
MuMudpt2->SetFillStyle(3001);
MuMudpt3->SetFillColor(ci);
MuMudpt4->SetFillColor(38);
MuMudpt4->Draw("hist");
MuMudpt3->Draw("same");
MuMudpt2->Draw("same");
MuMudpt1->Draw("same");
MuMudpt0->Draw("same");

Kinematic1->cd(4);
MuMudphi4->GetXaxis()->SetTitle("#mu#mu |#Delta #phi / #pi|");
MuMudphi4->GetYaxis()->SetTitle("# events / 0.02 ");
MuMudphi0->Sumw2();
MuMudphi0->SetLineWidth(2);
MuMudphi0->SetMarkerStyle(20);
MuMudphi1->SetFillColor(30);
MuMudphi2->SetFillColor(30);
MuMudphi2->SetFillStyle(3001);
MuMudphi3->SetFillColor(ci);
MuMudphi4->SetFillColor(38);
MuMudphi4->Draw("hist");
MuMudphi3->Draw("same");
MuMudphi2->Draw("same");
MuMudphi1->Draw("same");
MuMudphi0->Draw("same");

TCanvas *CASTOR = new TCanvas("CASTOR","CASTOR",800,550);
   CASTOR->SetFillColor(0);
   CASTOR->SetBorderMode(0);
   CASTOR->SetBorderSize(2);
   CASTOR->SetFrameBorderMode(0);

CASTOR->cd(1);
CastorSumE4->GetXaxis()->SetTitle("Castor Sum Energy [GeV]");
CastorSumE4->GetYaxis()->SetTitle("# events / 5 GeV");
CastorSumE0->Sumw2();
CastorSumE0->SetLineWidth(2);
CastorSumE0->SetMarkerStyle(20);
CastorSumE1->SetFillColor(30);
CastorSumE2->SetFillColor(30);
CastorSumE2->SetFillStyle(3001);
CastorSumE3->SetFillColor(ci);
CastorSumE4->SetFillColor(38);
CastorSumE4->Draw("hist");
CastorSumE3->Draw("same");
CastorSumE2->Draw("same");
CastorSumE1->Draw("same");
CastorSumE0->Draw("same");

TCanvas *ZDC = new TCanvas("ZDC","ZDC",800,550);
   ZDC->SetFillColor(0);
   ZDC->SetBorderMode(0);
   ZDC->SetBorderSize(2);
   ZDC->SetFrameBorderMode(0);

ZDC->Divide(4,2);
ZDC->cd(1);
ZDCemplus4->GetXaxis()->SetTitle("ZDC + em [GeV]");
ZDCemplus4->GetYaxis()->SetTitle("# events / 20 GeV");
ZDCemplus0->Sumw2();
ZDCemplus0->SetLineWidth(2);
ZDCemplus0->SetMarkerStyle(20);
ZDCemplus1->SetFillColor(30);
ZDCemplus2->SetFillColor(30);
ZDCemplus2->SetFillStyle(3001);
ZDCemplus3->SetFillColor(ci);
ZDCemplus4->SetFillColor(38);
ZDCemplus4->Draw("hist");
ZDCemplus3->Draw("same");
ZDCemplus2->Draw("same");
ZDCemplus1->Draw("same");
ZDCemplus0->Draw("same");

ZDC->cd(2);
ZDCemminus4->GetXaxis()->SetTitle("ZDC - em [GeV]");
ZDCemminus4->GetYaxis()->SetTitle("# events / 20 GeV");
ZDCemminus0->Sumw2();
ZDCemminus0->SetLineWidth(2);
ZDCemminus0->SetMarkerStyle(20);
ZDCemminus1->SetFillColor(30);
ZDCemminus2->SetFillColor(30);
ZDCemminus2->SetFillStyle(3001);
ZDCemminus3->SetFillColor(ci);
ZDCemminus4->SetFillColor(38);
ZDCemminus4->Draw("hist");
ZDCemminus3->Draw("same");
ZDCemminus2->Draw("same");
ZDCemminus1->Draw("same");
ZDCemminus0->Draw("same");

ZDC->cd(3);
ZDChadplus4->GetXaxis()->SetTitle("ZDC + had [GeV]");
ZDChadplus4->GetYaxis()->SetTitle("# events / 20 GeV");
ZDChadplus0->Sumw2();
ZDChadplus0->SetLineWidth(2);
ZDChadplus0->SetMarkerStyle(20);
ZDChadplus1->SetFillColor(30);
ZDChadplus2->SetFillColor(30);
ZDChadplus2->SetFillStyle(3001);
ZDChadplus3->SetFillColor(ci);
ZDChadplus4->SetFillColor(38);
ZDChadplus4->Draw("hist");
ZDChadplus3->Draw("same");
ZDChadplus2->Draw("same");
ZDChadplus1->Draw("same");
ZDChadplus0->Draw("same");

ZDC->cd(4);
ZDChadminus4->GetXaxis()->SetTitle("ZDC - had [GeV]");
ZDChadminus4->GetYaxis()->SetTitle("# events / 20 GeV");
ZDChadminus0->Sumw2();
ZDChadminus0->SetLineWidth(2);
ZDChadminus0->SetMarkerStyle(20);
ZDChadminus1->SetFillColor(30);
ZDChadminus2->SetFillColor(30);
ZDChadminus2->SetFillStyle(3001);
ZDChadminus3->SetFillColor(ci);
ZDChadminus4->SetFillColor(38);
ZDChadminus4->Draw("hist");
ZDChadminus3->Draw("same");
ZDChadminus2->Draw("same");
ZDChadminus1->Draw("same");
ZDChadminus0->Draw("same");

ZDC->cd(5);
ZDCtime4->GetXaxis()->SetTitle("ZDC hit time [???]");
ZDCtime4->GetYaxis()->SetTitle("# events / 1 ???");
ZDCtime0->Sumw2();
ZDCtime0->SetLineWidth(2);
ZDCtime0->SetMarkerStyle(20);
ZDCtime1->SetFillColor(30);
ZDCtime2->SetFillColor(30);
ZDCtime2->SetFillStyle(3001);
ZDCtime3->SetFillColor(ci);
ZDCtime4->SetFillColor(38);
ZDCtime4->Draw("hist");
ZDCtime3->Draw("same");
ZDCtime2->Draw("same");
ZDCtime1->Draw("same");
ZDCtime0->Draw("same");

ZDC->cd(7);
ZDCenergyEM4->GetXaxis()->SetTitle("ZDC hit EM energy [GeV ??]");
ZDCenergyEM4->GetYaxis()->SetTitle("# events / 30 GeV ???");
ZDCenergyEM0->Sumw2();
ZDCenergyEM0->SetLineWidth(2);
ZDCenergyEM0->SetMarkerStyle(20);
ZDCenergyEM1->SetFillColor(30);
ZDCenergyEM2->SetFillColor(30);
ZDCenergyEM2->SetFillStyle(3001);
ZDCenergyEM3->SetFillColor(ci);
ZDCenergyEM4->SetFillColor(38);
ZDCenergyEM4->Draw("hist");
ZDCenergyEM3->Draw("same");
ZDCenergyEM2->Draw("same");
ZDCenergyEM1->Draw("same");
ZDCenergyEM0->Draw("same");

ZDC->cd(8);
ZDCenergyHAD4->GetXaxis()->SetTitle("ZDC hit HAD energy [GeV ??]");
ZDCenergyHAD4->GetYaxis()->SetTitle("# events / 300 GeV ???");
ZDCenergyHAD0->Sumw2();
ZDCenergyHAD0->SetLineWidth(2);
ZDCenergyHAD0->SetMarkerStyle(20);
ZDCenergyHAD1->SetFillColor(30);
ZDCenergyHAD2->SetFillColor(30);
ZDCenergyHAD2->SetFillStyle(3001);
ZDCenergyHAD3->SetFillColor(ci);
ZDCenergyHAD4->SetFillColor(38);
ZDCenergyHAD4->Draw("hist");
ZDCenergyHAD3->Draw("same");
ZDCenergyHAD2->Draw("same");
ZDCenergyHAD1->Draw("same");
ZDCenergyHAD0->Draw("same");

}
cout << "END" << endl;   
}

