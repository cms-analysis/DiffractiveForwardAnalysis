void TestGhost_XXX_NUM_XXX()
{
gROOT->Reset();
gStyle->SetPalette(1);

#define pi 3.1413565359
//definition des fichiers + Tree
//  TFile *f0 = new TFile("zerobias_merge1.root"); // 
//  TTree *t0 = f0->Get("ntp1");

  TChain *t0 = new TChain("ntp1");
//  t0->Add("/storage/data/cms/store/user/schul/384_finalZB/final_ZBCom10_*.root");
  t0->Add("/storage/data/cms/store/user/schul/ZeroBias_prescale/prescale_ZBRunA_*.root");
  t0->Add("/storage/data/cms/store/user/schul/ZeroBias_prescale/prescale_MBRunA_*.root");
  t0->Add("/storage/data/cms/store/user/schul/ZeroBias_prescale/prescale_MBRunB_*.root");

// definitions des Trees pour la sauvegarde
  TH1F* hnVtx_loose_0 = new TH1F("nVtx_cond0","",21.,-1.,20.);

//Related to the closest track:
  TH1F* hDist_D0 = new TH1F("Dist_D0","",200,0.,0.2);
  TH1F* hDist_T0 = new TH1F("Dist_T0","",200,0.,0.2);
  TH1F* hDist_Z0 = new TH1F("Dist_Z0","",200,0.,0.2);

  TH1F* hpT  = new TH1F("pT0","",200,0.,20.);
  TH1F* hpZ  = new TH1F("pZ0","",200,-50.,50.);
  TH1F* heta = new TH1F("eta0","",120,-3.,3.);
  TH1F* hphi = new TH1F("phi0","",88,-1.1,1.1);

  TH1F* hchi2 = new TH1F("chi2_0","",200,0.,200.);
  TH1F* hProb = new TH1F("Prob_0","",100,0.,1.);
  TH1F* hPure = new TH1F("Pure_0","",5,-0.5,2.0);
  TH1F* hHits = new TH1F("Hits_0","",35,0.,35.);

// definitions des # d'entrée

//definition des variables
  Int_t techBit1[1][128], techBit2[1][128], techBit0[1][128];
// MuonID
// GenTrack
// RecoTrack
  Int_t var_nTrack1[1], var_nTrack2[1], var_nTrack0[1];
  Int_t var_nTrackQual1[1], var_nTrackQual2[1], var_nTrackQual0[1];
  Double_t var_TrackPt1[2000], var_TrackPt2[2000], var_TrackPt0[2000], var_TrackPz0[2000];
  Double_t var_TrackD1[2000],var_TrackD2[2000], var_TrackD0[2000];
  Double_t var_TrackQuality1[2000],var_TrackQuality2[2000], var_TrackQuality0[2000];
  Double_t var_TrackZ0[2000],  var_TrackX0[2000], var_TrackY0[2000];
  Double_t var_TrackChi2[2000], var_TrackNdof[2000];
  Int_t var_TrackHits[2000];
  Double_t var_TrackEta[2000], var_TrackPhi[2000];

// Prim Vtx
  Int_t var_nvtx1[1],var_nvtx2[1], var_nvtx0[1];
  Int_t var_vtxTrack1[30],var_vtxTrack2[30], var_vtxTrack0[30];
  Double_t var_vtxZ1[30],var_vtxZ2[30], var_vtxZ0[30];
  Double_t var_vtxX1[30],var_vtxX2[30], var_vtxX0[30];
  Double_t var_vtxY1[30],var_vtxY2[30], var_vtxY0[30];
  Double_t var_vertexChi2_1[30],var_vertexChi2_2[30], var_vertexChi2_0[30];
  Double_t var_vertexNdf1[30], var_vertexNdf2[30], var_vertexNdf0[30];

// CaloTowers 
  Int_t var_ncalo1[1], var_ncalo2[1], var_ncalo0[1];
  Int_t var_tower1[1], var_tower2[1], var_tower0[1];
  Int_t var_caloId1[2000],var_caloId2[2000], var_caloId0[2000];
  Double_t var_caloEn1[2000],var_caloEn2[2000], var_caloEn0[2000];
  Double_t var_caloEmE1[2000],var_caloEmE2[2000], var_caloEmE0[2000];
  Double_t var_caloHadE1[2000],var_caloHadE2[2000], var_caloHadE0[2000];
  Double_t var_caloTime1[2000], var_caloTime2[2000], var_caloTime0[2000];
  Double_t var_caloEta0[2000], var_caloX0[2000], var_caloY0[2000], var_caloZ0[2000], var_caloPhi0[2000];
// Bunch crossing
  Int_t var_bx0[1], var_run0[1], var_ls0[1];
// MuMu kinematics
// ZDC
  Int_t var_nZDC1[1], var_nZDC2[1], var_nZDC0[1];
  Int_t var_zdcsection1[5000],var_zdcsection2[5000],var_zdcsection0[5000];
  Int_t var_zdcside1[5000],var_zdcside2[5000],var_zdcside0[5000];
  Double_t var_zdcE1[5000], var_zdcE2[5000], var_zdcE0[5000];
  Double_t var_zdcEmMinus1[1], var_zdcEmMinus2[1],var_zdcEmMinus0[1], var_zdcHadMinus1[1], var_zdcHadMinus2[1], var_zdcHadMinus0[1];
  Double_t var_zdcEmPlus1[1], var_zdcEmPlus2[1], var_zdcEmPlus0[1],var_zdcHadPlus1[1], var_zdcHadPlus2[1], var_zdcHadPlus0[1];
  Double_t var_zdcTime1[5000], var_zdcTime2[5000], var_zdcTime0[5000];

//Castor
  Int_t var_nCastor1[1],var_nCastor2[1], var_nCastor0[1];
  Double_t var_CastorE1[1000],  var_CastorE2[1000], var_CastorE0[1000];
  Double_t var_CastorEta1[1000],  var_CastorEta2[1000], var_CastorEta0[1000];
  Double_t var_CastorPhi1[1000],  var_CastorPhi2[1000], var_CastorPhi0[1000];
  Double_t var_CastorRecHit1[1], var_CastorRecHit2[1], var_CastorRecHit0[1];

  t0->SetBranchAddress("L1TechnicalTriggers",techBit0);

  t0->SetBranchAddress("nTrackCand",var_nTrack0);
  t0->SetBranchAddress("nQualityTrackCand",var_nTrackQual0);
  t0->SetBranchAddress("TrackCand_pt",var_TrackPt0);
  t0->SetBranchAddress("TrackCand_pz",var_TrackPz0);
  t0->SetBranchAddress("TrackCand_purity",var_TrackQuality0);
  t0->SetBranchAddress("TrackCand_z",var_TrackZ0);
  t0->SetBranchAddress("TrackCand_x",var_TrackX0);
  t0->SetBranchAddress("TrackCand_y",var_TrackY0);
  t0->SetBranchAddress("TrackCand_chi2",var_TrackChi2);
  t0->SetBranchAddress("TrackCand_ndof",var_TrackNdof);
  t0->SetBranchAddress("TrackCand_nhits",var_TrackHits);
  t0->SetBranchAddress("TrackCand_eta",var_TrackEta);
  t0->SetBranchAddress("TrackCand_phi",var_TrackPhi);

  t0->SetBranchAddress("nVertexCand",var_nvtx0);
  t0->SetBranchAddress("VertexCand_z",var_vtxZ0);
  t0->SetBranchAddress("VertexCand_x",var_vtxX0);
  t0->SetBranchAddress("VertexCand_y",var_vtxY0);
  t0->SetBranchAddress("VertexCand_chi2",var_vertexChi2_0);
  t0->SetBranchAddress("VertexCand_ndof",var_vertexNdf0);
  t0->SetBranchAddress("VertexCand_tracks",var_vtxTrack0);

  t0->SetBranchAddress("nCaloCand",var_ncalo0);
  t0->SetBranchAddress("CaloTower_ID",var_caloId0);
  t0->SetBranchAddress("CaloTower_e",var_caloEn0);
  t0->SetBranchAddress("CaloTower_emE",var_caloEmE0);
  t0->SetBranchAddress("CaloTower_hadE",var_caloHadE0);
  t0->SetBranchAddress("CaloTower_t",var_caloTime0);
  t0->SetBranchAddress("CaloTower_eta",var_caloEta0);
  t0->SetBranchAddress("CaloTower_phi",var_caloPhi0);

  t0->SetBranchAddress("BX",var_bx0);
  t0->SetBranchAddress("Run",var_run0);
  t0->SetBranchAddress("LumiSection",var_ls0);

  t0->SetBranchAddress("nZDChitCand",var_nZDC0);
  t0->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus0);
  t0->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus0);
  t0->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus0);
  t0->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus0);
  t0->SetBranchAddress("ZDChit_time",var_zdcTime0);
  t0->SetBranchAddress("ZDChit_energy",var_zdcE0);
  t0->SetBranchAddress("ZDChit_section",var_zdcsection0);
  t0->SetBranchAddress("ZDChit_side",var_zdcside0);

  t0->SetBranchAddress("nCastorTowerCand",var_nCastor0);
  t0->SetBranchAddress("CastorTower_e",var_CastorE0);
  t0->SetBranchAddress("CastorTower_eta",var_CastorEta0);
  t0->SetBranchAddress("CastorTower_phi",var_CastorPhi0);
  t0->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit0);

  t0->SetBranchAddress("L1TechnicalTriggers",techBit0);
  const int NUM0 = t0->GetEntries();
  cout<<"NUM0="<<NUM0<<endl;
 
//Open txt files
/*  string filename = "out_BeamSpot"; //input
  ifstream infile(filename.c_str());
  if (! infile.is_open()) { cout << "\t ERROR: I can not open \"" << filename << "\"" << endl; return; }
  string temp_string, temp;
  istringstream curstring;*/

  const float fac_0=1.;

      double nTestVeto_loose[20], nTestVeto_tight[20];
      double nPassVeto_loose[20], nPassVeto_tight[20];
      for(int j=0; j<20; j++){nPassVeto_loose[j]=0; nPassVeto_tight[j]=0; nTestVeto_loose[j]=0; nTestVeto_tight[j]=0;}

  int filter0Gen(0);
  int filter0Track(0);
  int filter0Events(0);
  int n_cond0(0);
  int n_cond1(0);
  int n_cond2(0);
  for(Int_t i = XXX_START_XXX;i < XXX_STOP_XXX;i++){
    t0->GetEntry(i);
    if(var_run0[0]>=140058){
      int nPrimVtx = var_nvtx0[0];
      int nTrack=var_nTrack0[0];
      int nCalo=var_ncalo0[0];
      double nValidVtx_loose(0);
      if(!(techBit0[0][36]||techBit0[0][37]||techBit0[0][38]||techBit0[0][39]) && !((techBit0[0][42]&&!(techBit0[0][43]))||(techBit0[0][43]&&!(techBit0[0][42]))) ){
	if(techBit0[0][0]){ 
//		cout<<"beam-crossing"<<endl;
//		cout<<"----------"<<var_run0[0]<<endl;
		for(Int_t j=0; j<nPrimVtx; j++){if(var_vtxTrack0[j]>0 && TMath::Prob(var_vertexChi2_0[j],var_vertexNdf0[j]+0.5)>0.001) nValidVtx_loose++;}
		if(nValidVtx_loose==0){
		double vx(0),vy(0),vz(0),weight(0);
		if(var_run0[0]==140058){weight=(1./110975)*(3393.912/35322028.2939999923);vx=0.100794;vy= 0.0143116;vz=gRandom->Gaus(0.271343,6.06352);}
		else if(var_run0[0]==140059){weight=(1./1156814)*(23253.025/35322028.2939999923);vx=0.100335;vy= 0.014676;vz=gRandom->Gaus(0.0995124,6.42079);}
		else if(var_run0[0]==140124){weight=(1./1698837)*(44815.019/35322028.2939999923);vx=0.095623;vy= 0.015889;vz=gRandom->Gaus(-0.152702,6.36746);}
		else if(var_run0[0]==140158){weight=(1./391895)*(9983.773/35322028.2939999923);vx=0.0976755;vy= 0.0155846;vz=gRandom->Gaus(-0.587424,6.15816);}
		else if(var_run0[0]==140159){weight=(1./375558)*(8304.239/35322028.2939999923);vx=0.0969652;vy= 0.0166257;vz=gRandom->Gaus(0.0376294,6.06082);}
		else if(var_run0[0]==140160){weight=(1./216277)*(4157.372/35322028.2939999923);vx=0.0980627;vy= 0.0177718;vz=gRandom->Gaus(-0.49599,6.38909);}
		else if(var_run0[0]==140331){weight=(1./496714)*(13871.020/35322028.2939999923);vx=0.0966409;vy= 0.0135332;vz=gRandom->Gaus(-0.303373,6.31805);}
		else if(var_run0[0]==140359){weight=(1./36699)*(812.983/35322028.2939999923);vx=0.0940866;vy= 0.0146067;vz=gRandom->Gaus(-0.606706,5.95362);}
		else if(var_run0[0]==140361){weight=(1./149996)*(3216.423/35322028.2939999923);vx=0.0948746;vy= 0.0126356;vz=gRandom->Gaus(0.0671472,6.31849);}
		else if(var_run0[0]==140362){weight=(1./276827)*(5502.570/35322028.2939999923);vx=0.094178;vy= 0.0116153;vz=gRandom->Gaus(0.080781,6.46025);}
		else if(var_run0[0]==140379){weight=(1./116461)*(4482.460/35322028.2939999923);vx=0.0950191;vy= 0.0148368;vz=gRandom->Gaus(0.0353279,6.02474);}
		else if(var_run0[0]==140381){weight=(1./20514)*(749.960/35322028.2939999923);vx=0.0945141;vy= 0.0137429;vz=gRandom->Gaus(-0.552448,7.15243);}
		else if(var_run0[0]==140382){weight=(1./204860)*(7089.466/35322028.2939999923);vx=0.101838;vy= 0.0145831;vz=gRandom->Gaus(-0.138335,6.63679);}
		else if(var_run0[0]==140383){weight=(1./466126)*(13684.672/35322028.2939999923);vx=0.0938265;vy= 0.0144248;vz=gRandom->Gaus(0.305865,6.27479);}
		else if(var_run0[0]==140385){weight=(1./137196)*(3598.491/35322028.2939999923);vx=0.092982;vy= 0.0103742;vz=gRandom->Gaus(-0.434873,5.98944);}
		else if(var_run0[0]==140386){weight=(1./22320)*(569.441/35322028.2939999923);vx=0.091979;vy= 0.0139829;vz=gRandom->Gaus(-0.850734,6.33796);}
		else if(var_run0[0]==140387){weight=(1./128271)*(3110.306/35322028.2939999923);vx=0.0947584;vy= 0.0133499;vz=gRandom->Gaus(-0.421871,6.3771);}
		else if(var_run0[0]==140388){weight=(1./112853)*(2651.648/35322028.2939999923);vx=0.0926958;vy= 0.0137019;vz=gRandom->Gaus(-0.684303,6.3918);}
		else if(var_run0[0]==140399){weight=(1./703881)*(13139.036/35322028.2939999923);vx=0.0932278;vy= 0.0132279;vz=gRandom->Gaus(-0.106708,6.51808);}
		else if(var_run0[0]==140401){weight=(1./267918)*(4392.874/35322028.2939999923);vx=0.0931317;vy= 0.0145009;vz=gRandom->Gaus(0.186389,6.12158);}
		else if(var_run0[0]==141880){weight=(1./197073)*(4609.722/35322028.2939999923);vx=0.0945815;vy= 0.0103382;vz=gRandom->Gaus(0.0246312,6.22766);}
		else if(var_run0[0]==141881){weight=(1./696161)*(14363.950/35322028.2939999923);vx=0.095274;vy= 0.0111391;vz=gRandom->Gaus(-0.199211,6.47942);}
		else if(var_run0[0]==141956){weight=(1./221428)*(19230.703/35322028.2939999923);vx=0.0953798;vy= 0.0125119;vz=gRandom->Gaus(-0.0746152,6.42123);}
		else if(var_run0[0]==141957){weight=(1./55786)*(4383.057/35322028.2939999923);vx=0.0945655;vy= 0.0142757;vz=gRandom->Gaus(0.00655455,6.55368);}
		else if(var_run0[0]==141958){weight=(1./22321)*(1705.218/35322028.2939999923);vx=0.0929087;vy= 0.0126065;vz=gRandom->Gaus(-0.149173,6.28642);}
		else if(var_run0[0]==141959){weight=(1./53314)*(4026.959/35322028.2939999923);vx=0.0938471;vy= 0.0115851;vz=gRandom->Gaus(-0.175745,6.45066);}
		else if(var_run0[0]==141960){weight=(1./193785)*(12995.040/35322028.2939999923);vx=0.0944724;vy= 0.0129941;vz=gRandom->Gaus(-0.164615,6.46335);}
		else if(var_run0[0]==141961){weight=(1./28211)*(1741.197/35322028.2939999923);vx=0.0980082;vy= 0.0150895;vz=gRandom->Gaus(0.608491,5.91946);}
		else if(var_run0[0]==142035){weight=(1./292407)*(21984.860/35322028.2939999923);vx=0.0937506;vy= 0.0126625;vz=gRandom->Gaus(-0.0332531,6.38996);}
		else if(var_run0[0]==142036){weight=(1./185982)*(12289.327/35322028.2939999923);vx=0.0932561;vy= 0.0113626;vz=gRandom->Gaus(-0.163195,6.35545);}
		else if(var_run0[0]==142038){weight=(1./565676)*(33081.827/35322028.2939999923);vx=0.094352;vy= 0.0121395;vz=gRandom->Gaus(0.0304028,6.45642);}
		else if(var_run0[0]==142039){weight=(1./4055)*(207.039/35322028.2939999923);vx=0.0963541;vy= 0.00278066;vz=gRandom->Gaus(-1.49764,6.2135);}
		else if(var_run0[0]==142040){weight=(1./234628)*(12301.379/35322028.2939999923);vx=0.0926035;vy= 0.0115514;vz=gRandom->Gaus(-0.633006,6.5702);}
		else if(var_run0[0]==142076){weight=(1./98266)*(6503.550/35322028.2939999923);vx=0.0935682;vy= 0.0126095;vz=gRandom->Gaus(-0.188606,6.16845);}
		else if(var_run0[0]==142128){weight=(1./8869)*(865.923/35322028.2939999923);vx=0.0889196;vy= 0.0105045;vz=gRandom->Gaus(-0.390518,6.1126);}
		else if(var_run0[0]==142129){weight=(1./28954)*(2736.993/35322028.2939999923);vx=0.0925292;vy= 0.0119757;vz=gRandom->Gaus(-0.278862,6.28139);}
		else if(var_run0[0]==142130){weight=(1./545443)*(42616.540/35322028.2939999923);vx=0.0921308;vy= 0.0125107;vz=gRandom->Gaus(-0.128887,6.45175);}
		else if(var_run0[0]==142132){weight=(1./154954)*(10375.682/35322028.2939999923);vx=0.0916467;vy= 0.0118682;vz=gRandom->Gaus(-0.230859,6.23377);}
		else if(var_run0[0]==142135){weight=(1./102074)*(6547.628/35322028.2939999923);vx=0.092053;vy= 0.0114144;vz=gRandom->Gaus(-0.337647,6.62729);}
		else if(var_run0[0]==142136){weight=(1./151721)*(9458.010/35322028.2939999923);vx=0.0926426;vy= 0.0109235;vz=gRandom->Gaus(-0.310499,6.2499);}
		else if(var_run0[0]==142137){weight=(1./432783)*(24676.320/35322028.2939999923);vx=0.0923496;vy= 0.0120018;vz=gRandom->Gaus(-0.384746,6.49448);}
		else if(var_run0[0]==142189){weight=(1./119262)*(11612.672/35322028.2939999923);vx=0.0954219;vy= 0.0118383;vz=gRandom->Gaus(-0.782422,6.66343);}
		else if(var_run0[0]==142191){weight=(1./116751)*(10489.649/35322028.2939999923);vx=0.0948962;vy= 0.0103308;vz=gRandom->Gaus(-0.655677,6.57426);}
		else if(var_run0[0]==142264){weight=(1./25505)*(3532.925/35322028.2939999923);vx=0.0908284;vy= 0.0098164;vz=gRandom->Gaus(-0.00313011,6.35974);}
		else if(var_run0[0]==142265){weight=(1./84885)*(10816.338/35322028.2939999923);vx=0.0912857;vy= 0.0120297;vz=gRandom->Gaus(-0.363364,6.21843);}
		else if(var_run0[0]==142303){weight=(1./41931)*(5915.102/35322028.2939999923);vx=0.0936276;vy= 0.0123935;vz=gRandom->Gaus(-0.19073,6.3332);}
		else if(var_run0[0]==142304){weight=(1./4298)*(612.104/35322028.2939999923);vx=0.0953863;vy= 0.0124424;vz=gRandom->Gaus(0.805653,6.46557);}
		else if(var_run0[0]==142305){weight=(1./121260)*(15825.586/35322028.2939999923);vx=0.0930993;vy= 0.0123919;vz=gRandom->Gaus(-0.259715,6.18923);}
		else if(var_run0[0]==142308){weight=(1./19475)*(2345.748/35322028.2939999923);vx=0.0936005;vy= 0.0126497;vz=gRandom->Gaus(0.129302,6.52169);}
		else if(var_run0[0]==142309){weight=(1./50329)*(5847.004/35322028.2939999923);vx=0.094416;vy= 0.0119167;vz=gRandom->Gaus(-0.137821,6.18922);}
		else if(var_run0[0]==142311){weight=(1./301348)*(30114.738/35322028.2939999923);vx=0.0935093;vy= 0.0112015;vz=gRandom->Gaus(0.0334105,6.43299);}
		else if(var_run0[0]==142312){weight=(1./154097)*(13210.534/35322028.2939999923);vx=0.0929723;vy= 0.0114759;vz=gRandom->Gaus(-0.409077,6.484);}
		else if(var_run0[0]==142313){weight=(1./70715)*(5779.889/35322028.2939999923);vx=0.0948427;vy= 0.00903676;vz=gRandom->Gaus(-0.120377,6.23579);}
		else if(var_run0[0]==142413){weight=(1./41300)*(5491.433/35322028.2939999923);vx=0.0952755;vy= 0.0124298;vz=gRandom->Gaus(-0.246799,6.15094);}
		else if(var_run0[0]==142414){weight=(1./64369)*(7693.686/35322028.2939999923);vx=0.0927713;vy= 0.0139964;vz=gRandom->Gaus(-0.537765,6.66217);}
		else if(var_run0[0]==142419){weight=(1./113973)*(10875.088/35322028.2939999923);vx=0.0948953;vy= 0.0116192;vz=gRandom->Gaus(0.216078,6.65546);}
		else if(var_run0[0]==142422){weight=(1./628312)*(43168.626/35322028.2939999923);vx=0.094908;vy= 0.0130058;vz=gRandom->Gaus(-0.254085,6.48445);}
		else if(var_run0[0]==142513){weight=(1./45708)*(8432.289/35322028.2939999923);vx=0.0936079;vy= 0.0128486;vz=gRandom->Gaus(0.0553905,6.29823);}
		else if(var_run0[0]==142514){weight=(1./14977)*(2596.681/35322028.2939999923);vx=0.0928844;vy= 0.0117103;vz=gRandom->Gaus(0.471611,6.16774);}
		else if(var_run0[0]==142523){weight=(1./26419)*(4042.498/35322028.2939999923);vx=0.0949826;vy= 0.0135065;vz=gRandom->Gaus(-0.157981,6.08754);}
		else if(var_run0[0]==142524){weight=(1./149889)*(19399.445/35322028.2939999923);vx=0.0955073;vy= 0.013181;vz=gRandom->Gaus(0.0191425,6.16217);}
		else if(var_run0[0]==142525){weight=(1./84985)*(9230.657/35322028.2939999923);vx=0.0960522;vy= 0.0123776;vz=gRandom->Gaus(-0.327897,6.2091);}
		else if(var_run0[0]==142528){weight=(1./462617)*(41822.189/35322028.2939999923);vx=0.0964386;vy= 0.0132233;vz=gRandom->Gaus(-0.0868653,6.37082);}
		else if(var_run0[0]==142530){weight=(1./38479)*(3040.470/35322028.2939999923);vx=0.093721;vy= 0.0132565;vz=gRandom->Gaus(0.0801791,6.68744);}
		else if(var_run0[0]==142535){weight=(1./52981)*(4220.037/35322028.2939999923);vx=0.0960788;vy= 0.0107834;vz=gRandom->Gaus(-0.630023,6.30529);}
		else if(var_run0[0]==142537){weight=(1./47616)*(3707.233/35322028.2939999923);vx=0.0960362;vy= 0.0110998;vz=gRandom->Gaus(-0.198908,5.99704);}
		else if(var_run0[0]==142557){weight=(1./54776)*(8162.407/35322028.2939999923);vx=0.0972623;vy= 0.0142449;vz=gRandom->Gaus(0.0213277,6.11328);}
		else if(var_run0[0]==142558){weight=(1./120233)*(15924.613/35322028.2939999923);vx=0.0974421;vy= 0.0130626;vz=gRandom->Gaus(-0.226053,6.29172);}
		else if(var_run0[0]==142657){weight=(1./425)*(85.037/35322028.2939999923);vx=0.0885586;vy= 0.0294256;vz=gRandom->Gaus(-0.604557,3.23314);}
		else if(var_run0[0]==142658){weight=(1./1963)*(371.862/35322028.2939999923);vx=0.0924293;vy= 0.00755789;vz=gRandom->Gaus(-0.184352,5.82713);}
		else if(var_run0[0]==142659){weight=(1./28793)*(4967.979/35322028.2939999923);vx=0.0886328;vy= 0.0102111;vz=gRandom->Gaus(-0.365842,6.83124);}
		else if(var_run0[0]==142660){weight=(1./6993)*(1231.645/35322028.2939999923);vx=0.0883789;vy= 0.0122418;vz=gRandom->Gaus(-0.191347,6.07437);}
		else if(var_run0[0]==142661){weight=(1./7272)*(1270.152/35322028.2939999923);vx=0.0931001;vy= 0.0118804;vz=gRandom->Gaus(-0.437016,6.88346);}
		else if(var_run0[0]==142662){weight=(1./85315)*(13940.020/35322028.2939999923);vx=0.0911049;vy= 0.0155607;vz=gRandom->Gaus(-0.102537,6.77731);}
		else if(var_run0[0]==142663){weight=(1./77005)*(11595.491/35322028.2939999923);vx=0.089634;vy= 0.0251548;vz=gRandom->Gaus(0.0778907,6.7571);}
		else if(var_run0[0]==142664){weight=(1./36291)*(5184.042/35322028.2939999923);vx=0.0911891;vy= 0.0246174;vz=gRandom->Gaus(0.178342,7.0693);}
		else if(var_run0[0]==142928){weight=(1./286958)*(41225.161/35322028.2939999923);vx=0.0898088;vy= 0.0277331;vz=gRandom->Gaus(-0.21254,6.46752);}
		else if(var_run0[0]==142933){weight=(1./448940)*(44031.712/35322028.2939999923);vx=0.0897305;vy= 0.0278015;vz=gRandom->Gaus(-0.0292006,6.39151);}
		else if(var_run0[0]==142935){weight=(1./20705)*(1758.902/35322028.2939999923);vx=0.0889369;vy= 0.0274965;vz=gRandom->Gaus(-0.191199,6.78442);}
		else if(var_run0[0]==142936){weight=(1./6916)*(577.386/35322028.2939999923);vx=0.0915178;vy= 0.0266161;vz=gRandom->Gaus(-0.95817,6.59995);}
		else if(var_run0[0]==142953){weight=(1./111397)*(21007.826/35322028.2939999923);vx=0.0909008;vy= 0.0271429;vz=gRandom->Gaus(-0.203396,6.36024);}
		else if(var_run0[0]==142954){weight=(1./109030)*(18137.464/35322028.2939999923);vx=0.0909349;vy= 0.0267102;vz=gRandom->Gaus(0.0981885,6.3887);}
		else if(var_run0[0]==142970){weight=(1./64460)*(12565.900/35322028.2939999923);vx=0.089288;vy= 0.0232477;vz=gRandom->Gaus(-0.304972,6.68222);}
		else if(var_run0[0]==142971){weight=(1./788948)*(107097.262/35322028.2939999923);vx=0.0907004;vy= 0.0255772;vz=gRandom->Gaus(-0.120557,6.6521);}
		else if(var_run0[0]==143004){weight=(1./38072)*(6856.075/35322028.2939999923);vx=0.0910295;vy= 0.0272586;vz=gRandom->Gaus(-0.075168,6.32296);}
		else if(var_run0[0]==143005){weight=(1./108515)*(17889.789/35322028.2939999923);vx=0.0905392;vy= 0.0270044;vz=gRandom->Gaus(0.0446868,6.48798);}
		else if(var_run0[0]==143006){weight=(1./53272)*(8059.309/35322028.2939999923);vx=0.0908894;vy= 0.0276938;vz=gRandom->Gaus(0.157866,6.73279);}
		else if(var_run0[0]==143007){weight=(1./402306)*(49434.208/35322028.2939999923);vx=0.0901827;vy= 0.0271159;vz=gRandom->Gaus(-0.237678,6.56121);}
		else if(var_run0[0]==143008){weight=(1./20340)*(2034.899/35322028.2939999923);vx=0.0891607;vy= 0.0271899;vz=gRandom->Gaus(-0.107623,6.24831);}
		else if(var_run0[0]==143179){weight=(1./25946)*(3978.788/35322028.2939999923);vx=0.091262;vy= 0.028573;vz=gRandom->Gaus(-0.185253,6.35637);}
		else if(var_run0[0]==143181){weight=(1./492189)*(58549.114/35322028.2939999923);vx=0.0910366;vy= 0.026101;vz=gRandom->Gaus(-0.214549,6.36783);}
		else if(var_run0[0]==143187){weight=(1./155700)*(14349.063/35322028.2939999923);vx=0.0911821;vy= 0.0257043;vz=gRandom->Gaus(-0.26801,6.36053);}
		else if(var_run0[0]==143191){weight=(1./80176)*(6869.082/35322028.2939999923);vx=0.0917604;vy= 0.0260234;vz=gRandom->Gaus(-0.271733,6.55902);}
		else if(var_run0[0]==143192){weight=(1./22415)*(1859.083/35322028.2939999923);vx=0.0889398;vy= 0.0269094;vz=gRandom->Gaus(-0.0534999,6.63494);}
		else if(var_run0[0]==143193){weight=(1./35973)*(2930.126/35322028.2939999923);vx=0.0897632;vy= 0.0270503;vz=gRandom->Gaus(-0.0998894,6.41284);}
		else if(var_run0[0]==143318){weight=(1./23640)*(7079.072/35322028.2939999923);vx=0.0914484;vy= 0.0263239;vz=gRandom->Gaus(0.255529,6.5057);}
		else if(var_run0[0]==143319){weight=(1./22659)*(7355.080/35322028.2939999923);vx=0.0926588;vy= 0.0267072;vz=gRandom->Gaus(-0.0883608,6.41922);}
		else if(var_run0[0]==143320){weight=(1./78468)*(24112.801/35322028.2939999923);vx=0.0899219;vy= 0.0253984;vz=gRandom->Gaus(-0.246676,6.28565);}
		else if(var_run0[0]==143321){weight=(1./20980)*(6118.496/35322028.2939999923);vx=0.0906771;vy= 0.026495;vz=gRandom->Gaus(-0.361193,6.562);}
		else if(var_run0[0]==143322){weight=(1./67970)*(18110.010/35322028.2939999923);vx=0.0898169;vy= 0.0264218;vz=gRandom->Gaus(-0.321438,6.43413);}
		else if(var_run0[0]==143323){weight=(1./272123)*(62803.158/35322028.2939999923);vx=0.0900378;vy= 0.0263641;vz=gRandom->Gaus(-0.128706,6.51767);}
		else if(var_run0[0]==143326){weight=(1./125901)*(25240.384/35322028.2939999923);vx=0.0896332;vy= 0.0253725;vz=gRandom->Gaus(-0.674726,6.63082);}
		else if(var_run0[0]==143327){weight=(1./109157)*(20380.491/35322028.2939999923);vx=0.0906081;vy= 0.0271599;vz=gRandom->Gaus(-0.279194,6.54791);}
		else if(var_run0[0]==143328){weight=(1./195211)*(33413.441/35322028.2939999923);vx=0.0907095;vy= 0.0257698;vz=gRandom->Gaus(-0.404016,6.53572);}
		else if(var_run0[0]==143657){weight=(1./834587)*(269224.618/35322028.2939999923);vx=0.0905717;vy= 0.0270814;vz=gRandom->Gaus(-0.177333,6.57029);}
		else if(var_run0[0]==143665){weight=(1./71177)*(15576.111/35322028.2939999923);vx=0.0896075;vy= 0.0265908;vz=gRandom->Gaus(-0.223942,6.72243);}
		else if(var_run0[0]==143727){weight=(1./154412)*(62583.918/35322028.2939999923);vx=0.0912077;vy= 0.0268628;vz=gRandom->Gaus(-0.310465,6.3967);}
		else if(var_run0[0]==143731){weight=(1./10044)*(3637.397/35322028.2939999923);vx=0.0905155;vy= 0.0264595;vz=gRandom->Gaus(0.257504,6.68583);}
		else if(var_run0[0]==143827){weight=(1./279828)*(132241.198/35322028.2939999923);vx=0.089576;vy= 0.0277203;vz=gRandom->Gaus(-0.241826,6.52543);}
		else if(var_run0[0]==143833){weight=(1./306686)*(129787.949/35322028.2939999923);vx=0.0900497;vy= 0.0275122;vz=gRandom->Gaus(-0.124198,6.77496);}
		else if(var_run0[0]==143953){weight=(1./177885)*(106768.206/35322028.2939999923);vx=0.0906401;vy= 0.0263002;vz=gRandom->Gaus(-0.308051,6.53142);}
		else if(var_run0[0]==143955){weight=(1./25795)*(13033.752/35322028.2939999923);vx=0.0884157;vy= 0.0259701;vz=gRandom->Gaus(-0.555993,6.69627);}
		else if(var_run0[0]==143956){weight=(1./15241)*(7573.906/35322028.2939999923);vx=0.0915183;vy= 0.0250372;vz=gRandom->Gaus(-0.221492,6.54569);}
		else if(var_run0[0]==143957){weight=(1./55730)*(26516.592/35322028.2939999923);vx=0.0911123;vy= 0.02672;vz=gRandom->Gaus(-0.272577,6.63873);}
		else if(var_run0[0]==143959){weight=(1./47626)*(21580.020/35322028.2939999923);vx=0.0914659;vy= 0.0259829;vz=gRandom->Gaus(-0.329801,6.95139);}
		else if(var_run0[0]==143960){weight=(1./41523)*(17879.991/35322028.2939999923);vx=0.0909039;vy= 0.0260524;vz=gRandom->Gaus(-0.553233,7.01246);}
		else if(var_run0[0]==143961){weight=(1./122207)*(48844.926/35322028.2939999923);vx=0.0908619;vy= 0.0260308;vz=gRandom->Gaus(-0.31423,6.91023);}
		else if(var_run0[0]==143962){weight=(1./140651)*(51250.603/35322028.2939999923);vx=0.0897868;vy= 0.0267133;vz=gRandom->Gaus(-0.573746,7.1288);}
		else if(var_run0[0]==144010){weight=(1./28515)*(16817.848/35322028.2939999923);vx=0.0907405;vy= 0.0271876;vz=gRandom->Gaus(-0.467451,6.24458);}
		else if(var_run0[0]==144011){weight=(1./162465)*(86832.814/35322028.2939999923);vx=0.0904332;vy= 0.0269183;vz=gRandom->Gaus(-0.143259,6.5091);}
		else if(var_run0[0]==144086){weight=(1./121724)*(61382.679/35322028.2939999923);vx=0.0909809;vy= 0.0266411;vz=gRandom->Gaus(-0.196337,6.41842);}
		else if(var_run0[0]==144089){weight=(1./741126)*(249878.251/35322028.2939999923);vx=0.0909148;vy= 0.026542;vz=gRandom->Gaus(-0.262694,6.64063);}
		else if(var_run0[0]==144112){weight=(1./678091)*(289196.433/35322028.2939999923);vx=0.0898406;vy= 0.0243989;vz=gRandom->Gaus(-0.142117,6.56852);}
		else if(var_run0[0]==144114){weight=(1./34548)*(11552.512/35322028.2939999923);vx=0.0899083;vy= 0.0242026;vz=gRandom->Gaus(-0.00485743,6.58373);}
		else if(var_run0[0]==146428){weight=(1./5521)*(7395.316/35322028.2939999923);vx=0.0958546;vy= 0.0226069;vz=gRandom->Gaus(-0.0910449,5.66741);}
		else if(var_run0[0]==146430){weight=(1./5725)*(7129.284/35322028.2939999923);vx=0.0971417;vy= 0.0224034;vz=gRandom->Gaus(-0.0282573,5.87648);}
		else if(var_run0[0]==146431){weight=(1./1685)*(1980.971/35322028.2939999923);vx=0.09478;vy= 0.0245772;vz=gRandom->Gaus(0.0771557,6.07389);}
		else if(var_run0[0]==146436){weight=(1./35494)*(33426.223/35322028.2939999923);vx=0.0965146;vy= 0.022331;vz=gRandom->Gaus(-0.0526279,6.05393);}
		else if(var_run0[0]==146437){weight=(1./75347)*(53889.042/35322028.2939999923);vx=0.0956382;vy= 0.022949;vz=gRandom->Gaus(-0.084099,6.12018);}
		else if(var_run0[0]==146511){weight=(1./125078)*(260980.745/35322028.2939999923);vx=0.0956708;vy= 0.0210757;vz=gRandom->Gaus(0.161747,6.38767);}
		else if(var_run0[0]==146513){weight=(1./2493)*(3726.779/35322028.2939999923);vx=0.0983859;vy= 0.0215309;vz=gRandom->Gaus(0.0121315,6.41893);}
		else if(var_run0[0]==146514){weight=(1./191363)*(216421.987/35322028.2939999923);vx=0.0960731;vy= 0.0211199;vz=gRandom->Gaus(0.0439603,6.51007);}
		else if(var_run0[0]==146644){weight=(1./564044)*(1051881.804/35322028.2939999923);vx=0.0928772;vy= 0.0226648;vz=gRandom->Gaus(0.115525,5.91453);}
		else if(var_run0[0]==146804){weight=(1./140306)*(560556.097/35322028.2939999923);vx=0.0941828;vy= 0.0201881;vz=gRandom->Gaus(0.32071,6.16763);}
		else if(var_run0[0]==146807){weight=(1./76527)*(194015.961/35322028.2939999923);vx=0.0944907;vy= 0.0201597;vz=gRandom->Gaus(0.342937,6.3638);}
		else if(var_run0[0]==146944){weight=(1./125839)*(481088.796/35322028.2939999923);vx=0.0943141;vy= 0.0199461;vz=gRandom->Gaus(0.0295785,6.13587);}
		else if(var_run0[0]==147043){weight=(1./97753)*(251256.184/35322028.2939999923);vx=0.0964929;vy= 0.0213839;vz=gRandom->Gaus(0.0114308,5.63317);}
		else if(var_run0[0]==147048){weight=(1./99538)*(312050.409/35322028.2939999923);vx=0.0962723;vy= 0.0207036;vz=gRandom->Gaus(-0.0644547,5.33122);}
		else if(var_run0[0]==147114){weight=(1./120185)*(461879.570/35322028.2939999923);vx=0.0989798;vy= 0.0209683;vz=gRandom->Gaus(-0.14225,6.07775);}
		else if(var_run0[0]==147115){weight=(1./178556)*(481275.649/35322028.2939999923);vx=0.0988413;vy= 0.02093;vz=gRandom->Gaus(-0.180527,6.22696);}
		else if(var_run0[0]==147116){weight=(1./12538)*(29465.313/35322028.2939999923);vx=0.0996056;vy= 0.0212999;vz=gRandom->Gaus(-0.24068,6.33674);}
		else if(var_run0[0]==147196){weight=(1./11896)*(135377.274/35322028.2939999923);vx=0.0971283;vy=0.0205188;vz=gRandom->Gaus(-0.15672,5.88783);}
		else if(var_run0[0]==147214){weight=(1./18858)*(89278.344/35322028.2939999923);vx=0.0973382;vy=0.0206106;vz=gRandom->Gaus(-0.115323,6.35223);}
		else if(var_run0[0]==147216){weight=(1./15384)*(69158.433/35322028.2939999923);vx=0.0973536;vy=0.0206799;vz=gRandom->Gaus(-0.0602173,6.29366);}
		else if(var_run0[0]==147217){weight=(1./49125)*(206765.795/35322028.2939999923);vx=0.0975467;vy=0.0204553;vz=gRandom->Gaus(-0.141993,6.34426);}
		else if(var_run0[0]==147218){weight=(1./11550)*(45980.618/35322028.2939999923);vx=0.0978289;vy=0.0204433;vz=gRandom->Gaus(-0.232816,6.43283);}
		else if(var_run0[0]==147219){weight=(1./82627)*(296519.068/35322028.2939999923);vx=0.0977398;vy=0.0203929;vz=gRandom->Gaus(-0.130156,6.41196);}
		else if(var_run0[0]==147222){weight=(1./138507)*(371486.516/35322028.2939999923);vx=0.0978994;vy=0.0204988;vz=gRandom->Gaus(-0.111258,6.45418);}
		else if(var_run0[0]==147284){weight=(1./36684)*(383817.309/35322028.2939999923);vx=0.0957541;vy=0.0207941;vz=gRandom->Gaus(-0.0711814,6.10297);}
		else if(var_run0[0]==147390){weight=(1./131732)*(1054855.831/35322028.2939999923);vx=0.0963318;vy=0.0210244;vz=gRandom->Gaus(-0.25297,6.21035);}
		else if(var_run0[0]==147450){weight=(1./17127)*(167892.281/35322028.2939999923);vx=0.0957754;vy=0.0205546;vz=gRandom->Gaus(-0.142986,5.2659);}
		else if(var_run0[0]==147451){weight=(1./24392)*(233007.336/35322028.2939999923);vx=0.0959841;vy=0.0207997;vz=gRandom->Gaus(-0.100564,5.39198);}
		else if(var_run0[0]==147452){weight=(1./9302)*(67754.467/35322028.2939999923);vx=0.096096;vy=0.0208909;vz=gRandom->Gaus(-0.104289,5.74993);}
		else if(var_run0[0]==147453){weight=(1./31404)*(218624.969/35322028.2939999923);vx=0.0962718;vy=0.0208479;vz=gRandom->Gaus(-0.0206898,5.80622);}
		else if(var_run0[0]==147454){weight=(1./19215)*(126964.310/35322028.2939999923);vx=0.096188;vy=0.0208097;vz=gRandom->Gaus(-0.101956,5.88946);}
		else if(var_run0[0]==147754){weight=(1./75667)*(683515.893/35322028.2939999923);vx=0.0953413;vy=0.0206731;vz=gRandom->Gaus(0.0989693,5.7133);}
		else if(var_run0[0]==147755){weight=(1./33558)*(219110.476/35322028.2939999923);vx=0.095187;vy=0.0204986;vz=gRandom->Gaus(-0.122618,6.08225);}
		else if(var_run0[0]==147757){weight=(1./129222)*(446518.043/35322028.2939999923);vx=0.0953226;vy=0.0207534;vz=gRandom->Gaus(-0.154783,6.20663);}
		else if(var_run0[0]==147926){weight=(1./87566)*(937981.139/35322028.2939999923);vx=0.097043;vy=0.0208875;vz=gRandom->Gaus(0.396275,6.0927);}
		else if(var_run0[0]==147927){weight=(1./38137)*(254404.376/35322028.2939999923);vx=0.0971939;vy=0.0210484;vz=gRandom->Gaus(0.428636,6.30259);}
		else if(var_run0[0]==147929){weight=(1./166271)*(879561.939/35322028.2939999923);vx=0.0971233;vy=0.0209075;vz=gRandom->Gaus(0.208179,6.39423);}
		else if(var_run0[0]==148002){weight=(1./24214)*(323516.558/35322028.2939999923);vx=0.0960555;vy=0.0246036;vz=gRandom->Gaus(1.26235,5.3604);}
		else if(var_run0[0]==148029){weight=(1./95024)*(872764.458/35322028.2939999923);vx=0.0949575;vy=0.0195721;vz=gRandom->Gaus(0.88145,5.82658);}
		else if(var_run0[0]==148031){weight=(1./280800)*(950013.047/35322028.2939999923);vx=0.0952966;vy=0.0196396;vz=gRandom->Gaus(1.01894,6.08999);}
		else if(var_run0[0]==148032){weight=(1./121519)*(216311.823/35322028.2939999923);vx=0.095019;vy=0.0197582;vz=gRandom->Gaus(1.07986,6.22811);}
		else if(var_run0[0]==148058){weight=(1./12838)*(219753.667/35322028.2939999923);vx=0.0959153;vy=0.0214188;vz=gRandom->Gaus(0.626738,5.48652);}
		else if(var_run0[0]==148822){weight=(1./126563)*(1096016.253/35322028.2939999923);vx=0.0935971;vy=0.0216408;vz=gRandom->Gaus(0.587424,6.04144);}
		else if(var_run0[0]==148829){weight=(1./85788)*(650587.470/35322028.2939999923);vx=0.0936268;vy=0.0218682;vz=gRandom->Gaus(0.637791,6.23917);}
		else if(var_run0[0]==148860){weight=(1./9560)*(107338.259/35322028.2939999923);vx=0.0933432;vy=0.0211921;vz=gRandom->Gaus(0.686605,5.47467);}
		else if(var_run0[0]==148862){weight=(1./221601)*(2114190.211/35322028.2939999923);vx=0.0936483;vy=0.0211157;vz=gRandom->Gaus(0.692274,5.82953);}
		else if(var_run0[0]==148864){weight=(1./223111)*(1636848.589/35322028.2939999923);vx=0.0937442;vy=0.0211399;vz=gRandom->Gaus(0.942932,6.30535);}
		else if(var_run0[0]==148952){weight=(1./74356)*(793575.245/35322028.2939999923);vx=0.0932205;vy=0.0224833;vz=gRandom->Gaus(1.54875,5.38422);}
		else if(var_run0[0]==148953){weight=(1./39090)*(391755.240/35322028.2939999923);vx=0.0935521;vy=0.0224865;vz=gRandom->Gaus(1.57147,5.56769);}
		else if(var_run0[0]==149003){weight=(1./57921)*(604993.562/35322028.2939999923);vx=0.093383;vy=0.0218294;vz=gRandom->Gaus(1.47226,5.48786);}
		else if(var_run0[0]==149011){weight=(1./274237)*(2454242.751/35322028.2939999923);vx=0.0934937;vy=0.0219068;vz=gRandom->Gaus(1.562,5.85975);}
		else if(var_run0[0]==149058){weight=(1./24275)*(190179.491/35322028.2939999923);vx=0.0936254;vy=0.0219729;vz=gRandom->Gaus(1.66063,6.13628);}
		else if(var_run0[0]==149063){weight=(1./37954)*(290928.847/35322028.2939999923);vx=0.0936748;vy=0.0219186;vz=gRandom->Gaus(1.64352,6.16528);}
		else if(var_run0[0]==149181){weight=(1./650080)*(4973950.528/35322028.2939999923);vx=0.0935259;vy=0.0242632;vz=gRandom->Gaus(1.52393,6.10693);}
		else if(var_run0[0]==149182){weight=(1./157355)*(861507.489/35322028.2939999923);vx=0.0935598;vy=0.0240047;vz=gRandom->Gaus(1.22342,6.45628);}
		else if(var_run0[0]==149291){weight=(1./273461)*(1928712.013/35322028.2939999923);vx=0.0943346;vy=0.0236119;vz=gRandom->Gaus(1.07934,5.71377);}
		else if(var_run0[0]==149294){weight=(1./65824)*(341828.955/35322028.2939999923);vx=0.0947313;vy=0.0237154;vz=gRandom->Gaus(1.07326,5.90896);}

//		cout<<i<<"  -> run="<<var_run0[0]<<" ==> w="<<weight<<"; vx="<<vx<<"; vy="<<vy<<"; vz="<<vz<<endl;
	
		bool exclu(true);	
		double close_chi2(0.), close_ndof(0), close_nHits(0), close_pt(0), close_pz(0), close_eta(0), close_phi(0), close_distD(999.), close_distT(0), close_distZ(0);
		bool close_purity(false);
		for(int k=0; k<nTrack; k++){
			double distance_vtx_tk=sqrt(pow(var_TrackX0[k]-vx,2)+pow(var_TrackY0[k]-vy,2)+pow(var_TrackZ0[k]-vz,2));
			if(distance_vtx_tk<0.2 && distance_vtx_tk<close_distD){
				exclu=false;
				close_distD=distance_vtx_tk;
				close_chi2=var_TrackChi2[k]; close_ndof=var_TrackNdof[k]; close_nHits=var_TrackHits[k];
				close_pt=var_TrackPt0[k]; close_pz=var_TrackPz0[k]; close_eta=var_TrackEta[k]; close_phi=var_TrackPhi[k];
				close_distT=sqrt(pow(var_TrackX0[k]-vx,2)+pow(var_TrackY0[k]-vy,2)); close_distZ=sqrt(pow(var_TrackZ0[k]-vz,2));
				close_purity=var_TrackQuality0[k];
			}
		}
		if(exclu==true){nPassVeto_loose[nValidVtx_loose]+=weight;}
		if(exclu==false){
			hDist_D0->Fill(close_distD,weight) ;
			hDist_T0->Fill(close_distT,weight);
			hDist_Z0->Fill(close_distZ,weight);
			hpT->Fill(close_pt,weight);
			hpZ->Fill(close_pz,weight);
			heta->Fill(close_eta,weight);
			hphi->Fill(close_phi/pi,weight);
			hchi2->Fill(close_chi2,weight);
			hProb->Fill(TMath::Prob(close_chi2,close_ndof+0.5),weight);
			hPure->Fill(close_purity,weight);
			hHits->Fill(close_nHits,weight);
		}

		hnVtx_loose_0->Fill(nValidVtx_loose,weight);
		nTestVeto_loose[nValidVtx_loose]+=weight;
		}
	}
    }
  }
 }

cout<<"VERTEX Loose selection:"<<endl;
double total_loose(0);
double correction(0);
for(int j=0; j<20; j++){
	total_loose+=nTestVeto_loose[j];
}
for(int j=0; j<20; j++){
	if(nTestVeto_loose[j]!=0) cout<<j<<" extra vertex: "<<nPassVeto_loose[j]<<" / "<<nTestVeto_loose[j]<<"\t = "<<(double (nPassVeto_loose[j])/(nTestVeto_loose[j]))<<" +/- "<<((nTestVeto_loose[j]-nPassVeto_loose[j])*sqrt(nPassVeto_loose[j]))/(nTestVeto_loose[j]*nTestVeto_loose[j])<<"  \t  f= "<<(double (nTestVeto_loose[j])/(total_loose))<<endl;
	if(nTestVeto_loose[j]!=0) correction+=(double (nPassVeto_loose[j])/(nTestVeto_loose[j]))*(double (nTestVeto_loose[j])/(total_loose));
}

cout<<""<<endl;
cout<<"==> correction = "<<correction<<endl;

TCanvas *Vertex = new TCanvas("Vertex","Vertexrimeter",800,500);
   Vertex->SetFillColor(0);
   Vertex->SetBorderMode(0);
   Vertex->SetBorderSize(2);
   Vertex->SetFrameBorderMode(0);
Vertex->cd(1);
hnVtx_loose_0->GetXaxis()->SetTitle("n Primary vertex");
hnVtx_loose_0->GetYaxis()->SetTitle("# events");
hnVtx_loose_0->SetLineColor(2);
hnVtx_loose_0->Draw("hist");

  TFile *thefile1 = new TFile("global_XXX_NUM_XXX.root","RECREATE");
  thefile1->cd();
  hnVtx_loose_0->Write();
  hDist_D0->Write();
  hDist_T0->Write();
  hDist_Z0->Write();
  hpT->Write();
  hpZ->Write();
  heta->Write();
  hphi->Write();
  hchi2->Write();
  hProb->Write();
  hPure->Write();
  hHits->Write();

/*
ZDC->SaveAs("global_zdc_1.root");
Calo2->SaveAs("global_calo_1.root");
Calo3->SaveAs("global_etaphi_1.root");
Vertex->SaveAs("global_vertex_1.root");
Tracks->SaveAs("global_tracks_1.root");
*/
if(1==0){
// Save into histo
//Draw
//---------------------------
TCanvas *Calo = new TCanvas("Calo","Calorimeter",800,500);
   Calo->SetFillColor(0);
   Calo->SetBorderMode(0);
   Calo->SetBorderSize(2);
   Calo->SetFrameBorderMode(0);
Calo->Divide(3,2);
Calo->cd(4);
hnEB0->GetXaxis()->SetTitle("n EB towers (E>5GeV)");
hnEB0->GetYaxis()->SetTitle("# events");
hnEB1->SetLineColor(860);
hnEB0->SetLineColor(2);
hnEB0->Draw();
hnEB1->Draw("same");
hnEB2->Draw("same");

Calo->cd(5);
hnEE0->GetXaxis()->SetTitle("n EE towers (E>5GeV)");
hnEE0->GetYaxis()->SetTitle("# events");
hnEE1->SetLineColor(860);
hnEE0->SetLineColor(2);
hnEE0->Draw();
hnEE1->Draw("same");
hnEE2->Draw("same");

Calo->cd(1);
hnHB0->GetXaxis()->SetTitle("n HB towers (E>5GeV)");
hnHB0->GetYaxis()->SetTitle("# events");
hnHB1->SetLineColor(860);
hnHB0->SetLineColor(2);
hnHB0->Draw();
hnHB1->Draw("same");
hnHB2->Draw("same");

Calo->cd(2);
hnHE0->GetXaxis()->SetTitle("n HE towers (E>5GeV)");
hnHE0->GetYaxis()->SetTitle("# events");
hnHE1->SetLineColor(860);
hnHE0->SetLineColor(2);
hnHE0->Draw();
hnHE1->Draw("same");
hnHE2->Draw("same");

Calo->cd(3);
hnHF0->GetXaxis()->SetTitle("n HF towers (E>5GeV)");
hnHF0->GetYaxis()->SetTitle("# events (grouped by 5)");
hnHF1->SetLineColor(860);
hnHF0->SetLineColor(2);
hnHF0->Draw();
hnHF1->Draw("same");
hnHF2->Draw("same");

Calo->cd(6);
hnCastor0->GetXaxis()->SetTitle("n Castor towers (E>3GeV)");
hnCastor0->GetYaxis()->SetTitle("# events");
hnCastor1->SetLineColor(860);
hnCastor0->SetLineColor(2);
hnCastor0->Draw();
hnCastor1->Draw("same");
//Calo->SaveAs("global_calo.C");
/*
TCanvas *CASTOR = new TCanvas("CASTOR","CASTOR",800,550);
   CASTOR->SetFillColor(0);
   CASTOR->SetBorderMode(0);
   CASTOR->SetBorderSize(2);
   CASTOR->SetFrameBorderMode(0);

CASTOR->Divide(2,1);
CASTOR->cd(1);
CastorSumE0->GetXaxis()->SetTitle("Castor Sum Energy [GeV]");
CastorSumE0->GetYaxis()->SetTitle("# events / 5 GeV");
CastorSumE1->SetLineColor(860);
CastorSumE0->Draw();
CastorSumE2->Draw("same");

CASTOR->cd(2);
CastorEphi0->GetXaxis()->SetTitle("Castor #phi");
CastorEphi0->GetYaxis()->SetTitle("Castor Energy [GeV]");
CastorEphi0->Draw("colz");
*/
TCanvas *ZDC = new TCanvas("ZDC","ZDC",800,550);
   ZDC->SetFillColor(0);
   ZDC->SetBorderMode(0);
   ZDC->SetBorderSize(2);
   ZDC->SetFrameBorderMode(0);

ZDC->Divide(4,2);
ZDC->cd(1);
ZDCemplus0->GetXaxis()->SetTitle("ZDC + em [GeV]");
ZDCemplus0->GetYaxis()->SetTitle("# events / 20 GeV");
ZDCemplus1->SetLineColor(860);
ZDCemplus0->Draw();
ZDCemplus2->Draw("same");

ZDC->cd(2);
ZDCemminus0->GetXaxis()->SetTitle("ZDC - em [GeV]");
ZDCemminus0->GetYaxis()->SetTitle("# events / 20 GeV");
ZDCemminus1->SetLineColor(860);
ZDCemminus0->Draw();
ZDCemminus2->Draw("same");

ZDC->cd(3);
ZDChadplus0->GetXaxis()->SetTitle("ZDC + had [GeV]");
ZDChadplus0->GetYaxis()->SetTitle("# events / 20 GeV");
ZDChadplus1->SetLineColor(860);
ZDChadplus0->Draw();
ZDChadplus2->Draw("same");

ZDC->cd(4);
ZDChadminus0->GetXaxis()->SetTitle("ZDC - had [GeV]");
ZDChadminus0->GetYaxis()->SetTitle("# events / 20 GeV");
ZDChadminus1->SetLineColor(860);
ZDChadminus0->Draw();
ZDChadminus2->Draw("same");

ZDC->cd(5);
ZDCtime0->GetXaxis()->SetTitle("ZDC hit time [???]");
ZDCtime0->GetYaxis()->SetTitle("# events / 1 ???");
ZDCtime1->SetLineColor(860);
ZDCtime0->Draw();
ZDCtime2->Draw("same");

ZDC->cd(7);
ZDCenergyEM0->GetXaxis()->SetTitle("ZDC hit EM energy [GeV ??]");
ZDCenergyEM0->GetYaxis()->SetTitle("# events / 30 GeV ???");
ZDCenergyEM1->SetLineColor(860);
ZDCenergyEM0->Draw();
ZDCenergyEM2->Draw("same");

ZDC->cd(8);
ZDCenergyHAD0->GetXaxis()->SetTitle("ZDC hit HAD energy [GeV ??]");
ZDCenergyHAD0->GetYaxis()->SetTitle("# events / 300 GeV ???");
ZDCenergyHAD1->SetLineColor(860);
ZDCenergyHAD0->Draw();
ZDCenergyHAD2->Draw("same");
//ZDC->SaveAs("global_zdc.C");
}
cout << "END" << endl;   
}

