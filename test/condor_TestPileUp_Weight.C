void EFF_BX_XXX_NUM_XXX_XXX_LABEL_XXX() // 
{
gROOT->Reset();
gStyle->SetPalette(1);

#define pi 3.1413565359
//definition des fichiers + Tree
//  TFile *f0 = new TFile("zerobias_merge1.root"); // 
//  TTree *t0 = f0->Get("ntp1");

  TChain *t0 = new TChain("ntp1");
//  t0->Add("/storage/data/cms/store/user/schul/ZeroBias_prescale/prescale_ZBRunA_*.root");
//  t0->Add("/storage/data/cms/store/user/schul/ZeroBias_prescale/prescale_MBRunA_*.root");
  t0->Add("/storage/data/cms/store/user/schul/ZeroBias_prescale/prescale_MBRunB_*.root");

// definitions des Trees pour la sauvegarde
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

      double nTestVeto_bx[3263];
      double nPassVeto_bx[3263];
      for(int j=0; j<=3262; j++){nPassVeto_bx[j]=0; nTestVeto_bx[j]=0;}

        ifstream ifs1("/storage/data/cms/users/schul/Test_condor/CMSSW_3_8_6/src/PUvsRunEfficiency/distributionvertices_XXX_NUM_XXX.txt");

        Float_t tmpweight;
        Int_t tmpls, tmprun, tmpbx; 
        Int_t runweightind[824630];
        Int_t lsweightind[824630];
        Int_t bxweightind[824630];
        Float_t weightind[824630];

        Int_t k = 0;
        while(ifs1 >> tmprun >>  tmpls >> tmpbx >> tmpweight)
        {
                runweightind[k] = tmprun;
                lsweightind[k] = tmpls;
                bxweightind[k] = tmpbx;
                weightind[k] = tmpweight;
                k++;
        }
        ifs1.close();




  int filter0Gen(0);
  int filter0Track(0);
  int filter0Events(0);
  int n_cond0(0);
  int n_cond1(0);
  int n_cond2(0);
  for(Int_t i = XXX_START_XXX;i < XXX_STOP_XXX;i++){
    t0->GetEntry(i);
    if(var_run0[0]==XXX_NUM_XXX){
      int nPrimVtx = var_nvtx0[0];
      int nTrack=var_nTrack0[0];
      int nCalo=var_ncalo0[0];
      int nValidVtx_bx(0);
      if(!(techBit0[0][36]||techBit0[0][37]||techBit0[0][38]||techBit0[0][39]) && !((techBit0[0][42]&&!(techBit0[0][43]))||(techBit0[0][43]&&!(techBit0[0][42]))) ){
	if(techBit0[0][4]){ 
//		cout<<"beam-crossing"<<endl;
//		cout<<"----------"<<var_run0[0]<<endl;
/*		for(Int_t j=0; j<nPrimVtx; j++){if(var_vtxTrack0[j]>0 && TMath::Prob(var_vertexChi2_0[j],var_vertexNdf0[j]+0.5)>0.001) nValidVtx_bx++;}*/

		double vx(0),vy(0),vz(0),weight(0);
		if(var_run0[0]==140058){vx=0.100794;vy= 0.0143116;vz=gRandom->Gaus(0.271343,6.06352);}
		else if(var_run0[0]==140059){vx=0.100335;vy= 0.014676;vz=gRandom->Gaus(0.0995124,6.42079);}
		else if(var_run0[0]==140124){vx=0.095623;vy= 0.015889;vz=gRandom->Gaus(-0.152702,6.36746);}
		else if(var_run0[0]==140158){vx=0.0976755;vy= 0.0155846;vz=gRandom->Gaus(-0.587424,6.15816);}
		else if(var_run0[0]==140159){vx=0.0969652;vy= 0.0166257;vz=gRandom->Gaus(0.0376294,6.06082);}
		else if(var_run0[0]==140160){vx=0.0980627;vy= 0.0177718;vz=gRandom->Gaus(-0.49599,6.38909);}
		else if(var_run0[0]==140331){vx=0.0966409;vy= 0.0135332;vz=gRandom->Gaus(-0.303373,6.31805);}
		else if(var_run0[0]==140359){vx=0.0940866;vy= 0.0146067;vz=gRandom->Gaus(-0.606706,5.95362);}
		else if(var_run0[0]==140361){vx=0.0948746;vy= 0.0126356;vz=gRandom->Gaus(0.0671472,6.31849);}
		else if(var_run0[0]==140362){vx=0.094178;vy= 0.0116153;vz=gRandom->Gaus(0.080781,6.46025);}
		else if(var_run0[0]==140379){vx=0.0950191;vy= 0.0148368;vz=gRandom->Gaus(0.0353279,6.02474);}
		else if(var_run0[0]==140381){vx=0.0945141;vy= 0.0137429;vz=gRandom->Gaus(-0.552448,7.15243);}
		else if(var_run0[0]==140382){vx=0.101838;vy= 0.0145831;vz=gRandom->Gaus(-0.138335,6.63679);}
		else if(var_run0[0]==140383){vx=0.0938265;vy= 0.0144248;vz=gRandom->Gaus(0.305865,6.27479);}
		else if(var_run0[0]==140385){vx=0.092982;vy= 0.0103742;vz=gRandom->Gaus(-0.434873,5.98944);}
		else if(var_run0[0]==140386){vx=0.091979;vy= 0.0139829;vz=gRandom->Gaus(-0.850734,6.33796);}
		else if(var_run0[0]==140387){vx=0.0947584;vy= 0.0133499;vz=gRandom->Gaus(-0.421871,6.3771);}
		else if(var_run0[0]==140388){vx=0.0926958;vy= 0.0137019;vz=gRandom->Gaus(-0.684303,6.3918);}
		else if(var_run0[0]==140399){vx=0.0932278;vy= 0.0132279;vz=gRandom->Gaus(-0.106708,6.51808);}
		else if(var_run0[0]==140401){vx=0.0931317;vy= 0.0145009;vz=gRandom->Gaus(0.186389,6.12158);}
		else if(var_run0[0]==141880){vx=0.0945815;vy= 0.0103382;vz=gRandom->Gaus(0.0246312,6.22766);}
		else if(var_run0[0]==141881){vx=0.095274;vy= 0.0111391;vz=gRandom->Gaus(-0.199211,6.47942);}
		else if(var_run0[0]==141956){vx=0.0953798;vy= 0.0125119;vz=gRandom->Gaus(-0.0746152,6.42123);}
		else if(var_run0[0]==141957){vx=0.0945655;vy= 0.0142757;vz=gRandom->Gaus(0.00655455,6.55368);}
		else if(var_run0[0]==141958){vx=0.0929087;vy= 0.0126065;vz=gRandom->Gaus(-0.149173,6.28642);}
		else if(var_run0[0]==141959){vx=0.0938471;vy= 0.0115851;vz=gRandom->Gaus(-0.175745,6.45066);}
		else if(var_run0[0]==141960){vx=0.0944724;vy= 0.0129941;vz=gRandom->Gaus(-0.164615,6.46335);}
		else if(var_run0[0]==141961){vx=0.0980082;vy= 0.0150895;vz=gRandom->Gaus(0.608491,5.91946);}
		else if(var_run0[0]==142035){vx=0.0937506;vy= 0.0126625;vz=gRandom->Gaus(-0.0332531,6.38996);}
		else if(var_run0[0]==142036){vx=0.0932561;vy= 0.0113626;vz=gRandom->Gaus(-0.163195,6.35545);}
		else if(var_run0[0]==142038){vx=0.094352;vy= 0.0121395;vz=gRandom->Gaus(0.0304028,6.45642);}
		else if(var_run0[0]==142039){vx=0.0963541;vy= 0.00278066;vz=gRandom->Gaus(-1.49764,6.2135);}
		else if(var_run0[0]==142040){vx=0.0926035;vy= 0.0115514;vz=gRandom->Gaus(-0.633006,6.5702);}
		else if(var_run0[0]==142076){vx=0.0935682;vy= 0.0126095;vz=gRandom->Gaus(-0.188606,6.16845);}
		else if(var_run0[0]==142128){vx=0.0889196;vy= 0.0105045;vz=gRandom->Gaus(-0.390518,6.1126);}
		else if(var_run0[0]==142129){vx=0.0925292;vy= 0.0119757;vz=gRandom->Gaus(-0.278862,6.28139);}
		else if(var_run0[0]==142130){vx=0.0921308;vy= 0.0125107;vz=gRandom->Gaus(-0.128887,6.45175);}
		else if(var_run0[0]==142132){vx=0.0916467;vy= 0.0118682;vz=gRandom->Gaus(-0.230859,6.23377);}
		else if(var_run0[0]==142135){vx=0.092053;vy= 0.0114144;vz=gRandom->Gaus(-0.337647,6.62729);}
		else if(var_run0[0]==142136){vx=0.0926426;vy= 0.0109235;vz=gRandom->Gaus(-0.310499,6.2499);}
		else if(var_run0[0]==142137){vx=0.0923496;vy= 0.0120018;vz=gRandom->Gaus(-0.384746,6.49448);}
		else if(var_run0[0]==142189){vx=0.0954219;vy= 0.0118383;vz=gRandom->Gaus(-0.782422,6.66343);}
		else if(var_run0[0]==142191){vx=0.0948962;vy= 0.0103308;vz=gRandom->Gaus(-0.655677,6.57426);}
		else if(var_run0[0]==142264){vx=0.0908284;vy= 0.0098164;vz=gRandom->Gaus(-0.00313011,6.35974);}
		else if(var_run0[0]==142265){vx=0.0912857;vy= 0.0120297;vz=gRandom->Gaus(-0.363364,6.21843);}
		else if(var_run0[0]==142303){vx=0.0936276;vy= 0.0123935;vz=gRandom->Gaus(-0.19073,6.3332);}
		else if(var_run0[0]==142304){vx=0.0953863;vy= 0.0124424;vz=gRandom->Gaus(0.805653,6.46557);}
		else if(var_run0[0]==142305){vx=0.0930993;vy= 0.0123919;vz=gRandom->Gaus(-0.259715,6.18923);}
		else if(var_run0[0]==142308){vx=0.0936005;vy= 0.0126497;vz=gRandom->Gaus(0.129302,6.52169);}
		else if(var_run0[0]==142309){vx=0.094416;vy= 0.0119167;vz=gRandom->Gaus(-0.137821,6.18922);}
		else if(var_run0[0]==142311){vx=0.0935093;vy= 0.0112015;vz=gRandom->Gaus(0.0334105,6.43299);}
		else if(var_run0[0]==142312){vx=0.0929723;vy= 0.0114759;vz=gRandom->Gaus(-0.409077,6.484);}
		else if(var_run0[0]==142313){vx=0.0948427;vy= 0.00903676;vz=gRandom->Gaus(-0.120377,6.23579);}
		else if(var_run0[0]==142413){vx=0.0952755;vy= 0.0124298;vz=gRandom->Gaus(-0.246799,6.15094);}
		else if(var_run0[0]==142414){vx=0.0927713;vy= 0.0139964;vz=gRandom->Gaus(-0.537765,6.66217);}
		else if(var_run0[0]==142419){vx=0.0948953;vy= 0.0116192;vz=gRandom->Gaus(0.216078,6.65546);}
		else if(var_run0[0]==142422){vx=0.094908;vy= 0.0130058;vz=gRandom->Gaus(-0.254085,6.48445);}
		else if(var_run0[0]==142513){vx=0.0936079;vy= 0.0128486;vz=gRandom->Gaus(0.0553905,6.29823);}
		else if(var_run0[0]==142514){vx=0.0928844;vy= 0.0117103;vz=gRandom->Gaus(0.471611,6.16774);}
		else if(var_run0[0]==142523){vx=0.0949826;vy= 0.0135065;vz=gRandom->Gaus(-0.157981,6.08754);}
		else if(var_run0[0]==142524){vx=0.0955073;vy= 0.013181;vz=gRandom->Gaus(0.0191425,6.16217);}
		else if(var_run0[0]==142525){vx=0.0960522;vy= 0.0123776;vz=gRandom->Gaus(-0.327897,6.2091);}
		else if(var_run0[0]==142528){vx=0.0964386;vy= 0.0132233;vz=gRandom->Gaus(-0.0868653,6.37082);}
		else if(var_run0[0]==142530){vx=0.093721;vy= 0.0132565;vz=gRandom->Gaus(0.0801791,6.68744);}
		else if(var_run0[0]==142535){vx=0.0960788;vy= 0.0107834;vz=gRandom->Gaus(-0.630023,6.30529);}
		else if(var_run0[0]==142537){vx=0.0960362;vy= 0.0110998;vz=gRandom->Gaus(-0.198908,5.99704);}
		else if(var_run0[0]==142557){vx=0.0972623;vy= 0.0142449;vz=gRandom->Gaus(0.0213277,6.11328);}
		else if(var_run0[0]==142558){vx=0.0974421;vy= 0.0130626;vz=gRandom->Gaus(-0.226053,6.29172);}
		else if(var_run0[0]==142657){vx=0.0885586;vy= 0.0294256;vz=gRandom->Gaus(-0.604557,3.23314);}
		else if(var_run0[0]==142658){vx=0.0924293;vy= 0.00755789;vz=gRandom->Gaus(-0.184352,5.82713);}
		else if(var_run0[0]==142659){vx=0.0886328;vy= 0.0102111;vz=gRandom->Gaus(-0.365842,6.83124);}
		else if(var_run0[0]==142660){vx=0.0883789;vy= 0.0122418;vz=gRandom->Gaus(-0.191347,6.07437);}
		else if(var_run0[0]==142661){vx=0.0931001;vy= 0.0118804;vz=gRandom->Gaus(-0.437016,6.88346);}
		else if(var_run0[0]==142662){vx=0.0911049;vy= 0.0155607;vz=gRandom->Gaus(-0.102537,6.77731);}
		else if(var_run0[0]==142663){vx=0.089634;vy= 0.0251548;vz=gRandom->Gaus(0.0778907,6.7571);}
		else if(var_run0[0]==142664){vx=0.0911891;vy= 0.0246174;vz=gRandom->Gaus(0.178342,7.0693);}
		else if(var_run0[0]==142928){vx=0.0898088;vy= 0.0277331;vz=gRandom->Gaus(-0.21254,6.46752);}
		else if(var_run0[0]==142933){vx=0.0897305;vy= 0.0278015;vz=gRandom->Gaus(-0.0292006,6.39151);}
		else if(var_run0[0]==142935){vx=0.0889369;vy= 0.0274965;vz=gRandom->Gaus(-0.191199,6.78442);}
		else if(var_run0[0]==142936){vx=0.0915178;vy= 0.0266161;vz=gRandom->Gaus(-0.95817,6.59995);}
		else if(var_run0[0]==142953){vx=0.0909008;vy= 0.0271429;vz=gRandom->Gaus(-0.203396,6.36024);}
		else if(var_run0[0]==142954){vx=0.0909349;vy= 0.0267102;vz=gRandom->Gaus(0.0981885,6.3887);}
		else if(var_run0[0]==142970){vx=0.089288;vy= 0.0232477;vz=gRandom->Gaus(-0.304972,6.68222);}
		else if(var_run0[0]==142971){vx=0.0907004;vy= 0.0255772;vz=gRandom->Gaus(-0.120557,6.6521);}
		else if(var_run0[0]==143004){vx=0.0910295;vy= 0.0272586;vz=gRandom->Gaus(-0.075168,6.32296);}
		else if(var_run0[0]==143005){vx=0.0905392;vy= 0.0270044;vz=gRandom->Gaus(0.0446868,6.48798);}
		else if(var_run0[0]==143006){vx=0.0908894;vy= 0.0276938;vz=gRandom->Gaus(0.157866,6.73279);}
		else if(var_run0[0]==143007){vx=0.0901827;vy= 0.0271159;vz=gRandom->Gaus(-0.237678,6.56121);}
		else if(var_run0[0]==143008){vx=0.0891607;vy= 0.0271899;vz=gRandom->Gaus(-0.107623,6.24831);}
		else if(var_run0[0]==143179){vx=0.091262;vy= 0.028573;vz=gRandom->Gaus(-0.185253,6.35637);}
		else if(var_run0[0]==143181){vx=0.0910366;vy= 0.026101;vz=gRandom->Gaus(-0.214549,6.36783);}
		else if(var_run0[0]==143187){vx=0.0911821;vy= 0.0257043;vz=gRandom->Gaus(-0.26801,6.36053);}
		else if(var_run0[0]==143191){vx=0.0917604;vy= 0.0260234;vz=gRandom->Gaus(-0.271733,6.55902);}
		else if(var_run0[0]==143192){vx=0.0889398;vy= 0.0269094;vz=gRandom->Gaus(-0.0534999,6.63494);}
		else if(var_run0[0]==143193){vx=0.0897632;vy= 0.0270503;vz=gRandom->Gaus(-0.0998894,6.41284);}
		else if(var_run0[0]==143318){vx=0.0914484;vy= 0.0263239;vz=gRandom->Gaus(0.255529,6.5057);}
		else if(var_run0[0]==143319){vx=0.0926588;vy= 0.0267072;vz=gRandom->Gaus(-0.0883608,6.41922);}
		else if(var_run0[0]==143320){vx=0.0899219;vy= 0.0253984;vz=gRandom->Gaus(-0.246676,6.28565);}
		else if(var_run0[0]==143321){vx=0.0906771;vy= 0.026495;vz=gRandom->Gaus(-0.361193,6.562);}
		else if(var_run0[0]==143322){vx=0.0898169;vy= 0.0264218;vz=gRandom->Gaus(-0.321438,6.43413);}
		else if(var_run0[0]==143323){vx=0.0900378;vy= 0.0263641;vz=gRandom->Gaus(-0.128706,6.51767);}
		else if(var_run0[0]==143326){vx=0.0896332;vy= 0.0253725;vz=gRandom->Gaus(-0.674726,6.63082);}
		else if(var_run0[0]==143327){vx=0.0906081;vy= 0.0271599;vz=gRandom->Gaus(-0.279194,6.54791);}
		else if(var_run0[0]==143328){vx=0.0907095;vy= 0.0257698;vz=gRandom->Gaus(-0.404016,6.53572);}
		else if(var_run0[0]==143657){vx=0.0905717;vy= 0.0270814;vz=gRandom->Gaus(-0.177333,6.57029);}
		else if(var_run0[0]==143665){vx=0.0896075;vy= 0.0265908;vz=gRandom->Gaus(-0.223942,6.72243);}
		else if(var_run0[0]==143727){vx=0.0912077;vy= 0.0268628;vz=gRandom->Gaus(-0.310465,6.3967);}
		else if(var_run0[0]==143731){vx=0.0905155;vy= 0.0264595;vz=gRandom->Gaus(0.257504,6.68583);}
		else if(var_run0[0]==143827){vx=0.089576;vy= 0.0277203;vz=gRandom->Gaus(-0.241826,6.52543);}
		else if(var_run0[0]==143833){vx=0.0900497;vy= 0.0275122;vz=gRandom->Gaus(-0.124198,6.77496);}
		else if(var_run0[0]==143953){vx=0.0906401;vy= 0.0263002;vz=gRandom->Gaus(-0.308051,6.53142);}
		else if(var_run0[0]==143955){vx=0.0884157;vy= 0.0259701;vz=gRandom->Gaus(-0.555993,6.69627);}
		else if(var_run0[0]==143956){vx=0.0915183;vy= 0.0250372;vz=gRandom->Gaus(-0.221492,6.54569);}
		else if(var_run0[0]==143957){vx=0.0911123;vy= 0.02672;vz=gRandom->Gaus(-0.272577,6.63873);}
		else if(var_run0[0]==143959){vx=0.0914659;vy= 0.0259829;vz=gRandom->Gaus(-0.329801,6.95139);}
		else if(var_run0[0]==143960){vx=0.0909039;vy= 0.0260524;vz=gRandom->Gaus(-0.553233,7.01246);}
		else if(var_run0[0]==143961){vx=0.0908619;vy= 0.0260308;vz=gRandom->Gaus(-0.31423,6.91023);}
		else if(var_run0[0]==143962){vx=0.0897868;vy= 0.0267133;vz=gRandom->Gaus(-0.573746,7.1288);}
		else if(var_run0[0]==144010){vx=0.0907405;vy= 0.0271876;vz=gRandom->Gaus(-0.467451,6.24458);}
		else if(var_run0[0]==144011){vx=0.0904332;vy= 0.0269183;vz=gRandom->Gaus(-0.143259,6.5091);}
		else if(var_run0[0]==144086){vx=0.0909809;vy= 0.0266411;vz=gRandom->Gaus(-0.196337,6.41842);}
		else if(var_run0[0]==144089){vx=0.0909148;vy= 0.026542;vz=gRandom->Gaus(-0.262694,6.64063);}
		else if(var_run0[0]==144112){vx=0.0898406;vy= 0.0243989;vz=gRandom->Gaus(-0.142117,6.56852);}
		else if(var_run0[0]==144114){vx=0.0899083;vy= 0.0242026;vz=gRandom->Gaus(-0.00485743,6.58373);}
		else if(var_run0[0]==146428){vx=0.0958546;vy= 0.0226069;vz=gRandom->Gaus(-0.0910449,5.66741);}
		else if(var_run0[0]==146430){vx=0.0971417;vy= 0.0224034;vz=gRandom->Gaus(-0.0282573,5.87648);}
		else if(var_run0[0]==146431){vx=0.09478;vy= 0.0245772;vz=gRandom->Gaus(0.0771557,6.07389);}
		else if(var_run0[0]==146436){vx=0.0965146;vy= 0.022331;vz=gRandom->Gaus(-0.0526279,6.05393);}
		else if(var_run0[0]==146437){vx=0.0956382;vy= 0.022949;vz=gRandom->Gaus(-0.084099,6.12018);}
		else if(var_run0[0]==146511){vx=0.0956708;vy= 0.0210757;vz=gRandom->Gaus(0.161747,6.38767);}
		else if(var_run0[0]==146513){vx=0.0983859;vy= 0.0215309;vz=gRandom->Gaus(0.0121315,6.41893);}
		else if(var_run0[0]==146514){vx=0.0960731;vy= 0.0211199;vz=gRandom->Gaus(0.0439603,6.51007);}
		else if(var_run0[0]==146644){vx=0.0928772;vy= 0.0226648;vz=gRandom->Gaus(0.115525,5.91453);}
		else if(var_run0[0]==146804){vx=0.0941828;vy= 0.0201881;vz=gRandom->Gaus(0.32071,6.16763);}
		else if(var_run0[0]==146807){vx=0.0944907;vy= 0.0201597;vz=gRandom->Gaus(0.342937,6.3638);}
		else if(var_run0[0]==146944){vx=0.0943141;vy= 0.0199461;vz=gRandom->Gaus(0.0295785,6.13587);}
		else if(var_run0[0]==147043){vx=0.0964929;vy= 0.0213839;vz=gRandom->Gaus(0.0114308,5.63317);}
		else if(var_run0[0]==147048){vx=0.0962723;vy= 0.0207036;vz=gRandom->Gaus(-0.0644547,5.33122);}
		else if(var_run0[0]==147114){vx=0.0989798;vy= 0.0209683;vz=gRandom->Gaus(-0.14225,6.07775);}
		else if(var_run0[0]==147115){vx=0.0988413;vy= 0.02093;vz=gRandom->Gaus(-0.180527,6.22696);}
		else if(var_run0[0]==147116){vx=0.0996056;vy= 0.0212999;vz=gRandom->Gaus(-0.24068,6.33674);}
		else if(var_run0[0]==147196){vx=0.0971283;vy=0.0205188;vz=gRandom->Gaus(-0.15672,5.88783);}
		else if(var_run0[0]==147214){vx=0.0973382;vy=0.0206106;vz=gRandom->Gaus(-0.115323,6.35223);}
		else if(var_run0[0]==147216){vx=0.0973536;vy=0.0206799;vz=gRandom->Gaus(-0.0602173,6.29366);}
		else if(var_run0[0]==147217){vx=0.0975467;vy=0.0204553;vz=gRandom->Gaus(-0.141993,6.34426);}
		else if(var_run0[0]==147218){vx=0.0978289;vy=0.0204433;vz=gRandom->Gaus(-0.232816,6.43283);}
		else if(var_run0[0]==147219){vx=0.0977398;vy=0.0203929;vz=gRandom->Gaus(-0.130156,6.41196);}
		else if(var_run0[0]==147222){vx=0.0978994;vy=0.0204988;vz=gRandom->Gaus(-0.111258,6.45418);}
		else if(var_run0[0]==147284){vx=0.0957541;vy=0.0207941;vz=gRandom->Gaus(-0.0711814,6.10297);}
		else if(var_run0[0]==147390){vx=0.0963318;vy=0.0210244;vz=gRandom->Gaus(-0.25297,6.21035);}
		else if(var_run0[0]==147450){vx=0.0957754;vy=0.0205546;vz=gRandom->Gaus(-0.142986,5.2659);}
		else if(var_run0[0]==147451){vx=0.0959841;vy=0.0207997;vz=gRandom->Gaus(-0.100564,5.39198);}
		else if(var_run0[0]==147452){vx=0.096096;vy=0.0208909;vz=gRandom->Gaus(-0.104289,5.74993);}
		else if(var_run0[0]==147453){vx=0.0962718;vy=0.0208479;vz=gRandom->Gaus(-0.0206898,5.80622);}
		else if(var_run0[0]==147454){vx=0.096188;vy=0.0208097;vz=gRandom->Gaus(-0.101956,5.88946);}
		else if(var_run0[0]==147754){vx=0.0953413;vy=0.0206731;vz=gRandom->Gaus(0.0989693,5.7133);}
		else if(var_run0[0]==147755){vx=0.095187;vy=0.0204986;vz=gRandom->Gaus(-0.122618,6.08225);}
		else if(var_run0[0]==147757){vx=0.0953226;vy=0.0207534;vz=gRandom->Gaus(-0.154783,6.20663);}
		else if(var_run0[0]==147926){vx=0.097043;vy=0.0208875;vz=gRandom->Gaus(0.396275,6.0927);}
		else if(var_run0[0]==147927){vx=0.0971939;vy=0.0210484;vz=gRandom->Gaus(0.428636,6.30259);}
		else if(var_run0[0]==147929){vx=0.0971233;vy=0.0209075;vz=gRandom->Gaus(0.208179,6.39423);}
		else if(var_run0[0]==148002){vx=0.0960555;vy=0.0246036;vz=gRandom->Gaus(1.26235,5.3604);}
		else if(var_run0[0]==148029){vx=0.0949575;vy=0.0195721;vz=gRandom->Gaus(0.88145,5.82658);}
		else if(var_run0[0]==148031){vx=0.0952966;vy=0.0196396;vz=gRandom->Gaus(1.01894,6.08999);}
		else if(var_run0[0]==148032){vx=0.095019;vy=0.0197582;vz=gRandom->Gaus(1.07986,6.22811);}
		else if(var_run0[0]==148058){vx=0.0959153;vy=0.0214188;vz=gRandom->Gaus(0.626738,5.48652);}
		else if(var_run0[0]==148822){vx=0.0935971;vy=0.0216408;vz=gRandom->Gaus(0.587424,6.04144);}
		else if(var_run0[0]==148829){vx=0.0936268;vy=0.0218682;vz=gRandom->Gaus(0.637791,6.23917);}
		else if(var_run0[0]==148860){vx=0.0933432;vy=0.0211921;vz=gRandom->Gaus(0.686605,5.47467);}
		else if(var_run0[0]==148862){vx=0.0936483;vy=0.0211157;vz=gRandom->Gaus(0.692274,5.82953);}
		else if(var_run0[0]==148864){vx=0.0937442;vy=0.0211399;vz=gRandom->Gaus(0.942932,6.30535);}
		else if(var_run0[0]==148952){vx=0.0932205;vy=0.0224833;vz=gRandom->Gaus(1.54875,5.38422);}
		else if(var_run0[0]==148953){vx=0.0935521;vy=0.0224865;vz=gRandom->Gaus(1.57147,5.56769);}
		else if(var_run0[0]==149003){vx=0.093383;vy=0.0218294;vz=gRandom->Gaus(1.47226,5.48786);}
		else if(var_run0[0]==149011){vx=0.0934937;vy=0.0219068;vz=gRandom->Gaus(1.562,5.85975);}
		else if(var_run0[0]==149058){vx=0.0936254;vy=0.0219729;vz=gRandom->Gaus(1.66063,6.13628);}
		else if(var_run0[0]==149063){vx=0.0936748;vy=0.0219186;vz=gRandom->Gaus(1.64352,6.16528);}
		else if(var_run0[0]==149181){vx=0.0935259;vy=0.0242632;vz=gRandom->Gaus(1.52393,6.10693);}
		else if(var_run0[0]==149182){vx=0.0935598;vy=0.0240047;vz=gRandom->Gaus(1.22342,6.45628);}
		else if(var_run0[0]==149291){vx=0.0943346;vy=0.0236119;vz=gRandom->Gaus(1.07934,5.71377);}
		else if(var_run0[0]==149294){vx=0.0947313;vy=0.0237154;vz=gRandom->Gaus(1.07326,5.90896);}

		for(int m=0; m<k; m++){
			if(var_run0[0]==runweightind[m] && var_ls0[0]==lsweightind[m] && var_bx0[0]==bxweightind[m] ){weight=weightind[m]; continue;}
		}
		if(weight==0) continue;
//		cout<<i<<"  -> run="<<var_run0[0]<<" ==> w="<<weight<<"; vx="<<vx<<"; vy="<<vy<<"; vz="<<vz<<endl;
		bool exclu(true);	
		double close_distD(999.);
		bool close_purity(false);
		for(int l=0; l<nTrack; l++){
			double distance_vtx_tk=sqrt(pow(var_TrackX0[l]-vx,2)+pow(var_TrackY0[l]-vy,2)+pow(var_TrackZ0[l]-vz,2));
			if(distance_vtx_tk<0.2 && distance_vtx_tk<close_distD){
				exclu=false;
				close_distD=distance_vtx_tk;
			}
		}
		if(exclu==true){nPassVeto_bx[var_bx0[0]]+=weight;}
		nTestVeto_bx[var_bx0[0]]+=weight;
	}
    }
  }
 }

cout<<"VERTEX Loose selection:"<<endl;
double total_loose(0);
double correction(0);
for(int i=1; i<=3262; i++){
	if(nTestVeto_bx[i]!=0) cout<<"==> correctionBX "<<i<<" = "<<nPassVeto_bx[i]<<" / "<<nTestVeto_bx[i]<<endl;
}

cout << "END" << endl;   
}

