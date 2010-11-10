void DrawOneHistogram(THStack *mcstackhist, 
		      TH1F *datahist, 
		      TH1F *mchist1,
		      TH1F *mchist2,
		      TH1F *mchist3,
		      TH1F *mchist4,
		      TH1F *mchist5,			
		      TString xaxislabel, 
		      TString yaxislabel,
		      TString plottitle)
{
	TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas MuMu",800,500); 
   	Canvas1->SetFillColor(0); 
   	Canvas1->SetBorderMode(0); 
   	Canvas1->SetBorderSize(2); 
   	Canvas1->SetFrameBorderMode(0); 
	mcstackhist->SetTitle(0);

        ci = TColor::GetColor("#ffff99");
 
	datahist->Sumw2(); 
	datahist->SetStats(0);
	datahist->SetLineWidth(2); 
	datahist->SetMarkerStyle(20); 
	mchist1->SetFillColor(ci);mchist1->SetLineWidth(1);/*mchist1->SetLineColor(4);*/ 
	mchist2->SetFillColor(30);  mchist2->SetFillStyle(3001);
	mchist3->SetFillColor(30); 
	mchist4->SetFillColor(38); 
	mchist5->SetFillColor(903); 
	mcstackhist->Add(mchist3); 
	mcstackhist->Add(mchist2); 
	mcstackhist->Add(mchist1); 
	mcstackhist->Add(mchist4); 
	mcstackhist->Add(mchist5); 
	if(mcstackhist->GetMaximum() > 0) mcstackhist.SetMaximum(mcstackhist->GetMaximum() * 1.5);    
	else mcstackhist.SetMaximum(datahist->GetMaximum() * 2.0); 
	mcstackhist->Draw(); 
	mcstackhist->GetXaxis()->SetTitle(xaxislabel); 
	mcstackhist->GetYaxis()->SetTitle(yaxislabel); 
	datahist->Draw("same"); 
	Canvas1->SaveAs(plottitle);
	
}


void DrawOneHistogramBis(THStack *mcstackhist,
                      TH1F *datahist,
                      TH1F *mchist1,
                      TH1F *mchist2,
                      TH1F *mchist3,
                      TH1F *mchist4,
                      TH1F *mchist5,
                      TH1F *mchist6,
                      TString xaxislabel,
                      TString yaxislabel,
                      TString plottitle)
{
        TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas MuMu",800,500);
        Canvas1->SetFillColor(0);
        Canvas1->SetBorderMode(0);
        Canvas1->SetBorderSize(2);
        Canvas1->SetFrameBorderMode(0);
        mcstackhist->SetTitle(0);

        ci = TColor::GetColor("#ffff99");

        datahist->Sumw2();
        datahist->SetStats(0);
        datahist->SetLineWidth(2);
        datahist->SetMarkerStyle(20);
        mchist1->SetFillColor(800);mchist1->SetLineWidth(1);/*mchist1->SetLineColor(4);*/
        mchist2->SetFillColor(30);  /*mchist2->SetFillStyle(3001);*/
        mchist3->SetFillColor(419);
        mchist4->SetFillColor(38);
        mchist5->SetFillColor(903);
        mchist6->SetFillColor(2);
        mcstackhist->Add(mchist6);
        mcstackhist->Add(mchist3);
        mcstackhist->Add(mchist2);
        mcstackhist->Add(mchist1);
        mcstackhist->Add(mchist4);
        mcstackhist->Add(mchist5);
        if(mcstackhist->GetMaximum() > 0) mcstackhist.SetMaximum(mcstackhist->GetMaximum() * 1.5);
        else mcstackhist.SetMaximum(datahist->GetMaximum() * 2.0);
        mcstackhist->Draw();
        mcstackhist->GetXaxis()->SetTitle(xaxislabel);
        mcstackhist->GetYaxis()->SetTitle(yaxislabel);
        datahist->Draw("same");
        Canvas1->SaveAs(plottitle);

}



bool PassesTrigger(int triggerbit, int runNumber)
{
	bool pass = false;
	if(triggerbit == 1)
		pass = true;
/*	if(runNumber>=147146)
                pass = true;*/

	return pass;
}

bool PassesCosmicsCut()
{
}

bool PassesDphiCut(double deltaphioverpi)
{
	bool pass = false;
  	float dphicut = 0.9;     

	if(deltaphioverpi > dphicut)
		pass = true;
	return pass;
}

bool PassesDptCut(double deltapt)
{
	bool pass = false;
  	float dptcut = 1.0;
	
	if(deltapt < dptcut)
		pass = true;
	return pass;
}

bool PassesTowerCountVeto(int nEB, int nEE, int nHB, int nHE, int nHFp, int nHFm)
{
	bool pass = false;
	int ntowerthresh = 99999; //5;
	if((nEB+nEE+nHB+nHE+nHFp+nHFm) < ntowerthresh)
		pass = true;

	return pass;

}
 
bool FailsTrackDistanceVeto(double trackdistance, int trackquality)
{
	bool fail = false;
	float vtxtrackcountingcut = 0.2; 

	if(trackdistance < vtxtrackcountingcut)
		fail = true;

//	if((trackdistance < vtxtrackcountingcut) && (trackquality == 1))
//		fail = true;

	return fail;
}

double VertexSeparation(int nvtx, int* vtxtrks, double* vtxz, int* ismumuvtx)
{
	double closestvtx = 9999999.;
	if(nvtx == 1)
		return 9999999.;
	for(Int_t i = 0; i < nvtx; i++)
	{
		if(ismumuvtx[i] == 0)
			continue;
		for(Int_t j = 0; j < nvtx && (j!=i); j++)
		{
			if(fabs(vtxz[i]-vtxz[j]) < fabs(closestvtx))
				closestvtx = vtxz[i]-vtxz[j];
		} 
		if(i==0){
	                for(Int_t j = 1; j < nvtx && (j!=i); j++)
	                {
	                        if(fabs(vtxz[i]-vtxz[j]) < fabs(closestvtx))
	                                closestvtx = vtxz[i]-vtxz[j];
	                }
		}
	}
	return closestvtx;
}  

bool PassesZDCVeto(float em1, float em2, float had1, float had2)
{
	bool pass = false;
  	float ZDChadThresh =  999999999120.; //120.0;
  	float ZDCemThresh =   999999999916.; // 16.0;

	// Veto
	if(((em1<ZDCemThresh) && (em2<ZDCemThresh) && (had1<ZDChadThresh) && (had2<ZDChadThresh)))
		pass = true;

	// Anti-veto
//	if(((em1>ZDCemThresh) || (em2>ZDCemThresh) || (had1>ZDChadThresh) || (had2>ZDChadThresh)))
//		pass = true;

	// Anti-veto with double tags
	//	if((((em1>ZDCemThresh) || (had1>ZDChadThresh)) && ((em2>ZDCemThresh) || (had2>ZDChadThresh))))
	//		pass = true;

	// Anti-veto with single tags only
	//	if((((em1>ZDCemThresh) || (had1>ZDChadThresh)) && (em2<ZDCemThresh) && (had2<ZDChadThresh)) || 
	//	  (((em2>ZDCemThresh) || (had2>ZDChadThresh)) && (em1<ZDCemThresh) && (had1<ZDChadThresh)))
	//		pass = true;

	return pass;
}

bool PassesKinematicCuts(float mass, float ptPlus, float ptMinus, float etaPlus, float etaMinus)
{
	bool pass = true;
  	float lowermasscut1 = 0.0;
  	float uppermasscut1 = 8.5;
  	float lowermasscut2 = 11.5;
  	float uppermasscut2 = 999.0; 

        if((mass < lowermasscut1) || (mass > uppermasscut2)) 
		pass = false;
        if((mass > uppermasscut1) && (mass < lowermasscut2))  
		pass = false;

	if(ptPlus<3.0 || ptMinus<3.0)
		pass = false;

        if(fabs(etaPlus)>2.4 || fabs(etaMinus)>2.4)
                pass = false;

	return pass;
}

bool PassesMuonID(int trackerMuon1, int selectorMuon1, int globalMuon1, int trackerMuon2, int selectorMuon2, int globalMuon2, int nHitsMuon1, int nHitsMuon2)
{
	bool pass = false;
	int nhitsthresh = 10;

	bool muhit1pass = (nHitsMuon1 > nhitsthresh);
	bool muhit2pass = (nHitsMuon2 > nhitsthresh);
	bool muID1pass = (trackerMuon1 && selectorMuon1) || globalMuon1;
	bool muID2pass = (trackerMuon2 && selectorMuon2) || globalMuon2;  

	if(muhit1pass && muhit2pass && muID1pass && muID2pass)
	{
		pass = true;
	}
	return pass;
}

bool PassesVertexSelection(int vtxNTrack,double vtxChi2,double vtxNdf,double distance_vertex_z,double vtxZ,int isdimuonvtx) 
{
	bool pass = false;
 	float vtxseparationzthresh = 0.2; 
  	float PrimVertexZcut = 24.0; 
        if((vtxNTrack==2) 
	   && (isdimuonvtx == 1)	
	   && (TMath::Prob(vtxChi2,vtxNdf+0.5)>0.001) 
	   && (fabs(distance_vertex_z) > vtxseparationzthresh) 
	   && (fabs(vtxZ)<PrimVertexZcut)) 
		pass = true; 
	return pass;
}

void Cumul()
{

gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gROOT->SetStyle("Plain");
gStyle->SetStatStyle(0);
gROOT->SetTitle(0);
#define pi 3.14159265359

//definition des fichiers + Tree
  TFile *f0 = new TFile("cand_2tracks.root"); // 
//  TFile *f0 = new TFile("cand_2tracks.root"); // 
  TTree *t0 = f0->Get("ntp1");
  TFile *f1 = new TFile("El-El_highPt.root"); //
  TTree *t1 = f1->Get("ntp1");
  TFile *f2 = new TFile("Inel-El_highPt.root"); //
  TTree *t2 = f2->Get("ntp1");
  TFile *f3 = new TFile("Inel-Inel_highPt.root"); //
  TTree *t3 = f3->Get("ntp1");
  TFile *f4 = new TFile("Upsilon.root"); //
  TTree *t4 = f4->Get("ntp1");
  TFile *f5 = new TFile("Jpsi.root"); //
  TTree *t5 = f5->Get("ntp1");

  TFile *f6 = new TFile("DY_M2to10_10to120_20.root"); //
  TTree *t6 = f6->Get("ntp1");
//  TFile *f5 = new TFile("../Jpsi.root"); //
//  TTree *t5 = f5->Get("ntp1");

//  TFile *f6 = new TFile("Studies/bkg-2tracks.root"); //
//  TTree *t6 = f6->Get("ntp1");
  
// definitions : example for mass:
//  J/psi (5) -> ||    /\
//               ||   |  |  /\    /\  <<----Upsilon resonnances (4)
//               ||   |  | |  |  |  |  
//              |--------------------|
//              |    El-El (1)       |
//              |--------------------|
//              |    Inel-El (2)     |
//              |--------------------|
//              |    Inel-Inel (3)   |
//              |--------------------|
//		|    Inclusive (6)   -------|
//		|---------------------------|
//					+ Data (0)

  TH1F* hEB0 = new TH1F("neb_data","",16,-1.,15.);
  TH1F* hEE0 = new TH1F("nee_data","",16,-1.,15.);
  TH1F* hHB0 = new TH1F("nhb_data","",16,-1.,15.);
  TH1F* hHE0 = new TH1F("nhe_data","",16,-1.,15.);
  TH1F* hHFp0 = new TH1F("nhfp_data","",16,-1.,15.);
  TH1F* hHFm0 = new TH1F("ncastor_data","",16,-1.,15.);
  TH1F* nTower0 = new TH1F("nTower_data","",31,-1.,30.);
  TH1F* hEB1 = new TH1F("neb_cumulElEl","",16,-1.,15.);
  TH1F* hEE1 = new TH1F("nee_cumulElEl","",16,-1.,15.);
  TH1F* hHB1 = new TH1F("nhb_cumulElEl","",16,-1.,15.);
  TH1F* hHE1 = new TH1F("nhe_cumulElEl","",16,-1.,15.);
  TH1F* hHFp1 = new TH1F("nhfp_cumulElEl","",16,-1.,15.);
  TH1F* hHFm1 = new TH1F("ncastor_cumulElEl","",16,-1.,15.);
  TH1F* nTower1 = new TH1F("nTower_cumulElEl","",31,-1.,30.);
  TH1F* hEB2 = new TH1F("neb_cumulInelEl","",16,-1.,15.);
  TH1F* hEE2 = new TH1F("nee_cumulInelEl","",16,-1.,15.);
  TH1F* hHB2 = new TH1F("nhb_cumulInelEl","",16,-1.,15.);
  TH1F* hHE2 = new TH1F("nhe_cumulInelEl","",16,-1.,15.);
  TH1F* hHFp2 = new TH1F("nhfp_cumulInelEl","",16,-1.,15.);
  TH1F* hHFm2 = new TH1F("ncastor_cumulInelEl","",16,-1.,15.);
  TH1F* nTower2 = new TH1F("nTower_cumulInelEl","",31,-1.,30.);
  TH1F* hEB3 = new TH1F("neb_cumulInelInel","",16,-1.,15.);
  TH1F* hEE3 = new TH1F("nee_cumulInelInel","",16,-1.,15.);
  TH1F* hHB3 = new TH1F("nhb_cumulInelInel","",16,-1.,15.);
  TH1F* hHE3 = new TH1F("nhe_cumulInelInel","",16,-1.,15.);
  TH1F* hHFp3 = new TH1F("nhfp_cumulInelInel","",16,-1.,15.);
  TH1F* hHFm3 = new TH1F("ncastor_cumulInelInel","",16,-1.,15.);
  TH1F* nTower3 = new TH1F("nTower_cumulInelInel","",31,-1.,30.);
  TH1F* hEB4 = new TH1F("neb_cumulUps","",16,-1.,15.);
  TH1F* hEE4 = new TH1F("nee_cumulUps","",16,-1.,15.);
  TH1F* hHB4 = new TH1F("nhb_cumulUps","",16,-1.,15.);
  TH1F* hHE4 = new TH1F("nhe_cumulUps","",16,-1.,15.);
  TH1F* hHFp4 = new TH1F("nhfp_cumulUps","",16,-1.,15.);
  TH1F* hHFm4 = new TH1F("ncastor_cumulUps","",16,-1.,15.);
  TH1F* nTower4 = new TH1F("nTower_cumulUps","",31,-1.,30.);
  TH1F* hEB5 = new TH1F("neb_cumulJpsi","",16,-1.,15.);
  TH1F* hEE5 = new TH1F("nee_cumulJpsi","",16,-1.,15.);
  TH1F* hHB5 = new TH1F("nhb_cumulJpsi","",16,-1.,15.);
  TH1F* hHE5 = new TH1F("nhe_cumulJpsi","",16,-1.,15.);
  TH1F* hHFp5 = new TH1F("nhfp_cumulJpsi","",16,-1.,15.);
  TH1F* hHFm5 = new TH1F("ncastor_cumulJpsi","",16,-1.,15.);
  TH1F* nTower5 = new TH1F("nTower_cumulJpsi","",31,-1.,30.);
  TH1F* hEB6 = new TH1F("neb_cumulInclu","",16,-1.,15.);
  TH1F* hEE6 = new TH1F("nee_cumulInclu","",16,-1.,15.);
  TH1F* hHB6 = new TH1F("nhb_cumulInclu","",16,-1.,15.);
  TH1F* hHE6 = new TH1F("nhe_cumulInclu","",16,-1.,15.);
  TH1F* hHFp6 = new TH1F("nhfp_cumulInclu","",16,-1.,15.);
  TH1F* hHFm6 = new TH1F("ncastor_cumulInclu","",16,-1.,15.);
  TH1F* nTower6 = new TH1F("nTower_cumulInclu","",31,-1.,30.);
  THStack *sEB = new THStack("sEB","stack EB");
  THStack *sEE = new THStack("sEE","stack EE");
  THStack *sHB = new THStack("sHB","stack HB");
  THStack *sHE = new THStack("sHE","stack HE");
  THStack *sHFp = new THStack("sHFp","stack HFp");
  THStack *sHFm = new THStack("sHFm","stack Hfm");
  THStack *sTower = new THStack("sTower","stack nTower");

  TH1F* nTrack0 = new TH1F("ntrack_data","",15,-1.,14.);
  TH1F* nTrack1 = new TH1F("ntrack_cumulElEl","",15,-1.,14.);
  TH1F* nTrack2 = new TH1F("ntrack_cumulInelEl","",15,-1.,14.);
  TH1F* nTrack3 = new TH1F("ntrack_cumulInelInel","",15,-1.,14.);
  TH1F* nTrack4 = new TH1F("ntrack_cumulUps","",15,-1.,14.);
  TH1F* nTrack5 = new TH1F("ntrack_cumulJpsi","",15,-1.,14.);
  TH1F* nTrack6 = new TH1F("ntrack_cumulInclu","",15,-1.,14.);
  THStack *sTrack = new THStack("sTrack","stack Tracks");

  // 120.,0.,40.
  // 10.,0.,100.	
  // 26.,-4.0,100.0 
  TH1F* MuMuMass0 = new TH1F("mass_data","",36.,-1.,89.);
  TH1F* MuMuMass1 = new TH1F("mass_cumulElEl","",36.,-1.,89.);
  TH1F* MuMuMass2 = new TH1F("mass_cumulInelEl","",36.,-1.,89.);
  TH1F* MuMuMass3 = new TH1F("mass_cumulInelInel","",36.,-1.,89.);
  TH1F* MuMuMass4 = new TH1F("mass_cumulUps","",36.,-1.,89.);
  TH1F* MuMuMass5 = new TH1F("mass_cumulJpsi","",36.,-1.,89.);
  TH1F* MuMuMass6 = new TH1F("mass_cumulInclu","",36.,-1.,89.);
  THStack *sMuMuMass = new THStack("sMuMuMass","stack Mass");

  TH1F* MuMuMassUps0 = new TH1F("massUps_data","",40.,8.,12.);
  TH1F* MuMuMassUps1 = new TH1F("massUps_cumulElEl","",40.,8.,12.);
  TH1F* MuMuMassUps2 = new TH1F("massUps_cumulInelEl","",40.,8.,12.);
  TH1F* MuMuMassUps3 = new TH1F("massUps_cumulInelInel","",40.,8.,12.);
  TH1F* MuMuMassUps4 = new TH1F("massUps_cumulUps","",40.,8.,12.);
  TH1F* MuMuMassUps5 = new TH1F("massUps_cumulJpsi","",40.,8.,12.);
  TH1F* MuMuMassUps6 = new TH1F("massUps_cumulInclu","",40.,8.,12.);
  THStack *sMuMuMassUps = new THStack("sMuMuMassUps","stack MassUps");

  TH1F* MuMuMassJpsi0 = new TH1F("massJpsi_data","",40.,2.0,4.0);
  TH1F* MuMuMassJpsi1 = new TH1F("massJpsi_cumulElEl","",40.,2.0,4.0);
  TH1F* MuMuMassJpsi2 = new TH1F("massJpsi_cumulInelEl","",40.,2.0,4.0);
  TH1F* MuMuMassJpsi3 = new TH1F("massJpsi_cumulInelInel","",40.,2.0,4.0);
  TH1F* MuMuMassJpsi4 = new TH1F("massJpsi_cumulUps","",40.,2.0,4.0);
  TH1F* MuMuMassJpsi5 = new TH1F("massJpsi_cumulJpsi","",40.,2.0,4.0);
  TH1F* MuMuMassJpsi6 = new TH1F("massJpsi_cumulInclu","",40.,2.0,4.0);
  THStack *sMuMuMassJpsi = new THStack("sMuMuMassJpsi","stack MassJpsi");

  // 40,-0.5,1.5
  TH1F* MuMudpt0 = new TH1F("dpt_data","",40,-0.5,1.5);
  TH1F* MuMudpt1 = new TH1F("dpt_cumulElEl","",40,-0.5,1.5);
  TH1F* MuMudpt2 = new TH1F("dpt_cumulInelEl","",40,-0.5,1.5);
  TH1F* MuMudpt3 = new TH1F("dpt_cumulInelInel","",40,-0.5,1.5);
  TH1F* MuMudpt4 = new TH1F("dpt_cumulUps","",40,-0.5,1.5);
  TH1F* MuMudpt5 = new TH1F("dpt_cumulJpsi","",40,-0.5,1.5);
  TH1F* MuMudpt6 = new TH1F("dpt_cumulInclu","",40,-0.5,1.5);
  THStack *sMuMudpt = new THStack("sMuMudpt","stack dpt");

  // 60,-0.1,1.1
  TH1F* MuMudphi0 = new TH1F("dphi_data","",14,0.88,1.02);
  TH1F* MuMudphi1 = new TH1F("dphi_cumulElEl","",14,0.88,1.02);
  TH1F* MuMudphi2 = new TH1F("dphi_cumulInelEl","",14,0.88,1.02);
  TH1F* MuMudphi3 = new TH1F("dphi_cumulInelInel","",14,0.88,1.02);
  TH1F* MuMudphi4 = new TH1F("dphi_cumulUps","",14,0.88,1.02);
  TH1F* MuMudphi5 = new TH1F("dphi_cumulJpsi","",14,0.88,1.02);
  TH1F* MuMudphi6 = new TH1F("dphi_cumulInclu","",14,0.88,1.02);
  THStack *sMuMudphi = new THStack("sMuMudphi","stack dphi");

  TH1F* MuMuSymdphi0 = new TH1F("Symdphi_data","",40,-0.1,0.1); 
  TH1F* MuMuSymdphi1 = new TH1F("Symdphi_cumulElEl","",40,-0.1,0.1); 
  TH1F* MuMuSymdphi2 = new TH1F("Symdphi_cumulInelEl","",40,-0.1,0.1); 
  TH1F* MuMuSymdphi3 = new TH1F("Symdphi_cumulInelInel","",40,-0.1,0.1); 
  TH1F* MuMuSymdphi4 = new TH1F("Symdphi_cumulUps","",40,-0.1,0.1); 
  TH1F* MuMuSymdphi5 = new TH1F("Symdphi_cumulJpsi","",40,-0.1,0.1); 
  TH1F* MuMuSymdphi6 = new TH1F("Symdphi_cumulInclu","",40,-0.1,0.1); 
  THStack *sMuMuSymdphi = new THStack("sMuMuSymdphi","stack Symdphi"); 


  // 28.,-1.,6.
  TH1F* MuMudeta0 = new TH1F("deta_data","",17,-0.2,3.2);
  TH1F* MuMudeta1 = new TH1F("deta_cumulElEl","",17,-0.2,3.2);
  TH1F* MuMudeta2 = new TH1F("deta_cumulInelEl","",17,-0.2,3.2);
  TH1F* MuMudeta3 = new TH1F("deta_cumulInelInel","",17,-0.2,3.2);
  TH1F* MuMudeta4 = new TH1F("deta_cumulUps","",17,-0.2,3.2);
  TH1F* MuMudeta5 = new TH1F("deta_cumulJpsi","",17,-0.2,3.2);
  TH1F* MuMudeta6 = new TH1F("deta_cumulInclu","",17,-0.2,3.2);
  THStack *sMuMudeta = new THStack("sMuMudeta","stack deta");

  TH1F* MuMuvtxXY0 = new TH1F("vtxXY_data","",30,-0.1,0.2);
  TH1F* MuMuvtxXY1 = new TH1F("vtxXY_cumulElEl","",30,-0.1,0.2);
  TH1F* MuMuvtxXY2 = new TH1F("vtxXY_cumulInelEl","",30,-0.1,0.2);
  TH1F* MuMuvtxXY3 = new TH1F("vtxXY_cumulInelInel","",30,-0.1,0.2);
  TH1F* MuMuvtxXY4 = new TH1F("vtxXY_cumulUps","",30,-0.1,0.2);
  TH1F* MuMuvtxXY5 = new TH1F("vtxXY_cumulJpsi","",30,-0.1,0.2);
  TH1F* MuMuvtxXY6 = new TH1F("vtxXY_cumulInclu","",30,-0.1,0.2);
  THStack *sMuMuvtxXY = new THStack("sMuMuvtxXY","stack vtxXY");

  TH1F* Tdist0 = new TH1F("tDist_data","",10,-0.2,1.8); 
  TH1F* Tdist1 = new TH1F("tDist_cumulElEl","",10,-0.2,1.8); 
  TH1F* Tdist2 = new TH1F("tDist_cumulInelEl","",10,-0.2,1.8); 
  TH1F* Tdist3 = new TH1F("tDist_cumulInelInel","",10,-0.2,1.8); 
  TH1F* Tdist4 = new TH1F("tDist_cumulUps","",10,-0.2,1.8); 
  TH1F* Tdist5 = new TH1F("tDist_cumulJpsi","",10,-0.2,1.8); 
  TH1F* Tdist6 = new TH1F("tDist_cumulInclu","",10,-0.2,1.8); 
  THStack *sTdist = new THStack("sTdist","stack Tdist");

  TH1F* etaPair0 = new TH1F("etaPair_data","",26,-2.6,2.6);
  TH1F* etaPair1 = new TH1F("etaPair_cumulElEl","",26,-2.6,2.6);
  TH1F* etaPair2 = new TH1F("etaPair_cumulInelEl","",26,-2.6,2.6);
  TH1F* etaPair3 = new TH1F("etaPair_cumulInelInel","",26,-2.6,2.6);
  TH1F* etaPair4 = new TH1F("etaPair_cumulUps","",26,-2.6,2.6);
  TH1F* etaPair5 = new TH1F("etaPair_cumulJpsi","",26,-2.6,2.6);
  TH1F* etaPair6 = new TH1F("etaPair_cumulInclu","",26,-2.6,2.6);
  THStack *setaPair = new THStack("setaPair","stack eta Pair");

  TH1F* pTPair0 = new TH1F("pTPair_data","",20,-1.,4.);
  TH1F* pTPair1 = new TH1F("pTPair_cumulElEl","",20,-1.,4.);
  TH1F* pTPair2 = new TH1F("pTPair_cumulInelEl","",20,-1.,4.);
  TH1F* pTPair3 = new TH1F("pTPair_cumulInelInel","",20,-1.,4.);
  TH1F* pTPair4 = new TH1F("pTPair_cumulUps","",20,-1.,4.);
  TH1F* pTPair5 = new TH1F("pTPair_cumulJpsi","",20,-1.,4.);
  TH1F* pTPair6 = new TH1F("pTPair_cumulInclu","",20,-1.,4.);
  THStack *spTPair = new THStack("spTPair","stack pT Pair");


  TH1F* phiSingleP0 = new TH1F("phiSingleP_data","",12,-1.2,1.2);
  TH1F* phiSingleP1 = new TH1F("phiSingleP_cumulElEl","",12,-1.2,1.2);
  TH1F* phiSingleP2 = new TH1F("phiSingleP_cumulInelEl","",12,-1.2,1.2);
  TH1F* phiSingleP3 = new TH1F("phiSingleP_cumulInelInel","",12,-1.2,1.2);
  TH1F* phiSingleP4 = new TH1F("phiSingleP_cumulUps","",12,-1.2,1.2);
  TH1F* phiSingleP5 = new TH1F("phiSingleP_cumulJpsi","",12,-1.2,1.2);
  TH1F* phiSingleP6 = new TH1F("phiSingleP_cumulInclu","",12,-1.2,1.2);
  THStack *sphiSingleP = new THStack("sphiSingleP","stack eta Single +");

  TH1F* phiSingleM0 = new TH1F("phiSingleM_data","",12,-1.2,1.2);
  TH1F* phiSingleM1 = new TH1F("phiSingleM_cumulElEl","",12,-1.2,1.2);
  TH1F* phiSingleM2 = new TH1F("phiSingleM_cumulInelEl","",12,-1.2,1.2);
  TH1F* phiSingleM3 = new TH1F("phiSingleM_cumulInelInel","",12,-1.2,1.2);
  TH1F* phiSingleM4 = new TH1F("phiSingleM_cumulUps","",12,-1.2,1.2);
  TH1F* phiSingleM5 = new TH1F("phiSingleM_cumulJpsi","",12,-1.2,1.2);
  TH1F* phiSingleM6 = new TH1F("phiSingleM_cumulInclu","",12,-1.2,1.2);
  THStack *sphiSingleM = new THStack("sphiSingleM","stack eta Single -");


  TH1F* etaSingleP0 = new TH1F("etaSingleP_data","",26,-2.6,2.6);
  TH1F* etaSingleP1 = new TH1F("etaSingleP_cumulElEl","",26,-2.6,2.6);
  TH1F* etaSingleP2 = new TH1F("etaSingleP_cumulInelEl","",26,-2.6,2.6);
  TH1F* etaSingleP3 = new TH1F("etaSingleP_cumulInelInel","",26,-2.6,2.6);
  TH1F* etaSingleP4 = new TH1F("etaSingleP_cumulUps","",26,-2.6,2.6);
  TH1F* etaSingleP5 = new TH1F("etaSingleP_cumulJpsi","",26,-2.6,2.6);
  TH1F* etaSingleP6 = new TH1F("etaSingleP_cumulInclu","",26,-2.6,2.6);
  THStack *setaSingleP = new THStack("setaSingleP","stack eta SingleP");

  TH1F* etaSingleM0 = new TH1F("etaSingleM_data","",26,-2.6,2.6);
  TH1F* etaSingleM1 = new TH1F("etaSingleM_cumulElEl","",26,-2.6,2.6);
  TH1F* etaSingleM2 = new TH1F("etaSingleM_cumulInelEl","",26,-2.6,2.6);
  TH1F* etaSingleM3 = new TH1F("etaSingleM_cumulInelInel","",26,-2.6,2.6);
  TH1F* etaSingleM4 = new TH1F("etaSingleM_cumulUps","",26,-2.6,2.6);
  TH1F* etaSingleM5 = new TH1F("etaSingleM_cumulJpsi","",26,-2.6,2.6);
  TH1F* etaSingleM6 = new TH1F("etaSingleM_cumulInclu","",26,-2.6,2.6);
  THStack *setaSingleM = new THStack("setaSingleM","stack eta SingleM");


  TH1F* pTSingleP0 = new TH1F("pTSingleP_data","",41,-1.,40.);
  TH1F* pTSingleP1 = new TH1F("pTSingleP_cumulElEl","",41,-1.,40.);
  TH1F* pTSingleP2 = new TH1F("pTSingleP_cumulInelEl","",41,-1.,40.);
  TH1F* pTSingleP3 = new TH1F("pTSingleP_cumulInelInel","",41,-1.,40.);
  TH1F* pTSingleP4 = new TH1F("pTSingleP_cumulUps","",41,-1.,40.);
  TH1F* pTSingleP5 = new TH1F("pTSingleP_cumulJpsi","",41,-1.,40.);
  TH1F* pTSingleP6 = new TH1F("pTSingleP_cumulInclu","",41,-1.,40.);
  THStack *spTSingleP = new THStack("spTSingleP","stack pT muon +");

  TH1F* pTSingleM0 = new TH1F("pTSingleM_data","",41,-1.,40.);
  TH1F* pTSingleM1 = new TH1F("pTSingleM_cumulElEl","",41,-1.,40.);
  TH1F* pTSingleM2 = new TH1F("pTSingleM_cumulInelEl","",41,-1.,40.);
  TH1F* pTSingleM3 = new TH1F("pTSingleM_cumulInelInel","",41,-1.,40.);
  TH1F* pTSingleM4 = new TH1F("pTSingleM_cumulUps","",41,-1.,40.);
  TH1F* pTSingleM5 = new TH1F("pTSingleM_cumulJpsi","",41,-1.,40.);
  TH1F* pTSingleM6 = new TH1F("pTSingleM_cumulInclu","",41,-1.,40.);
  THStack *spTSingleM = new THStack("spTSingleM","stack pT muon -");


  TH1F* ZDCemplus0 = new TH1F("zdcEm+_data","",72.,-100.,3500.);
  TH1F* ZDCemplus1 = new TH1F("zdcEm+_cumulElEl","",72.,-100.,3500.);
  TH1F* ZDCemplus2 = new TH1F("zdcEm+_cumulInelEl","",72.,-100.,3500.);
  TH1F* ZDCemplus3 = new TH1F("zdcEm+_cumulInelInel","",72.,-100.,3500.);
  TH1F* ZDCemplus4 = new TH1F("zdcEm+_cumulUps","",72.,-100.,3500.);
  TH1F* ZDCemplus5 = new TH1F("zdcEm+_cumulJpsi","",72.,-100.,3500.);
  TH1F* ZDCemplus6 = new TH1F("zdcEm+_cumulInclu","",72.,-100.,3500.);
  THStack *sZDCemplus = new THStack("sZDCemplus","stack ZDC EM+");

  TH1F* ZDChadplus0 = new TH1F("zdcHad+_data","",92.,-1000.,45000.);
  TH1F* ZDChadplus1 = new TH1F("zdcHad+_cumulElEl","",92.,-1000.,45000.);
  TH1F* ZDChadplus2 = new TH1F("zdcHad+_cumulInelEl","",92.,-1000.,45000.);
  TH1F* ZDChadplus3 = new TH1F("zdcHad+_cumulInelInel","",92.,-1000.,45000.);
  TH1F* ZDChadplus4 = new TH1F("zdcHad+_cumulUps","",92.,-1000.,45000.);
  TH1F* ZDChadplus5 = new TH1F("zdcHad+_cumulJpsi","",92.,-1000.,45000.);
  TH1F* ZDChadplus6 = new TH1F("zdcHad+_cumulInclu","",92.,-1000.,45000.);
  THStack *sZDChadplus = new THStack("sZDChadplus","stack ZDC HAD+");

  TH1F* ZDCemminus0 = new TH1F("zdcEm-_data","",72.,-100.,3500.);
  TH1F* ZDCemminus1 = new TH1F("zdcEm-_cumulElEl","",72.,-100.,3500.);
  TH1F* ZDCemminus2 = new TH1F("zdcEm-_cumulInelEl","",72.,-100.,3500.);
  TH1F* ZDCemminus3 = new TH1F("zdcEm-_cumulInelInel","",72.,-100.,3500.);
  TH1F* ZDCemminus4 = new TH1F("zdcEm-_cumulUps","",72.,-100.,3500.);
  TH1F* ZDCemminus5 = new TH1F("zdcEm-_cumulJpsi","",72.,-100.,3500.);
  TH1F* ZDCemminus6 = new TH1F("zdcEm-_cumulInclu","",72.,-100.,3500.);
  THStack *sZDCemminus = new THStack("sZDCemminus","stack ZDC EM-");

  TH1F* ZDChadminus0 = new TH1F("zdcHad-_data","",92.,-1000.,45000.);
  TH1F* ZDChadminus1 = new TH1F("zdcHad-_cumulElEl","",92.,-1000.,45000.);
  TH1F* ZDChadminus2 = new TH1F("zdcHad-_cumulInelEl","",92.,-1000.,45000.);
  TH1F* ZDChadminus3 = new TH1F("zdcHad-_cumulInelInel","",92.,-1000.,45000.);
  TH1F* ZDChadminus4 = new TH1F("zdcHad-_cumulUps","",92.,-1000.,45000.);
  TH1F* ZDChadminus5 = new TH1F("zdcHad-_cumulJpsi","",92.,-1000.,45000.);
  TH1F* ZDChadminus6 = new TH1F("zdcHad-_cumulInclu","",92.,-1000.,45000.);
  THStack *sZDChadminus = new THStack("sZDChadminus","stack ZDC HAD-");

  TH1F* ZDCtime0 =  new TH1F("zdcTime_data","",80.,-20.,60.);
  TH1F* ZDCtime1 =  new TH1F("zdcTime_cumulElEl","",80.,-20.,60.);
  TH1F* ZDCtime2 =  new TH1F("zdcTime_cumulInelEl","",80.,-20.,60.);
  TH1F* ZDCtime3 =  new TH1F("zdcTime_cumulInelInel","",80.,-20.,60.);
  TH1F* ZDCtime4 =  new TH1F("zdcTime_cumulUps","",80.,-20.,60.);
  TH1F* ZDCtime5 =  new TH1F("zdcTime_cumulJpsi","",80.,-20.,60.);
  TH1F* ZDCtime6 =  new TH1F("zdcTime_cumulInclu","",80.,-20.,60.);
  THStack* sZDCtime = new THStack("sZDCtime","stack ZDC time");

  TH1F* ZDCenergyEM0 =  new TH1F("zdcEnEM_data","",100.,-500.,2500.);
  TH1F* ZDCenergyEM1 =  new TH1F("zdcEnEM_cumulElEl","",100.,-500.,2500.);
  TH1F* ZDCenergyEM2 =  new TH1F("zdcEnEM_cumulInelEl","",100.,-500.,2500.);
  TH1F* ZDCenergyEM3 =  new TH1F("zdcEnEM_cumulInelInel","",100.,-500.,2500.);
  TH1F* ZDCenergyEM4 =  new TH1F("zdcEnEM_cumulUps","",100.,-500.,2500.);
  TH1F* ZDCenergyEM5 =  new TH1F("zdcEnEM_cumulJpsi","",100.,-500.,2500.);
  TH1F* ZDCenergyEM6 =  new TH1F("zdcEnEM_cumulInclu","",100.,-500.,2500.);
  THStack* sZDCenergyEM = new THStack("sZDCenergyEM","stack ZDC energy EM");

  TH1F* ZDCenergyHAD0 =  new TH1F("zdcEnHAD_data","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD1 =  new TH1F("zdcEnHAD_cumulElEl","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD2 =  new TH1F("zdcEnHAD_cumulInelEl","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD3 =  new TH1F("zdcEnHAD_cumulInelInel","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD4 =  new TH1F("zdcEnHAD_cumulUps","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD5 =  new TH1F("zdcEnHAD_cumulJpsi","",100.,-5000.,25000.);
  TH1F* ZDCenergyHAD6 =  new TH1F("zdcEnHAD_cumulInclu","",100.,-5000.,25000.);
  THStack* sZDCenergyHAD = new THStack("sZDCenergyHAD","stack ZDC energy HAD");

  TH1F* CastorSumE0 = new TH1F("castorE_data","",140.,-500.,10300.);
  TH1F* CastorSumE1 = new TH1F("castorE_cumulElEl","",140.,-500.,10300.);
  TH1F* CastorSumE2 = new TH1F("castorE_cumulInelEl","",140.,-500.,10300.);
  TH1F* CastorSumE3 = new TH1F("castorE_cumulInelInel","",140.,-500.,10300.);
  TH1F* CastorSumE4 = new TH1F("castorE_cumulUps","",140.,-500.,10300.);
  TH1F* CastorSumE5 = new TH1F("castorE_cumulJpsi","",140.,-500.,10300.);
  TH1F* CastorSumE6 = new TH1F("castorE_cumulInclu","",140.,-500.,10300.);
  THStack* sCastorSumE = new THStack("sCastorSumE","stack Castor energy");

  TH1F* VtxT0 = new TH1F("vtx_data","",300,0.,0.6);
  TH1F* VtxT1 = new TH1F("vtx_cumulElEl","",300,0.,0.6);
  TH1F* VtxT2 = new TH1F("vtx_cumulInelEl","",300,0.,0.6);
  TH1F* VtxT3 = new TH1F("vtx_cumulInelInel","",300,0.,0.6);
  TH1F* VtxT4 = new TH1F("vtx_cumulUps","",300,0.,0.6);
  TH1F* VtxT5 = new TH1F("vtx_cumulJpsi","",300,0.,0.6);
  TH1F* VtxT6 = new TH1F("vtx_cumulInclu","",300,0.,0.6);
  THStack* sVtxT = new THStack("sVtxT","stack vtx Transverse");


// definitions des # d'entrées
  const int NUM0 = t0->GetEntries();
  const int NUM1 = t1->GetEntries();
  const int NUM2 = t2->GetEntries();
  const int NUM3 = t3->GetEntries();
  const int NUM4 = t4->GetEntries();
  const int NUM5 = t5->GetEntries();
  const int NUM6 = t6->GetEntries();

//  const float integrated_lumi = 299.30569*0.5857; //in nb-1
//  const float integrated_lumi = 2872.246*0.5; // in nb-1	
//  const float integrated_lumi = 2872.246;
  const float integrated_lumi = 34686.819159*0.84;
//  const float integrated_lumi = 3044.0 * 0.5;
//  const float integrated_lumi = 4426.5;
//  const float integrated_lumi = 4426.5 - 3044.0;

  const float doublemuopenfractionallumi = 1.0;
  const float doublemuopentightfractionallumi = 0.0;
  const float HBThresh = 1.25; 
  const float HEThresh = 1.80; 
  const float EEThresh = 2.40; 
  const float EBThresh = 0.6; 
  const float HFPlusThresh = 4.5; 
  const float HFMinusThresh = 4.0; 
  const float ZDChadThresh = 120.0; 
  const float ZDCemThresh = 16.0;                
  const float dRcone = 0.3;

  const float fac_lumi0 = 1.0;
  const float fac_lumi1 = 1.0850e-6*integrated_lumi;			//  1.0820e-6*integrated_lumi;
  const float fac_lumi2 = 3.0525e-6*integrated_lumi;			//  3.05250e-6*integrated_lumi;
  const float fac_lumi3 = 4.7400e-6*integrated_lumi;			//  4.740e-6*integrated_lumi;
  const float fac_lumi4 = 1.350e-6*integrated_lumi;			//  1.350e-6*integrated_lumi;
  const float fac_lumi5 = 3.02430e-4*integrated_lumi;

  const float fac_lumiBkg[21]={0.,1.,6.7736e-5*integrated_lumi,3.,4.,5.,6.,7.,8.,9.,1.194e-6*integrated_lumi,11.,12.,13.,14.,15.,16.,17.,18.,19.,5.68e-7*integrated_lumi};
//  for(Int_t i=0; i<=20; i++){
//	cout<<"fac("<<i<<")="<<fac_lumiBkg[i]<<endl;
//  }
//definition des variables
  Int_t hlt_d1[1], hlt_d2[1], hlt_d0[1], hlt_d3[1], hlt_d4[1], hlt_d5[1],  hlt_d6[1], hlt_d0bis[1];
  Int_t techBit1[1][128], techBit2[1][128], techBit0[1][128], techBit3[1][128], techBit4[1][128], techBit5[1][128];
// MuonID
  Int_t var_idA1[10], var_idA2[10], var_idA0[10],var_idA3[10],var_idA4[10],var_idA5[10],var_idA6[10],var_idB1[10], var_idB2[10], var_idB0[10],var_idB3[10],var_idB4[10],var_idB5[10],var_idB6[10],var_idC1[10], var_idC2[10],var_idC0[10],var_idC3[10],var_idC4[10],var_idC5[10],var_idC6[10],var_idD1[10], var_idD2[10], var_idD0[10], var_idD3[10],var_idD4[10],var_idD5[10],var_idD6[10],var_idE1[10], var_idE2[10], var_idE0[10] , var_idE3[10], var_idE4[10], var_idE5[10], var_idE6[10];
  Int_t var_nMuon[1];
// RecoTrack
  Int_t var_nTrack1[1], var_nTrack2[1], var_nTrack0[1],var_nTrack3[1],var_nTrack4[1],var_nTrack5[1],var_nTrack6[1];
  Int_t var_nTrackQual1[1], var_nTrackQual2[1], var_nTrackQual0[1], var_nTrackQual3[1], var_nTrackQual4[1], var_nTrackQual5[1], var_nTrackQual6[1];
  Double_t var_TrackPt1[2000], var_TrackPt2[2000], var_TrackPt0[2000], var_TrackPt3[2000], var_TrackPt4[2000], var_TrackPt5[2000], var_TrackPt6[2000];
  Double_t var_TrackD1[2000],var_TrackD2[2000], var_TrackD0[2000], var_TrackD3[2000], var_TrackD4[2000], var_TrackD5[2000], var_TrackD6[2000];
  Double_t var_TrackQuality1[2000],var_TrackQuality2[2000], var_TrackQuality0[2000], var_TrackQuality3[2000], var_TrackQuality4[2000], var_TrackQuality5[2000], var_TrackQuality6[2000];;
// Prim Vtx
  Int_t var_nvtx1[1],var_nvtx2[1], var_nvtx0[1], var_nvtx3[1], var_nvtx4[1], var_nvtx5[1], var_nvtx6[1];
  Int_t var_vtxTrack1[10],var_vtxTrack2[10], var_vtxTrack0[10],  var_vtxTrack3[10],  var_vtxTrack4[10], var_vtxTrack5[10], var_vtxTrack6[10];
  Int_t var_vtxmumu1[10],var_vtxmumu2[10], var_vtxmumu0[10],  var_vtxmumu3[10],  var_vtxmumu4[10], var_vtxmumu5[10], var_vtxmumu6[10]; 

  Double_t var_vtxZ1[10],var_vtxZ2[10], var_vtxZ0[10], var_vtxZ3[10], var_vtxZ4[10], var_vtxZ5[10], var_vtxZ6[10];
  Double_t var_vtxX1[10],var_vtxX2[10], var_vtxX0[10], var_vtxX3[10], var_vtxX4[10], var_vtxX5[10], var_vtxX6[10];
  Double_t var_vtxY1[10],var_vtxY2[10], var_vtxY0[10], var_vtxY3[10], var_vtxY4[10], var_vtxY5[10], var_vtxY6[10];
  Double_t var_MuMuvtxX1[10],var_MuMuvtxX2[10], var_MuMuvtxX0[10], var_MuMuvtxX3[10], var_MuMuvtxX4[10], var_MuMuvtxX5[10], var_MuMuvtxX6[10];
  Double_t var_MuMuvtxY1[10],var_MuMuvtxY2[10], var_MuMuvtxY0[10], var_MuMuvtxY3[10], var_MuMuvtxY4[10], var_MuMuvtxY5[10], var_MuMuvtxY6[10];
  Double_t var_MuMuvtxZ1[10],var_MuMuvtxZ2[10], var_MuMuvtxZ0[10], var_MuMuvtxZ3[10], var_MuMuvtxZ4[10], var_MuMuvtxZ5[10], var_MuMuvtxZ6[10];;
  Int_t var_MuMuvtxValid1[10],var_MuMuvtxValid2[10], var_MuMuvtxValid0[10], var_MuMuvtxValid3[10], var_MuMuvtxValid4[10], var_MuMuvtxValid5[10], var_MuMuvtxValid6[10];
  Double_t var_vertexChi2_1[10],var_vertexChi2_2[10], var_vertexChi2_0[10], var_vertexChi2_3[10], var_vertexChi2_4[10], var_vertexChi2_5[10], var_vertexChi2_6[10];
  Double_t var_vertexNdf1[10], var_vertexNdf2[10], var_vertexNdf0[10],  var_vertexNdf3[10], var_vertexNdf4[10], var_vertexNdf5[10], var_vertexNdf6[10];

// CaloTowers 
  Int_t var_ncalo1[1], var_ncalo2[1], var_ncalo0[1], var_ncalo3[1], var_ncalo4[1], var_ncalo5[1], var_ncalo6[1];
  Int_t var_tower1[1], var_tower2[1], var_tower0[1], var_tower3[1], var_tower4[1], var_tower5[1], var_tower6[1];;
  Int_t var_caloId1[2000],var_caloId2[2000], var_caloId0[2000], var_caloId3[2000], var_caloId4[2000], var_caloId5[2000], var_caloId6[2000];
  Double_t var_caloEn1[2000],var_caloEn2[2000], var_caloEn0[2000],  var_caloEn3[2000],  var_caloEn4[2000],  var_caloEn5[2000],  var_caloEn6[2000];
  Double_t var_caloEmE1[2000],var_caloEmE2[2000], var_caloEmE0[2000], var_caloEmE3[2000],var_caloEmE4[2000],var_caloEmE5[2000], var_caloEmE6[2000];
  Double_t var_caloHadE1[2000],var_caloHadE2[2000], var_caloHadE0[2000], var_caloHadE3[2000],var_caloHadE4[2000],var_caloHadE5[2000], var_caloHadE6[2000];
  Double_t var_caloTime1[2000], var_caloTime2[2000], var_caloTime0[2000], var_caloTime3[2000], var_caloTime4[2000], var_caloTime5[2000], var_caloTime6[2000];
  Double_t var_caloEta0[2000],var_caloEta1[2000],var_caloEta2[2000],var_caloEta3[2000],var_caloEta4[2000],var_caloEta5[2000],var_caloEta6[2000];
  Double_t var_caloPhi0[2000],var_caloPhi1[2000],var_caloPhi2[2000],var_caloPhi3[2000],var_caloPhi4[2000],var_caloPhi5[2000],var_caloPhi6[2000];
  Double_t var_etmiss1[1], var_etmiss2[1], var_etmiss0[1], var_etmiss3[1], var_etmiss4[1], var_etmiss5[1], var_etmiss6[1];
  Double_t var_calodR1[2000], var_calodR2[2000], var_calodR0[2000], var_calodR3[2000], var_calodR4[2000], var_calodR5[2000], var_calodR6[2000];
  Double_t var_caloZ1[2000], var_caloZ2[2000], var_caloZ0[2000], var_caloZ3[2000], var_caloZ4[2000], var_caloZ5[2000], var_caloZ6[2000];

// Bunch crossing
  Int_t var_bx0[1], var_run0[1], var_run6[1], var_ls0[1], var_event0[1], var_event3[1], var_event4[1], var_event5[1], var_event6[1];
// MuMu kinematics
  Double_t var_mass1[5], var_mass2[5] ,var_mass0[2], var_mass3[2], var_mass4[2], var_mass5[2], var_mass6[2];
  Double_t var_dpt1[5], var_dpt2[5], var_dpt0[5], var_dpt3[5], var_dpt4[5], var_dpt5[5], var_dpt6[5];
  Double_t var_dphi1[5], var_dphi2[5], var_dphi0[5], var_dphi3[5], var_dphi4[5], var_dphi5[5], var_dphi6[5];
  Int_t var_global1[10], var_global2[10], var_global0[10], var_global3[10], var_global4[10], var_global5[10],var_global6[10],var_tracker1[10], var_tracker2[10],var_tracker0[10],var_tracker3[10],var_tracker4[10],var_tracker5[10],var_tracker6[10],var_standalone1[10], var_standalone2[10], var_standalone0[10], var_standalone3[10],var_standalone4[10],var_standalone5[10],var_standalone6[10];
  Double_t var_pt1[10], var_pt2[10], var_pt0[10],var_pt3[10],var_pt4[10],var_pt5[10],var_pt6[10],var_pz1[10], var_pz2[10], var_pz0[10],var_pz3[10],var_pz4[10],var_pz5[10],var_pz6[10],var_phi1[10], var_phi2[10], var_phi0[10],var_phi3[10], var_phi4[10],var_phi5[10],var_phi6[10],var_eta1[10], var_eta2[10], var_eta0[10], var_eta3[10], var_eta4[10], var_eta5[10], var_eta6[10];
  Int_t var_nhitsTrack1[10], var_nhitsTrack2[10], var_nhitsTrack0[10],  var_nhitsTrack3[10], var_nhitsTrack4[10], var_nhitsTrack5[10], var_nhitsTrack6[10];;
  Int_t var_Pair0[2], var_Pair1[2],var_Pair2[2],var_Pair3[2],var_Pair4[2],var_Pair5[2],var_Pair6[2];
  Double_t var_p0[10],var_p1[10],var_p2[10],var_p3[10],var_p4[10],var_p5[10],var_p6[10];
  Double_t var_px0[10],var_px1[10],var_px2[10],var_px3[10],var_px4[10],var_px5[10],var_px6[10];
  Double_t var_py0[10],var_py1[10],var_py2[10],var_py3[10],var_py4[10],var_py5[10],var_py6[10];
  Double_t var_pz0[10],var_pz1[10],var_pz2[10],var_pz3[10],var_pz4[10],var_pz5[10],var_pz6[10];
  Int_t var_charge0[10],var_charge1[10],var_charge2[10],var_charge3[10],var_charge4[10],var_charge5[10],var_charge6[10];

// ZDC
  Int_t var_nZDC1[1], var_nZDC2[1], var_nZDC0[1], var_nZDC3[1], var_nZDC4[1], var_nZDC5[1],var_nZDC6[1];
  Int_t var_zdcsection1[5000],var_zdcsection2[5000],var_zdcsection0[5000],var_zdcsection3[5000],var_zdcsection4[5000],var_zdcsection5[5000],var_zdcsection6[5000];
  Double_t var_zdcE1[5000], var_zdcE2[5000], var_zdcE0[5000], var_zdcE3[5000], var_zdcE4[5000],var_zdcE5[5000],var_zdcE6[5000];
  Double_t var_zdcEmMinus1[1], var_zdcEmMinus2[1],var_zdcEmMinus0[1],var_zdcEmMinus3[1], var_zdcEmMinus4[1],var_zdcEmMinus5[1],var_zdcEmMinus6[1],var_zdcHadMinus1[1], var_zdcHadMinus2[1], var_zdcHadMinus0[1], var_zdcHadMinus3[1], var_zdcHadMinus4[1],var_zdcHadMinus5[1],var_zdcHadMinus6[1];
  Double_t var_zdcEmPlus1[1], var_zdcEmPlus2[1], var_zdcEmPlus0[1], var_zdcEmPlus3[1],var_zdcEmPlus4[1],var_zdcEmPlus5[1],var_zdcEmPlus6[1],var_zdcHadPlus1[1], var_zdcHadPlus2[1], var_zdcHadPlus0[1], var_zdcHadPlus3[1], var_zdcHadPlus4[1], var_zdcHadPlus5[1], var_zdcHadPlus6[1];
  Double_t var_zdcTime1[5000], var_zdcTime2[5000], var_zdcTime0[5000], var_zdcTime3[5000], var_zdcTime4[5000], var_zdcTime5[5000], var_zdcTime6[5000];

//Castor
  Int_t var_nCastor1[1],var_nCastor2[1], var_nCastor0[1], var_nCastor3[1], var_nCastor4[1], var_nCastor5[1], var_nCastor6[1];
  Double_t var_CastorE1[1000],  var_CastorE2[1000], var_CastorE0[1000], var_CastorE3[1000], var_CastorE4[1000], var_CastorE5[1000], var_CastorE6[1000];
  Double_t var_CastorEta1[1000],  var_CastorEta2[1000], var_CastorEta0[1000], var_CastorEta3[1000], var_CastorEta4[1000], var_CastorEta5[1000], var_CastorEta6[1000];
  Double_t var_CastorPhi1[1000],  var_CastorPhi2[1000], var_CastorPhi0[1000], var_CastorPhi3[1000], var_CastorPhi4[1000], var_CastorPhi5[1000], var_CastorPhi6[1000];
  Double_t var_CastorRecHit1[1], var_CastorRecHit2[1], var_CastorRecHit0[1], var_CastorRecHit3[1], var_CastorRecHit4[1], var_CastorRecHit5[1], var_CastorRecHit6[1];

//Efficiency
  Double_t var_eff1[10],var_eff2[10],var_eff3[10],var_eff4[10],var_eff5[10],var_eff6[10];

  TString hlttrigger0 = "HLT_DoubleMu0";
  TString hlttrigger = "HLT_L1DoubleMuOpen";	
  TString hlttrigger2 = "HLT_L1DoubleMuOpen_Tight";
  TString hlttrigger3 = "HLT_DoubleMu3";

  t1->SetBranchAddress(hlttrigger3,hlt_d1);
  t2->SetBranchAddress(hlttrigger3,hlt_d2);
  t0->SetBranchAddress(hlttrigger3,hlt_d0);
  t3->SetBranchAddress(hlttrigger3,hlt_d3);
  t4->SetBranchAddress(hlttrigger3,hlt_d4);
  t5->SetBranchAddress(hlttrigger3,hlt_d5);
  t6->SetBranchAddress(hlttrigger3,hlt_d6);

  t0->SetBranchAddress("L1TechnicalTriggers",techBit0);
  t1->SetBranchAddress("L1TechnicalTriggers",techBit1);
  t2->SetBranchAddress("L1TechnicalTriggers",techBit2);
  t3->SetBranchAddress("L1TechnicalTriggers",techBit3);
  t4->SetBranchAddress("L1TechnicalTriggers",techBit4);
  t5->SetBranchAddress("L1TechnicalTriggers",techBit5);

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
  t5->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA5);
  t5->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB5);
  t5->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC5);
  t5->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD5);
  t5->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE5);
  t6->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA6);
  t6->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB6);
  t6->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC6);
  t6->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD6);
  t6->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE6);


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
  t4->SetBranchAddress("nTrackCand",var_nTrack4);
  t4->SetBranchAddress("nQualityTrackCand",var_nTrackQual4);
  t4->SetBranchAddress("TrackCand_pt",var_TrackPt4);
  t4->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD4);
  t4->SetBranchAddress("TrackCand_purity",var_TrackQuality4);
  t5->SetBranchAddress("nTrackCand",var_nTrack5);
  t5->SetBranchAddress("nQualityTrackCand",var_nTrackQual5);
  t5->SetBranchAddress("TrackCand_pt",var_TrackPt5);
  t5->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD5);
  t5->SetBranchAddress("TrackCand_purity",var_TrackQuality5);
  t6->SetBranchAddress("nTrackCand",var_nTrack6);
  t6->SetBranchAddress("nQualityTrackCand",var_nTrackQual6);
  t6->SetBranchAddress("TrackCand_pt",var_TrackPt6);
  t6->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD6);
  t6->SetBranchAddress("TrackCand_purity",var_TrackQuality6);

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
  t5->SetBranchAddress("nPrimVertexCand",var_nvtx5);
  t5->SetBranchAddress("PrimVertexCand_z",var_vtxZ5);
  t5->SetBranchAddress("PrimVertexCand_x",var_vtxX5);
  t5->SetBranchAddress("PrimVertexCand_y",var_vtxY5);
  t5->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_5);
  t5->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf5);
  t5->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack5);
  t6->SetBranchAddress("nPrimVertexCand",var_nvtx6);
  t6->SetBranchAddress("PrimVertexCand_z",var_vtxZ6);
  t6->SetBranchAddress("PrimVertexCand_x",var_vtxX6);
  t6->SetBranchAddress("PrimVertexCand_y",var_vtxY6);
  t6->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_6);
  t6->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf6);
  t6->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack6);
  t0->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu0); 
  t1->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu1);  
  t2->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu2);  
  t3->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu3);  
  t4->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu4);  
  t5->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu5);  
  t6->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu6);   

  t1->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX1);
  t2->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX2);
  t1->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY1);
  t2->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY2);
  t1->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ1);
  t2->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ2);
  t0->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX0);
  t0->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY0);
  t0->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ0);
  t0->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid0);
  t1->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid1);
  t2->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid2);
  t3->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX3);
  t3->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY3);
  t3->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ3);
  t3->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid3);
  t4->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX4);
  t4->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY4);
  t4->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ4);
  t4->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid4);
  t5->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX5);
  t5->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY5);
  t5->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ5);
  t5->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid5);
  t6->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX6);
  t6->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY6);
  t6->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ6);
  t6->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid6);

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
  t1->SetBranchAddress("CaloTower_eta",var_caloEta1);
  t1->SetBranchAddress("CaloTower_phi",var_caloPhi1);
  t2->SetBranchAddress("CaloTower_eta",var_caloEta2);
  t2->SetBranchAddress("CaloTower_phi",var_caloPhi2);
  t1->SetBranchAddress("Etmiss",var_etmiss1);
  t2->SetBranchAddress("Etmiss",var_etmiss2);
  t1->SetBranchAddress("nExtraCaloTowersE5",var_tower1);
  t2->SetBranchAddress("nExtraCaloTowersE5",var_tower2);
  t0->SetBranchAddress("nCaloCand",var_ncalo0);
  t0->SetBranchAddress("CaloTower_ID",var_caloId0);
  t0->SetBranchAddress("CaloTower_e",var_caloEn0);
  t0->SetBranchAddress("CaloTower_t",var_caloTime0);
  t0->SetBranchAddress("CaloTower_dr",var_calodR0);
  t0->SetBranchAddress("CaloTower_eta",var_caloEta0);
  t0->SetBranchAddress("CaloTower_phi",var_caloPhi0);
  t0->SetBranchAddress("Etmiss",var_etmiss0);
  t0->SetBranchAddress("nExtraCaloTowersE5",var_tower0);
  t3->SetBranchAddress("nCaloCand",var_ncalo3);
  t3->SetBranchAddress("CaloTower_ID",var_caloId3);
  t3->SetBranchAddress("CaloTower_e",var_caloEn3);
  t3->SetBranchAddress("CaloTower_t",var_caloTime3);
  t3->SetBranchAddress("CaloTower_dr",var_calodR3);
  t3->SetBranchAddress("CaloTower_eta",var_caloEta3);
  t3->SetBranchAddress("CaloTower_phi",var_caloPhi3);
  t3->SetBranchAddress("Etmiss",var_etmiss3);
  t3->SetBranchAddress("nExtraCaloTowersE5",var_tower3);
  t4->SetBranchAddress("nCaloCand",var_ncalo4);
  t4->SetBranchAddress("CaloTower_ID",var_caloId4);
  t4->SetBranchAddress("CaloTower_e",var_caloEn4);
  t4->SetBranchAddress("CaloTower_t",var_caloTime4);
  t4->SetBranchAddress("CaloTower_dr",var_calodR4);
  t4->SetBranchAddress("CaloTower_eta",var_caloEta4);
  t4->SetBranchAddress("CaloTower_phi",var_caloPhi4);
  t4->SetBranchAddress("Etmiss",var_etmiss4);
  t4->SetBranchAddress("nExtraCaloTowersE5",var_tower4);
  t5->SetBranchAddress("nCaloCand",var_ncalo5);
  t5->SetBranchAddress("CaloTower_ID",var_caloId5);
  t5->SetBranchAddress("CaloTower_e",var_caloEn5);
  t5->SetBranchAddress("CaloTower_t",var_caloTime5);
  t5->SetBranchAddress("CaloTower_dr",var_calodR5);
  t5->SetBranchAddress("CaloTower_eta",var_caloEta5);
  t5->SetBranchAddress("CaloTower_phi",var_caloPhi5);
  t5->SetBranchAddress("Etmiss",var_etmiss5);
  t5->SetBranchAddress("nExtraCaloTowersE5",var_tower5);
  t6->SetBranchAddress("nCaloCand",var_ncalo6);
  t6->SetBranchAddress("CaloTower_ID",var_caloId6);
  t6->SetBranchAddress("CaloTower_e",var_caloEn6);
  t6->SetBranchAddress("CaloTower_t",var_caloTime6);
  t6->SetBranchAddress("CaloTower_dr",var_calodR6);
  t6->SetBranchAddress("CaloTower_eta",var_caloEta6);
  t6->SetBranchAddress("CaloTower_phi",var_caloPhi6);
  t6->SetBranchAddress("Etmiss",var_etmiss6);
  t6->SetBranchAddress("nExtraCaloTowersE5",var_tower6);


  t0->SetBranchAddress("BX",var_bx0);
  t0->SetBranchAddress("Run",var_run0);
  t6->SetBranchAddress("Run",var_run6);
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
  t5->SetBranchAddress("MuMu_mass",var_mass5);
  t5->SetBranchAddress("MuMu_dpt",var_dpt5);
  t5->SetBranchAddress("MuMu_dphi",var_dphi5);
  t6->SetBranchAddress("MuMu_mass",var_mass6);
  t6->SetBranchAddress("MuMu_dpt",var_dpt6);
  t6->SetBranchAddress("MuMu_dphi",var_dphi6);

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
  t5->SetBranchAddress("MuonCand_isglobal",var_global5);
  t5->SetBranchAddress("MuonCand_istracker",var_tracker5);
  t5->SetBranchAddress("MuonCand_isstandalone",var_standalone5);
  t5->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack5);
  t6->SetBranchAddress("MuonCand_isglobal",var_global6);
  t6->SetBranchAddress("MuonCand_istracker",var_tracker6);
  t6->SetBranchAddress("MuonCand_isstandalone",var_standalone6);
  t6->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack6);
  t0->SetBranchAddress("MuonCand_charge",var_charge0);
  t1->SetBranchAddress("MuonCand_charge",var_charge1);
  t2->SetBranchAddress("MuonCand_charge",var_charge2);
  t3->SetBranchAddress("MuonCand_charge",var_charge3);
  t4->SetBranchAddress("MuonCand_charge",var_charge4);
  t5->SetBranchAddress("MuonCand_charge",var_charge5);
  t6->SetBranchAddress("MuonCand_charge",var_charge6);


  t1->SetBranchAddress("MuonCand_pt",var_pt1);
  t2->SetBranchAddress("MuonCand_pt",var_pt2);
  t1->SetBranchAddress("MuonCand_pz",var_pz1);
  t2->SetBranchAddress("MuonCand_pz",var_pz2);
  t1->SetBranchAddress("MuonCand_phi",var_phi1);
  t2->SetBranchAddress("MuonCand_phi",var_phi2);
  t1->SetBranchAddress("MuonCand_eta",var_eta1);
  t2->SetBranchAddress("MuonCand_eta",var_eta2);
  t1->SetBranchAddress("MuonCand_p",var_p1);
  t1->SetBranchAddress("MuonCand_px",var_px1);
  t1->SetBranchAddress("MuonCand_py",var_py1);
  t1->SetBranchAddress("MuonCand_pz",var_pz1);
  t2->SetBranchAddress("MuonCand_p",var_p2);
  t2->SetBranchAddress("MuonCand_px",var_px2);
  t2->SetBranchAddress("MuonCand_py",var_py2);
  t2->SetBranchAddress("MuonCand_pz",var_pz2);
  t0->SetBranchAddress("MuonCand_pt",var_pt0);
  t0->SetBranchAddress("MuonCand_pz",var_pz0);
  t0->SetBranchAddress("MuonCand_phi",var_phi0);
  t0->SetBranchAddress("MuonCand_eta",var_eta0);
  t0->SetBranchAddress("MuonCand_p",var_p0);
  t0->SetBranchAddress("MuonCand_px",var_px0);
  t0->SetBranchAddress("MuonCand_py",var_py0);
  t0->SetBranchAddress("MuonCand_pz",var_pz0);
  t3->SetBranchAddress("MuonCand_pt",var_pt3);
  t3->SetBranchAddress("MuonCand_pz",var_pz3);
  t3->SetBranchAddress("MuonCand_phi",var_phi3);
  t3->SetBranchAddress("MuonCand_eta",var_eta3);
  t3->SetBranchAddress("MuonCand_p",var_p3);
  t3->SetBranchAddress("MuonCand_px",var_px3);
  t3->SetBranchAddress("MuonCand_py",var_py3);
  t3->SetBranchAddress("MuonCand_pz",var_pz3);
  t4->SetBranchAddress("MuonCand_pt",var_pt4);
  t4->SetBranchAddress("MuonCand_pz",var_pz4);
  t4->SetBranchAddress("MuonCand_phi",var_phi4);
  t4->SetBranchAddress("MuonCand_eta",var_eta4);
  t4->SetBranchAddress("MuonCand_p",var_p4);
  t4->SetBranchAddress("MuonCand_px",var_px4);
  t4->SetBranchAddress("MuonCand_py",var_py4);
  t4->SetBranchAddress("MuonCand_pz",var_pz4);
  t5->SetBranchAddress("MuonCand_pt",var_pt5);
  t5->SetBranchAddress("MuonCand_pz",var_pz5);
  t5->SetBranchAddress("MuonCand_phi",var_phi5);
  t5->SetBranchAddress("MuonCand_eta",var_eta5);
  t5->SetBranchAddress("MuonCand_p",var_p5);
  t5->SetBranchAddress("MuonCand_px",var_px5);
  t5->SetBranchAddress("MuonCand_py",var_py5);
  t5->SetBranchAddress("MuonCand_pz",var_pz5);
  t6->SetBranchAddress("MuonCand_pt",var_pt6);
  t6->SetBranchAddress("MuonCand_pz",var_pz6);
  t6->SetBranchAddress("MuonCand_phi",var_phi6);
  t6->SetBranchAddress("MuonCand_eta",var_eta6);
  t6->SetBranchAddress("MuonCand_p",var_p6);
  t6->SetBranchAddress("MuonCand_px",var_px6);
  t6->SetBranchAddress("MuonCand_py",var_py6);
  t6->SetBranchAddress("MuonCand_pz",var_pz6);

  t0->SetBranchAddress("MuonPairCand",var_Pair0);
  t1->SetBranchAddress("MuonPairCand",var_Pair1);
  t2->SetBranchAddress("MuonPairCand",var_Pair2);
  t3->SetBranchAddress("MuonPairCand",var_Pair3);
  t4->SetBranchAddress("MuonPairCand",var_Pair4);
  t5->SetBranchAddress("MuonPairCand",var_Pair5);
  t6->SetBranchAddress("MuonPairCand",var_Pair6);

  t1->SetBranchAddress("MuonCand_efficiency",var_eff1);
  t2->SetBranchAddress("MuonCand_efficiency",var_eff2);
  t3->SetBranchAddress("MuonCand_efficiency",var_eff3);
  t4->SetBranchAddress("MuonCand_efficiency",var_eff4);
  t5->SetBranchAddress("MuonCand_efficiency",var_eff5);
  t6->SetBranchAddress("MuonCand_efficiency",var_eff6);

  t0->SetBranchAddress("nZDChitCand",var_nZDC0);
  t1->SetBranchAddress("nZDChitCand",var_nZDC1);
  t2->SetBranchAddress("nZDChitCand",var_nZDC2);
  t3->SetBranchAddress("nZDChitCand",var_nZDC3);
  t4->SetBranchAddress("nZDChitCand",var_nZDC4);
  t5->SetBranchAddress("nZDChitCand",var_nZDC5);
  t6->SetBranchAddress("nZDChitCand",var_nZDC6);
  t0->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus0);
  t1->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus1);
  t2->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus2);
  t3->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus3);
  t4->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus4);
  t5->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus5);
  t6->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus6);
  t0->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus0);
  t1->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus1);
  t2->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus2);
  t3->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus3);
  t4->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus4);
  t5->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus5);
  t6->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus6);
  t0->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus0);
  t1->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus1);
  t2->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus2);
  t3->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus3);
  t4->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus4);
  t5->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus5);
  t6->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus6);
  t0->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus0);
  t1->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus1);
  t2->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus2);
  t3->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus3);
  t4->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus4);
  t5->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus5);
  t6->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus6);
  t0->SetBranchAddress("ZDChit_time",var_zdcTime0);
  t1->SetBranchAddress("ZDChit_time",var_zdcTime1);
  t2->SetBranchAddress("ZDChit_time",var_zdcTime2);
  t3->SetBranchAddress("ZDChit_time",var_zdcTime3);
  t4->SetBranchAddress("ZDChit_time",var_zdcTime4);
  t5->SetBranchAddress("ZDChit_time",var_zdcTime5);
  t6->SetBranchAddress("ZDChit_time",var_zdcTime6);
  t0->SetBranchAddress("ZDChit_energy",var_zdcE0);
  t1->SetBranchAddress("ZDChit_energy",var_zdcE1);
  t2->SetBranchAddress("ZDChit_energy",var_zdcE2);
  t3->SetBranchAddress("ZDChit_energy",var_zdcE3);
  t4->SetBranchAddress("ZDChit_energy",var_zdcE4);
  t5->SetBranchAddress("ZDChit_energy",var_zdcE5);
  t6->SetBranchAddress("ZDChit_energy",var_zdcE6);
  t0->SetBranchAddress("ZDChit_section",var_zdcsection0);
  t1->SetBranchAddress("ZDChit_section",var_zdcsection1);
  t2->SetBranchAddress("ZDChit_section",var_zdcsection2);
  t3->SetBranchAddress("ZDChit_section",var_zdcsection3);
  t4->SetBranchAddress("ZDChit_section",var_zdcsection4);
  t5->SetBranchAddress("ZDChit_section",var_zdcsection5);
  t6->SetBranchAddress("ZDChit_section",var_zdcsection6);

  t0->SetBranchAddress("nCastorTowerCand",var_nCastor0);
  t1->SetBranchAddress("nCastorTowerCand",var_nCastor1);
  t2->SetBranchAddress("nCastorTowerCand",var_nCastor2);
  t3->SetBranchAddress("nCastorTowerCand",var_nCastor3);
  t4->SetBranchAddress("nCastorTowerCand",var_nCastor4);
  t5->SetBranchAddress("nCastorTowerCand",var_nCastor5);
  t6->SetBranchAddress("nCastorTowerCand",var_nCastor6);
  t0->SetBranchAddress("CastorTower_e",var_CastorE0);
  t1->SetBranchAddress("CastorTower_e",var_CastorE1);
  t2->SetBranchAddress("CastorTower_e",var_CastorE2);
  t3->SetBranchAddress("CastorTower_e",var_CastorE3);
  t4->SetBranchAddress("CastorTower_e",var_CastorE4);
  t5->SetBranchAddress("CastorTower_e",var_CastorE5);
  t6->SetBranchAddress("CastorTower_e",var_CastorE6);
  t0->SetBranchAddress("CastorTower_eta",var_CastorEta0);
  t1->SetBranchAddress("CastorTower_eta",var_CastorEta1);
  t2->SetBranchAddress("CastorTower_eta",var_CastorEta2);
  t3->SetBranchAddress("CastorTower_eta",var_CastorEta3);
  t4->SetBranchAddress("CastorTower_eta",var_CastorEta4);
  t5->SetBranchAddress("CastorTower_eta",var_CastorEta5);
  t6->SetBranchAddress("CastorTower_eta",var_CastorEta6);
  t0->SetBranchAddress("CastorTower_phi",var_CastorPhi0);
  t1->SetBranchAddress("CastorTower_phi",var_CastorPhi1);
  t2->SetBranchAddress("CastorTower_phi",var_CastorPhi2);
  t3->SetBranchAddress("CastorTower_phi",var_CastorPhi3);
  t4->SetBranchAddress("CastorTower_phi",var_CastorPhi4);
  t5->SetBranchAddress("CastorTower_phi",var_CastorPhi5);
  t6->SetBranchAddress("CastorTower_phi",var_CastorPhi6);
  t0->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit0);
  t1->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit1);
  t2->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit2);
  t3->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit3);
  t4->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit4);
  t5->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit5);
  t6->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit6);

  t0->SetBranchAddress("L1TechnicalTriggers",techBit0);
  t1->SetBranchAddress("L1TechnicalTriggers",techBit1);
  t2->SetBranchAddress("L1TechnicalTriggers",techBit2);
  t3->SetBranchAddress("L1TechnicalTriggers",techBit3);
  t4->SetBranchAddress("L1TechnicalTriggers",techBit4);
  t5->SetBranchAddress("L1TechnicalTriggers",techBit5);

  int filter0Gen(0);
  int filter0Track(0);
  int filter0Events(0);
  for(Int_t i = 0;i < NUM0;i++){
      	t0->GetEntry(i);
	int pair1 = var_Pair0[0]; int pair2 = var_Pair0[1];
        int muID1 = var_idA0[pair1];
        int muID2 = var_idA0[pair2];
	int muAng1 = var_idB0[pair1];
        int muAng2 = var_idB0[pair2];
        int hlt_pass = hlt_d0[0];
	int nPrimVtx = var_nvtx0[0];
	int nTrack=var_nTrack0[0];
	int nTrackQual=var_nTrackQual0[0]; 
	int nCalo=var_ncalo0[0];
	int label_vertex(99);
	double openangle(-1.);
//	cout<<"--------------------"<<var_event0[0]<<"-----------------------"<<endl;

        if(PassesTrigger(hlt_pass,var_run0[0]) == false) 
                continue;  

        if(PassesKinematicCuts(var_mass0[0],var_pt0[pair1],var_pt0[pair2],var_eta0[pair1],var_eta0[pair2]) == false)  
                continue; 

	if(nPrimVtx>=1){
	  double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack0,var_vtxZ0,var_vtxmumu0);          

          for(Int_t j=0; j<nPrimVtx; j++){
		if(PassesVertexSelection(var_vtxTrack0[j],var_vertexChi2_0[j],var_vertexNdf0[j],distance_vertex_z,var_vtxZ0[j],var_vtxmumu0[j])
                   && (techBit0[0][0]==1))   
			{label_vertex=j;}
          }
        }
	if(label_vertex!=99
	   && sqrt(pow(var_vtxX0[label_vertex],2)+pow(var_vtxY0[label_vertex],2))< 0.2
           && PassesMuonID(var_tracker0[pair1], muAng1, var_global0[pair1], var_tracker0[pair2], muAng2, var_global0[pair2], var_nhitsTrack0[pair1], var_nhitsTrack0[pair2])
 	   && PassesDptCut(var_dpt0[0]) 
	   && PassesDphiCut(var_dphi0[0]/pi)
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
		if(FailsTrackDistanceVeto(var_TrackD0[j], var_TrackQuality0[j])){
                  filter0Track++;
  		  nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0), nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR0[k]>dRcone){
                        //remove known noisy towers
                                if(var_caloId0[k]==4 && var_run0[0]>=139779 && var_run0[0]<=140159
                                                     && (double)(var_caloEta0[k])>1.1749 && (double)(var_caloEta0[k])<1.175
                                                     && (double)(var_caloPhi0[k])>2.661  && (double)(var_caloPhi0[k])<2.662) {continue;}

                                if(var_caloId0[k]==2 && var_caloEn0[k]>EBThresh) nEB++;
                                if(var_caloId0[k]==3 && var_caloEn0[k]>EEThresh) nEE++;
                                if(var_caloId0[k]==4 && var_caloEn0[k]>HBThresh) nHB++;
                                if(var_caloId0[k]==5 && var_caloEn0[k]>HEThresh) nHE++;
                                if(var_caloId0[k]==1 && var_caloEta0[k]>0 && var_caloEn0[k]>HFPlusThresh) nHFp++;
                                if(var_caloId0[k]==1 && var_caloEta0[k]<0 && var_caloEn0[k]>HFMinusThresh) nHFm++;

                                if(var_caloId0[k]==6){
                                        if(var_caloHadE0[k]>HBThresh) nHB++;
                                        if(var_caloEmE0[k]>EBThresh) nEB++;}
                                if(var_caloId0[k]==7){
                                        if(var_caloHadE0[k]>HEThresh) nHE++;
                                        if(var_caloEmE0[k]>EBThresh) nEB++;}
                                if(var_caloId0[k]==8){
                                        if(var_caloHadE0[k]>HEThresh) nHE++;
                                        if(var_caloEmE0[k]>EEThresh) nEE++;}
	       }
	    }

	if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus0[0],var_zdcEmPlus0[0],var_zdcHadMinus0[0],var_zdcHadPlus0[0]))
	{
		nTower0->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi0); 
	}

	if(nTrackExclu<1 
		&&  PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm)
		&&  PassesZDCVeto(var_zdcEmMinus0[0],var_zdcEmPlus0[0],var_zdcHadMinus0[0],var_zdcHadPlus0[0]))
	{ 
          filter0Events++;

	  hEB0->Fill(nEB,fac_lumi0);
          hEE0->Fill(nEE,fac_lumi0);
          hHB0->Fill(nHB,fac_lumi0);
          hHE0->Fill(nHE,fac_lumi0);
          hHFp0->Fill(nHFp,fac_lumi0);
          hHFm0->Fill(nHFm,fac_lumi0);
	  nTrack0->Fill(nTrackExclu,fac_lumi0);

	  MuMuMass0->Fill(var_mass0[0],fac_lumi0);
          MuMuMassUps0->Fill(var_mass0[0],fac_lumi0);
	  MuMuMassJpsi0->Fill(var_mass0[0],fac_lumi0);
          MuMudpt0->Fill(var_dpt0[0],fac_lumi0);
          MuMudphi0->Fill(var_dphi0[0]/pi,fac_lumi0);
	  double symdphi0 = 1 - fabs(var_phi0[pair1]-var_phi0[pair2])/pi;   
	  MuMuSymdphi0->Fill(symdphi0,fac_lumi0);
	//          MuMudeta0->Fill(fabs(var_eta0[pair1]+var_eta0[pair2]),fac_lumi0);
	  MuMuvtxXY0->Fill(sqrt(var_MuMuvtxX0[0]*var_MuMuvtxX0[0]+var_MuMuvtxY0[0]*var_MuMuvtxY0[0]),fac_lumi0);

	  TLorentzVector mu1, mu2, dimuon;
	  mu1.SetPtEtaPhiM(var_pt0[pair1],var_eta0[pair1],var_phi0[pair1],0.1057);	
          mu2.SetPtEtaPhiM(var_pt0[pair2],var_eta0[pair2],var_phi0[pair2],0.1057);   
	  dimuon = mu1 + mu2;
	  Tdist0->Fill(dimuon.Pt()*dimuon.Pt()); 
	  MuMudeta0->Fill(mu1.Angle(mu2.Vect()),fac_lumi0);
	  openangle = mu1.Angle(mu2.Vect());
	  cout<<"candidate  Run "<<var_run0[0]<<"  LS "<<var_ls0[0]<<"\tEvt "<<var_event0[0]<<"\t mass="<<var_mass0[0]<<" GeV"<<" Opening angle="<<openangle<<endl; 

	  ZDCemminus0->Fill(var_zdcEmMinus0[0],fac_lumi0); ZDCemplus0->Fill(var_zdcEmPlus0[0],fac_lumi0);
	  ZDChadminus0->Fill(var_zdcHadMinus0[0],fac_lumi0); ZDChadplus0->Fill(var_zdcHadPlus0[0],fac_lumi0);
	  for(Int_t l=0; l<var_nZDC0[0]; l++){
	     if(var_zdcsection0[l]==1 && var_zdcE0[l]>ZDCemThresh){ 
                ZDCtime0->Fill(var_zdcTime0[l],fac_lumi0); ZDCenergyEM0->Fill(var_zdcE0[l],fac_lumi0);
	     }
	     if(var_zdcsection0[l]==2 && var_zdcE0[l]>ZDChadThresh){
		ZDCtime0->Fill(var_zdcTime0[l],fac_lumi0); ZDCenergyHAD0->Fill(var_zdcE0[l],fac_lumi0);
	     }
	  }

	  CastorSumE0->Fill(var_CastorRecHit0[0],fac_lumi0);

	  double vertexT=sqrt(pow(var_vtxX0[label_vertex],2)+pow(var_vtxY0[label_vertex],2));
	  VtxT0->Fill(vertexT,fac_lumi0);

          double eta_pair=0.5*TMath::Log((double)((var_p0[pair1]+var_p0[pair2]+var_pz0[pair1]+var_pz0[pair2])/(var_p0[pair1]+var_p0[pair2]-var_pz0[pair1]-var_pz0[pair2])));
          double rap_pair=0.5*TMath::Log((double)((sqrt(var_p0[pair1]*var_p0[pair1]+0.1057*0.1057)+sqrt(var_p0[pair2]*var_p0[pair2]+0.1057*0.1057)+var_pz0[pair1]+var_pz0[pair2])/(sqrt(var_p0[pair1]*var_p0[pair1]+0.1057*0.1057)+sqrt(var_p0[pair2]*var_p0[pair2]+0.1057*0.1057)-var_pz0[pair1]-var_pz0[pair2])));
	  double pt_pair =sqrt((var_px0[pair1]+var_px0[pair2])*(var_px0[pair1]+var_px0[pair2])+(var_py0[pair1]+var_py0[pair2])*(var_py0[pair1]+var_py0[pair2]));
          etaPair0->Fill(rap_pair,fac_lumi0);
	  pTPair0->Fill(pt_pair,fac_lumi0);

	  if(var_charge0[pair1]>0) {phiSingleP0->Fill(var_phi0[pair1]/pi,fac_lumi0);etaSingleP0->Fill(var_eta0[pair1],fac_lumi0);pTSingleP0->Fill(var_pt0[pair1],fac_lumi0);}   
	  else {phiSingleM0->Fill(var_phi0[pair1]/pi,fac_lumi0);etaSingleM0->Fill(var_eta0[pair1],fac_lumi0);pTSingleM0->Fill(var_pt0[pair1],fac_lumi0);}

	  if(var_charge0[pair2]>0) {phiSingleP0->Fill(var_phi0[pair2]/pi,fac_lumi0);etaSingleP0->Fill(var_eta0[pair2],fac_lumi0);pTSingleP0->Fill(var_pt0[pair2],fac_lumi0);}   
	  else {phiSingleM0->Fill(var_phi0[pair2]/pi,fac_lumi0);etaSingleM0->Fill(var_eta0[pair2],fac_lumi0);pTSingleM0->Fill(var_pt0[pair2],fac_lumi0);}

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"Data :"<<endl;
cout<<"  # Dimuon events = "<<filter0Events<<endl;

  int filter1Gen(0);
 int filter1Track(0);
  double filter1Events(0.);
  for(Int_t i = 0;i < NUM1;i++){
      	t1->GetEntry(i);
	int pair1 = var_Pair1[0]; int pair2 = var_Pair1[1];
        int muID1 = var_idA1[pair1];
        int muID2 = var_idA1[pair2];
	int muAng1 = var_idB1[pair1];
        int muAng2 = var_idB1[pair2];
	double effcorrection1 = (var_eff1[pair1]*var_eff1[pair2]*doublemuopenfractionallumi+doublemuopentightfractionallumi);
        int hlt_pass = hlt_d1[0];
	int nPrimVtx = var_nvtx1[0];
	int nTrack=var_nTrack1[0];
	int nTrackQual=var_nTrackQual1[0]; 
	int nCalo=var_ncalo1[0];
	int label_vertex(99);
//	cout<<"--------------------"<<i<<"-----------------------"<<endl;

	if(PassesTrigger(hlt_pass,1) == false)
		continue; 

        if(PassesKinematicCuts(var_mass1[0],var_pt1[pair1],var_pt1[pair2],var_eta1[pair1],var_eta1[pair2]) == false)  
		continue; 

	if(nPrimVtx>=1){
          double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack1,var_vtxZ1,var_vtxmumu1);           
	  for(Int_t j=0; j<nPrimVtx; j++){
                if(PassesVertexSelection(var_vtxTrack1[j],var_vertexChi2_1[j],var_vertexNdf1[j],distance_vertex_z,var_vtxZ1[j],var_vtxmumu1[j]) 
		   && !(techBit1[0][0]==1))   
			{label_vertex=j;}
	  }
	}

	if(label_vertex!=99
	   && PassesMuonID(var_tracker1[pair1], muAng1, var_global1[pair1], var_tracker1[pair2], muAng2, var_global1[pair2], var_nhitsTrack1[pair1], var_nhitsTrack1[pair2]) 
           && PassesDptCut(var_dpt1[0])  
           && PassesDphiCut(var_dphi1[0]/pi) 
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
		if(FailsTrackDistanceVeto(var_TrackD1[j], var_TrackQuality1[j])){
                  filter1Track++;
		  nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR1[k]>dRcone){
                                if(var_caloId1[k]==2 && var_caloEn1[k]>EBThresh) nEB++;
                                if(var_caloId1[k]==3 && var_caloEn1[k]>EEThresh) nEE++;
                                if(var_caloId1[k]==4 && var_caloEn1[k]>HBThresh) nHB++;
                                if(var_caloId1[k]==5 && var_caloEn1[k]>HEThresh) nHE++;
                                if(var_caloId1[k]==1 && var_caloEta1[k]>0 && var_caloEn1[k]>HFPlusThresh) nHFp++;
                                if(var_caloId1[k]==1 && var_caloEta1[k]<0 && var_caloEn1[k]>HFMinusThresh) nHFm++;

                                if(var_caloId1[k]==6){
                                        if(var_caloHadE1[k]>HBThresh) nHB++;
                                        if(var_caloEmE1[k]>EBThresh) nEB++;}
                                if(var_caloId1[k]==7){
                                        if(var_caloHadE1[k]>HEThresh) nHE++;
                                        if(var_caloEmE1[k]>EBThresh) nEB++;}
                                if(var_caloId1[k]==8){
                                        if(var_caloHadE1[k]>HEThresh) nHE++;
                                        if(var_caloEmE1[k]>EEThresh) nEE++;}
	       }
	    }

	if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus1[0],var_zdcEmPlus1[0],var_zdcHadMinus1[0],var_zdcHadPlus1[0]))
	{
		nTower1->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi1*effcorrection1); 
	}

	if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm) 
		&& PassesZDCVeto(var_zdcEmMinus1[0],var_zdcEmPlus1[0],var_zdcHadMinus1[0],var_zdcHadPlus1[0]))
	){ 
          filter1Events+=fac_lumi1*effcorrection1;

	  hEB1->Fill(nEB,fac_lumi1*effcorrection1);
          hEE1->Fill(nEE,fac_lumi1*effcorrection1);
          hHB1->Fill(nHB,fac_lumi1*effcorrection1);
          hHE1->Fill(nHE,fac_lumi1*effcorrection1);
          hHFp1->Fill(nHFp,fac_lumi1*effcorrection1);
          hHFm1->Fill(nHFm,fac_lumi1*effcorrection1);
	  nTrack1->Fill(nTrackExclu,fac_lumi1*effcorrection1);

	  MuMuMass1->Fill(var_mass1[0],fac_lumi1*effcorrection1);
          MuMuMassUps1->Fill(var_mass1[0],fac_lumi1*effcorrection1);
          MuMuMassJpsi1->Fill(var_mass1[0],fac_lumi1*effcorrection1);
          MuMudpt1->Fill(var_dpt1[0],fac_lumi1*effcorrection1);
          MuMudphi1->Fill(var_dphi1[0]/pi,fac_lumi1*effcorrection1);
          double symdphi1 = 1 - fabs(var_phi1[pair1]-var_phi1[pair2])/pi;    
          MuMuSymdphi1->Fill(symdphi1,fac_lumi1*effcorrection1); 

	//          MuMudeta1->Fill(fabs(var_eta1[pair1]+var_eta1[pair2]),fac_lumi1*effcorrection1);
          MuMuvtxXY1->Fill(sqrt(var_MuMuvtxX1[0]*var_MuMuvtxX1[0]+var_MuMuvtxY1[0]*var_MuMuvtxY1[0]),fac_lumi1*effcorrection1);

          TLorentzVector mu11, mu12, dimuon1;   
          mu11.SetPtEtaPhiM(var_pt1[pair1],var_eta1[pair1],var_phi1[pair1],0.1057);     
          mu12.SetPtEtaPhiM(var_pt1[pair2],var_eta1[pair2],var_phi1[pair2],0.1057);      
          dimuon1 = mu11 + mu12;   
          Tdist1->Fill(dimuon1.Pt()*dimuon1.Pt(),fac_lumi1*effcorrection1);    
          MuMudeta1->Fill(mu11.Angle(mu12.Vect()),fac_lumi1*effcorrection1); 

	  ZDCemminus1->Fill(var_zdcEmMinus1[0],fac_lumi1*effcorrection1);
	  ZDCemplus1->Fill(var_zdcEmPlus1[0],fac_lumi1*effcorrection1);
	  ZDChadminus1->Fill(var_zdcHadMinus1[0],fac_lumi1*effcorrection1);
	  ZDChadplus1->Fill(var_zdcHadPlus1[0],fac_lumi1*effcorrection1); 

	  for(Int_t l=0; l<var_nZDC1[0]; l++){
	     if(var_zdcsection1[l]==1 && var_zdcE1[l]>ZDCemThresh){ 
                ZDCtime1->Fill(var_zdcTime1[l],fac_lumi1*effcorrection1); 
                ZDCenergyEM1->Fill(var_zdcE1[l],fac_lumi1*effcorrection1);
	     }
	     if(var_zdcsection1[l]==2 && var_zdcE1[l]>ZDChadThresh){
		ZDCtime1->Fill(var_zdcTime1[l],fac_lumi1*effcorrection1);
                ZDCenergyHAD1->Fill(var_zdcE1[l],fac_lumi1*effcorrection1);
	     }
	  }

	  CastorSumE1->Fill(var_CastorRecHit1[0],fac_lumi1*effcorrection1);

          double vertexT=sqrt(pow(var_vtxX1[label_vertex],2)+pow(var_vtxY1[label_vertex],2))-0.3633;
          VtxT1->Fill(vertexT,fac_lumi1*effcorrection1);

          double eta_pair=0.5*TMath::Log((double)((var_p1[pair1]+var_p1[pair2]+var_pz1[pair1]+var_pz1[pair2])/(var_p1[pair1]+var_p1[pair2]-var_pz1[pair1]-var_pz1[pair2])));
          double rap_pair=0.5*TMath::Log((double)((sqrt(var_p1[pair1]*var_p1[pair1]+0.1057*0.1057)+sqrt(var_p1[pair2]*var_p1[pair2]+0.1057*0.1057)+var_pz1[pair1]+var_pz1[pair2])/(sqrt(var_p1[pair1]*var_p1[pair1]+0.1057*0.1057)+sqrt(var_p1[pair2]*var_p1[pair2]+0.1057*0.1057)-var_pz1[pair1]-var_pz1[pair2])));
          double pt_pair =sqrt((var_px1[pair1]+var_px1[pair2])*(var_px1[pair1]+var_px1[pair2])+(var_py1[pair1]+var_py1[pair2])*(var_py1[pair1]+var_py1[pair2]));
          etaPair1->Fill(rap_pair,fac_lumi1*effcorrection1);
          pTPair1->Fill(pt_pair,fac_lumi1*effcorrection1);

          if(var_charge1[pair1]>0) {phiSingleP1->Fill(var_phi1[pair1]/pi,fac_lumi1*effcorrection1);etaSingleP1->Fill(var_eta1[pair1],fac_lumi1*effcorrection1);pTSingleP1->Fill(var_pt1[pair1],fac_lumi1*effcorrection1);}   
	  else {phiSingleM1->Fill(var_phi1[pair1]/pi,fac_lumi1*effcorrection1);etaSingleM1->Fill(var_eta1[pair1],fac_lumi1*effcorrection1);pTSingleM1->Fill(var_pt1[pair1],fac_lumi1*effcorrection1);}
          if(var_charge1[pair2]>0) {phiSingleP1->Fill(var_phi1[pair2]/pi,fac_lumi1*effcorrection1);etaSingleP1->Fill(var_eta1[pair2],fac_lumi1*effcorrection1);pTSingleP1->Fill(var_pt1[pair2],fac_lumi1*effcorrection1);}   
	  else {phiSingleM1->Fill(var_phi1[pair2]/pi,fac_lumi1*effcorrection1);etaSingleM1->Fill(var_eta1[pair2],fac_lumi1*effcorrection1);pTSingleM1->Fill(var_pt1[pair2],fac_lumi1*effcorrection1);}

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"ElEl :"<<endl;
cout<<"  # Dimuon events = "<<filter1Events<<endl;

  int filter2Gen(0);
  int filter2Track(0);
  double filter2Events(0);
  for(Int_t i = 0;i < NUM2;i++){
      	t2->GetEntry(i);
	int pair2 = var_Pair2[0]; int pair2 = var_Pair2[1];
        int muID1 = var_idA2[pair1];
        int muID2 = var_idA2[pair2];
	int muAng1 = var_idB2[pair1];
        int muAng2 = var_idB2[pair2];
	double effcorrection2 = (var_eff2[pair1]*var_eff2[pair2]*doublemuopenfractionallumi+doublemuopentightfractionallumi);
        int hlt_pass = hlt_d2[0];
	int nPrimVtx = var_nvtx2[0];
	int nTrack=var_nTrack2[0];
	int nTrackQual=var_nTrackQual2[0]; 
	int nCalo=var_ncalo2[0];
	int label_vertex(99);
//	cout<<"--------------------"<<var_event2[0]<<"-----------------------"<<endl;

       if(PassesTrigger(hlt_pass,1) == false) 
                continue;  

        if(PassesKinematicCuts(var_mass2[0],var_pt2[pair1],var_pt2[pair2],var_eta2[pair1],var_eta2[pair2]) == false)   
                continue;  

        if(nPrimVtx>=1){
          double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack2,var_vtxZ2,var_vtxmumu2);           
          for(Int_t j=0; j<nPrimVtx; j++){
                if(PassesVertexSelection(var_vtxTrack2[j],var_vertexChi2_2[j],var_vertexNdf2[j],distance_vertex_z,var_vtxZ2[j],var_vtxmumu2[j]) 
                   && !(techBit2[0][0]==1))   
			{label_vertex=j;}
          }
        }
	if(label_vertex!=99
           && PassesMuonID(var_tracker2[pair1], muAng1, var_global2[pair1], var_tracker2[pair2], muAng2, var_global2[pair2], var_nhitsTrack2[pair1], var_nhitsTrack2[pair2])  
           && PassesDptCut(var_dpt2[0])   
           && PassesDphiCut(var_dphi2[0]/pi)  
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(FailsTrackDistanceVeto(var_TrackD2[j], var_TrackQuality2[j])){ 
                  filter2Track++;
		  nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR2[k]>dRcone){
                                if(var_caloId2[k]==2 && var_caloEn2[k]>EBThresh) nEB++;
                                if(var_caloId2[k]==3 && var_caloEn2[k]>EEThresh) nEE++;
                                if(var_caloId2[k]==4 && var_caloEn2[k]>HBThresh) nHB++;
                                if(var_caloId2[k]==5 && var_caloEn2[k]>HEThresh) nHE++;
                                if(var_caloId2[k]==1 && var_caloEta2[k]>0 && var_caloEn2[k]>HFPlusThresh) nHFp++;
                                if(var_caloId2[k]==1 && var_caloEta2[k]<0 && var_caloEn2[k]>HFMinusThresh) nHFm++;

                                if(var_caloId2[k]==6){
                                        if(var_caloHadE2[k]>HBThresh) nHB++;
                                        if(var_caloEmE2[k]>EBThresh) nEB++;}
                                if(var_caloId2[k]==7){
                                        if(var_caloHadE2[k]>HEThresh) nHE++;
                                        if(var_caloEmE2[k]>EBThresh) nEB++;}
                                if(var_caloId2[k]==8){
                                        if(var_caloHadE2[k]>HEThresh) nHE++;
                                        if(var_caloEmE2[k]>EEThresh) nEE++;}
	       }
	    }

	if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus2[0],var_zdcEmPlus2[0],var_zdcHadMinus2[0],var_zdcHadPlus2[0]))
	{
		nTower2->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi2*effcorrection2);
	}

	if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm) 
		&& PassesZDCVeto(var_zdcEmMinus2[0],var_zdcEmPlus2[0],var_zdcHadMinus2[0],var_zdcHadPlus2[0]))
	){ 
          filter2Events+=fac_lumi2*effcorrection2;

	  hEB2->Fill(nEB,fac_lumi2*effcorrection2);
          hEE2->Fill(nEE,fac_lumi2*effcorrection2);
          hHB2->Fill(nHB,fac_lumi2*effcorrection2);
          hHE2->Fill(nHE,fac_lumi2*effcorrection2);
          hHFp2->Fill(nHFp,fac_lumi2*effcorrection2);
          hHFm2->Fill(nHFm,fac_lumi2*effcorrection2); 
//	  nTower2->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi2*effcorrection2);
	  nTrack2->Fill(nTrackExclu,fac_lumi2*effcorrection2);

	  MuMuMass2->Fill(var_mass2[0],fac_lumi2*effcorrection2); 
          MuMuMassUps2->Fill(var_mass2[0],fac_lumi2*effcorrection2);
          MuMuMassJpsi2->Fill(var_mass2[0],fac_lumi2*effcorrection2);
          MuMudpt2->Fill(var_dpt2[0],fac_lumi2*effcorrection2); 
          MuMudphi2->Fill(var_dphi2[0]/pi,fac_lumi2*effcorrection2); 
          double symdphi2 = 1 - fabs(var_phi2[pair1]-var_phi2[pair2])/pi;    
          MuMuSymdphi2->Fill(symdphi2,fac_lumi2*effcorrection2); 
	//          MuMudeta2->Fill(fabs(var_eta2[pair1]+var_eta2[pair2]),fac_lumi2*effcorrection2);
	  MuMuvtxXY2->Fill(sqrt(var_MuMuvtxX2[0]*var_MuMuvtxX2[0]+var_MuMuvtxY2[0]*var_MuMuvtxY2[0]),fac_lumi2*effcorrection2);

          TLorentzVector mu21, mu22, dimuon2;    
          mu21.SetPtEtaPhiM(var_pt2[pair1],var_eta2[pair1],var_phi2[pair1],0.1057);      
          mu22.SetPtEtaPhiM(var_pt2[pair2],var_eta2[pair2],var_phi2[pair2],0.1057);       
          dimuon2 = mu21 + mu22;    
          Tdist2->Fill(dimuon2.Pt()*dimuon2.Pt(),fac_lumi2*effcorrection2);     
          MuMudeta2->Fill(mu21.Angle(mu22.Vect()),fac_lumi2*effcorrection2); 

	  ZDCemminus2->Fill(var_zdcEmMinus2[0],fac_lumi2*effcorrection2); 
	  ZDCemplus2->Fill(var_zdcEmPlus2[0],fac_lumi2*effcorrection2);
	  ZDChadminus2->Fill(var_zdcHadMinus2[0],fac_lumi2*effcorrection2); 
	  ZDChadplus2->Fill(var_zdcHadPlus2[0],fac_lumi2*effcorrection2); 

	  for(Int_t l=0; l<var_nZDC2[0]; l++){
	     if(var_zdcsection2[l]==1 && var_zdcE2[l]>ZDCemThresh){ 
                ZDCtime2->Fill(var_zdcTime2[l],fac_lumi2*effcorrection2);
                ZDCenergyEM2->Fill(var_zdcE2[l],fac_lumi2*effcorrection2);
	     }
	     if(var_zdcsection2[l]==2 && var_zdcE2[l]>ZDChadThresh){
		ZDCtime2->Fill(var_zdcTime2[l],fac_lumi2*effcorrection2);
                ZDCenergyHAD2->Fill(var_zdcE2[l],fac_lumi2*effcorrection2);
	     }
	  }

	  CastorSumE2->Fill(var_CastorRecHit2[0],fac_lumi2*effcorrection2);

          double vertexT=sqrt(pow(var_vtxX2[label_vertex],2)+pow(var_vtxY2[label_vertex],2))-0.3633;
          VtxT2->Fill(vertexT,fac_lumi2*effcorrection2);

          double eta_pair=0.5*TMath::Log((double)((var_p2[pair1]+var_p2[pair2]+var_pz2[pair1]+var_pz2[pair2])/(var_p2[pair1]+var_p2[pair2]-var_pz2[pair1]-var_pz2[pair2])));
          double rap_pair=0.5*TMath::Log((double)((sqrt(var_p2[pair1]*var_p2[pair1]+0.1057*0.1057)+sqrt(var_p2[pair2]*var_p2[pair2]+0.1057*0.1057)+var_pz2[pair1]+var_pz2[pair2])/(sqrt(var_p2[pair1]*var_p2[pair1]+0.1057*0.1057)+sqrt(var_p2[pair2]*var_p2[pair2]+0.1057*0.1057)-var_pz2[pair1]-var_pz2[pair2])));
          double pt_pair =sqrt((var_px2[pair1]+var_px2[pair2])*(var_px2[pair1]+var_px2[pair2])+(var_py2[pair1]+var_py2[pair2])*(var_py2[pair1]+var_py2[pair2]));
          etaPair2->Fill(rap_pair,fac_lumi2*effcorrection2);
          pTPair2->Fill(pt_pair,fac_lumi2*effcorrection2);

          if(var_charge2[pair1]>0) {phiSingleP2->Fill(var_phi2[pair1]/pi,fac_lumi2*effcorrection2);etaSingleP2->Fill(var_eta2[pair1],fac_lumi2*effcorrection2);pTSingleP2->Fill(var_pt2[pair1],fac_lumi2*effcorrection2);}
          else {phiSingleM2->Fill(var_phi2[pair1]/pi,fac_lumi2*effcorrection2);etaSingleM2->Fill(var_eta2[pair1],fac_lumi2*effcorrection2);pTSingleM2->Fill(var_pt2[pair1],fac_lumi2*effcorrection2);}
          if(var_charge2[pair2]>0) {phiSingleP2->Fill(var_phi2[pair2]/pi,fac_lumi2*effcorrection2);etaSingleP2->Fill(var_eta2[pair2],fac_lumi2*effcorrection2);pTSingleP2->Fill(var_pt2[pair2],fac_lumi2*effcorrection2);}
          else {phiSingleM2->Fill(var_phi2[pair2]/pi,fac_lumi2*effcorrection2);etaSingleM2->Fill(var_eta2[pair2],fac_lumi2*effcorrection2);pTSingleM2->Fill(var_pt2[pair2],fac_lumi2*effcorrection2);}

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"InelEl :"<<endl;
cout<<"  # Dimuon events = "<<filter2Events<<endl;

  int filter3Gen(0);
  int filter3Track(0);
  double filter3Events(0);
  for(Int_t i = 0;i < NUM3;i++){
      	t3->GetEntry(i);
	int pair1 = var_Pair3[0]; int pair2 = var_Pair3[1];
        int muID1 = var_idA3[pair1];
        int muID2 = var_idA3[pair2];
	int muAng1 = var_idB3[pair1];
        int muAng2 = var_idB3[pair2];
	double effcorrection3 = (var_eff3[pair1]*var_eff3[pair2]*doublemuopenfractionallumi+doublemuopentightfractionallumi);
        int hlt_pass = hlt_d3[0];
	int nPrimVtx = var_nvtx3[0];
	int nTrack=var_nTrack3[0];
	int nTrackQual=var_nTrackQual3[0]; 
	int nCalo=var_ncalo3[0];
	int label_vertex(99);
//	cout<<"--------------------"<<var_event1[0]<<"-----------------------"<<endl;

       if(PassesTrigger(hlt_pass,1) == false) 
                continue;  

        if(PassesKinematicCuts(var_mass3[0],var_pt3[pair1],var_pt3[pair2],var_eta3[pair1],var_eta3[pair2]) == false)   
                continue;  

        if(nPrimVtx>=1){
          double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack3,var_vtxZ3,var_vtxmumu3);           
          for(Int_t j=0; j<nPrimVtx; j++){
                if(PassesVertexSelection(var_vtxTrack3[j],var_vertexChi2_3[j],var_vertexNdf3[j],distance_vertex_z,var_vtxZ3[j],var_vtxmumu3[j]) 
                   && !(techBit3[0][0]==1))   
			{label_vertex=j;}
          }
        }

	if(label_vertex!=99
           && PassesMuonID(var_tracker3[pair1], muAng1, var_global3[pair1], var_tracker3[pair2], muAng2, var_global3[pair2], var_nhitsTrack3[pair1], var_nhitsTrack3[pair2])   
	   && PassesDptCut(var_dpt3[0])   
           && PassesDphiCut(var_dphi3[0]/pi)  
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(FailsTrackDistanceVeto(var_TrackD3[j], var_TrackQuality3[j])){  
                  filter3Track++;
		  nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR3[k]>dRcone){
                                if(var_caloId3[k]==2 && var_caloEn3[k]>EBThresh) nEB++;
                                if(var_caloId3[k]==3 && var_caloEn3[k]>EEThresh) nEE++;
                                if(var_caloId3[k]==4 && var_caloEn3[k]>HBThresh) nHB++;
                                if(var_caloId3[k]==5 && var_caloEn3[k]>HEThresh) nHE++;
                                if(var_caloId3[k]==1 && var_caloEta3[k]>0 && var_caloEn3[k]>HFPlusThresh) nHFp++;
                                if(var_caloId3[k]==1 && var_caloEta3[k]<0 && var_caloEn3[k]>HFMinusThresh) nHFm++;

                                if(var_caloId3[k]==6){
                                        if(var_caloHadE3[k]>HBThresh) nHB++;
                                        if(var_caloEmE3[k]>EBThresh) nEB++;}
                                if(var_caloId3[k]==7){
                                        if(var_caloHadE3[k]>HEThresh) nHE++;
                                        if(var_caloEmE3[k]>EBThresh) nEB++;}
                                if(var_caloId3[k]==8){
                                        if(var_caloHadE3[k]>HEThresh) nHE++;
                                        if(var_caloEmE3[k]>EEThresh) nEE++;}
	       }
	    }

	if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus3[0],var_zdcEmPlus3[0],var_zdcHadMinus3[0],var_zdcHadPlus3[0]))
	{
		nTower3->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi3*effcorrection3);
	}

	if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm) 
		&& PassesZDCVeto(var_zdcEmMinus3[0],var_zdcEmPlus3[0],var_zdcHadMinus3[0],var_zdcHadPlus3[0]))	
	){ 
          filter3Events+=fac_lumi3*effcorrection3;

	  hEB3->Fill(nEB,fac_lumi3*effcorrection3); 
          hEE3->Fill(nEE,fac_lumi3*effcorrection3);
          hHB3->Fill(nHB,fac_lumi3*effcorrection3); 
          hHE3->Fill(nHE,fac_lumi3*effcorrection3);
          hHFp3->Fill(nHFp,fac_lumi3*effcorrection3);
          hHFm3->Fill(nHFm,fac_lumi3*effcorrection3);
	  nTrack3->Fill(nTrackExclu,fac_lumi3*effcorrection3);

	  MuMuMass3->Fill(var_mass3[0],fac_lumi3*effcorrection3);
          MuMuMassUps3->Fill(var_mass3[0],fac_lumi3*effcorrection3);
          MuMuMassJpsi3->Fill(var_mass3[0],fac_lumi3*effcorrection3);
          MuMudpt3->Fill(var_dpt3[0],fac_lumi3*effcorrection3);
          MuMudphi3->Fill(var_dphi3[0]/pi,fac_lumi3*effcorrection3);
          double symdphi3 = 1 - fabs(var_phi3[pair1]-var_phi3[pair2])/pi;    
          MuMuSymdphi3->Fill(symdphi3,fac_lumi3*effcorrection3); 
	//          MuMudeta3->Fill(fabs(var_eta3[pair1]+var_eta3[pair2]),fac_lumi3*effcorrection3);
	  MuMuvtxXY3->Fill(sqrt(var_MuMuvtxX3[0]*var_MuMuvtxX3[0]+var_MuMuvtxY3[0]*var_MuMuvtxY3[0]),fac_lumi3*effcorrection3);

	  ZDCemminus3->Fill(var_zdcEmMinus3[0],fac_lumi3*effcorrection3);
	  ZDCemplus3->Fill(var_zdcEmPlus3[0],fac_lumi3*effcorrection3);
	  ZDChadminus3->Fill(var_zdcHadMinus3[0],fac_lumi3*effcorrection3);
	  ZDChadplus3->Fill(var_zdcHadPlus3[0],fac_lumi3*effcorrection3);

          TLorentzVector mu31, mu32, dimuon3; 
          mu31.SetPtEtaPhiM(var_pt3[pair1],var_eta3[pair1],var_phi3[pair1],0.1057);   
          mu32.SetPtEtaPhiM(var_pt3[pair2],var_eta3[pair2],var_phi3[pair2],0.1057);    
          dimuon3 = mu31 + mu32; 
          Tdist3->Fill(dimuon3.Pt()*dimuon3.Pt(),fac_lumi3*effcorrection3);  
          MuMudeta3->Fill(mu31.Angle(mu32.Vect()),fac_lumi3*effcorrection3);  

	  for(Int_t l=0; l<var_nZDC3[0]; l++){
	     if(var_zdcsection3[l]==1 && var_zdcE3[l]>ZDCemThresh){ 
                ZDCtime3->Fill(var_zdcTime3[l],fac_lumi3*effcorrection3); 
                ZDCenergyEM3->Fill(var_zdcE3[l],fac_lumi3*effcorrection3);
	     }
	     if(var_zdcsection3[l]==2 && var_zdcE3[l]>ZDChadThresh){
		ZDCtime3->Fill(var_zdcTime3[l],fac_lumi3*effcorrection3);
                ZDCenergyHAD3->Fill(var_zdcE3[l],fac_lumi3*effcorrection3);
	     }
	  }

	  CastorSumE3->Fill(var_CastorRecHit3[0],fac_lumi3*effcorrection3);

          double vertexT=sqrt(pow(var_vtxX3[label_vertex],2)+pow(var_vtxY3[label_vertex],2))-0.3633;
          VtxT3->Fill(vertexT,fac_lumi3*effcorrection3);

          double eta_pair=0.5*TMath::Log((double)((var_p3[pair1]+var_p3[pair2]+var_pz3[pair1]+var_pz3[pair2])/(var_p3[pair1]+var_p3[pair2]-var_pz3[pair1]-var_pz3[pair2])));
          double rap_pair=0.5*TMath::Log((double)((sqrt(var_p3[pair1]*var_p3[pair1]+0.1057*0.1057)+sqrt(var_p3[pair2]*var_p3[pair2]+0.1057*0.1057)+var_pz3[pair1]+var_pz3[pair2])/(sqrt(var_p3[pair1]*var_p3[pair1]+0.1057*0.1057)+sqrt(var_p3[pair2]*var_p3[pair2]+0.1057*0.1057)-var_pz3[pair1]-var_pz3[pair2])));
          double pt_pair =sqrt((var_px3[pair1]+var_px3[pair2])*(var_px3[pair1]+var_px3[pair2])+(var_py3[pair1]+var_py3[pair2])*(var_py3[pair1]+var_py3[pair2]));
          etaPair3->Fill(rap_pair,fac_lumi3*effcorrection3);
          pTPair3->Fill(pt_pair,fac_lumi3*effcorrection3);

          if(var_charge3[pair1]>0) {phiSingleP3->Fill(var_phi3[pair1]/pi,fac_lumi3*effcorrection3);etaSingleP3->Fill(var_eta3[pair1],fac_lumi3*effcorrection3);pTSingleP3->Fill(var_pt3[pair1],fac_lumi3*effcorrection3);}
          else {phiSingleM3->Fill(var_phi3[pair1]/pi,fac_lumi3*effcorrection3);etaSingleM3->Fill(var_eta3[pair1],fac_lumi3*effcorrection3);pTSingleM3->Fill(var_pt3[pair1],fac_lumi3*effcorrection3);}
          if(var_charge3[pair2]>0) {phiSingleP3->Fill(var_phi3[pair2]/pi,fac_lumi3*effcorrection3);etaSingleP3->Fill(var_eta3[pair2],fac_lumi3*effcorrection3);pTSingleP3->Fill(var_pt3[pair2],fac_lumi3*effcorrection3);}
          else {phiSingleM3->Fill(var_phi3[pair2]/pi,fac_lumi3*effcorrection3);etaSingleM3->Fill(var_eta3[pair2],fac_lumi3*effcorrection3);pTSingleM3->Fill(var_pt3[pair2],fac_lumi3*effcorrection3);}

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"InelInel :"<<endl;
cout<<"  # Dimuon events = "<<filter3Events<<endl;


  int filter4Gen(0);
  int filter4Track(0);
  double filter4Events(0);
  for(Int_t i = 0;i < NUM4;i++){
      	t4->GetEntry(i);
	int pair1 = var_Pair4[0]; int pair2 = var_Pair4[1];
        int muID1 = var_idA4[pair1];
        int muID2 = var_idA4[pair2];
	int muAng1 = var_idB4[pair1];
        int muAng2 = var_idB4[pair2];
	double effcorrection4 = (var_eff4[pair1]*var_eff4[pair2]*doublemuopenfractionallumi+doublemuopentightfractionallumi);
        int hlt_pass = hlt_d4[0];
	int nPrimVtx = var_nvtx4[0];
	int nTrack=var_nTrack4[0];
	int nTrackQual=var_nTrackQual4[0]; 
	int nCalo=var_ncalo4[0];
	int label_vertex(99);
//	cout<<"--------------------"<<var_event4[0]<<"-----------------------"<<endl;

       if(PassesTrigger(hlt_pass,1) == false) 
                continue;  

	if(PassesKinematicCuts(var_mass4[0],var_pt4[pair1],var_pt4[pair2],var_eta4[pair1],var_eta4[pair2]) == false) 
		continue;

        if(nPrimVtx>=1){
          double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack4,var_vtxZ4,var_vtxmumu4);           
          for(Int_t j=0; j<nPrimVtx; j++){
                if(PassesVertexSelection(var_vtxTrack4[j],var_vertexChi2_4[j],var_vertexNdf4[j],distance_vertex_z,var_vtxZ4[j],var_vtxmumu4[j])
                   && !(techBit4[0][0]==1))   
			{label_vertex=j;}
          }
        }

	if(label_vertex!=99
           && PassesMuonID(var_tracker4[pair1], muAng1, var_global4[pair1], var_tracker4[pair2], muAng2, var_global4[pair2], var_nhitsTrack4[pair1], var_nhitsTrack4[pair2])    
           && PassesDptCut(var_dpt4[0])   
           && PassesDphiCut(var_dphi4[0]/pi)  
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
                if(FailsTrackDistanceVeto(var_TrackD4[j], var_TrackQuality4[j])){   
                  filter4Track++;
		  nTrackExclu++;
		}
            }

    	    int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
	    for(Int_t k=0; k<nCalo; k++){
	       if(var_calodR4[k]>dRcone){
                                if(var_caloId4[k]==2 && var_caloEn4[k]>EBThresh) nEB++;
                                if(var_caloId4[k]==3 && var_caloEn4[k]>EEThresh) nEE++;
                                if(var_caloId4[k]==4 && var_caloEn4[k]>HBThresh) nHB++;
                                if(var_caloId4[k]==5 && var_caloEn4[k]>HEThresh) nHE++;
                                if(var_caloId4[k]==1 && var_caloEta4[k]>0 && var_caloEn4[k]>HFPlusThresh) nHFp++;
                                if(var_caloId4[k]==1 && var_caloEta4[k]<0 && var_caloEn4[k]>HFMinusThresh) nHFm++;

                                if(var_caloId4[k]==6){
                                        if(var_caloHadE4[k]>HBThresh) nHB++;
                                        if(var_caloEmE4[k]>EBThresh) nEB++;}
                                if(var_caloId4[k]==7){
                                        if(var_caloHadE4[k]>HEThresh) nHE++;
                                        if(var_caloEmE4[k]>EBThresh) nEB++;}
                                if(var_caloId4[k]==8){
                                        if(var_caloHadE4[k]>HEThresh) nHE++;
                                        if(var_caloEmE4[k]>EEThresh) nEE++;}
	       }
	    }

        if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus4[0],var_zdcEmPlus4[0],var_zdcHadMinus4[0],var_zdcHadPlus4[0]))
	{
	        nTower4->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi4*effcorrection4);
	}

	if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm) 
		&& PassesZDCVeto(var_zdcEmMinus4[0],var_zdcEmPlus4[0],var_zdcHadMinus4[0],var_zdcHadPlus4[0]))
	){ 
          filter4Events+=fac_lumi4*effcorrection4;

	  hEB4->Fill(nEB,fac_lumi4*effcorrection4);
          hEE4->Fill(nEE,fac_lumi4*effcorrection4);
          hHB4->Fill(nHB,fac_lumi4*effcorrection4);
          hHE4->Fill(nHE,fac_lumi4*effcorrection4);
          hHFp4->Fill(nHFp,fac_lumi4*effcorrection4);
          hHFm4->Fill(nHFm,fac_lumi4*effcorrection4);
//	  nTower4->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi4*effcorrection4);
	  nTrack4->Fill(nTrackExclu,fac_lumi4*effcorrection4);

	  MuMuMass4->Fill(var_mass4[0],fac_lumi4*effcorrection4);
          MuMuMassUps4->Fill(var_mass4[0],fac_lumi4*effcorrection4);
          MuMuMassJpsi4->Fill(var_mass4[0],fac_lumi4*effcorrection4);
          MuMudpt4->Fill(var_dpt4[0],fac_lumi4*effcorrection4);
          MuMudphi4->Fill(var_dphi4[0]/pi,fac_lumi4*effcorrection4);
          double symdphi4 = 1 - fabs(var_phi4[pair1]-var_phi4[pair2])/pi;    
          MuMuSymdphi4->Fill(symdphi4,fac_lumi4*effcorrection4); 
	//          MuMudeta4->Fill(fabs(var_eta4[pair1]+var_eta4[pair2]),fac_lumi4*effcorrection4);
	  MuMuvtxXY4->Fill(sqrt(var_MuMuvtxX4[0]*var_MuMuvtxX4[0]+var_MuMuvtxY4[0]*var_MuMuvtxY4[0]),fac_lumi4*effcorrection4);

	  ZDCemminus4->Fill(var_zdcEmMinus4[0],fac_lumi4*effcorrection4); 
	  ZDChadminus4->Fill(var_zdcHadMinus4[0],fac_lumi4*effcorrection4);
          ZDCemminus4->Fill(var_zdcEmMinus4[0],fac_lumi4*effcorrection4); 
          ZDChadminus4->Fill(var_zdcHadMinus4[0],fac_lumi4*effcorrection4);

          TLorentzVector mu41, mu42, dimuon4;  
          mu41.SetPtEtaPhiM(var_pt4[pair1],var_eta4[pair1],var_phi4[pair1],0.1057);    
          mu42.SetPtEtaPhiM(var_pt4[pair2],var_eta4[pair2],var_phi4[pair2],0.1057);     
          dimuon4 = mu41 + mu42;  
          Tdist4->Fill(dimuon4.Pt()*dimuon4.Pt(),fac_lumi4*effcorrection4);   
          MuMudeta4->Fill(mu41.Angle(mu42.Vect()),fac_lumi4*effcorrection4);  

	  for(Int_t l=0; l<var_nZDC4[0]; l++){
	     if(var_zdcsection4[l]==1 && var_zdcE4[l]>ZDCemThresh){ 
                ZDCtime4->Fill(var_zdcTime4[l],fac_lumi4*effcorrection4);
                ZDCenergyEM4->Fill(var_zdcE4[l],fac_lumi4*effcorrection4);
	     }
	     if(var_zdcsection4[l]==2 && var_zdcE4[l]>ZDChadThresh){
		ZDCtime4->Fill(var_zdcTime4[l],fac_lumi4*effcorrection4);
                ZDCenergyHAD4->Fill(var_zdcE4[l],fac_lumi4*effcorrection4);
	     }
	  }

	  CastorSumE4->Fill(var_CastorRecHit4[0],fac_lumi4*effcorrection4);

          double vertexT=sqrt(pow(var_vtxX4[label_vertex],2)+pow(var_vtxY4[label_vertex],2))-0.3633;
          VtxT4->Fill(vertexT,fac_lumi4*effcorrection4);

          double eta_pair=0.5*TMath::Log((double)((var_p4[pair1]+var_p4[pair2]+var_pz4[pair1]+var_pz4[pair2])/(var_p4[pair1]+var_p4[pair2]-var_pz4[pair1]-var_pz4[pair2])));
          double rap_pair=0.5*TMath::Log((double)((sqrt(var_p4[pair1]*var_p4[pair1]+0.1057*0.1057)+sqrt(var_p4[pair2]*var_p4[pair2]+0.1057*0.1057)+var_pz4[pair1]+var_pz4[pair2])/(sqrt(var_p4[pair1]*var_p4[pair1]+0.1057*0.1057)+sqrt(var_p4[pair2]*var_p4[pair2]+0.1057*0.1057)-var_pz4[pair1]-var_pz4[pair2])));
          double pt_pair =sqrt((var_px4[pair1]+var_px4[pair2])*(var_px4[pair1]+var_px4[pair2])+(var_py4[pair1]+var_py4[pair2])*(var_py4[pair1]+var_py4[pair2]));
          etaPair4->Fill(rap_pair,fac_lumi4*effcorrection4);
          pTPair4->Fill(pt_pair,fac_lumi4*effcorrection4);

          if(var_charge4[pair1]>0) {phiSingleP4->Fill(var_phi4[pair1]/pi,fac_lumi4*effcorrection4);etaSingleP4->Fill(var_eta4[pair1],fac_lumi4*effcorrection4);pTSingleP4->Fill(var_pt4[pair1],fac_lumi4*effcorrection4);}
          else {phiSingleM4->Fill(var_phi4[pair1]/pi,fac_lumi4*effcorrection4);etaSingleM4->Fill(var_eta4[pair1],fac_lumi4*effcorrection4);pTSingleM4->Fill(var_pt4[pair1],fac_lumi4*effcorrection4);}
          if(var_charge4[pair2]>0) {phiSingleP4->Fill(var_phi4[pair2]/pi,fac_lumi4*effcorrection4);etaSingleP4->Fill(var_eta4[pair2],fac_lumi4*effcorrection4);pTSingleP4->Fill(var_pt4[pair2],fac_lumi4*effcorrection4);}
          else {phiSingleM4->Fill(var_phi4[pair2]/pi,fac_lumi4*effcorrection4);etaSingleM4->Fill(var_eta4[pair2],fac_lumi4*effcorrection4);pTSingleM4->Fill(var_pt4[pair2],fac_lumi4*effcorrection4);}

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"Upsilon :"<<endl;
cout<<"  # Dimuon events = "<<filter4Events<<endl;


  int filter5Gen(0);
  int filter5Track(0);
  double filter5Events(0);
  for(Int_t i = 0;i < NUM5;i++){
        t5->GetEntry(i);
        int pair1 = var_Pair5[0]; int pair2 = var_Pair5[1];
        int muID1 = var_idA5[pair1];
        int muID2 = var_idA5[pair2];
        int muAng1 = var_idB5[pair1];
        int muAng2 = var_idB5[pair2];
	double effcorrection5 = (var_eff5[pair1]*var_eff5[pair2]*doublemuopenfractionallumi+doublemuopentightfractionallumi);
        int hlt_pass = hlt_d5[0];
        int nPrimVtx = var_nvtx5[0];
        int nTrack=var_nTrack5[0];
        int nTrackQual=var_nTrackQual5[0];
        int nCalo=var_ncalo5[0];
        int label_vertex(99);
//      cout<<"--------------------"<<var_event5[0]<<"-----------------------"<<endl;

       if(PassesTrigger(hlt_pass,1) == false) 
                continue;  

        if(PassesKinematicCuts(var_mass5[0],var_pt5[pair1],var_pt5[pair2],var_eta5[pair1],var_eta5[pair2]) == false)  
                continue; 

        if(nPrimVtx>=1){
          double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack5,var_vtxZ5,var_vtxmumu5);           
          for(Int_t j=0; j<nPrimVtx; j++){
                if(PassesVertexSelection(var_vtxTrack5[j],var_vertexChi2_5[j],var_vertexNdf5[j],distance_vertex_z,var_vtxZ5[j],var_vtxmumu5[j])
                   && !(techBit5[0][0]==1))   
			{label_vertex=j;}
          }
        }

        if(label_vertex!=99
           && PassesMuonID(var_tracker5[pair1], muAng1, var_global5[pair1], var_tracker5[pair2], muAng2, var_global5[pair2], var_nhitsTrack5[pair1], var_nhitsTrack5[pair2])     
           && PassesDptCut(var_dpt5[0])   
           && PassesDphiCut(var_dphi5[0]/pi)  
                ) {
            int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){
                if(FailsTrackDistanceVeto(var_TrackD5[j], var_TrackQuality5[j])){   
                  filter5Track++;
		  nTrackExclu++;
                }
            }

            int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
            for(Int_t k=0; k<nCalo; k++){
               if(var_calodR5[k]>dRcone){
                                if(var_caloId5[k]==2 && var_caloEn5[k]>EBThresh) nEB++;
                                if(var_caloId5[k]==3 && var_caloEn5[k]>EEThresh) nEE++;
                                if(var_caloId5[k]==4 && var_caloEn5[k]>HBThresh) nHB++;
                                if(var_caloId5[k]==5 && var_caloEn5[k]>HEThresh) nHE++;
                                if(var_caloId5[k]==1 && var_caloEta5[k]>0 && var_caloEn5[k]>HFPlusThresh) nHFp++;
                                if(var_caloId5[k]==1 && var_caloEta5[k]<0 && var_caloEn5[k]>HFMinusThresh) nHFm++;

                                if(var_caloId5[k]==6){
                                        if(var_caloHadE5[k]>HBThresh) nHB++;
                                        if(var_caloEmE5[k]>EBThresh) nEB++;}
                                if(var_caloId5[k]==7){
                                        if(var_caloHadE5[k]>HEThresh) nHE++;
                                        if(var_caloEmE5[k]>EBThresh) nEB++;}
                                if(var_caloId5[k]==8){
                                        if(var_caloHadE5[k]>HEThresh) nHE++;
                                        if(var_caloEmE5[k]>EEThresh) nEE++;}
                      }
            }

        if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus5[0],var_zdcEmPlus5[0],var_zdcHadMinus5[0],var_zdcHadPlus5[0]))
	{
	          nTower5->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi5*effcorrection5);
	}

        if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm)
		&& PassesZDCVeto(var_zdcEmMinus5[0],var_zdcEmPlus5[0],var_zdcHadMinus5[0],var_zdcHadPlus5[0]))
	){
          filter5Events+=fac_lumi5*effcorrection5;

          hEB5->Fill(nEB,fac_lumi5*effcorrection5);
          hEE5->Fill(nEE,fac_lumi5*effcorrection5);
          hHB5->Fill(nHB,fac_lumi5*effcorrection5);
          hHE5->Fill(nHE,fac_lumi5*effcorrection5);
          hHFp5->Fill(nHFp,fac_lumi5*effcorrection5);
          hHFm5->Fill(nHFm,fac_lumi5*effcorrection5);
//          nTower5->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumi5*effcorrection5);
          nTrack5->Fill(nTrackExclu,fac_lumi5*effcorrection5);

          MuMuMass5->Fill(var_mass5[0],fac_lumi5*effcorrection5);
          MuMuMassUps5->Fill(var_mass5[0],fac_lumi5*effcorrection5);
          MuMuMassJpsi5->Fill(var_mass5[0],fac_lumi5*effcorrection5);
          MuMudpt5->Fill(var_dpt5[0],fac_lumi5*effcorrection5);
          MuMudphi5->Fill(var_dphi5[0]/pi,fac_lumi5*effcorrection5);
          double symdphi5 = 1 - fabs(var_phi5[pair1]-var_phi5[pair2])/pi;    
          MuMuSymdphi5->Fill(symdphi5,fac_lumi5*effcorrection5); 
	//          MuMudeta5->Fill(fabs(var_eta5[pair1]+var_eta5[pair2]),fac_lumi5*effcorrection5);
	  MuMuvtxXY5->Fill(sqrt(var_MuMuvtxX5[0]*var_MuMuvtxX5[0]+var_MuMuvtxY5[0]*var_MuMuvtxY5[0]),fac_lumi5*effcorrection5);

          TLorentzVector mu51, mu52, dimuon5;   
          mu51.SetPtEtaPhiM(var_pt5[pair1],var_eta5[pair1],var_phi5[pair1],0.1057);     
          mu52.SetPtEtaPhiM(var_pt5[pair2],var_eta5[pair2],var_phi5[pair2],0.1057);      
          dimuon5 = mu51 + mu52;   
          Tdist5->Fill(dimuon5.Pt()*dimuon5.Pt(),fac_lumi5*effcorrection5);     
          MuMudeta5->Fill(mu51.Angle(mu52.Vect()),fac_lumi5*effcorrection5);  

          ZDCemminus5->Fill(var_zdcEmMinus5[0],fac_lumi5*effcorrection5); ZDCemplus5->Fill(var_zdcEmPlus5[0],fac_lumi5*effcorrection5);
          ZDChadminus5->Fill(var_zdcHadMinus5[0],fac_lumi5*effcorrection5); ZDChadplus5->Fill(var_zdcHadPlus5[0],fac_lumi5*effcorrection5);
          for(Int_t l=0; l<var_nZDC5[0]; l++){
             if(var_zdcsection5[l]==1 && var_zdcE5[l]>ZDCemThresh){
                ZDCtime5->Fill(var_zdcTime5[l],fac_lumi5*effcorrection5); ZDCenergyEM5->Fill(var_zdcE5[l],fac_lumi5*effcorrection5);
             }
             if(var_zdcsection5[l]==2 && var_zdcE5[l]>ZDChadThresh){
                ZDCtime5->Fill(var_zdcTime5[l],fac_lumi5*effcorrection5); ZDCenergyHAD5->Fill(var_zdcE5[l],fac_lumi5*effcorrection5);
             }
          }

          CastorSumE5->Fill(var_CastorRecHit5[0],fac_lumi5*effcorrection5);

          double vertexT=sqrt(pow(var_vtxX5[label_vertex],2)+pow(var_vtxY5[label_vertex],2))-0.3633;
          VtxT5->Fill(vertexT,fac_lumi5*effcorrection5);

          double eta_pair=0.5*TMath::Log((double)((var_p5[pair1]+var_p5[pair2]+var_pz5[pair1]+var_pz5[pair2])/(var_p5[pair1]+var_p5[pair2]-var_pz5[pair1]-var_pz5[pair2])));
	  double rap_pair=0.5*TMath::Log((double)((sqrt(var_p5[pair1]*var_p5[pair1]+0.1057*0.1057)+sqrt(var_p5[pair2]*var_p5[pair2]+0.1057*0.1057)+var_pz5[pair1]+var_pz5[pair2])/(sqrt(var_p5[pair1]*var_p5[pair1]+0.1057*0.1057)+sqrt(var_p5[pair2]*var_p5[pair2]+0.1057*0.1057)-var_pz5[pair1]-var_pz5[pair2])));
          double pt_pair =sqrt((var_px5[pair1]+var_px5[pair2])*(var_px5[pair1]+var_px5[pair2])+(var_py5[pair1]+var_py5[pair2])*(var_py5[pair1]+var_py5[pair2]));
          etaPair5->Fill(rap_pair,fac_lumi5*effcorrection5);
          pTPair5->Fill(pt_pair,fac_lumi5*effcorrection5);

          if(var_charge5[pair1]>0) {phiSingleP5->Fill(var_phi5[pair1]/pi,fac_lumi5*effcorrection5);etaSingleP5->Fill(var_eta5[pair1],fac_lumi5*effcorrection5);pTSingleP5->Fill(var_pt5[pair1],fac_lumi5*effcorrection5);}
          else {phiSingleM5->Fill(var_phi5[pair1]/pi,fac_lumi5*effcorrection5);etaSingleM5->Fill(var_eta5[pair1],fac_lumi5*effcorrection5);pTSingleM5->Fill(var_pt5[pair1],fac_lumi5*effcorrection5);}
          if(var_charge5[pair2]>0) {phiSingleP5->Fill(var_phi5[pair2]/pi,fac_lumi5*effcorrection5);etaSingleP5->Fill(var_eta5[pair2],fac_lumi5*effcorrection5);pTSingleP5->Fill(var_pt5[pair2],fac_lumi5*effcorrection5);}
          else {phiSingleM5->Fill(var_phi5[pair2]/pi,fac_lumi5*effcorrection5);etaSingleM5->Fill(var_eta5[pair2],fac_lumi5*effcorrection5);pTSingleM5->Fill(var_pt5[pair2],fac_lumi5*effcorrection5);}


          } // if nTrack&nCalo if relevant
        }
  }
cout<<"Jpsi :"<<endl;
cout<<"  # Dimuon events = "<<filter5Events<<endl;


  int filter6Gen(0);
  int filter6Track(0);
  double filter6Events_norm(0);
  for(Int_t i = 0;i < NUM6;i++){
        t6->GetEntry(i);
        int pair1 = var_Pair6[0]; int pair2 = var_Pair6[1];
        int muID1 = var_idA6[pair1];
        int muID2 = var_idA6[pair2];
        int muAng1 = var_idB6[pair1];
        int muAng2 = var_idB6[pair2];
        int hlt_pass = hlt_d6[0];
        int nPrimVtx = var_nvtx6[0];
        int nTrack=var_nTrack6[0];
        int nTrackQual=var_nTrackQual6[0];
        int nCalo=var_ncalo6[0];
        int label_vertex(99);
	int bkgNum=var_run6[0];
        double effcorrection6 = (var_eff6[pair1]*var_eff6[pair2]*doublemuopenfractionallumi+doublemuopentightfractionallumi);
//      cout<<"--------------------"<<var_event1[0]<<"-----------------------"<<endl;

       if(PassesTrigger(hlt_pass,1) == false)
                continue;

        if(PassesKinematicCuts(var_mass6[0],var_pt6[pair1],var_pt6[pair2],var_eta6[pair1],var_eta6[pair2]) == false)
                continue;

        if(nPrimVtx>=1){
          double distance_vertex_z = VertexSeparation(nPrimVtx,var_vtxTrack6,var_vtxZ6,var_vtxmumu6);
          for(Int_t j=0; j<nPrimVtx; j++){
                if(PassesVertexSelection(var_vtxTrack6[j],var_vertexChi2_6[j],var_vertexNdf6[j],distance_vertex_z,var_vtxZ6[j],var_vtxmumu6[j]))
                   {label_vertex=j;}
          }
        }
        if(label_vertex!=99
           && PassesMuonID(var_tracker6[pair1], muAng1, var_global6[pair1], var_tracker6[pair2], muAng2, var_global6[pair2], var_nhitsTrack6[pair1], var_nhitsTrack6[pair2])      
	   && PassesDptCut(var_dpt6[0])  
           && PassesDphiCut(var_dphi6[0]/pi) 
                ) {
            int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){
                if(FailsTrackDistanceVeto(var_TrackD6[j], var_TrackQuality6[j])){   
                  filter6Track++;
		  nTrackExclu++;
                }
            }

            int nEB(0),nEE(0),nHB(0),nHE(0),nHFp(0),nHFm(0);
            for(Int_t k=0; k<nCalo; k++){
               if(var_calodR6[k]>0.3){
                                if(var_caloId6[k]==2 && var_caloEn6[k]>EBThresh) nEB++;
                                if(var_caloId6[k]==3 && var_caloEn6[k]>EEThresh) nEE++;
                                if(var_caloId6[k]==4 && var_caloEn6[k]>HBThresh) nHB++;
                                if(var_caloId6[k]==5 && var_caloEn6[k]>HEThresh) nHE++;
                                if(var_caloId6[k]==1 && var_caloEta6[k]>0 && var_caloEn6[k]>HFPlusThresh) nHFp++;
                                if(var_caloId6[k]==1 && var_caloEta6[k]<0 && var_caloEn6[k]>HFMinusThresh) nHFm++;

                                if(var_caloId6[k]==6){
                                        if(var_caloHadE6[k]>1.25) nHB++;
                                        if(var_caloEmE6[k]>0.60) nEB++;}
                                if(var_caloId6[k]==7){
                                        if(var_caloHadE6[k]>1.80) nHE++;
                                        if(var_caloEmE6[k]>0.60) nEB++;}
                                if(var_caloId6[k]==8){
                                        if(var_caloHadE6[k]>1.80) nHE++;
                                        if(var_caloEmE6[k]>2.40) nEE++;}
               }
            }

	if(nTrackExclu<1 && PassesZDCVeto(var_zdcEmMinus6[0],var_zdcEmPlus6[0],var_zdcHadMinus6[0],var_zdcHadPlus6[0]))
	{
		nTower6->Fill(nEB+nEE+nHB+nHE+nHFp+nHFm,fac_lumiBkg[bkgNum]*effcorrection6);
	}

        if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm)
		&& PassesZDCVeto(var_zdcEmMinus6[0],var_zdcEmPlus6[0],var_zdcHadMinus6[0],var_zdcHadPlus6[0])
	){
          filter6Events_norm+=fac_lumiBkg[bkgNum]*effcorrection6;

          hEB6->Fill(nEB,fac_lumiBkg[bkgNum]*effcorrection6); 
          hEE6->Fill(nEE,fac_lumiBkg[bkgNum]*effcorrection6); 
          hHB6->Fill(nHB,fac_lumiBkg[bkgNum]*effcorrection6);
          hHE6->Fill(nHE,fac_lumiBkg[bkgNum]*effcorrection6);
          hHFp6->Fill(nHFp,fac_lumiBkg[bkgNum]*effcorrection6);
          hHFm6->Fill(nHFm,fac_lumiBkg[bkgNum]*effcorrection6);
          nTrack6->Fill(nTrackExclu,fac_lumiBkg[bkgNum]*effcorrection6);

          MuMuMass6->Fill(var_mass6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          MuMuMassUps6->Fill(var_mass6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          MuMuMassJpsi6->Fill(var_mass6[0],fac_lumiBkg[bkgNum]*effcorrection6); 
          MuMudpt6->Fill(var_dpt6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          MuMudphi6->Fill(var_dphi6[0]/pi,fac_lumiBkg[bkgNum]*effcorrection6);
          double symdphi6 = 1 - fabs(var_phi6[pair1]-var_phi6[pair2])/pi;
          MuMuSymdphi6->Fill(symdphi6,fac_lumiBkg[bkgNum]*effcorrection6);
          MuMudeta6->Fill(fabs(var_eta6[pair1]+var_eta6[pair2]),fac_lumiBkg[bkgNum]*effcorrection6);
          MuMuvtxXY6->Fill(sqrt(var_MuMuvtxX6[0]*var_MuMuvtxX6[0]+var_MuMuvtxY6[0]*var_MuMuvtxY6[0]),fac_lumiBkg[bkgNum]*effcorrection6);

          TLorentzVector mu61, mu62, dimuon6;   
          mu61.SetPtEtaPhiM(var_pt6[pair1],var_eta6[pair1],var_phi6[pair1],0.1057);     
          mu62.SetPtEtaPhiM(var_pt6[pair2],var_eta6[pair2],var_phi6[pair2],0.1057);      
          dimuon6 = mu61 + mu62;   
          Tdist6->Fill(dimuon6.Pt()*dimuon6.Pt(),fac_lumiBkg[bkgNum]*effcorrection6);    

          ZDCemminus6->Fill(var_zdcEmMinus6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          ZDCemplus6->Fill(var_zdcEmPlus6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          ZDChadminus6->Fill(var_zdcHadMinus6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          ZDChadplus6->Fill(var_zdcHadPlus6[0],fac_lumiBkg[bkgNum]*effcorrection6);
          for(Int_t l=0; l<var_nZDC6[0]; l++){
             if(var_zdcsection6[l]==1 && var_zdcE6[l]>ZDCemThresh){
                ZDCtime6->Fill(var_zdcTime6[l],fac_lumiBkg[bkgNum]*effcorrection6);
                ZDCenergyEM6->Fill(var_zdcE6[l],fac_lumiBkg[bkgNum]*effcorrection6);
             }
             if(var_zdcsection6[l]==2 && var_zdcE6[l]>ZDChadThresh){
                ZDCtime6->Fill(var_zdcTime6[l],fac_lumiBkg[bkgNum]*effcorrection6);
                ZDCenergyHAD6->Fill(var_zdcE6[l],fac_lumiBkg[bkgNum]*effcorrection6);
             }
          }

          CastorSumE6->Fill(var_CastorRecHit6[0],fac_lumiBkg[bkgNum]*effcorrection6);

          double vertexT=sqrt(pow(var_vtxX6[label_vertex],2)+pow(var_vtxY6[label_vertex],2))-0.3633;
          VtxT6->Fill(vertexT,fac_lumiBkg[bkgNum]*effcorrection6);

          double eta_pair=0.5*TMath::Log((double)((var_p6[pair1]+var_p6[pair2]+var_pz6[pair1]+var_pz6[pair2])/(var_p6[pair1]+var_p6[pair2]-var_pz6[pair1]-var_pz6[pair2])));
          double rap_pair=0.5*TMath::Log((double)((sqrt(var_p6[pair1]*var_p6[pair1]+0.1057*0.1057)+sqrt(var_p6[pair2]*var_p6[pair2]+0.1057*0.1057)+var_pz6[pair1]+var_pz6[pair2])/(sqrt(var_p6[pair1]*var_p6[pair1]+0.1057*0.1057)+sqrt(var_p6[pair2]*var_p6[pair2]+0.1057*0.1057)-var_pz6[pair1]-var_pz6[pair2])));
          double pt_pair =sqrt((var_px6[pair1]+var_px6[pair2])*(var_px6[pair1]+var_px6[pair2])+(var_py6[pair1]+var_py6[pair2])*(var_py6[pair1]+var_py6[pair2]));
          etaPair6->Fill(rap_pair,fac_lumiBkg[bkgNum]*effcorrection6);
          pTPair6->Fill(pt_pair,fac_lumiBkg[bkgNum]*effcorrection6);

          if(var_charge6[pair1]>0) {phiSingleP6->Fill(var_phi6[pair1]/pi,fac_lumiBkg[bkgNum]*effcorrection6);etaSingleP6->Fill(var_eta6[pair1],fac_lumiBkg[bkgNum]*effcorrection6);pTSingleP6->Fill(var_pt6[pair1],fac_lumiBkg[bkgNum]*effcorrection6);}
          else {phiSingleM6->Fill(var_phi6[pair1]/pi,fac_lumiBkg[bkgNum]*effcorrection6);etaSingleM6->Fill(var_eta6[pair1],fac_lumiBkg[bkgNum]*effcorrection6);pTSingleM6->Fill(var_pt6[pair1],fac_lumiBkg[bkgNum]*effcorrection6);}
          if(var_charge6[pair2]>0) {phiSingleP6->Fill(var_phi6[pair2]/pi,fac_lumiBkg[bkgNum]*effcorrection6);etaSingleP6->Fill(var_eta6[pair2],fac_lumiBkg[bkgNum]*effcorrection6);pTSingleP6->Fill(var_pt6[pair2],fac_lumiBkg[bkgNum]*effcorrection6);}
          else {phiSingleM6->Fill(var_phi6[pair2]/pi,fac_lumiBkg[bkgNum]*effcorrection6);etaSingleM6->Fill(var_eta6[pair2],fac_lumiBkg[bkgNum]*effcorrection6);pTSingleM6->Fill(var_pt6[pair2],fac_lumiBkg[bkgNum]*effcorrection6);}

          } // if nTrack&nCalo if relevant

        }
  }
cout<<"Inclusive :"<<endl;
cout<<"  # Dimuon events = "<<filter6Events_norm<<endl;



ci = TColor::GetColor("#ffff99");

DrawOneHistogramBis(sMuMuMass,MuMuMass0,MuMuMass1,MuMuMass2,MuMuMass3,MuMuMass4,MuMuMass5,MuMuMass6,"#mu#mu mass [GeV]","Events/0.5 GeV","MuMuMass_34pb-1_trackExclu2mm.png");
DrawOneHistogramBis(sMuMudeta,MuMudeta0,MuMudeta1,MuMudeta2,MuMudeta3,MuMudeta4,MuMudeta5,MuMudeta6,"#mu#mu 3D opening angle","Events/0.1","MuMu3dangle_34pb-1_trackExclu2mm.png"); 
DrawOneHistogramBis(sMuMudpt,MuMudpt0,MuMudpt1,MuMudpt2,MuMudpt3,MuMudpt4,MuMudpt5,MuMudpt6,"#mu#mu |#Delta p_{T}| [GeV]","Events/0.1 GeV","MuMudpt_34pb-1_trackExclu2mm.png"); 
DrawOneHistogramBis(sMuMudphi,MuMudphi0,MuMudphi1,MuMudphi2,MuMudphi3,MuMudphi4,MuMudphi5,MuMudphi6,"#mu#mu |#Delta #phi / #pi|","Events/0.01","MuMudphi_34pb-1_trackExclu2mm.png"); 
DrawOneHistogramBis(sMuMuSymdphi,MuMuSymdphi0,MuMuSymdphi1,MuMuSymdphi2,MuMuSymdphi3,MuMuSymdphi4,MuMuSymdphi5,MuMuSymdphi6,"1 - |#phi(#mu^{-})-#phi(#mu^{+})| / #pi","Events/0.01","MuMuSymdphi_34pb-1_trackExclu2mm.png");  
DrawOneHistogramBis(setaPair,etaPair0,etaPair1,etaPair2,etaPair3,etaPair4,etaPair5,etaPair6,"#eta(#mu#mu)","Events/0.5","etaPair_34pb-1_trackExclu2mm.png");  
DrawOneHistogramBis(spTPair,pTPair0,pTPair1,pTPair2,pTPair3,pTPair4,pTPair5,pTPair6,"p_{T}(#mu#mu) [GeV]","Events/0.5 GeV","pTPair_34pb-1_trackExclu2mm.png");   
DrawOneHistogramBis(setaSingleP,etaSingleP0,etaSingleP1,etaSingleP2,etaSingleP3,etaSingleP4,etaSingleP5,etaSingleP6,"#eta(#mu^{+})","Events/0.5","etaSingleP_34pb-1_trackExclu2mm.png");  
DrawOneHistogramBis(setaSingleM,etaSingleM0,etaSingleM1,etaSingleM2,etaSingleM3,etaSingleM4,etaSingleM5,etaSingleM6,"#eta(#mu^{-})","Events/0.5","etaSingleM_34pb-1_trackExclu2mm.png"); 
DrawOneHistogramBis(spTSingleP,pTSingleP0,pTSingleP1,pTSingleP2,pTSingleP3,pTSingleP4,pTSingleP5,pTSingleP6,"p_{T}(#mu^{+}) [GeV]","Events/0.5 GeV","pTSingleP_34pb-1_trackExclu2mm.png");  
DrawOneHistogramBis(spTSingleM,pTSingleM0,pTSingleM1,pTSingleM2,pTSingleM3,pTSingleM4,pTSingleM5,pTSingleM6,"p_{T}(#mu^{-}) [GeV]","Events/0.5 GeV","pTSingleM_34pb-1_trackExclu2mm.png");
DrawOneHistogramBis(sphiSingleP,phiSingleP0,phiSingleP1,phiSingleP2,phiSingleP3,phiSingleP4,phiSingleP5,phiSingleP6,"#phi(#mu^{+})","Events/0.5","phiSingleP_34pb-1_trackExclu2mm.png");    
DrawOneHistogramBis(sphiSingleM,phiSingleM0,phiSingleM1,phiSingleM2,phiSingleM3,phiSingleM4,phiSingleM5,phiSingleM6,"#phi(#mu^{-})","Events/0.5","phiSingleM_34pb-1_trackExclu2mm.png");
//DrawOneHistogram(sMuMuMassUps,MuMuMassUps0,MuMuMassUps1,MuMuMassUps2,MuMuMassUps3,MuMuMassUps4,MuMuMassUps5,"#mu#mu mass","Events/0.1 GeV","MuMuMassUps_34pb-1_trackExclu2mm.png");
//DrawOneHistogramBis(sMuMuMassJpsi,MuMuMassJpsi0,MuMuMassJpsi1,MuMuMassJpsi2,MuMuMassJpsi3,MuMuMassJpsi4,MuMuMassJpsi5,MuMuMassJpsi6,"#mu#mu mass","Events/0.4","MuMuMassJpsi_34pb-1_trackExclu2mm.png");
//DrawOneHistogram(sTrack,nTrack0,nTrack1,nTrack2,nTrack3,nTrack4,nTrack5,"# track (|d|<2mm)","Events/0.1 GeV","nTrack_34pb-1_trackExclu2mm.png");
DrawOneHistogramBis(sVtxT,VtxT0,VtxT1,VtxT2,VtxT3,VtxT4,VtxT5,VtxT6,"#mu#mu transverse vtx [cm]","Events/0.1 cm","VtxT_34pb-1_trackExclu2mm.png");


// Save into histo
//Draw
if(0){
TCanvas *PtSpectrumEl = new TCanvas("PtSpectrumEl","Tracks",800,500);
   PtSpectrumEl->SetFillColor(0);
   PtSpectrumEl->SetBorderMode(0);
   PtSpectrumEl->SetBorderSize(2);
   PtSpectrumEl->SetFrameBorderMode(0);

//   leg->Draw();

PtSpectrumEl->cd(1);
nTrack0->Sumw2();
nTrack0->SetLineWidth(2);
nTrack0->SetMarkerStyle(20);
nTrack1->SetFillColor(ci);
nTrack2->SetFillColor(30);
nTrack2->SetFillStyle(3001);
nTrack3->SetFillColor(30);
nTrack4->SetFillColor(38);
nTrack5->SetFillColor(903);
sTrack->Add(nTrack3);
sTrack->Add(nTrack2);
sTrack->Add(nTrack1);
sTrack->Add(nTrack4);
sTrack->Add(nTrack5);
sTrack.SetMaximum(sTrack->GetMaximum() * 1.5);
sTrack->Draw();
sTrack->GetXaxis()->SetTitle("n Extra Tracks");
sTrack->GetYaxis()->SetTitle("# events");
nTrack0->Draw("same");
}

if(0){
//------------
TCanvas *Vertex = new TCanvas("Vertex","Vertex",800,500);
   Vertex->SetFillColor(0);
   Vertex->SetBorderMode(0);
   Vertex->SetBorderSize(2);
   Vertex->SetFrameBorderMode(0);

Vertex->cd(1);
MuMuvtxXY0->Sumw2();
MuMuvtxXY0->SetLineWidth(2);
MuMuvtxXY0->SetMarkerStyle(20);
MuMuvtxXY1->SetFillColor(ci);
MuMuvtxXY2->SetFillColor(30);
MuMuvtxXY2->SetFillStyle(3001);
MuMuvtxXY3->SetFillColor(30);
MuMuvtxXY4->SetFillColor(38);
MuMuvtxXY5->SetFillColor(903);
sMuMuvtxXY->Add(MuMuvtxXY3);
sMuMuvtxXY->Add(MuMuvtxXY2);
sMuMuvtxXY->Add(MuMuvtxXY1);
sMuMuvtxXY->Add(MuMuvtxXY4);
sMuMuvtxXY->Add(MuMuvtxXY5);
sMuMuvtxXY->Draw();
sMuMuvtxXY->GetXaxis()->SetTitle("transverse vertex position [cm]");
sMuMuvtxXY->GetYaxis()->SetTitle("# events");
MuMuvtxXY0->Draw("same");
}

if(0){
//---------------------------
TCanvas *Calo = new TCanvas("Calo","Calorimeter",800,500);
   Calo->SetFillColor(0);
   Calo->SetBorderMode(0);
   Calo->SetBorderSize(2);
   Calo->SetFrameBorderMode(0);
Calo->Divide(3,2);
Calo->cd(4);
hEB0->Sumw2();
hEB0->SetLineWidth(2);
hEB0->SetMarkerStyle(20);
hEB1->SetFillColor(ci);
hEB2->SetFillColor(30);
hEB2->SetFillStyle(3001);
hEB3->SetFillColor(30);
hEB4->SetFillColor(38);
hEB5->SetFillColor(903);
sEB->Add(hEB3);
sEB->Add(hEB2);
sEB->Add(hEB1);
sEB->Add(hEB4);
sEB->Add(hEB5);
sEB.SetMaximum(sEB->GetMaximum() * 1.5); 
sEB->Draw();
sEB->GetXaxis()->SetTitle("n EB towers (E>0.60GeV)");
sEB->GetYaxis()->SetTitle("# events");
hEB0->Draw("same");

Calo->cd(5);
hEE0->Sumw2();
hEE0->SetLineWidth(2);
hEE0->SetMarkerStyle(20);
hEE1->SetFillColor(ci);
hEE2->SetFillColor(30);
hEE2->SetFillStyle(3001);
hEE3->SetFillColor(30);
hEE4->SetFillColor(38);
hEE5->SetFillColor(903);
sEE->Add(hEE3);
sEE->Add(hEE2);
sEE->Add(hEE1);
sEE->Add(hEE4);
sEE->Add(hEE5);
sEE.SetMaximum(sEE->GetMaximum() * 1.5); 
sEE->Draw();
sEE->GetXaxis()->SetTitle("n EE towers (E>2.40GeV)");
sEE->GetYaxis()->SetTitle("# events");
hEE0->Draw("same");

Calo->cd(1);
hHB0->Sumw2();
hHB0->SetLineWidth(2);
hHB0->SetMarkerStyle(20);
hHB1->SetFillColor(ci);
hHB2->SetFillColor(30);
hHB2->SetFillStyle(3001);
hHB3->SetFillColor(30);
hHB4->SetFillColor(38);
hHB5->SetFillColor(903);
hHB5->SetFillColor(903);
sHB->Add(hHB3);
sHB->Add(hHB2);
sHB->Add(hHB1);
sHB->Add(hHB4);
sHB->Add(hHB5);
sHB.SetMaximum(sHB->GetMaximum() * 1.5); 
sHB->Draw();
sHB->GetXaxis()->SetTitle("n HB towers (E>1.25GeV)");
sHB->GetYaxis()->SetTitle("# events");
hHB0->Draw("same");

Calo->cd(2);
hHE0->Sumw2();
hHE0->SetLineWidth(2);
hHE0->SetMarkerStyle(20);
hHE1->SetFillColor(ci);
hHE2->SetFillColor(30);
hHE2->SetFillStyle(3001);
hHE3->SetFillColor(30);
hHE4->SetFillColor(38);
hHE5->SetFillColor(903);
sHE->Add(hHE3);
sHE->Add(hHE2);
sHE->Add(hHE1);
sHE->Add(hHE4);
sHE->Add(hHE5);
sHE.SetMaximum(sHE->GetMaximum() * 1.5); 
sHE->Draw();
sHE->GetXaxis()->SetTitle("n HE towers (E>1.90GeV)");
sHE->GetYaxis()->SetTitle("# events");
hHE0->Draw("same");

Calo->cd(3);
hHFp0->Sumw2();
hHFp0->SetLineWidth(2);
hHFp0->SetMarkerStyle(20);
hHFp1->SetFillColor(ci);
hHFp2->SetFillColor(30);
hHFp2->SetFillStyle(3001);
hHFp3->SetFillColor(30);
hHFp4->SetFillColor(38);
hHFp5->SetFillColor(903);
sHFp->Add(hHFp3);
sHFp->Add(hHFp2);
sHFp->Add(hHFp1);
sHFp->Add(hHFp4);
sHFp->Add(hHFp5);
sHFp.SetMaximum(sHFp->GetMaximum() * 1.5); 
sHFp->Draw();
sHFp->GetXaxis()->SetTitle("n HF+ towers (E>4.5GeV)");
sHFp->GetYaxis()->SetTitle("# events");
hHFp0->Draw("same");

Calo->cd(6);
hHFm0->Sumw2();
hHFm0->SetLineWidth(2);
hHFm0->SetMarkerStyle(20);
hHFm1->SetFillColor(ci);
hHFm2->SetFillColor(30);
hHFm2->SetFillStyle(3001);
hHFm3->SetFillColor(30);
hHFm4->SetFillColor(38);
hHFm5->SetFillColor(903);
sHFm->Add(hHFm3);
sHFm->Add(hHFm2);
sHFm->Add(hHFm1);
sHFm->Add(hHFm4);
sHFm->Add(hHFm5);
sHFm.SetMaximum(sHFm->GetMaximum() * 1.5); 
sHFm->Draw();
//Calo->SaveAs("Calo_851nb.eps"); 
sHFm->GetXaxis()->SetTitle("n HF- towers (E>4.0GeV)");
sHFm->GetYaxis()->SetTitle("# events");
hHFm0->Draw("same");

TCanvas *Calo2 = new TCanvas("Calo2","Calorimeter 2",800,500);
   Calo2->SetFillColor(0);
   Calo2->SetBorderMode(0);
   Calo2->SetBorderSize(2);
   Calo2->SetFrameBorderMode(0);
Calo2->cd(1);
nTower0->Sumw2();
nTower0->SetLineWidth(2);
nTower0->SetMarkerStyle(20);
nTower1->SetFillColor(ci);
nTower2->SetFillColor(30);
nTower2->SetFillStyle(3001);
nTower3->SetFillColor(30);
nTower4->SetFillColor(38);
nTower5->SetFillColor(903);
sTower->Add(nTower3);
sTower->Add(nTower2);
sTower->Add(nTower1);
sTower->Add(nTower4);
sTower->Add(nTower5);
sTower.SetMaximum(sTower->GetMaximum() * 1.5); 
sTower->Draw();
sTower->GetXaxis()->SetTitle("n towers (E>E_{Threshold})");
sTower->GetYaxis()->SetTitle("# events");
nTower0->Draw("same");

//Calo2->SaveAs("Calo2_851nb.eps");  
}

if(0){
TCanvas *Upsilon_Jpsi = new TCanvas("Upsilon_Jpsi","Upsilon J/psi",800,500);
   Upsilon_Jpsi->SetFillColor(0);
   Upsilon_Jpsi->SetBorderMode(0);
   Upsilon_Jpsi->SetBorderSize(2);
   Upsilon_Jpsi->SetFrameBorderMode(0);

Upsilon_Jpsi->Divide(2,1);
Upsilon_Jpsi->cd(1);
MuMuMassUps0->Sumw2();
MuMuMassUps0->SetLineWidth(2);
MuMuMassUps0->SetMarkerStyle(20); 
MuMuMassUps1->SetFillColor(0);MuMuMassUps1->SetLineWidth(3);MuMuMassUps1->SetLineColor(4); 
MuMuMassUps2->SetFillColor(2); 
MuMuMassUps3->SetFillColor(5); 
MuMuMassUps4->SetFillColor(38);
MuMuMassUps5->SetFillColor(903);
sMuMuMassUps->Add(MuMuMassUps3);
sMuMuMassUps->Add(MuMuMassUps2);
sMuMuMassUps->Add(MuMuMassUps1);
sMuMuMassUps->Add(MuMuMassUps4);
sMuMuMassUps->Add(MuMuMassUps5);
sMuMuMassUps.SetMaximum(sMuMuMassUps->GetMaximum() * 5); 
sMuMuMassUps->Draw();
sMuMuMassUps->GetXaxis()->SetTitle("#mu#mu mass [GeV]");
sMuMuMassUps->GetYaxis()->SetTitle("# events / 0.05 GeV");
MuMuMassUps0->Draw("esame");

Upsilon_Jpsi->cd(2);
MuMuMassJpsi0->Sumw2();
MuMuMassJpsi0->SetLineWidth(2);
MuMuMassJpsi0->SetMarkerStyle(20); 
MuMuMassJpsi1->SetFillColor(0);MuMuMassJpsi1->SetLineWidth(3);MuMuMassJpsi1->SetLineColor(4); 
MuMuMassJpsi2->SetFillColor(2); 
//MuMuMassJpsi2->SetFillStyle(3001); 
MuMuMassJpsi3->SetFillColor(5); 
MuMuMassJpsi4->SetFillColor(38);
MuMuMassJpsi5->SetFillColor(903);
sMuMuMassJpsi->Add(MuMuMassJpsi3);
sMuMuMassJpsi->Add(MuMuMassJpsi2);
sMuMuMassJpsi->Add(MuMuMassJpsi1);
sMuMuMassJpsi->Add(MuMuMassJpsi4);
sMuMuMassJpsi->Add(MuMuMassJpsi5);
sMuMuMassJpsi.SetMaximum(sMuMuMassJpsi->GetMaximum() * 5);  
sMuMuMassJpsi->Draw();
sMuMuMassJpsi->GetXaxis()->SetTitle("#mu#mu mass [GeV]");
sMuMuMassJpsi->GetYaxis()->SetTitle("# events / 0.025 GeV");
MuMuMassJpsi0->Draw("esame");
//Upsilon_Jpsi->SaveAs("Upsilon_Jpsi_851nb.eps");  
}

if(0){
TCanvas *Kinematic1 = new TCanvas("Kinematic1","Kinematic MuMu",800,500);
   Kinematic1->SetFillColor(0);
   Kinematic1->SetBorderMode(0);
   Kinematic1->SetBorderSize(2);
   Kinematic1->SetFrameBorderMode(0);

Kinematic1->Divide(2,1);
Kinematic1->cd(1);
MuMuMass0->Sumw2();
MuMuMass0->SetLineWidth(2);
MuMuMass0->SetMarkerStyle(20);
MuMuMass1->SetFillColor(0);MuMuMass1->SetLineWidth(3);MuMuMass1->SetLineColor(4);
MuMuMass2->SetFillColor(2);
//MuMuMass2->SetFillStyle(3001);
MuMuMass3->SetFillColor(5);
MuMuMass4->SetFillColor(38);
MuMuMass5->SetFillColor(903);
sMuMuMass->Add(MuMuMass3);
sMuMuMass->Add(MuMuMass2);
sMuMuMass->Add(MuMuMass1);
sMuMuMass->Add(MuMuMass4);
sMuMuMass->Add(MuMuMass5);
if(sMuMuMass->GetMaximum() > 0) sMuMuMass.SetMaximum(sMuMuMass->GetMaximum() * 1.5);   
else sMuMuMass.SetMaximum(MuMuMass0->GetMaximum() * 2.0);
//sMuMuMass.SetMinimum(0.01);
//sMuMuMass.SetMaximum(50); 
sMuMuMass->Draw();
sMuMuMass->GetXaxis()->SetTitle("#mu#mu mass [GeV]");
sMuMuMass->GetYaxis()->SetTitle("# events / 1 GeV");
MuMuMass0->Draw("same");

Kinematic1->cd(2);
MuMudeta0->Sumw2();
MuMudeta0->SetLineWidth(2);
MuMudeta0->SetMarkerStyle(20);
MuMudeta1->SetFillColor(0);MuMudeta1->SetLineWidth(3);MuMudeta1->SetLineColor(4);
MuMudeta2->SetFillColor(2);
//MuMudeta2->SetFillStyle(3001);
MuMudeta3->SetFillColor(5);
MuMudeta4->SetFillColor(38);
MuMudeta5->SetFillColor(903);
sMuMudeta->Add(MuMudeta3);
sMuMudeta->Add(MuMudeta2);
sMuMudeta->Add(MuMudeta1);
sMuMudeta->Add(MuMudeta4);
sMuMudeta->Add(MuMudeta5);
if(sMuMudeta->GetMaximum() > 0) sMuMudeta.SetMaximum(sMuMudeta->GetMaximum() * 2.0);    
else sMuMudeta.SetMaximum(MuMudeta0->GetMaximum() * 2.0);  
sMuMudeta->Draw();
//sMuMudeta->GetXaxis()->SetTitle("#eta(#mu^{+}) + #eta(#mu^{-})");
sMuMudeta->GetXaxis()->SetTitle("3D opening angle");
sMuMudeta->GetYaxis()->SetTitle("# events / 0.5 ");
MuMudeta0->Draw("same");
}

if(0){
TCanvas *Kinematic2 = new TCanvas("Kinematic2","Kinematic MuMu",800,500);
   Kinematic2->SetFillColor(0);
   Kinematic2->SetBorderMode(0);
   Kinematic2->SetBorderSize(2);
   Kinematic2->SetFrameBorderMode(0);

Kinematic2->Divide(2,1);
Kinematic2->cd(1);
MuMudpt0->Sumw2();
MuMudpt0->SetLineWidth(2);
MuMudpt0->SetMarkerStyle(20);
MuMudpt1->SetFillColor(0);MuMudpt1->SetLineColor(4);MuMudpt1->SetLineWidth(3);
MuMudpt2->SetFillColor(2);
//MuMudpt2->SetFillStyle(3001);
MuMudpt3->SetFillColor(5);
MuMudpt4->SetFillColor(38);
MuMudpt5->SetFillColor(903);
sMuMudpt->Add(MuMudpt3);
sMuMudpt->Add(MuMudpt2);
sMuMudpt->Add(MuMudpt1);
sMuMudpt->Add(MuMudpt4);
sMuMudpt->Add(MuMudpt5);
if(sMuMudpt->GetMaximum() > 0) sMuMudpt.SetMaximum(sMuMudpt->GetMaximum() * 1.5);   
else sMuMudpt.SetMaximum(MuMudpt0->GetMaximum() * 2.0); 
sMuMudpt->Draw();
sMuMudpt->GetXaxis()->SetTitle("#mu#mu |#Delta p_{T}| [GeV]");
sMuMudpt->GetYaxis()->SetTitle("# events / 0.25 GeV ");
MuMudpt0->Draw("same");

Kinematic2->cd(2);
MuMudphi0->Sumw2();
MuMudphi0->SetLineWidth(2);
MuMudphi0->SetMarkerStyle(20);
MuMudphi1->SetFillColor(0);MuMudphi1->SetLineColor(4); MuMudphi1->SetLineWidth(3);
MuMudphi2->SetFillColor(2);
//MuMudphi2->SetFillStyle(3001);
MuMudphi3->SetFillColor(5);
MuMudphi4->SetFillColor(38);
MuMudphi5->SetFillColor(903);
sMuMudphi->Add(MuMudphi3);
sMuMudphi->Add(MuMudphi2);
sMuMudphi->Add(MuMudphi1);
sMuMudphi->Add(MuMudphi4);
sMuMudphi->Add(MuMudphi5);
if(sMuMudphi->GetMaximum() > 0) sMuMudphi.SetMaximum(sMuMudphi->GetMaximum() * 1.5);  
else sMuMudphi.SetMaximum(MuMudphi0->GetMaximum() * 2.0);
sMuMudphi->Draw();
sMuMudphi->GetXaxis()->SetTitle("#mu#mu |#Delta #phi / #pi|");
sMuMudphi->GetYaxis()->SetTitle("# events / 0.02 ");
MuMudphi0->Draw("same");
//Kinematic1->SaveAs("Kinematic1_851nb.eps");
}

if(0){
TCanvas *Kinematic3 = new TCanvas("Kinematic3","Kinematic MuMu",800,500);
   Kinematic3->SetFillColor(0);
   Kinematic3->SetBorderMode(0);
   Kinematic3->SetBorderSize(2);
   Kinematic3->SetFrameBorderMode(0);
   
Kinematic3->Divide(3,1);
Kinematic3->cd(1);
etaPair0->Sumw2();
etaPair0->SetLineWidth(2);
etaPair0->SetMarkerStyle(20);
etaPair1->SetFillColor(0);etaPair1->SetLineColor(4);etaPair1->SetLineWidth(3); 
etaPair2->SetFillColor(2);
//etaPair2->SetFillStyle(3001);
etaPair3->SetFillColor(5);
etaPair4->SetFillColor(38);
etaPair5->SetFillColor(903);
setaPair->Add(etaPair3);
setaPair->Add(etaPair2);
setaPair->Add(etaPair1);
setaPair->Add(etaPair4);
setaPair->Add(etaPair5);
if(setaPair->GetMaximum() > 0) setaPair.SetMaximum(setaPair->GetMaximum() * 2.5);  
else setaPair.SetMaximum(etaPair0->GetMaximum() * 3.0);
setaPair->Draw();
setaPair->GetXaxis()->SetTitle("#mu#mu #eta");
setaPair->GetYaxis()->SetTitle("# events / 0.5");
etaPair0->Draw("same");

Kinematic3->cd(2);
pTPair0->Sumw2();
pTPair0->SetLineWidth(2);
pTPair0->SetMarkerStyle(20);
pTPair1->SetFillColor(0);pTPair1->SetLineColor(4); pTPair1->SetLineWidth(3);
pTPair2->SetFillColor(2);
//pTPair2->SetFillStyle(3001);
pTPair3->SetFillColor(5);
pTPair4->SetFillColor(38);
pTPair5->SetFillColor(903);
spTPair->Add(pTPair3);
spTPair->Add(pTPair2);
spTPair->Add(pTPair1);
spTPair->Add(pTPair4);
spTPair->Add(pTPair5);
if(spTPair->GetMaximum() > 0) spTPair.SetMaximum(spTPair->GetMaximum() * 1.5);   
else spTPair.SetMaximum(spTPair->GetMaximum() * 2.0); 
spTPair->Draw();
spTPair->GetXaxis()->SetTitle("#mu#mu p_{T} [GeV]");
spTPair->GetYaxis()->SetTitle("# events / 0.5 GeV");
pTPair0->Draw("same");

Kinematic3->cd(3);
MuMuSymdphi0->Sumw2(); 
MuMuSymdphi0->SetLineWidth(2); 
MuMuSymdphi0->SetMarkerStyle(20); 
MuMuSymdphi1->SetFillColor(0);MuMuSymdphi1->SetLineWidth(3);MuMuSymdphi1->SetLineColor(4); 
MuMuSymdphi2->SetFillColor(2); 
//MuMuSymdphi2->SetFillStyle(3001); 
MuMuSymdphi3->SetFillColor(5); 
MuMuSymdphi4->SetFillColor(38); 
MuMuSymdphi5->SetFillColor(903); 
sMuMuSymdphi->Add(MuMuSymdphi3); 
sMuMuSymdphi->Add(MuMuSymdphi2); 
sMuMuSymdphi->Add(MuMuSymdphi1); 
sMuMuSymdphi->Add(MuMuSymdphi4); 
sMuMuSymdphi->Add(MuMuSymdphi5); 
if(sMuMuSymdphi->GetMaximum() > 0) sMuMuSymdphi.SetMaximum(sMuMuSymdphi->GetMaximum() * 1.5);   
else sMuMuSymdphi.SetMaximum(MuMuSymdphi0->GetMaximum() * 2.0); 
sMuMuSymdphi->Draw(); 
sMuMuSymdphi->GetXaxis()->SetTitle("#mu#mu 1 - |#Delta #phi / #pi|"); 
sMuMuSymdphi->GetYaxis()->SetTitle("# events / 0.02 "); 
MuMuSymdphi0->Draw("same"); 

}

if(0){
TCanvas *Kinematic4 = new TCanvas("Kinematic4","Kinematic single Muon",800,500);
   Kinematic4->SetFillColor(0);
   Kinematic4->SetBorderMode(0);
   Kinematic4->SetBorderSize(2);
   Kinematic4->SetFrameBorderMode(0);

Kinematic4->Divide(3,1);
Kinematic4->cd(1);
etaSingle0->Sumw2();
etaSingle0->SetLineWidth(2);
etaSingle0->SetMarkerStyle(20);
etaSingle1->SetFillColor(ci);
etaSingle2->SetFillColor(30);
etaSingle2->SetFillStyle(3001);
etaSingle3->SetFillColor(30);
etaSingle4->SetFillColor(38);
etaSingle5->SetFillColor(903);
setaSingle->Add(etaSingle3);
setaSingle->Add(etaSingle2);
setaSingle->Add(etaSingle1);
setaSingle->Add(etaSingle4);
setaSingle->Add(etaSingle5);
setaSingle.SetMaximum(setaSingle->GetMaximum() * 1.5);  
setaSingle->Draw();
setaSingle->GetXaxis()->SetTitle("#mu #eta");
setaSingle->GetYaxis()->SetTitle("# events / 0.5");
etaSingle0->Draw("same");

Kinematic4->cd(2);
phiSingle0->Sumw2();
phiSingle0->SetLineWidth(2);
phiSingle0->SetMarkerStyle(20);
phiSingle1->SetFillColor(ci);
phiSingle2->SetFillColor(30);
phiSingle2->SetFillStyle(3001);
phiSingle3->SetFillColor(30);
phiSingle4->SetFillColor(38);
phiSingle5->SetFillColor(903);
sphiSingle->Add(phiSingle3);
sphiSingle->Add(phiSingle2);
sphiSingle->Add(phiSingle1);
sphiSingle->Add(phiSingle4);
sphiSingle->Add(phiSingle5);
sphiSingle.SetMaximum(sphiSingle->GetMaximum() * 1.5);   
sphiSingle->Draw();
sphiSingle->GetXaxis()->SetTitle("#mu #phi");
sphiSingle->GetYaxis()->SetTitle("# events / 0.5");
phiSingle0->Draw("same");

Kinematic4->cd(3);
pTSingle0->Sumw2();
pTSingle0->SetLineWidth(2);
pTSingle0->SetMarkerStyle(20);
pTSingle1->SetFillColor(ci);
pTSingle2->SetFillColor(30);
pTSingle2->SetFillStyle(3001);
pTSingle3->SetFillColor(30);
pTSingle4->SetFillColor(38);
pTSingle5->SetFillColor(903);
spTSingle->Add(pTSingle3);
spTSingle->Add(pTSingle2);
spTSingle->Add(pTSingle1);
spTSingle->Add(pTSingle4);
spTSingle->Add(pTSingle5);
spTSingle.SetMaximum(spTSingle->GetMaximum() * 1.5);   
spTSingle->Draw();
spTSingle->GetXaxis()->SetTitle("#mu p_{T} [GeV]");
spTSingle->GetYaxis()->SetTitle("# events / 0.5 GeV");
pTSingle0->Draw("same");
}

//Kinematic1->SaveAs("Kinematic1_851nb.eps");

if(0){
TCanvas *CASTOR = new TCanvas("CASTOR","CASTOR",800,550);
   CASTOR->SetFillColor(0);
   CASTOR->SetBorderMode(0);
   CASTOR->SetBorderSize(2);
   CASTOR->SetFrameBorderMode(0);

CASTOR->cd(1);
CastorSumE0->Sumw2();
CastorSumE0->SetLineWidth(2);
CastorSumE0->SetMarkerStyle(20);
CastorSumE1->SetFillColor(ci);
CastorSumE2->SetFillColor(30);
CastorSumE2->SetFillStyle(3001);
CastorSumE3->SetFillColor(30);
CastorSumE4->SetFillColor(38);
CastorSumE5->SetFillColor(903);
sCastorSumE->Add(CastorSumE3);
sCastorSumE->Add(CastorSumE2);
sCastorSumE->Add(CastorSumE1);
sCastorSumE->Add(CastorSumE4);
sCastorSumE->Add(CastorSumE5);
sCastorSumE.SetMaximum(sCastorSumE->GetMaximum() * 1.5);   
sCastorSumE->Draw();
sCastorSumE->GetXaxis()->SetTitle("Castor Sum Energy [GeV]");
sCastorSumE->GetYaxis()->SetTitle("# events / 5 GeV");
CastorSumE0->Draw("same");
}

if(0){
TCanvas *ZDC = new TCanvas("ZDC","ZDC",800,550);
   ZDC->SetFillColor(0);
   ZDC->SetBorderMode(0);
   ZDC->SetBorderSize(2);
   ZDC->SetFrameBorderMode(0);

ZDC->Divide(4,2);
ZDC->cd(1);
ZDCemplus0->Sumw2();
ZDCemplus0->SetLineWidth(2);
ZDCemplus0->SetMarkerStyle(20);
ZDCemplus1->SetFillColor(ci);
ZDCemplus2->SetFillColor(30);
ZDCemplus2->SetFillStyle(3001);
ZDCemplus3->SetFillColor(30);
ZDCemplus4->SetFillColor(38);
ZDCemplus5->SetFillColor(903);
sZDCemplus->Add(ZDCemplus3);
sZDCemplus->Add(ZDCemplus2);
sZDCemplus->Add(ZDCemplus1);
sZDCemplus->Add(ZDCemplus4);
sZDCemplus->Add(ZDCemplus5);
sZDCemplus.SetMaximum(sZDCemplus->GetMaximum() * 1.5);   
sZDCemplus->Draw();
sZDCemplus->GetXaxis()->SetTitle("ZDC + em [GeV]");
sZDCemplus->GetYaxis()->SetTitle("# events / 20 GeV");
ZDCemplus0->Draw("same");

ZDC->cd(2);
ZDCemminus0->Sumw2();
ZDCemminus0->SetLineWidth(2);
ZDCemminus0->SetMarkerStyle(20);
ZDCemminus1->SetFillColor(ci);
ZDCemminus2->SetFillColor(30);
ZDCemminus2->SetFillStyle(3001);
ZDCemminus3->SetFillColor(30);
ZDCemminus4->SetFillColor(38);
ZDCemminus5->SetFillColor(903);
sZDCemminus->Add(ZDCemminus3);
sZDCemminus->Add(ZDCemminus2);
sZDCemminus->Add(ZDCemminus1);
sZDCemminus->Add(ZDCemminus4);
sZDCemminus->Add(ZDCemminus5);
sZDCemminus.SetMaximum(sZDCemminus->GetMaximum() * 1.5);    
sZDCemminus->Draw();
sZDCemminus->GetXaxis()->SetTitle("ZDC - em [GeV]");
sZDCemminus->GetYaxis()->SetTitle("# events / 20 GeV");
ZDCemminus0->Draw("same");

ZDC->cd(3);
ZDChadplus0->Sumw2();
ZDChadplus0->SetLineWidth(2);
ZDChadplus0->SetMarkerStyle(20);
ZDChadplus1->SetFillColor(ci);
ZDChadplus2->SetFillColor(30);
ZDChadplus2->SetFillStyle(3001);
ZDChadplus3->SetFillColor(30);
ZDChadplus4->SetFillColor(38);
ZDChadplus5->SetFillColor(903);
sZDChadplus->Add(ZDChadplus3);
sZDChadplus->Add(ZDChadplus2);
sZDChadplus->Add(ZDChadplus1);
sZDChadplus->Add(ZDChadplus4);
sZDChadplus->Add(ZDChadplus5);
sZDChadplus.SetMaximum(sZDChadplus->GetMaximum() * 1.5);    
sZDChadplus->Draw();
sZDChadplus->GetXaxis()->SetTitle("ZDC + had [GeV]");
sZDChadplus->GetYaxis()->SetTitle("# events / 20 GeV");
ZDChadplus0->Draw("same");

ZDC->cd(4);
ZDChadminus0->Sumw2();
ZDChadminus0->SetLineWidth(2);
ZDChadminus0->SetMarkerStyle(20);
ZDChadminus1->SetFillColor(ci);
ZDChadminus2->SetFillColor(30);
ZDChadminus2->SetFillStyle(3001);
ZDChadminus3->SetFillColor(30);
ZDChadminus4->SetFillColor(38);
ZDChadminus5->SetFillColor(903);
sZDChadminus->Add(ZDChadminus3);
sZDChadminus->Add(ZDChadminus2);
sZDChadminus->Add(ZDChadminus1);
sZDChadminus->Add(ZDChadminus4);
sZDChadminus->Add(ZDChadminus5);
sZDCemminus.SetMaximum(sZDCemminus->GetMaximum() * 1.5);    
sZDChadminus->Draw();
sZDChadminus->GetXaxis()->SetTitle("ZDC - had [GeV]");
sZDChadminus->GetYaxis()->SetTitle("# events / 20 GeV");
ZDChadminus0->Draw("same");

ZDC->cd(5);
ZDCtime0->Sumw2();
ZDCtime0->SetLineWidth(2);
ZDCtime0->SetMarkerStyle(20);
ZDCtime1->SetFillColor(ci);
ZDCtime2->SetFillColor(30);
ZDCtime2->SetFillStyle(3001);
ZDCtime3->SetFillColor(30);
ZDCtime4->SetFillColor(38);
ZDCtime5->SetFillColor(903);
sZDCtime->Add(ZDCtime3);
sZDCtime->Add(ZDCtime2);
sZDCtime->Add(ZDCtime1);
sZDCtime->Add(ZDCtime4);
sZDCtime->Add(ZDCtime5);
sZDCtime.SetMaximum(sZDCtime->GetMaximum() * 1.5);    
sZDCtime->Draw();
sZDCtime->GetXaxis()->SetTitle("ZDC hit time [???]");
sZDCtime->GetYaxis()->SetTitle("# events / 1 ???");
ZDCtime0->Draw("same");

ZDC->cd(7);
ZDCenergyEM0->Sumw2();
ZDCenergyEM0->SetLineWidth(2);
ZDCenergyEM0->SetMarkerStyle(20);
ZDCenergyEM1->SetFillColor(ci);
ZDCenergyEM2->SetFillColor(30);
ZDCenergyEM2->SetFillStyle(3001);
ZDCenergyEM3->SetFillColor(30);
ZDCenergyEM4->SetFillColor(38);
ZDCenergyEM5->SetFillColor(903);
sZDCenergyEM->Add(ZDCenergyEM3);
sZDCenergyEM->Add(ZDCenergyEM2);
sZDCenergyEM->Add(ZDCenergyEM1);
sZDCenergyEM->Add(ZDCenergyEM4);
sZDCenergyEM->Add(ZDCenergyEM5);
sZDCenergyEM.SetMaximum(sZDCenergyEM->GetMaximum() * 1.5);    
sZDCenergyEM->Draw();
sZDCenergyEM->GetXaxis()->SetTitle("ZDC hit EM energy [GeV ??]");
sZDCenergyEM->GetYaxis()->SetTitle("# events / 30 GeV ???");
ZDCenergyEM0->Draw("same");

ZDC->cd(8);
ZDCenergyHAD0->Sumw2();
ZDCenergyHAD0->SetLineWidth(2);
ZDCenergyHAD0->SetMarkerStyle(20);
ZDCenergyHAD1->SetFillColor(ci);
ZDCenergyHAD2->SetFillColor(30);
ZDCenergyHAD2->SetFillStyle(3001);
ZDCenergyHAD3->SetFillColor(30);
ZDCenergyHAD4->SetFillColor(38);
ZDCenergyHAD5->SetFillColor(903);
sZDCenergyHAD->Add(ZDCenergyHAD3);
sZDCenergyHAD->Add(ZDCenergyHAD2);
sZDCenergyHAD->Add(ZDCenergyHAD1);
sZDCenergyHAD->Add(ZDCenergyHAD4);
sZDCenergyHAD->Add(ZDCenergyHAD5);
sZDCenergyHAD.SetMaximum(sZDCenergyHAD->GetMaximum() * 1.5);    
sZDCenergyHAD->Draw();
sZDCenergyHAD->GetXaxis()->SetTitle("ZDC hit HAD energy [GeV ??]");
sZDCenergyHAD->GetYaxis()->SetTitle("# events / 300 GeV ???");
ZDCenergyHAD0->Draw("same");
}

if(0){
TCanvas *Tdist = new TCanvas("Tdist","Tdist",800,550); 
   Tdist->SetFillColor(0); 
   Tdist->SetBorderMode(0); 
   Tdist->SetBorderSize(2); 
   Tdist->SetFrameBorderMode(0); 
 
Tdist->cd(1); 
Tdist0->Sumw2(); 
Tdist0->SetLineWidth(2); 
Tdist0->SetMarkerStyle(20); 
Tdist1->SetFillColor(ci); 
Tdist2->SetFillColor(30); 
Tdist2->SetFillStyle(3001); 
Tdist3->SetFillColor(30); 
Tdist4->SetFillColor(38); 
Tdist5->SetFillColor(903); 
//Tdist5->SetMaximum(Tdist0->GetMaximum() * 1.5);
sTdist->Add(Tdist3);
sTdist->Add(Tdist2);
sTdist->Add(Tdist1);
sTdist->Add(Tdist4);
sTdist->Add(Tdist5);
sTdist->Draw();
sTdist->GetXaxis()->SetTitle("#mu#mu p_{T}^{2} [GeV]"); 
sTdist->GetYaxis()->SetTitle("# events / 5 GeV"); 
Tdist0->Draw("same"); 
//Tdist->SaveAs("Tdist_851nb.eps");
}

	cout << "END" << endl;   
}

