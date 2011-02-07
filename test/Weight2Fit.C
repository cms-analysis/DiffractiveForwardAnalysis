void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

// For the Global title:

//  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

// personnal additions
  tdrStyle->SetPalette(1);
  tdrStyle->cd();
}



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
	TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas MuMu",600,600); 
   	Canvas1->SetFillColor(0); 
   	Canvas1->SetBorderMode(0); 
   	Canvas1->SetBorderSize(2); 
   	Canvas1->SetFrameBorderMode(0);

  Canvas1->SetDefH(600); //Height of canvas
  Canvas1->SetDefW(600); //Width of canvas
  Canvas1->SetDefX(0);   //POsition on screen
  Canvas1->SetDefY(0);


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
	int ntowerthresh = 9999;//5;
	if((nEB+nEE+nHB+nHE+nHFp+nHFm) < ntowerthresh)
		pass = true;

	return pass;

}
 
bool FailsTrackDistanceVeto(double trackdistance, int trackquality, int number_hits)
{
	bool fail = false;
	float vtxtrackcountingcut = 0.2; 

	if(trackdistance < vtxtrackcountingcut)
		fail = true;

//	if((trackdistance < vtxtrackcountingcut) && (trackquality == 1) && (number_hits > 10))
//		fail = true;

	return fail;
}

double VertexSeparation(int nvtx, int* vtxtrks, double* vtxz, int* ismumuvtx)
{
	double closestvtx = 9999.;
	if(nvtx == 1)
		return 9999.;
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
  	float ZDChadThresh =  9999120.; //120.0;
  	float ZDCemThresh =   9999918.; // 16.0;

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
  	float lowermasscut1 = 10.0;
  	float uppermasscut1 = 8.5;
  	float lowermasscut2 = 11.5;
  	float uppermasscut2 = 25.0; 

        if((mass < lowermasscut1) || (mass > uppermasscut2)) 
		pass = false;
        if((mass > uppermasscut1) && (mass < lowermasscut2))  
		pass = false;

	if(ptPlus<4.0 || ptMinus<4.0)
		pass = false;

        if(fabs(etaPlus)>2.1 || fabs(etaMinus)>2.1)
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

Double_t selectHist(Double_t *x, Double_t *par){
//	TFile *fMC = new TFile("MonteCarloHisto.root","READ"); // 
	Double_t xx=x[0];
//	cout<<"check param "<<par[0]<<" "<<par[1]<<" "<<par[2]<<"  for x="<<xx<<endl;
//-------------
	histo_inel=(TH1F*)fit_Inel;
	Int_t bin = histo_inel->GetXaxis()->FindBin(xx);
	Double_t pt_pair=histo_inel->GetBinContent(bin);
	Double_t dumping=TMath::Exp(-(par[0])*pt_pair*pt_pair);
	Double_t value=pt_pair*par[2]*par[1]*dumping;
//	cout<<" -->add "<<pt_pair<<" x "<<par[1]<<"(n) x "<<dumping<<"(d) x "<<par[2]<<"(l)"<<endl;
//----------------------
	histo_others=(TH1F*)fit_others;
	bin = histo_others->GetXaxis()->FindBin(xx);
	value += histo_others->GetBinContent(bin)*par[2];
//	cout<<" -->add "<<histo_others->GetBinContent(bin)<<" x "<<par[2]<<"(l)"<<endl;
//	cout<<" ==>f(x)="<<value<<endl;
	return value;
}

void Weight2Fit()
{
setTDRStyle();
/*gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gROOT->SetStyle("Plain");
gStyle->SetStatStyle(0);
gROOT->SetTitle(0);*/
#define pi 3.14159265359

//definition des fichiers + Tree
  TFile *f0 = new TFile("../cand_2tracks.root"); // 
//  TFile *f0 = new TFile("cand_2tracks.root"); // 
  TTree *t0 = f0->Get("ntp1");
  TFile *f1 = new TFile("/home/fynu/schul/scratch/data_analyses/TagAndProbe/CMSSW_3_8_5/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/El-El.root"); //
  TTree *t1 = f1->Get("ntp1");
  TFile *f2 = new TFile("/home/fynu/schul/scratch/data_analyses/TagAndProbe/CMSSW_3_8_5/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/Inel-El.root"); //
  TTree *t2 = f2->Get("ntp1");
  TFile *f3 = new TFile("/home/fynu/schul/scratch/data_analyses/TagAndProbe/CMSSW_3_8_5/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/Inel-Inel.root"); //
  TTree *t3 = f3->Get("ntp1");
  TFile *f4 = new TFile("../Upsilon.root"); //
  TTree *t4 = f4->Get("ntp1");
  TFile *f5 = new TFile("../Jpsi.root"); //
  TTree *t5 = f5->Get("ntp1");

  TFile *f6 = new TFile("/home/fynu/schul/scratch/data_analyses/TagAndProbe/CMSSW_3_8_5/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/DY_tp385.root"); //
//  TFile *f6 = new TFile("/storage/data/cms/store/user/schul/384_novQCD/QCD2MU_merge.root");
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
  TH1F* MuMuMass0 = new TH1F("mass_data","",81,-1.,80.);
  TH1F* MuMuMass1 = new TH1F("mass_cumulElEl","",81,-1.,80.);
  TH1F* MuMuMass2 = new TH1F("mass_cumulInelEl","",81,-1.,80.);
  TH1F* MuMuMass3 = new TH1F("mass_cumulInelInel","",81,-1.,80.);
  TH1F* MuMuMass4 = new TH1F("mass_cumulUps","",81,-1.,80.);
  TH1F* MuMuMass5 = new TH1F("mass_cumulJpsi","",81,-1.,80.);
  TH1F* MuMuMass6 = new TH1F("mass_cumulInclu","",81,-1.,80.);
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
  TH1F* MuMudpt0 = new TH1F("dpt_data","",60,-0.5,5.5);
  TH1F* MuMudpt1 = new TH1F("dpt_cumulElEl","",60,-0.5,5.5);
  TH1F* MuMudpt2 = new TH1F("dpt_cumulInelEl","",60,-0.5,5.5);
  TH1F* MuMudpt3 = new TH1F("dpt_cumulInelInel","",60,-0.5,5.5);
  TH1F* MuMudpt4 = new TH1F("dpt_cumulUps","",60,-0.5,5.5);
  TH1F* MuMudpt5 = new TH1F("dpt_cumulJpsi","",60,-0.5,5.5);
  TH1F* MuMudpt6 = new TH1F("dpt_cumulInclu","",60,-0.5,5.5);
  THStack *sMuMudpt = new THStack("sMuMudpt","stack dpt");

  // 60,-0.1,1.1
  TH1F* MuMudphi0 = new TH1F("dphi_data","",44,-0.02,0.42);
  TH1F* MuMudphi1 = new TH1F("dphi_cumulElEl","",44,-0.02,0.42);
  TH1F* MuMudphi2 = new TH1F("dphi_cumulInelEl","",44,-0.02,0.42);
  TH1F* MuMudphi3 = new TH1F("dphi_cumulInelInel","",44,-0.02,0.42);
  TH1F* MuMudphi4 = new TH1F("dphi_cumulUps","",44,-0.02,0.42);
  TH1F* MuMudphi5 = new TH1F("dphi_cumulJpsi","",44,-0.02,0.42);
  TH1F* MuMudphi6 = new TH1F("dphi_cumulInclu","",44,-0.02,0.42);
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
  TH1F* MuMu3DAng0 = new TH1F("a3DAng_data","",12,-0.1,1.1);
  TH1F* MuMu3DAng1 = new TH1F("a3DAng_cumulElEl","",12,-0.1,1.1);
  TH1F* MuMu3DAng2 = new TH1F("a3DAng_cumulInelEl","",12,-0.1,1.1);
  TH1F* MuMu3DAng3 = new TH1F("a3DAng_cumulInelInel","",12,-0.1,1.1);
  TH1F* MuMu3DAng4 = new TH1F("a3DAng_cumulUps","",12,-0.1,1.1);
  TH1F* MuMu3DAng5 = new TH1F("a3DAng_cumulJpsi","",12,-0.1,1.1);
  TH1F* MuMu3DAng6 = new TH1F("a3DAng_cumulInclu","",12,-0.1,1.1);
  THStack *sMuMu3DAng = new THStack("sMuMu3DAng","stack 3DAng");


  TH1F* MuMuDeta0 = new TH1F("aDeta_data","",27,-0.2,5.2);
  TH1F* MuMuDeta1 = new TH1F("aDeta_cumulElEl","",27,-0.2,5.2);
  TH1F* MuMuDeta2 = new TH1F("aDeta_cumulInelEl","",27,-0.2,5.2);
  TH1F* MuMuDeta3 = new TH1F("aDeta_cumulInelInel","",27,-0.2,5.2);
  TH1F* MuMuDeta4 = new TH1F("aDeta_cumulUps","",27,-0.2,5.2);
  TH1F* MuMuDeta5 = new TH1F("aDeta_cumulJpsi","",27,-0.2,5.2);
  TH1F* MuMuDeta6 = new TH1F("aDeta_cumulInclu","",27,-0.2,5.2);
  THStack *sMuMuDeta = new THStack("sMuMuDeta","stack Deta");


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

  TH1F* etaPair0 = new TH1F("etaPair_data","",12,-3.,3.);
  TH1F* etaPair1 = new TH1F("etaPair_cumulElEl","",12,-3.,3.);
  TH1F* etaPair2 = new TH1F("etaPair_cumulInelEl","",12,-3.,3.);
  TH1F* etaPair3 = new TH1F("etaPair_cumulInelInel","",12,-3.,3.);
  TH1F* etaPair4 = new TH1F("etaPair_cumulUps","",12,-3.,3.);
  TH1F* etaPair5 = new TH1F("etaPair_cumulJpsi","",12,-3.,3.);
  TH1F* etaPair6 = new TH1F("etaPair_cumulInclu","",12,-3.,3.);
  THStack *setaPair = new THStack("setaPair","stack eta Pair");

  TH1F* pTPair0 = new TH1F("pTPair_data","",20,0.,3.);
  TH1F* pTPair1 = new TH1F("pTPair_cumulElEl","",20,0.,3.);
  TH1F* pTPair2 = new TH1F("pTPair_cumulInelEl","",20,0.,3.);
  TH1F* pTPair3 = new TH1F("pTPair_cumulInelInel","",20,0.,3.);
  TH1F* pTPair4 = new TH1F("pTPair_cumulUps","",20,0.,3.);
  TH1F* pTPair5 = new TH1F("pTPair_cumulJpsi","",20,0.,3.);
  TH1F* pTPair6 = new TH1F("pTPair_cumulInclu","",20,0.,3.);
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


  TH1F* etaSingleP0 = new TH1F("etaSingleP_data","",12,-3.,3.);
  TH1F* etaSingleP1 = new TH1F("etaSingleP_cumulElEl","",12,-3.,3.);
  TH1F* etaSingleP2 = new TH1F("etaSingleP_cumulInelEl","",12,-3.,3.);
  TH1F* etaSingleP3 = new TH1F("etaSingleP_cumulInelInel","",12,-3.,3.);
  TH1F* etaSingleP4 = new TH1F("etaSingleP_cumulUps","",12,-3.,3.);
  TH1F* etaSingleP5 = new TH1F("etaSingleP_cumulJpsi","",12,-3.,3.);
  TH1F* etaSingleP6 = new TH1F("etaSingleP_cumulInclu","",12,-3.,3.);
  THStack *setaSingleP = new THStack("setaSingleP","stack eta SingleP");

  TH1F* etaSingleM0 = new TH1F("etaSingleM_data","",12,-3.,3.);
  TH1F* etaSingleM1 = new TH1F("etaSingleM_cumulElEl","",12,-3.,3.);
  TH1F* etaSingleM2 = new TH1F("etaSingleM_cumulInelEl","",12,-3.,3.);
  TH1F* etaSingleM3 = new TH1F("etaSingleM_cumulInelInel","",12,-3.,3.);
  TH1F* etaSingleM4 = new TH1F("etaSingleM_cumulUps","",12,-3.,3.);
  TH1F* etaSingleM5 = new TH1F("etaSingleM_cumulJpsi","",12,-3.,3.);
  TH1F* etaSingleM6 = new TH1F("etaSingleM_cumulInclu","",12,-3.,3.);
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

  TH1F* VtxT0 = new TH1F("vtx_data","",50,0.,0.2);
  TH1F* VtxT1 = new TH1F("vtx_cumulElEl","",50,0.,0.2);
  TH1F* VtxT2 = new TH1F("vtx_cumulInelEl","",50,0.,0.2);
  TH1F* VtxT3 = new TH1F("vtx_cumulInelInel","",50,0.,0.2);
  TH1F* VtxT4 = new TH1F("vtx_cumulUps","",50,0.,0.2);
  TH1F* VtxT5 = new TH1F("vtx_cumulJpsi","",50,0.,0.2);
  TH1F* VtxT6 = new TH1F("vtx_cumulInclu","",50,0.,0.2);
  THStack* sVtxT = new THStack("sVtxT","stack vtx Transverse");

  TH1F* VtxZ0 = new TH1F("vtxZ_data","",20,-25.,25.);
  TH1F* VtxZ1 = new TH1F("vtxZ_cumulElEl","",20,-25.,25.);
  TH1F* VtxZ2 = new TH1F("vtxZ_cumulInelEl","",20,-25.,25.);
  TH1F* VtxZ3 = new TH1F("vtxZ_cumulInelInel","",20,-25.,25.);
  TH1F* VtxZ4 = new TH1F("vtxZ_cumulUps","",20,-25.,25.);
  TH1F* VtxZ5 = new TH1F("vtxZ_cumulJpsi","",20,-25.,25.);
  TH1F* VtxZ6 = new TH1F("vtxZ_cumulInclu","",20,-25.,25.);
  THStack* sVtxZ = new THStack("sVtxZ","stack vtx Z");


/*
  TH2F* correl0 = new TH2F("correl_data","",80,0.,1.,48,0.,12.); //80,-0.5,3.5,68,-0.02,0.32
  TH2F* correl1 = new TH2F("correl_signal","",80,0.,1.,48,0.,12.);
  TH2F* correl2 = new TH2F("correl_inelel","",80,0.,1.,48,0.,12.);
  TH2F* correl3 = new TH2F("correl_inelinel","",80,0.,1.,48,0.,12.);
  TH2F* correl6 = new TH2F("correl_drell","",80,0.,1.,48,0.,12.);
*/
  TH1F* correl0 = new TH1F("correl_data","",20,0.,3.);
  TH1F* correl1 = new TH1F("correl_1","",20,0.,3.);
  TH1F* correl2 = new TH1F("correl_2","",20,0.,3.);
  TH1F* correl3 = new TH1F("correl_3","",20,0.,3.);
  TH1F* correl6 = new TH1F("correl_6","",20,0.,3.);

/*
  TH2F* hdata = new TH2F("hdata","",80,0.,1.,48,0.,12.);
  TH2F* hEl = new TH2F("hEl","",80,0.,1.,48,0.,12.);
*/
  TH1F* hdata = new TH1F("hdata","",20,0.,3.);
  TH1F* hEl = new TH1F("hEl","",20,0.,3.);


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
  const float integrated_lumi = 35437.511542*0.892;
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
  const float fac_lumi2 = 3.0543555e-7*integrated_lumi;                    //  3.05250e-6*integrated_lumi;
  const float fac_lumi3 = 3.6488765e-7*integrated_lumi;                    //  4.740e-6*integrated_lumi;
  const float fac_lumi4 = 1.350e-6*integrated_lumi;			//  1.350e-6*integrated_lumi;
  const float fac_lumi5 = 3.02430e-4*integrated_lumi;

  const float fac_lumiBkg[21]={0.,1.,3.4204003e-5*integrated_lumi,3.,4.,5.,6.,7.,8.,1.8058285e-3*integrated_lumi,1.1935327e-6*integrated_lumi,11.,12.,13.,14.,15.,16.,17.,18.,19.,5.6770716e-7*integrated_lumi};
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
  Int_t var_nHits0[2000], var_nHits1[2000],var_nHits2[2000],var_nHits3[2000],var_nHits4[2000],var_nHits5[2000],var_nHits6[2000];

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

  Double_t var_MassYp2[1], var_MassYm2[1], var_MassYp3[1], var_MassYm3[1];

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

  t0->SetBranchAddress("TrackCand_nhits",var_nHits0);
  t1->SetBranchAddress("TrackCand_nhits",var_nHits1);
  t2->SetBranchAddress("TrackCand_nhits",var_nHits2);
  t3->SetBranchAddress("TrackCand_nhits",var_nHits3);
  t4->SetBranchAddress("TrackCand_nhits",var_nHits4);
  t5->SetBranchAddress("TrackCand_nhits",var_nHits5);
  t6->SetBranchAddress("TrackCand_nhits",var_nHits6);

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

  t2->SetBranchAddress("MassYPlus",var_MassYp2);
  t2->SetBranchAddress("MassYMinus",var_MassYm2);
  t3->SetBranchAddress("MassYPlus",var_MassYp3);
  t3->SetBranchAddress("MassYMinus",var_MassYm3);

  int filter0Gen(0);
  int filter0Track(0);
  int filter0Events(0);
  for(Int_t i = 0;i < NUM0;i++){
      	t0->GetEntry(i);
	if(var_charge0[var_Pair0[0]]>0){ int pair1 = var_Pair0[0]; int pair2 = var_Pair0[1];}
	else{int pair1 = var_Pair0[1]; int pair2 = var_Pair0[0];}
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

        double pt_pair =sqrt((var_px0[pair1]+var_px0[pair2])*(var_px0[pair1]+var_px0[pair2])+(var_py0[pair1]+var_py0[pair2])*(var_py0[pair1]+var_py0[pair2]));

	if(label_vertex!=99
	   && sqrt(pow(var_vtxX0[label_vertex],2)+pow(var_vtxY0[label_vertex],2))< 0.15
	   && sqrt(pow(var_vtxX0[label_vertex],2)+pow(var_vtxY0[label_vertex],2))> 0.05
           && PassesMuonID(var_tracker0[pair1], muAng1, var_global0[pair1], var_tracker0[pair2], muAng2, var_global0[pair2], var_nhitsTrack0[pair1], var_nhitsTrack0[pair2])
 	   && PassesDptCut(var_dpt0[0]) 
	   && PassesDphiCut(var_dphi0[0]/pi)
//	   && pt_pair<0.8
		) {
	    int nTrackExclu(0);
            for(Int_t j=0; j<nTrack; j++){  
		if(FailsTrackDistanceVeto(var_TrackD0[j], var_TrackQuality0[j], var_nHits0[j])){
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

	if(nTrackExclu<1 
		&&  PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm)
		&&  PassesZDCVeto(var_zdcEmMinus0[0],var_zdcEmPlus0[0],var_zdcHadMinus0[0],var_zdcHadPlus0[0]))
	{ 
          filter0Events++;

          MuMudpt0->Fill(var_dpt0[0],fac_lumi0);
          MuMudphi0->Fill(1-(var_dphi0[0]/pi),fac_lumi0);
	  pTPair0->Fill(pt_pair/**pt_pair*/,fac_lumi0);
	  hdata->Fill(pt_pair/**pt_pair*/);
          correl0->Fill(pt_pair/**pt_pair*/);
		if(pt_pair>2){
			cout<<"----------"<<endl;
			cout<<"pt pair="<<pt_pair<<endl;
			cout<<"m="<<var_mass0[0]<<"pt="<<var_pt0[pair1]<<" & "<<var_pt0[pair2]<<", eta="<<var_eta0[pair1]<<" & "<<var_eta0[pair2]<<endl;
			cout<<"dpt="<<var_dpt0[0]<<", dphi/pi="<<var_dphi0[0]/pi<<endl;
		}
	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"Data :"<<endl;
cout<<"  # Dimuon events = "<<filter0Events<<endl;

 correl0->Sumw2();

 TFile *fMC = new TFile("MonteCarloHisto.root","READ"); // 
 fit_Inel=(TH1F*)fMC->Get("fit_Inel");
 fit_others=(TH1F*)fMC->Get("fit_others");

 TF1 *ffit = new TF1("ffit",selectHist,0.,3.,3); // xmin=0, xmax=3, param=3: par[0/1/2]=a/norm/lumi 
 ffit->SetParameter(0,0.047);
 ffit->SetParameter(1,1.0);
 ffit->SetParameter(2,0.75);
 ffit->SetParLimits(0,-0.4,0.4);
 ffit->SetParLimits(1,0.1,3.5);
 ffit->SetParLimits(2,0.1,2.0);
 correl0->Fit("ffit","LMV");
cout<<"-----------"<<endl;
 correl0->Fit("ffit","LM");
cout<<"------I-----"<<endl;
 correl0->Fit("ffit","LM");
cout<<"------IE-----"<<endl;
 correl0->Fit("ffit","LEM");
cout<<"------IEM-----"<<endl;
 correl0->Fit("ffit","LME");
cout<<endl;
cout<<endl;
cout<<endl;
  cout << "chi^2 / nDF = " << ffit->GetChisquare() <<  "\t" << ffit->GetNDF()
       << "\t = " << ffit->GetChisquare()/ffit->GetNDF() << "\t ["
       << 1-sqrt(2)/sqrt(ffit->GetNDF()) << ";" << 1+sqrt(2)/sqrt(ffit->GetNDF()) << "]\n";
  cout << endl;

  float g2 = ffit->GetParameter(2), g2e = ffit->GetParError(2);
  float g1 = ffit->GetParameter(1), g1e = ffit->GetParError(1);
  float g0 = ffit->GetParameter(0), g0e = ffit->GetParError(0);

  cout<<"best-fit values: a="<<g0<<"+/-"<<g0e<<endl;
  cout<<"                 n="<<g1<<"+/-"<<g1e<<endl;
  cout<<"                 L="<<g2<<"+/-"<<g2e<<endl;
  cout << endl;

  cout << endl;
  matrix = new double[9];
  gMinuit->mnemat (matrix,3);

 cout<<"Correlation matrix:"<<endl;
 cout<<"("<<matrix[0]<<"\t"<<matrix[1]<<"\t"<<matrix[2]<<")"<<endl;
 cout<<"("<<matrix[3]<<"\t"<<matrix[4]<<"\t"<<matrix[5]<<")"<<endl;
 cout<<"("<<matrix[6]<<"\t"<<matrix[7]<<"\t"<<matrix[8]<<")"<<endl;
  cout << endl;

   TVirtualFitter *fitter = TVirtualFitter::Fitter(ffit);
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx;
   fitter->GetStats(amin,edm,errdef,nvpar,nparx);
   printf("amin=%g, chi2=%g\n", amin, ffit->GetChisquare());

/*
   TCanvas *c2 = new TCanvas("c2","contours",10,10,800,600);
   c2->Divide(3,2);
   c2->cd(1);
   //get first contour for parameter 0 versus parameter 1
   gMinuit->SetErrorDef(1);
   TGraph *gr01 = (TGraph*)gMinuit->Contour(100,0,1);
   gr01->Draw("alp");
   c2->cd(2);
   TGraph *gr02 = (TGraph*)gMinuit->Contour(100,0,2);
   gr02->Draw("alp");
   c2->cd(3);
   TGraph *gr12 = (TGraph*)gMinuit->Contour(100,1,2);
   gr12->Draw("alp");
   // 1sigma and 2 sigma contours:
   c2->cd(4);
   gMinuit->SetErrorDef(4); //note 4 and not 2!
   TGraph *gr01_2 = (TGraph*)gMinuit->Contour(100,0,1);
   gr01_2->SetFillColor(42);
   gr01_2->Draw("alf");
   gMinuit->SetErrorDef(1);
   TGraph *gr01_1 = (TGraph*)gMinuit->Contour(100,0,1);
   gr01_1->SetFillColor(38);
   gr01_1->Draw("lf");
   // 1sigma and 2 sigma contours:
   c2->cd(5);
   gMinuit->SetErrorDef(4); //note 4 and not 2!
   TGraph *gr02_2 = (TGraph*)gMinuit->Contour(100,0,2);
   gr02_2->SetFillColor(42);
   gr02_2->Draw("alf");
   gMinuit->SetErrorDef(1);
   TGraph *gr02_1 = (TGraph*)gMinuit->Contour(100,0,2);
   gr02_1->SetFillColor(38);
   gr02_1->Draw("lf");
   // 1sigma and 2 sigma contours:
   c2->cd(6);
   gMinuit->SetErrorDef(4); //note 4 and not 2!
   TGraph *gr12_2 = (TGraph*)gMinuit->Contour(100,1,2);
   gr12_2->SetFillColor(42);
   gr12_2->Draw("alf");
   gMinuit->SetErrorDef(1);
   TGraph *gr12_1 = (TGraph*)gMinuit->Contour(100,1,2);
   gr12_1->SetFillColor(38);
   gr12_1->Draw("lf");
*/


// Save into histo
//Draw
   TCanvas *csummary = new TCanvas("csummary","Summary");
   csummary->SetFillColor(0);
   csummary->SetBorderMode(0);
   csummary->SetBorderSize(2);
   csummary->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

  csummary->cd();
  correl0->Draw();
  ffit->SetNpx(2000);
  ffit->Draw("same");  
  csummary->SaveAs("Test.C");
	cout << "END" << endl;   
}

