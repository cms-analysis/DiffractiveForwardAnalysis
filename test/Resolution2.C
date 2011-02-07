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
  	float dphicut = 0.000009;

	if(deltaphioverpi > dphicut)
		pass = true;
	return pass;
}

bool PassesDptCut(double deltapt)
{
	bool pass = false;
  	float dptcut = 100000.0;
	
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

void Resolution2()
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
  TFile *f1 = new TFile("/home/fynu/schul/scratch/data_analyses/TagAndProbe/CMSSW_3_8_5/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/forReso.root"); //
  TTree *t1 = f1->Get("ntp1");

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

  TH1F* pTResolution_rel[20];
  TH1F* pTResolution_abs[20];
  char ntr[15];
  char nta[15];
  for(int k=0; k<20; k++){
    sprintf(ntr,"rel_%d",k);
    sprintf(nta,"abs_%d",k);    
	 pTResolution_rel[k]=new TH1F(ntr,"",100,0.,2.);
         pTResolution_abs[k]=new TH1F(nta,"",100,-1.,1.);
  }
  
cout<<"pass ini"<<endl;
// definitions des # d'entrée
  const int NUM1 = t1->GetEntries();
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

//GEn Muons
  Int_t nGenMuon[1];
  Double_t GenMuon_eta[10], GenMuon_py[10],GenMuon_px[10],GenMuon_pz[10],GenMuon_phi[10],GenMuon_charge[10];
  Double_t GenMuMu_pt[1];

  TString hlttrigger0 = "HLT_DoubleMu0";
  TString hlttrigger = "HLT_L1DoubleMuOpen";	
  TString hlttrigger2 = "HLT_L1DoubleMuOpen_Tight";
  TString hlttrigger3 = "HLT_DoubleMu3";

  t1->SetBranchAddress(hlttrigger3,hlt_d1);
    t1->SetBranchAddress("L1TechnicalTriggers",techBit1);
  t1->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA1);
  t1->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB1);
  t1->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC1);
  t1->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD1);
  t1->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE1);
    t1->SetBranchAddress("nTrackCand",var_nTrack1);
  t1->SetBranchAddress("nQualityTrackCand",var_nTrackQual1);
  t1->SetBranchAddress("TrackCand_pt",var_TrackPt1);
  t1->SetBranchAddress("TrackCand_vtxdxyz",var_TrackD1);
  t1->SetBranchAddress("TrackCand_purity",var_TrackQuality1);
  t1->SetBranchAddress("nPrimVertexCand",var_nvtx1);
  t1->SetBranchAddress("PrimVertexCand_z",var_vtxZ1);
  t1->SetBranchAddress("PrimVertexCand_x",var_vtxX1);
  t1->SetBranchAddress("PrimVertexCand_y",var_vtxY1);
  t1->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2_1);
  t1->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf1);
  t1->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack1);
  t1->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxmumu1);  
  t1->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX1);
  t1->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY1);
  t1->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ1);
  t1->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid1);
    t1->SetBranchAddress("nCaloCand",var_ncalo1);
  t1->SetBranchAddress("CaloTower_ID",var_caloId1);
  t1->SetBranchAddress("CaloTower_e",var_caloEn1);
  t1->SetBranchAddress("CaloTower_t",var_caloTime1);
  t1->SetBranchAddress("CaloTower_dr",var_calodR1);
  t1->SetBranchAddress("CaloTower_eta",var_caloEta1);
  t1->SetBranchAddress("CaloTower_phi",var_caloPhi1);
  t1->SetBranchAddress("Etmiss",var_etmiss1);
  t1->SetBranchAddress("nExtraCaloTowersE5",var_tower1);
    t1->SetBranchAddress("MuMu_mass",var_mass1);
  t1->SetBranchAddress("MuMu_dpt",var_dpt1);
  t1->SetBranchAddress("MuMu_dphi",var_dphi1);
    t1->SetBranchAddress("MuonCand_isglobal",var_global1);
  t1->SetBranchAddress("MuonCand_istracker",var_tracker1);
  t1->SetBranchAddress("MuonCand_isstandalone",var_standalone1);
  t1->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack1);
  t1->SetBranchAddress("MuonCand_charge",var_charge1);
  t1->SetBranchAddress("MuonCand_pt",var_pt1);
  t1->SetBranchAddress("MuonCand_pz",var_pz1);
  t1->SetBranchAddress("MuonCand_phi",var_phi1);
  t1->SetBranchAddress("MuonCand_eta",var_eta1);
  t1->SetBranchAddress("MuonCand_p",var_p1);
  t1->SetBranchAddress("MuonCand_px",var_px1);
  t1->SetBranchAddress("MuonCand_py",var_py1);
  t1->SetBranchAddress("MuonCand_pz",var_pz1);
  t1->SetBranchAddress("MuonPairCand",var_Pair1);
  t1->SetBranchAddress("MuonCand_efficiency",var_eff1);
  t1->SetBranchAddress("nZDChitCand",var_nZDC1);
  t1->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus1);
  t1->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus1);
  t1->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus1);
  t1->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus1);
  t1->SetBranchAddress("ZDChit_time",var_zdcTime1);
  t1->SetBranchAddress("ZDChit_energy",var_zdcE1);
  t1->SetBranchAddress("ZDChit_section",var_zdcsection1);
  t1->SetBranchAddress("nCastorTowerCand",var_nCastor1);
  t1->SetBranchAddress("CastorTower_e",var_CastorE1);
  t1->SetBranchAddress("CastorTower_eta",var_CastorEta1);
  t1->SetBranchAddress("CastorTower_phi",var_CastorPhi1);
  t1->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit1);
  t1->SetBranchAddress("L1TechnicalTriggers",techBit1);

  t1->SetBranchAddress("nGenMuonCand",nGenMuon);
  t1->SetBranchAddress("GenMuonCand_px",GenMuon_px);
  t1->SetBranchAddress("GenMuonCand_py",GenMuon_py);
  t1->SetBranchAddress("GenMuonCand_pz",GenMuon_pz);
  t1->SetBranchAddress("GenMuonCand_phi",GenMuon_phi);
  t1->SetBranchAddress("GenMuonCand_eta",GenMuon_eta);
  t1->SetBranchAddress("GenMuonCand_charge",GenMuon_charge);
  t1->SetBranchAddress("GenMuMu_pt",GenMuMu_pt);

  const float HBThresh = 1.25;
  const float HEThresh = 1.80;
  const float EEThresh = 2.40;
  const float EBThresh = 0.6;
  const float HFPlusThresh = 4.5;
  const float HFMinusThresh = 4.0;
  const float ZDChadThresh = 120.0;
  const float ZDCemThresh = 16.0;
  const float dRcone = 0.3;

  int filter1Gen(0);
 int filter1Track(0);
  double filter1Events(0.);
  for(Int_t i = 0;i < NUM1;i++){
      	t1->GetEntry(i);
//--weight--
        double weight=1.0;
        if(i<=6587) weight=3.5;
//----------
        if(var_charge1[var_Pair1[0]]>0){ int pair1 = var_Pair1[0]; int pair2 = var_Pair1[1];}
        else{int pair1 = var_Pair1[1]; int pair2 = var_Pair1[0];}
        int muID1 = var_idA1[pair1];
        int muID2 = var_idA1[pair2];
	int muAng1 = var_idB1[pair1];
        int muAng2 = var_idB1[pair2];
	double effcorrection1 = (var_eff1[pair1]*var_eff1[pair2]);
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

//          double pt_pair =sqrt((var_px1[pair1]+var_px1[pair2])*(var_px1[pair1]+var_px1[pair2])+(var_py1[pair1]+var_py1[pair2])*(var_py1[pair1]+var_py1[pair2]));
	if(label_vertex!=99
	   && PassesMuonID(var_tracker1[pair1], muAng1, var_global1[pair1], var_tracker1[pair2], muAng2, var_global1[pair2], var_nhitsTrack1[pair1], var_nhitsTrack1[pair2]) 
           && PassesDptCut(var_dpt1[0])  
           && PassesDphiCut(var_dphi1[0]/pi) 
//	   && pt_pair < 0.8
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

	if(nTrackExclu<1 
		&& PassesTowerCountVeto(nEB,nEE,nHB,nHE,nHFp,nHFm) 
		&& PassesZDCVeto(var_zdcEmMinus1[0],var_zdcEmPlus1[0],var_zdcHadMinus1[0],var_zdcHadPlus1[0])
		&& nGenMuon[0]==2
	){ 
          filter1Events+=effcorrection1;
	  double reso_pt0(-1), reso_pt1(-1), reso_pt01(-1);
          double gen_pt0=sqrt(GenMuon_px[0]*GenMuon_px[0]+GenMuon_py[0]*GenMuon_py[0]);
          double gen_pt1=sqrt(GenMuon_px[1]*GenMuon_px[1]+GenMuon_py[1]*GenMuon_py[1]);

/*          TLorentzVector gen_mu1, gen_mu2, gen_dimuon;
          gen_mu1.SetPtEtaPhiM(gen_pt0,GenMuon_eta[0],GenMuon_phi[0],0.1057);
          gen_mu2.SetPtEtaPhiM(gen_pt1,GenMuon_eta[1],GenMuon_phi[1],0.1057);
          gen_dimuon = gen_mu1 + gen_mu2;*/
          double gen_pt2=GenMuMu_pt[0];//gen_dimuon.Pt();

          TLorentzVector mu1, mu2, dimuon;
          mu1.SetPtEtaPhiM(var_pt1[pair1],var_eta1[pair1],var_phi1[pair1],0.1057);
          mu2.SetPtEtaPhiM(var_pt1[pair2],var_eta1[pair2],var_phi1[pair2],0.1057);
          dimuon = mu1 + mu2;
          double pt2=dimuon.Pt();

//---------------For reco/gen plots-----------------------
	  for(int k=0; k<20; k++){
		if(gen_pt2>=k*0.15 && gen_pt2<(k+1)*0.15) {pTResolution_abs[k]->Fill((double)(pt2-gen_pt2),weight); pTResolution_rel[k]->Fill((double)(pt2/(gen_pt2)),weight);}
	  }
//----------------------------------------------------------

	  } // if nTrack&nCalo if relevant
        }
  }
cout<<"ElEl :"<<endl;
cout<<"  # Dimuon events = "<<filter1Events<<endl;

TCanvas *Greso = new TCanvas("Gauss reso","corel",800,500);
   Greso->SetFillColor(0);
   Greso->SetBorderMode(0);
   Greso->SetBorderSize(2);
   Greso->SetFrameBorderMode(0);
Greso->Divide(4,5);
double pt_v[20], pt_e[20];
double rms_v[20], rms_e[20];
for(int k=1; k<=20; k++){
	Greso->cd(k);
	pTResolution_abs[k-1]->Draw();
	pTResolution_abs[k-1]->Fit("gaus","ILME");
	rms_v[k-1]=gaus->GetParameter(2);
        rms_e[k-1]=gaus->GetParError(2);
	pt_v[k-1]=(k-1)*0.15+0.075;
	pt_e[k-1]=0.075;
}
Greso->Update();
Greso->SaveAs("abs.root");

TCanvas *Greso2 = new TCanvas("Gauss reso2","corel",800,500);
   Greso2->SetFillColor(0);
   Greso2->SetBorderMode(0);
   Greso2->SetBorderSize(2);
   Greso2->SetFrameBorderMode(0);
Greso2->Divide(4,5);
double pt_v2[20], pt_e2[20];
double rms_v2[20], rms_e2[20];
for(int k=1; k<=20; k++){
        Greso2->cd(k);
        pTResolution_rel[k-1]->Draw();
        pTResolution_rel[k-1]->Fit("gaus","ILME");
        rms_v2[k-1]=gaus->GetParameter(2);
        rms_e2[k-1]=gaus->GetParError(2);
        pt_v2[k-1]=(k-1)*0.15+0.075;
        pt_e2[k-1]=0.075;
}
Greso2->Update();
Greso2->SaveAs("rel.root");

TGraphErrors *gAbs = new TGraphErrors(20,pt_v,rms_v,pt_e,rms_e);
TCanvas *Greso3 = new TCanvas("Gauss reso3","corel",800,500);
   Greso3->SetFillColor(0);
   Greso3->SetBorderMode(0);
   Greso3->SetBorderSize(2);
   Greso3->SetFrameBorderMode(0);
gAbs->Draw("APL");
Greso3->SaveAs("reso_abs.root");

TGraphErrors *gAbs2 = new TGraphErrors(20,pt_v2,rms_v2,pt_e2,rms_e2);
TCanvas *Greso4 = new TCanvas("Gauss reso4","corel",800,500);
   Greso4->SetFillColor(0);
   Greso4->SetBorderMode(0);
   Greso4->SetBorderSize(2);
   Greso4->SetFrameBorderMode(0);
gAbs2->Draw("APL");
Greso4->SaveAs("reso_rel.root");

	cout << "END" << endl;   
}

