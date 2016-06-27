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
  tdrStyle->cd();

}


void DileptonFitter()
{
  setTDRStyle();
  //  using namespace RooFit;

  TString flatfile = "upsilondata.txt";
  //  TString flatfile = "upsilonmc.txt";

  RooRealVar* mMuMu = new RooRealVar("mMuMu","Dimuon mass",8.0,12.0,"GeV");
  RooRealVar* mMuMuFull = new RooRealVar("mMuMuFull","Dimuon mass",6.0,50.0,"GeV");
  
  RooRealVar* mdimucoef1 = new RooRealVar("Two-photon p0","Two-photon p0",0.0,-50000.0,50000.0);
  RooRealVar* mdimucoef2 = new RooRealVar("Two-photon p1","Two-photon p1",0.0,-20000.0,20000.0);

  RooRealVar* ups1smdimucoef1 = new RooRealVar("Upsilon(1S) mass","Upsilon(1S) mass",9.6,9.0,10.0);
  RooRealVar* ups1smdimucoef2 = new RooRealVar("Upsilon(1S) width","Upsilon(1S) width",1.0,0.0,3.0);
  //  RooRealVar* ups2smdimucoef1 = new RooRealVar("Upsilon(2S) mass","Upsilin(2S) mass",10.0,9.7,10.3);
  //  RooRealVar* ups2smdimucoef2 = new RooRealVar("Upsilon(2S) width","Upsilon(2S) width",1.0,0.0,3.0);
  //  RooRealVar* ups3smdimucoef1 = new RooRealVar("Upsilon(3S) mass","Upsilon(3S) mass",10.4,10.2,10.7);
  //  RooRealVar* ups3smdimucoef2 = new RooRealVar("Upsilon(3S) width","Upsilon(3S) width",1.0,0.0,3.0);
  RooRealVar* ups2smdimucoef1 = new RooRealVar("Upsilon(2S) mass","Upsilin(2S) mass",10.0207);
  RooRealVar* ups2smdimucoef2 = new RooRealVar("Upsilon(2S) width","Upsilon(2S) width",0.0914582);
  RooRealVar* ups3smdimucoef1 = new RooRealVar("Upsilon(3S) mass","Upsilon(3S) mass",10.3498);
  RooRealVar* ups3smdimucoef2 = new RooRealVar("Upsilon(3S) width","Upsilon(3S) width",0.117446);

  //  RooRealVar* nmm = new RooRealVar("N(two-photon events)","number of signal events",0);
  RooRealVar* nmm = new RooRealVar("N(two-photon events)","number of signal events",1000.0,0.0,20000.0); 
  RooRealVar* nu1s = new RooRealVar("N(Upsilon(1S))","number of Upsilon 1S",50.0,0.0,20000.0); 
  RooRealVar* nu2s = new RooRealVar("N(Upsilon(2S))","number of Upsilon 2S",100.0,0.0,20000.0);
  RooRealVar* nu3s = new RooRealVar("N(Upsilon(3S))","number of Upsilon 3S",100.0,0.0,20000.0);

  RooPolynomial* mm = new RooPolynomial("mm shape","mm shape", *mMuMu, RooArgList(*mdimucoef1, *mdimucoef2),2);
  RooGaussian* ups1 = new RooGaussian("ups1s sigpeak","ups1s sigpeak", *mMuMu,*ups1smdimucoef1,*ups1smdimucoef2);
  RooGaussian* ups2 = new RooGaussian("ups2s sigpeak","ups2s sigpeak", *mMuMu,*ups2smdimucoef1,*ups2smdimucoef2);
  RooGaussian* ups3 = new RooGaussian("ups3s sigpeak","ups3s sigpeak", *mMuMu,*ups3smdimucoef1,*ups3smdimucoef2);

  RooAddPdf* totshape = new RooAddPdf("totshape","total PDF",RooArgList(*ups1,*ups2,*ups3,*mm),RooArgList(*nu1s,*nu2s,*nu3s,*nmm));

  RooDataSet *datatmp = 0;
  datatmp = RooDataSet::read(flatfile,RooArgList(*mMuMu));
  datatmpfull = RooDataSet::read(flatfile,RooArgList(*mMuMuFull));

  //  RooFitResult* fitres = 0;
  //  fitres = totshape->fitTo(*datatmp,"er","","8,12");
  //  fitres = totshape->fitTo(*datatmp,"er");
  //  fitres->Print();
  RooFitResult *fitres = totshape->fitTo(*datatmp,RooFit::FitOptions("MHTER")); 

  TCanvas *c = new TCanvas("Dilepton signal","Dilepton signal",800,400);
  //  c->Divide(2,1);
  c->cd();
  RooPlot* xframe = mMuMu->frame() ;
  //  RooPlot* yframe = mMuMuFull->frame();
  //  c->cd(1);
  xframe->SetMaximum(15);
  datatmp->plotOn(xframe,RooFit::Binning(30));
  //  datatmp->plotOn(xframe,RooFit::Binning(100)); 
  //  totshape->plotOn(xframe);
  //  mm-plotOn(xframe,LineColor(kRed),LineStyle(kDashed));
  //  mm->plotOn(xframe);
  totshape->plotOn(xframe,RooFit::Components(*ups1),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  totshape->plotOn(xframe,RooFit::Components(*ups2),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  totshape->plotOn(xframe,RooFit::Components(*ups3),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  totshape->plotOn(xframe,RooFit::Components(*mm),RooFit::LineColor(kBlue),RooFit::LineStyle(kDashed));
  totshape->plotOn(xframe);
  totshape->paramOn(xframe);

  xframe->Draw() ;
  //  c->cd(2);
  //  datatmpfull->plotOn(yframe);
  //  yframe->Draw();
  c->SaveAs("SigPlusBGFit.gif");
}
