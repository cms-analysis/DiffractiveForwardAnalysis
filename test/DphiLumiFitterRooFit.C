/*
 * DphiLumiFitterRooFit.C
 * Fit the delta-phi distribution for exclusive dimuons.
 */


/* Here we just setup the CMS plotting style */
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

/* Here's the part that actually does the fit */
void DphiLumiFitterRooFit()
{
  setTDRStyle();

  //  TString flatfile = "inelastic.txt";
  //  TString flatfile = "el-el.txt";
  //  TString flatfile = "lumifitdata_34pb.txt";
  TString flatfile = "lumifitdata_34pb_pt3gev.txt";

  /* Variables to read from the input file: dpT, dphi, and an integer saying which type of event this is (for plotting) */
  RooRealVar* dpt = new RooRealVar("#Delta p_{T} (#mu #mu)","#Delta p_{T}",-0.5,0.5,"GeV");
  RooRealVar* dphi = new RooRealVar("#Delta #phi (#mu #mu)","1 - |#Delta #phi (#mu #mu)| / #pi",-0.1,0.1,"");

  RooRealVar* eventtype = new RooRealVar("eventtype","Event class",0,5,"class");
  RooRealVar* signaltype = new RooRealVar("signaltype","Elastic signal class",1,1,"class");
  RooRealVar* bkgtype = new RooRealVar("bkgtype","Background class",2,4,"class");
  RooRealVar* datatype = new RooRealVar("datatype","Data",5,5,"class");

  /* Construct the signal PDF as 2 Gaussians with mean zero and floating widths */
  RooRealVar* elmdimucoef1 = new RooRealVar("dphi peak","dphi peak",0.0);
  RooRealVar* elmdimucoef2 = new RooRealVar("#Delta #phi width 1","dphi width",0.00480892);
  //  RooRealVar* elmdimucoef2 = new RooRealVar("#Delta #phi width 1","dphi width",0.005,0.0002,0.5);
  RooBreitWigner* elmm = new RooBreitWigner("elastics1","elastics1", *dphi,*elmdimucoef1,*elmdimucoef2);
  RooRealVar* nemm = new RooRealVar("N(elastic #gamma #gamma #rightarrow #mu #mu)","number of signal",700.0,0.0,20000.0);
  //  RooRealVar* nemm = new RooRealVar("N(elastic #gamma #gamma #rightarrow #mu #mu)","number of signal",0);
  //  RooAddPdf *totsig = new RooAddPdf("totsig","total signal PDF",RooArgList(*elmm,*elmm2),RooArgList(*nemm,*nemm));
  RooAddPdf *totsig = new RooAddPdf("totsig","total signal PDF",RooArgList(*elmm),RooArgList(*nemm));

  /* Construct the background PDF as 1 Gaussian with mean zero and floating width*/
  RooRealVar* inelmdimucoef1 = new RooRealVar("Two-photon p0","Two-photon p0",0.0);
  RooRealVar* inelmdimucoef2 = new RooRealVar("background #Delta #phi width ","inelastic dphi width",0.0558694);
  //  RooRealVar* inelmdimucoef2 = new RooRealVar("background #Delta #phi width ","inelastic dphi width",0.053,0.01,0.1);
  RooRealVar* ninmm = new RooRealVar("N(p-dissociation)","number of bkg",100.0,50.0,20000.0);
  //  RooRealVar* ninmm = new RooRealVar("N(background)","number of bkg",0.0);
  RooBreitWigner *inelmm = new RooBreitWigner("inelastics","inelastics",*dphi,*inelmdimucoef1,*inelmdimucoef2);

  /* Alternative - use polynomial shape for the background */
  //  RooPolynomial *inelmm = new RooPolynomial("inelastics","inelastics",*dphi,RooArgList(*inelmdimucoef3),1);

  /* Now make the total signal+background shape as a RooAddPdf */
  RooAddPdf* totshape = new RooAddPdf("totshape","total PDF",RooArgList(*totsig,*inelmm),RooArgList(*nemm,*ninmm));

  /* Now read the dataset from the text file */
  RooDataSet *datatmp = 0;
  datatmp = RooDataSet::read(flatfile,RooArgList(*dphi,*dpt,*eventtype));
  Int_t nentries = datatmp->numEntries();

  /* Now record the true number of elastic signal events being fit */
  RooDataSet *datasig = 0;
  datasig = RooDataSet::read(flatfile,RooArgList(*dpt,*dphi,*signaltype));
  Int_t ntrue = datasig->numEntries();

  RooDataSet *databkg = 0;
  databkg = RooDataSet::read(flatfile,RooArgList(*dpt,*dphi,*bkgtype));

  /* Now do the fit with the "er" option (unbinned extended maximum likelihood fit)!*/
  //  RooFitResult* fitres = 0;
  //  fitres = totshape->fitTo(*datatmp,"er");
  RooFitResult *fitres = totshape->fitTo(*datatmp,RooFit::FitOptions("MHTER"));
  //  fitres->Print();

  /* Now draw the data points and PDFs for signal, background, and total */
  TCanvas *c = new TCanvas("Dilepton signal","Dilepton signal",800,400);
  c->cd();
  RooPlot* xframe = dphi->frame() ;
  xframe->SetTitle(0);
  xframe->SetMaximum(50.0);
  //  datatmp->plotOn(xframe,RooFit::Binning(100));
  datatmp->plotOn(xframe,RooFit::Binning(100));

  totshape->plotOn(xframe,RooFit::Components(*totsig),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  totshape->plotOn(xframe,RooFit::Components(*inelmm),RooFit::LineColor(kBlue),RooFit::LineStyle(kDashed));
  totshape->plotOn(xframe);
  totshape->paramOn(xframe);

  /* Print out the chi2/dof from the plot. */
  Double_t chi2 = xframe->chiSquare(5);
  cout << "Chi^2 = " << chi2 << endl;
  cout << "N(sig input) = " << ntrue << endl;
  cout << "N(total input) = " << nentries << endl;

  xframe->Draw();

  /* Here we do a toy MC study. Take the PDF, generate pseudo-experiments with Poisson fluctuations, and redo the fit.
   * If the fit is unbiased, the pull distribution Nfit-Ntrue/Efit should be a gaussian with mean 0 and width 1.
   */
  if(1) 
    {
      RooMCStudy *toymc = new RooMCStudy(*totshape,RooArgSet(*dphi),RooFit::Binned(kTRUE),RooFit::FitModel(*totshape),RooFit::Silence(),RooFit::Extended(),RooFit::FitOptions(RooFit::Extended()));

      //      RooMCStudy *toymc = new RooMCStudy(*totshape,*totshape,*dphi,"","er");
      toymc->generateAndFit(1000);
      
      TCanvas *c2 = new TCanvas("c2","c2");
      c2->Divide(3,1);
      c2->cd(1);
      RooPlot* pullframe = toymc->plotPull(*nemm,-5.0,5.0,100,kTRUE);
      pullframe->Draw();
      c2->cd(2);
      RooPlot* nsigframe = nemm->frame(0,200);
      toymc->plotParamOn(nsigframe);
      nsigframe->Draw();
      c2->cd(3);
      RooPlot* errsigframe = toymc->plotError(*nemm);
      errsigframe->Draw();
    }
}

