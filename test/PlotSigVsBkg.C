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

void PlotSigVsBkg(Int_t thevar = 6,Int_t theleptonmode, bool save=false)
{
 setTDRStyle();
  TH1F *htmp[8]; 
  if(theleptonmode == 1)
    {
      htmp[0] = GetMuMuHist(thevar,1,save);
      htmp[1] = GetMuMuHist(thevar,2,save);
      htmp[2] = GetMuMuHist(thevar,3,save);
      htmp[3] = GetMuMuHist(thevar,4,save);
      htmp[4] = GetMuMuHist(thevar,5,save);
      htmp[5] = GetMuMuHist(thevar,6,save);
      htmp[6] = GetMuMuHist(thevar,7,save);
      htmp[7] = GetMuMuHist(thevar,8,save);
      htmp[8] = GetMuMuHist(thevar,9,save);

      // Draw stacked histograms for mass plots
      if(thevar > 4 && thevar < 7)
	{
	  htmp[0]->Clone("htot");
	  htot->Add(htmp[1]);
	  htot->Add(htmp[2]);
	  htot->Add(htmp[3]);
	  htot->Add(htmp[4]);
	  htot->Add(htmp[5]);
	  htot->Add(htmp[6]);
	  htot->Add(htmp[7]);
	  htot->Add(htmp[8]);
		htot->SetLineWidth(4);
	  
	  htmp[1]->Clone("hdytot");
	  hdytot->Add(htmp[2]);
	  hdytot->Add(htmp[3]);
	  hdytot->Add(htmp[4]);
	  hdytot->Add(htmp[5]);
	  hdytot->Add(htmp[6]);
	  hdytot->Add(htmp[7]);
	  
	  htmp[5]->Clone("hups1tot");
	  hups1tot->Add(htmp[4]);
	  hups1tot->Add(htmp[6]);
	  hups1tot->Add(htmp[7]);
	  
	  htmp[6]->Clone("hups2tot");
	  hups2tot->Add(htmp[4]);
	  hups2tot->Add(htmp[7]);

	  htmp[7]->Clone("hups3tot");
	  hups3tot->Add(htmp[4]);
	  
	  htot->Draw("hist");
	  hdytot->Draw("histsame");
	  htmp[8]->Draw("histsame");
	  hups1tot->Draw("histsame");
	  hups2tot->Draw("histsame");
	  hups3tot->Draw("histsame");
	  htmp[4]->Draw("histsame");
	  
	  TLegend *l1 = new TLegend(0.55,0.72,0.91,0.91);
	  l1->AddEntry(htot,"#gamma #gamma #rightarrow #mu^{+}#mu^{-} plus backgrounds","lf");
	  l1->AddEntry(hdytot,"Drell-Yan","lf");
	  l1->AddEntry(hups1tot,"#Upsilon (1S)","lf");
	  l1->AddEntry(hups2tot,"#Upsilon (2S)","lf");
	  l1->AddEntry(hups3tot,"#Upsilon (3S)","lf");
	  l1->AddEntry(htmp[4],"W^{+}W^{-}","lf");	  
	  l1->AddEntry(htmp[8],"Singly inelastic #gamma #gamma #rightarrow #mu^{+}#mu^{-}","lf");
	  l1->SetFillColor(0);
	  l1->Draw("same");
	}
      else if(thevar == 8)
	{
          htmp[1]->Clone("htot");
          htot->Add(htmp[2]);
          htot->Add(htmp[3]);
	  htot->Add(htmp[4]);
          htot->Add(htmp[5]);
          htot->Add(htmp[6]);
          htot->Add(htmp[7]);

	  htot->SetLineColor(1);
	  htot->SetMarkerStyle(20);
	  htot->SetLineWidth(3);
	  htot->Draw("e");
	  htmp[0]->Draw("same");
	}
      // Otherwise draw signal and the sum of backgrounds overlaid      
      else
	{	// for elastic mumu
	  TH1F *hsig = htmp[0];
	  TH1F *hbkg = htmp[1];
	  hbkg->Add(htmp[2]);
	  hbkg->Add(htmp[3]);
	  hbkg->Add(htmp[4]);
	  if(thevar == 3)
	    hbkg->Scale(0.08);
	  else if(thevar == 1)
	    hbkg->Scale(0.05);
	  else if(thevar == 7)
	    hbkg->Scale(0.03);
	  else
	    hbkg->Scale(0.02);
	  hbkg->SetMaximum(hsig->GetMaximum() + (0.2 * hsig->GetMaximum()));
	  hbkg->Draw("hist");
	  hsig->Draw("histsame");

	  TLegend *l1 = new TLegend(0.55,0.72,0.91,0.91);
	  l1->AddEntry(hsig,"#gamma #gamma #rightarrow #mu^{+}#mu^{-}","lf");
	  l1->AddEntry(hbkg,"Backgrounds","lf");
	  l1->SetFillColor(kWhite);
	  l1->SetBorderSize(0);
	  l1->Draw("same");

		// for inelastics mumu
/*	  TH1F *hsig = htmp[8];
		hsig->SetFillColor(15);
		hsig->SetFillStyle(3003);
          TLegend *l1 = new TLegend(0.55,0.72,0.91,0.91);
          l1->AddEntry(hsig,"#gamma #gamma #rightarrow #mu^{+}#mu^{-} inel","lf");
          l1->SetFillColor(kWhite);
          l1->SetBorderSize(0);
          l1->Draw("same");	*/

	}
    }

  if(theleptonmode == 2)
    {
      htmp[0] = GetEEHist(thevar,1,save);
      htmp[2] = GetEEHist(thevar,3,save);
      htmp[3] = GetEEHist(thevar,4,save);
      htmp[4] = GetEEHist(thevar,5,save);

      // Draw stacked histograms for mass plots
      if(thevar > 4)
	{	  
	  htmp[0]->Clone("htot");
	  htot->Add(htmp[2]);
	  htot->Add(htmp[3]);
	  htot->Add(htmp[4]);
	  
	  htmp[2]->Clone("hdytot");
	  hdytot->Add(htmp[3]);
	  
	  htot->Draw("hist");
		htot->SetLineWidth(4);
	  hdytot->Draw("histsame");
		hdytot->SetFillColor(10);
	  htmp[4]->Draw("histsame");

	  TLegend *l1 = new TLegend(0.55,0.72,0.91,0.91);
	  l1->AddEntry(htot,"#gamma #gamma #rightarrow e^{+} e^{-} plus backgrounds","lf");
	  l1->AddEntry(hdytot,"Drell-Yan","lf");
	  l1->AddEntry(htmp[4],"W^{+}W^{-}","lf");	  
	  l1->SetFillColor(0);
	  l1->Draw("same");
	}
      // Otherwise draw signal and the sum of backgrounds overlaid
      else
	{
	  TH1F *hsig = htmp[0];
	  TH1F *hbkg = htmp[2];
	  hbkg->Add(htmp[3]);
	  hbkg->Add(htmp[4]);
	  if(thevar == 4 || thevar == 3)
	    hbkg->Scale(0.05);
	  else
	    hbkg->Scale(0.2);
	  hbkg->SetMaximum(hsig->GetMaximum() + (0.2 * hsig->GetMaximum()));
	  hbkg->Draw("hist");
	  hsig->Draw("histsame");

	  TLegend *l1 = new TLegend(0.55,0.72,0.91,0.91);
	  l1->AddEntry(hsig,"#gamma #gamma #rightarrow e^{+}e^{-}","lf");
	  l1->AddEntry(hbkg,"Backgrounds","lf");
	  l1->SetFillColor(kWhite);
	  l1->SetBorderSize(0);
	  l1->Draw("same");
	}
    }

	gPad->SetLogy();
}

// Return histogrammed quantities for mu+mu- samples
TH1F *GetMuMuHist(Int_t plotvar = 6,Int_t physsample = 1, bool save = false)
{
  Double_t ecalocut = 5.0;
  Double_t etcalocut = 0.2; 
  Double_t deltarcut = 0.3;
  Double_t ntrackcut = 3;
  Double_t dptcut = 2.0;
  Double_t dphicut = 2.9;
  Double_t xsec = 0.0;
  Double_t upsloveto = 9.0;
  Double_t upshiveto = 11.0;
  TString st = "";
  Double_t lumi = 100.0;
  Int_t linecolor = 1;
  Int_t linewidth = 0;
  Int_t fillcolor = kWhite; //15;
  Int_t fillstyle = 1001;

  switch(physsample) {
  case 1:
    xsec = 74.65 * lumi / 100000.0;
    st = "gamgammumu.lpair.anal.root";
    linecolor = 1;
    linewidth = 3;
    fillcolor = 15;
    fillstyle = 3003;
    break;
  case 2:
    xsec = 0.5 * 37820.0 * lumi / 10000.0;
    st = "dymumu.610.anal.root";
    fillcolor = 10;
    break;
  case 3:
    xsec = 0.5 * 5951.0 * lumi / 38800.0;
    st = "dymumu.1040.anal.root";
    fillcolor = 10;
    break;
  case 4:
    xsec = 0.5 * 1797.0 * lumi / 22000.0;
    st = "dymumu.40.anal.root";
    fillcolor = 10;
    break;
  case 5:
    xsec = 7.5 * 0.9 * lumi / 1000.0;
    st = "wwmumu.anal.root";
    fillcolor = 9;
    break;
  case 6:
    xsec = 39.0 * lumi / 10000.0;
    st = "upsilonmumu.starlight.anal.root";
    fillcolor = 3;
    break;
  case 7:
    xsec = 13.0 * lumi / 3359.0;
    st = "upsilonmumu.starlight2s.anal.root";
    fillcolor = 6;
    break;
  case 8:
    xsec = 10.0 * lumi / 2636.0;
    st = "upsilonmumu.starlight3s.anal.root";
    fillcolor = 1;
    break;
  case 9:
    xsec = 76.2 * lumi / 20000.0;
    st = "gamgammumu.lpairinelastic.anal.root";
    fillcolor = 5;
    break;
  default:
    break;
  }

  TFile *f1 = new TFile(st);

  TTree *tr1 = f1->Get("ntp1");

  TH1F *hmll;
  TH1F *hmll2;
  TH1F *hdphi;
  TH1F *hdpt;
  TH1F *hnextracaloe;
  TH1F *hnextracaloet;
  TH1F *hntrack;
  TH1F *hncalofinal;

  hmll = new TH1F("mll","mll",100,0,200);
  hmll2 = new TH1F("mll2","mll2",100,0,200);
  hdphi = new TH1F("hdphi","hdphi",192,0.0,3.2);
  hdpt = new TH1F("hdpt","hdpt",240,-20,20);
  hnextracaloe = new TH1F("hnextracaloe","hnextracaloe",50,0,50);
  hnextracaloet = new TH1F("hnextracaloet","hnextracaloet",50,0,50);
  hntrack = new TH1F("hntrack","hntrack",50,0,50);
  hncalofinal = new TH1F("hncalofinal","hncalofinal",50,0,50);

  TH2F *hntrackncalo = new TH2F("hntrackncalo","hntrackncalo",500,0,100,500,0,100);
  TH2F *hdr = new TH2F("hdr","hdr",400,-6.0,6.0,400,-6.0,6.0);  
  TH1F *hjete = new TH1F("hjete","hjete",50,0,200);

  Double_t mumumass;
  Double_t mumudphi;
  Double_t mupt[2];
  Double_t mueta[2];
  Double_t muphi[2];
  Int_t nmuons;
  Int_t njets;
  Int_t ntracks;
  Int_t ncalo;
  Double_t deltara;
  Double_t deltarb;

  Double_t caloe[1000];
  Double_t caloet[1000];
  Double_t caloeta[1000];
  Double_t calophi[1000];

  tr1->SetBranchAddress("MuonCand_pt",&mupt);
  tr1->SetBranchAddress("MuMu_mass",&mumumass);
  tr1->SetBranchAddress("MuMu_dphi",&mumudphi);
  tr1->SetBranchAddress("MuonCand_eta",&mueta);
  tr1->SetBranchAddress("MuonCand_phi",&muphi);
  tr1->SetBranchAddress("nJetCand",&njets);
  tr1->SetBranchAddress("nTrackCand",&ntracks);
  tr1->SetBranchAddress("nMuonCand",&nmuons);
  tr1->SetBranchAddress("CaloTower_e",&caloe);
  tr1->SetBranchAddress("CaloTower_et",&caloet);
  tr1->SetBranchAddress("CaloTower_eta",&caloeta);
  tr1->SetBranchAddress("CaloTower_phi",&calophi);
  tr1->SetBranchAddress("nCaloCand",&ncalo);

  Int_t ent1 = tr1->GetEntries();
  Int_t nisocalo = 0;
  Int_t ntotsig = tr1->GetEntries();
  Int_t ndphipass, ndptpass, ncalopass, ntrackpass = 0;
  Int_t nupspass = 0;
  Int_t ntotpass = 0;

  for(Int_t i = 0;i < ent1;i++)
    {
      nisocalo = 0;
      tr1->GetEntry(i);
      hdphi->Fill(mumudphi);
      hdpt->Fill((mupt[0]-mupt[1]));
      hmll->Fill(mumumass);
      hntrack->Fill(ntracks);

      for(Int_t x = 0;x < ncalo;x++)
	{
	  deltara = sqrt((caloeta[x]-mueta[0])*(caloeta[x]-mueta[0]) + (calophi[x]-muphi[0])*(calophi[x]-muphi[0]));
	  deltarb = sqrt((caloeta[x]-mueta[1])*(caloeta[x]-mueta[1]) + (calophi[x]-muphi[1])*(calophi[x]-muphi[1]));
	  if(deltara > 0.3 && deltarb > 0.3 && caloe[x] > ecalocut)
	  //  if(deltara > 0.3 && deltarb > 0.3 && caloet[x] > etcalocut) ///////// change here
	    nisocalo++;

	  if(caloe[x] > ecalocut)
	    {
	      if(deltara < deltarb)
		{
		  hdr->Fill((caloeta[x]-mueta[0]),(calophi[x]-muphi[0]));
		}
	      else
		{
		  hdr->Fill((caloeta[x]-mueta[1]),(calophi[x]-muphi[1]));
		}
	    }
	}

      hnextracaloe->Fill(nisocalo);

      if(mumudphi > dphicut)
	{
	  ndphipass++;
	  if(fabs(mupt[0]-mupt[1]) < dptcut)
	    {
	      ndptpass++;
	      if(nisocalo < 5) 
		{
		  ncalopass++;
		  if(ntracks < ntrackcut)
		    {
		      ntrackpass++;
		      ntotpass++;
		      if(mumumass < upsloveto || mumumass > upshiveto)
			{
			  nupspass++;
			  hmll2->Fill(mumumass);
			}
		    }
		}
	      if(ntracks < 200)
		{
		  if(mumumass < upsloveto || mumumass > upshiveto)
		    hncalofinal->Fill(nisocalo);
		}
	    }
	}
    }

  Double_t eff = (Double_t)ntrackpass/(Double_t)ntotsig;
  cout << "Efficiency/rejection for file: " << st << endl;
  cout << "\td(phi) eff = " << (Double_t)ndphipass/(Double_t)ntotsig << endl;
  cout << "\td(pt) eff = " << (Double_t)ndptpass/(Double_t)ntotsig << endl;
  cout << "\tcalo excl eff = " << (Double_t)ncalopass/(Double_t)ntotsig << endl;
  cout << "\ttrack excl eff = " << eff << endl;
  cout << "\tUpsilon veto eff = " << (Double_t)nupspass/(Double_t)ntotsig << endl;
  //  cout << "Eff = " << eff << endl;
  cout << "N(sig) = " << (Double_t)nupspass * xsec << " +- " << sqrt(nupspass) * xsec << endl;


  if(0){
  gStyle->SetPalette(51);
   hdr->SetStats(0);
   hdr->SetTitle(0);
   hdr->SetXTitle("#Delta #eta (#mu - Calo Tower)");
   hdr->SetYTitle("#Delta #phi (#mu - Calo Tower)");
   hdr->Scale(5951.0/1797.0);
   hdr->Add(hdrdylo);
  }

  if(0){
  gStyle->SetPalette(51);
  hntrackncalo->Scale((1797.0/5951.0)*(10000.0/15000.0));  
  hntrackncalo->Add(hntrackncalo);
  hntrackncalo->SetXTitle("N(extra calo towers)");
  hntrackncalo->SetYTitle("N(tracks)");
  hntrackncalo->SetStats(0);
  hntrackncalo->SetTitle(0);
  hntrackncalo->Draw("");
  }

  if(plotvar == 1){
  hntrack->SetFillColor(4);
  hntrack->SetMaximum(20000);
  hntrack->Scale(1);
  hntrack->SetStats(0);
  hntrack->SetTitle(0);
  hntrack->SetXTitle("Track multiplicity");
  hntrack->Sumw2();
  hntrack->Scale(xsec);
  hntrack->Scale(100.0);
  hntrack->SetLineColor(linecolor);
  hntrack->SetFillColor(fillcolor);
  hntrack->SetFillStyle(fillstyle);
  hntrack->SetLineWidth(linewidth);
  hntrack->Draw("hist");
  return(hntrack);
  }

  if(plotvar == 2){
  hnextracaloe->SetFillColor(fillcolor);
  hnextracaloe->SetFillStyle(fillstyle);
  hnextracaloe->SetMaximum(10000);
  hnextracaloe->Scale(1);
  hnextracaloe->SetStats(0);
  hnextracaloe->SetTitle(0);
  hnextracaloe->SetXTitle("Tower multiplicity");
  hnextracaloe->Sumw2();
  hnextracaloe->Scale(xsec);
  hnextracaloe->SetLineColor(linecolor);
  hnextracaloe->SetLineWidth(linewidth);
  hnextracaloe->Scale(200.0);
  hnextracaloe->Draw("hist");
  return(hnextracaloe);
  }

  if(plotvar == 3){
  hdpt->SetFillColor(fillcolor);
  hdpt->SetFillStyle(fillstyle);
  hdpt->SetMaximum(6000);
  hdpt->SetStats(0);
  hdpt->SetTitle(0);
  hdpt->SetXTitle("|#Delta p_{T}(#mu #mu)|");
  hdpt->Draw("hist");
  hdpt->Sumw2();
  hdpt->Scale(xsec);
  hdpt->Scale(40);
  hdpt->SetLineColor(linecolor);
  hdpt->SetLineWidth(linewidth);
  hdpt->Draw("hist");
  return(hdpt);
  }

  if(plotvar == 4){
  hdphi->SetFillColor(fillcolor);
  hdphi->SetFillStyle(fillstyle);
  hdphi->SetMaximum(8000);
  hdphi->SetStats(0);
  hdphi->SetTitle(0);
  hdphi->SetXTitle("|#Delta #phi(#mu #mu)|");
  hdphi->Sumw2();
  hdphi->Scale(xsec);
  hdphi->Scale(40);
  hdphi->SetLineColor(linecolor);
  hdphi->SetLineWidth(linewidth);
  hdphi->Draw("hist");
  return(hdphi);
  }

  if(plotvar == 5){
  hmll2->Sumw2();
  hmll2->Scale(xsec);
  hmll2->SetFillColor(fillcolor);
  hmll2->SetFillStyle(fillstyle);
  hmll2->SetLineColor(linecolor);
  hmll2->SetLineWidth(linewidth);
  hmll2->SetMinimum(0.001);
  hmll2->SetMaximum(50.0*lumi);
  hmll2->SetTitle(0);
  hmll2->SetStats(0);
  hmll2->SetXTitle("m (#mu #mu)");
  hmll2->Draw("hist");
  return(hmll2);
  }

  if(plotvar == 6){
  hmll->Sumw2();
  hmll->Scale(xsec);
  hmll->SetMinimum(0.001);
  hmll->SetMaximum(20.0*lumi);
  hmll->SetTitle(0);
  hmll->SetStats(0);
  hmll->SetXTitle("m (#mu #mu)");
  hmll->SetFillColor(fillcolor);
  hmll->SetFillStyle(fillstyle);
  hmll->SetLineColor(linecolor);
  hmll->SetLineWidth(linewidth);
  hmll->Draw("hist");
  return(hmll);
  }

  if(plotvar == 7){
  hnextracaloet->SetFillColor(fillcolor);
  hnextracaloet->SetFillStyle(fillstyle);
  hnextracaloet->SetMaximum(10000);
  hnextracaloet->Scale(1);
  hnextracaloet->SetStats(0);
  hnextracaloet->SetTitle(0);
  hnextracaloet->SetXTitle("Tower multiplicity");
  hnextracaloet->Sumw2();
  hnextracaloet->Scale(xsec);
  hnextracaloet->SetLineColor(linecolor);
  hnextracaloet->SetLineWidth(linewidth);
  hnextracaloet->Scale(200.0);
  hnextracaloet->Draw("hist");
  return(hnextracaloet);
  }

  if(plotvar == 8){
    hncalofinal->SetFillColor(fillcolor);
    hncalofinal->SetFillStyle(fillstyle);
    hncalofinal->SetMaximum(10000);
    hncalofinal->Scale(xsec);
    hncalofinal->SetStats(0);
    hncalofinal->SetTitle(0);
    hncalofinal->SetXTitle("Tower multiplicity");
    hncalofinal->Sumw2();
    hncalofinal->Scale(xsec);
    hncalofinal->SetLineColor(linecolor);
    hncalofinal->SetLineWidth(linewidth);
    hncalofinal->Scale(200.0);
    hncalofinal->Draw("hist");
    return(hncalofinal);
  }

}

// Return histogrammed quantities for e+e- events
TH1F *GetEEHist(Int_t plotvar = 6,Int_t physsample = 1,bool save=false)
{
  Double_t ecalocut = 5.0;
  Double_t etcalocut = 2.0;
  Double_t deltarcut = 0.3;
  Double_t ntrackcut = 3;
  Double_t dptcut = 5.0;
  Double_t dphicut = 2.7;
  Double_t xsec = 0.0;
  TString st = "";
  Double_t lumi = 100.0;
  Int_t linecolor = 1;
  Int_t linewidth = 0;
  Int_t fillcolor = kWhite; //15;
  Int_t fillstyle = 1001;

  switch(physsample) {
  case 1:
    xsec = 10.35 * lumi / 100000.0;
    st = "gamgamee.lpair.anal.root";
    linecolor = 1;
    linewidth = 3;
    fillcolor = 15;
    fillstyle = 3003;
    break;
  case 2:
    xsec = 0.5 * 37820.0 * lumi / 10000.0;
    st = "dyee.610.anal.root";
    fillcolor = 10;
    break;
  case 3:
    xsec = 0.5 * 5951.0 * lumi / 38800.0;
    st = "dyee.1040.anal.root";
    fillcolor = 10;
    break;
  case 4:
    xsec = 0.5 * 1797.0 * lumi / 22000.0;
    st = "dyee.40.anal.root";
    fillcolor = 10;
    break;
  case 5:
    xsec = 7.5 * 0.9 * lumi / 1000.0;
    st = "wwee.anal.root";
    fillcolor = 9;
    break;
  case 6:
    break;
  case 7:
    break;
  case 8:
    break;
  default:
    break;
  }

  TFile *f1 = new TFile(st);

  TTree *tr1 = f1->Get("ntp1");

  TH1F *hmll;
  TH1F *hmll2;
  TH1F *hdphi;
  TH1F *hdpt;
  TH1F *hnextracaloe;
  TH1F *hnextracaloet;
  TH1F *hntrack;

  hmll = new TH1F("mll","mll",100,0,200);
  hmll2 = new TH1F("mll2","mll2",100,0,200);
  hdphi = new TH1F("hdphi","hdphi",196,0.0,3.2);
  hdpt = new TH1F("hdpt","hdpt",240,-20,20);
  hnextracaloe = new TH1F("hnextracaloe","hnextracaloe",50,0,50);
  hnextracaloet = new TH1F("hnextracaloet","hnextracaloet",50,0,50);
  hntrack = new TH1F("hntrack","hntrack",50,0,50);

  TH2F *hntrackncalo = new TH2F("hntrackncalo","hntrackncalo",500,0,100,500,0,100);
  TH2F *hdr = new TH2F("hdr","hdr",400,-6.0,6.0,400,-6.0,6.0);  
  TH1F *hjete = new TH1F("hjete","hjete",50,0,200);

  Double_t elelmass;
  Double_t eleldphi;
  Double_t elept[2];
  Double_t eleeta[2];
  Double_t elephi[2];
  Int_t nelectrons;
  Int_t njets;
  Int_t ntracks;
  Int_t ncalo;
  Double_t deltara;
  Double_t deltarb;

  Double_t caloe[1000];
  Double_t caloeta[1000];
  Double_t calophi[1000];

  tr1->SetBranchAddress("EleCand_et",&elept);
  tr1->SetBranchAddress("ElEl_mass",&elelmass);
  tr1->SetBranchAddress("ElEl_dphi",&eleldphi);
  tr1->SetBranchAddress("EleCand_eta",&eleeta);
  tr1->SetBranchAddress("EleCand_phi",&elephi);
  tr1->SetBranchAddress("nJetCand",&njets);
  tr1->SetBranchAddress("nTrackCand",&ntracks);
  tr1->SetBranchAddress("nEleCand",&nelectrons);
  tr1->SetBranchAddress("CaloTower_e",&caloe);
  tr1->SetBranchAddress("CaloTower_eta",&caloeta);
  tr1->SetBranchAddress("CaloTower_phi",&calophi);
  tr1->SetBranchAddress("nCaloCand",&ncalo);

  Int_t ent1 = tr1->GetEntries();
  Int_t nisocalo = 0;
  Int_t ntotsig = tr1->GetEntries();
  Int_t ndphipass, ndptpass, ncalopass, ntrackpass = 0;
  Int_t ntotpass = 0;

  for(Int_t i = 0;i < ent1;i++)
    {
      nisocalo = 0;
      tr1->GetEntry(i);
      hdphi->Fill(eleldphi);
      hdpt->Fill((elept[0]-elept[1]));
      hmll->Fill(elelmass);
      hntrack->Fill(ntracks);

      for(Int_t x = 0;x < ncalo;x++)
	{
	  deltara = sqrt((caloeta[x]-eleeta[0])*(caloeta[x]-eleeta[0]) + (calophi[x]-elephi[0])*(calophi[x]-elephi[0]));
	  deltarb = sqrt((caloeta[x]-eleeta[1])*(caloeta[x]-eleeta[1]) + (calophi[x]-elephi[1])*(calophi[x]-elephi[1]));
	  if(deltara > 0.3 && deltarb > 0.3 && caloe[x] > ecalocut)
	    nisocalo++;

	  if(caloe[x] > ecalocut)
	    {
	      if(deltara < deltarb)
		{
		  hdr->Fill((caloeta[x]-eleeta[0]),(calophi[x]-elephi[0]));
		}
	      else
		{
		  hdr->Fill((caloeta[x]-eleeta[1]),(calophi[x]-elephi[1]));
		}
	    }
	}

      hnextracaloe->Fill(nisocalo);

      if(eleldphi > dphicut)
	{
	  ndphipass++;
	  if(fabs(elept[0]-elept[1]) < dptcut)
	    {
	      ndptpass++;
	      if(nisocalo < 5)
		{
		  ncalopass++;
		  if(ntracks < ntrackcut)
		    {
		      ntotpass++;
		      hmll2->Fill(elelmass);
		    }
		}
	    }
	}
    }

  Double_t eff = (Double_t)ntotpass/(Double_t)ntotsig;
  cout << "Efficiency/rejection for file: " << st << endl;
  cout << "\td(phi) eff = " << (Double_t)ndphipass/(Double_t)ntotsig << endl;
  cout << "\td(pt) eff = " << (Double_t)ndptpass/(Double_t)ntotsig << endl;
  cout << "\tcalo excl eff = " << (Double_t)ncalopass/(Double_t)ntotsig << endl;
  cout << "\tTotal eff = " << eff << endl;
  //  cout << "Eff = " << eff << endl;
  cout << "N(sig) = " << (Double_t)ntotpass * xsec << " +- " << sqrt(ntotpass) * xsec << endl;

  if(0){
  gStyle->SetPalette(51);
   hdr->SetStats(0);
   hdr->SetTitle(0);
   hdr->SetXTitle("#Delta #eta (#ele - Calo Tower)");
   hdr->SetYTitle("#Delta #phi (#ele - Calo Tower)");
   hdr->Scale(5951.0/1797.0);
   hdr->Add(hdrdylo);
  }

  if(0){
  gStyle->SetPalette(51);
  hntrackncalo->Scale((1797.0/5951.0)*(10000.0/15000.0));  
  hntrackncalo->Add(hntrackncalo);
  hntrackncalo->SetXTitle("N(extra calo towers)");
  hntrackncalo->SetYTitle("N(tracks)");
  hntrackncalo->SetStats(0);
  hntrackncalo->SetTitle(0);
  hntrackncalo->Draw("");
  }

  if(plotvar == 1){
  hntrack->SetFillColor(4);
  hntrack->SetMaximum(20000);
  hntrack->Scale(1);
  hntrack->SetStats(0);
  hntrack->SetTitle(0);
  hntrack->SetXTitle("Track multiplicity");
  hntrack->Sumw2();
  hntrack->Scale(xsec);
  hntrack->Scale(100.0);
  hntrack->SetLineColor(linecolor);
  hntrack->SetFillColor(fillcolor);
  hntrack->SetFillStyle(fillstyle);
  hntrack->SetLineWidth(linewidth);
  hntrack->Draw("hist");
  return(hntrack);
  }

  if(plotvar == 2){
  hnextracaloe->SetFillColor(fillcolor);
  hnextracaloe->SetFillStyle(fillstyle);
  hnextracaloe->SetMaximum(10000);
  hnextracaloe->Scale(1);
  hnextracaloe->SetStats(0);
  hnextracaloe->SetTitle(0);
  hnextracaloe->SetXTitle("Tower multiplicity");
  hnextracaloe->Sumw2();
  hnextracaloe->Scale(xsec);
  hnextracaloe->SetLineColor(linecolor);
  hnextracaloe->SetLineWidth(linewidth);
  hnextracaloe->Scale(200.0);
  hnextracaloe->Draw("hist");
  return(hnextracaloe);
  }

  if(plotvar == 3){
  hdpt->SetFillColor(fillcolor);
  hdpt->SetFillStyle(fillstyle);
  hdpt->SetMaximum(6000);
  hdpt->SetStats(0);
  hdpt->SetTitle(0);
  hdpt->SetXTitle("|#Delta E_{T}(e e)|");
  hdpt->Draw("hist");
  hdpt->Sumw2();
  hdpt->Scale(xsec);
  hdpt->Scale(40);
  hdpt->SetLineColor(linecolor);
  hdpt->SetLineWidth(linewidth);
  hdpt->Draw("hist");
  return(hdpt);
  }

  if(plotvar == 4){
  hdphi->SetFillColor(fillcolor);
  hdphi->SetFillStyle(fillstyle);
  hdphi->SetMaximum(8000);
  hdphi->SetStats(0);
  hdphi->SetTitle(0);
  hdphi->SetXTitle("|#Delta #phi(e e)|");
  hdphi->Sumw2();
  hdphi->Scale(xsec);
  hdphi->Scale(40);
  hdphi->SetLineColor(linecolor);
  hdphi->SetLineWidth(linewidth);
  hdphi->Draw("hist");
  return(hdphi);
  }

  if(plotvar == 5){
  hmll2->Sumw2();
  hmll2->Scale(xsec);
  hmll2->SetFillColor(fillcolor);
  hmll2->SetFillStyle(fillstyle);
  hmll2->SetLineColor(linecolor);
  hmll2->SetLineWidth(linewidth);
  hmll2->SetMinimum(0.001);
  hmll2->SetMaximum(5.0*lumi);
  hmll2->SetTitle(0);
  hmll2->SetStats(0);
  hmll2->SetXTitle("m (e e)");
  hmll2->Draw("hist");
  return(hmll2);
  }

  if(plotvar == 6){
  hmll->Sumw2();
  hmll->Scale(xsec);
  hmll->SetMinimum(0.001);
  hmll->SetMaximum(2.0*lumi);
  hmll->SetTitle(0);
  hmll->SetStats(0);
  hmll->SetXTitle("m (e e)");
  hmll->SetFillColor(fillcolor);
  hmll->SetFillStyle(fillstyle);
  hmll->SetLineColor(linecolor);
  hmll->SetLineWidth(linewidth);
  hmll->Draw("hist");
  return(hmll);
  }
}
