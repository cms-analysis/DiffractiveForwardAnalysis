#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

/*
 * Template fit yields
 * 0.02-0.03:    20.1(45) 0.2(56)
 * 0.03-0.04:    21.8(45) 20.0(56)
 * 0.04-0.06:    11.0(45) 12.5(56) 
 * > 0.06:       19.5(45) 10.7(56) 
 * > 0.04:       32.4(45) 25.3(56)
 * > 0.02:       63.3(45)
 * > 0.03:       40.5(56)
 * > 0.04 mult:  10.4(45) 7.5(56)
 * > 0.02,0.03:  99.7(45+56)
 * > 0.04 mult:  17.9(45+56) 
 */

void FitBinnedPlots(TString hist = "hressum", Int_t rebinfact = 5, Int_t year = 2017, Int_t mode = 1, Int_t pot = 220, Float_t fitrange = 0.5, Int_t minimizer = 0)
{
  TFile *f0,*f1,*f2;
  if(year == 2017)
    {
      if(pot == 220)
	{
          f0 = TFile::Open("MoreDimuons2017BCDEFSingleTrackPixelsLegacyFinalFromDBWithMass.root");
	  f1 = TFile::Open("../macros/MoreDimuons2017BCDEFSingleTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracks.root");
	  f2 = TFile::Open("MoreDimuons2017preTS2MCallanglesSingleTrack220LegacyFinalFromDBWithMass.root");
	}
      if(pot == 210)
	{
	  f0 = TFile::Open("MoreDimuons2017BCDEFSingleTrackStripsLegacyFinalFromDBWithMass.root");
	  f1 = TFile::Open("MoreDimuons2017BCDEFSingleTrackNearLegacyFinalFromDBWithMassAntiAcop1to5tracks.root");
          f2 = TFile::Open("MoreDimuons2017MCapreTS2llanglesSingleTrack210LegacyFinalFromDBWithMass.root");
	}
    }
  if(year == 2018)
    {
      if(pot == 220)
	{
	  f0 = TFile::Open("MoreDimuons2018ABCDSingleTrackPixelsFarLegacyFinalFromDBWithMass.root");
	  f1 = TFile::Open("MoreDimuons2018ABCDSingleTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracks.root");
	  f2 = TFile::Open("MoreDimuons2018MCallanglesWithMultiLegacyFinalFromDBWithMass.root");
	}
      if(pot == 210)
	{
	  f0 = TFile::Open("MoreDimuons2018ApartialBCDSingleTrackPixelsNearLegacyFinalFromDBWithMass.root");
	  f1 = TFile::Open("MoreDimuons2018ApartialBCDSingleTrackNearLegacyFinalFromDBWithMassAntiAcop1to5tracks.root");
	  f2 = TFile::Open("MoreDimuons2018MCallanglesSingleTrack210LegacyFinalFromDBWithMass.root");
	}
    }
  if(year == 20172018)
    {
      //      f0 = TFile::Open("MoreDimuons2017BCDEF2018ApartialBCDWithMultiTrackFarLegacyFinalFromDBWithMassMoreBins.root");
      //      f0 = TFile::Open("MoreDimuons2018ApartialBCDWithMultiTrackNearLegacyFinalFromDBWithMassMoreBins.root");
      //      f1 = TFile::Open("MoreDimuons2017BCDEF2018BCDWithMultiTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracksMoreBins.root");
      if(pot == 220)
	{
	  f0 = TFile::Open("MoreDimuons2017BCDEF2018ABCDSingleTrackPixelsFarLegacyFinalFromDBWithMass.root");
	  f1 = TFile::Open("MoreDimuons2017BCDEF2018ABCDSingleTrackFarLegacyFinalFromDBWithMassAntiAcop1to5tracks.root");
	  //	  f2 = TFile::Open("../../../../../CMSSW_10_6_0/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/systematicsmacros/MoreDimuons20172018MCallanglesSingleTrack220LegacyFinalFromDBWithMass.root");
	  //	  f2 = TFile::Open("MoreDimuons20172018MCallanglesWithMultiLegacyFinalFromDBWithMass.root");
	  f2 = TFile::Open("MoreDimuons20172018MCallanglesWithMultiLegacyFinalFromDBWithMassMoreBins.root");
	}
      if(pot == 210)
	{
	  f0 = TFile::Open("MoreDimuons2017BCDEF2018ApartialBCDSingleTrack210LegacyFinalFromDBWithMass.root");
	  f1 = TFile::Open("MoreDimuons2017BCDEF2018ApartialBCDSingleTrack210LegacyFinalFromDBWithMassAntiAcop1to5tracks.root");
          f2 = TFile::Open("MoreDimuons20172018MCallanglesSingleTrack210LegacyFinalFromDBWithMass.root");
	}
    }

  //  TH1F *hb145 = (TH1F *)f0->Get("hbin145");
  TH1F *hall;
  TH1F *hbkg;
  TH1F *hsig;

  Float_t siglo = -0.2;
  Float_t sighi = 0.2;
  Float_t bkglo = -1.5;
  Float_t bkghi = -0.5;
  Float_t fitmax = 25;
  Float_t submax = 25;

  if(mode == 1)
    {
      hall = (TH1F *)f0->Get(hist);
      hbkg = (TH1F *)f1->Get(hist);
      hsig = (TH1F *)f2->Get(hist);
    }

  TH1F *hsumall145 = (TH1F *)f0->Get("hbin145"); TH1F *hsumall245 = (TH1F *)f0->Get("hbin245"); 
  TH1F *hsumall345 = (TH1F *)f0->Get("hbin345"); TH1F *hsumall445 = (TH1F *)f0->Get("hbin445");
  TH1F *hsumall156 = (TH1F *)f0->Get("hbin156"); TH1F *hsumall256 = (TH1F *)f0->Get("hbin256"); 
  TH1F *hsumall356 = (TH1F *)f0->Get("hbin356"); TH1F *hsumall456 = (TH1F *)f0->Get("hbin456");

  TH1F *hsumbkg145 = (TH1F *)f1->Get("hbin145"); TH1F *hsumbkg245 = (TH1F *)f1->Get("hbin245"); 
  TH1F *hsumbkg345 = (TH1F *)f1->Get("hbin345"); TH1F *hsumbkg445 = (TH1F *)f1->Get("hbin445");
  TH1F *hsumbkg156 = (TH1F *)f1->Get("hbin156"); TH1F *hsumbkg256 = (TH1F *)f1->Get("hbin256"); 
  TH1F *hsumbkg356 = (TH1F *)f1->Get("hbin356"); TH1F *hsumbkg456 = (TH1F *)f1->Get("hbin456");

  TH1F *hsumsig145 = (TH1F *)f2->Get("hbin145"); TH1F *hsumsig245 = (TH1F *)f2->Get("hbin245"); 
  TH1F *hsumsig345 = (TH1F *)f2->Get("hbin345"); TH1F *hsumsig445 = (TH1F *)f2->Get("hbin445");
  TH1F *hsumsig156 = (TH1F *)f2->Get("hbin156"); TH1F *hsumsig256 = (TH1F *)f2->Get("hbin256"); 
  TH1F *hsumsig356 = (TH1F *)f2->Get("hbin356"); TH1F *hsumsig456 = (TH1F *)f2->Get("hbin456");

  TH1F *hsummultall45 = (TH1F *)f0->Get("hbin3mult45"); 
  TH1F *hsummultbkg45 = (TH1F *)f1->Get("hbin3mult45");
  TH1F *hsummultsig45 = (TH1F *)f2->Get("hbin3mult45");
  TH1F *hsummultall56 = (TH1F *)f0->Get("hbin3mult56");
  TH1F *hsummultbkg56 = (TH1F *)f1->Get("hbin3mult56");
  TH1F *hsummultsig56 = (TH1F *)f2->Get("hbin3mult56");

  TH1F *hsummultall345 = (TH1F *)f0->Get("hbin3mult45");
  TH1F *hsummultbkg345 = (TH1F *)f1->Get("hbin3mult45");
  TH1F *hsummultsig345 = (TH1F *)f2->Get("hbin3mult45");
  TH1F *hsummultall445 = (TH1F *)f0->Get("hbin4mult45");
  TH1F *hsummultbkg445 = (TH1F *)f1->Get("hbin4mult45");
  TH1F *hsummultsig445 = (TH1F *)f2->Get("hbin4mult45");

  TH1F *hsummultall356 = (TH1F *)f0->Get("hbin3mult56");
  TH1F *hsummultbkg356 = (TH1F *)f1->Get("hbin3mult56");
  TH1F *hsummultsig356 = (TH1F *)f2->Get("hbin3mult56");
  TH1F *hsummultall456 = (TH1F *)f0->Get("hbin4mult56");
  TH1F *hsummultbkg456 = (TH1F *)f1->Get("hbin4mult56");
  TH1F *hsummultsig456 = (TH1F *)f2->Get("hbin4mult56");


  if(mode == 2)
    {
      hall = (TH1F *)hsumall145->Clone("hall"); hall->Add(hsumall245); hall->Add(hsumall345); hall->Add(hsumall445);
      hbkg = (TH1F *)hsumbkg145->Clone("hbkg"); hbkg->Add(hsumbkg245); hbkg->Add(hsumbkg345); hbkg->Add(hsumbkg445);
      hsig = (TH1F *)hsumsig145->Clone("hsig"); hsig->Add(hsumsig245); hsig->Add(hsumsig345); hsig->Add(hsumsig445);
      fitmax = 50; submax = 35;
      if(year == 20172018)
        {
          fitmax = 80;
          submax = 60;
        }
    }
  if(mode == 3)
    {
      hall = (TH1F *)hsumall256->Clone("hall"); hall->Add(hsumall356); hall->Add(hsumall456);
      hbkg = (TH1F *)hsumbkg256->Clone("hbkg"); hbkg->Add(hsumbkg356); hbkg->Add(hsumbkg456);
      hsig = (TH1F *)hsumsig256->Clone("hsig"); hsig->Add(hsumsig356); hsig->Add(hsumsig456);
      fitmax = 50; submax = 35;
      if(year == 20172018)
        {
          fitmax = 80;
          submax = 60;
        }
    }
  if(mode == 4)
    {
      if(pot == 220)
	{
	  hall = (TH1F *)hsumall145->Clone("hall"); 
	  hall->Add(hsumall245); hall->Add(hsumall345); hall->Add(hsumall445); hall->Add(hsumall256); hall->Add(hsumall356); hall->Add(hsumall456);
	  hbkg = (TH1F *)hsumbkg145->Clone("hbkg"); 
	  hbkg->Add(hsumbkg245); hbkg->Add(hsumbkg345); hbkg->Add(hsumbkg445); hbkg->Add(hsumbkg256); hbkg->Add(hsumbkg356); hbkg->Add(hsumbkg456);
	  hsig = (TH1F *)hsumsig145->Clone("hsig"); 
	  hsig->Add(hsumsig245); hsig->Add(hsumsig345); hsig->Add(hsumsig445); hsig->Add(hsumsig256); hsig->Add(hsumsig356); hsig->Add(hsumsig456);
	  fitmax = 50; submax = 30;
	}
      if(pot == 210) // Start from larger xi acceptance
	{
	  hall = (TH1F *)hsumall245->Clone("hall");
	  hall->Add(hsumall345); hall->Add(hsumall445); hall->Add(hsumall356); hall->Add(hsumall456);
	  hbkg = (TH1F *)hsumbkg245->Clone("hbkg");
          hbkg->Add(hsumbkg345); hbkg->Add(hsumbkg445); hbkg->Add(hsumbkg356); hbkg->Add(hsumbkg456);
          hsig = (TH1F *)hsumsig245->Clone("hsig");
          hsig->Add(hsumsig345); hsig->Add(hsumsig445); hsig->Add(hsumsig356); hsig->Add(hsumsig456);
	  fitmax = 50; submax = 30;
	}
      if(year == 20172018)
	{
	  fitmax = 80;
	  submax = 60;
	}
    }
  if(mode == 5)
    {
      hall = (TH1F *)hsummultall345->Clone("hall"); hall->Add(hsummultall356); hall->Add(hsummultall445); hall->Add(hsummultall456);
      hbkg = (TH1F *)hsummultbkg345->Clone("hbkg"); hbkg->Add(hsummultbkg356); hbkg->Add(hsummultbkg445); hbkg->Add(hsummultbkg456);
      hsig = (TH1F *)hsummultsig345->Clone("hsig"); hsig->Add(hsummultsig356); hsig->Add(hsummultsig445); hsig->Add(hsummultsig456);
      fitmax = 50; submax = 40;
      if(year == 20172018)
        {
          fitmax = 80;
          submax = 60;
        }

    }
  if(mode == 6)
    {
      hall = (TH1F *)hsumall345->Clone("hall"); hall->Add(hsumall445);
      hbkg = (TH1F *)hsumbkg345->Clone("hbkg"); hbkg->Add(hsumbkg445);
      hsig = (TH1F *)hsumsig345->Clone("hsig"); hsig->Add(hsumsig445);
    }
  if(mode == 7)
    {
      hall = (TH1F *)hsumall356->Clone("hall"); hall->Add(hsumall456);
      hbkg = (TH1F *)hsumbkg356->Clone("hbkg"); hbkg->Add(hsumbkg456);
      hsig = (TH1F *)hsumsig356->Clone("hsig"); hsig->Add(hsumsig456);
    }
  if(mode == 8)
    {
      hall = (TH1F *)hsummultall345->Clone("hall"); hall->Add(hsummultall445);
      hbkg = (TH1F *)hsummultbkg345->Clone("hbkg"); hbkg->Add(hsummultbkg445);
      hsig = (TH1F *)hsummultsig345->Clone("hsig"); hsig->Add(hsummultsig445);

      fitmax = 35; submax = 20;
      if(year == 20172018)
        {
          fitmax = 80;
          submax = 60;
        }
    }
  if(mode == 9)
    {
      hall = (TH1F *)hsummultall356->Clone("hall"); hall->Add(hsummultall456);
      hbkg = (TH1F *)hsummultbkg356->Clone("hbkg"); hbkg->Add(hsummultbkg456);
      hsig = (TH1F *)hsummultsig356->Clone("hsig"); hsig->Add(hsummultsig456);

      fitmax = 35; submax = 20;
      if(year == 20172018)
        {
          fitmax = 80;
          submax = 60;
        }
    }

  if(year == 20172018)
    {
      fitmax = 50;
      submax = 50;
    }


  TString thetitle = "";
  if(hist == "hbin145" || hist == "hbin156")
    thetitle = "0.02 < #xi(#mu#mu) < 0.03";
  if(hist == "hbin245" || hist == "hbin256")
    thetitle = "0.03 < #xi(#mu#mu) < 0.04";
  if(hist == "hbin345" || hist == "hbin356")
    thetitle = "0.04 < #xi(#mu#mu) < 0.06";
    //    thetitle = "#xi(#mu#mu) > 0.04";
  if(hist == "hbin445" || hist == "hbin456")
    thetitle = "#xi(#mu#mu) > 0.06";
  if(mode == 2)
    thetitle = "#xi(#mu#mu) > 0.02";
  if(mode == 3)
    thetitle = "#xi(#mu#mu) > 0.03";
  if(mode == 4)
    {
      if(pot == 220)
	thetitle = "#xi(#mu#mu) > 0.02, 0.03";
      if(pot == 210)
	thetitle = "#xi(#mu#mu) > 0.03, 0.04";
    }
  if(mode == 5)
    thetitle = "#xi(#mu#mu) > 0.04";
  if(mode == 6 || mode == 7)
    thetitle = "#xi(#mu#mu) > 0.04";
  if(mode == 8)
    thetitle = "#xi(#mu#mu) > 0.04";
  if(mode == 9)
    thetitle = "#xi(#mu#mu) > 0.04";
  if(hist == "hbin3mult45" || hist == "hbin3mult56")
    thetitle = "0.04 < #xi(#mu#mu) < 0.06";
  if(hist == "hbin4mult45" || hist == "hbin4mult56")
    thetitle = "#xi(#mu#mu) > 0.06";
  
  if(year == 2017)
    thetitle += ", 2017, L=36.9 fb^{-1}";
  if(year == 2018)
    thetitle += ", 2018, L=52.2 fb^{-1}";
  if(year == 20172018)
    thetitle += ", 2017+2018, L=89.1 fb^{-1}";

  hall->Rebin(rebinfact); hbkg->Rebin(rebinfact); hsig->Rebin(rebinfact);
  hall->GetXaxis()->SetRangeUser(-3,1); hbkg->GetXaxis()->SetRangeUser(-3,1); hsig->GetXaxis()->SetRangeUser(-3,1);

  hall->SetBinErrorOption(TH1::kPoisson);
  //  hbkg->SetBinErrorOption(TH1::kPoisson);
  //  hsig->SetBinErrorOption(TH1::kPoisson);

  TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array                                                            
  mc->Add(hbkg);
  mc->Add(hsig);
  
  TFractionFitter *fit = new TFractionFitter(hall,mc);
  TVirtualFitter* vFit = (TVirtualFitter* )fit->GetFitter();
  fit->Constrain(0,0.0,1.0);
  fit->Constrain(1,0.0,1.0);
  fit->SetRangeX(hall->FindBin(-3),hall->FindBin(1));
  
  //  vFit->SetParameter(0,"bkg",0.0,1.0,0.0,1000.0);
  //  vFit->SetParameter(1,"sig",0.15,1.0,0.10,1000.0);
  fit->Fit();
      
  TH1* result = fit->GetPlot();
  result->SetLineStyle(3);
  
  double theFrac[2], err[2];
  fit->GetResult(0, theFrac[0], err[0]);
  fit->GetResult(1, theFrac[1], err[1]);


  if(1)
    {
      Float_t bkgshapenorm = hbkg->Integral(hbkg->FindBin(-5),hbkg->FindBin(-0.5)) + hbkg->Integral(hbkg->FindBin(0.5),hbkg->FindBin(1));
      Float_t datnorm = hall->Integral(hall->FindBin(-5),hall->FindBin(-0.5)) + hall->Integral(hall->FindBin(0.5),hall->FindBin(1));
      Float_t allnorm = hall->Integral(hall->FindBin(-5),hall->FindBin(1));
      hbkg->Scale(datnorm/bkgshapenorm);
      Float_t scaledbkgnorm = hbkg->Integral(hbkg->FindBin(-5),hbkg->FindBin(1));
      Float_t signorm = allnorm - scaledbkgnorm;
      hsig->Scale(signorm / hsig->Integral(hsig->FindBin(-5),hsig->FindBin(1)));
    }

  //  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  //  c2->Divide(1,2);
  TCanvas *c2 = new TCanvas("c2","c2",1200,400);
  c2->Divide(2,1);
  c2->cd(1);
  hbkg->Scale(theFrac[0] * hall->GetSumOfWeights() / hbkg->GetSumOfWeights());
  hsig->Scale(theFrac[1] * hall->GetSumOfWeights() / hsig->GetSumOfWeights());

  /* All floating */
  TH1F *hsplusb = (TH1F *)hsig->Clone("hsplusb");
  
  hsplusb->Add(hbkg);  // single-inel + double-inel                                                                                   
                                    
  hsplusb->SetMaximum(fitmax);
  hsplusb->SetFillColor(0);
  hsplusb->SetLineColor(4);
  hsplusb->SetLineWidth(3);
  hsplusb->SetLineColor(4);
  //  hsplusb->SetFillColor(4); 
  //  hsplusb->SetTitle(thetitle);
  hsplusb->SetTitle(0);
  hsplusb->GetXaxis()->SetTitleSize(0.05);
  hsplusb->GetYaxis()->SetTitleOffset(0.9);
  hsplusb->SetXTitle("1 - #xi(p)/#xi(#mu#mu)");
  hsplusb->SetYTitle("Events");
  hsplusb->Draw("hist");
  hbkg->SetLineColor(2); hbkg->SetMarkerStyle(0);
  hbkg->SetFillColor(2); hbkg->Draw("histsame");
  hall->SetLineColor(1); hall->SetLineWidth(3); hall->SetMarkerColor(1); hall->SetMarkerStyle(20); hall->SetLineWidth(3);

  TH1F *hsplusberr = (TH1F *)hsplusb->Clone("hsplusberr");
  hsplusberr->SetFillColor(13);
  hsplusberr->SetMarkerSize(0);
  hsplusberr->Draw("E2same");
  // JH - PPD comments
  hsplusb->Draw("histsame");
  hall->Draw("E0same");

  hsplusberr->SetFillColor(13);


  TLegend *le1 = new TLegend(0.15,0.3,0.65,0.85);
  le1->AddEntry(hall,"Data");
  le1->AddEntry(hbkg,"Background shape");
  le1->AddEntry(hsplusberr,"Signal shape (pp#rightarrowp#mu#mup MC) + background");
  le1->SetLineWidth(0);
  le1->Draw("same");

  TLatex latex;                                                                                                                                                                           
  float l = gPad->GetLeftMargin();                                                                                                                                                        
  float t = gPad->GetTopMargin();                                                                                                                                                         
  latex.SetTextAlign(11);                                                                                                                                                                 
  latex.SetTextSize(0.6*t);                                                                                                                                                               
  //  TString cmsText     = "#font[61]{CMS-TOTEM} #scale[0.76]{#font[52]{Preliminary}}";
  //  latex.DrawLatexNDC(l+0.5,1-t-0.7*t,cmsText);                                                                                                 
  TString cmsText     = "#font[61]{CMS-TOTEM}";
  latex.DrawLatexNDC(l+0.5,1-t-0.7*t,cmsText);                                      

  TLatex latexb;
  float lb = gPad->GetLeftMargin();
  float tb = gPad->GetTopMargin();
  latexb.SetTextAlign(11);
  latexb.SetTextSize(0.6*t);
  TString cmsTextb     = "#font[42]{92.3 fb^{-1} (13 TeV)}";
  //  latexb.DrawLatexNDC(lb+0.6,1-tb+0.2*tb,cmsTextb);
  latexb.DrawLatexNDC(lb+0.5,1-tb+0.2*tb,cmsTextb);

  /*
  TLatex latex;
  l = gPad->GetLeftMargin();
  t = gPad->GetTopMargin();
  latex.SetTextAlign(11);
  latex.SetTextSize(0.6*t);
  TString cmsText     = "#font[61]{CMS} #scale[0.76]{#font[52]{Preliminary 2017+2018, Sector 45 (#xi > 0.02) + Sector 56 (#xi > 0.03), 220m FAR}}";
  latex.DrawLatexNDC(l,1-t+0.2*t,cmsText);
  */

  c2->cd(2);
  TH1F *hsub = (TH1F *)hall->Clone("hsub");
  hbkg->Sumw2();
  hsub->Add(hbkg,-1);
  //  hsub->Add(hsig,-1);hsub->Add(hbkg,-1);
  hsub->SetBinErrorOption(TH1::kPoisson);

  hsig->GetXaxis()->SetRangeUser(-1,1);
  //  hsig->GetXaxis()->SetRangeUser(-3,1);
  hsig->SetMinimum(-3); hsig->SetMaximum(hsub->GetMaximum()*2.0); //hsig->SetMaximum(submax);
  hsig->SetStats(0);
  hsig->GetXaxis()->SetTitleSize(0.05);
  hsig->GetYaxis()->SetTitleOffset(0.9);
  hsig->SetXTitle("1 - #xi(p)/#xi(#mu#mu)");
  hsig->SetFillColor(0);
  hsig->SetLineColor(4);
  hsig->SetLineWidth(3);
  hsig->SetLineColor(4);
  //  hsig->SetFillColor(4);
  hsig->SetTitle(0);
  hsig->SetYTitle("Data - B");
  hsig->Draw("hist");

  TH1F *hsigerr = (TH1F *)hsig->Clone("hsigerr");
  hsigerr->SetFillColor(13); 
  hsigerr->SetMarkerSize(0);
  hsigerr->SetMarkerColor(0);
  hsigerr->Draw("E2same");
  hsig->Draw("histsame");

  // JH - testing Sept 8, 2020. ROOT doesn't handle asymmetric errors in the histogram subtraction. 
  // So test a "by hand" quadrature sum of the upper and lower errors instead. 
  //  hsub->Draw("E0same");
  TGraphAsymmErrors *hsubasym = new TGraphAsymmErrors(hsub);
  for(Int_t z = 0; z < hall->GetNbinsX(); z++)
    {
      float xup1 = hall->GetBinErrorUp(z+1);
      float xup2 = hbkg->GetBinError(z+1);
      float xupdiff = TMath::Sqrt(xup1*xup1 + xup2*xup2);
      float xlo1 = hall->GetBinErrorLow(z+1);
      float xlo2 = hbkg->GetBinError(z+1);
      float xlodiff = TMath::Sqrt(xlo1*xlo1 + xlo2*xlo2);
      hsubasym->SetPointEYhigh(z,xupdiff);
      hsubasym->SetPointEYlow(z,xlodiff);
      hsubasym->SetPointEXhigh(z,0);
      hsubasym->SetPointEXlow(z,0);
    }
  hsubasym->SetMarkerStyle(20);
  hsubasym->SetMarkerColor(0);
  hsubasym->Draw("E0same");

  Float_t nbkgpeak = hbkg->Integral(hbkg->FindBin(-0.15),hbkg->FindBin(0.15));

  TLegend *le2 = new TLegend(0.15,0.3,0.45,0.85);
  le2->AddEntry(hsub,"Background-subtracted data");
  le2->AddEntry(hsigerr,"Signal shape (pp#rightarrowp#mu#mup MC)");
  le2->SetLineWidth(0);
  le2->Draw("same");

  TLatex latex2;
  float l2 = gPad->GetLeftMargin();
  float t2 = gPad->GetTopMargin();
  latex2.SetTextAlign(11);
  latex2.SetTextSize(0.6*t);
  //  TString cmsText2     = "#font[61]{CMS-TOTEM} #scale[0.76]{#font[52]{Preliminary}}";
  TString cmsText2     = "#font[61]{CMS-TOTEM}";
  //  latex2.DrawLatexNDC(l2+0.5,1-t2-0.7*t,cmsText2);
  latex2.DrawLatexNDC(l2+0.5,1-t2-0.7*t,cmsText2);

  TLatex latex2b;
  float l2b = gPad->GetLeftMargin();
  float t2b = gPad->GetTopMargin();
  latex2b.SetTextAlign(11);
  latex2b.SetTextSize(0.6*t);
  TString cmsText2b     = "#font[42]{92.3 fb^{-1} (13 TeV)}";
  //  latex2b.DrawLatexNDC(l2b+0.6,1-t2b+0.2*t2b,cmsText2b);
  latex2b.DrawLatexNDC(l2b+0.5,1-t2b+0.2*t2b,cmsText2b);

  // Systematic error bands!
  if(1)
    {
      TFile *f2s = TFile::Open("../reminiaodtests/MoreDimuons20172018MCallanglesWithMultiLegacyFinal_SystJanShiftDown.root");
      TFile *f3s = TFile::Open("../reminiaodtests/MoreDimuons20172018MCallanglesWithMultiLegacyFinal_SystJanShiftUp.root");
      //      TH1F *h2 = (TH1F *)f2->Get("hressummult");
      //      TH1F *h3 = (TH1F *)f3->Get("hressummult");

      TH1F *hsummultsigshiftd345 = (TH1F *)f2s->Get("hbin3mult45");
      TH1F *hsummultsigshiftd445 = (TH1F *)f2s->Get("hbin4mult45");
      TH1F *hsummultsigshiftd356 = (TH1F *)f2s->Get("hbin3mult56");
      TH1F *hsummultsigshiftd456 = (TH1F *)f2s->Get("hbin4mult56");
      TH1F *hsummultsigshiftu345 = (TH1F *)f3s->Get("hbin3mult45");
      TH1F *hsummultsigshiftu445 = (TH1F *)f3s->Get("hbin4mult45");
      TH1F *hsummultsigshiftu356 = (TH1F *)f3s->Get("hbin3mult56");
      TH1F *hsummultsigshiftu456 = (TH1F *)f3s->Get("hbin4mult56");

      TH1F *h2 = (TH1F *)hsummultsigshiftd345->Clone("h2"); h2->Add(hsummultsigshiftd445);
      h2->Add(hsummultsigshiftd356); h2->Add(hsummultsigshiftd456);

      TH1F *h3 = (TH1F *)hsummultsigshiftu345->Clone("h3"); h3->Add(hsummultsigshiftu445);
      h3->Add(hsummultsigshiftu356); h3->Add(hsummultsigshiftu456);

      h2->Rebin(rebinfact); h3->Rebin(rebinfact);

      h2->Scale(theFrac[1] * hall->GetSumOfWeights() / h2->GetSumOfWeights());
      h3->Scale(theFrac[1] * hall->GetSumOfWeights() / h3->GetSumOfWeights());


      Int_t nbins = hsig->GetNbinsX();
      Float_t binup[nbins];
      Float_t bindown[nbins];
      Float_t binup2[nbins];
      Float_t bindown2[nbins];
      Float_t binx[nbins];
      Float_t biny[nbins];
      Float_t binxup[nbins];
      Float_t binxdown[nbins];


      TH1F *hband = (TH1F *)hsig->Clone("hband");

      for(Int_t n = 1; n < hsig->GetNbinsX(); n++)
	{
	  Float_t ndef = hsig->GetBinContent(n);
	  Float_t nup = h3->GetBinContent(n);
	  Float_t ndown = h2->GetBinContent(n);

	  Float_t highshift = TMath::Max(nup,ndown);
	  Float_t lowshift = TMath::Min(nup,ndown);

	  //	  Float_t high = TMath::Max(highshift,ndef);
	  //	  Float_t low = TMath::Min(lowshift,ndef);
	  Float_t high = TMath::Max(ndown,ndef);
	  Float_t low = TMath::Min(ndown,ndef);
	  Float_t high2 = TMath::Max(nup,ndef);
	  Float_t low2 = TMath::Min(nup,ndef);

	  binup[n-1] = high-ndef;
	  bindown[n-1] = ndef-low;
          binup2[n-1] = high2-ndef;
          bindown2[n-1] = ndef-low2;
	  
	  binxdown[n-1] = (hsig->GetBinWidth(n)/2.0);
	  binxup[n-1] = (hsig->GetBinWidth(n)/2.0);
	  binx[n-1] = hsig->GetBinCenter(n);
	  biny[n-1] = hsig->GetBinContent(n);
	}
      TGraphAsymmErrors *gr = new TGraphAsymmErrors(nbins,binx,biny,binxdown,binxup,bindown,binup);
      TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(nbins,binx,biny,binxdown,binxup,bindown2,binup2);

      gr->SetMarkerStyle(0); gr->SetFillColor(13); gr->SetFillStyle(3004); gr->SetLineStyle(0); //gr->Draw("PE2same");
      gr2->SetMarkerStyle(0); gr2->SetFillColor(13); gr2->SetFillStyle(3005); gr2->SetLineStyle(0); //gr2->Draw("PE2same");
      gr->Draw("E2same");
      gr2->Draw("E2same");

      hsigerr->Draw("E2same");
      hsig->Draw("histsame");
      //      hsub->Draw("E0same");
      hsubasym->Draw("E0same");

      le2->AddEntry(gr,"#xi systematic uncertainty (shift down)");
      le2->AddEntry(gr2,"#xi systematic uncertainty (shift up)");
      le2->Draw("same");

      //      hband->SetFillColor(13); hband->SetMarkerSize(0); hband->Draw("E2same");
      //      h1->SetLineWidth(3); h1->SetLineColor(4); h1->Draw("histsame");
    }

  // -0.4, 0.4
  TString fitopt = "EMVS";
  if(minimizer == 1)
    fitopt = "LEMVS";

  //  TF1 *gausfitdat = new TF1("gausfitdat","gaus(0)",-fitrange,fitrange);
  //  gausfitdat->SetLineColor(2); gausfitdat->SetLineWidth(3); gausfitdat->SetLineStyle(2);
  //  //  TFitResultPtr fitdat1 = hsub->Fit("gausfitdat","LEMVS","",-fitrange,fitrange);
  //  TFitResultPtr fitdat1 = hsub->Fit("gausfitdat",fitopt,"",-fitrange,fitrange);
  //  TF1 *gausfitmc = new TF1("gausfitmc","gaus(0)",-fitrange,fitrange);
  //  gausfitmc->SetLineColor(3); gausfitmc->SetLineWidth(3);
  //  TFitResultPtr fitmc1 = hsigerr->Fit("gausfitmc","LEMVS","",-fitrange,fitrange);
  //  std::cout << "# " << hist << ", " << rebinfact << ", " << year << ", " << mode << ", " << pot << ", " << fitrange << std::endl;
  //  std::cout << "# Data mean and error, MC mean and error , Data sigma and error, MC sigma and error" << std::endl;
  //  std::cout << gausfitdat->GetParameter(1) << " "
  //	    << gausfitdat->GetParError(1) << " "
  //            << gausfitmc->GetParameter(1) << " "
  //            << gausfitmc->GetParError(1) << " "
  //	    << gausfitdat->GetParameter(2) << " "
  //	    << gausfitdat->GetParError(2) << " "
  //            << gausfitmc->GetParameter(2) << " "
  //            << gausfitmc->GetParError(2) << std::endl;

  //  std::cout << "N(data, gaussian integral) = " << gausfitdat->Integral(-0.5,0.5) << std::endl;
  std::cout << "N(sig) = " << hsig->GetSumOfWeights() << std::endl;
  std::cout << "N(bkg, -0.15-0.15) = " << nbkgpeak << std::endl;
  std::cout << "RMS(sig) = " << hsig->GetRMS() << std::endl;
  std::cout << "chi2_BC = " << fit->GetChisquare() << std::endl;
  std::cout << "NDF = " << fit->GetNDF() << std::endl;
  std::cout << "p(chi2) = " << fit->GetProb() << std::endl;

  //  std::cout << "# chi2_BC/NDF (data) = " << gausfitdat->GetChisquare() << "/" 
  //	    << gausfitdat->GetNDF() << std::endl;
  //  std::cout << "# p(chi2, data) = " << gausfitdat->GetProb() << std::endl;
  //  std::cout << "# chi2_BC/NDF (MC) = " << gausfitmc->GetChisquare() << "/" 
  //            << gausfitmc->GetNDF() << std::endl;
  //  std::cout << "# p(chi2, MC) = " << gausfitmc->GetProb() << std::endl;

  //  Float_t nbkgpeak = hbkg->Integral(hbkg->FindBin(-0.1),hbkg->FindBin(0.1));                                                                                                      
  Float_t ndatapeak = hall->Integral(hbkg->FindBin(-0.1),hbkg->FindBin(0.1));
  Float_t nbkgpeak2 = hbkg->Integral(hbkg->FindBin(-0.05),hbkg->FindBin(0.2));
  Float_t ndatapeak2 = hall->Integral(hbkg->FindBin(-0.05),hbkg->FindBin(0.2));

  std::cout << "Data in 45 signal region = " << ndatapeak << ", bkg in signal region = " << nbkgpeak << std::endl;
  std::cout << "Data in 56 signal region = " << ndatapeak2 << ", bkg in signal region = " << nbkgpeak2 << std::endl;

  if(mode == 4)
    CMSTOTEM_lumi((TPad*)c2->GetPad(2),0,0,"single-RP");
  if(mode == 5)
    {
      //      CMSTOTEM_lumi((TPad*)c2->GetPad(1),0,0,"multi-RP");
      //      CMSTOTEM_lumi((TPad*)c2->GetPad(2),0,0,"multi-RP");
    }
}
