/// Final file for the evaluation of the luminosity from the fit
/// X.Rouby 10/7/8, after PhD predefense. 
/// This is a modified version of Luminosity2.C
/// Based on the exclusive dilepton analysis CMS PAS DIF-07-001
/// Only dimuon analysis 
///
/// Usage : Luminosity(100,4,0.072,0.016,5,2.4,1) 
///         after compilation in ROOT/CINT:	 .L Luminosity3.C++

#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"
#include "TPaveText.h"
#include "FuncDef.h"
#include "TF1.h"
#include <cmath>

double xsection = 74.65;	// pb
double efftrigger = 0.10;	//
double effselsignal =0;		//

// Return histogrammed quantities for mu+mu- samples
TH1F *GetMuMuHist(const double lumi = 100., const int plotvar = 6, const int physsample = 1, const float ptcut=3., const float etacut=2.4, const int maxevent=-1, const bool nozdcveto=true, const double phicut = 2.9) // phicut can change if restricted fit range !
// plotvar 
///     3) : dpht2
///     4) : dphi2
///     5) : mll2
///     9) : sumpt
{
  // analysis cuts
  Double_t ecalocut  = 5.0;	// 5.0
  //Double_t etcalocut = 0.2; 	// 0.2
  Double_t deltarcut = 0.3;	// 0.3
  Double_t ntrackcut = 3;	// 3
  Double_t ncalocut  = 5;	// 5
  Double_t dptcut    = 2.0;	// 2.0
  Double_t dphicut   = phicut;	// 2.9 // dphicut = 3.0913272;
  Double_t upsloveto = 9; 	// 9
  Double_t upshiveto = 11; 	// 10.8; // 11

  // initialization
  TString st = "";
  Double_t xsec = 0.0;		
  Int_t linecolor = 1;
  Int_t linewidth = 0;
  Int_t fillcolor = kWhite; //15;
  Int_t fillstyle = 0;

  switch(physsample) {
  case 1:
    xsec = xsection * lumi / 100000.0;
    st = "gamgammumu.lpair.anal.root";
    linecolor = 1; linewidth = 3; fillcolor = 15; fillstyle = 1; 
    break;
  case 2:
    xsec = 0.5 * 37820.0 * lumi / 96500.0;
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
    xsec = 39.0 * lumi / 10000.0;	// cross-section(pb) * 100 (pb-1) / samplesize
    st = "upsilonmumu.starlight.anal.root";
    fillcolor = 3;
    break;
  case 7:
    xsec = 13.0 * lumi / 3359.0;	// cross-section(pb) * 100 (pb-1) / samplesize
    st = "upsilonmumu.starlight2s.anal.root";
    fillcolor = 6;
    break;
  case 8:
    xsec = 10.0 * lumi / 2636.0;	// cross-section(pb) * 100 (pb-1) / samplesize
    st = "upsilonmumu.starlight3s.anal.root";
    fillcolor = 1;
    break;
  case 9:
    xsec = 76.2 * lumi / 20000.0 ;	// cross-section(pb) * 100 (pb-1) / samplesize
    st = "gamgammumu.lpairinelasticcteq.anal.root";
    fillcolor = 5;
    break;
  default:
    break;
  }

  if(maxevent>0) xsec = 1.;

  TFile *f1 = new TFile(st);
  TTree *tr1 = (TTree*) f1->Get("ntp1");

  TH1F *hmll2 = new TH1F("mll2","mll2",200,0,200); // 400 or 50
  TH1F *hdphi2 = new TH1F("hdphi2","hdphi",80,-0.1,0.1); // (phi1 - phi2)/pi
  TH1F *hdpt2 = new TH1F("hdpt2","hdpt",20,-2,2); //"hdpt2","hdpt",50,-5,5);
  TH1F *sumpt = new TH1F("sumpt","sumpt",9,0,180);

  Double_t mumumass, mumudphi, mupt[2], mueta[2], muphi[2];
  Int_t nmuons, njets, ntracks, ncalo, zdchit, hitincastor;
  Double_t deltara, deltarb, caloe[1000], caloet[1000], caloeta[1000], calophi[1000];

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

  if(physsample == 9) {
    tr1->SetBranchAddress("HitInZDC",&zdchit);
    tr1->SetBranchAddress("HitInCastor",&hitincastor);
  } else { zdchit =0; hitincastor =0; }

  int ntotsig = tr1->GetEntries(), step =0;
  int nkinepass = 0, nisocalo = 0, ndphipass =0 , ndptpass=0 , ncalopass=0, ntrackpass = 0, nupspass = 0;

  TRandom rnd1;
  bool firsttime =false; // for the cout -- should be true



  /////// EVENT LOOP ///////////
  for(int i = 0;i < ntotsig;i++) {
    tr1->GetEntry(i);
    step++ ; // counter for the evaluation of the efficiencies

      // Stupid trick - ignore half of events with Castor hits to simulate presence of only one Castor
      double castorrnd1 = rnd1.Uniform();
      if((hitincastor >= 1) && (castorrnd1 < 0.5)) hitincastor = 0;


      // computes the number of calo towers
      nisocalo = 0;
      for(Int_t x = 0;x < ncalo;x++) {
	  deltara = sqrt((caloeta[x]-mueta[0])*(caloeta[x]-mueta[0]) + (calophi[x]-muphi[0])*(calophi[x]-muphi[0]));
	  deltarb = sqrt((caloeta[x]-mueta[1])*(caloeta[x]-mueta[1]) + (calophi[x]-muphi[1])*(calophi[x]-muphi[1]));
	  if(deltara > deltarcut && deltarb > deltarcut && caloe[x] > ecalocut) nisocalo++; // caloe or caloet
      }


  // ** 1) pt and eta cuts
      if( (fabs(mueta[0]) <= etacut) && (fabs(mueta[1]) <= etacut) && (mupt[0] > ptcut) && (mupt[1] > ptcut)   ) {
	nkinepass++;
  // ** 2) Delta phi cut
	if(mumudphi > dphicut) { if (firsttime) {cout << "dphicut = " << dphicut <<endl; firsttime=false;}
	  ndphipass++;
  // ** 3) Delta pt cut
	  if(fabs(mupt[0]-mupt[1]) < dptcut) {	//ptcut
	      ndptpass++;
  // ** 4) calo exclusivity cut
		if(nisocalo < ncalocut) {  
		  ncalopass++;
  // ** 5) tracking exclusivity cut
		  if(ntracks < ntrackcut) {  //trackcut
		      ntrackpass++;
  // ** 6) Upsilon veto 
	 	      if(mumumass < upsloveto || mumumass > upshiveto) { 
  // ** 7) CASTOR/ZDC veto
		        if ( nozdcveto ||  ((zdchit == 0)  && (hitincastor == 0)) ) {
			  nupspass++;
			  hmll2->Fill(mumumass);
			  sumpt->Fill(mupt[0]+mupt[1]);
			  hdpt2->Fill((mupt[0]-mupt[1]));
			  hdphi2->Fill( (fabs(varphi(muphi[0])) -fabs(varphi(muphi[1])))/pi   );
			} // zdc/castor veto
		     } // upsilon veto
		  } // trackcut
	  	} //calocut
	  } // ptcut
	} // phicut
      } //eta+pt cut
  if (nupspass==maxevent) break;
  } ////////// loop on events



	cout << "Efficiency/rejection for file: " << st << endl;
/*	cout << "\tnkine  eff = " 	<< (double)nkinepass/step << endl;
	cout << "\td(phi) eff = " 	<< (double)ndphipass/step << endl;
	cout << "\td(pt) eff = " 	<< (double)ndptpass/step << endl;
	cout << "\tcalo excl eff = " 	<< (double)ncalopass/step << endl;
	cout << "\ttrack excl eff = " << (double)ntrackpass/step << endl;
	cout << "\tUpsilon veto eff = " << (double)nupspass/step << endl;*/
	if(physsample==1 )  { effselsignal = (double) nupspass/step; cout << "effselsignal in GetMuMu = " << effselsignal << endl;}
//	cout << "Efficiency = " 	<< (double) nupspass/step <<endl;
  
  cout << "N(sig) = " << (double)nupspass * xsec << " +- " << sqrt(nupspass) * xsec << " after " << lumi << " pb-1 ";
  if (maxevent<0) cout << "(errors from MC sample size)\n\n"; else cout << "(errors from Poisson statistics)\n\n";

  TCanvas * c2 = new TCanvas("c2","GetMuMu");
  c2->cd();
  if(plotvar == 3){ // dimuon Delta pt after selection
	hdpt2->SetFillColor(fillcolor);
	hdpt2->SetFillStyle(fillstyle);
	hdpt2->SetMaximum(6000);
	hdpt2->SetStats(0);
	hdpt2->SetTitle(0);
	hdpt2->SetXTitle("#Delta p_{T} (GeV)");
	char title[500]; sprintf(title,"Events/(%g GeV)",hdpt2->GetBinWidth(2));
	hdpt2->SetYTitle(title);
	hdpt2->Draw("hist");
	hdpt2->Sumw2();
	//  if(size<0)hdpt2->Scale(xsec);
	hdpt2->SetLineColor(linecolor);
	hdpt2->SetLineWidth(linewidth);
	hdpt2->Draw("hist");
	return(hdpt2);
  } else if(plotvar == 4){ //dphi2 = |phi1| - |phi2| after selection
	hdphi2->SetFillColor(fillcolor);
	hdphi2->SetFillStyle(fillstyle);
	hdphi2->SetStats(0);
	hdphi2->SetTitle(0);
	hdphi2->SetXTitle("(|#phi_{1}| - |#phi_{2}|)/#pi");
	hdphi2->GetXaxis()->SetNdivisions(509);
	char text[500]; sprintf(text,"Events/(%g)",hdphi2->GetBinWidth(1));
	hdphi2->GetYaxis()->SetTitle(text);
	hdphi2->Sumw2();
	hdphi2->Scale(xsec);
	hdphi2->SetLineColor(linecolor);
	hdphi2->SetLineWidth(linewidth);
	hdphi2->Draw("hist");
	return(hdphi2);
  } else if(plotvar == 5){  // dimuon invariant mass after analysis
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
	hmll2->SetXTitle("M(#mu#mu)");
	hmll2->Draw("hist");
	return(hmll2);
  } else if(plotvar == 9){ // pt sum
	sumpt->Sumw2();
	sumpt->Scale(xsec);
	sumpt->SetMinimum(0.001);
	sumpt->SetMaximum(20.0*lumi);
	sumpt->SetTitle(0);
	sumpt->SetStats(0);
	sumpt->SetXTitle("#Sigma pt (#mu #mu)");
	sumpt->SetFillColor(fillcolor);
	sumpt->SetFillStyle(fillstyle);
	sumpt->SetLineColor(linecolor);
	sumpt->SetLineWidth(linewidth);
	sumpt->Draw("hist");
	return(sumpt);
  } else {
	cout << "'plotvar' not equal to 3, 4, 5 or 9; no draw !\n";
	return new TH1F("temp","",200,0,200);
  }
} // end of GetMuMuHist



///////////////////////// Plots for µ+µ-, after selection and 100pb-1
// Luminosity(4,0.1,0.016,5,2.4) ; nozdcveto ==true if you do NOT want any veto with ZDC/CASTOR
void Luminosity(const double lumi = 100., const int plottype=-1, const float bin = 2., const float binsig=0.016, const float ptcut=3., const float etacut=2.4, const bool nozdcveto=true) {
  // plottype ==
  // 	3 = Delta pt
  // 	4 = Delta phi  Luminosity(4,0,0.1,0.016,5,2.4) 
  // 	5 = invariant mass

  if(plottype < 0) {
	cout << "Luminosity(plotvar,...)\n";
	cout << "\tplotvar= 3 for Delta pt\n";
	cout << "\t         4 for Delta phi\n";
	cout << "\t         5 for dimuon invariant mass M\n";
	cout << "\t         9 for dimuon pt\n";
	return;
  }

  setTDRStyle();

  TH1F  * hsignalmc,	// signalmc = all elastic events, from MC sample 
	* hinelmc, 	// inelmc = all inelastic events, from MC sample
	* hsignal,	// signal = expected elastic events after a given luminosity
	* hinel, 	// same
	* hsignalrestricted, // same but for a narrow range of delta phi
	* hinelrestricted;   // same but for a narrow range of delta phi


  /// *********** FIRST PART : estimate the expected number of events after the selection
  ///
  /// use hsignalmc and hinelmc to do so, and fills the nsigtrue and nbkgtrue variables

  cout << endl << "**************\n";
  cout << "Histograms with the number of events corresponding to the MC sample size\n";
  hsignalmc = (TH1F*) (GetMuMuHist(lumi,plottype,1,ptcut,etacut,-1,nozdcveto))->Clone("hsignalmc");
  hinelmc   = (TH1F*) (GetMuMuHist(lumi,plottype,9,ptcut,etacut,-1,nozdcveto))->Clone("hinelmc");
  float nsigmc = hsignalmc->GetEntries();
  float nsigtrue = nsigmc * xsection * lumi / 100000.0;
  cout << "Histo signal     " << nsigmc  << "\t rescaled: " << nsigtrue << endl;
  // TH1F * hups1S, * hups2S, * hups3S, * hups;
  // hups1S  = (TH1F*) (GetMuMuHist(lumi,plottype,6,ptcut,etacut,-1,nozdcveto))->Clone("hups1S"); 
  // hups2S  = (TH1F*) (GetMuMuHist(lumi,plottype,7,ptcut,etacut,-1,nozdcveto))->Clone("hups2S");
  // hups3S  = (TH1F*) (GetMuMuHist(lumi,plottype,8,ptcut,etacut,-1,nozdcveto))->Clone("hups3S");
  // hups    = (TH1F*) hups1S->Clone("hups");
  // hups->Add(hups2S);
  // hups->Add(hups3S);    */
  float nbkgmc = hinelmc->GetEntries();
  float nbkgtrue = nbkgmc * 76.2 * lumi / 20000.0;
  cout << "Histo inelastic  " << nbkgmc << "\t rescaled: " << nbkgtrue << endl;


  /// *********** SECOND PART : create the histograms with the true number of events, after selection, assuming a given luminosity
  ///
  /// this creates hsignal and hinel
 
  cout << endl << "**************\n";
  cout << "Histograms with the true number of events and not more\n";
  hsignal = (TH1F*) ( GetMuMuHist(lumi,plottype,1,ptcut,etacut,(int)round(nsigtrue),nozdcveto) )->Clone("hsignal");
  hinel   = (TH1F*) ( GetMuMuHist(lumi,plottype,9,ptcut,etacut,(int)round(nbkgtrue),nozdcveto) )->Clone("hinel");

	// ********* Histo cosmetics 
	  hsignalmc->SetFillColor(14); 
	  hsignalmc->SetLineColor(14); 
	  hsignalmc->SetFillStyle(1001);

	// hups->SetLineColor(kBlack);
	// hups->SetLineWidth(2);

	  hinel->SetLineColor(kBlack); 
	  hinel->SetLineWidth(3);

  TLegend * leg =0; 
  TLatex * mytex = new TLatex(1,1,"");
  THStack * s = new THStack("s","");
  TCanvas * c1 = new TCanvas("c1","");
  c1->Divide(2,2); 
  c1->cd(1);
  s->Add(hsignalmc,"hist");
  //s->Add(hups,   "hist");
  s->Add(hinel,  "hist");
  s->Draw("nostack");
  char temp[500];

/*
  ///// JUST TO CHECK THE function used for estimating the backgrounds
  cout << "\n\n***************************** INELASTICS\n";
  TCanvas *cinel = new TCanvas("cinel","inelastics");
  cinel->cd();
  hinel->Draw("hist");
  hinel->SetMinimum(0);
  TF1 * b1 = new TF1("b1","[1]/(pow(x,2)+pow([0],2)/4.)",-bin,bin);
  hinel->Fit("b1","RQ");
  hinel->Fit("b1","RQ");
  hinel->Fit("b1","RIQ");
  hinel->Fit("b1","RI");
  cout << "hinel *****\n";
  cout << b1->GetExpFormula() << endl;
  cout << "chi^2 / nDF = " << b1->GetChisquare() <<  "\t" << b1->GetNDF() << "\t = " << b1->GetChisquare()/b1->GetNDF() << "\t [" << 1-sqrt(2)/sqrt(b1->GetNDF()) << ";" << 1+sqrt(2)/sqrt(b1->GetNDF()) << "]\n";
  cout << endl;

  TPaveText * pave1 = new TPaveText(0.16,0.68,0.47,0.92,"NDC");
  pave1->SetBorderSize(0);
  pave1->SetFillColor(0);
  pave1->SetTextFont(42);
  pave1->AddText("#frac{A_{1}}{x^{2} + #Gamma_{1}^{2}/4}");
  pave1->AddText("Inelastics");
  pave1->AddText("CMSSW 1 6 0");
  pave1->Draw();
  TPaveText * pave1a = new TPaveText(0.62,0.79,0.91,0.92,"NDC");
  pave1a->SetBorderSize(0);
  pave1a->SetFillColor(0);
  if (plottype==4) {
    sprintf(temp,"#Gamma_{1}=(%.3f #pm %.3f)#times 10^{-3}",1000*b1->GetParameter(0),1000*b1->GetParError(0));
    pave1a->AddText(temp);
    sprintf(temp,"A_{1}=(%.3f #pm %.3f)#times 10^{-3}",1000*b1->GetParameter(1),1000*b1->GetParError(1));
    pave1a->AddText(temp);
    sprintf(temp,"#chi^{2}/n=%.3f/%d=%.3f",b1->GetChisquare(),b1->GetNDF(),b1->GetChisquare()/b1->GetNDF());
    pave1a->AddText(temp);
  } else { 
    sprintf(temp,"#Gamma_{1}=%.3f #pm %.3f",b1->GetParameter(0),b1->GetParError(0));
    pave1a->AddText(temp);
    sprintf(temp,"A_{1}=%.3f #pm %.3f",b1->GetParameter(1),b1->GetParError(1));
    pave1a->AddText(temp);
    sprintf(temp,"#chi^{2}/n=%.3f/%d=%.3f",b1->GetChisquare(),b1->GetNDF(),b1->GetChisquare()/b1->GetNDF());
    pave1a->AddText(temp);
  }
  pave1a->Draw();
*/

  /// *********** THIRD PART : fits the signal MC sample
  ///
  /// use hsignalmc and fits it with the "b2" function in the narrow [-binsig,binsig] range

  cout << "\n\n******************************* ELASTICS\n";
  TCanvas *csignal = new TCanvas("csignal","elastics");
  csignal->cd();
  hsignalmc->Draw("hist");
  TF1 * b2 = new TF1("b2","[1]/(pow(x,2)+pow([0],2)/4.) + [3]/(pow(x,2)+pow([2],2)/4.)",-binsig,binsig);

  hsignalmc->Fit("b2","RQ");
  hsignalmc->Fit("b2","RQ");
  if (binsig < 0.03) {
    hsignalmc->Fit("b2","RQ");
    hsignalmc->Fit("b2","RIQ");
    hsignalmc->Fit("b2","RIEQ");
    hsignalmc->Fit("b2","RIE");
  } else hsignalmc->Fit("b2","RE");
  cout << "hsignalmc *****\n";
  cout << b2->GetExpFormula() << endl;
  cout << "chi^2 / nDF = " << b2->GetChisquare() <<  "\t" << b2->GetNDF() 
       << "\t = " << b2->GetChisquare()/b2->GetNDF() << "\t [" 
       << 1-sqrt(2)/sqrt(b2->GetNDF()) << ";" << 1+sqrt(2)/sqrt(b2->GetNDF()) << "]\n";
  cout << endl;

  TPaveText * pave2 = new TPaveText(0.16,0.68,0.47,0.92,"NDC");
  pave2->SetBorderSize(0);
  pave2->SetFillColor(0);
  pave2->AddText("#frac{A_{1}}{x^{2} + #Gamma_{1}^{2}/4} + #frac{A_{2}}{x^{2} + #Gamma_{2}^{2}/4}");
  pave2->AddText("Elastics");
  pave2->AddText("CMSSW 1 6 0");
  pave2->Draw();
  TPaveText * pave2a = new TPaveText(0.62,0.66,0.91,0.92,"NDC");
  pave2a->SetBorderSize(0);
  pave2a->SetFillColor(0);

  // parameters from the fit of MC signal
  float g1 = b2->GetParameter(0), g1e = b2->GetParError(0);
  float a1 = b2->GetParameter(1), a1e = b2->GetParError(1);
  float g2 = b2->GetParameter(2), g2e = b2->GetParError(2);
  float a2 = b2->GetParameter(3), a2e = b2->GetParError(3);
  float chi2 = b2->GetChisquare(), ndf=b2->GetNDF();

  if (plottype==4) {
    sprintf(temp,"#Gamma_{1}=(%.3f #pm %.3f)#times 10^{-3}",1000*g1,1000*g1e);
    pave2a->AddText(temp);
    sprintf(temp,"A_{1}=(%.3f #pm %.3f)#times 10^{-3}",1000*a1,1000*a1e);
    pave2a->AddText(temp);
    sprintf(temp,"#Gamma_{2}=(%.3f #pm %.3f)#times 10^{-3}",1000*g2,1000*g2e);
    pave2a->AddText(temp);
    sprintf(temp,"A_{2}=(%.3f #pm %.3f)#times 10^{-3}",1000*a2,1000*a2e);
    pave2a->AddText(temp);
    sprintf(temp,"#chi^{2}/n=%.3f/%d=%.3f",chi2,(int)ndf,chi2/ndf);
    pave2a->AddText(temp);
  } else {
    sprintf(temp,"#Gamma_{1}=%.3f #pm %.3f",g1,g1e);
    pave2a->AddText(temp);
    sprintf(temp,"A_{1}=%.3f #pm %.3f",a1,a1e);
    pave2a->AddText(temp);
    sprintf(temp,"#Gamma_{2}=%.3f #pm %.3f",g2,g2e);
    pave2a->AddText(temp);
    sprintf(temp,"A_{2}=%.3f #pm %.3f",a2,a2e);
    pave2a->AddText(temp);
    sprintf(temp,"#chi^{2}/n=%.3f/%d=%.3f",chi2,(int)ndf,chi2/ndf);
    pave2a->AddText(temp);
  }
  pave2a->Draw();

  if(strstr(gROOT->GetVersion(),"5.18")) {
           cout << "Integrale du signal dans [-0.02, 0.02]" << b2->Integral(-0.02,0.02)/hsignalmc->GetBinWidth(1) << " +- "
                << b2->IntegralError(-0.02,0.02)/hsignalmc->GetBinWidth(1) << endl; 
	   cout << "Valid only if fit range is [0.02;0.02]\n";
  } else {
           cout << "Integrale du signal dans  [-0.02, 0.02] " << b2->Integral(-0.02,0.02)/hsignalmc->GetBinWidth(1) << endl;
           cout << "Pour l'erreur sur l'intégrale, utiliser root 5.18\n";
  }


  /// *********** FOURTH PART : fits the observed data
  ///
  /// use hsum and fits it with the "b3" function over the wide [-bin,bin] range

  cout << "\n\n************************************ SUM \n";
  TCanvas *csum = new TCanvas("csum","elastics + inelastics");
  csum->cd();
  TH1F* hsum = (TH1F*) hsignal->Clone("hsum"); //hsignal and not hsignalmc !
    hsum->Add(hinel);
    hsum->GetXaxis()->SetNdivisions(507);
    hsum->GetYaxis()->SetNdivisions(507);
    hsum->Draw("hist");
    hsum->Rebin(2);

  // copies the result from the fit to the signal (a1,g1,a2,g2)
  char formuleb3[500];
  sprintf(formuleb3,"[2]*(%g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.)) + [1]/(pow(x,2)+pow([0],2)/4.)",a1,g1,a2,g2);
  // sprintf(formuleb3,"[2]*(%g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.)) + gaus",a1,g1,a2,g2);
  TF1 * b3 = new TF1("b3",formuleb3,-bin,bin); 

  b3->SetParLimits(0,0,20000); // gamma >0
  b3->SetParLimits(1,0,20000); // gamma >0
  b3->SetParLimits(2,0,20000);
  hsum->Fit("b3","RQ");
  hsum->Fit("b3","RQ");
  hsum->Fit("b3","RQ");
  hsum->Fit("b3","RIQ");
  hsum->Fit("b3","RIEQ");
  hsum->Fit("b3","RIE");
  cout << "hsum *****\n";
  cout << b3->GetExpFormula() << endl;
  cout << "chi^2 / nDF = " << b3->GetChisquare() <<  "\t" << b3->GetNDF() << "\t = " << b3->GetChisquare()/b3->GetNDF() << "\t [" << 1-sqrt(2)/sqrt(b3->GetNDF()) << ";" << 1+sqrt(2)/sqrt(b3->GetNDF()) << "]\n";
  cout << endl;

  TPaveText * pave3 = new TPaveText(0.16,0.68,0.47,0.92,"NDC");
  pave3->SetTextFont(42);
  pave3->SetTextAlign(12);
  pave3->SetBorderSize(0);
  pave3->SetFillColor(0);
  pave3->AddText("A S_{sig} + #frac{B}{x^{2} + #Gamma^{2}/4}");
  pave3->AddText("All events");
  pave3->AddText("CMSSW 1 6 0");
  pave3->Draw();
  TPaveText * pave3a = new TPaveText(0.62,0.66,0.91,0.92,"NDC");
  pave3a->SetBorderSize(0);
  pave3a->SetFillColor(0);

  float g3 = b3->GetParameter(0), g3e = b3->GetParError(0);
  float B3 = b3->GetParameter(1), B3e = b3->GetParError(1);
  float a3 = b3->GetParameter(2), a3e = b3->GetParError(2);
  if (plottype==4) {
    pave3a->SetTextFont(42);
    pave3a->SetTextAlign(12);
    sprintf(temp,"#Gamma=(%.3f #pm %.3f)#times 10^{-3}",1000*g3,1000*g3e);
    pave3a->AddText(temp);
    sprintf(temp,"A=(%.3f #pm %.3f)",a3,a3e);
    pave3a->AddText(temp);
    sprintf(temp,"B=(%.3f #pm %.3f)#times 10^{-3}",1000*B3,1000*B3e);
    pave3a->AddText(temp);
    sprintf(temp,"#chi^{2}/n=%.3f/%d=%.3f",b3->GetChisquare(),b3->GetNDF(),b3->GetChisquare()/b3->GetNDF());
    pave3a->AddText(temp);
  } else {
    sprintf(temp,"#Gamma_{1}=%.3f #pm %.3f",g3,g3e);
    pave3a->AddText(temp);
    sprintf(temp,"A_{1}=%.3f #pm %.3f",a3,a3e);
    pave3a->AddText(temp);
    sprintf(temp,"#Gamma_{2}=%.3f #pm %.3f",b3->GetParameter(2),b3->GetParError(2));
    pave3a->AddText(temp);
    sprintf(temp,"A_{2}=%.3f #pm %.3f",b3->GetParameter(3),b3->GetParError(3));
    pave3a->AddText(temp);
    sprintf(temp,"#chi^{2}/n=%.3f/%d=%.3f",b3->GetChisquare(),b3->GetNDF(),b3->GetChisquare()/b3->GetNDF());
    pave3a->AddText(temp);
  }
  pave3a->Draw();



// copy here the result of fit to b2
  char formuleb4[500]; // exact copy of b2 (= MC signal)
  sprintf(formuleb4,"%g*(%g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.))",a3,a1,g1,a2,g2);
  TF1 *b4= new TF1("b4",formuleb4,-binsig,binsig);
  cout << "b3= " << b3->GetExpFormula() << endl; // fit of hsum
  cout << "b4= " << b4->GetExpFormula() << endl; // fit of hsignal

  //char formuleb5[500]; // modified copy of b3, for systematic error studies
  //sprintf(formuleb5,"%g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.)",a1,g1,a2,g2,a3,0.095*g3);
  //TF1 *b5= new TF1("b5",formuleb5,-binsig,binsig);
  //cout << "b5= " << b5->GetExpFormula() << endl;

   int N_bkg = (int) round( (b3->Integral(-binsig,binsig) - b4->Integral(-binsig,binsig))/hsum->GetBinWidth(1) ); // is a double as a result of a fit
   if(strstr(gROOT->GetVersion(),"5.18")) {
	   cout << "Integrale de la somme dans [" << -binsig << "," << binsig << "] " << b3->Integral(-binsig,binsig)/hsum->GetBinWidth(1) << " +- " << b3->IntegralError(-binsig,binsig)/hsum->GetBinWidth(1) << endl; 
	   cout << "Integrale du bkg dans [" << -binsig << "," << binsig << "] " << N_bkg << " +- " << b3->IntegralError(-binsig,binsig)/hsum->GetBinWidth(1) << endl; 
	 //  cout << "Integrale du bkg dans [" << -binsig << "," << binsig << "] " << (b5->Integral(-binsig,binsig) - b4->Integral(-binsig,binsig))/hsum->GetBinWidth(1) << " for the modified distribution (syst. study)\n"; 
   } else {
	   cout << "Integrale de la somme dans [" << -binsig << "," << binsig << "] " << b3->Integral(-binsig,binsig)/hsum->GetBinWidth(1) << endl; 
	   cout << "Integrale du bkg dans [" << -binsig << "," << binsig << "] " << (b3->Integral(-binsig,binsig) - b4->Integral(-binsig,binsig))/hsum->GetBinWidth(1) << endl; 
	   cout << "Pour l'erreur sur l'intégrale, utiliser root 5.18\n";
   }


  ///////// ZOOMED version + fits
  // b1 = fit of inelastics only
  // b2 = fit of signal only (elastic) -- MC
  // b3 = fit of sum (elastic + inelastic) -- true number
  // b4 = exact copy of b2, for integration
  // b5 = modified copy of b3, for systematic error studies -- do not reuse !
  // b6 = exact copy of b3, for zoomed version
  // b7 = function of \Gamma and A only, for drawing the estimated bkg function
  char formuleb6[500]; // exact copy of b3, for the zoomed version
  sprintf(formuleb6,"%g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.) + %g/(pow(x,2)+pow(%g,2)/4.)",a1,g1,a2,g2,a3,g3);
  TF1 *b6= new TF1("b6",formuleb6,-binsig,binsig);
//  cout << "b6= " << b6->GetExpFormula() << endl;

  char formuleb7[500]; // copy of the bkg part of b3, for zoomed version
  sprintf(formuleb7,"%g/(pow(x,2)+pow(%g,2)/4.)",a3,g3);
  TF1 *b7= new TF1("b7",formuleb7,-binsig,binsig);
//  cout << "b7= " << b7->GetExpFormula() << endl;

  TCanvas *csumzoomed = new TCanvas("csumzoomed","elastics + inelastics (zoom)");
  csumzoomed->cd();
  TH1F* hsumzoomed = (TH1F*) hsignal->Clone("hsum");
    hsumzoomed->GetXaxis()->SetNdivisions(506);
    hsumzoomed->GetYaxis()->SetNdivisions(507);
    hsumzoomed->Draw();
    hsumzoomed->GetXaxis()->SetRangeUser(-binsig,binsig);
    hsumzoomed->SetMaximum(300);
    hsumzoomed->SetMinimum(1);
    b6->SetLineWidth(3);
    b6->SetLineColor(kBlack);
    b7->SetLineColor(kBlack);
    b4->SetLineColor(kBlack);
    b7->SetLineStyle(7);
    b6->Draw("same"); // sum
    b7->Draw("same"); // bkg
    b4->Draw("same"); // signal
    gPad->SetLogy();
    gPad->Update();

    TLegend * legzoom = new TLegend(0.66,0.66,0.92,0.92,"");
    legzoom->SetBorderSize(0);
    legzoom->SetFillStyle(0);
    legzoom->SetTextFont(42);
    legzoom->SetTextSizePixels(23);
    legzoom->AddEntry(hsumzoomed,"All events","p");
    legzoom->AddEntry(b4,"S_{fit}","l");
    legzoom->AddEntry(b7,"B_{fit}","l");
    legzoom->AddEntry(b6,"(B+S)_{fit}","l");
    legzoom->Draw();

    TPaveText * pavezoom = new TPaveText(0.16,0.85,0.47,0.90,"NDC");
    pavezoom->SetBorderSize(0);
    pavezoom->SetFillStyle(0);
    pavezoom->SetTextFont(42);
    pavezoom->AddText("CMSSW 1 6 0");
    pavezoom->Draw();







    /////// OLD plots 
    c1->cd(2);
    TH1F* sum = (TH1F*) hsignal->Clone("sum");
    sum->Add(hinel);
    sum->GetXaxis()->SetNdivisions(507);
    sum->GetYaxis()->SetNdivisions(507);


    /******* checking that the fit matches a posteriori the inelastics *********/
    c1->cd(3);
    TH1F* inelastics = (TH1F*) hinel->Clone("inelastics");
    inelastics->Draw();
    inelastics->GetXaxis()->SetNdivisions(507);
    inelastics->GetYaxis()->SetNdivisions(507);


	  float signal_true = hsignal->Integral();
//	  float signal_fit = 0; // fitsignal->Integral(-bin,bin)/sum->GetBinWidth(1);
//	  float error_signal_fit = fabs(signal_true - signal_fit)/signal_true;
//	  float inel_true = hinel->Integral();
//	  float inel_fit = 0; //fitinel->Integral(-bin,bin)/sum->GetBinWidth(1);
//	  float error_inel_fit = fabs(inel_true - inel_fit)/inel_true;
	  float sum_true = sum->Integral();
	  float sum_fit = 0; //fit3->Integral(-bin,bin)/sum->GetBinWidth(1);
//	  float error_sum_fit =  fabs(sum_true - sum_fit)/sum_true;

	//cout << "Signal : events from fit = " << signal_fit << " and N_true = " << signal_true << " (err= " << 100*error_signal_fit << "%)\n"; 
	//cout << "Inelastics : events from fit= " << inel_fit << " and N_true = " << inel_true << " (err= " << 100*error_inel_fit << "%)\n"; 
	//cout << "Sum = " << sum_fit << " and N_true = " << sum_true << " (err= " << 100*error_sum_fit << "%)\n";
//	  cout << "Signal true = " << signal_true << "; Sum true = " << sum_true << " in [ " << bin<< "," << bin << "]\n";
     if(plottype==3 || plottype==4) {
          if (plottype==3) cout << "Here: " << hinel->GetBinWidth(1) << " GeV per bin" << endl;
          else cout << "Here: " << hinel->GetBinWidth(1) << " per bin" << endl;
     } // if(plottype==3 || plottype==4) 





  // drawing old plots
  c1->cd(4);
  if (plottype==5) { 
	s->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
	s->SetMaximum(200); s->GetYaxis()->SetTitle("Events/(0.5 GeV)"); 
	s->GetXaxis()->SetRangeUser(0,50);
  	leg = new TLegend(0.40,0.65,0.9,0.9);
	mytex->DrawTextNDC(0.6,0.2,"CMSSW 1_6_0");
	//cout << "Invariant mass\n";
  } else if(plottype==3) {
	s->GetXaxis()->SetTitle("#Deltap_{T} (GeV)");
	s->SetMaximum(190);
	sum->SetMaximum(190);
	//s->GetXaxis()->SetRangeUser(-30,30);
	char title[500]; sprintf(title,"Events/(%g GeV)",sum->GetBinWidth(1));
	s->GetYaxis()->SetTitle(title);
	s->GetYaxis()->SetTitleOffset(1.05);
	s->GetXaxis()->SetNdivisions(507);
	s->GetYaxis()->SetNdivisions(507);
  	leg = new TLegend(0.01,0.35,0.95,0.9);

	c1->cd(1); gPad->Update();
	mytex->SetTextSizePixels(18);
	mytex->DrawTextNDC(0.18,0.86,"CMSSW 1_6_0");

        c1->cd(2); gPad->Update();
	mytex->SetTextSizePixels(18);
        mytex->DrawTextNDC(0.18,0.86,"N observed");
        char text[500];
        sprintf(text,"MC = %d",(int)round(sum_true));
        mytex->DrawTextNDC(0.18,0.76,text);
        sprintf(text,"Fit = %d",(int)round(sum_fit));
        mytex->DrawTextNDC(0.18,0.66,text);

        c1->cd(3); gPad->Update();
	mytex->SetTextSizePixels(18);
        mytex->DrawTextNDC(0.18,0.86,"Inelastics");
	//cout << "Delta pT\n";
  } else if(plottype==4) {
	s->GetXaxis()->SetTitle("( |#phi_{1}| - |#phi_{2}| ) / #pi");
	s->SetMaximum(200);
	s->GetXaxis()->SetRangeUser(-0.1,0.1);
	s->GetXaxis()->SetNdivisions(509);
	sum->SetMaximum(200);
	char title[500]; sprintf(title,"Events/(%g)",sum->GetBinWidth(1));
	s->GetYaxis()->SetTitle(title);
  	leg = new TLegend(0.01,0.35,0.95,0.9);
	c1->cd(1); gPad->Update();
	mytex->SetTextSizePixels(18);
	mytex->DrawTextNDC(0.18,0.86,"CMSSW 1_6_0");
	c1->cd(2); gPad->Update();
	mytex->SetTextSizePixels(18);
	mytex->DrawTextNDC(0.18,0.86,"N observed");
	char text[500]; 
	sprintf(text,"MC = %d",(int)round(sum_true));
	mytex->DrawTextNDC(0.18,0.76,text);
	sprintf(text,"Fit = %d",(int)round(sum_fit));
	mytex->DrawTextNDC(0.18,0.66,text);
	c1->cd(3); gPad->Update();
	mytex->DrawTextNDC(0.18,0.86,"Inelastics");
	mytex->SetTextSizePixels(18);
	//cout << "Delta phi\n";
  }
  s->GetYaxis()->SetTitleOffset(1.17);

  c1->cd(4);
  leg->AddEntry(hsignal,"#gamma #gamma #rightarrow #mu^{+} #mu^{-}","f");
  leg->AddEntry(hinel,"Singly inelastic","l");
//  char text[500]; sprintf(text,"Fit result (N_{inel}=%d)",(int)round(inel_fit));
  leg->Draw();
  leg->SetFillColor(0);
  leg->SetFillStyle(1);
  leg->SetBorderSize(0);
  delete c1;

  TCanvas *ctemp = new TCanvas("ctemp");
  ctemp->cd();
  // for the evaluation of the signal efficiency, with the narrow range of fit
  // i.e. put a smaller cut on D|phi|
  hsignalrestricted = (TH1F*) ( GetMuMuHist(lumi,plottype,1,ptcut,etacut,-1,nozdcveto, pi*(1-binsig)) )->Clone("hsignal");
  hinelrestricted   = (TH1F*) ( GetMuMuHist(lumi,plottype,9,ptcut,etacut,-1,nozdcveto, pi*(1-binsig)) )->Clone("hsignal");
  // this will modify the content of the effselsignal variable !
  // not clean coding, but this is life...
  int N_obs = (int) round(hsignalrestricted->Integral()) + round(hinelrestricted->Integral());
  delete ctemp;

  cout << "True luminosity = " << lumi << endl;
  cout << "Castor/Zdc veto = " ;  if(nozdcveto) cout << "off\n"; else cout << "on\n";
  cout << "N_obs = " << N_obs << "\t N_bkg = " << N_bkg << "\t eff_sel = " << effselsignal << "\t eff_trigger = " << efftrigger << "\t xsection = " << xsection << endl;
  double lumi_estimate = (double) (N_obs - N_bkg) / ( effselsignal * efftrigger  * xsection) ;
  cout << "Luminosity estimate = " << lumi_estimate << endl;
}

 /************ end of muon analysis *************/
