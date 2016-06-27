Double_t fitfunc(Double_t *v, Double_t *par)
{
  Double_t arg = 0;
  arg = v[0];

  Double_t fitval = par[0] + TMath::Exp(-arg * par[1]);
  //  Double_t fitval = TMath::Exp(-arg*par[1] + par[0]);
  //  Double_t fitval = par[0] + par[1]*arg + par[2]*arg*arg;
  //  Double_t fitval = par[0] + (par[1]*arg) + (par[2]*arg*arg) + (par[3]*arg*arg);

  return fitval;
}

TF1 * BackgroundTowerFitter(Int_t seed = 32423)
{
  TFile *f1 = new TFile("dymumu.610.anal.root");
  TFile *f2 = new TFile("dymumu.1040.anal.root");
  TFile *f3 = new TFile("dymumu.40.anal.root");
  TFile *f4 = new TFile("gamgammumu.lpair.anal.root");
  TFile *f5 = new TFile("gamgammumu.lpairinelastic.anal.root");
  TFile *f6 = new TFile("upsilonmumu.starlight.anal.root");
  TFile *f7 = new TFile("upsilonmumu.starlight2s.anal.root");
  TFile *f8 = new TFile("upsilonmumu.starlight3s.anal.root");
  TFile *f9 = new TFile("mumu.stew.histos.root");
  TFile *f10 = new TFile("mumu.gumbo.histos.root");
  TFile *f11 = new TFile("mumu.chowder.histos.root");

  TTree *t1 = f1->Get("ntp1");
  TTree *t2 = f2->Get("ntp1");
  TTree *t3 = f3->Get("ntp1");
  TTree *t4 = f4->Get("ntp1");
  TTree *t5 = f5->Get("ntp1");
  TTree *t6 = f6->Get("ntp1");
  TTree *t7 = f7->Get("ntp1");
  TTree *t8 = f8->Get("ntp1");

  Int_t ntrack1, ntrack2, ntrack3, ntrack4, ntrack5, ntrack6, ntrack7, ntrack8;
  Int_t ntower1, ntower2, ntower3, ntower4, ntower5, ntower6, ntower7, ntower8;
  Double_t mmmass1, mmmass2, mmmass3, mmmass4, mmmass5, mmmass6, mmmass7, mmmass8;
  Double_t mmdpt1[2], mmdpt2[2], mmdpt3[2], mmdpt4[2];
  Double_t mmdpt5[2], mmdpt6[2], mmdpt7[2], mmdpt8[2], mmdpt9[2];
  Double_t mmdphi1, mmdphi2, mmdphi3, mmdphi4, mmdphi5, mmdphi6, mmdphi7, mmdphi8, mmdphi9;
  Int_t hitinzdc5, hitincastor5;

  t1->SetBranchAddress("nTrackCand",&ntrack1);
  t1->SetBranchAddress("nExtraCaloTowersE5",&ntower1);
  t2->SetBranchAddress("nTrackCand",&ntrack2);
  t2->SetBranchAddress("nExtraCaloTowersE5",&ntower2);
  t3->SetBranchAddress("nTrackCand",&ntrack3);
  t3->SetBranchAddress("nExtraCaloTowersE5",&ntower3);
  t4->SetBranchAddress("nTrackCand",&ntrack4);
  t4->SetBranchAddress("nExtraCaloTowersE5",&ntower4);
  t5->SetBranchAddress("nTrackCand",&ntrack5);
  t5->SetBranchAddress("nExtraCaloTowersE5",&ntower5);
  t6->SetBranchAddress("nTrackCand",&ntrack6);
  t6->SetBranchAddress("nExtraCaloTowersE5",&ntower6);
  t7->SetBranchAddress("nTrackCand",&ntrack7);
  t7->SetBranchAddress("nExtraCaloTowersE5",&ntower7);
  t8->SetBranchAddress("nTrackCand",&ntrack8);
  t8->SetBranchAddress("nExtraCaloTowersE5",&ntower8);

  t1->SetBranchAddress("MuonCand_pt",&mmdpt1);
  t1->SetBranchAddress("MuMu_dphi",&mmdphi1);
  t2->SetBranchAddress("MuonCand_pt",&mmdpt2);
  t2->SetBranchAddress("MuMu_dphi",&mmdphi2);
  t3->SetBranchAddress("MuonCand_pt",&mmdpt3);
  t3->SetBranchAddress("MuMu_dphi",&mmdphi3);
  t4->SetBranchAddress("MuonCand_pt",&mmdpt4);
  t4->SetBranchAddress("MuMu_dphi",&mmdphi4);
  t5->SetBranchAddress("MuonCand_pt",&mmdpt5);
  t5->SetBranchAddress("MuMu_dphi",&mmdphi5);
  t6->SetBranchAddress("MuonCand_pt",&mmdpt6);
  t6->SetBranchAddress("MuMu_dphi",&mmdphi6);
  t7->SetBranchAddress("MuonCand_pt",&mmdpt7);
  t7->SetBranchAddress("MuMu_dphi",&mmdphi7);
  t8->SetBranchAddress("MuonCand_pt",&mmdpt8);
  t8->SetBranchAddress("MuMu_dphi",&mmdphi8);


  t1->SetBranchAddress("MuMu_mass",&mmmass1);
  t2->SetBranchAddress("MuMu_mass",&mmmass2);
  t3->SetBranchAddress("MuMu_mass",&mmmass3);
  t4->SetBranchAddress("MuMu_mass",&mmmass4);
  t5->SetBranchAddress("MuMu_mass",&mmmass5);
  t6->SetBranchAddress("MuMu_mass",&mmmass6);
  t7->SetBranchAddress("MuMu_mass",&mmmass7);
  t8->SetBranchAddress("MuMu_mass",&mmmass8);

  t5->SetBranchAddress("HitInZDC",&hitinzdc5);
  t5->SetBranchAddress("HitInCastor",&hitincastor5);

  TH1F *h1 = new TH1F("h1","h1",25,0,25);
  TH1F *h2 = new TH1F("h2","h2",25,0,25);
  TH1F *h3 = new TH1F("h3","h3",25,0,25);
  TH1F *h4 = new TH1F("h4","h4",25,0,25);
  TH1F *h5 = new TH1F("h5","h5",25,0,25);
  TH1F *h6 = new TH1F("h6","h6",25,0,25);
  TH1F *h7 = new TH1F("h7","h7",25,0,25);
  TH1F *h8 = new TH1F("h8","h8",25,0,25);
  TH1F *h9 = new TH1F("h9","h9",25,0,25);

  TH1F *hm1 = new TH1F("hm1","hm1",20,8,12);

//   t1->Draw("nExtraCaloTowersE5 >> h1","nTrackCand < 5","e");
//   t2->Draw("nExtraCaloTowersE5 >> h2","nTrackCand < 5","e");
//   t3->Draw("nExtraCaloTowersE5 >> h3","nTrackCand < 5","e");
//   t4->Draw("nExtraCaloTowersE5 >> h4","nTrackCand < 5","e");
//   t5->Draw("nExtraCaloTowersE5 >> h5","nTrackCand < 5","e");
//   t6->Draw("nExtraCaloTowersE5 >> h6","nTrackCand < 5","e");
//   t7->Draw("nExtraCaloTowersE5 >> h7","nTrackCand < 5","e");
//   t8->Draw("nExtraCaloTowersE5 >> h8","nTrackCand < 5","e");

  TRandom *rnd = new TRandom(seed);

  double lumi = 10.0;
  double xsec1 = 0.5 * 37820.0 * lumi * (t1->GetEntries()/96500.0);
  double xsec2 = 0.5 * 5951.0 * lumi  * (t2->GetEntries()/ 38800.0);
  double xsec3 = 0.5 * 1797.0 * lumi * (t3->GetEntries()/ 22000.0);
  double xsec4 = 74.65 * lumi * (t4->GetEntries()/ 100000.0);
  double xsec5 = 76.2 * lumi * (t5->GetEntries()/ 20000.0);
  double xsec6 = 39.0 * lumi * (t6->GetEntries()/ 10000.0);
  double xsec7 = 13.0 * lumi * (t7->GetEntries()/ 3359.0);
  double xsec8 = 10.0 * lumi * (t8->GetEntries()/ 2636.0);

  h9 = (TH1F *)f9->Get("hcalofit");
  h10 = (TH1F *)f10->Get("hcalofit");
  h11 = (TH1F *)f11->Get("hcalofit");
  //  h9 = (TH1F *)f9->Get("hnextracaloe");
  //  h10 = (TH1F *)f10->Get("hnextracaloe");
  h9->Scale(0.1);h10->Scale(0.1);h11->Scale(0.1);

  Int_t count1=0; Int_t count2=0; Int_t count3=0; Int_t count4=0; Int_t count5=0; Int_t count6=0; Int_t count7=0; Int_t count8=0;
  Int_t peakingcount = 0; Int_t gamgamcount = 0; Int_t upsiloncount = 0;
  Int_t sigregioncount = 0;
  Int_t trkcut = 4;
  Double_t dphicut = 2.9;
  Double_t dptcut = 2.0;

  for(Int_t i1 = 0;i1 < t1->GetEntries();i1++)
    {
      t1->GetEntry(i1);
      if(rnd->Uniform() > 0.5)
        continue;

      if(ntrack1 < trkcut && count1 < xsec1 && fabs(mmdpt1[0]-mmdpt1[1]) < dptcut && mmdphi1 > dphicut)
	{
	  h1->Fill(ntower1);
	  if(ntower1 < 5)
	    {
	      sigregioncount++;
	      hm1->Fill(mmmass1);
	    }
	}
      count1++;
    }
  for(Int_t i2 = 0;i2 < t2->GetEntries();i2++)
    {
      t2->GetEntry(i2);
      if(rnd->Uniform() > 0.5)
	continue;

      if(ntrack2 < trkcut && count2 < xsec2  && fabs(mmdpt2[0]-mmdpt2[1]) < dptcut && mmdphi2 > dphicut)
	{
	  h2->Fill(ntower2);
          if(ntower2 < 5)
	    {
	      sigregioncount++;
              hm1->Fill(mmmass2);
            }
        }

      count2++;
    }
  for(Int_t i3 = 0;i3 < t3->GetEntries();i3++)
    {
      t3->GetEntry(i3);
      if(rnd->Uniform() > 0.5)
        continue;

      if(ntrack3 < trkcut && count3 < xsec3  && fabs(mmdpt3[0]-mmdpt3[1]) < dptcut && mmdphi3 > dphicut)
	{
	  h3->Fill(ntower3);
          if(ntower3 < 5)
	    {
	      sigregioncount++;
              hm1->Fill(mmmass3);
            }
        }

      count3++;
    }
  for(Int_t i4 = 0;i4 < t4->GetEntries();i4++)
    {
      t4->GetEntry(i4);
      if(rnd->Uniform() > 0.5)
        continue;

      if(ntrack4 < trkcut && count4 < xsec4  && fabs(mmdpt4[0]-mmdpt4[1]) < dptcut && mmdphi4 > dphicut)
	{
	  h4->Fill(ntower4);
          if(ntower4 < 5)
	    {
	      sigregioncount++;
              hm1->Fill(mmmass4);
            }
        }

      count4++;
    }
  for(Int_t i5 = 0;i5 < t5->GetEntries();i5++)
    {
      t5->GetEntry(i5);
      if(rnd->Uniform() > 0.5)
        continue;

      if(ntrack5 < trkcut && count5 < xsec5 && hitinzdc5 == 0  && fabs(mmdpt5[0]-mmdpt5[1]) < dptcut && mmdphi5 > dphicut)
	{
	  h5->Fill(ntower5);
	  if(ntower5 < 5)
	    {
              hm1->Fill(mmmass5);

	      peakingcount++;
	      sigregioncount++;
	    }
	}
      count5++;
    }
  for(Int_t i6 = 0;i6 < t6->GetEntries();i6++)
    {
      t6->GetEntry(i6);
      if(rnd->Uniform() > 0.5)
        continue;

      if(ntrack6 < trkcut && count6 < xsec6  && fabs(mmdpt6[0]-mmdpt6[1]) < dptcut && mmdphi6 > dphicut)
	{
	  h6->Fill(ntower6);
          if(ntower6 < 5)
	    {
              hm1->Fill(mmmass6);
	      sigregioncount++;
	    }
        }

      count6++;
    }
  for(Int_t i7 = 0;i7 < t7->GetEntries();i7++)
    {
      t7->GetEntry(i7);
      if(rnd->Uniform() > 0.5)
        continue;

      if(ntrack7 < trkcut && count7 < xsec7  && fabs(mmdpt7[0]-mmdpt7[1]) < dptcut && mmdphi7 > dphicut)
	{
	  h7->Fill(ntower7);
          if(ntower7 < 5)
	    {
	      sigregioncount++;
              hm1->Fill(mmmass7);
            }
        }

      count7++;
    }
  for(Int_t i8 = 0;i8 < t8->GetEntries();i8++)
    {
      t8->GetEntry(i8);
      if(rnd->Uniform() > 0.5)
        continue;

      //      cout << ntower8 << ", " << ntrack8 << ", Xsec = " << xsec8 << ", i8 = " << i8 << endl;
      if(ntrack8 < trkcut && count8 < xsec8  && fabs(mmdpt8[0]-mmdpt8[1]) < dptcut && mmdphi8 > dphicut)
	{
	  h8->Fill(ntower8);
          if(ntower8 < 5)
	    {
	      sigregioncount++;
              hm1->Fill(mmmass8);
            }
	}

      count8++;
    }


  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h4->Sumw2();
  h5->Sumw2();
  h6->Sumw2();
  h7->Sumw2();
  h8->Sumw2();
  //  h9->Sumw2();


//   h1->Scale(xsec1);
//   h2->Scale(xsec2);
//   h3->Scale(xsec3);
//   h4->Scale(xsec4);
//   h5->Scale(xsec5);
//   h6->Scale(xsec6);
//   h7->Scale(xsec7);
//   h8->Scale(xsec8);

//  h1->Add(h2); h1->Add(h3); h1->Add(h4); h1->Add(h5); h1->Add(h6); h1->Add(h7); h1->Add(h8); h1->Add(h9);
  h10->Add(h4); h10->Add(h5); h10->Add(h6); h10->Add(h7); h10->Add(h8); h10->Add(h9); h10->Add(h11);
  h4->Add(h5); h4->Add(h6); h4->Add(h7); h4->Add(h8);
  h4->SetFillColor(4);

  h6->Add(h5);h6->Add(h7);h6->Add(h8);

  h6->SetFillColor(3);
  h7->SetFillColor(6);
  h8->SetFillColor(1);
  h9->SetFillColor(2);

  h5->SetFillColor(5);

  //  h1->SetLineColor(2);
  h10->SetFillColor(0);
  h10->SetXTitle("N(Extra towers)");
  h10->SetYTitle("Events");
  h10->Draw("e");
  
  TF1 *func = new TF1("thefit",fitfunc,-10,10,2);
  func->SetParameters(0.05,-1.0);
  func->SetParNames("a","b");  
  h10->Fit("thefit","","",5,25);
  Double_t par1 = func->GetParameter("a");
  Double_t par2 = func->GetParameter("b");
  //  Double_t par3 = func->GetParameter("c");
  //  Double_t par4 = func->GetParameter("d");
  TString strpar1 = Form("%f", par1);
  TString strpar2 = Form("%f", par2);
  //  TString strpar3 = Form("%f", par3);
  //  TString strpar4 = Form("%f", par4);

  //  cout << strpar1 << ", " << strpar2 << endl;

  TF1 *func2 = new TF1("func1",strpar1 + " + TMath::Exp(-x * " + strpar2 + ")",0,5);
  //  TF1 *func2 = new TF1("func1",strpar1 + " + (x * " + strpar2 + ") + (x * x * " + strpar3 + ")",0,5);
  //  TF1 *func2 = new TF1("func1",strpar1 + " + (x*" + strpar2 + ") + (x*x*" + strpar3 + ") + (x*x*x*" + strpar4 + ")",0,5); 
  func->SetLineColor(2);
  func->SetLineWidth(3);
  func2->SetLineColor(2);
  func2->SetLineStyle(2);
  func2->SetLineWidth(3);
  func2->Draw("same");
  cout << "Integral = " << func->Integral(0,5) << endl;
  //  func->IntegralError(0,5);
  cout << "Inelastic = " << peakingcount << endl;
  cout << "Events in signal region = " << sigregioncount << endl;
  //  f1->Close();
  //  f2->Close();
  //  f3->Close();
  //  f4->Close();
  //  f5->Close();
  //  f6->Close();
  //  f7->Close();
  //  f8->Close();

  h10->Draw("e");
  h4->Draw("histsame");
  h6->Draw("histsame");
  h5->Draw("histsame");
  //  h6->Draw("histsame");
  //  h7->Draw("histsame");
  //  h8->Draw("histsame");

  h2->SetMarkerStyle(20);
  h2->SetLineWidth(3);
  h10->SetMarkerStyle(20);
  h10->SetLineWidth(3);

  h10->Draw("esame");
  func->Draw("same");
  func2->Draw("same");

  TLegend *l1 = new TLegend(0.5,0.5,0.8,0.8);
  l1->AddEntry(h10,"\"Data\" (Signal plus bkg.)");
  l1->AddEntry(h4,"Elastic #gamma #gamma #rightarrow #mu^{+} #mu^{-}");
  l1->AddEntry(h6,"Elastic #Upsilon #rightarrow #mu^{+} #mu^{-}");
  l1->AddEntry(h5,"Singly inelastic #gamma #gamma #rightarrow #mu^{+} #mu^{-}");
  l1->AddEntry(func,"Fit");
  l1->Draw("same");

  ofstream ofs("fitfileunweightedtowers.txt");
  for(Int_t z = 1;z < 26;z++)
    {
      Int_t nt = h1->GetBinContent(z);
      for(Int_t y = 0;y < nt;y++)
	ofs << z-1 << endl;
    }

  //  return (func->Integral(0,5));
  //  hm1->Draw("e");
  return func;
}

void FitLooper()
{
  TH1F *h1 = new TH1F("h1","h1",10,0,10);
  Double_t nbkg = 0.0;

  TRandom *sd = new TRandom(65539);
  for(Int_t k = 0;k < 100;k++)
    {
      Int_t theseed = sd->Uniform(1,1000000);
      
      nbkg = UnweightedTowerFitter(theseed);
      h1->Fill(nbkg);
    }
  h1->Draw("hist");
}
