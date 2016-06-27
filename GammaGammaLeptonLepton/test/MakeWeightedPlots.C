void MakeWeightedPlots()
{
  //  TFile *f = new TFile("ee.stew.weighted.root");
  //  TTree *tr1 = f->Get("ntp1");
  TChain *tr1 = new TChain("ntp1");

  if(0){
    //   tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaEE-Gumbo/ee.recolevel_1.root");
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaEE-Stew/ee.recolevel_1.root");  
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaEE-Stew/ee.recolevel_2.root"); 
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaEE-Stew/ee.recolevel_3.root"); 
  }

  if(0){
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Gumbo/mumu.recolevel_1.root"); 
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Gumbo/mumu.recolevel_2.root");  
  }

  if(0){ 
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Chowder/mumu.recolevel_1.root");  
    tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Chowder/mumu.recolevel_2.root");   
  } 

  if(1){
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_1.root");
  cout << "Added 1" << endl;
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_3.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_4.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_5.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_6.root");
  cout << "Added 6" << endl;
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_7.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_8.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_9.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_10.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_11.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_12.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_13.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_14.root");
  cout << "Added 14" << endl;
  }
  if(1){
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_15.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_16.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_17.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_18.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_19.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_20.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_21.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_22.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_23.root");
  cout << "Added 23" << endl;
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_24.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_25.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_26.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_27.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_28.root");
  cout << "Added 28" << endl;
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_29.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_30.root");
  tr1->Add("rfio:/castor/cern.ch/user/j/jjhollar/gammagammaMuMu-Stew/mumu.recolevel_31.root");
  cout << "Added all" << endl;
  }

  TH1F *hm = new TH1F("hm","hm",50,0.0,200.0);
  TH1F *hm2 = new TH1F("hm2","hm2",50,0,200.0);
  TH1F *hmres = new TH1F("hmres","hmres",600,0,12.0);
  TH1F *hmz = new TH1F("hmz","hmz",100,50,150);
  TH1F *hmhi = new TH1F("hmhi","hmhi",100,150,1000);
  TH1F *hmsideband = new TH1F("hmsideband","hmsideband",50,0.0,200.0);
  TH1F *hdphi = new TH1F("hdphi","hdphi",192,0,3.2);
  TH1F *hdpt = new TH1F("hdpt","hdpt",240,-20,20);
  TH1F *hnextracaloe = new TH1F("hnextracaloe","hnextracaloe",50,0,50);
  TH1F *hnextracaloe4 = new TH1F("hnextracaloe4","hnextracaloe4",50,0,50);
  TH1F *hnextracaloe3 = new TH1F("hnextracaloe3","hnextracaloe3",50,0,50); 
  TH1F *hnextracaloe2 = new TH1F("hnextracaloe2","hnextracaloe2",50,0,50); 

  TH1F *hdphinoexcl = new TH1F("hdphinoexcl","hdphinoexcl",50,2.9,3.2);
  TH1F *hdptnoexcl = new TH1F("hdptnoexcl","hdptnoexcl",50,-2,2);

  TH1F *hnextracaloetpt2 = new TH1F("hnextracaloetpt2","hnextracaloetpt2",50,0,50); 
  TH1F *hnextracaloetpt5 = new TH1F("hnextracaloetpt5","hnextracaloetpt5",50,0,50); 
  TH1F *hnextracaloet1 = new TH1F("hnextracaloet1","hnextracaloet1",50,0,50);  
  TH1F *hnextracaloet2 = new TH1F("hnextracaloet2","hnextracaloet2",50,0,50);  

  TH1F *hntrack = new TH1F("hntrack","hntrack",50,0,50);
  TH1F *hcalofit = new TH1F("hcalofit","hcalofit",25,0,25);
  TH1F *htrackfit = new TH1F("htrackfit","htrackfit",50,0,50);

  Double_t mumumass;
  Double_t evweight;
  Int_t ntowers;

  Double_t mumumass;
  Double_t mumudphi;
  Double_t mupt[2];
  Double_t mueta[2];
  Double_t muphi[2];
  Int_t ntracks;
  Int_t ncalo, ncaloe4, ncaloe3, ncaloe2;
  Int_t ncaloetpt2, ncaloetpt5, ncaloet1, ncaloet2;

  tr1->SetBranchAddress("MuonCand_pt",&mupt);
  tr1->SetBranchAddress("MuMu_mass",&mumumass);
  tr1->SetBranchAddress("MuMu_dphi",&mumudphi);
  tr1->SetBranchAddress("nTrackCand",&ntracks);
  tr1->SetBranchAddress("nExtraCaloTowersE5",&ncalo);
  tr1->SetBranchAddress("nExtraCaloTowersE4",&ncaloe4);
  tr1->SetBranchAddress("nExtraCaloTowersE3",&ncaloe3); 
  tr1->SetBranchAddress("nExtraCaloTowersE2",&ncaloe2); 
  tr1->SetBranchAddress("nExtraCaloTowersEt0pt2",&ncaloetpt2); 
  tr1->SetBranchAddress("nExtraCaloTowersEt0pt5",&ncaloetpt5); 
  tr1->SetBranchAddress("nExtraCaloTowersEt1",&ncaloet1);  
  tr1->SetBranchAddress("nExtraCaloTowersEt2",&ncaloet2);  
  tr1->SetBranchAddress("evweight",&evweight);

  //  tr1->SetBranchAddress("EleCand_et",&mupt); 
  //  tr1->SetBranchAddress("ElEl_mass",&mumumass); 
  //  tr1->SetBranchAddress("ElEl_dphi",&mumudphi); 


  hm->Sumw2(); 
  hm2->Sumw2(); 
  hmres->Sumw2(); 
  hdphi->Sumw2(); 
  hdpt->Sumw2(); 
  hntrack->Sumw2(); 
  hnextracaloe->Sumw2(); 
  hnextracaloe4->Sumw2();  
  hnextracaloe3->Sumw2();  
  hnextracaloe2->Sumw2();  
  hnextracaloetpt2->Sumw2();
  hnextracaloetpt5->Sumw2(); 
  hnextracaloet1->Sumw2(); 
  hnextracaloet2->Sumw2(); 
  hcalofit->Sumw2(); 
  hmz->Sumw2(); 
  hmhi->Sumw2(); 

  Double_t dptcut = 5.0;
  Double_t dphicut = 2.7;

  Int_t nent = tr1->GetEntries();
  cout << "Nent = " << nent << endl;
  for(Int_t j = 0;j < nent;j++)
    {
      tr1->GetEntry(j);

      //      if(j%1000 == 0)
      cout << "Entry " << j << endl;

      hm->Fill(mumumass,evweight);
      hdpt->Fill(mupt[0]-mupt[1],evweight);
      hdphi->Fill(mumudphi,evweight);
      hmres->Fill(mumumass,evweight);
      hmhi->Fill(mumumass,evweight);
      hmz->Fill(mumumass,evweight);
      hntrack->Fill(ntracks,evweight);
      hnextracaloe->Fill(ncalo,evweight);
      hnextracaloe4->Fill(ncaloe4,evweight);
      hnextracaloe3->Fill(ncaloe3,evweight);
      hnextracaloe2->Fill(ncaloe2,evweight);

      hnextracaloetpt2->Fill(ncaloetpt2,evweight);
      hnextracaloetpt5->Fill(ncaloetpt5,evweight); 
      hnextracaloet1->Fill(ncaloet1,evweight); 
      hnextracaloet2->Fill(ncaloet2,evweight); 

      if((mupt[0]-mupt[1]) < dptcut && mumudphi > dphicut && ntracks < 3 && ncalo > 5)
	hmsideband->Fill(mumumass,evweight);

      if((mupt[0]-mupt[1]) < dptcut && mumudphi > dphicut)
	{
	  hdphinoexcl->Fill(mumudphi,evweight);
	  hdptnoexcl->Fill(mupt[0]-mupt[1],evweight);
	}

      if((mupt[0]-mupt[1]) < dptcut && mumudphi > dphicut && ntracks < 3)
	hcalofit->Fill(ncalo,evweight);
      
      if((mupt[0]-mupt[1]) < dptcut && mumudphi > dphicut && ncalo < 5 && ntracks < 3)
	hm2->Fill(mumumass,evweight);

      if((mupt[0]-mupt[1]) < dptcut && mumudphi > dphicut && ncalo < 5) 
        htrackfit->Fill(ntracks,evweight); 
    }

  
  hm->SetFillColor(10);
  hm2->SetFillColor(10);
  hmres->SetFillColor(10);
  hdphi->SetFillColor(10);
  hdpt->SetFillColor(10);
  hntrack->SetFillColor(10);
  hnextracaloe->SetFillColor(10);
  hnextracaloe4->SetFillColor(10); 
  hnextracaloe3->SetFillColor(10); 
  hnextracaloe2->SetFillColor(10); 
  hnextracaloetpt2->SetFillColor(10); 
  hnextracaloetpt5->SetFillColor(10);  
  hnextracaloet1->SetFillColor(10);  
  hnextracaloet2->SetFillColor(10);
  hcalofit->SetFillColor(8);
  hmz->SetFillColor(8);
  hmhi->SetFillColor(8);
  htrackfit->SetFillColor(8);
  hdphinoexcl->SetFillColor(10);
  hdptnoexcl->SetFillColor(10);
  hmsideband->SetFillColor(10);

  TFile *outfile = new TFile("/tmp/jjhollar/mumu.stew.histos1.root","RECREATE");

  hm->Write();
  hm2->Write();
  hmres->Write();
  hmz->Write();
  hmhi->Write();
  hdpt->Write();
  hdphi->Write();
  hntrack->Write();
  hnextracaloe->Write();
  hnextracaloe4->Write(); 
  hnextracaloe3->Write(); 
  hnextracaloe2->Write(); 
  hnextracaloetpt2->Write();
  hnextracaloetpt5->Write(); 
  hnextracaloet1->Write(); 
  hnextracaloet2->Write(); 
  hcalofit->Write();
  htrackfit->Write();
  hdphinoexcl->Write();
  hdptnoexcl->Write();
  hmsideband->Write();

  //  tr1->Write();  
  outfile->Write();

  if(0){
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,4);
  c1->cd(1);
  hm->Draw("hist");
  c1->cd(2);
  hmres->Draw("hist");
  c1->cd(3);
  hdpt->Draw("hist");
  c1->cd(4);
  hdphi->Draw("hist");
  c1->cd(5);
  hntrack->Draw("hist");
  c1->cd(6);
  hnextracaloe->Draw("hist");
  c1->cd(7);
  hm2->Draw("hist");
  }
}
