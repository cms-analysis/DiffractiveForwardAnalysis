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
#include "TF1.h"
#include <cmath>

void SelectorQED_MC_Data() //NO PU, no use of Calo Tower, cumulative TGraph
{
gROOT->Reset();
//definition des fichiers + Tree
  TFile *f0 = new TFile("../Root_files/ggmm_0PU_Ideal_COMPLETE.root"); // signal A156
  TTree *t0 = f0->Get("ntp1");
  TFile *f1 = new TFile("../Root_files/igmm_0PU_Ideal_COMPLETE.root"); // inelastic A156
  TTree *t1 = f1->Get("ntp1");
  TFile *f7 = new TFile("../Root_files/ExclusiveUpsilon_0PU_Ideal_COMPLETE.root"); // Exclusive Upsilon
  TTree *t7 = f7->Get("ntp1"); // << NOT USED ANYMORE
  TFile *f9 = new TFile("../Root_files/DiffDY_0PU_Ideal_COMPLETE.root"); // Diff Zmumu
  TTree *t9 = f9->Get("ntp1");
  TFile *f11 = new TFile("../Root_files/gampupsilon1smumu.pat.root"); // Exclusive Upsilon
  TTree *t11 = f11->Get("ntp1");
  TFile *f12 = new TFile("../Root_files/gampupsilon2smumu.pat.root"); // Exclusive Upsilon
  TTree *t12 = f12->Get("ntp1");
  TFile *f13 = new TFile("../Root_files/gampupsilon3smumu.pat.root"); // Exclusive Upsilon
  TTree *t13 = f13->Get("ntp1");



  
// definitions des Trees pour la sauvegarde
// MASS

// phi1-phi2 / pi

  TH1F *hmc_phi = new TH1F("hmc_phi","",20.,-0.1,0.1);
  TH1F *hinel_phi = new TH1F("hinel_phi","",20.,-0.1,0.1);
  TH1F *hdata_phi = new TH1F("hdata_phi","",20.,-0.1,0.1);

  TH1F *hmc_pt = new TH1F("hmc_pt","",20.,0.,2.);
  TH1F *hinel_pt = new TH1F("hinel_pt","",20.,0.,2.);
  TH1F *hdata_pt = new TH1F("hdata_pt","",20.,0.,2.);

  TH1F *hmc_zdc = new TH1F("hmc_zdc","",11.,-1.,10.);
  TH1F *hinel_zdc = new TH1F("hinel_zdc","",11.,-1.,10.);
  TH1F *hdata_zdc = new TH1F("hdata_zdc","",11.,-1.,10.);

  TH1F *hmc_castor = new TH1F("hmc_castor","",21.,-1.,20.);
  TH1F *hinel_castor = new TH1F("hinel_castor","",21.,-1.,20.);
  TH1F *hdata_castor = new TH1F("hdata_castor","",21.,-1.,20.);

  TH1F *hmc_mass = new TH1F("hmc_mass","",50.,0.,100.);
  TH1F *hinel_mass = new TH1F("hinel_mass","",50.,0.,100.);
  TH1F *hdata_mass = new TH1F("hdata_mass","",50.,0.,100.);

//Other histo (filled after excl. cut, before mass cut)

// definitions des # d'entrées
  const int NUM0 = t0->GetEntries();
  const int NUM1 = t1->GetEntries();
  const int NUM7 = t7->GetEntries();
  const int NUM9 = t9->GetEntries();
  const int NUM11 = t11->GetEntries();
  const int NUM12 = t12->GetEntries();
  const int NUM13 = t13->GetEntries();


  const double xsection = 64.1;
 
  const double Pi = 3.1415926536;
  const float binsig=0.040;

  const int N=0;
  const double lumi_int = 20; // en pb-1 

  double faclumi0 = 0.003205*lumi_int;
  double faclumi1 = 0.004040*lumi_int;
  double faclumi7 = 0.00203112*lumi_int;
  double faclumi9 = 0.02601479*lumi_int; // to be tuned
  double faclumi11 = 0.008233333*lumi_int;
  double faclumi12 = 0.002766666*lumi_int;
  double faclumi13 = 0.00247*lumi_int;

  const double mass_upsilon_max = 11.5; //GeV
  const double mass_upsilon_min = 8.5; //GeV
  const double dpt_max = 1.5;  //0.65 GeV
  const double dphi_min = 2.9;  //rad
  const double vtxT_max = 999.;//0.15; // cm
  const double closeTrk_min = 998.; //cm
  const int nextraCaloTower_max = 5; 
  const int nhitZDC_max = 999;

  const double pt_purity = 3.5;  //GeV
  const double eta_purity = 2.3;  //rad

  const bool purity = false;
  const bool draw_mc = true;

 cout << " ---- Exclusivity cuts ---- " << endl;
 cout << " d(pT)  < " << dpt_max << "GeV" << endl;
 cout << " d(phi) > " << dphi_min << "rad" << endl;
 cout << " vtxT   < " << vtxT_max << "cm" << endl;
 cout << " close track > " << closeTrk_min << "cm from µµ vertex" << endl;
 cout << " # Calo Tower < " << nextraCaloTower_max << endl;
 cout << " # ZDC Hits   < " << nhitZDC_max << endl;
 cout << " excluded mass = [" << mass_upsilon_min << " ; " << mass_upsilon_max << "] GeV"<< endl;
 cout << " " << endl;
 cout << " ----    Purity cuts   ---- " << endl;
 cout << " pT(µ)  < " << pt_purity << "GeV" << endl;
 cout << " |eta|  > " << eta_purity << endl;
 if(purity==true){cout << " Applied = true"<< endl;}
 else{cout << " Applied = false"<< endl;}
 cout << " ---- ---------------- ---- " << endl;


//definition des variables
//HLT Double Mu3 only
    Int_t hlt_d0[1], hlt_d1[1], hlt_d7[1], hlt_d9[1],  hlt_d11[1], hlt_d12[1], hlt_d13[1];

//Mass
  Double_t var_a0[1], var_a1[1], var_a7[1], var_a9[1], var_a11[1], var_a12[1],var_a13[1],;

// dpT
  Double_t var_b0[1], var_b1[1], var_b7[1], var_b9[1], var_b11[1], var_b12[1], var_b13[1];

// dphi
  Double_t var_c0[1], var_c1[1], var_c7[1], var_c9[1], var_c11[1], var_c12[1], var_c13[1];

// vtxT
  Double_t var_d0[1], var_d1[1], var_d7[1], var_d9[1], var_d11[1], var_d12[1], var_d13[1];

// MuonID
     Int_t var_e0[2], var_e1[2], var_e7[2], var_e9[2], var_e11[2], var_e12[2], var_e13[2];

// Closest extra track
  Double_t var_f0[1], var_f1[1], var_f7[1], var_f9[1], var_f11[1], var_f12[1], var_f13[1];
 
// NextraCaloE5
     Int_t var_g0[1], var_g1[1], var_g7[1], var_g9[1], var_g11[1], var_g12[1], var_g13[1];

// Hit ZDC
     Int_t var_k0[1], var_k1[1], var_k7[1], var_k9[1], var_k11[1], var_k12[1], var_k13[1];
// Hit Castor
     Int_t var_kk0[1], var_kk1[1], var_kk7[1], var_kk9[1], var_kk11[1], var_kk12[1], var_kk13[1];

// Muon Phi
  Double_t var_h0[2], var_h1[2], var_h7[2], var_h9[2], var_h11[2], var_h12[2], var_h13[2];
// Muon pT
  Double_t var_i0[2], var_i1[2], var_i7[2], var_i9[2], var_i11[2], var_i12[2], var_i13[2];
// Muon eta
  Double_t var_j0[2], var_j1[2], var_j7[2], var_j9[2], var_j11[2], var_j12[2], var_j13[2];
 
  t0->SetBranchAddress("HLT2MuonNonIso",hlt_d0);
  t1->SetBranchAddress("HLT2MuonNonIso",hlt_d1);
  t7->SetBranchAddress("HLT2MuonNonIso",hlt_d7);
  t9->SetBranchAddress("HLT2MuonNonIso",hlt_d9);
  t11->SetBranchAddress("HLT_DoubleMu3",hlt_d11);
  t12->SetBranchAddress("HLT_DoubleMu3",hlt_d12);
  t13->SetBranchAddress("HLT_DoubleMu3",hlt_d13);

  t0->SetBranchAddress("MuMu_mass",var_a0);
  t1->SetBranchAddress("MuMu_mass",var_a1);
  t7->SetBranchAddress("MuMu_mass",var_a7);
  t9->SetBranchAddress("MuMu_mass",var_a9);
  t11->SetBranchAddress("MuMu_mass",var_a11);
  t12->SetBranchAddress("MuMu_mass",var_a12);
  t13->SetBranchAddress("MuMu_mass",var_a13);

  t0->SetBranchAddress("MuMu_dpt",var_b0);
  t1->SetBranchAddress("MuMu_dpt",var_b1);
  t7->SetBranchAddress("MuMu_dpt",var_b7);
  t9->SetBranchAddress("MuMu_dpt",var_b9);
  t11->SetBranchAddress("MuMu_dpt",var_b11);
  t12->SetBranchAddress("MuMu_dpt",var_b12);
  t13->SetBranchAddress("MuMu_dpt",var_b13);

  t0->SetBranchAddress("MuMu_dphi",var_c0);
  t1->SetBranchAddress("MuMu_dphi",var_c1);
  t7->SetBranchAddress("MuMu_dphi",var_c7);
  t9->SetBranchAddress("MuMu_dphi",var_c9);
  t11->SetBranchAddress("MuMu_dphi",var_c11);
  t12->SetBranchAddress("MuMu_dphi",var_c12);
  t13->SetBranchAddress("MuMu_dphi",var_c13);

  t0->SetBranchAddress("MuMu_vtxT",var_d0);
  t1->SetBranchAddress("MuMu_vtxT",var_d1);
  t7->SetBranchAddress("MuMu_vtxT",var_d7);
  t9->SetBranchAddress("MuMu_vtxT",var_d9);
  t11->SetBranchAddress("MuMu_vtxT",var_d11);
  t12->SetBranchAddress("MuMu_vtxT",var_d12);
  t13->SetBranchAddress("MuMu_vtxT",var_d13);

  t0->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e0);
  t1->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e1);
  t7->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e7);
  t9->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e9);
  t11->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e11);
  t12->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e12);
  t13->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_e13);

  t0->SetBranchAddress("MuonCand_phi",var_h0);
  t1->SetBranchAddress("MuonCand_phi",var_h1);
  t7->SetBranchAddress("MuonCand_phi",var_h7);
  t9->SetBranchAddress("MuonCand_phi",var_h9);
  t11->SetBranchAddress("MuonCand_phi",var_h11);
  t12->SetBranchAddress("MuonCand_phi",var_h12);
  t13->SetBranchAddress("MuonCand_phi",var_h13);

  t0->SetBranchAddress("MuonCand_pt",var_i0);
  t1->SetBranchAddress("MuonCand_pt",var_i1);
  t7->SetBranchAddress("MuonCand_pt",var_i7);
  t9->SetBranchAddress("MuonCand_pt",var_i9);
  t11->SetBranchAddress("MuonCand_pt",var_i11);
  t12->SetBranchAddress("MuonCand_pt",var_i12);
  t13->SetBranchAddress("MuonCand_pt",var_i13);

  t0->SetBranchAddress("MuonCand_eta",var_j0);
  t1->SetBranchAddress("MuonCand_eta",var_j1);
  t7->SetBranchAddress("MuonCand_eta",var_j7);
  t9->SetBranchAddress("MuonCand_eta",var_j9);
  t11->SetBranchAddress("MuonCand_eta",var_j11);
  t12->SetBranchAddress("MuonCand_eta",var_j12);
  t13->SetBranchAddress("MuonCand_eta",var_j13);


  t0->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f0);
  t1->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f1);
  t7->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f7);
  t9->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f9);
  t11->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f11);
  t12->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f12);
  t13->SetBranchAddress("ClosestExtraTrack_vtxdxyz",var_f13);


  t0->SetBranchAddress("nExtraCaloTowersE5",var_g0);
  t1->SetBranchAddress("nExtraCaloTowersE5",var_g1);
  t7->SetBranchAddress("nExtraCaloTowersE5",var_g7);
  t9->SetBranchAddress("nExtraCaloTowersE5",var_g9);
  t11->SetBranchAddress("nExtraCaloTowersE5",var_g11);
  t12->SetBranchAddress("nExtraCaloTowersE5",var_g12);
  t13->SetBranchAddress("nExtraCaloTowersE5",var_g13);


  t0->SetBranchAddress("HitInZDC",var_k0);
  t1->SetBranchAddress("HitInZDC",var_k1);
  t7->SetBranchAddress("HitInZDC",var_k7);
  t9->SetBranchAddress("HitInZDC",var_k9);
  t11->SetBranchAddress("HitInZDC",var_k11);
  t12->SetBranchAddress("HitInZDC",var_k12);
  t13->SetBranchAddress("HitInZDC",var_k13);

  t0->SetBranchAddress("HitInCastor",var_kk0);
  t1->SetBranchAddress("HitInCastor",var_kk1);
  t7->SetBranchAddress("HitInCastor",var_kk7);
  t9->SetBranchAddress("HitInCastor",var_kk9);
  t11->SetBranchAddress("HitInCastor",var_kk11);
  t12->SetBranchAddress("HitInCastor",var_kk12);
  t13->SetBranchAddress("HitInCastor",var_kk13);




// Boucle #0
  cout << " " << endl;
  cout << " - selected el-el" << endl;
  Int_t hlt_0(0);
  Int_t pass_0(0);
  Int_t excl_0(0);
  Int_t mass_0(0);
  Int_t pure_0(0);
  Int_t i_true(0);
  Int_t i_exp(0);
  for(Int_t i = 0;i < NUM0;i++){
        t0->GetEntry(i);
        double muID1 = var_e0[0];
        double muID2 = var_e0[1];
	double mass  = var_a0[0];
        int hlt_pass = hlt_d0[0];
        if(hlt_pass==1){
         hlt_0++;  
	 if(muID1==1 && muID2==1){
	   pass_0++;

	   if(var_b0[0]<dpt_max&&var_c0[0]>dphi_min&&var_d0[0]<vtxT_max&&var_f0[0]>closeTrk_min&&var_g0[0]<nextraCaloTower_max&&var_k0[0]<nhitZDC_max){
	      excl_0++;
	      if(!(mass<=mass_upsilon_max&&mass>=mass_upsilon_min)){
		mass_0++;
		if(purity==false){
			double phi0 = var_h0[0] > 0 ? var_h0[0] : -var_h0[0];
        	        double phi1 = (var_h0[1]+Pi) > 0 ? (var_h0[1]+Pi) : -(var_h0[1]+Pi);
			phi1 = phi1 < Pi ? phi1 : fabs(phi1-2*Pi);
                        double dpt0 = var_b0[0];

	                hmc_phi->Fill((phi0 - phi1)/Pi,faclumi0);
                        hmc_pt->Fill(dpt0,faclumi0);
                        hmc_zdc->Fill(var_k0[0],faclumi0);
                        hmc_castor->Fill(var_kk0[0],faclumi0);
                        hmc_mass->Fill(mass,faclumi0);

			if(4.926*lumi_int <= 1537){
	         		if(i_true < 4.926*lumi_int) {hdata_phi->Fill((phi0 - phi1)/Pi,1.); i_true++;
                                                            hdata_pt->Fill(dpt0,1.);
                                                            hdata_zdc->Fill(var_k0[0],1.);
                                                            hdata_castor->Fill(var_kk0[0],1.);
                                                            hdata_mass->Fill(mass,1.);}
			}
			else cout <<"!!!!!!!!"<<endl ;
		}
		else{ //if purity = true
                 cout << "WARNING: PURITY ON!!!" << endl;	
		}
             }
	   }
         }
	}
  }
  cout << "" << i_true << endl;
  double eff0_0=(hlt_0*100.)/(NUM0);
  double eff1_0=(pass_0*100.)/(hlt_0);
  double eff2_0=(excl_0*100.)/(pass_0);
  double eff3_0=(mass_0*100.)/(excl_0);
  double eff4_0=(pure_0*100.)/(mass_0);
  cout << "hlt  eff = "<<hlt_0 <<"/"<<NUM0<<" = \t"<< eff0_0 <<"%"<< endl;
  cout << "muID eff = "<<pass_0<<"/"<<hlt_0 <<" = \t"<< eff1_0 <<"%"<< endl;
  cout << "excl eff = "<<excl_0<<"/"<<pass_0<<" = \t"<< eff2_0 <<"%"<< endl;
  cout << "mass eff = "<<mass_0<<"/"<<excl_0<<" = \t"<< eff3_0 <<"%"<< endl;
  if(purity==true){ cout << "pure eff = "<<pure_0<<"/"<<mass_0<<" = \t"<< eff4_0 <<"%"<< endl;
 		    double effselsignal = (pure_0*1.)/(hlt_0);}
  else{double effselsignal = (mass_0*1.)/(hlt_0);}
  double efftrigger = (hlt_0*1.)/20000;
  


// Boucle #1
  cout << " " << endl;
  cout << " - selected inel-el" << endl;
  Int_t hlt_1(0);
  Int_t pass_1(0);
  Int_t excl_1(0);
  Int_t mass_1(0);
  Int_t pure_1(0);
  Int_t i_true(0);
  Int_t i_exp(0);
  for(Int_t i = 0;i < NUM1;i++){
      	t1->GetEntry(i);
        double muID1 = var_e1[0];
        double muID2 = var_e1[1];
        double mass  = var_a1[0];
        int hlt_pass = hlt_d1[0];
        if(hlt_pass==1){
          hlt_1++;
          if(muID1==1 && muID2==1) {
		pass_1++;
           if(var_b1[0]<dpt_max&&var_c1[0]>dphi_min&&var_d1[0]<vtxT_max&&var_f1[0]>closeTrk_min&&var_g1[0]<nextraCaloTower_max&&var_k1[0]<nhitZDC_max){
	     excl_1++;
	     if(!(mass<=mass_upsilon_max&&mass>=mass_upsilon_min)){
                mass_1++;
		if(purity == false){
	                double phi0 = var_h1[0] > 0 ? var_h1[0] : -var_h1[0];
	                double phi1 = (var_h1[1]+Pi) > 0 ? (var_h1[1]+Pi) : -(var_h1[1]+Pi);
	                phi1 = phi1 < Pi ? phi1 : fabs(phi1-2*Pi);
                        double dpt1 = var_b1[0];

	                hmc_phi->Fill((phi0 - phi1)/Pi,faclumi1);
                        hmc_pt->Fill(dpt1,faclumi1);
                        hmc_zdc->Fill(var_k1[0],faclumi1);
                        hmc_castor->Fill(var_kk1[0],faclumi1);
                        hmc_mass->Fill(mass,faclumi1);

                        hinel_phi->Fill((phi0 - phi1)/Pi,faclumi1);
                        hinel_pt->Fill(dpt1,faclumi1);
                        hinel_zdc->Fill(var_k1[0],faclumi1);
                        hinel_castor->Fill(var_kk1[0],faclumi1);
                        hinel_mass->Fill(mass,faclumi1);

	                if(4.29*lumi_int <= 1062){
				if(i_true < 4.29*lumi_int) {hdata_phi->Fill((phi0 - phi1)/Pi,1.); i_true++;
							    hdata_pt->Fill(dpt1,1.);
				                            hdata_zdc->Fill(var_k1[0],1.);
				                            hdata_castor->Fill(var_kk1[0],1.);
                                                            hdata_mass->Fill(mass,1.);}
                              i_exp++;
			}
			else {cout << "!!!!!!!!!!!" << endl;}
		}
		else{ // if purity is true
                  cout << "WARNING: PURITY IS ON !!!" << endl;
		}
	     }	
	   }
	}
     }
  }
  double eff0_1=(hlt_1*100.)/(NUM1);
  double eff1_1=(pass_1*100.)/(hlt_1);
  double eff2_1=(excl_1*100.)/(pass_1);
  double eff3_1=(mass_1*100.)/(excl_1);
  double eff4_1=(pure_1*100.)/(mass_1);
  cout << "hlt eff  = "<<hlt_1 <<"/"<<NUM1<<" = \t"<< eff0_1 <<"%"<< endl;
  cout << "muID eff = "<<pass_1<<"/"<<hlt_1 <<" = \t"<< eff1_1 <<"%"<< endl;
  cout << "excl eff = "<<excl_1<<"/"<<pass_1<<" = \t"<< eff2_1 <<"%"<< endl;
  cout << "mass eff = "<<mass_1<<"/"<<excl_1<<" = \t"<< eff3_1 <<"%"<< endl;
  if(purity==true) cout << "pure eff = "<<pure_1<<"/"<<mass_1<<" = \t"<< eff4_1 <<"%"<< endl;

// Boucle #11
  cout << " " << endl;
  cout << " - Exclusive Upsilon 1S" << endl; 
  Int_t hlt_11(0);
  Int_t pass_11(0);
  Int_t excl_11(0);
  Int_t mass_11(0);
  Int_t pure_11(0);
  for(Int_t i = 0;i < NUM11;i++){
        t11->GetEntry(i);
        double muID1 = var_e11[0];
        double muID2 = var_e11[1];
        double mass  = var_a11[0];
        int hlt_pass = hlt_d11[0];
        if(hlt_pass==1){
          hlt_11++;
          if(muID1==1 && muID2==1) {
                pass_11++;
           if(var_b11[0]<dpt_max&&var_c11[0]>dphi_min&&var_d11[0]<vtxT_max&&var_f11[0]>closeTrk_min&&var_g11[0]<nextraCaloTower_max&&var_k11[0]<nhitZDC_max){
	     excl_11++;
             if(!(mass<=mass_upsilon_max&&mass>=mass_upsilon_min)){
                mass_11++;
		if(purity==true){
                     if(var_i11[0]>pt_purity && var_i11[1]>pt_purity &&
                     fabs(var_j11[0])<eta_purity && fabs(var_j11[1])<eta_purity)
                        pure_11++;
		}
	     }
	   }
        }
     }
  }
  double eff0_11=(hlt_11*100)/(NUM11);
  double eff1_11=(pass_11*100.)/(hlt_11);
  double eff2_11=(excl_11*100.)/(pass_11);
  double eff3_11=(mass_11*100.)/(excl_11);
//  double eff4_11=(pure_11*100.)/(mass_11);
  cout << "hlt  eff = "<<hlt_11 <<"/"<<NUM11<<" = \t"<< eff0_11 <<"%"<< endl;
  cout << "muID eff = "<<pass_11<<"/"<<hlt_11 <<" = \t"<< eff1_11 <<"%"<< endl;
  cout << "excl eff = "<<excl_11<<"/"<<pass_11<<" = \t"<< eff2_11 <<"%"<< endl;
  cout << "mass eff = "<<mass_11<<"/"<<excl_11<<" = \t"<< eff3_11 <<"%"<< endl;
//  if(purity==true) cout << "pure eff = "<<pure_11<<"/"<<mass_11<<" = \t"<< eff4_11 <<"%"<< endl;

// Boucle #12
  cout << " " << endl;
  cout << " - Exclusive Upsilon 2S" << endl; 
  Int_t hlt_12(0);
  Int_t pass_12(0);
  Int_t excl_12(0);
  Int_t mass_12(0);
  Int_t pure_12(0);
  for(Int_t i = 0;i < NUM12;i++){
        t12->GetEntry(i);
        double muID1 = var_e12[0];
        double muID2 = var_e12[1];
        double mass  = var_a12[0];
        int hlt_pass = hlt_d12[0];
        if(hlt_pass==1){
          hlt_12++;
          if(muID1==1 && muID2==1) {
                pass_12++;
           if(var_b12[0]<dpt_max&&var_c12[0]>dphi_min&&var_d12[0]<vtxT_max&&var_f12[0]>closeTrk_min&&var_g12[0]<nextraCaloTower_max&&var_k12[0]<nhitZDC_max){
	     excl_12++;
             if(!(mass<=mass_upsilon_max&&mass>=mass_upsilon_min)){
                mass_12++;
		if(purity==true){
                     if(var_i12[0]>pt_purity && var_i12[1]>pt_purity &&
                     fabs(var_j12[0])<eta_purity && fabs(var_j12[1])<eta_purity)
                        pure_12++;
		}
	     }
	   }
        }
     }
  }
  double eff0_12=(hlt_12*100)/(NUM12);
  double eff1_12=(pass_12*100.)/(hlt_12);
  double eff2_12=(excl_12*100.)/(pass_12);
  double eff3_12=(mass_12*100.)/(excl_12);
//  double eff4_12=(pure_12*100.)/(mass_12);
  cout << "hlt  eff = "<<hlt_12 <<"/"<<NUM12<<" = \t"<< eff0_12 <<"%"<< endl;
  cout << "muID eff = "<<pass_12<<"/"<<hlt_12 <<" = \t"<< eff1_12 <<"%"<< endl;
  cout << "excl eff = "<<excl_12<<"/"<<pass_12<<" = \t"<< eff2_12 <<"%"<< endl;
  cout << "mass eff = "<<mass_12<<"/"<<excl_12<<" = \t"<< eff3_12 <<"%"<< endl;
//  if(purity==true) cout << "pure eff = "<<pure_12<<"/"<<mass_12<<" = \t"<< eff4_12 <<"%"<< endl;

// Boucle #13
  cout << " " << endl;
  cout << " - Exclusive Upsilon 3S" << endl; 
  Int_t hlt_13(0);
  Int_t pass_13(0);
  Int_t excl_13(0);
  Int_t mass_13(0);
  Int_t pure_13(0);
  for(Int_t i = 0;i < NUM13;i++){
        t13->GetEntry(i);
        double muID1 = var_e13[0];
        double muID2 = var_e13[1];
        double mass  = var_a13[0];
        int hlt_pass = hlt_d13[0];
        if(hlt_pass==1){
          hlt_13++;
          if(muID1==1 && muID2==1) {
                pass_13++;
           if(var_b13[0]<dpt_max&&var_c13[0]>dphi_min&&var_d13[0]<vtxT_max&&var_f13[0]>closeTrk_min&&var_g13[0]<nextraCaloTower_max&&var_k13[0]<nhitZDC_max){
	     excl_13++;
             if(!(mass<=mass_upsilon_max&&mass>=mass_upsilon_min)){
                mass_13++;
		if(purity==true){
                     if(var_i13[0]>pt_purity && var_i13[1]>pt_purity &&
                     fabs(var_j13[0])<eta_purity && fabs(var_j13[1])<eta_purity)
                        pure_13++;
		}
	     }
	   }
        }
     }
  }
  double eff0_13=(hlt_13*100)/(NUM13);
  double eff1_13=(pass_13*100.)/(hlt_13);
  double eff2_13=(excl_13*100.)/(pass_13);
  double eff3_13=(mass_13*100.)/(excl_13);
//  double eff4_13=(pure_13*100.)/(mass_13);
  cout << "hlt  eff = "<<hlt_13 <<"/"<<NUM13<<" = \t"<< eff0_13 <<"%"<< endl;
  cout << "muID eff = "<<pass_13<<"/"<<hlt_13 <<" = \t"<< eff1_13 <<"%"<< endl;
  cout << "excl eff = "<<excl_13<<"/"<<pass_13<<" = \t"<< eff2_13 <<"%"<< endl;
  cout << "mass eff = "<<mass_13<<"/"<<excl_13<<" = \t"<< eff3_13 <<"%"<< endl;
//  if(purity==true) cout << "pure eff = "<<pure_13<<"/"<<mass_13<<" = \t"<< eff4_13 <<"%"<< endl;


/*
// Boucle #9
  cout << " " << endl;
  cout << " - Diffractive Z/gamma*" << endl;
  Int_t pass_9(0);
  Int_t excl_9(0);
  Int_t mass_9(0);
  Int_t pure_9(0);
  Int_t i_true(0);
  for(Int_t i = 0;i < NUM9;i++){
        t9->GetEntry(i);
        double muID1 = var_e9[0];
        double muID2 = var_e9[1];
        double mass  = var_a9[0];
        if(muID1==1 && muID2==1) {
                pass_9++;
           if(var_b9[0]<dpt_max&&var_c9[0]>dphi_min&&var_d9[0]<vtxT_max&&var_f9[0]>closeTrk_min&&var_g9[0]<nextraCaloTower_max&&var_k9[0]<nhitZDC_max){
	     excl_9++;
             if(!(mass<=mass_upsilon_max&&mass>=mass_upsilon_min)){
		mass_9++;
		if(purity==false){
	                double phi0 = var_h9[0] > 0 ? var_h9[0] : -var_h9[0];
        	        double phi1 = (var_h9[1]+Pi) > 0 ? (var_h9[1]+Pi) : -(var_h9[1]+Pi);
	                phi1 = phi1 < Pi ? phi1 : fabs(phi1-2*Pi);

	                hbkgmc->Fill((phi0 - phi1)/Pi,faclumi9);

        	        if(0.12*lumi_int <= 9){
                	        if(i_true < 0.12*lumi_int) {hbkg->Fill((phi0 - phi1)/Pi,1.); i_true++;
							    hdy->Fill((phi0 - phi1)/Pi,1.);}
	                }
        	        else {hbkg->Fill((phi0 - phi1)/Pi,faclumi9);
			      hdy->Fill((phi0 - phi1)/Pi,faclumi9);}
		}
		else{ // if purity = true
		   if(var_i9[0]>pt_purity && var_i9[1]>pt_purity &&
                     fabs(var_j9[0])<eta_purity && fabs(var_j9[1])<eta_purity){
			pure_9++;
                        double phi0 = var_h9[0] > 0 ? var_h9[0] : -var_h9[0];
                        double phi1 = (var_h9[1]+Pi) > 0 ? (var_h9[1]+Pi) : -(var_h9[1]+Pi);
                        phi1 = phi1 < Pi ? phi1 : fabs(phi1-2*Pi);

                        hbkgmc->Fill((phi0 - phi1)/Pi,faclumi9);

                        if(0.05*lumi_int <= 4){
                                if(i_true < 0.05*lumi_int) {hbkg->Fill((phi0 - phi1)/Pi,1.); i_true++;
							    hdy->Fill((phi0 - phi1)/Pi,1.);}
                        }
                        else {hbkg->Fill((phi0 - phi1)/Pi,faclumi9);
			      hdy->Fill((phi0 - phi1)/Pi,faclumi9);}
		   }
		}
	    }
          }
        }
  }
  double eff1_9=(pass_9*100.)/(NUM9);
  double eff2_9=(excl_9*100.)/(pass_9);
  double eff3_9=(mass_9*100.)/(excl_9);
  double eff4_9=(pure_9*100.)/(mass_9);
  cout << "muID eff = "<<pass_9<<"/"<<NUM9<<" = \t"<< eff1_9 <<"%"<< endl;
  cout << "excl eff = "<<excl_9<<"/"<<pass_9<<" = \t"<< eff2_9 <<"%"<< endl;
  cout << "mass eff = "<<mass_9<<"/"<<excl_9<<" = \t"<< eff3_9 <<"%"<< endl;
  if(purity==true) cout << "pure eff = "<<pure_9<<"/"<<mass_9<<" = \t"<< eff4_9 <<"%"<< endl;
  cout << "\n \n \n" << endl;
*/

// Save into histo

//Draw
TFile *file = new TFile("Output_SelectorCumul.root","recreate");
if(draw_mc==true){
//-------dPt
TCanvas *c_dpt = new TCanvas("dPt","delta(pT)");
   c_dpt->SetFillColor(0);
   c_dpt->SetBorderMode(0);
   c_dpt->SetBorderSize(2);
   c_dpt->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

   TPaveText *qed = new TPaveText(1.010887,56.2797,1.933908,62.03707,"br");
   qed->SetBorderSize(0);
   qed->SetFillColor(0);
   char text[500]; sprintf(text,"QED signal: m_{#mu#mu} < %g GeV and m_{#mu#mu} > %g GeV",mass_upsilon_min,mass_upsilon_max);
   qed->AddText(text);
 
   c_dpt->cd();
   hdata_pt->SetLineWidth(3);
   hdata_pt->SetMarkerStyle(20);
   hdata_pt->SetMarkerSize(1.4);
   hdata_pt->SetXTitle("|#Delta(p_{T})| [GeV]");
   hdata_pt->SetYTitle("Events/0.1 GeV");
   hdata_pt->GetXaxis()->SetNdivisions(509);
   hdata_pt->Sumw2();
   hdata_pt->Draw();

   hmc_pt->SetFillColor(1);
   hmc_pt->SetFillStyle(3007);
   hmc_pt->SetLineWidth(2);
   hmc_pt->Draw("same");

   hinel_pt->SetFillColor(1);
   hinel_pt->SetFillStyle(3005);
   hinel_pt->SetLineWidth(2);
   hinel_pt->Draw("same");

   qed->Draw();

   TLegend * legdpt = new TLegend(0.56,0.66,0.88,0.94,"");
    legdpt->SetBorderSize(0);
    legdpt->SetFillStyle(0);
    legdpt->SetTextFont(42);
    legdpt->SetTextSizePixels(23);
//    legdpt->AddEntry(hdpt_signal,"el-el #gamma#gamma#rightarrow#mu#mu","f");
//    legdpt->AddEntry(hdpt_inel,"inel-el #gamma#gamma#rightarrow#mu#mu","l");
//    legdpt->AddEntry(hdpt_ups,"exclusive #Upsilon#rightarrow#mu#mu","l");
/*    legdpt->AddEntry(hdpt_diffdy,"diff. Z/#gamma*#rightarrow#mu#mu","l");*/
    legdpt->Draw();

//-------dPhi
TCanvas *c_dphi = new TCanvas("dphi","delta(phi)");
   c_dphi->SetFillColor(0);
   c_dphi->SetBorderMode(0);
   c_dphi->SetBorderSize(2);
   c_dphi->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

   c_dphi->cd();
   hdata_phi->SetLineWidth(3);
   hdata_phi->SetMarkerStyle(20);
   hdata_phi->SetMarkerSize(1.4);
   hdata_phi->SetXTitle("|#Delta(#phi)|/#pi [GeV]");
   hdata_phi->SetYTitle("Events/0.002");
   hdata_phi->GetXaxis()->SetNdivisions(509);
   hdata_phi->Sumw2();
   hdata_phi->Draw();

   hmc_phi->SetFillColor(1);
   hmc_phi->SetFillStyle(3007);
   hmc_phi->SetLineWidth(2);
   hmc_phi->Draw("same");

   hinel_phi->SetFillColor(1);
   hinel_phi->SetFillStyle(3005);
   hinel_phi->SetLineWidth(2);
   hinel_phi->Draw("same");
   qed->Draw();

//-------ZDC
TCanvas *c_zdc = new TCanvas("ZDC","ZDC hits");
   c_zdc->SetFillColor(0);
   c_zdc->SetBorderMode(0);
   c_zdc->SetBorderSize(2);
   c_zdc->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

   c_zdc->cd();
   hdata_zdc->SetLineWidth(3);
   hdata_zdc->SetMarkerStyle(20);
   hdata_zdc->SetMarkerSize(1.4);
   hdata_zdc->SetXTitle("n ZDC hits");
   hdata_zdc->SetYTitle("Events/hit");
   hdata_zdc->GetXaxis()->SetNdivisions(509);
   hdata_zdc->Sumw2();
   hdata_zdc->Draw();

   hmc_zdc->SetFillColor(1);
   hmc_zdc->SetFillStyle(3007);
   hmc_zdc->SetLineWidth(2);
   hmc_zdc->Draw("same");

   hinel_zdc->SetFillColor(1);
   hinel_zdc->SetFillStyle(3005);
   hinel_zdc->SetLineWidth(2);
   hinel_zdc->Draw("same");
   qed->Draw();

//-------castor
TCanvas *c_castor = new TCanvas("Castor","Castor h");
   c_castor->SetFillColor(0);
   c_castor->SetBorderMode(0);
   c_castor->SetBorderSize(2);
   c_castor->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

   c_castor->cd();
   hdata_castor->SetLineWidth(3);
   hdata_castor->SetMarkerStyle(20);
   hdata_castor->SetMarkerSize(1.4);
   hdata_castor->SetXTitle("n ZDC hits");
   hdata_castor->SetYTitle("Events/hit");
   hdata_castor->GetXaxis()->SetNdivisions(509);

   hdata_castor->Sumw2();
   hdata_castor->Draw();

   hmc_castor->SetFillColor(1);
   hmc_castor->SetFillStyle(3007);
   hmc_castor->SetLineWidth(2);
   hmc_castor->Draw("same");

   hinel_castor->SetFillColor(1);
   hinel_castor->SetFillStyle(3005);
   hinel_castor->SetLineWidth(2);
   hinel_castor->Draw("same");
   qed->Draw();

//----Mass
/*
TCanvas *c_mass = new TCanvas("mass","invariant mu mu mass");
   c_mass->SetFillColor(0);
   c_mass->SetBorderMode(0);
   c_mass->SetBorderSize(2);
   c_mass->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

   c_mass->cd();
   hdata_mass->SetLineWidth(3);
   hdata_mass->SetMarkerStyle(20);
   hdata_mass->SetMarkerSize(1.4);
   hdata_mass->SetXTitle("m_{#mu#mu} [GeV]");
   hdata_mass->SetYTitle("Events/2 GeV");
   hdata_mass->GetXaxis()->SetNdivisions(509);
   hdata_mass->Sumw2();
   hdata_mass->Draw();

   hmc_mass->SetFillColor(1);
   hmc_mass->SetFillStyle(3007);
   hmc_mass->SetLineWidth(2);
   hmc_mass->Draw("same");

   hinel_mass->SetFillColor(1);
   hinel_mass->SetFillStyle(3005);
   hinel_mass->SetLineWidth(2);
   hinel_mass->Draw("same");
   qed->Draw();
*/
}

//-----FITTING

cout << "END" << endl;   
}

