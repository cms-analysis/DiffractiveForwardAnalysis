#define Dimuons2017Macro_cxx
#include "Dimuons2017Macro.h"
#include <TGraph.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TH1.h>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"


#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>

#include "RoccoR.cc"


#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "./align/alignment_classes.h"
#include "./align/fill_info.h"
#include "./align/fill_info.cc"

using namespace std;
int getXangle(int run,int lumi, const char* filename);


void Dimuons2017Macro::Loop()
{
int preTS2=1;


   if (fChain == 0) return;

   TH1F *hacop = new TH1F("hacop","hacop",1000,0,0.01);
   TH2F *hcorr45 = new TH2F("hcorr45","hcorr45",500,0,0.25,500,0,0.25);
   TH2F *hcorr56 = new TH2F("hcorr56","hcorr56",500,0,0.25,500,0,0.25);

   TH1F *hres45 = new TH1F("hres45","hres45",300,-5,1);
   TH1F *hres56 = new TH1F("hres56","hres56",300,-5,1);
   TH1F *hressum = new TH1F("hressum","hressum",300,-5,1);

   TH1F *hbin145 = new TH1F("hbin145","hbin145",300,-5,1);
   TH1F *hbin156 = new TH1F("hbin156","hbin156",300,-5,1);
   TH1F *hbin1sum = new TH1F("hbin1sum","hbin1sum",300,-5,1);
   TH1F *hbin245 = new TH1F("hbin245","hbin245",300,-5,1);
   TH1F *hbin256 = new TH1F("hbin256","hbin256",300,-5,1);
   TH1F *hbin2sum = new TH1F("hbin2sum","hbin2sum",300,-5,1);
   TH1F *hbin345 = new TH1F("hbin345","hbin345",300,-5,1);
   TH1F *hbin356 = new TH1F("hbin356","hbin356",300,-5,1);
   TH1F *hbin3sum = new TH1F("hbin3sum","hbin3sum",300,-5,1);

   TH1F *hres45mult = new TH1F("hres45mult","hres45mult",300,-5,1);
   TH1F *hres56mult = new TH1F("hres56mult","hres56mult",300,-5,1);
   TH1F *hressummult = new TH1F("hressummult","hressummult",300,-5,1);

RoccoR  rc("RoccoR2017.txt");


/////////////////////////////////////////

// b2
 TFile *f_in_B2_120 = new TFile("./inp/updatedshiftoptics/xi_as_a_function_of_x_graph_b2_120_murad.root");
 TFile *f_in_B2_130 = new TFile("./inp/updatedshiftoptics/xi_as_a_function_of_x_graph_b2_130_murad.root");
 TFile *f_in_B2_140 = new TFile("./inp/updatedshiftoptics/xi_as_a_function_of_x_graph_b2_140_murad.root");

 TGraph *m_s_x_to_xi_120_l   = (TGraph *) f_in_B2_120->Get("XRPH_B6L5_B2");
 TGraph *m_s_x_to_xi_130_l   = (TGraph *) f_in_B2_130->Get("XRPH_B6L5_B2");
 TGraph *m_s_x_to_xi_140_l   = (TGraph *) f_in_B2_140->Get("XRPH_B6L5_B2");


// b1
 TFile *f_in_B1_120 = new TFile("./inp/updatedshiftoptics/xi_as_a_function_of_x_graph_b1_120_murad.root");
 TFile *f_in_B1_130 = new TFile("./inp/updatedshiftoptics/xi_as_a_function_of_x_graph_b1_130_murad.root");
 TFile *f_in_B1_140 = new TFile("./inp/updatedshiftoptics/xi_as_a_function_of_x_graph_b1_140_murad.root");

 TGraph *m_s_x_to_xi_120_r   = (TGraph *) f_in_B1_120->Get("XRPH_B6R5_B1");
 TGraph *m_s_x_to_xi_130_r   = (TGraph *) f_in_B1_130->Get("XRPH_B6R5_B1"); 
 TGraph *m_s_x_to_xi_140_r   = (TGraph *) f_in_B1_140->Get("XRPH_B6R5_B1");

/////////////////////////////////////////

double angle=0;


TGraph *graph[2];
graph[0] = new TGraph();
graph[1] = new TGraph();


 AlignmentResults *alignments = NULL; 
AlignmentResultsCollection alignmentCollection;
if(preTS2==1) alignmentCollection.Load("./align/collect_alignments_2018_10_26.4.out");
if(preTS2==0) alignmentCollection.Load("./align/collect_alignments_2018_10_25.5.out"); // post-TS2
InitFillInfoCollection();


   // 0.02-0.05, 0.05-0.075, and 0.075-0.25
   
   Int_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry % 10000 == 0)
	std::cout << "Entry " << jentry << "/" << nentries << std::endl;


//////////////////////////////////////
double dataSFa=0;
double dataSFb=0;
double corrMass=0;
dataSFa = rc.kScaleDT(MuonCand_charge[0], MuonCand_pt[0],MuonCand_eta[0],MuonCand_phi[0], 0, 0);
dataSFb = rc.kScaleDT(MuonCand_charge[1], MuonCand_pt[1],MuonCand_eta[1],MuonCand_phi[1], 0, 0);
double pt0cor=MuonCand_pt[0]*dataSFa;
double pt1cor=MuonCand_pt[1]*dataSFb;
TLorentzVector l1;
TLorentzVector l2;

TLorentzVector mutmp1;
TLorentzVector mutmp2;
mutmp1.SetPtEtaPhiE(MuonCand_pt[0], MuonCand_eta[0], MuonCand_phi[0], MuonCand_e[0]);
mutmp2.SetPtEtaPhiE(MuonCand_pt[1], MuonCand_eta[1], MuonCand_phi[1], MuonCand_e[1]);

TLorentzVector pair;
double MASS_MU=0.1057;

  l1.SetXYZM(mutmp1.Px()*dataSFa, 

		 mutmp1.Py()*dataSFa,

		 mutmp1.Pz()*dataSFa,

		 MASS_MU); 

      l2.SetXYZM(mutmp2.Px()*dataSFb,

		 mutmp2.Py()*dataSFb,

		 mutmp2.Pz()*dataSFb,

		 MASS_MU); 
pair = l1+l2;
corrMass =  pair.M();





/////////////////////////////////////////
      if(/*Pair_mass[0]*/corrMass<110.)
	continue;
      if(/*MuonCand_pt[0]<50 || MuonCand_pt[1]<50*/pt0cor<50. || pt1cor<50.)
	continue;
      if(MuonCand_istight[0]!=1 || MuonCand_istight[1]!=1)
	continue;
      if(MuonCand_charge[0] == MuonCand_charge[1])
	continue;
      if(1-fabs(Pair_dphi[0])/3.14159 > 0.009)
	continue;
      if(Pair_extratracks0p5mm[0] > 0)
	continue;

/////////////////////



FillInfo fillInfo;
unsigned int ret = fillInfoCollection.FindByRun(Run, fillInfo);
const auto alignment_it = alignmentCollection.find(fillInfo.alignmentTag);
//	alignments = &alignment_it->second;
  //                      aligned=1;



/////////////////////////////////////////////////
double shx_l=0;
double shx_r=0;
 unsigned int rpDecId = 23;
 shx_l = alignment_it->second.Apply(rpDecId,1,0);
  rpDecId = 123;
 shx_r = alignment_it->second.Apply(rpDecId,1,0);
double shx_l_s=0;
double shx_r_s=0;
  rpDecId = 3;
 shx_l_s = alignment_it->second.Apply(rpDecId,1,0);
  rpDecId = 103;
 shx_r_s = alignment_it->second.Apply(rpDecId,1,0);



double shy_l=0;
double shy_r=0;
// unsigned int
 rpDecId = 23;
 shy_l = alignment_it->second.Apply(rpDecId,0,1);
  rpDecId = 123;
 shy_r = alignment_it->second.Apply(rpDecId,0,1);
double shy_l_s=0;
double shy_r_s=0;
  rpDecId = 3;
 shy_l_s = alignment_it->second.Apply(rpDecId,0,1);
  rpDecId = 103;
 shy_r_s = alignment_it->second.Apply(rpDecId,0,1);




///////////////////////


      hacop->Fill(1-fabs(Pair_dphi[0])/3.14159);

      Float_t mumuxisol1 = (1.0/13000.0)*((pt0cor*TMath::Exp(MuonCand_eta[0]) + pt1cor*TMath::Exp(MuonCand_eta[1])));
      Float_t mumuxisol2 = (1.0/13000.0)*((pt0cor*TMath::Exp(-1*MuonCand_eta[0]) + pt1cor*TMath::Exp(-1*MuonCand_eta[1])));


      Float_t protxi45 = 0.0;
      Float_t protxi56 = 0.0;
      Float_t mumuxi45 = 0.0;
      Float_t mumuxi56 = 0.0; double prot_x_p=0;
      Int_t ntrk45 = 0;
      Int_t ntrk56 = 0;

double prot_x_l=0;
double prot_x_r=0;

      for(Int_t p = 0; p < nLocalProtCand; p++)
	{//if(ProtCand_ismultirp[p])continue;
	 




if(LocalProtCand_rpid[p] == 23)
	    { prot_x_l=LocalProtCand_x[p]+shx_l/1000.;
	      ntrk45++;
	    }
	  if(LocalProtCand_rpid[p] == 123)
	    { prot_x_r=LocalProtCand_x[p]+shx_r/1000.;
	      ntrk56++;
	    }



/*
          if(ProtCand_arm[p] == 0 && ProtCand_ismultirp[p] == 1)
            {
              protxi45 = ProtCand_xi[p];
              hres45mult->Fill(1 - (protxi45/mumuxisol1));
	      hressummult->Fill(1 - (protxi45/mumuxisol1));
            }

	  if(ProtCand_arm[p] == 1 && ProtCand_ismultirp[p] == 1)
	    {
              protxi56 = ProtCand_xi[p];
              hres56mult->Fill(1 - (protxi56/mumuxisol2));
              hressummult->Fill(1 - (protxi56/mumuxisol2));
	    }
*/
	}


/////////////////////////////////////////////////////////

if(preTS2==1){
if(ntrk45==1 || ntrk56==1){angle=getXangle(Run,  LumiSection , "./inp/xangle_tillTS2_stableonly_cleanup.csv");}
}


if(preTS2==0){
if(ntrk45==1 || ntrk56==1){angle=getXangle(Run,  LumiSection , "./inp/xangle_afterTS2_cleanup.csv");}
}


             if(ntrk45 ==1)
	    { prot_x_p=prot_x_l;

	      protxi45 = fabs(m_s_x_to_xi_120_l->Eval( prot_x_p) + ( ((120.-angle)/(120.-140.))*(m_s_x_to_xi_140_l->Eval( prot_x_p )-m_s_x_to_xi_120_l->Eval( prot_x_p ) )) ); //ProtCand_xi[p];
	      //protxi45 = ProtCand_xi[p];
	      hcorr45->Fill(protxi45,mumuxisol1);
	      hres45->Fill(1 - (protxi45/mumuxisol1));
	      hressum->Fill(1 - (protxi45/mumuxisol1));
              graph[0]->SetPoint(graph[0]->GetN(), protxi45 , mumuxisol1);
              if(mumuxisol1 >= 0.02 && mumuxisol1 < 0.04)
                hbin145->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.04 && mumuxisol1 < 0.06)
                hbin245->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.07)
                hbin345->Fill(1 - (protxi45/mumuxisol1));
                   
	      ntrk45++;
	    }
	  if(ntrk56 ==1)
	    { prot_x_p=prot_x_r;


	      protxi56 = fabs(m_s_x_to_xi_120_r->Eval( prot_x_p) + ( ((120.-angle)/(120.-140.))*(m_s_x_to_xi_140_r->Eval( prot_x_p )-m_s_x_to_xi_120_r->Eval( prot_x_p ) )) ); //ProtCand_xi[p];
	      
              hcorr56->Fill(protxi56,mumuxisol2);
	      hres56->Fill(1 - (protxi56/mumuxisol2));
	      hressum->Fill(1 - (protxi56/mumuxisol2));


              graph[1]->SetPoint(graph[1]->GetN(), protxi56 , mumuxisol2);

	      if(mumuxisol2 >= 0.02 && mumuxisol2 < 0.04)
		hbin156->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.04 && mumuxisol2 < 0.07)
                hbin256->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.07)
                hbin356->Fill(1 - (protxi56/mumuxisol2));

	      ntrk56++;
	    }

   }//loop

   TFile *fx = new TFile("MoreDimuons2017AllWithMultSingleTrack_BCD_single_RC.root","RECREATE");
graph[0]->SetName("S45_xiRP_xiMUMU");
// graph[i]->Draw("AP");
 graph[0]->Write();
graph[1]->SetName("S56_xiRP_xiMUMU");
// graph[i]->Draw("AP");
 graph[1]->Write();


   hacop->Write();
   hcorr45->Write();
   hres45->Write();
   hcorr56->Write();
   hres56->Write();
   hressum->Write();
   hres45mult->Write();
   hres56mult->Write();
   hressummult->Write();

   hbin145->Write();
   hbin245->Write();
   hbin345->Write();

   hbin156->Write();
   hbin256->Write();
   hbin356->Write();

   fx->Write();
}


int run()
 {
 Dimuons2017Macro m;
 m.Loop();
        
 return 0;
 }
          





int getXangle(int run, int lumi, const char* filename)
 {
 TString drun;
 TString dfill;
 TString dlumi;
 TString temp;

int Xangle=-1;

 TString runs;runs.Form("%d", run);
 TString lumis;lumis.Form("%d", lumi);


 ifstream F;
 F.open((const char*)filename);
 
 //F.clear();
 //F.seekg(0, ios::beg);
 int counter=0;
  if(F){ 
  while (!F.eof())
    {


     F>>drun;
     F>>dfill;
     F>>dlumi;
     F>>temp;
     F>>Xangle;
//if(runfill==runlumi && )

if( runs == drun &&  lumis==dlumi )
  {
  //  cout << "Read: "<< run<<" " << lumi<<"    -----> "<<Xangle << endl;   
break; 
    //return Xangle;
  }
//    prevrunfill=runfill;
//    prevXangle=Xangle;

    }
   }//endif
  else cout << "[!] gerXangle():Error reading from file" << endl;

  if(F.eof())
   {
    cout << "[getXangle() warning:] No Xangle data for this run found!" <<endl;
    F.close();
    return -1;
   }

  else  {F.close();
         return Xangle;
        }  
     
 }



