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
#include "TGraph2D.h"
#include "TMath.h"
#include "TF2.h"
#include "TF1.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include <cmath>

void LoopOnRun(int runno, int lslow, int lshigh, TTree *t1, TTree *t2, ofstream& outputfile)
{
  const int NUM1 =100000 /*t1->GetEntries()*/;
  const int NUM2 = t2->GetEntries();
  
  //  cout<<"#data = "<<NUM2<<"\t #MC = "<<NUM1<<endl;
  
  //Define variables  

  Int_t var_tracks1[1], var_calo1_e1[1], var_calo1_e3[1],var_calo1_e5[1], var_nCalo1[1], var_ID1[3500], var_vertex1Track[5], var_tracksQCD1[1], var_vertex1Ndf[5];
  Double_t var_sumcalo1[1], var_vertex1Z[5], var_e1[3500], var_eta1[3500], var_time1[3500], HF_plus_time1[1], HF_minus_time1[1], var_eEM1[3500], var_eHAD1[3500], var_vertex1X[5],var_vertex1Y[5], var_vertex1Chi2[5]/*, var_vertex1Ndf[5]*/;
  Int_t techBit2[1][128], techBit1[1][128], var_genProcess[1];
  Int_t var_tracks2[1], var_calo2_e1[1], var_calo2_e3[1],var_calo2_e5[1], var_nCalo2[1], var_ID2[3500], var_vertex2Track[5], var_Run[1], var_LS[1], var_vertex2Ndf[5], var_tracks2QCD[1];
  Double_t var_sumcalo2[1], var_vertex2Z[5], var_e2[3500], var_eta2[3500], var_time2[3500], HF_plus_time2[1], HF_minus_time2[1], var_eEM2[3500], var_eHAD2[3500], var_vertex2X[5],var_vertex2Y[5], var_vertex2Chi2[5];
  Int_t var_nvertex1[1], var_nvertex2[1], var_bx[1];
  Int_t var_hltMinBiasBSC1[1],var_hltMinBiasBSC2[1];

  t1->SetBranchAddress("nExtraCaloTowersE5hf",var_calo1_e5);
  t1->SetBranchAddress("nExtraTrackCand",var_tracks1);
  t1->SetBranchAddress("nTrackCandPassQCDCuts",var_tracksQCD1);
  t1->SetBranchAddress("SumCalo_e",var_sumcalo1);
  t1->SetBranchAddress("VertexCand_z",var_vertex1Z);
  t1->SetBranchAddress("VertexCand_x",var_vertex1X);
  t1->SetBranchAddress("VertexCand_y",var_vertex1Y);
  t1->SetBranchAddress("VertexCand_tracks",var_vertex1Track);
  t1->SetBranchAddress("VertexCand_chi2",var_vertex1Chi2);
  t1->SetBranchAddress("VertexCand_ndof",var_vertex1Ndf);
  t1->SetBranchAddress("nCaloCand",var_nCalo1);
  t1->SetBranchAddress("CaloTower_ID",var_ID1);
  t1->SetBranchAddress("CaloTower_e",var_e1);
  t1->SetBranchAddress("CaloTower_eta",var_eta1);
  t1->SetBranchAddress("CaloTower_t",var_time1);
  t1->SetBranchAddress("L1TechnicalTriggers",techBit1);
  t1->SetBranchAddress("HF_plus_time",HF_plus_time1);
  t1->SetBranchAddress("HF_minus_time",HF_minus_time1);
  t1->SetBranchAddress("CaloTower_emE",var_eEM1);
  t1->SetBranchAddress("CaloTower_hadE",var_eHAD1);
  t1->SetBranchAddress("nVertexCand",var_nvertex1);
  t1->SetBranchAddress("GenProcessId",var_genProcess);
  t1->SetBranchAddress("HLTMinBiasBSC",var_hltMinBiasBSC1);


  t2->SetBranchAddress("nExtraCaloTowersE5hf",var_calo2_e5);
  t2->SetBranchAddress("nExtraTrackCand",var_tracks2);  
  t2->SetBranchAddress("nTrackCandPassQCDCuts",var_tracks2QCD);
  t2->SetBranchAddress("SumCalo_e",var_sumcalo2);
  t2->SetBranchAddress("VertexCand_z",var_vertex2Z);
  t2->SetBranchAddress("VertexCand_x",var_vertex2X);
  t2->SetBranchAddress("VertexCand_y",var_vertex2Y);
  t2->SetBranchAddress("VertexCand_tracks",var_vertex2Track);
  t2->SetBranchAddress("VertexCand_chi2",var_vertex2Chi2);
  t2->SetBranchAddress("VertexCand_ndof",var_vertex2Ndf);
  t2->SetBranchAddress("nCaloCand",var_nCalo2);
  t2->SetBranchAddress("CaloTower_ID",var_ID2);
  t2->SetBranchAddress("CaloTower_e",var_e2);
  t2->SetBranchAddress("CaloTower_eta",var_eta2);
  t2->SetBranchAddress("CaloTower_t",var_time2);
  t2->SetBranchAddress("HF_plus_time",HF_plus_time2);
  t2->SetBranchAddress("HF_minus_time",HF_minus_time2);
  t2->SetBranchAddress("L1TechnicalTriggers",techBit2);
  t2->SetBranchAddress("CaloTower_emE",var_eEM2);
  t2->SetBranchAddress("CaloTower_hadE",var_eHAD2);
  t2->SetBranchAddress("nVertexCand",var_nvertex2);
  t2->SetBranchAddress("Run",var_Run);
  t2->SetBranchAddress("LumiSection",var_LS);
  t2->SetBranchAddress("BX",var_bx);
  t2->SetBranchAddress("HLTMinBiasBSC",var_hltMinBiasBSC2);

  TH1F *h1 = new TH1F("h1","h1",200,-20,20);

  double energy_threshold(1.);
  int veto_mc_sd(0) , veto_mc_dd(0), veto_mc_nd(0), veto_data(0), veto_all_mc(0);
  int select_data(0), select_mc_sd(0) , select_mc_dd(0), select_mc_nd(0), select_all_mc(0);
  int trigger_data(0), trigger_mc_sd(0) , trigger_mc_dd(0), trigger_mc_nd(0), trigger_all_mc(0);
  int all_mc_sd(0) , all_mc_dd(0), all_mc_nd(0), all_data(0), all_mc(0);

  cout<<"Starting luminosity counting for Run " << runno << ", LumiSections " << lslow << " to " << lshigh << ""<<endl;
  cout << "Getting MC efficiencies" << endl;
  for(Int_t i = 0;i < NUM1;i++){
        t1->GetEntry(i);
	all_mc++;
         if(var_genProcess[0]==92||var_genProcess[0]==93) all_mc_sd++;
         else if(var_genProcess[0]==94) all_mc_dd++;
         else all_mc_nd++;

        bool passed(false);
//        for(int v=0; v<var_nvertex1[0]; v++){
//                if(var_vertex1Track[v]>0 && var_vertex1Z[v]<14.1161 && var_vertex1Z[v]>-14.09539) passed=true;
//        }
//	if((techBit1[0][40]==1 || techBit1[0][41]==1)){   
	if(var_hltMinBiasBSC1[0] == 1){
	  trigger_all_mc++;
	  if(var_genProcess[0]==92||var_genProcess[0]==93) trigger_mc_sd++;
	  else if(var_genProcess[0]==94) trigger_mc_dd++;
	  else trigger_mc_nd++;
	  
	   if(var_vertex1Track[0]>0 && var_vertex1Z[0]<14.1161 && var_vertex1Z[0]>-14.09539) passed=true;
	   if(var_vertex1Track[0]>2 && TMath::Prob(var_vertex1Chi2[0],var_vertex1Ndf[0])<=0.01) passed=false;

	   if(passed==true){
	     select_all_mc++;
		if(var_genProcess[0]==92||var_genProcess[0]==93) select_mc_sd++;
                else if(var_genProcess[0]==94) select_mc_dd++;
                else select_mc_nd++;

		if((var_tracksQCD1[0]>5 || var_calo1_e5[0]>5)){
		  veto_all_mc++;
	                if(var_genProcess[0]==92||var_genProcess[0]==93) veto_mc_sd++;
	                else if(var_genProcess[0]==94) veto_mc_dd++;
	                else veto_mc_nd++;
		}
	   }
	}
  }
  //cout<<"MC efficiencies:"<<endl;
  //cout<<"n single diff = "<<all_mc_sd<<"\t +trigger = "<<(double)(trigger_mc_sd*100)/all_mc_sd<<"%\t +vertex = "<<(double)(select_mc_sd*100)/all_mc_sd<<"%\t +veto-diff = "<<(double)(veto_mc_sd*100)/all_mc_sd<<"%"<<endl;
  //cout<<"n double diff = "<<all_mc_dd<<"\t +trigger = "<<(double)(trigger_mc_dd*100)/all_mc_dd<<"%\t +vertex = "<<(double)(select_mc_dd*100)/all_mc_dd<<"%\t +veto-diff = "<<(double)(veto_mc_dd*100)/all_mc_dd<<"%"<<endl;
  //cout<<"n non diff    = "<<all_mc_nd<<"\t +trigger = "<<(double)(trigger_mc_nd*100)/all_mc_nd<<"%\t +vertex = "<<(double)(select_mc_nd*100)/all_mc_nd<<"%\t +veto-diff = "<<(double)(veto_mc_nd*100)/all_mc_nd<<"%"<<endl;

  cout << "Fitting vertex z distribution" << endl;
  int previous_run(0);
  //  cout << "Num2 = " << NUM2 << endl;
  // First get the width of the vertex z distribution
  for(Int_t i = 0 ;i < NUM2; i++){
    t2->GetEntry(i);
    if(var_Run[0]==runno && var_LS[0]>=lslow && var_LS[0]<=lshigh){
      if(var_vertex2Track[0]>0) {
	if(!(var_vertex2Track[0]>2 && TMath::Prob(var_vertex2Chi2[0],var_vertex2Ndf[0])<=0.001)){
	  h1->Fill(var_vertex2Z[0]);
	}
      }      
    }
  }

  TF1 *fitf = new TF1("f1","gaus",-20,20);
  //  TF1 *fitf = (TF1*)h1->GetFunction("gaus");
  h1->Fit(fitf);
  double vertexzmean = fitf->GetParameter(1);
  double vertexzsigma = fitf->GetParameter(2);
  double vertexzlow = vertexzmean - (2.0 * vertexzsigma);
  double vertexzhigh = vertexzmean + (2.0 * vertexzsigma);

  cout << "Counting events in data" << endl;
  for(Int_t i = 0 ;i < NUM2; i++){
     t2->GetEntry(i);
     //        if(var_Run[0]!=previous_run) cout<<"i = "<<i<<"\t Run# "<<var_Run[0]<<endl;
        previous_run=var_Run[0];

     bool passed(false);
     if(var_Run[0]==runno && var_LS[0]>=lslow && var_LS[0]<=lshigh){
	all_data++;
//	for(int v=0; v<var_nvertex2[0]; v++){
//		if(var_vertex2Track[v]>0 && var_vertex2Z[v]<vertexzhigh && var_vertex2Z[v]>vertexzlow) passed=true;
//	}
//        if((techBit2[0][40]==1||techBit2[0][41]==1) && techBit2[0][0]==0/*(var_bx[0]==51 || var_bx[0]==2724)*/){
        if((var_hltMinBiasBSC2[0] == 1) && (techBit2[0][0]==1)){
	  trigger_data++;
	  if(var_vertex2Track[0]>0 && var_vertex2Z[0]<vertexzhigh && var_vertex2Z[0]>vertexzlow) passed=true;
	  if(var_vertex2Track[0]>2 && TMath::Prob(var_vertex2Chi2[0],var_vertex2Ndf[0])<=0.001) passed=false;

	  if(passed==true){
	    select_data++;
	    if(var_tracks2QCD[0]>5 || var_calo2_e5[0]>5) veto_data++;
	  }
	}
     }
  }


  double select_all_mcpermb = 52.1 * (double)(select_all_mc)/all_mc;
  double veto_all_mcpermb = 52.1 * (double)(veto_all_mc)/all_mc;
  
  cout<<endl<<endl;
  cout << "*************************************************************************" << endl;
  cout<<"Run " << runno << " (LumiSections " << lslow << "-" << lshigh << "):"<<endl;
  cout<<"\tMC efficiencies:"<<endl;
  cout<<"\t\tn single diff = "<<all_mc_sd<<"\t +trigger = "<<(double)(trigger_mc_sd*100)/all_mc_sd<<"%\t +vertex = "<<(double)(select_mc_sd*100)/all_mc_sd<<"%\t+veto-diff = "<<(double)(veto_mc_sd*100)/all_mc_sd<<"%"<<endl;
  cout<<"\t\tn double diff = "<<all_mc_dd<<"\t +trigger = "<<(double)(trigger_mc_dd*100)/all_mc_dd<<"%\t +vertex = "<<(double)(select_mc_dd*100)/all_mc_dd<<"%\t+veto-diff = "<<(double)(veto_mc_dd*100)/all_mc_dd<<"%"<<endl;
  cout<<"\t\tn non diff    = "<<all_mc_nd<<"\t +trigger = "<<(double)(trigger_mc_nd*100)/all_mc_nd<<"%\t +vertex = "<<(double)(select_mc_nd*100)/all_mc_nd<<"%\t+veto-diff = "<<(double)(veto_mc_nd*100)/all_mc_nd<<"%"<<endl;
  cout<<"\t\tn total    = "<<all_mc<<"\t +trigger = "<<(double)(trigger_all_mc*100)/all_mc<<"%\t +vertex = "<<(double) (select_all_mc*100)/all_mc<<"%\t+veto-diff = "<<(double)(veto_all_mc*100)/all_mc<<"%"<<endl;
  cout<<"\t\tTotal events per mb-1 (no veto): " << select_all_mcpermb << endl;
  cout<<"\t\tTotal events per mb-1 (veto): " << veto_all_mcpermb << endl;
  
  cout << "\tData:" << endl;
  cout<<"\t\tn data = "<<all_data<<"\t +trigger = "<<trigger_data<<"\t +vertex = "<<select_data<<"\t +veto-diff = "<<veto_data<<endl;
  double lumi = (double) select_data / select_all_mcpermb; 
  double lumiveto = (double) veto_data / veto_all_mcpermb;
 
  cout << "\t\tL (no veto) = " << lumi << " mb-1" << endl;
  cout << "\t\tL (veto) = " << lumiveto << " mb-1" << endl;
  cout << "*************************************************************************" << endl;

  outputfile << "| " << runno << " | " << lslow << "-" << lshigh
	     << " | " << select_data << " | " << select_all_mcpermb << " | " << lumi 
	     << " | " << veto_data << " | " << veto_all_mcpermb << " | " << lumiveto << " |\n";  

}

void MinbiasLuminosities()
{
  gROOT->Reset();
  //definition des fichiers + Tree
  TFile *f1 = new TFile("/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_5/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/runs_rereco_all/minbiasmc.335patch2rerunpixtrgnewvtxcuts.900gev.root"); // mc
  TTree *t1 = f1->Get("ntp1");
  TChain *t2 = new TChain("ntp1");
  t2->Add("/storage/data/cms/store/user/schul/out_rereco14Dec_*");

  cout << "Data entries = " << t2->GetEntries() << endl
       << "MC entries = " << t1->GetEntries() << endl;

  TString fil = ("LumiTable.twiki");
  ofstream outFile(fil.Data());
  if (!outFile){cout<<"Error opening output file"<< endl;}
  outFile.setf(ios::floatfield,ios::fixed);
  outFile<<setprecision(2);
  outFile << "| *Run* | *LumiSections* | *Data events* | *MC events per mb-1* | *Lumi* | *Data events (veto)* | *MC events per mb-1 (veto)* | *Lumi (veto)* |\n";  

  /*
   * Start analyzing runs and LumiSection ranges for 2009 data. 
   * Runs with no good (PhysicsEnabled) LumiSections are commented out.
   */

  //  LoopOnRun(123586, 0, 0, t1, t2, outFile);
  LoopOnRun(123587, 5, 0, t1, t2, outFile);
  LoopOnRun(123591, 71, 72, t1, t2, outFile);
  LoopOnRun(123592, 2, 14, t1, t2, outFile);
  LoopOnRun(123596, 2, 144, t1, t2, outFile);
  //  LoopOnRun(123600, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123603, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123606, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123608, 0, 0, t1, t2, outFile);
  LoopOnRun(123615, 70, 91, t1, t2, outFile);
  LoopOnRun(123732, 56, 112, t1, t2, outFile);
  //  LoopOnRun(123734, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123737, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123740, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123791, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123801, 0, 0, t1, t2, outFile);
  LoopOnRun(123815, 7, 16, t1, t2, outFile);
  LoopOnRun(123818, 2, 42, t1, t2, outFile);
  //  LoopOnRun(123893, 0, 0, t1, t2, outFile);
  LoopOnRun(123906, 17, 28, t1, t2, outFile);
  LoopOnRun(123908, 2, 13, t1, t2, outFile);
  LoopOnRun(123909, 12, 29, t1, t2, outFile);
  //  LoopOnRun(123910, 0, 0, t1, t2, outFile);
  //  LoopOnRun(123967, 0, 0, t1, t2, outFile);
  LoopOnRun(123970, 7, 20, t1, t2, outFile);
  //  LoopOnRun(123976, 0, 0, t1, t2, outFile);
  LoopOnRun(123977, 1, 68, t1, t2, outFile);
  LoopOnRun(123978, 2, 25, t1, t2, outFile);
  LoopOnRun(123985, 1, 85, t1, t2, outFile);
  LoopOnRun(123987, 1, 21, t1, t2, outFile);
  //  LoopOnRun(123997, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124006, 0, 0, t1, t2, outFile);
  LoopOnRun(124009, 1, 68, t1, t2, outFile);
  //  LoopOnRun(124017, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124019, 0, 0, t1, t2, outFile);
  LoopOnRun(124020, 12, 94, t1, t2, outFile);
  LoopOnRun(124022, 66, 179, t1, t2, outFile);
  LoopOnRun(124023, 38, 96, t1, t2, outFile);
  LoopOnRun(124024, 2, 83, t1, t2, outFile);
  LoopOnRun(124025, 3, 13, t1, t2, outFile);
  //  LoopOnRun(124026, 0, 0, t1, t2, outFile);
  LoopOnRun(124027, 23, 39, t1, t2, outFile);
  //  LoopOnRun(124028, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124029, 0, 0, t1, t2, outFile);
  LoopOnRun(124030, 1, 32, t1, t2, outFile);
  //  LoopOnRun(124108, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124112, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124115, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124119, 0, 0, t1, t2, outFile);
  LoopOnRun(124120, 1, 59, t1, t2, outFile);
  //  LoopOnRun(124228, 0, 0, t1, t2, outFile);
  LoopOnRun(124230, 23, 68, t1, t2, outFile);
  //  LoopOnRun(124239, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124240, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124265, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124267, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124270, 0, 0, t1, t2, outFile);
  LoopOnRun(124275, 3, 63, t1, t2, outFile);
  //  LoopOnRun(124301, 0, 0, t1, t2, outFile);
  //  LoopOnRun(124310, 0, 0, t1, t2, outFile);
}
