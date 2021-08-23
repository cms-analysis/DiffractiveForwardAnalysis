//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 13 14:47:30 2019 by ROOT version 6.12/07
// from TChain ggll_aod/ntp1/
//////////////////////////////////////////////////////////

#ifndef Dimuons2018Macro_h
#define Dimuons2018Macro_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class Dimuons2018Macro {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          LumiSection;
   UInt_t          BX;
   UInt_t          EventNum;
   Float_t        CrossingAngle;
   UInt_t          nHLT;
   Int_t           HLT_Accept[3];   //[nHLT]
   Int_t           HLT_Prescl[3];   //[nHLT]
   vector<string>  *HLT_Name;
   UInt_t          nMuonCand;
   Double_t        MuonCand_pt[10];   //[nMuonCand]
   Double_t        MuonCand_eta[10];   //[nMuonCand]
   Double_t        MuonCand_phi[10];   //[nMuonCand]
   Double_t        MuonCand_e[10];   //[nMuonCand]
   Int_t           MuonCand_charge[10];   //[nMuonCand]
   Double_t        MuonCand_vtxx[10];   //[nMuonCand]
   Double_t        MuonCand_vtxy[10];   //[nMuonCand]
   Double_t        MuonCand_vtxz[10];   //[nMuonCand]
   Double_t        MuonCand_dxy[10];   //[nMuonCand]
   Int_t           MuonCand_nstatseg[10];   //[nMuonCand]
   Int_t           MuonCand_ntrklayers[10];   //[nMuonCand]
   Int_t           MuonCand_npxlhits[10];   //[nMuonCand]
   Int_t           MuonCand_isglobal[10];   //[nMuonCand]
   Int_t           MuonCand_istracker[10];   //[nMuonCand]
   Int_t           MuonCand_isstandalone[10];   //[nMuonCand]
   Int_t           MuonCand_ispfmuon[10];   //[nMuonCand]
   Int_t           MuonCand_istight[10];   //[nMuonCand]
   Int_t           MuonCandTrack_nmuchits[10];   //[nMuonCand]
   Double_t        MuonCandTrack_chisq[10];   //[nMuonCand]
   Double_t        MuonCand_innerTrackPt[10];   //[nMuonCand]
   Double_t        MuonCand_innerTrackEta[10];   //[nMuonCand]
   Double_t        MuonCand_innerTrackPhi[10];   //[nMuonCand]
   Double_t        MuonCand_innerTrackVtxz[10];   //[nMuonCand]
   UInt_t          nPrimVertexCand;
   Double_t        PrimVertexCand_x[97];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_y[97];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_z[97];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_chi2[97];   //[nPrimVertexCand]
   UInt_t          PrimVertexCand_ndof[97];   //[nPrimVertexCand]
   UInt_t          PrimVertexCand_tracks[97];   //[nPrimVertexCand]
   UInt_t          nPair;
   Int_t           Pair_lepton1[10];   //[nPair]
   Int_t           Pair_lepton2[10];   //[nPair]
   Double_t        Pair_mass[10];   //[nPair]
   Double_t        Pair_pt[10];   //[nPair]
   Double_t        Pair_eta[10];   //[nPair]
   Double_t        Pair_phi[10];   //[nPair]
   Double_t        Pair_dpt[10];   //[nPair]
   Double_t        Pair_dphi[10];   //[nPair]
   Double_t        Pair_3Dangle[10];   //[nPair]
   UInt_t          Pair_extratracks0p5mm[10];   //[nPair]
   UInt_t          Pair_extratracks1mm[10];   //[nPair]
   UInt_t          Pair_extratracks2mm[10];   //[nPair]
   UInt_t          Pair_extratracks3mm[10];   //[nPair]
   UInt_t          Pair_extratracks4mm[10];   //[nPair]
   UInt_t          Pair_extratracks5mm[10];   //[nPair]
   UInt_t          Pair_extratracks1cm[10];   //[nPair]
   UInt_t          Pair_extratracks2cm[10];   //[nPair]
   UInt_t          Pair_extratracks3cm[10];   //[nPair]
   UInt_t          Pair_extratracks4cm[10];   //[nPair]
   UInt_t          Pair_extratracks5cm[10];   //[nPair]
   UInt_t          Pair_extratracks10cm[10];   //[nPair]
   Double_t        KalmanVertexCand_x[10];   //[nPair]
   Double_t        KalmanVertexCand_y[10];   //[nPair]
   Double_t        KalmanVertexCand_z[10];   //[nPair]
   UInt_t          nLocalProtCand;
   Double_t        LocalProtCand_x[100];   //[nLocalProtCand]
   Double_t        LocalProtCand_y[100];   //[nLocalProtCand]
   Double_t        LocalProtCand_t[100];   //[nLocalProtCand]
   Double_t        LocalProtCand_xSigma[100];   //[nLocalProtCand]
   Double_t        LocalProtCand_ySigma[100];   //[nLocalProtCand]
   Double_t        LocalProtCand_tSigma[100];   //[nLocalProtCand]
   Int_t           LocalProtCand_arm[100];   //[nLocalProtCand]
   Int_t           LocalProtCand_station[100];   //[nLocalProtCand]
   Int_t           LocalProtCand_pot[100];   //[nLocalProtCand]
   Int_t           LocalProtCand_rpid[100];   //[nLocalProtCand]
   UInt_t          nRecoProtCand;
   Double_t        ProtCand_xi[50];   //[nRecoProtCand]
   Double_t        ProtCand_t[50];   //[nRecoProtCand]
   Double_t        ProtCand_ThX[50];   //[nRecoProtCand]
   Double_t        ProtCand_ThY[50];   //[nRecoProtCand]
   Double_t        ProtCand_ystar[50];   //[nRecoProtCand]                                                                                                               
   Int_t           ProtCand_rpid[50];   //[nRecoProtCand]
   Int_t           ProtCand_arm[50];   //[nRecoProtCand]
   Int_t           ProtCand_ismultirp[50];   //[nRecoProtCand]
   Int_t           ProtCand_trackpixshift1[50]; //[nRecoProtCand]
   Double_t        ProtCand_trackx1[50]; //[nRecoProtCand]                                                                                                
   Double_t        ProtCand_tracky1[50]; //[nRecoProtCand]                                                                                                
   Double_t        ProtCand_trackx2[50]; //[nRecoProtCand]                                                                                                
   Double_t        ProtCand_tracky2[50]; //[nRecoProtCand]                                                                                                
   UInt_t          nExtraTracks;
   Int_t           ExtraTrack_pair[2000];   //[nExtraTracks]
   Int_t           ExtraTrack_purity[2000];   //[nExtraTracks]
   UInt_t          ExtraTrack_nhits[2000];   //[nExtraTracks]
   Int_t           ExtraTrack_charge[2000];   //[nExtraTracks]
   UInt_t          ExtraTrack_ndof[2000];   //[nExtraTracks]
   Double_t        ExtraTrack_px[2000];   //[nExtraTracks]
   Double_t        ExtraTrack_py[2000];   //[nExtraTracks]
   UInt_t          nGenMuonCand;
   Double_t        GenMuonCand_pt[10];   //[nMuonCand]                                                                                                                                     
   Double_t        GenMuonCand_eta[10];   //[nMuonCand]                                                                                                                                     

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_BX;   //!
   TBranch        *b_EventNum;   //!
   TBranch        *b_CrossingAngle;  //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT_Accept;   //!
   TBranch        *b_HLT_Prescl;   //!
   TBranch        *b_HLT_Name;   //!
   TBranch        *b_nMuonCand;   //!
   TBranch        *b_MuonCand_pt;   //!
   TBranch        *b_MuonCand_eta;   //!
   TBranch        *b_MuonCand_phi;   //!
   TBranch        *b_MuonCand_e;   //!
   TBranch        *b_MuonCand_charge;   //!
   TBranch        *b_MuonCand_vtxx;   //!
   TBranch        *b_MuonCand_vtxy;   //!
   TBranch        *b_MuonCand_vtxz;   //!
   TBranch        *b_MuonCand_dxy;   //!
   TBranch        *b_MuonCand_nstatseg;   //!
   TBranch        *b_MuonCand_ntrklayers;   //!
   TBranch        *b_MuonCand_npxlhits;   //!
   TBranch        *b_MuonCand_isglobal;   //!
   TBranch        *b_MuonCand_istracker;   //!
   TBranch        *b_MuonCand_isstandalone;   //!
   TBranch        *b_MuonCand_ispfmuon;   //!
   TBranch        *b_MuonCand_istight;   //!
   TBranch        *b_MuonCandTrack_nmuchits;   //!
   TBranch        *b_MuonCandTrack_chisq;   //!
   TBranch        *b_MuonCand_innerTrackPt;   //!
   TBranch        *b_MuonCand_innerTrackEta;   //!
   TBranch        *b_MuonCand_innerTrackPhi;   //!
   TBranch        *b_MuonCand_innerTrackVtxz;   //!
   TBranch        *b_nPrimVertexCand;   //!
   TBranch        *b_PrimVertexCand_x;   //!
   TBranch        *b_PrimVertexCand_y;   //!
   TBranch        *b_PrimVertexCand_z;   //!
   TBranch        *b_PrimVertexCand_chi2;   //!
   TBranch        *b_PrimVertexCand_ndof;   //!
   TBranch        *b_PrimVertexCand_tracks;   //!
   TBranch        *b_nPair;   //!
   TBranch        *b_Pair_lepton1;   //!
   TBranch        *b_Pair_lepton2;   //!
   TBranch        *b_Pair_mass;   //!
   TBranch        *b_Pair_pt;   //!
   TBranch        *b_Pair_eta;   //!
   TBranch        *b_Pair_phi;   //!
   TBranch        *b_Pair_dpt;   //!
   TBranch        *b_Pair_dphi;   //!
   TBranch        *b_Pair_3Dangle;   //!
   TBranch        *b_Pair_extratracks0p5mm;   //!
   TBranch        *b_Pair_extratracks1mm;   //!
   TBranch        *b_Pair_extratracks2mm;   //!
   TBranch        *b_Pair_extratracks3mm;   //!
   TBranch        *b_Pair_extratracks4mm;   //!
   TBranch        *b_Pair_extratracks5mm;   //!
   TBranch        *b_Pair_extratracks1cm;   //!
   TBranch        *b_Pair_extratracks2cm;   //!
   TBranch        *b_Pair_extratracks3cm;   //!
   TBranch        *b_Pair_extratracks4cm;   //!
   TBranch        *b_Pair_extratracks5cm;   //!
   TBranch        *b_Pair_extratracks10cm;   //!
   TBranch        *b_KalmanVertexCand_x;   //!
   TBranch        *b_KalmanVertexCand_y;   //!
   TBranch        *b_KalmanVertexCand_z;   //!
   TBranch        *b_nLocalProtCand;   //!
   TBranch        *b_LocalProtCand_x;   //!
   TBranch        *b_LocalProtCand_y;   //!
   TBranch        *b_LocalProtCand_t;   //!
   TBranch        *b_LocalProtCand_xSigma;   //!
   TBranch        *b_LocalProtCand_ySigma;   //!
   TBranch        *b_LocalProtCand_tSigma;   //!
   TBranch        *b_LocalProtCand_arm;   //!
   TBranch        *b_LocalProtCand_station;   //!
   TBranch        *b_LocalProtCand_pot;   //!
   TBranch        *b_LocalProtCand_rpid;   //!
   TBranch        *b_nRecoProtCand;   //!
   TBranch        *b_ProtCand_xi;   //!
   TBranch        *b_ProtCand_t;   //!
   TBranch        *b_ProtCand_ThX;   //!
   TBranch        *b_ProtCand_ThY;   //!
   TBranch        *b_ProtCand_ystar;   //!                                                                                                                               
   TBranch        *b_ProtCand_rpid;   //!
   TBranch        *b_ProtCand_arm;   //!
   TBranch        *b_ProtCand_ismultirp;   //!
   TBranch        *b_ProtCand_trackpixshift1;   //!
   TBranch        *b_ProtCand_trackx1; //!                                                                                                               
   TBranch        *b_ProtCand_tracky1; //!                                                                                                               
   TBranch        *b_ProtCand_trackx2; //!                                                                                                               
   TBranch        *b_ProtCand_tracky2; //!                                                                                                               
   TBranch        *b_nExtraTracks;   //!
   TBranch        *b_ExtraTrack_pair;   //!
   TBranch        *b_ExtraTrack_purity;   //!
   TBranch        *b_ExtraTrack_nhits;   //!
   TBranch        *b_ExtraTrack_charge;   //!
   TBranch        *b_ExtraTrack_ndof;   //!
   TBranch        *b_ExtraTrack_px;   //!
   TBranch        *b_ExtraTrack_py;   //!
   TBranch        *b_nGenMuonCand;   //!                                                                                                                                                     
   TBranch        *b_GenMuonCand_pt;   //!                                                                                                                                                  
   TBranch        *b_GenMuonCand_eta;   //!                                                                                                               

   // Efficiency correction histograms                                                                                                                   
   TH2F *hpixeffA45210, *hpixeffB145210, *hpixeffB245210, *hpixeffC45210, *hpixeffD145210, *hpixeffD245210;
   TH2F *hpixeffA56210, *hpixeffB156210, *hpixeffB256210, *hpixeffC56210, *hpixeffD156210, *hpixeffD256210;
   TH2F *hpixeffA45220, *hpixeffB145220, *hpixeffB245220, *hpixeffC45220, *hpixeffD145220, *hpixeffD245220;
   TH2F *hpixeffA56220, *hpixeffB156220, *hpixeffB256220, *hpixeffC56220, *hpixeffD156220, *hpixeffD256220;
                         
   // Systematics                                                                                                     
   TFile *fsyst;
   TGraphErrors *grsyst45;
   TGraphErrors *grsyst56;
         

   Dimuons2018Macro(TTree *tree=0);
   virtual ~Dimuons2018Macro();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t multi=0, Int_t mc=0, Int_t sb=1, Int_t yr=2017, Int_t nearfar=0);
   virtual bool     FiducalCuts(Float_t trackx210, Float_t tracky210, Float_t trackx220, Float_t tracky220, Int_t arm, Int_t run);
   virtual Float_t  MultiRPEffCorr(Float_t trackx210, Float_t tracky210, Float_t trackx220, Float_t tracky220, Int_t arm, Int_t run);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Dimuons2018Macro_cxx
Dimuons2018Macro::Dimuons2018Macro(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/2017/ProtonReconstruction_DoubleMu2017F_merge.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/2017/ProtonReconstruction_DoubleMu2017F_merge.root");
      }
      f->GetObject("ggll_aod/ntp1",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ggll_aod/ntp1","");

      // 2017, Ultra-Legacy
      /*
      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2017B_finalUL_merge.root/ggll_aod/ntp1");
      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2017C_finalUL_merge.root/ggll_aod/ntp1");
      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2017D_finalUL_merge.root/ggll_aod/ntp1");
      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2017E_finalUL_merge.root/ggll_aod/ntp1");
      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2017F_finalUL_merge.root/ggll_aod/ntp1");
      */


      //      chain->Add("output_xangleall2017MC.root");
      //      chain->Add("output_xangleall2017MC_noPUprotons.root");

      // Final
      //      chain->Add("../test/output_xangleall2017MCpreTS2.root");
      //      chain->Add("../test/output_xangleall2018MC.root");
      //      chain->Add("../test/output_xangleall2017postTS2MC.root");


      // 2018, final legacy tests
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018A_finalUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018B_finalUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018C_finalUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018D_finalUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("output_xangleall2018MC.root");

      // 2018D, re-MiniAOD tests
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018A_reMiniAODtestUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018B_reMiniAODtestUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018C_reMiniAODtestUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018Dpart1_reMiniAODtestUL_merge.root/ggll_aod/ntp1");
      //      chain->Add("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/FinalUltraLegacy/ProtonReconstruction_DoubleMu2018Dpart2_reMiniAODtestUL_merge.root/ggll_aod/ntp1");

      // final DPS MC
      chain->Add("../reminiaodtests/output_MC_2018.root/ggll_aod/ntp1");

	tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);

   // Systematics                                                                                                     
   // Old version for DPS
   //   fsyst = TFile::Open("/tmp/jjhollar/xi_uncertainty.root");
   fsyst = TFile::Open("reco_charactersitics_version1.root");

   TFile *fpixeff = TFile::Open("pixelEfficiencies_radiation.root");
   hpixeffA45210 = (TH2F *)fpixeff->Get("Pixel/2018/2018A/h45_210_2018A_all_2D");
   hpixeffB145210 = (TH2F *)fpixeff->Get("Pixel/2018/2018B1/h45_210_2018B1_all_2D");
   hpixeffB245210 = (TH2F *)fpixeff->Get("Pixel/2018/2018B2/h45_210_2018B2_all_2D");
   hpixeffC45210 = (TH2F *)fpixeff->Get("Pixel/2018/2018C/h45_210_2018C_all_2D");
   hpixeffD145210 = (TH2F *)fpixeff->Get("Pixel/2018/2018D1/h45_210_2018D1_all_2D");
   hpixeffD245210 = (TH2F *)fpixeff->Get("Pixel/2018/2018D2/h45_210_2018D2_all_2D");

   hpixeffA56210 = (TH2F *)fpixeff->Get("Pixel/2018/2018A/h56_210_2018A_all_2D");
   hpixeffB156210 = (TH2F *)fpixeff->Get("Pixel/2018/2018B1/h56_210_2018B1_all_2D");
   hpixeffB256210 = (TH2F *)fpixeff->Get("Pixel/2018/2018B2/h56_210_2018B2_all_2D");
   hpixeffC56210 = (TH2F *)fpixeff->Get("Pixel/2018/2018C/h56_210_2018C_all_2D");
   hpixeffD156210 = (TH2F *)fpixeff->Get("Pixel/2018/2018D1/h56_210_2018D1_all_2D");
   hpixeffD256210 = (TH2F *)fpixeff->Get("Pixel/2018/2018D2/h56_210_2018D2_all_2D");

   hpixeffA45220 = (TH2F *)fpixeff->Get("Pixel/2018/2018A/h45_220_2018A_all_2D");
   hpixeffB145220 = (TH2F *)fpixeff->Get("Pixel/2018/2018B1/h45_220_2018B1_all_2D");
   hpixeffB245220 = (TH2F *)fpixeff->Get("Pixel/2018/2018B2/h45_220_2018B2_all_2D");
   hpixeffC45220 = (TH2F *)fpixeff->Get("Pixel/2018/2018C/h45_220_2018C_all_2D");
   hpixeffD145220 = (TH2F *)fpixeff->Get("Pixel/2018/2018D1/h45_220_2018D1_all_2D");
   hpixeffD245220 = (TH2F *)fpixeff->Get("Pixel/2018/2018D2/h45_220_2018D2_all_2D");

   hpixeffA56220 = (TH2F *)fpixeff->Get("Pixel/2018/2018A/h56_220_2018A_all_2D");
   hpixeffB156220 = (TH2F *)fpixeff->Get("Pixel/2018/2018B1/h56_220_2018B1_all_2D");
   hpixeffB256220 = (TH2F *)fpixeff->Get("Pixel/2018/2018B2/h56_220_2018B2_all_2D");
   hpixeffC56220 = (TH2F *)fpixeff->Get("Pixel/2018/2018C/h56_220_2018C_all_2D");
   hpixeffD156220 = (TH2F *)fpixeff->Get("Pixel/2018/2018D1/h56_220_2018D1_all_2D");
   hpixeffD256220 = (TH2F *)fpixeff->Get("Pixel/2018/2018D2/h56_220_2018D2_all_2D");
}

Dimuons2018Macro::~Dimuons2018Macro()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Dimuons2018Macro::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Dimuons2018Macro::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Dimuons2018Macro::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   HLT_Name = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   fChain->SetBranchAddress("BX", &BX, &b_BX);
   fChain->SetBranchAddress("EventNum", &EventNum, &b_EventNum);
   fChain->SetBranchAddress("CrossingAngle", &CrossingAngle, &b_CrossingAngle);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT_Accept", HLT_Accept, &b_HLT_Accept);
   fChain->SetBranchAddress("HLT_Prescl", HLT_Prescl, &b_HLT_Prescl);
   fChain->SetBranchAddress("HLT_Name", &HLT_Name, &b_HLT_Name);
   fChain->SetBranchAddress("nMuonCand", &nMuonCand, &b_nMuonCand);
   fChain->SetBranchAddress("MuonCand_pt", MuonCand_pt, &b_MuonCand_pt);
   fChain->SetBranchAddress("MuonCand_eta", MuonCand_eta, &b_MuonCand_eta);
   fChain->SetBranchAddress("MuonCand_phi", MuonCand_phi, &b_MuonCand_phi);
   fChain->SetBranchAddress("MuonCand_e", MuonCand_e, &b_MuonCand_e);
   fChain->SetBranchAddress("MuonCand_charge", MuonCand_charge, &b_MuonCand_charge);
   fChain->SetBranchAddress("MuonCand_vtxx", MuonCand_vtxx, &b_MuonCand_vtxx);
   fChain->SetBranchAddress("MuonCand_vtxy", MuonCand_vtxy, &b_MuonCand_vtxy);
   fChain->SetBranchAddress("MuonCand_vtxz", MuonCand_vtxz, &b_MuonCand_vtxz);
   fChain->SetBranchAddress("MuonCand_dxy", MuonCand_dxy, &b_MuonCand_dxy);
   fChain->SetBranchAddress("MuonCand_nstatseg", MuonCand_nstatseg, &b_MuonCand_nstatseg);
   fChain->SetBranchAddress("MuonCand_ntrklayers", MuonCand_ntrklayers, &b_MuonCand_ntrklayers);
   fChain->SetBranchAddress("MuonCand_npxlhits", MuonCand_npxlhits, &b_MuonCand_npxlhits);
   fChain->SetBranchAddress("MuonCand_isglobal", MuonCand_isglobal, &b_MuonCand_isglobal);
   fChain->SetBranchAddress("MuonCand_istracker", MuonCand_istracker, &b_MuonCand_istracker);
   fChain->SetBranchAddress("MuonCand_isstandalone", MuonCand_isstandalone, &b_MuonCand_isstandalone);
   fChain->SetBranchAddress("MuonCand_ispfmuon", MuonCand_ispfmuon, &b_MuonCand_ispfmuon);
   fChain->SetBranchAddress("MuonCand_istight", MuonCand_istight, &b_MuonCand_istight);
   fChain->SetBranchAddress("MuonCandTrack_nmuchits", MuonCandTrack_nmuchits, &b_MuonCandTrack_nmuchits);
   fChain->SetBranchAddress("MuonCandTrack_chisq", MuonCandTrack_chisq, &b_MuonCandTrack_chisq);
   fChain->SetBranchAddress("MuonCand_innerTrackPt", MuonCand_innerTrackPt, &b_MuonCand_innerTrackPt);
   fChain->SetBranchAddress("MuonCand_innerTrackEta", MuonCand_innerTrackEta, &b_MuonCand_innerTrackEta);
   fChain->SetBranchAddress("MuonCand_innerTrackPhi", MuonCand_innerTrackPhi, &b_MuonCand_innerTrackPhi);
   fChain->SetBranchAddress("MuonCand_innerTrackVtxz", MuonCand_innerTrackVtxz, &b_MuonCand_innerTrackVtxz);
   fChain->SetBranchAddress("nPrimVertexCand", &nPrimVertexCand, &b_nPrimVertexCand);
   fChain->SetBranchAddress("PrimVertexCand_x", PrimVertexCand_x, &b_PrimVertexCand_x);
   fChain->SetBranchAddress("PrimVertexCand_y", PrimVertexCand_y, &b_PrimVertexCand_y);
   fChain->SetBranchAddress("PrimVertexCand_z", PrimVertexCand_z, &b_PrimVertexCand_z);
   fChain->SetBranchAddress("PrimVertexCand_chi2", PrimVertexCand_chi2, &b_PrimVertexCand_chi2);
   fChain->SetBranchAddress("PrimVertexCand_ndof", PrimVertexCand_ndof, &b_PrimVertexCand_ndof);
   fChain->SetBranchAddress("PrimVertexCand_tracks", PrimVertexCand_tracks, &b_PrimVertexCand_tracks);
   fChain->SetBranchAddress("nPair", &nPair, &b_nPair);
   fChain->SetBranchAddress("Pair_lepton1", Pair_lepton1, &b_Pair_lepton1);
   fChain->SetBranchAddress("Pair_lepton2", Pair_lepton2, &b_Pair_lepton2);
   fChain->SetBranchAddress("Pair_mass", Pair_mass, &b_Pair_mass);
   fChain->SetBranchAddress("Pair_pt", Pair_pt, &b_Pair_pt);
   fChain->SetBranchAddress("Pair_eta", Pair_eta, &b_Pair_eta);
   fChain->SetBranchAddress("Pair_phi", Pair_phi, &b_Pair_phi);
   fChain->SetBranchAddress("Pair_dpt", Pair_dpt, &b_Pair_dpt);
   fChain->SetBranchAddress("Pair_dphi", Pair_dphi, &b_Pair_dphi);
   fChain->SetBranchAddress("Pair_3Dangle", Pair_3Dangle, &b_Pair_3Dangle);
   fChain->SetBranchAddress("Pair_extratracks0p5mm", Pair_extratracks0p5mm, &b_Pair_extratracks0p5mm);
   fChain->SetBranchAddress("Pair_extratracks1mm", Pair_extratracks1mm, &b_Pair_extratracks1mm);
   fChain->SetBranchAddress("Pair_extratracks2mm", Pair_extratracks2mm, &b_Pair_extratracks2mm);
   fChain->SetBranchAddress("Pair_extratracks3mm", Pair_extratracks3mm, &b_Pair_extratracks3mm);
   fChain->SetBranchAddress("Pair_extratracks4mm", Pair_extratracks4mm, &b_Pair_extratracks4mm);
   fChain->SetBranchAddress("Pair_extratracks5mm", Pair_extratracks5mm, &b_Pair_extratracks5mm);
   fChain->SetBranchAddress("Pair_extratracks1cm", Pair_extratracks1cm, &b_Pair_extratracks1cm);
   fChain->SetBranchAddress("Pair_extratracks2cm", Pair_extratracks2cm, &b_Pair_extratracks2cm);
   fChain->SetBranchAddress("Pair_extratracks3cm", Pair_extratracks3cm, &b_Pair_extratracks3cm);
   fChain->SetBranchAddress("Pair_extratracks4cm", Pair_extratracks4cm, &b_Pair_extratracks4cm);
   fChain->SetBranchAddress("Pair_extratracks5cm", Pair_extratracks5cm, &b_Pair_extratracks5cm);
   fChain->SetBranchAddress("Pair_extratracks10cm", Pair_extratracks10cm, &b_Pair_extratracks10cm);
   fChain->SetBranchAddress("KalmanVertexCand_x", KalmanVertexCand_x, &b_KalmanVertexCand_x);
   fChain->SetBranchAddress("KalmanVertexCand_y", KalmanVertexCand_y, &b_KalmanVertexCand_y);
   fChain->SetBranchAddress("KalmanVertexCand_z", KalmanVertexCand_z, &b_KalmanVertexCand_z);
   fChain->SetBranchAddress("nLocalProtCand", &nLocalProtCand, &b_nLocalProtCand);
   fChain->SetBranchAddress("LocalProtCand_x", LocalProtCand_x, &b_LocalProtCand_x);
   fChain->SetBranchAddress("LocalProtCand_y", LocalProtCand_y, &b_LocalProtCand_y);
   fChain->SetBranchAddress("LocalProtCand_t", LocalProtCand_t, &b_LocalProtCand_t);
   fChain->SetBranchAddress("LocalProtCand_xSigma", LocalProtCand_xSigma, &b_LocalProtCand_xSigma);
   fChain->SetBranchAddress("LocalProtCand_ySigma", LocalProtCand_ySigma, &b_LocalProtCand_ySigma);
   fChain->SetBranchAddress("LocalProtCand_tSigma", LocalProtCand_tSigma, &b_LocalProtCand_tSigma);
   fChain->SetBranchAddress("LocalProtCand_arm", LocalProtCand_arm, &b_LocalProtCand_arm);
   fChain->SetBranchAddress("LocalProtCand_station", LocalProtCand_station, &b_LocalProtCand_station);
   fChain->SetBranchAddress("LocalProtCand_pot", LocalProtCand_pot, &b_LocalProtCand_pot);
   fChain->SetBranchAddress("LocalProtCand_rpid", LocalProtCand_rpid, &b_LocalProtCand_rpid);
   fChain->SetBranchAddress("nRecoProtCand", &nRecoProtCand, &b_nRecoProtCand);
   fChain->SetBranchAddress("ProtCand_xi", ProtCand_xi, &b_ProtCand_xi);
   fChain->SetBranchAddress("ProtCand_t", ProtCand_t, &b_ProtCand_t);
   fChain->SetBranchAddress("ProtCand_ThX", ProtCand_ThX, &b_ProtCand_ThX);
   fChain->SetBranchAddress("ProtCand_ThY", ProtCand_ThY, &b_ProtCand_ThY);
   fChain->SetBranchAddress("ProtCand_ystar", ProtCand_ystar, &b_ProtCand_ystar);
   fChain->SetBranchAddress("ProtCand_rpid", ProtCand_rpid, &b_ProtCand_rpid);
   fChain->SetBranchAddress("ProtCand_arm", ProtCand_arm, &b_ProtCand_arm);
   fChain->SetBranchAddress("ProtCand_ismultirp", ProtCand_ismultirp, &b_ProtCand_ismultirp);
   fChain->SetBranchAddress("ProtCand_trackpixshift1", ProtCand_trackpixshift1, &b_ProtCand_trackpixshift1);
   fChain->SetBranchAddress("ProtCand_trackx1", ProtCand_trackx1, &b_ProtCand_trackx1);
   fChain->SetBranchAddress("ProtCand_tracky1", ProtCand_tracky1, &b_ProtCand_tracky1);
   fChain->SetBranchAddress("ProtCand_trackx2", ProtCand_trackx2, &b_ProtCand_trackx2);
   fChain->SetBranchAddress("ProtCand_tracky2", ProtCand_tracky2, &b_ProtCand_tracky2);
   fChain->SetBranchAddress("nExtraTracks", &nExtraTracks, &b_nExtraTracks);
   fChain->SetBranchAddress("ExtraTrack_pair", ExtraTrack_pair, &b_ExtraTrack_pair);
   fChain->SetBranchAddress("ExtraTrack_purity", ExtraTrack_purity, &b_ExtraTrack_purity);
   fChain->SetBranchAddress("ExtraTrack_nhits", ExtraTrack_nhits, &b_ExtraTrack_nhits);
   fChain->SetBranchAddress("ExtraTrack_charge", ExtraTrack_charge, &b_ExtraTrack_charge);
   fChain->SetBranchAddress("ExtraTrack_ndof", ExtraTrack_ndof, &b_ExtraTrack_ndof);
   fChain->SetBranchAddress("ExtraTrack_px", ExtraTrack_px, &b_ExtraTrack_px);
   fChain->SetBranchAddress("ExtraTrack_py", ExtraTrack_py, &b_ExtraTrack_py);
   fChain->SetBranchAddress("nGenMuonCand", &nGenMuonCand, &b_nGenMuonCand);
   fChain->SetBranchAddress("GenMuonCand_pt", GenMuonCand_pt, &b_GenMuonCand_pt);
   fChain->SetBranchAddress("GenMuonCand_eta", GenMuonCand_eta, &b_GenMuonCand_eta);

   Notify();
}

Bool_t Dimuons2018Macro::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Dimuons2018Macro::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Dimuons2018Macro::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Dimuons2018Macro_cxx
