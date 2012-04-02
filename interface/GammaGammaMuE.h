#ifndef DiffractiveForwardAnalysis_GammaGammaMuE
#define DiffractiveForwardAnalysis_GammaGammaMuE

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h> 
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include <FWCore/ParameterSet/interface/ParameterDescriptionNode.h> 
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/Utilities/interface/InputTag.h" 

#include "DataFormats/Common/interface/TriggerResults.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" 
#include "FWCore/Common/interface/TriggerNames.h"

//#include "MuonAnalysis/TagAndProbe/interface/MuonPerformanceReadback.h" 
//#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h"  

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <fstream>
#include <string>

class GammaGammaMuE : public edm::EDAnalyzer {
 public:
  explicit GammaGammaMuE(const edm::ParameterSet&);
  ~GammaGammaMuE();
  
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  edm::InputTag recTrackLabel;
  edm::InputTag recVertexLabel;
  edm::InputTag theGLBMuonLabel;
  edm::InputTag thePixelGsfELabel;
  edm::InputTag theJetLabel;
  edm::InputTag theMetLabel;
  edm::InputTag thePhotonLabel;
  edm::InputTag theCaloTowLabel;
  edm::InputTag recCastorTowerLabel;   
  edm::InputTag recZDCRecHitsLabel;
  edm::InputTag recCastorRecHitsLabel;
  std::string mcPileupFile;
  std::string mcPileupPath;
  std::string dataPileupFile;
  std::string dataPileupPath;
  std::string hltMenuLabel;

  double mudptmax;
  double mudphimin;
  double drisocalo; 
  double minmuevtxd;
  bool keepsamesign;

  std::string rootfilename;

  TFile *thefile;
  TTree *thetree;
  std::ofstream outdebug;

  int BX;
  int Run;
  int LumiSection;
  int EventNum;
  double AvgInstDelLumi;  
  double BunchInstLumi[3]; 
 
  int L1TechnicalTriggers[128]; 

  int nEvt;
  int nMuonCand;
  int MUONMAX;// used to set maximum of arrays
  double MuonCand_px[10];
  double MuonCand_py[10];
  double MuonCand_pz[10];
  double MuonCand_vtxx[10]; 
  double MuonCand_vtxy[10]; 
  double MuonCand_vtxz[10]; 
  double MuonCand_p[10];
  double MuonCand_eta[10];
  double MuonCand_pt[10];
  double MuonCand_phi[10];
  double MuonCand_e[10];
  double MuonCandTrack_p[100];
  double MuonCand_efficiency[10];
  int MuonCand_charge[10];
  int MuonCand_tmlsloosemuonid[10];
  int MuonCand_tmlsOptLowPtloosemuonid[10];
  int MuonCand_tm2dloosemuid[10];
  int MuonCand_arbmuid[10];
  int MuonCand_tmlsAngloosemuonid[10];
  int MuonCand_tmlsAngtightmuonid[10]; 
  int MuonCand_tmosAngloosemuonid[10]; 
  int MuonCand_tmosAngtightmuonid[10]; 
  int MuonCand_gmPromptTight[10];
  int MuonCand_isglobal[10];
  int MuonCand_istracker[10];
  int MuonCand_isstandalone[10];
  double MuonCand_ecalisor3[10]; 
  double MuonCand_hcalisor3[10]; 
  double MuonCand_trkisor3[10];
  double MuonCand_hoisor3[10]; 
  double MuonCand_ecalisor5[10];  
  double MuonCand_hcalisor5[10];  
  double MuonCand_trkisor5[10]; 
  double MuonCand_hoisor5[10];
  double MuonCand_timein[10]; 	 
  double MuonCand_timeout[10]; 	 
  double MuonCand_timeouterr[10]; 	 
  double MuonCand_timeinerr[10]; 	 
  int MuonCand_validtrackhits[10];
  int MuonCand_validhits[10];
  int MuonCand_validpixelhits[10];
  int MuonCand_validmuonhits[10];
  int MuonCand_matches[10];
  double MuonCand_normchi2[10];
  double MuonCand_normtrackchi2[10];
  double MuonCand_dB[10];
  int MuonCand_nlayers[10]; 
  int MuonCand_tightID[10]; 
  int MuonCand_PF[10];

  int nEleCand;
  int ELEMAX;// used to set maximum of arrays
  double EleCand_px[10];
  double EleCand_py[10];
  double EleCand_pz[10];
  double EleCand_p[10];
  double EleCand_e[10];
  double EleCand_et[10];
  double EleCand_phi[10];
  double EleCand_eta[10];
  double EleCandTrack_p[100];
  double EleCandTrack_pt[100]; 
  double EleCandTrack_eta[100]; 
  double EleCandTrack_phi[100];
  double EleCandTrack_vtxz[100]; 
  int EleCand_charge[10];
  int EleCand_looseid[10];
  double EleCand_likelihoodid[10];
  int EleCand_robustid[10];
  double EleCand_vtxx[10];  
  double EleCand_vtxy[10];  
  double EleCand_vtxz[10];  
  double EleCand_deltaPhi[10];
  double EleCand_deltaEta[10];
  double EleCand_HoverE[10];
  double EleCand_trackiso[10];
  double EleCand_ecaliso[10];
  double EleCand_hcaliso[10];
  double EleCand_sigmaIetaIeta[10];
  double EleCand_convDist[10];
  double EleCand_convDcot[10];
  int EleCand_ecalDriven[10]; 
  int EleCand_wp80[10];

  int MuEPairCand[2];

  int nHLTMu10Ele10MuonCand;
  double HLT_Mu10Ele10_MuonCand_pt[10];
  double HLT_Mu10Ele10_MuonCand_eta[10];
  double HLT_Mu10Ele10_MuonCand_phi[10];
  int HLT_Mu10Ele10_MuonCand_charge[10];

  int nHLTMu8Ele17MuonCand; 
  double HLT_Mu8Ele17_MuonCand_pt[10]; 
  double HLT_Mu8Ele17_MuonCand_eta[10]; 
  double HLT_Mu8Ele17_MuonCand_phi[10]; 
  int HLT_Mu8Ele17_MuonCand_charge[10]; 

  int nHLTMu17Ele8MuonCand;  
  double HLT_Mu17Ele8_MuonCand_pt[10];  
  double HLT_Mu17Ele8_MuonCand_eta[10];  
  double HLT_Mu17Ele8_MuonCand_phi[10];  
  int HLT_Mu17Ele8_MuonCand_charge[10];  

  double MuE_mass;
  double MuE_dphi;
  double MuE_dpt;
  double MuE_pt;
  double MuE_phi;
  double MuE_3Dangle;
  double MuE_Kalmanvtxx;
  double MuE_Kalmanvtxy;
  double MuE_Kalmanvtxz;
  double MuE_KalmanvtxT;
  double MuE_Kalmanvtxchi2dof;
  int MuE_Kalmanvtxisvalid;
  int MuE_extratracks1mm;
  int MuE_extratracks2mm; 
  int MuE_extratracks3mm;
  int MuE_extratracks4mm; 
  int MuE_extratracks5mm;
  int MuE_extratracks1cm;
  int MuE_extratracks3cm;
  int MuE_extratracks5cm;
  int MuE_extratracks10cm;

  int nJetCand;
  int JETMAX;// used to set maximum of arrays
  double JetCand_px[30];
  double JetCand_py[30];
  double JetCand_pz[30];
  double JetCand_e[30];
  double JetCand_eta[30];
  double JetCand_phi[30];
  double HighestJet_e;
  double HighestJet_eta;
  double HighestJet_phi;
  double SumJet_e;

  int HitInZDC; 
  int HitInCastor; 

  int nGenMuonCand;
  int GENMUONMAX;
  double GenMuonCand_px[10]; 
  double GenMuonCand_py[10]; 
  double GenMuonCand_pz[10]; 
  double GenMuonCand_pt[10];  
  double GenMuonCand_eta[10];  
  int nGenEleCand; 
  int GENELEMAX; 
  double GenEleCand_px[10];  
  double GenEleCand_py[10];  
  double GenEleCand_pz[10];  
  double GenEleCand_pt[10];   
  double GenEleCand_eta[10];   

  double GenMuE_eta; 
  double GenMuE_pt; 



  double Etmiss;
  double Etmiss_phi;
  double Etmiss_x;
  double Etmiss_y;
  double Etmiss_z;
  double Etmiss_significance;

  int nCaloCand;

  double CaloTower_e[1000];
  double CaloTower_et[1000];
  double CaloTower_eta[1000];
  double CaloTower_phi[1000];
  double CaloTower_dr[1000];
  double CaloTower_emE[1000];
  double CaloTower_hadE[1000];
  double CaloTower_outE[1000];
  int CaloTower_ID[1000];
  double CaloTower_x[1000];
  double CaloTower_y[1000];
  double CaloTower_z[1000];
  double CaloTower_t[1000];
  int CaloTower_badhcalcells[1000];
  int CaloTower_problemhcalcells[1000];
  int CaloTower_badecalcells[1000];
  int CaloTower_problemecalcells[1000];

  double HighestCaloTower_e;
  double HighestCaloTower_eta;
  double HighestCaloTower_phi;
  double HighestCaloTower_dr;
  double HighestEtCaloTower_et;
  double HighestEtCaloTower_eta;
  double HighestEtCaloTower_phi;
  double HighestEtCaloTower_dr;
  double SumCalo_e;

  int nCastorTowerCand;   
  int nCastorTowerCandE3;
  double CastorTower_e[1000];   
  double CastorTower_eta[1000];    
  double CastorTower_phi[1000];   
  double CastorTower_emratio[1000];   
  double HighestCastorTowerFwd_e;   
  double HighestCastorTowerBwd_e;   
  double SumCastorFwd_e; 
  double SumCastorBwd_e; 

  int nZDChitCand;
  int ZDChit_section[5000];
  double ZDChit_energy[5000];
  double ZDChit_time[5000];
  int ZDChit_side[5000];
  double ZDCsumEMplus;
  double ZDCsumHADplus;
  double ZDCsumEMminus;
  double ZDCsumHADminus;

  double CASTORsumRecHitsE;

  int nPrimVertexCand;
  double PrimVertexCand_x[40];
  double PrimVertexCand_y[40];
  double PrimVertexCand_z[40];
  int PrimVertexCand_tracks[40];
  double PrimVertexCand_chi2[40];
  double PrimVertexCand_ndof[40];
  int PrimVertexCand_mueTwoTracks[40];
  int PrimVertexCand_mueExactlyTwoTracks[40];
  int PrimVertexCand_mueTwoTracksMuIndex[40];
  int PrimVertexCand_mueTwoTracksEleIndex[40]; 

    
  int nTrackCand;
  int nQualityTrackCand;
  int TRACKMAX;
  double TrackCand_purity[500];
  int TrackCand_nhits[500];
  double TrackCand_chi2[500];
  double TrackCand_ndof[500];
  double TrackCand_px[500];
  double TrackCand_py[500];
  double TrackCand_pz[500];
  double TrackCand_p[500];
  double TrackCand_eta[500];
  double TrackCand_pt[500];
  double TrackCand_phi[500];
  double TrackCand_vtxdxyz[500];
  double TrackCand_vtxT[500];
  double TrackCand_vtxZ[500];
  double TrackCand_X[500];
  double TrackCand_Y[500];
  double TrackCand_Z[500];
  int TrackCand_charge[500];
  double ClosestExtraTrack_vtxdxyz;
  double ClosestHighPurityExtraTrack_vtxdxyz;

  int nPFPhotonCand;  
  int PHOTONMAX;
  double PFPhotonCand_pt[50];
  double PFPhotonCand_eta[50];
  double PFPhotonCand_phi[50];
  double PFPhotonCand_drtrue[50];

  double evweight;
  
  int HLT_Mu8Ele17;
  int HLT_Mu10Ele10;
  int HLT_Mu17Ele8;
  int HLT_Mu8Ele17_Prescl; 
  int HLT_Mu10Ele10_Prescl; 
  int HLT_Mu17Ele8_Prescl; 

  double LowPt_pt[10];
  double LowPt_eta[10];


  double nTruePUafterPUWeight;
  double nTruePUafterPUWeightBXM1;
  double nTruePUafterPUWeightBXP1;
  double nTruePUafterPUWeightBX0;
  double Weight3D;


  int nTruePUforPUWeight, nTruePUforPUWeightBXM1, nTruePUforPUWeightBXP1, nTruePUforPUWeightBX0;
  double PUWeightTrue;
  HLTConfigProvider hltConfig_;  

  AcceptanceTableHelper helper420beam1;   
  AcceptanceTableHelper helper420beam2;   
  AcceptanceTableHelper helper220beam1;   
  AcceptanceTableHelper helper220beam2;   
  AcceptanceTableHelper helper420a220beam1;   
  AcceptanceTableHelper helper420a220beam2;   

  edm::TriggerNames trigNames ;
};
#endif

/*  LocalWords:  vtxz
 */
