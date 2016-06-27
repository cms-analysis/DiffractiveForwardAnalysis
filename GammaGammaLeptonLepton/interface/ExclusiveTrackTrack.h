#ifndef DiffractiveForwardAnalysis_ExclusiveTrackTrack
#define DiffractiveForwardAnalysis_ExclusiveTrackTrack

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/Utilities/interface/InputTag.h" 

#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Common/interface/TriggerNames.h" 

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLorentzVector.h>

class ExclusiveTrackTrack : public edm::EDAnalyzer {
 public:
  explicit ExclusiveTrackTrack(const edm::ParameterSet&);
  ~ExclusiveTrackTrack();
  
  
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  edm::InputTag recTrackLabel;
  edm::InputTag theCaloTowLabel;
  edm::InputTag recCastorTowerLabel;
  edm::InputTag recZDCRecHitsLabel;
  edm::InputTag recCastorRecHitsLabel;

  double drisocalo;
  std::string rootfilename;
  bool fillallmc;
  bool requiretwotracks;

  TFile *thefile;
  TTree *thetree;

  int GenProcessId;
  int GenHasRho;

  int nTrackCand;
  int TRACKMAX;// used to set maximum of arrays
  double TrackCand_px[100];
  double TrackCand_py[100];
  double TrackCand_pz[100];
  double TrackCand_p[100];
  double TrackCand_e[100];
  double TrackCand_pt[100];
  double TrackCand_phi[100];
  double TrackCand_eta[100];
  double TrackCand_chi2[100];
  double TrackCand_ndof[100];
  int TrackCand_purity[100]; 
  int TrackCand_charge[100];
  //You can access dE/dx Estimation of your track with: 
  double TrackCand_dEdx[100];
  int TrackCand_nSaturated[100];
  int TrackCand_nMeasurements[100];

  double TrTr_mass;
  double TrTr_dphi;
  double TrTr_dpt;
  double TrTr_pt;
  double TrTr_eta;
  double TrTr_mass_khyp;
  double TrTr_pt_khyp;
  double TrTr_eta_khyp;

  double TrTr_Kalmanvtxx; 
  double TrTr_Kalmanvtxy; 
  double TrTr_Kalmanvtxz; 
  double TrTr_Kalmanvtxchi2dof; 
  int TrTr_Kalmanvtxisvalid; 

  int nCaloCand;
  int nExtraCaloTowersE0pt6eb, nExtraCaloTowersE2pt45ee, nExtraCaloTowersE1pt25hb, nExtraCaloTowersE1pt9he, nExtraCaloTowersE4pt5hfp, nExtraCaloTowersE4pt0hfm; 
  int nExtraCaloTowersE1, nExtraCaloTowersE2, nExtraCaloTowersE3, nExtraCaloTowersE4, nExtraCaloTowersE5, nExtraCaloTowersE6, nExtraCaloTowersE7, nExtraCaloTowersE8, nExtraCaloTowersE9;  
  int nExtraCaloTowersEt0pt1, nExtraCaloTowersEt0pt2, nExtraCaloTowersEt0pt5, nExtraCaloTowersEt1, nExtraCaloTowersEt2, nExtraCaloTowersEt3, nExtraCaloTowersEt4;  

  double CaloTower_e[1000];
  double CaloTower_et[1000];
  double CaloTower_eta[1000];
  double CaloTower_phi[1000];
  double CaloTower_dr[1000];
  double CaloTower_emE[1000]; 
  double CaloTower_hadE[1000]; 
  double CaloTower_outE[1000]; 
  int CaloTower_ID[1000]; 
  double HighestCaloTower_e;
  double HighestCaloTower_eta;
  double HighestCaloTower_phi;
  double HighestCaloTower_dr;
  double HighestEtCaloTower_et;
  double HighestEtCaloTower_eta;
  double HighestEtCaloTower_phi;
  double HighestEtCaloTower_dr;
  double SumCalo_e;
  double SumHFPlus_e;
  double SumHFMinus_e;

  int nCastorTowerCand;
  double CastorTower_e[1000];
  double CastorTower_eta[1000];
  double CastorTower_phi[1000];
  double CastorTower_emratio[1000];
  double HighestCastorTowerFwd_e;
  double HighestCastorTowerBwd_e;
  double SumCastorFwd_e;
  double SumCastorBwd_e;

  int nZDChitCand;
  int ZDChit_section[500];
  double ZDChit_energy[500];
  double ZDChit_time[500];
  int ZDChit_side[500];
  double ZDCsumEMplus;
  double ZDCsumHADplus;
  double ZDCsumEMminus;
  double ZDCsumHADminus;
  double CASTORsumRecHitsE;

  int nVertexCand;
  double VertexCand_x[10];
  double VertexCand_y[10];
  double VertexCand_z[10];
  int VertexCand_tracks[10];
  double VertexCand_chi2[10];
  double VertexCand_ndof[10];

  int nMCPar; 
  int MCPARMAX;// used to set maximum of arrays 
  int MCPar_status[500]; 
  double MCPar_px[500]; 
  double MCPar_py[500]; 
  double MCPar_pz[500]; 
  double MCPar_phi[500]; 
  double MCPar_eta[500]; 
  double MCPar_mass[500]; 
  int MCPar_pdgid[500]; 

  int BX;
  int Run;
  int LumiSection;

  int L1TechnicalTriggers[128];
  int HLTZeroBiasPixelSingleTrack;
  int HLTZeroBias;
  int HLT_L1_BscMinBiasOR_BptxPlusORMinus;
  int HLTPhysicsDeclared;

  edm::TriggerNames trigNames ;
};
#endif
