#ifndef DiffractiveForwardAnalysis_ExclusiveTrackTrack
#define DiffractiveForwardAnalysis_ExclusiveTrackTrack

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Framework/interface/TriggerNames.h" 

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLorentzVector.h>

class ExclusiveTrackTrack : public edm::EDAnalyzer {
 public:
  explicit ExclusiveTrackTrack(const edm::ParameterSet&);
  ~ExclusiveTrackTrack();
  
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  edm::InputTag recTrackLabel;
  edm::InputTag theCaloTowLabel;
  edm::InputTag recCastorTowerLabel;

  double drisocalo;

  std::string rootfilename;

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
  int TrackCand_charge[100];

  double TrTr_mass;
  double TrTr_dphi;
  double TrTr_dpt;
  double TrTr_pt;
  double TrTr_vtxx; 
  double TrTr_vtxy; 
  double TrTr_vtxz; 
  double TrTr_vtxchi2dof; 
  int TrTr_vtxisvalid; 

  int nCaloCand;
  int nExtraCaloTowersE1, nExtraCaloTowersE2, nExtraCaloTowersE3, nExtraCaloTowersE4, nExtraCaloTowersE5, nExtraCaloTowersE6, nExtraCaloTowersE7, nExtraCaloTowersE8, nExtraCaloTowersE9;  
  int nExtraCaloTowersEt0pt1, nExtraCaloTowersEt0pt2, nExtraCaloTowersEt0pt5, nExtraCaloTowersEt1, nExtraCaloTowersEt2, nExtraCaloTowersEt3, nExtraCaloTowersEt4;  
  int nExtraCaloTowersE0hf, nExtraCaloTowersE1hf, nExtraCaloTowersE2hf, nExtraCaloTowersE3hf, nExtraCaloTowersE4hf, nExtraCaloTowersE5hf;
  int nExtraCaloTowersE0hfp, nExtraCaloTowersE1hfp, nExtraCaloTowersE2hfp, nExtraCaloTowersE3hfp, nExtraCaloTowersE4hfp, nExtraCaloTowersE5hfp;
  int nExtraCaloTowersE0hfm, nExtraCaloTowersE1hfm, nExtraCaloTowersE2hfm, nExtraCaloTowersE3hfm, nExtraCaloTowersE4hfm, nExtraCaloTowersE5hfm;
  int nExtraCaloTowersE1he, nExtraCaloTowersE2he, nExtraCaloTowersE3he;
  int nExtraCaloTowersE2hb, nExtraCaloTowersE3hb, nExtraCaloTowersE4hb;
  double CaloTower_e[1000];
  double CaloTower_et[1000];
  double CaloTower_eta[1000];
  double CaloTower_phi[1000];
  double CaloTower_dr[1000];
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

  int nVertexCand;
  double VertexCand_x[10];
  double VertexCand_y[10];
  double VertexCand_z[10];
  int VertexCand_tracks[10];
  double VertexCand_chi2[10];
  double VertexCand_ndof[10];

  int BX;
  int Run;
  int LumiSection;

  int L1TechnicalTriggers[128];
  int HLTMinBiasPixelSingleTrack;
  int HLTPhysicsDeclared;

  edm::TriggerNames trigNames ;
};
#endif
