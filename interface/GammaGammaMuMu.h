#ifndef DiffractiveForwardAnalysis_GammaGammaMuMu
#define DiffractiveForwardAnalysis_GammaGammaMuMu

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

class GammaGammaMuMu : public edm::EDAnalyzer {
 public:
  explicit GammaGammaMuMu(const edm::ParameterSet&);
  ~GammaGammaMuMu();
  
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  edm::InputTag recTrackLabel;
  edm::InputTag theGLBMuonLabel;
  edm::InputTag thePixelGsfELabel;
  edm::InputTag theJetLabel;
  edm::InputTag theMetLabel;
  edm::InputTag thePhotonLabel;
  edm::InputTag theCaloTowLabel;

  double mudptmax;
  double mudphimin;
  int njetsmax;
  double highestjetemax;
  double sumjetemax;

  std::string rootfilename;

  TFile *thefile;
  TTree *thetree;

  int nMuonCand;
  int MUONMAX;// used to set maximum of arrays
  double MuonCand_px[4];
  double MuonCand_py[4];
  double MuonCand_pz[4];
  double MuonCand_p[4];
  double MuonCand_eta[4];
  double MuonCand_pt[4];
  double MuonCand_phi[4];
  double MuonCand_e[4];
  int MuonCand_charge[4];
  double MuMu_mass;
  double MuMu_dphi;

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

  double Etmiss;

  int nCaloCand;
  double CaloTower_e[1000];
  double CaloTower_eta[1000];
  double CaloTower_phi[1000];
  double HighestCaloTower_e;
  double HighestCaloTower_eta;
  double HighestCaloTower_phi;
  double SumCalo_e;

  int nTrackCand;
  double TrackCand_px[100];
  double TrackCand_py[100];
  double TrackCand_pz[100];
  double TrackCand_p[100];
  double TrackCand_eta[100];
  double TrackCand_pt[100];
  double TrackCand_phi[100];
  int TrackCand_charge[100];

  double eventWeight;

};
#endif
