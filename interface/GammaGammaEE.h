#ifndef DiffractiveForwardAnalysis_GammaGammaEE
#define DiffractiveForwardAnalysis_GammaGammaEE

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

class GammaGammaEE : public edm::EDAnalyzer {
 public:
  explicit GammaGammaEE(const edm::ParameterSet&);
  ~GammaGammaEE();
  
  
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

  double eldetmax;
  double eldphimin;
  int njetsmax;
  double highestjetemax;
  double sumjetemax;

  std::string rootfilename;

  TFile *thefile;
  TTree *thetree;

  int nEleCand;
  int ELEMAX;// used to set maximum of arrays
  double EleCand_px[4];
  double EleCand_py[4];
  double EleCand_pz[4];
  double EleCand_p[4];
  double EleCand_e[4];
  double EleCand_et[4];
  double EleCand_phi[4];
  double EleCand_eta[4];
  double EleCand_charge[4];
  double ElEl_mass;
  double ElEl_dphi;

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
  double CaloTower_e[100];
  double CaloTower_eta[100];
  double CaloTower_phi[100];
  double HighestCaloTower_e;
  double HighestCaloTower_eta;
  double HighestCaloTower_phi;
  double SumCalo_e;

  int nTrackCand;
  double TrackCand_px[4];
  double TrackCand_py[4];
  double TrackCand_pz[4];
  double TrackCand_p[4];
  double TrackCand_eta[4];
  double TrackCand_pt[4];
  double TrackCand_phi[4];
  int TrackCand_charge[4];

  double eventWeight;

};
#endif
