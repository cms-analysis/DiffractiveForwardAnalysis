#ifndef DiffractiveForwardAnalysis_GammaGammaMuMuMC
#define DiffractiveForwardAnalysis_GammaGammaMuMuMC

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/InputTag.h"
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

class GammaGammaMuMuMC : public edm::EDAnalyzer {
 public:
  explicit GammaGammaMuMuMC(const edm::ParameterSet&);
  ~GammaGammaMuMuMC();
  
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  std::string rootfilename;

  TFile *thefile;
  TTree *thetree;

  int nEvt;

  int nMCPar;
  int MCPARMAX;// used to set maximum of arrays
  int MCPar_status[50];
  double MCPar_px[50];
  double MCPar_py[50];
  double MCPar_pz[50];
  double MCPar_e[50];
  double MCPar_phi[50];
  double MCPar_eta[50];
  double MCPar_mass[50];
  int MCPar_pdgid[50];

  int nGenMuonCand;
  int MUONMAX;// used to set maximum of arrays
  double GenMuonCand_px[4];
  double GenMuonCand_py[4];
  double GenMuonCand_pz[4];
  double GenMuonCand_p[4];
  double GenMuonCand_eta[4];
  double GenMuonCand_pt[4];
  double GenMuonCand_phi[4];
  double GenMuonCand_e[4];
  int GenMuonCand_charge[4];
  double GenMuMu_mass;
  double GenMuMu_dphi;

  int nGenProtCand;
  int PROTMAX;// used to set maximum of arrays
  double GenProtCand_px[30];
  double GenProtCand_py[30];
  double GenProtCand_pz[30];
  double GenProtCand_e[30];

  double eventWeight;

};
#endif
