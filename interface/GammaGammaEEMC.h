#ifndef DiffractiveForwardAnalysis_GammaGammaEEMC
#define DiffractiveForwardAnalysis_GammaGammaEEMC

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

class GammaGammaEEMC : public edm::EDAnalyzer {
 public:
  explicit GammaGammaEEMC(const edm::ParameterSet&);
  ~GammaGammaEEMC();
  
  
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

  int nGenEleCand;
  int ELEMAX;// used to set maximum of arrays
  double GenEleCand_px[4];
  double GenEleCand_py[4];
  double GenEleCand_pz[4];
  double GenEleCand_p[4];
  double GenEleCand_e[4];
  double GenEleCand_pt[4];
  double GenEleCand_phi[4];
  double GenEleCand_eta[4];
  double GenEleCand_charge[4];
  double GenElEl_mass;
  double GenElEl_dphi;

  int nGenProtCand;
  int PROTMAX;// used to set maximum of arrays
  double GenProtCand_px[30];
  double GenProtCand_py[30];
  double GenProtCand_pz[30];
  double GenProtCand_e[30];

  double eventWeight;

};
#endif
