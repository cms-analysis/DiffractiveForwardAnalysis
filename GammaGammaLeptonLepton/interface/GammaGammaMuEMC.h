#ifndef DiffractiveForwardAnalysis_GammaGammaMuEMC
#define DiffractiveForwardAnalysis_GammaGammaMuEMC

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
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <fstream>
#include <string>

class GammaGammaMuEMC : public edm::EDAnalyzer {
 public:
  explicit GammaGammaMuEMC(const edm::ParameterSet&);
  ~GammaGammaMuEMC();
  
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions); 
  
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  std::string rootfilename;
  bool fillallmc;

  TFile *thefile;
  TTree *thetree;

  int nEvt, nWeird;

  int nMCPar;
  int MCPARMAX;// used to set maximum of arrays
  int MCPar_status[10000];
  double MCPar_px[10000];
  double MCPar_py[10000];
  double MCPar_pz[10000];
  double MCPar_e[10000];
  double MCPar_phi[10000];
  double MCPar_eta[10000];
  double MCPar_mass[10000];
  int MCPar_pdgid[10000];

  int PassPartonLevel_JetVeto;

  int HitInZDC;
  int HitInCastor;

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
  int GenMuonCand_status[4];
  int GenMuonCand_charge[4];

  int nGenEleCand; 
  int ELEMAX;// used to set maximum of arrays 
  double GenEleCand_px[4]; 
  double GenEleCand_py[4]; 
  double GenEleCand_pz[4]; 
  double GenEleCand_p[4]; 
  double GenEleCand_eta[4]; 
  double GenEleCand_pt[4]; 
  double GenEleCand_phi[4]; 
  double GenEleCand_e[4]; 
  int GenEleCand_status[4];
  int GenEleCand_charge[4]; 

  int nGenPhotCand;  
  int PHOTMAX;// used to set maximum of arrays  
  double GenPhotCand_px[4];  
  double GenPhotCand_py[4];  
  double GenPhotCand_pz[4];  
  double GenPhotCand_p[4];  
  double GenPhotCand_eta[4];  
  double GenPhotCand_pt[4];  
  double GenPhotCand_phi[4];  
  double GenPhotCand_e[4];  
  int GenPhotCand_charge[4];  

  double GenGamGam_mass;
  double GenMuE_mass;
  double GenMuE_dphi;
  double GenMuE_pt;

  double GenMET;

  int nGenJets;
  double GenJet_pt[1000];
  double GenJet_eta[1000]; 
  double GenJet_phi[1000]; 
  double GenJet_dR[1000];

  int nGenProtCand;
  int PROTMAX;// used to set maximum of arrays
  double GenProtCand_px[30];
  double GenProtCand_py[30];
  double GenProtCand_pz[30];
  double GenProtCand_e[30];

  double eventWeight;

};
#endif
