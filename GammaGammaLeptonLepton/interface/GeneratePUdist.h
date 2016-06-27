#ifndef GENERATEPUDIST_H
#define GENERATEPUDIST_H

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h> 
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include <FWCore/ParameterSet/interface/ParameterDescriptionNode.h> 
#include <FWCore/Framework/interface/Event.h>

#include "FWCore/Utilities/interface/InputTag.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TFile.h>
#include <TH1.h>
//
// class declaration
//

class GeneratePUdist : public edm::EDAnalyzer {
 public:
  explicit GeneratePUdist(const edm::ParameterSet&);
  ~GeneratePUdist();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  
  edm::InputTag recVertexLabel;
  
  TH1D *TNPUInTime_;
  TH1D *TNPUTrue_;
  TH1D *TNVTX_;
  TH1D *TNPVIcount_;
  TH1D *TPU_;
};

#endif
