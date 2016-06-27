// -*- C++ -*-
//
// Package:    GeneratePUdist
// Class:      GeneratePUdist
// 
/**\class GeneratePUdist GeneratePUdist.cc DiffractiveForwardAnalysis/GeneratePUdist/src/GeneratePUdist.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme,40 4-B20,+41227671567,
//         Created:  Wed Jan 25 18:52:33 CET 2012
// $Id: GeneratePUdist.cc,v 1.1 2012/04/10 14:23:35 lforthom Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <TFile.h>
#include <TH1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "FWCore/Framework/interface/ESHandle.h" 

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"  
#include "DataFormats/Common/interface/Ref.h"   

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GeneratePUdist.h"


using namespace std;
using namespace edm;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GeneratePUdist::GeneratePUdist(const edm::ParameterSet& pset)

{
  //now do what ever initialization is needed

  // setup file service for PU information
  string filename;


  recVertexLabel  = pset.getParameter<edm::InputTag>("RecoVertexLabel");
  //filename        = pset.getParameter<std::string>("outfilename");


  //cout << "Producing output : Pileup distribution file (" << filename << ")" << endl;
  edm::Service<TFileService> fs;
  // setup histograms
  TNPUInTime_ = fs->make<TH1D>("TNPUInTime","No. in-time pileup interactions",40,0.,40.);
  TNPUTrue_ = fs->make<TH1D>("TNPUTrue","True pileup interactions",40,0.,40.);
  TNPVIcount_ = fs->make<TH1D>("TNPVIcount","No. PVI",40,0.,40.);
  TNVTX_ = fs->make<TH1D>("TNVTX","No. reconstructed vertices",40,0.,40.);
  TPU_ = fs->make<TH1D>("TPU","Debug",40,0.,40.);

}


GeneratePUdist::~GeneratePUdist()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GeneratePUdist::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
   event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
   
   std::vector<PileupSummaryInfo>::const_iterator PVI;

   int npv = 0;
   int ni = 0;
   int npvtrue = -1;
   int npvm1true = 0;
   int npv0true = 0;
   int npvp1true = 0;
   double nm1 = -999.;
   double n0 = -999.;
   double np1 = -999.;
   double npT = -999.;

   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
     int BX = PVI->getBunchCrossing();
     if(abs(BX)<=1) {
       if(BX == -1) {
	 npvm1true++;
	 nm1 = PVI->getPU_NumInteractions();
       }
       if(BX == 0) {
	 npv0true++;
	 npT = PVI->getTrueNumInteractions();
	 n0 = PVI->getPU_NumInteractions();      
       }
       if(BX == 1) {
	 npvp1true++;
	 np1 = PVI->getPU_NumInteractions();
       }
       npv += PVI->getPU_NumInteractions();
       ni++;
     }
     npvtrue = PVI->getTrueNumInteractions();
   }

   edm::Handle< std::vector<reco::Vertex> > vertices_h;
   event.getByLabel(recVertexLabel, vertices_h);
   float NVtx = vertices_h->size();

   TNPVIcount_->Fill(PupInfo->size());
   TNVTX_->Fill(float(NVtx));
   TNPUTrue_->Fill(npT);
   TNPUInTime_->Fill(n0);
   if (npv>=0) TPU_->Fill(float(npv)/float(ni));
}


// ------------ method called once each job just before starting event loop  ------------
void 
GeneratePUdist::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GeneratePUdist::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GeneratePUdist::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GeneratePUdist::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GeneratePUdist::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GeneratePUdist::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GeneratePUdist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GeneratePUdist);
