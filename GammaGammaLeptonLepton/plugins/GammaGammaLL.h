#ifndef DiffractiveForwardAnalysis_GammaGammaLL_h
#define DiffractiveForwardAnalysis_GammaGammaLL_h

// system include files
#include <fstream>
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// L1 collections
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"

// HLT information
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// Generator level collection
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// Pileup
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "DiffractiveForwardAnalysis/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Electrons collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//// tweaked electron ID
//#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ElectronID.h"

// Photons collection
#include "DataFormats/PatCandidates/interface/Photon.h"

// Particle flow collection
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// Jets/MET collection
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

// CT-PPS objects
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/HLTMatcher.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// LHC fill information
//#include "DataFormats/Common/interface/ConditionsInEdm.h" // L1 method
//#include "CondFormats/RunInfo/interface/FillInfo.h"
//#include "CondFormats/DataRecord/interface/FillInfoRcd.h" // db method

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"

#include <map>

typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

//
// class declaration
//

class GammaGammaLL : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit GammaGammaLL( const edm::ParameterSet& );
      ~GammaGammaLL();

      static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

   private:
      virtual void beginJob() ;
      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const& );
      virtual void endRun(edm::Run const&, edm::EventSetup const& );

      virtual void lookAtTriggers( const edm::Event&, const edm::EventSetup& );

      virtual void analyzeMCEventContent( const edm::Event& );
      void fetchElectrons( const edm::Event& );
      void fetchMuons( const edm::Event& );
      void fetchPhotons( const edm::Event& );
      void fetchProtons( const edm::Event& );
      void fetchJets( const edm::Event& );
      void fetchVertices( const edm::Event& );

      void newVertexInfoRetrieval( const edm::Event& );
      bool newTracksInfoRetrieval( int, int );

      // ----------member data ---------------------------
      ggll::TreeType leptonsType_;

      TTree* tree_;
      ggll::AnalysisEvent evt_;

      bool fetchMuons_, fetchElectrons_, fetchProtons_;
      bool foundPairInEvent_;

      // Input tags
      std::string hltMenuLabel_;
      std::vector<std::string> triggersList_;

      edm::EDGetTokenT<pat::TriggerEvent> triggerEventToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > recoVertexToken_;
      edm::EDGetTokenT<edm::View<reco::Track> > recoTrackToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
      edm::EDGetTokenT<edm::View<pat::Muon> > muonToken_;
      edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;
      /*edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_, eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_, phoTightIdMapToken_;*/
      edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
      edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
      edm::EDGetTokenT<edm::View<pat::MET> > metToken_;
      edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack> > totemRPHitToken_;
      edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;

      bool runOnMC_, printCandidates_;
      double minPtMC_, minEtaMC_;
      double sqrts_;
      unsigned int maxExTrkVtx_;

      // Trigger information
      ggll::HLTMatcher hlts_;
      HLTConfigProvider hltConfig_;
      HLTPrescaleProvider hltPrescale_;

      // E/gamma identification
      edm::ParameterSet eleIdLabelSet_;
      std::string eleMediumIdLabel_, eleTightIdLabel_;
      edm::ParameterSet phoIdLabelSet_;
      std::string phoMediumIdLabel_, phoTightIdLabel_;

      // Pileup information
      edm::LumiReWeighting lumiWeights_;
      std::string mcPileupFile_, dataPileupFile_;
      std::string mcPileupPath_, dataPileupPath_;

      edm::Handle<edm::View<reco::Track> > trackColl_;
      edm::ESHandle<TransientTrackBuilder> KalVtx_;
      std::map<int,TLorentzVector> muonsMomenta_, electronsMomenta_;
      std::map<unsigned int,reco::TransientTrack> muonTransientTracks_, eleTransientTracks_;

      unsigned int nCandidates_;
};

#endif
