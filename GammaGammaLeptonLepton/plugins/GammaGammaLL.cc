// -*- C++ -*-
//
// Package:    GammaGammaLL
// Class:      GammaGammaLL
//
/**\class GammaGammaLL GammaGammaLL.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/plugins/GammaGammaLL.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme,40 4-B20,+41227671567,
//         Created:  Thu Sep 13 15:17:14 CET 2012
// $Id: GammaGammaLL.cc,v 1.3 2013/04/28 08:40:45 lforthom Exp $
//
//

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

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

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

// Particle flow collection
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


// PPS objects
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"


#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/HLTMatcher.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/PrimaryVertexSelector.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// LHC fill information
//#include "DataFormats/Common/interface/ConditionsInEdm.h" // L1 method
//#include "CondFormats/RunInfo/interface/FillInfo.h"
//#include "CondFormats/DataRecord/interface/FillInfoRcd.h" // db method
#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

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

class GammaGammaLL : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources>
{
  public:
    explicit GammaGammaLL( const edm::ParameterSet& );
    ~GammaGammaLL() {}

    static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

  private:
    void beginJob() override;
    void beginRun( const edm::Run&, const edm::EventSetup& ) override;
    void analyze( const edm::Event&, const edm::EventSetup& ) override;
    void endRun( const edm::Run&, const edm::EventSetup& ) override {}
    void endJob() override;

    void lookAtTriggers( const edm::Event&, const edm::EventSetup& );
    void analyzeMCEventContent( const edm::Event& );

    void fetchElectrons( const edm::Event& );
    void fetchMuons( const edm::Event& );
    void fetchProtons( const edm::Event& );
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
    edm::InputTag triggerResults_;
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
    edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
    edm::EDGetTokenT<edm::View<CTPPSLocalTrackLite> > ppsLocalTrackToken_;
    edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
    edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;

    bool runOnMC_, printCandidates_;
    double minPtMC_, minEtaMC_;
    double sqrts_;
    bool saveExtraTracks_;
    std::string year_; 

    // Trigger information
    ggll::HLTMatcher hlts_;
    HLTConfigProvider hltConfig_;
    HLTPrescaleProvider hltPrescale_;

    unsigned int maxExTrkVtx_;
    double minMpair_, maxMpair_;

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

const unsigned int ggll::AnalysisEvent::MAX_ET;

GammaGammaLL::GammaGammaLL( const edm::ParameterSet& iConfig ) :
  tree_( 0 ),
  fetchMuons_( false ), fetchElectrons_( false ),
  fetchProtons_       ( iConfig.getParameter<bool>( "fetchProtons" ) ),
  triggerResults_     ( iConfig.getParameter<edm::InputTag>            ( "triggerResults" ) ),
  triggersList_       ( iConfig.getParameter<std::vector<std::string> >( "triggersList" ) ),
  triggerResultsToken_( consumes<edm::TriggerResults>                  ( iConfig.getParameter<edm::InputTag>( "triggerResults" ) ) ),
  pileupToken_        ( consumes<edm::View<PileupSummaryInfo> >        ( iConfig.getParameter<edm::InputTag>( "pileupInfo" ) ) ),
  recoVertexToken_    ( consumes<edm::View<reco::Vertex> >             ( iConfig.getParameter<edm::InputTag>( "vertexTag" ) ) ),
  recoTrackToken_     ( consumes<edm::View<reco::Track> >              ( iConfig.getParameter<edm::InputTag>( "trackTag" ) ) ),
  muonToken_          ( consumes<edm::View<pat::Muon> >                ( iConfig.getParameter<edm::InputTag>( "muonTag" ) ) ),
  eleToken_           ( consumes<edm::View<pat::Electron> >            ( iConfig.getParameter<edm::InputTag>( "electronTag" ) ) ),
  /*eleMediumIdMapToken_( consumes<edm::ValueMap<bool> >                 ( iConfig.getParameter<edm::InputTag>( "eleMediumIdMap" ) ) ),
  eleTightIdMapToken_ ( consumes<edm::ValueMap<bool> >                 ( iConfig.getParameter<edm::InputTag>( "eleTightIdMap" ) ) ),
  phoMediumIdMapToken_( consumes<edm::ValueMap<bool> >                 ( iConfig.getParameter<edm::InputTag>( "phoMediumIdMap" ) ) ),
  phoTightIdMapToken_ ( consumes<edm::ValueMap<bool> >                 ( iConfig.getParameter<edm::InputTag>( "phoTightIdMap" ) ) ),*/
  fixedGridRhoFastjetAllToken_( consumes<double>                       ( iConfig.getParameter<edm::InputTag>( "fixedGridRhoFastjetAllLabel" ) ) ),
  ppsLocalTrackToken_ ( consumes<edm::View<CTPPSLocalTrackLite> >      ( iConfig.getParameter<edm::InputTag>( "ppsLocalTrackTag" ) ) ),
  recoProtonsSingleRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonSingleRPTag" ) ) ),
  recoProtonsMultiRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonMultiRPTag" ) ) ),
  runOnMC_            ( iConfig.getParameter<bool>( "runOnMC" ) ),
  printCandidates_    ( iConfig.getParameter<bool>( "printCandidates" ) ),
  sqrts_              ( iConfig.getParameter<double>( "sqrtS" ) ),
  saveExtraTracks_    ( iConfig.getParameter<bool>( "saveExtraTracks" ) ),
  year_               ( iConfig.getParameter<std::string>( "year" ) ),
  hlts_               ( triggersList_ ),
  hltPrescale_        ( iConfig, consumesCollector(), *this ),
  // Central selection
  maxExTrkVtx_        ( iConfig.getUntrackedParameter<unsigned int>( "maxExtraTracks", ggll::AnalysisEvent::MAX_ET ) ),
  minMpair_           ( iConfig.getUntrackedParameter<double>( "minMpair", -1. ) ),
  maxMpair_           ( iConfig.getUntrackedParameter<double>( "maxMpair", -1. ) ),
  // Pileup input tags
  mcPileupFile_       ( iConfig.getParameter<std::string>( "mcpufile" ) ),
  dataPileupFile_     ( iConfig.getParameter<std::string>( "datapufile" ) ),
  mcPileupPath_       ( iConfig.getParameter<std::string>( "mcpupath" ) ),
  dataPileupPath_     ( iConfig.getParameter<std::string>( "datapupath" ) ),
  nCandidates_( 0 )
{
  // Generator level
  if ( runOnMC_ ) {
    genToken_ = consumes<edm::View<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>( "genParticleTag" ) );
    minPtMC_ = iConfig.getUntrackedParameter<double>( "MCAcceptPtCut", 10. );
    minEtaMC_ = iConfig.getUntrackedParameter<double>( "MCAcceptEtaCut", 2.5 );
  }

  // Leptons input tags
  const std::string ltype = iConfig.getParameter<std::string>( "leptonsType" );
  if ( ltype == "ElectronMuon" ) {
    leptonsType_ = ggll::ElectronMuon;
    fetchElectrons_ = true;
    fetchMuons_ = true;
  }
  else if ( ltype == "Muon" ) {
    leptonsType_ = ggll::DiMuon;
    fetchMuons_ = true;
  }
  else if ( ltype == "Electron" ) {
    leptonsType_ = ggll::DiElectron;
    fetchElectrons_ = true;
  }
  else throw cms::Exception( "GammaGammaLL" ) << "'LeptonsType' parameter should either be:\n"
                                            << "   * 'ElectronMuon'       (for mixed leptons pair)\n"
                                            << "   * 'Electron' or 'Muon' (for same-flavour leptons)";

  if ( fetchElectrons_ ) {
    // electron identification variables
    const edm::ParameterSet eleIdLabelSet = iConfig.getParameter<edm::ParameterSet>( "eleIdLabels" );
    eleMediumIdLabel_ = eleIdLabelSet.getParameter<edm::InputTag>( "mediumLabel" ).encode();
    eleTightIdLabel_ = eleIdLabelSet.getParameter<edm::InputTag>( "tightLabel" ).encode();
  }

  // Pileup reweighting utilities
  if ( runOnMC_ )
    lumiWeights_ = edm::LumiReWeighting( mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_ );

  // Book the output tree
  usesResource( "TFileService" );
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>( "ntp1", "ntp1" );
}

void
GammaGammaLL::lookAtTriggers( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // Get the trigger information from the event
  evt_.nHLT = triggersList_.size();
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken( triggerResultsToken_, hltResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);

  std::ostringstream os;
  os << "Prescale set: " << hltPrescale_.prescaleSet(iEvent, iSetup) << "\n"
     << "Trigger names: " << std::endl;
  for ( unsigned int i = 0; i < trigNames.size(); ++i ) {
    os << "* " << trigNames.triggerNames().at( i ) << std::endl;

    // ensure trigger matches the interesting ones
    const int trigNum = hlts_.TriggerNum( trigNames.triggerNames().at( i ) );
    if ( trigNum < 0 )
      continue;

    evt_.HLT_Accept[trigNum] = hltResults->accept( i );

    if ( !evt_.HLT_Accept[trigNum] )
      continue;
    // extract prescale value for this path
    if ( !iEvent.isRealData() ) {
      evt_.HLT_Prescl[trigNum] = 1.;
      continue;
    } //FIXME
    int prescale_set = hltPrescale_.prescaleSet( iEvent, iSetup );
    if ( prescale_set >= 0 )
      evt_.HLT_Prescl[trigNum] = hltConfig_.prescaleValue( prescale_set, trigNames.triggerNames().at( i ) );
    /*std::pair<int,int> prescales = hltPrescale_.prescaleValues( iEvent, iSetup, trigNames.triggerNames().at( i ) );
    //std::cout << "trigger path " << trigNames.triggerNames().at( i ) << " has L1/HLT prescales: " << prescales.first << "/" << prescales.second
    //  << "   event r/l/e: " << iEvent.id().run() << "/" << iEvent.luminosityBlock() << "/" << iEvent.id().event() << std::endl;
    evt_.HLT_Prescl[trigNum] = prescales.second;*/
  }
  //std::cout
  LogDebug( "GammaGammaLL" )
    << os.str();
}

// ------------ method called for each event  ------------
void
GammaGammaLL::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // Kalman filtering
  iSetup.get<TransientTrackRecord>().get( "TransientTrackBuilder", KalVtx_ );

  // First initialization of the variables
  evt_.clear();

  evt_.Weight = 1.;

  // Run and BX information
  evt_.BX = iEvent.bunchCrossing();
  evt_.Run = iEvent.id().run();
  evt_.LumiSection = iEvent.luminosityBlock();
  evt_.EventNum = iEvent.id().event();

  // High level trigger information retrieval
  // JH - for Summer17 MC this crashes...
  if ( !runOnMC_ ) {
    lookAtTriggers( iEvent, iSetup);
  }

  LogDebug( "GammaGammaLL" ) << "Passed trigger filtering stage";

  // Crossing angle information from DB - available for all years in legacy re-RECO
  if(year_ == "2018" || year_ == "2017" || year_ == "2016")
    {
      edm::ESHandle<LHCInfo> pSetup;
      const string label = "";
      iSetup.get<LHCInfoRcd>().get(label, pSetup);
      
      // re-initialise algorithm upon crossing-angle change                                                                                            
      const LHCInfo* pInfo = pSetup.product();
      evt_.CrossingAngle = pInfo->crossingAngle();
    }
  else
    evt_.CrossingAngle = -999;

  // Generator level information
  if ( runOnMC_ ) {
    analyzeMCEventContent( iEvent );

    // Pileup information
    edm::Handle<edm::View<PileupSummaryInfo> > pu_info;
    iEvent.getByToken( pileupToken_, pu_info);

    int npv0true = 0;
    for ( unsigned int i = 0; i < pu_info->size(); ++i ) {
      const edm::Ptr<PileupSummaryInfo> PVI = pu_info->ptrAt( i);

      const int beamXing = PVI->getBunchCrossing(),
                npvtrue = PVI->getTrueNumInteractions();

      if ( beamXing == 0) npv0true += npvtrue;
    }

    evt_.Weight = lumiWeights_.weight( npv0true );
    LogDebug( "GammaGammaLL" ) << "Passed Pileup retrieval stage";
  }

  muonsMomenta_.clear();
  electronsMomenta_.clear();
  muonTransientTracks_.clear();
  eleTransientTracks_.clear();

  if ( fetchMuons_ ) fetchMuons( iEvent );
  if ( fetchElectrons_ ) fetchElectrons( iEvent );

  newVertexInfoRetrieval( iEvent );

  if ( !foundPairInEvent_ ) {
    //LogDebug( "GammaGammaLL" ) << "No pair retrieved in event";
    return; // avoid to unpack RP/jet/MET if no dilepton candidate has been found
  }

  if ( fetchProtons_ ) {
    fetchProtons( iEvent );
    if ( evt_.nLocalProtCand < 1 )
      return;
  }

  if ( printCandidates_ )
    std::cout << "Event " << evt_.Run << ":" << evt_.EventNum << " has " << evt_.nPair << " leptons pair(s) candidate(s) (vertex mult. : " << evt_.nPrimVertexCand << " )" << std::endl;

  tree_->Fill();
}

void
GammaGammaLL::analyzeMCEventContent( const edm::Event& iEvent )
{
  edm::Handle<edm::View<reco::GenParticle> > genPartColl;
  iEvent.getByToken( genToken_, genPartColl );

  for ( unsigned int i = 0; i < genPartColl->size(); ++i ) {
    const edm::Ptr<reco::GenParticle> genPart = genPartColl->ptrAt( i );

    if ( !genPart->isPromptFinalState() ) continue;

    // check the particles out of acceptance
    if ( genPart->pt() < minPtMC_ || ( minEtaMC_ != -1. && fabs( genPart->eta() ) > minEtaMC_ ) ) {
      if ( fabs( genPart->pdgId() ) == 13 ) evt_.nGenMuonCandOutOfAccept++;
      if ( fabs( genPart->pdgId() ) == 11 ) evt_.nGenEleCandOutOfAccept++;
      if ( fabs( genPart->pdgId() ) == 22 ) evt_.nGenPhotCandOutOfAccept++;
      if ( genPart->pdgId() != 2212 ) continue; // we keep the forward protons
    }

    // generated outgoing protons
    if ( genPart->pdgId() == 2212 && evt_.nGenProtCand < ggll::AnalysisEvent::MAX_GENPRO ) {
      evt_.GenProtCand_pt[evt_.nGenProtCand] = genPart->pt();
      evt_.GenProtCand_eta[evt_.nGenProtCand] = genPart->eta();
      evt_.GenProtCand_phi[evt_.nGenProtCand] = genPart->phi();
      evt_.GenProtCand_e[evt_.nGenProtCand] = genPart->energy();
      evt_.GenProtCand_status[evt_.nGenProtCand] = genPart->status();
      evt_.nGenProtCand++;
    }

    // generated central dimuon system
    if ( fabs( genPart->pdgId() ) == 13 && evt_.nGenMuonCand < ggll::AnalysisEvent::MAX_GENMU ) {
      evt_.GenMuonCand_pt[evt_.nGenMuonCand] = genPart->pt();
      evt_.GenMuonCand_eta[evt_.nGenMuonCand] = genPart->eta();
      evt_.GenMuonCand_phi[evt_.nGenMuonCand] = genPart->phi();
      evt_.GenMuonCand_e[evt_.nGenMuonCand] = genPart->energy();
      evt_.nGenMuonCand++;
    }

    // generated central dimuon system
    if ( fabs( genPart->pdgId() ) == 11 && evt_.nGenEleCand < ggll::AnalysisEvent::MAX_GENELE ) {
      evt_.GenEleCand_pt[evt_.nGenEleCand] = genPart->pt();
      evt_.GenEleCand_eta[evt_.nGenEleCand] = genPart->eta();
      evt_.GenEleCand_phi[evt_.nGenEleCand] = genPart->phi();
      evt_.GenEleCand_e[evt_.nGenEleCand] = genPart->energy();
      evt_.nGenEleCand++;
    }

    if ( genPart->pdgId() == 2212 && fabs( genPart->pz() ) > 3000. ) {
      // Kinematic quantities computation
      // xi = fractional momentum loss
      /*if ( genPart->pz() > 0. ) xi = 1.-2.*genPart->pz()/sqrts_;
      else xi = 1.+2.*genPart->pz()/sqrts_;
      t = -( std::pow( genPart->pt(), 2 )+std::pow( genPart->mass()*xi, 2 ) )/( 1.-xi );*/
    }
  }

  bool foundGenCandPairInEvent = false;

  TLorentzVector l1, l2;
  // electron+muon
  if ( leptonsType_ == ggll::ElectronMuon && evt_.nGenMuonCand == 1 && evt_.nGenEleCand == 1 ) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE( evt_.GenMuonCand_pt[0], evt_.GenMuonCand_eta[0], evt_.GenMuonCand_phi[0], evt_.GenMuonCand_e[0] );
    l2.SetPtEtaPhiE( evt_.GenEleCand_pt[0], evt_.GenEleCand_eta[0], evt_.GenEleCand_phi[0], evt_.GenEleCand_e[0] );
    foundGenCandPairInEvent = true;
  }
  // dielectron
  else if ( leptonsType_ == ggll::DiElectron && evt_.nGenEleCand == 2 ) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE( evt_.GenEleCand_pt[0], evt_.GenEleCand_eta[0], evt_.GenEleCand_phi[0], evt_.GenEleCand_e[0] );
    l2.SetPtEtaPhiE( evt_.GenEleCand_pt[1], evt_.GenEleCand_eta[1], evt_.GenEleCand_phi[1], evt_.GenEleCand_e[1] );
    foundGenCandPairInEvent = true;
  }
  // dimuon
  else if ( leptonsType_ == ggll::DiMuon && evt_.nGenMuonCand == 2 ) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE( evt_.GenMuonCand_pt[0], evt_.GenMuonCand_eta[0], evt_.GenMuonCand_phi[0], evt_.GenMuonCand_e[0] );
    l2.SetPtEtaPhiE( evt_.GenMuonCand_pt[1], evt_.GenMuonCand_eta[1], evt_.GenMuonCand_phi[1], evt_.GenMuonCand_e[1] );
    foundGenCandPairInEvent = true;
  }
  if ( foundGenCandPairInEvent ) {
    const TLorentzVector pair = l1+l2;
    evt_.GenPair_mass = pair.M();
    evt_.GenPair_pt = pair.Pt();
    evt_.GenPair_eta = pair.Eta();
    evt_.GenPair_phi = pair.Phi();

    double dphi = fabs( l1.Phi()-l2.Phi() );
    // dphi lies in [-pi, pi]
    while ( dphi < -M_PI ) dphi += 2.*M_PI;
    while ( dphi >  M_PI ) dphi -= 2.*M_PI;
    evt_.GenPair_dphi = dphi;

    evt_.GenPair_dpt = fabs( l1.Pt()-l2.Pt() );
    evt_.GenPair_3Dangle = l1.Angle( l2.Vect() )/M_PI;
  }
}

void
GammaGammaLL::fetchMuons( const edm::Event& iEvent )
{
  // Get the muons collection from the event
  // PAT muons
  edm::Handle<edm::View<pat::Muon> > muonColl;
  edm::View<pat::Muon>::const_iterator muon;

  iEvent.getByToken( muonToken_, muonColl );

  for ( unsigned int i = 0; i < muonColl->size() && evt_.nMuonCand < ggll::AnalysisEvent::MAX_MUONS; ++i ) {
    const edm::Ptr<pat::Muon> muon = muonColl->ptrAt( i);

    // JH
    if( muon->pt() < 20) continue;
    if( muon->isGlobalMuon() == 0 || muon->isTrackerMuon() == 0) continue;

    evt_.MuonCand_pt[evt_.nMuonCand] = muon->pt();
    evt_.MuonCand_eta[evt_.nMuonCand] = muon->eta();
    evt_.MuonCand_phi[evt_.nMuonCand] = muon->phi();
    evt_.MuonCand_e[evt_.nMuonCand] = muon->energy();
    evt_.MuonCand_charge[evt_.nMuonCand] = muon->charge();
    evt_.MuonCand_dxy[evt_.nMuonCand] = muon->dB();
    evt_.MuonCand_nstatseg[evt_.nMuonCand] = muon->numberOfMatchedStations();

    evt_.MuonCand_vtxx[evt_.nMuonCand] = muon->vertex().x();
    evt_.MuonCand_vtxy[evt_.nMuonCand] = muon->vertex().y();
    evt_.MuonCand_vtxz[evt_.nMuonCand] = muon->vertex().z();

    evt_.MuonCand_isglobal[evt_.nMuonCand] = muon->isGlobalMuon();
    evt_.MuonCand_istracker[evt_.nMuonCand] = muon->isTrackerMuon();
    evt_.MuonCand_isstandalone[evt_.nMuonCand] = muon->isStandAloneMuon();
    evt_.MuonCand_ispfmuon[evt_.nMuonCand] = muon->isPFMuon();

    TLorentzVector leptonptmp;
    leptonptmp.SetPtEtaPhiM( muon->pt(), muon->eta(), muon->phi(), muon->mass() );

    if ( evt_.MuonCand_istracker[evt_.nMuonCand] ) {
      evt_.MuonCand_npxlhits[evt_.nMuonCand] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();
      evt_.MuonCand_ntrklayers[evt_.nMuonCand] = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      leptonptmp.SetPtEtaPhiM( muon->innerTrack()->pt(), muon->innerTrack()->eta(), muon->innerTrack()->phi(), muon->mass() );
    }
    muonsMomenta_.insert( std::pair<int,TLorentzVector>( evt_.nMuonCand, leptonptmp ) );

    if ( evt_.MuonCand_isglobal[evt_.nMuonCand] && evt_.MuonCand_istracker[evt_.nMuonCand] ) {
      evt_.MuonCand_innerTrackPt[evt_.nMuonCand] = muon->innerTrack()->pt();
      evt_.MuonCand_innerTrackEta[evt_.nMuonCand] = muon->innerTrack()->eta();
      evt_.MuonCand_innerTrackPhi[evt_.nMuonCand] = muon->innerTrack()->phi();
      evt_.MuonCand_innerTrackVtxz[evt_.nMuonCand] = muon->innerTrack()->vertex().z();
      evt_.MuonCandTrack_nmuchits[evt_.nMuonCand] = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
      evt_.MuonCandTrack_chisq[evt_.nMuonCand] = muon->globalTrack()->normalizedChi2();
      const bool istight = ( evt_.MuonCand_ispfmuon[evt_.nMuonCand]
                        && ( evt_.MuonCandTrack_chisq[evt_.nMuonCand] < 10. )
                        && ( evt_.MuonCandTrack_nmuchits[evt_.nMuonCand] >= 1 )
                        && ( evt_.MuonCand_nstatseg[evt_.nMuonCand] >= 2 )
                        && ( evt_.MuonCand_dxy[evt_.nMuonCand] < 0.2 )
                        && ( evt_.MuonCand_npxlhits[evt_.nMuonCand] > 0 )
                        && ( evt_.MuonCand_ntrklayers[evt_.nMuonCand] > 5 ) );
      evt_.MuonCand_istight[evt_.nMuonCand] = istight;

      muonTransientTracks_[evt_.nMuonCand] = KalVtx_->build( *muon->innerTrack() );
    }

    evt_.nMuonCand++;
  }
  LogDebug( "GammaGammaLL" ) << "Passed Muon retrieval stage. Got " << evt_.nMuonCand << " muon(s)";
}

void
GammaGammaLL::fetchElectrons( const edm::Event& iEvent )
{
  // Get the electrons collection from the event
  edm::Handle<edm::View<pat::Electron> > eleColl;
  iEvent.getByToken( eleToken_, eleColl );

  /*edm::Handle<edm::ValueMap<bool> > medium_id_decisions, tight_id_decisions;
  iEvent.getByToken( eleMediumIdMapToken_, medium_id_decisions );
  iEvent.getByToken( eleTightIdMapToken_, tight_id_decisions );*/

  edm::Handle<double> rhoJECJets;
  iEvent.getByToken( fixedGridRhoFastjetAllToken_, rhoJECJets );//kt6PFJets

  for ( unsigned int i = 0; i < eleColl->size(); ++i ) {
    const edm::Ptr<pat::Electron> electron = eleColl->ptrAt( i );

    evt_.EleCand_et[evt_.nEleCand] = electron->et();
    evt_.EleCand_eta[evt_.nEleCand] = electron->eta();
    evt_.EleCand_phi[evt_.nEleCand] = electron->phi();
    evt_.EleCand_e[evt_.nEleCand] = electron->energy();
    evt_.EleCand_charge[evt_.nEleCand] = electron->charge();

    evt_.EleCand_vtxx[evt_.nEleCand] = electron->vertex().x();
    evt_.EleCand_vtxy[evt_.nEleCand] = electron->vertex().y();
    evt_.EleCand_vtxz[evt_.nEleCand] = electron->vertex().z();

    TLorentzVector leptonptmp;
    leptonptmp.SetPtEtaPhiM( electron->et(), electron->eta(), electron->phi(), electron->mass() );

    if ( electron->closestCtfTrackRef().isNonnull() ) { // Only for pat::Electron
      evt_.EleCand_innerTrackPt[evt_.nEleCand] = electron->closestCtfTrackRef()->pt();
      evt_.EleCand_innerTrackEta[evt_.nEleCand] = electron->closestCtfTrackRef()->eta();
      evt_.EleCand_innerTrackPhi[evt_.nEleCand] = electron->closestCtfTrackRef()->phi();
      evt_.EleCand_innerTrackVtxz[evt_.nEleCand] = electron->closestCtfTrackRef()->vertex().z();
      leptonptmp.SetPtEtaPhiM( evt_.EleCand_innerTrackPt[evt_.nEleCand], evt_.EleCand_innerTrackEta[evt_.nEleCand], evt_.EleCand_innerTrackPhi[evt_.nEleCand], electron->mass() );
      //eleTransientTrack[evt_.nEleCand] = KalVtx_->build( *electron->closestCtfTrackRef() );
    }
    eleTransientTracks_[evt_.nEleCand] = KalVtx_->build( *electron->gsfTrack() );

    electronsMomenta_.insert( std::pair<int,TLorentzVector>( evt_.nEleCand, leptonptmp ) );

    evt_.EleCand_deltaPhi[evt_.nEleCand] = electron->deltaPhiSuperClusterTrackAtVtx();
    evt_.EleCand_deltaEta[evt_.nEleCand] = electron->deltaEtaSuperClusterTrackAtVtx();
    evt_.EleCand_HoverE[evt_.nEleCand] = electron->hcalOverEcal();
    evt_.EleCand_trackiso[evt_.nEleCand] = electron->dr03TkSumPt() / electron->et();
    evt_.EleCand_ecaliso[evt_.nEleCand] = electron->dr03EcalRecHitSumEt() / electron->et();
    evt_.EleCand_hcaliso[evt_.nEleCand] = electron->dr03HcalTowerSumEt() / electron->et();
    evt_.EleCand_sigmaIetaIeta[evt_.nEleCand] = electron->sigmaIetaIeta();
    evt_.EleCand_convDist[evt_.nEleCand] = fabs( electron->convDist() );
    evt_.EleCand_convDcot[evt_.nEleCand] = fabs( electron->convDcot() );
    evt_.EleCand_ecalDriven[evt_.nEleCand] = electron->ecalDrivenSeed();

    const std::vector<pat::Electron::IdPair> ids = electron->electronIDs();
    for ( unsigned int j = 0; j < ids.size(); ++j ) {
      pat::Electron::IdPair idp = ids.at( j );
      //FIXME make me private attributes
      if ( eleMediumIdLabel_.find( idp.first ) != std::string::npos ) evt_.EleCand_mediumID[evt_.nEleCand] = idp.second;
      if ( eleTightIdLabel_.find( idp.first ) != std::string::npos ) evt_.EleCand_tightID[evt_.nEleCand] = idp.second;
    }

    //edm::RefToBase<pat::Electron> electronRef( eleColl->refAt( i ) );
    //evt_.EleCand_mediumID[evt_.nEleCand] = medium_id_decisions->operator[]( electronRef ),
    //evt_.EleCand_tightID[evt_.nEleCand] = tight_id_decisions->operator[]( electronRef ),

    evt_.nEleCand++;
  }
  LogDebug( "GammaGammaLL" ) << "Passed Electron retrieval stage. Got " << evt_.nEleCand << " electron(s)";
}

void
GammaGammaLL::fetchProtons( const edm::Event& iEvent )
{
  // Forward proton tracks
  edm::Handle<edm::View<CTPPSLocalTrackLite> > ppslocaltracks;
  iEvent.getByToken( ppsLocalTrackToken_, ppslocaltracks );

  evt_.nLocalProtCand = 0;
  for ( const auto trk : *ppslocaltracks ) {
    const CTPPSDetId det_id( trk.getRPId() );
    if ( evt_.nLocalProtCand == ggll::AnalysisEvent::MAX_LOCALPCAND-1 )
      throw cms::Exception( "GammaGammaLL" ) << "maximum number of local tracks in RPs is reached!\n"
        << "increase MAX_LOCALPCAND=" << ggll::AnalysisEvent::MAX_LOCALPCAND << " in GammaGammaLL.cc";

    // transform the raw, 32-bit unsigned integer detId into the TOTEM "decimal" notation                                          
    const unsigned int raw_id = 100*det_id.arm()+10*det_id.station()+det_id.rp();

    evt_.LocalProtCand_x[evt_.nLocalProtCand] = ( trk.getX() )/1.e3;
    evt_.LocalProtCand_y[evt_.nLocalProtCand] = ( trk.getY() )/1.e3;
    evt_.LocalProtCand_xSigma[evt_.nLocalProtCand] = ( trk.getXUnc() )/1.e3;
    evt_.LocalProtCand_ySigma[evt_.nLocalProtCand] = ( trk.getYUnc() )/1.e3;
    evt_.LocalProtCand_arm[evt_.nLocalProtCand] = det_id.arm();
    evt_.LocalProtCand_station[evt_.nLocalProtCand] = det_id.station();
    evt_.LocalProtCand_pot[evt_.nLocalProtCand] = det_id.rp();
    evt_.LocalProtCand_rpid[evt_.nLocalProtCand] = raw_id;
    evt_.nLocalProtCand++;
    LogDebug( "GammaGammaLL" ) << "Proton track candidate with origin: ( " << trk.getX() << ", " << trk.getY() << " ) extracted!";
  }
  LogDebug( "GammaGammaLL" ) << "Passed TOTEM RP info retrieval stage. Got " << evt_.nLocalProtCand << " local track(s)";

  // Full reco protons - 2 separate collections for legacy re-reco
  edm::Handle<vector<reco::ForwardProton>> recoMultiRPProtons;
  iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
  edm::Handle<vector<reco::ForwardProton>> recoSingleRPProtons;
  iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);

  int ismultirp = -999;
  unsigned int decRPId = -999;
  unsigned int armId = -999;
  float th_y = -999;
  float th_x = -999;
  float t = -999;
  float xi = -999;

  evt_.nRecoProtCand = 0;

  // Single-RP algorithm
  for (const auto & proton : *recoSingleRPProtons)
    {
      if (proton.validFit())
        {
          th_y = proton.thetaY();
	  th_x = proton.thetaX();
          xi = proton.xi();
          t = proton.t();

	  CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->getRPId());
	  decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
          ismultirp = 0;

          evt_.ProtCand_xi[evt_.nRecoProtCand] = xi;
          evt_.ProtCand_t[evt_.nRecoProtCand] = t;
          evt_.ProtCand_ThX[evt_.nRecoProtCand] = th_x;
          evt_.ProtCand_ThY[evt_.nRecoProtCand] = th_y;
          evt_.ProtCand_rpid[evt_.nRecoProtCand] = decRPId;
          evt_.ProtCand_arm[evt_.nRecoProtCand] = armId;
          evt_.ProtCand_ismultirp[evt_.nRecoProtCand] = ismultirp;
          evt_.nRecoProtCand++;
	}
    }

  // Multi-RP algorithm
  for (const auto & proton : *recoMultiRPProtons)
    {
      if (proton.validFit())
	{
	  th_y = proton.thetaY();
	  th_x = proton.thetaX();
	  xi = proton.xi();
	  t = proton.t();

	  CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->getRPId());
	  armId = rpId.arm();
	  ismultirp = 1;

	  evt_.ProtCand_xi[evt_.nRecoProtCand] = xi;
	  evt_.ProtCand_t[evt_.nRecoProtCand] = t;
          evt_.ProtCand_ThX[evt_.nRecoProtCand] = th_x;
          evt_.ProtCand_ThY[evt_.nRecoProtCand] = th_y;
          evt_.ProtCand_rpid[evt_.nRecoProtCand] = decRPId;
	  evt_.ProtCand_arm[evt_.nRecoProtCand] = armId;
	  evt_.ProtCand_ismultirp[evt_.nRecoProtCand] = ismultirp;
	  evt_.nRecoProtCand++;
	}
    }

}

void
GammaGammaLL::fetchVertices( const edm::Event& iEvent )
{
  // Get the vertex collection from the event
  edm::Handle<edm::View<reco::Vertex> > recoVertexColl;
  iEvent.getByToken( recoVertexToken_, recoVertexColl );

  for ( unsigned int i = 0; i < recoVertexColl->size() && evt_.nPrimVertexCand < ggll::AnalysisEvent::MAX_VTX; ++i ) {
    const edm::Ptr<reco::Vertex> vertex = recoVertexColl->ptrAt( i);

    evt_.PrimVertexCand_x[evt_.nPrimVertexCand] = vertex->x();
    evt_.PrimVertexCand_y[evt_.nPrimVertexCand] = vertex->y();
    evt_.PrimVertexCand_z[evt_.nPrimVertexCand] = vertex->z();
    evt_.PrimVertexCand_tracks[evt_.nPrimVertexCand] = vertex->nTracks();
    evt_.PrimVertexCand_chi2[evt_.nPrimVertexCand] = vertex->chi2();
    evt_.PrimVertexCand_ndof[evt_.nPrimVertexCand] = vertex->ndof();
    evt_.nPrimVertexCand++;
  }
}

void
GammaGammaLL::newVertexInfoRetrieval( const edm::Event& iEvent )
{
  iEvent.getByToken( recoTrackToken_, trackColl_ );

  foundPairInEvent_ = false;
  switch ( leptonsType_ ) {
    case ggll::ElectronMuon: {
      for ( unsigned int i = 0; i < evt_.nMuonCand; ++i ) {
        for ( unsigned int j = 0; j < evt_.nEleCand; ++j ) {
          if ( evt_.MuonCand_charge[i]*evt_.EleCand_charge[j] > 0 ) continue;
          if ( newTracksInfoRetrieval( i, j ) ) foundPairInEvent_ = true;
        }
      }
    } break;
    case ggll::DiElectron: {
      for ( unsigned int i = 0; i < evt_.nEleCand; ++i ) {
        for ( unsigned int j = i+1; j < evt_.nEleCand; ++j ) {
          if ( evt_.EleCand_charge[i]*evt_.EleCand_charge[j] > 0 ) continue;
          if ( newTracksInfoRetrieval( i, j ) ) foundPairInEvent_ = true;
        }
      }
    } break;
    case ggll::DiMuon: {
      for ( unsigned int i = 0; i < evt_.nMuonCand; ++i ) {
       for ( unsigned int j = i+1; j < evt_.nMuonCand; ++j ) {
          if ( evt_.MuonCand_charge[i]*evt_.MuonCand_charge[j] > 0 ) continue;
          if ( newTracksInfoRetrieval( i, j ) ) foundPairInEvent_ = true;
        }
      }
    } break;
    default:
      throw cms::Exception( "GammaGammaLL" ) << "invalid leptons type to be retrieved: " << leptonsType_;
  }
  fetchVertices( iEvent );
}

bool
GammaGammaLL::newTracksInfoRetrieval( int l1id, int l2id )
{
  double l1cand_pt, l2cand_pt;
  std::vector<reco::TransientTrack> translepttrks;

  switch ( leptonsType_ ) {
    case ggll::ElectronMuon: {
      translepttrks.push_back(muonTransientTracks_[l1id] );
      translepttrks.push_back(eleTransientTracks_[l2id] );
      l1cand_pt = evt_.MuonCand_innerTrackPt[l1id];
      l2cand_pt = evt_.EleCand_innerTrackPt[l2id];
    } break;
    case ggll::DiElectron: {
      translepttrks.push_back(eleTransientTracks_[l1id] );
      translepttrks.push_back(eleTransientTracks_[l2id] );
      l1cand_pt = evt_.EleCand_innerTrackPt[l1id];
      l2cand_pt = evt_.EleCand_innerTrackPt[l2id];
    } break;
    case ggll::DiMuon: {
      translepttrks.push_back(muonTransientTracks_[l1id] );
      translepttrks.push_back(muonTransientTracks_[l2id] );
      l1cand_pt = evt_.MuonCand_innerTrackPt[l1id];
      l2cand_pt = evt_.MuonCand_innerTrackPt[l2id];
    } break;
    default: throw cms::Exception( "GammaGammaLL" ) << "Invalid leptons type: " << leptonsType_;
  }

  if ( translepttrks.size() < 2 ) return false; // just in case...

  const bool use_smoothing = true;
  KalmanVertexFitter fitter( use_smoothing );
  TransientVertex dileptonVertex;
  try {
    dileptonVertex = fitter.vertex(translepttrks);
  } catch (...) { return false; }

  if ( !dileptonVertex.isValid() ) return false; // only keep the pairs with valid vertex

  evt_.KalmanVertexCand_x[evt_.nPair] = dileptonVertex.position().x();
  evt_.KalmanVertexCand_y[evt_.nPair] = dileptonVertex.position().y();
  evt_.KalmanVertexCand_z[evt_.nPair] = dileptonVertex.position().z();

  // Count nearby tracks
  int closesttrkid = -1, closesthighpuritytrkid = -1;
  double closesttrkdxyz = 999., closesthighpuritytrkdxyz = 999.;

  for ( unsigned int i = 0; i < trackColl_->size() && evt_.nExtraTracks < ggll::AnalysisEvent::MAX_ET; ++i ) {
    const edm::Ptr<reco::Track> track = trackColl_->ptrAt( i);

    if ( track->pt() == l1cand_pt || track->pt() == l2cand_pt) continue;

    const double vtxdst = sqrt( std::pow( ( track->vertex().x()-evt_.KalmanVertexCand_x[evt_.nPair] ), 2 )
                              + std::pow( ( track->vertex().y()-evt_.KalmanVertexCand_y[evt_.nPair] ), 2 )
                              + std::pow( ( track->vertex().z()-evt_.KalmanVertexCand_z[evt_.nPair] ), 2 ) );
    if ( vtxdst < 0.05 ) evt_.Pair_extratracks0p5mm[evt_.nPair]++;
    if ( vtxdst < 0.1 ) evt_.Pair_extratracks1mm[evt_.nPair]++;
    if ( vtxdst < 0.2 ) evt_.Pair_extratracks2mm[evt_.nPair]++;
    if ( vtxdst < 0.3 ) evt_.Pair_extratracks3mm[evt_.nPair]++;
    if ( vtxdst < 0.4 ) evt_.Pair_extratracks4mm[evt_.nPair]++;
    if ( vtxdst < 0.5 ) evt_.Pair_extratracks5mm[evt_.nPair]++;
    if ( vtxdst < 1.0 ) evt_.Pair_extratracks1cm[evt_.nPair]++;
    if ( vtxdst < 2.0 ) evt_.Pair_extratracks2cm[evt_.nPair]++;
    if ( vtxdst < 3.0 ) evt_.Pair_extratracks3cm[evt_.nPair]++;
    if ( vtxdst < 4.0 ) evt_.Pair_extratracks4cm[evt_.nPair]++;
    if ( vtxdst < 5.0 ) evt_.Pair_extratracks5cm[evt_.nPair]++;
    if ( vtxdst < 10. ) evt_.Pair_extratracks10cm[evt_.nPair]++;
    if ( vtxdst < closesttrkdxyz ) {
      closesttrkid = i;
      closesttrkdxyz = vtxdst;
    }
    if ( track->quality( reco::TrackBase::highPurity ) == 1 && vtxdst < closesthighpuritytrkdxyz ) {
      closesthighpuritytrkid = i;
      closesthighpuritytrkdxyz = vtxdst;
    }

    // Save track properties if within 5mm
    if ( (saveExtraTracks_ == true) && (vtxdst < 0.5) ) {
      evt_.ExtraTrack_pair[evt_.nExtraTracks] = evt_.nPair;
      evt_.ExtraTrack_purity[evt_.nExtraTracks] = track->quality( reco::TrackBase::highPurity );
      evt_.ExtraTrack_nhits[evt_.nExtraTracks] = track->numberOfValidHits();

      evt_.ExtraTrack_px[evt_.nExtraTracks] = track->px();
      evt_.ExtraTrack_py[evt_.nExtraTracks] = track->py();
      evt_.ExtraTrack_pz[evt_.nExtraTracks] = track->pz();
      evt_.ExtraTrack_charge[evt_.nExtraTracks] = track->charge();
      evt_.ExtraTrack_chi2[evt_.nExtraTracks] = track->chi2();
      evt_.ExtraTrack_ndof[evt_.nExtraTracks] = track->ndof();
      evt_.ExtraTrack_vtxdxyz[evt_.nExtraTracks] = vtxdst;
      evt_.ExtraTrack_vtxT[evt_.nExtraTracks] = std::hypot( track->vertex().x()-evt_.KalmanVertexCand_x[evt_.nPair],
                                                            track->vertex().y()-evt_.KalmanVertexCand_y[evt_.nPair] );
      evt_.ExtraTrack_vtxZ[evt_.nExtraTracks] = fabs( track->vertex().z()-evt_.KalmanVertexCand_z[evt_.nPair] );
      evt_.ExtraTrack_x[evt_.nExtraTracks] = track->vertex().x();
      evt_.ExtraTrack_y[evt_.nExtraTracks] = track->vertex().y();
      evt_.ExtraTrack_z[evt_.nExtraTracks] = track->vertex().z();

      evt_.nExtraTracks++;
    }
  }
  evt_.ClosestExtraTrack_vtxdxyz[evt_.nPair] = closesttrkdxyz;
  evt_.ClosestExtraTrack_id[evt_.nPair] = closesttrkid;
  evt_.ClosestHighPurityExtraTrack_vtxdxyz[evt_.nPair] = closesthighpuritytrkdxyz;
  evt_.ClosestHighPurityExtraTrack_id[evt_.nPair] = closesthighpuritytrkid;

  TLorentzVector l1, l2;
  switch (leptonsType_ ) {
    case ggll::ElectronMuon: {
      l1.SetPtEtaPhiE( evt_.MuonCand_pt[l1id], evt_.MuonCand_eta[l1id], evt_.MuonCand_phi[l1id], evt_.MuonCand_e[l1id] );
      l2.SetPtEtaPhiE( evt_.EleCand_et[l2id], evt_.EleCand_eta[l2id], evt_.EleCand_phi[l2id], evt_.EleCand_e[l2id] );
    } break;
    case ggll::DiElectron: {
      l1.SetPtEtaPhiE( evt_.EleCand_et[l1id], evt_.EleCand_eta[l1id], evt_.EleCand_phi[l1id], evt_.EleCand_e[l1id] );
      l2.SetPtEtaPhiE( evt_.EleCand_et[l2id], evt_.EleCand_eta[l2id], evt_.EleCand_phi[l2id], evt_.EleCand_e[l2id] );
    } break;
    case ggll::DiMuon: {
      l1.SetPtEtaPhiE( evt_.MuonCand_pt[l1id], evt_.MuonCand_eta[l1id], evt_.MuonCand_phi[l1id], evt_.MuonCand_e[l1id] );
      l2.SetPtEtaPhiE( evt_.MuonCand_pt[l2id], evt_.MuonCand_eta[l2id], evt_.MuonCand_phi[l2id], evt_.MuonCand_e[l2id] );
    } break;
    default: throw cms::Exception( "GammaGammaLL" ) << "Invalid leptons type: " << leptonsType_;
  }

  const TLorentzVector pair = l1+l2;

  if ( minMpair_ > 0. && pair.M() < minMpair_ )
    return false;
  if ( maxMpair_ > 0. && pair.M() > maxMpair_ )
    return false;

  evt_.Pair_pt[evt_.nPair] = pair.Pt();
  evt_.Pair_mass[evt_.nPair] = pair.M();
  evt_.Pair_phi[evt_.nPair] = pair.Phi();
  evt_.Pair_eta[evt_.nPair] = pair.Eta();
  evt_.Pair_lepton1[evt_.nPair] = l1id;
  evt_.Pair_lepton2[evt_.nPair] = l2id;

  double dphi = fabs( l1.Phi()-l2.Phi() );
  // dphi lies in [-pi, pi]
  while ( dphi < -M_PI ) dphi += 2.*M_PI;
  while ( dphi >  M_PI ) dphi -= 2.*M_PI;
  evt_.Pair_dphi[evt_.nPair] = dphi;

  evt_.Pair_dpt[evt_.nPair] = fabs( l1.Pt()-l2.Pt() );
  evt_.Pair_3Dangle[evt_.nPair] = l1.Angle( l2.Vect() )/M_PI;

  evt_.nPair++;
  nCandidates_++;


  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void
GammaGammaLL::beginJob()
{
  // Filling the ntuple
  evt_.attach( tree_, leptonsType_, runOnMC_, saveExtraTracks_ );
  *evt_.HLT_Name = triggersList_;
}

// ------------ method called once each job just after ending the event loop  ------------
void
GammaGammaLL::endJob()
{
  std::cout << "==> Number of candidates in the dataset : " << nCandidates_ << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void
GammaGammaLL::beginRun( const edm::Run& iRun, const edm::EventSetup& iSetup )
{
  bool changed = true;
  if ( !hltPrescale_.init( iRun, iSetup, triggerResults_.process(), changed ) )
    edm::LogError( "GammaGammaLL" ) << "prescales extraction failure with process name " << triggerResults_.process();

  // Initialise HLTConfigProvider
  hltConfig_ = hltPrescale_.hltConfigProvider();
  if ( hltConfig_.size() == 0 )
    edm::LogError( "GammaGammaLL" ) << "HLT config size error";
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GammaGammaLL::fillDescriptions( edm::ConfigurationDescriptions& descriptions ) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( GammaGammaLL );

