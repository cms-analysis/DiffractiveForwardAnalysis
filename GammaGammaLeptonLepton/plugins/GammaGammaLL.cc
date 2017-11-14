// -*- C++ -*-
//
// Package:    GammaGammaLL
// Class:      GammaGammaLL
//
/**\class GammaGammaLL GammaGammaLL.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/src/GammaGammaLL.cc

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

#include "GammaGammaLL.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/PrimaryVertexSelector.h"

GammaGammaLL::GammaGammaLL( const edm::ParameterSet& iConfig ) :
  tree_( 0 ),
  fetchMuons_( false ), fetchElectrons_( false ),
  fetchProtons_       ( iConfig.getParameter<bool>( "fetchProtons" ) ),
  hltMenuLabel_       ( iConfig.getParameter<std::string>( "HLTMenuTag" ) ),
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
  jetToken_           ( consumes<edm::View<pat::Jet> >                 ( iConfig.getParameter<edm::InputTag>( "jetTag" ) ) ),
  fixedGridRhoFastjetAllToken_( consumes<double>                       ( iConfig.getParameter<edm::InputTag>( "fixedGridRhoFastjetAllLabel" ) ) ),
  metToken_           ( consumes<edm::View<pat::MET> >                 ( iConfig.getParameter<edm::InputTag>( "metTag" ) ) ),
  totemRPHitToken_    ( consumes<edm::DetSetVector<TotemRPLocalTrack> >( iConfig.getParameter<edm::InputTag>( "totemRPLocalTrackTag" ) ) ),
  photonToken_        ( consumes<edm::View<pat::Photon> >              ( iConfig.getParameter<edm::InputTag>( "photonTag" ) ) ),
  runOnMC_            ( iConfig.getParameter<bool>( "runOnMC" ) ),
  printCandidates_    ( iConfig.getParameter<bool>( "printCandidates" ) ),
  sqrts_              ( iConfig.getParameter<double>( "sqrtS" ) ),
  maxExTrkVtx_        ( iConfig.getUntrackedParameter<unsigned int>( "maxExtraTracks", 1000 ) ),
  hlts_               ( triggersList_ ),
  hltPrescale_        ( iConfig, consumesCollector(), *this ),
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
  // photon identification variables
  const edm::ParameterSet phoIdLabelSet = iConfig.getParameter<edm::ParameterSet>( "phoIdLabels" );
  phoMediumIdLabel_ = phoIdLabelSet.getParameter<edm::InputTag>( "mediumLabel" ).encode();
  phoTightIdLabel_ = phoIdLabelSet.getParameter<edm::InputTag>( "tightLabel" ).encode();

  // Pileup reweighting utilities
  if ( runOnMC_ ) {
    lumiWeights_ = edm::LumiReWeighting(mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_ );
  }

  // Book the output tree
  usesResource( "TFileService" );
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>( "ntp1", "ntp1" );
}

GammaGammaLL::~GammaGammaLL()
{}

//
// member functions
//

void
GammaGammaLL::lookAtTriggers( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // Get the trigger information from the event
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken( triggerResultsToken_, hltResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);

  std::ostringstream os;
  os << "Trigger names: " << std::endl;

  evt_.HLT_Accept.reserve( trigNames.size() );
  evt_.HLT_Prescl.reserve( trigNames.size() );

  for ( unsigned int i = 0; i < trigNames.size(); ++i ) {
    os << "--> " << trigNames.triggerNames().at( i ) << std::endl;

    const int trigNum = hlts_.TriggerNum( trigNames.triggerNames().at( i ) );
    if ( trigNum < 0 ) continue; // Trigger didn't match the interesting ones

    evt_.HLT_Accept[trigNum] = hltResults->accept( i );

    // extract prescale value for this path
    if ( runOnMC_ ) { evt_.HLT_Prescl[trigNum] = 1.; continue; } //FIXME
    int prescale_set = hltPrescale_.prescaleSet( iEvent, iSetup);
    evt_.HLT_Prescl[trigNum] = ( prescale_set < 0 ) ? 0. : hltConfig_.prescaleValue( prescale_set, trigNames.triggerNames().at( i ) ); //FIXME
  }
  LogDebug( "GammaGammaLL" ) << os.str();
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
  lookAtTriggers( iEvent, iSetup);

  LogDebug( "GammaGammaLL" ) << "Passed trigger filtering stage";

  // Generator level information
  if ( runOnMC_ ) {
    analyzeMCEventContent( iEvent );

    // Pileup information
    edm::Handle<edm::View<PileupSummaryInfo> > pu_info;
    iEvent.getByToken( pileupToken_, pu_info);

    int npv0true = 0;
    for ( unsigned int i = 0; i < pu_info->size(); ++i ) {
      const edm::Ptr<PileupSummaryInfo> PVI = pu_info->ptrAt( i );

      const int beamXing = PVI->getBunchCrossing(),
                npvtrue = PVI->getTrueNumInteractions();

      if ( beamXing == 0 ) npv0true += npvtrue;
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

  fetchPhotons( iEvent );

  newVertexInfoRetrieval( iEvent );

  if ( !foundPairInEvent_ ) {
    //LogDebug( "GammaGammaLL" ) << "No pair retrieved in event";
    return; // avoid to unpack RP/jet/MET if no dilepton candidate has been found
  }

  if ( fetchProtons_ ) fetchProtons( iEvent );

  fetchJets( iEvent );

  // Missing ET
  edm::Handle<edm::View<pat::MET> > MET;
  iEvent.getByToken( metToken_, MET);
  const edm::View<pat::MET>* metColl = MET.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();

  evt_.Etmiss = met->et();
  evt_.Etmiss_phi = met->phi();
  evt_.Etmiss_significance = met->significance();

  LogDebug( "GammaGammaLL" ) << "Passed MET retrieval stage";

  if ( printCandidates_ ) {
    std::cout << "Event " << evt_.Run << ":" << evt_.EventNum << " has " << evt_.nPair << " leptons pair(s) candidate(s) (vertex mult. : " << evt_.nPrimVertexCand << " )" << std::endl;
  }

  tree_->Fill();
}

void
GammaGammaLL::analyzeMCEventContent( const edm::Event& iEvent )
{
  edm::Handle<edm::View<reco::GenParticle> > genPartColl;
  iEvent.getByToken( genToken_, genPartColl );

  unsigned int ngpro = 0, ngmu = 0, ngele = 0, ngpho = 0;
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
    if ( genPart->pdgId() == 2212 ) {
      evt_.GenProtCand_pt[ngpro] = genPart->pt();
      evt_.GenProtCand_eta[ngpro] = genPart->eta();
      evt_.GenProtCand_phi[ngpro] = genPart->phi();
      evt_.GenProtCand_e[ngpro] = genPart->energy();
      evt_.GenProtCand_status[ngpro] = genPart->status();
      ngpro++;
    }

    // generated central dimuon system
    if ( fabs( genPart->pdgId() ) == 13 ) {
      evt_.GenMuonCand_pt[ngmu] = genPart->pt();
      evt_.GenMuonCand_eta[ngmu] = genPart->eta();
      evt_.GenMuonCand_phi[ngmu] = genPart->phi();
      evt_.GenMuonCand_e[ngmu] = genPart->energy();
      ngmu++;
    }

    // generated central dimuon system
    if ( fabs( genPart->pdgId() ) == 11 ) {
      evt_.GenEleCand_pt[ngele] = genPart->pt();
      evt_.GenEleCand_eta[ngele] = genPart->eta();
      evt_.GenEleCand_phi[ngele] = genPart->phi();
      evt_.GenEleCand_e[ngele] = genPart->energy();
      ngele++;
    }

    // generated inner photon line
    if ( genPart->pdgId() == 22 ) {
      evt_.GenPhotCand_e[ngpho] = genPart->energy();
      evt_.GenPhotCand_pt[ngpho] = genPart->pt();
      evt_.GenPhotCand_eta[ngpho] = genPart->eta();
      evt_.GenPhotCand_phi[ngpho] = genPart->phi();
      ngpho++;
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
  if ( leptonsType_ == ggll::ElectronMuon && evt_.GenMuonCand_pt.size() == 1 && evt_.GenEleCand_pt.size() == 1 ) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE( evt_.GenMuonCand_pt[0], evt_.GenMuonCand_eta[0], evt_.GenMuonCand_phi[0], evt_.GenMuonCand_e[0] );
    l2.SetPtEtaPhiE( evt_.GenEleCand_pt[0], evt_.GenEleCand_eta[0], evt_.GenEleCand_phi[0], evt_.GenEleCand_e[0] );
    foundGenCandPairInEvent = true;
  }
  // dielectron
  else if ( leptonsType_ == ggll::DiElectron && evt_.GenEleCand_pt.size() == 2 ) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE( evt_.GenEleCand_pt[0], evt_.GenEleCand_eta[0], evt_.GenEleCand_phi[0], evt_.GenEleCand_e[0] );
    l2.SetPtEtaPhiE( evt_.GenEleCand_pt[1], evt_.GenEleCand_eta[1], evt_.GenEleCand_phi[1], evt_.GenEleCand_e[1] );
    foundGenCandPairInEvent = true;
  }
  // dimuon
  else if ( leptonsType_ == ggll::DiMuon && evt_.GenMuonCand_pt.size() == 2 ) { // FIXME maybe a bit tight according to the newer PU conditions?
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
    while ( dphi < -M_PI ) { dphi += 2.*M_PI; }
    while ( dphi >  M_PI ) { dphi -= 2.*M_PI; }
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

  unsigned int nmu = 0;
  for ( unsigned int i = 0; i < muonColl->size(); ++i ) {
    const edm::Ptr<pat::Muon> muon = muonColl->ptrAt( i );

    evt_.MuonCand_pt[nmu] = muon->pt();
    evt_.MuonCand_eta[nmu] = muon->eta();
    evt_.MuonCand_phi[nmu] = muon->phi();
    evt_.MuonCand_e[nmu] = muon->energy();
    evt_.MuonCand_charge[nmu] = muon->charge();
    evt_.MuonCand_dxy[nmu] = muon->dB();
    evt_.MuonCand_nstatseg[nmu] = muon->numberOfMatchedStations();

    evt_.MuonCand_vtxx[nmu] = muon->vertex().x();
    evt_.MuonCand_vtxy[nmu] = muon->vertex().y();
    evt_.MuonCand_vtxz[nmu] = muon->vertex().z();

    evt_.MuonCand_isglobal[nmu] = muon->isGlobalMuon();
    evt_.MuonCand_istracker[nmu] = muon->isTrackerMuon();
    evt_.MuonCand_isstandalone[nmu] = muon->isStandAloneMuon();
    evt_.MuonCand_ispfmuon[nmu] = muon->isPFMuon();

    if ( evt_.MuonCand_istracker[nmu] ) {
      evt_.MuonCand_npxlhits[nmu] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();
      evt_.MuonCand_ntrklayers[nmu] = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      muonsMomenta_[nmu].SetPtEtaPhiM( muon->innerTrack()->pt(), muon->innerTrack()->eta(), muon->innerTrack()->phi(), muon->mass() );
    }
    else {
      evt_.MuonCand_npxlhits[nmu] = evt_.MuonCand_ntrklayers[nmu] = -1;
      muonsMomenta_[nmu].SetPtEtaPhiM( muon->pt(), muon->eta(), muon->phi(), muon->mass() );
    }

    if ( evt_.MuonCand_isglobal[nmu] && evt_.MuonCand_istracker[nmu] ) {
      evt_.MuonCand_innerTrackPt[nmu] = muon->innerTrack()->pt();
      evt_.MuonCand_innerTrackEta[nmu] = muon->innerTrack()->eta();
      evt_.MuonCand_innerTrackPhi[nmu] = muon->innerTrack()->phi();
      evt_.MuonCand_innerTrackVtxz[nmu] = muon->innerTrack()->vertex().z();
      evt_.MuonCandTrack_nmuchits[nmu] = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
      evt_.MuonCandTrack_chisq[nmu] = muon->globalTrack()->normalizedChi2();
      const bool istight = ( evt_.MuonCand_ispfmuon[nmu]
                        && ( evt_.MuonCandTrack_chisq[nmu] < 10. )
                        && ( evt_.MuonCandTrack_nmuchits[nmu] >= 1 )
                        && ( evt_.MuonCand_nstatseg[nmu] >= 2 )
                        && ( evt_.MuonCand_dxy[nmu] < 0.2 )
                        && ( evt_.MuonCand_npxlhits[nmu] > 0 )
                        && ( evt_.MuonCand_ntrklayers[nmu] > 5 ) );
      evt_.MuonCand_istight[nmu] = istight;

      muonTransientTracks_[nmu] = KalVtx_->build( *muon->innerTrack() );
    }
    else {
      evt_.MuonCand_innerTrackPt[nmu] = evt_.MuonCand_innerTrackEta[nmu] = evt_.MuonCand_innerTrackPhi[nmu] = evt_.MuonCand_innerTrackVtxz[nmu] = 0.;
      evt_.MuonCandTrack_nmuchits[nmu] = evt_.MuonCand_istight[nmu] = 0;
      evt_.MuonCandTrack_chisq[nmu] = 0.;
      evt_.MuonCand_istight[nmu] = 0;
    }

    nmu++;
  }
  LogDebug( "GammaGammaLL" ) << "Passed Muon retrieval stage. Got " << nmu << " muon(s)";
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

  unsigned int nele = 0;
  for ( unsigned int i = 0; i < eleColl->size(); ++i ) {
    const edm::Ptr<pat::Electron> electron = eleColl->ptrAt( i );

    evt_.EleCand_et[nele] = electron->et();
    evt_.EleCand_eta[nele] = electron->eta();
    evt_.EleCand_phi[nele] = electron->phi();
    evt_.EleCand_e[nele] = electron->energy();
    evt_.EleCand_charge[nele] = electron->charge();

    evt_.EleCand_vtxx[nele] = electron->vertex().x();
    evt_.EleCand_vtxy[nele] = electron->vertex().y();
    evt_.EleCand_vtxz[nele] = electron->vertex().z();

    if ( electron->closestCtfTrackRef().isNonnull() ) { // Only for pat::Electron
      evt_.EleCand_innerTrackPt[nele] = electron->closestCtfTrackRef()->pt();
      evt_.EleCand_innerTrackEta[nele] = electron->closestCtfTrackRef()->eta();
      evt_.EleCand_innerTrackPhi[nele] = electron->closestCtfTrackRef()->phi();
      evt_.EleCand_innerTrackVtxz[nele] = electron->closestCtfTrackRef()->vertex().z();
      electronsMomenta_[nele].SetPtEtaPhiM( evt_.EleCand_innerTrackPt[nele], evt_.EleCand_innerTrackEta[nele], evt_.EleCand_innerTrackPhi[nele], electron->mass() );
      //eleTransientTrack[nele] = KalVtx_->build( *electron->closestCtfTrackRef() );
    }
    else {
      evt_.EleCand_innerTrackPt[nele] = evt_.EleCand_innerTrackEta[nele] = evt_.EleCand_innerTrackPhi[nele] = evt_.EleCand_innerTrackVtxz[nele] = 0.;
      electronsMomenta_[nele].SetPtEtaPhiM( electron->et(), electron->eta(), electron->phi(), electron->mass() );
    }

    eleTransientTracks_[nele] = KalVtx_->build( *electron->gsfTrack() );

    evt_.EleCand_deltaPhi[nele] = electron->deltaPhiSuperClusterTrackAtVtx();
    evt_.EleCand_deltaEta[nele] = electron->deltaEtaSuperClusterTrackAtVtx();
    evt_.EleCand_HoverE[nele] = electron->hcalOverEcal();
    evt_.EleCand_trackiso[nele] = electron->dr03TkSumPt() / electron->et();
    evt_.EleCand_ecaliso[nele] = electron->dr03EcalRecHitSumEt() / electron->et();
    evt_.EleCand_hcaliso[nele] = electron->dr03HcalTowerSumEt() / electron->et();
    evt_.EleCand_sigmaIetaIeta[nele] = electron->sigmaIetaIeta();
    evt_.EleCand_convDist[nele] = fabs( electron->convDist() );
    evt_.EleCand_convDcot[nele] = fabs( electron->convDcot() );
    evt_.EleCand_ecalDriven[nele] = electron->ecalDrivenSeed();

    evt_.EleCand_mediumID[nele] = evt_.EleCand_tightID[nele] = 0;
    const std::vector<pat::Electron::IdPair> ids = electron->electronIDs();
    for ( unsigned int j = 0; j < ids.size(); ++j ) {
      pat::Electron::IdPair idp = ids.at( j );
      //FIXME make me private attributes
      if ( eleMediumIdLabel_.find( idp.first ) != std::string::npos ) evt_.EleCand_mediumID[nele] = idp.second;
      if ( eleTightIdLabel_.find( idp.first ) != std::string::npos ) evt_.EleCand_tightID[nele] = idp.second;
    }

    //edm::RefToBase<pat::Electron> electronRef( eleColl->refAt( i ) );
    //evt_.EleCand_mediumID[nele] = medium_id_decisions->operator[]( electronRef ),
    //evt_.EleCand_tightID[nele] = tight_id_decisions->operator[]( electronRef ),

    nele++;
  }
  LogDebug( "GammaGammaLL" ) << "Passed Electron retrieval stage. Got " << nele << " electron(s)";
}

void
GammaGammaLL::fetchPhotons( const edm::Event& iEvent )
{
  // Get the photons collection from the event
  edm::Handle<edm::View<pat::Photon> > photonColl;
  iEvent.getByToken( photonToken_, photonColl );

  // identification
  /*edm::Handle<edm::ValueMap<bool> > medium_id_decisions, tight_id_decisions;
  iEvent.getByToken( phoMediumIdMapToken_, medium_id_decisions );
  iEvent.getByToken( phoTightIdMapToken_, tight_id_decisions );*/

  unsigned int np = 0;
  for ( unsigned int i = 0; i < photonColl->size(); ++i ) {
    const edm::Ptr<pat::Photon> photon = photonColl->ptrAt( i );

    evt_.PhotonCand_pt[np] = photon->pt();
    evt_.PhotonCand_eta[np] = photon->eta();
    evt_.PhotonCand_phi[np] = photon->phi();
    evt_.PhotonCand_e[np] = photon->energy();
    evt_.PhotonCand_r9[np] = photon->r9();

    if ( runOnMC_ ) {
      double photdr = 999., photdeta = 999., photdphi = 999.;
      double endphotdr = 999., endphotdeta = 999., endphotdphi = 999.;
      for ( unsigned int j = 0; j < evt_.GenPhotCand_pt.size(); ++j ) { // matching with the 'true' photon object from MC
        photdeta = ( evt_.PhotonCand_eta[np]-evt_.GenPhotCand_eta[j] );
        photdphi = ( evt_.PhotonCand_phi[np]-evt_.GenPhotCand_phi[j] );
        photdr = sqrt( photdeta*photdeta + photdphi*photdphi );
        if ( photdr < endphotdr ) {
          endphotdr = photdr;
          endphotdeta = photdeta;
          endphotdphi = photdphi;
        }
      }
      evt_.PhotonCand_detatrue[np] = endphotdeta;
      evt_.PhotonCand_dphitrue[np] = endphotdphi;
      evt_.PhotonCand_drtrue[np] = endphotdr;
    }

    evt_.PhotonCand_mediumID[np] = evt_.PhotonCand_tightID[np] = 0;
    const std::vector<pat::Photon::IdPair> ids = photon->photonIDs();
    for ( unsigned int j = 0; j < ids.size(); ++j ) {
      pat::Photon::IdPair idp = ids.at( j );
      //FIXME make me private attributes
      if ( phoMediumIdLabel_.find( idp.first ) != std::string::npos ) evt_.PhotonCand_mediumID[np] = idp.second;
      if ( phoTightIdLabel_.find( idp.first ) != std::string::npos ) evt_.PhotonCand_tightID[np] = idp.second;
    }

    //edm::RefToBase<pat::Photon> photonRef = photonColl->refAt( i );
    //const edm::Ptr<reco::Photon> photonRef = photonColl->ptrAt( i );
    //evt_.PhotonCand_mediumID[np] = medium_id_decisions->operator[]( photonRef );
    //evt_.PhotonCand_tightID[np] = tight_id_decisions->operator[]( photonRef );

    np++;
  }
  LogDebug( "GammaGammaLL" ) << "Passed photons retrieval stage. Got " << np << " photon(s)";
}

void
GammaGammaLL::fetchProtons( const edm::Event& iEvent )
{
  // Forward proton tracks
  edm::Handle<edm::DetSetVector<TotemRPLocalTrack> > rplocaltracks;
  iEvent.getByToken( totemRPHitToken_, rplocaltracks);

  unsigned int nlpc = 0;
  for ( const auto rplocaltrack : *rplocaltracks ) {
    const CTPPSDetId det_id( rplocaltrack.detId() );
    const unsigned short arm = det_id.arm(), // 0->L, 1->R (  2/  3->L, 102/103->R)
                         pot = det_id.rp();  // 2->N, 3->F (  2/102->N,   3/103->F)

    for ( const auto track : rplocaltrack ) {
      if ( !track.isValid() ) continue;
      // express distances in metres
      evt_.LocalProtCand_x[nlpc] = track.getX0() * 1.e-3;
      evt_.LocalProtCand_y[nlpc] = track.getY0() * 1.e-3;
      evt_.LocalProtCand_z[nlpc] = track.getZ0() * 1.e-3;
      evt_.LocalProtCand_xSigma[nlpc] = track.getX0Sigma() * 1.e-3;
      evt_.LocalProtCand_ySigma[nlpc] = track.getY0Sigma() * 1.e-3;
      evt_.LocalProtCand_Tx[nlpc] = track.getTx();
      evt_.LocalProtCand_Ty[nlpc] = track.getTy();
      evt_.LocalProtCand_TxSigma[nlpc] = track.getTxSigma();
      evt_.LocalProtCand_TySigma[nlpc] = track.getTySigma();
      evt_.LocalProtCand_arm[nlpc] = arm;
      evt_.LocalProtCand_pot[nlpc] = pot;
      nlpc++;
      LogDebug( "GammaGammaLL" ) << "Proton track candidate with origin: ( " << track.getX0() << ", " << track.getY0() << ", " << track.getZ0() << " ) extracted!";
    }
  }
  LogDebug( "GammaGammaLL" ) << "Passed TOTEM RP info retrieval stage. Got " << nlpc << " local track(s)";
}

void
GammaGammaLL::fetchJets( const edm::Event& iEvent )
{
  // Get the Jet collection from the event
  edm::Handle<edm::View<pat::Jet> > jetColl; // PAT
  iEvent.getByToken( jetToken_, jetColl );

  double totalJetEnergy = 0., HEJet_pt = 0., HEJet_eta = 0., HEJet_phi = 0., HEJet_e = 0.;

  unsigned int njc = 0;
  for ( unsigned int i = 0; i < jetColl->size(); ++i ) {
    const edm::Ptr<pat::Jet> jet = jetColl->ptrAt( i );

    evt_.JetCand_e[njc] = jet->energy();
    evt_.JetCand_pt[njc] = jet->pt();
    evt_.JetCand_eta[njc] = jet->eta();
    evt_.JetCand_phi[njc] = jet->phi();
    totalJetEnergy += evt_.JetCand_e[njc];
    // Find kinematics quantities associated to the highest energy jet
    if ( evt_.JetCand_e[njc] > HEJet_e ) {
      HEJet_e = evt_.JetCand_e[njc];
      HEJet_pt = evt_.JetCand_pt[njc];
      HEJet_eta = evt_.JetCand_eta[njc];
      HEJet_phi = evt_.JetCand_phi[njc];
    }
    njc++;
  }
  evt_.HighestJet_pt = HEJet_pt;
  evt_.HighestJet_eta = HEJet_eta;
  evt_.HighestJet_phi = HEJet_phi;
  evt_.HighestJet_e = HEJet_e;
  evt_.SumJet_e = totalJetEnergy;

  LogDebug( "GammaGammaLL" ) << "Passed Loop on jets";
}

void
GammaGammaLL::fetchVertices( const edm::Event& iEvent )
{
  // Get the vertex collection from the event
  edm::Handle<edm::View<reco::Vertex> > recoVertexColl;
  iEvent.getByToken( recoVertexToken_, recoVertexColl );

  unsigned int nvtx = 0;
  for ( unsigned int i = 0; i < recoVertexColl->size(); ++i ) {
    const edm::Ptr<reco::Vertex> vertex = recoVertexColl->ptrAt( i );

    evt_.PrimVertexCand_id[nvtx] = nvtx;
    evt_.PrimVertexCand_x[nvtx] = vertex->x();
    evt_.PrimVertexCand_y[nvtx] = vertex->y();
    evt_.PrimVertexCand_z[nvtx] = vertex->z();
    evt_.PrimVertexCand_tracks[nvtx] = vertex->nTracks();
    evt_.PrimVertexCand_chi2[nvtx] = vertex->chi2();
    evt_.PrimVertexCand_ndof[nvtx] = vertex->ndof();
    nvtx++;
  }
}

void
GammaGammaLL::newVertexInfoRetrieval( const edm::Event& iEvent )
{
  iEvent.getByToken( recoTrackToken_, trackColl_ );

  foundPairInEvent_ = false;
  switch ( leptonsType_ ) {
    case ggll::ElectronMuon: {
      for ( unsigned int i = 0; i < evt_.MuonCand_pt.size(); ++i ) {
        for ( unsigned int j = 0; j < evt_.EleCand_et.size(); ++j ) {
          if ( evt_.MuonCand_charge[i]*evt_.EleCand_charge[j] > 0 ) continue;
          if ( newTracksInfoRetrieval( i, j ) ) foundPairInEvent_ = true;
        }
      }
    } break;
    case ggll::DiElectron: {
      for ( unsigned int i = 0; i < evt_.EleCand_et.size(); ++i ) {
        for ( unsigned int j = i+1; j < evt_.EleCand_et.size(); ++j ) {
          if ( evt_.EleCand_charge[i]*evt_.EleCand_charge[j] > 0 ) continue;
          if ( newTracksInfoRetrieval( i, j ) ) foundPairInEvent_ = true;
        }
      }
    } break;
    case ggll::DiMuon: {
      for ( unsigned int i = 0; i < evt_.MuonCand_pt.size(); ++i ) {
        for ( unsigned int j = i+1; j < evt_.MuonCand_pt.size(); ++j ) {
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

  const unsigned int npair = evt_.Pair_mass.size();

  switch ( leptonsType_ ) {
    case ggll::ElectronMuon: {
      translepttrks.push_back( muonTransientTracks_[l1id] );
      translepttrks.push_back( eleTransientTracks_[l2id] );
      l1cand_pt = evt_.MuonCand_innerTrackPt[l1id];
      l2cand_pt = evt_.EleCand_innerTrackPt[l2id];
    } break;
    case ggll::DiElectron: {
      translepttrks.push_back( eleTransientTracks_[l1id] );
      translepttrks.push_back( eleTransientTracks_[l2id] );
      l1cand_pt = evt_.EleCand_innerTrackPt[l1id];
      l2cand_pt = evt_.EleCand_innerTrackPt[l2id];
    } break;
    case ggll::DiMuon: {
      translepttrks.push_back( muonTransientTracks_[l1id] );
      translepttrks.push_back( muonTransientTracks_[l2id] );
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
    dileptonVertex = fitter.vertex( translepttrks );
  } catch (...) { return false; }

  if ( !dileptonVertex.isValid() ) return false; // only keep the pairs with valid vertex

  evt_.KalmanVertexCand_x[npair] = dileptonVertex.position().x();
  evt_.KalmanVertexCand_y[npair] = dileptonVertex.position().y();
  evt_.KalmanVertexCand_z[npair] = dileptonVertex.position().z();

  // Count nearby tracks
  int closesttrkid = -1, closesthighpuritytrkid = -1;
  double closesttrkdxyz = 999., closesthighpuritytrkdxyz = 999.;

  for ( unsigned int i = 0; i < trackColl_->size(); ++i ) {
    const edm::Ptr<reco::Track> track = trackColl_->ptrAt( i );

    if ( track->pt() == l1cand_pt || track->pt() == l2cand_pt) continue;

    const double vtxdst = sqrt( std::pow( ( track->vertex().x()-evt_.KalmanVertexCand_x[npair] ), 2 )
                              + std::pow( ( track->vertex().y()-evt_.KalmanVertexCand_y[npair] ), 2 )
                              + std::pow( ( track->vertex().z()-evt_.KalmanVertexCand_z[npair] ), 2 ) );
    if ( vtxdst < 0.05 ) evt_.Pair_extratracks0p5mm[npair]++;
    if ( vtxdst < 0.1 ) evt_.Pair_extratracks1mm[npair]++;
    if ( vtxdst < 0.2 ) evt_.Pair_extratracks2mm[npair]++;
    if ( vtxdst < 0.3 ) evt_.Pair_extratracks3mm[npair]++;
    if ( vtxdst < 0.4 ) evt_.Pair_extratracks4mm[npair]++;
    if ( vtxdst < 0.5 ) evt_.Pair_extratracks5mm[npair]++;
    if ( vtxdst < 1.0 ) evt_.Pair_extratracks1cm[npair]++;
    if ( vtxdst < 2.0 ) evt_.Pair_extratracks2cm[npair]++;
    if ( vtxdst < 3.0 ) evt_.Pair_extratracks3cm[npair]++;
    if ( vtxdst < 4.0 ) evt_.Pair_extratracks4cm[npair]++;
    if ( vtxdst < 5.0 ) evt_.Pair_extratracks5cm[npair]++;
    if ( vtxdst < 10. ) evt_.Pair_extratracks10cm[npair]++;
    if ( vtxdst < closesttrkdxyz ) {
      closesttrkid = i;
      closesttrkdxyz = vtxdst;
    }
    if ( track->quality( reco::TrackBase::highPurity ) == 1 && vtxdst < closesthighpuritytrkdxyz ) {
      closesthighpuritytrkid = i;
      closesthighpuritytrkdxyz = vtxdst;
    }

    // Save track properties if within 5mm
    if ( vtxdst < 0.5 ) {
      evt_.ExtraTrack_pair[evt_.nExtraTracks] = npair;
      evt_.ExtraTrack_purity[evt_.nExtraTracks] = track->quality( reco::TrackBase::highPurity );
      evt_.ExtraTrack_nhits[evt_.nExtraTracks] = track->numberOfValidHits();

      evt_.ExtraTrack_px[evt_.nExtraTracks] = track->px();
      evt_.ExtraTrack_py[evt_.nExtraTracks] = track->py();
      evt_.ExtraTrack_pz[evt_.nExtraTracks] = track->pz();
      evt_.ExtraTrack_charge[evt_.nExtraTracks] = track->charge();
      evt_.ExtraTrack_chi2[evt_.nExtraTracks] = track->chi2();
      evt_.ExtraTrack_ndof[evt_.nExtraTracks] = track->ndof();
      evt_.ExtraTrack_vtxdxyz[evt_.nExtraTracks] = vtxdst;
      evt_.ExtraTrack_vtxT[evt_.nExtraTracks] = sqrt( std::pow( track->vertex().x()-evt_.KalmanVertexCand_x[npair], 2)
                                                    + std::pow( track->vertex().y()-evt_.KalmanVertexCand_y[npair], 2 ) );
      evt_.ExtraTrack_vtxZ[evt_.nExtraTracks] = fabs( track->vertex().z()-evt_.KalmanVertexCand_z[npair] );
      evt_.ExtraTrack_x[evt_.nExtraTracks] = track->vertex().x();
      evt_.ExtraTrack_y[evt_.nExtraTracks] = track->vertex().y();
      evt_.ExtraTrack_z[evt_.nExtraTracks] = track->vertex().z();

      evt_.nExtraTracks++;
    }
  }
  evt_.ClosestExtraTrack_vtxdxyz[npair] = closesttrkdxyz;
  evt_.ClosestExtraTrack_id[npair] = closesttrkid;
  evt_.ClosestHighPurityExtraTrack_vtxdxyz[npair] = closesthighpuritytrkdxyz;
  evt_.ClosestHighPurityExtraTrack_id[npair] = closesthighpuritytrkid;

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

  const TLorentzVector pair(l1+l2);

  evt_.Pair_pt[npair] = pair.Pt();
  evt_.Pair_mass[npair] = pair.M();
  evt_.Pair_phi[npair] = pair.Phi();
  evt_.Pair_eta[npair] = pair.Eta();
  evt_.Pair_lepton1[npair] = l1id;
  evt_.Pair_lepton2[npair] = l2id;

  double dphi = fabs( l1.Phi()-l2.Phi() );
  // dphi lies in [-pi, pi]
  while ( dphi < -M_PI ) { dphi += 2.*M_PI; }
  while ( dphi >  M_PI ) { dphi -= 2.*M_PI; }
  evt_.Pair_dphi[npair] = dphi;

  evt_.Pair_dpt[npair] = fabs( l1.Pt()-l2.Pt() );
  evt_.Pair_3Dangle[npair] = l1.Angle( l2.Vect() )/M_PI;

  for ( unsigned int j = 0; j < evt_.nPhotonCand; ++j ) {
    TLorentzVector pho;
    pho.SetPtEtaPhiE( evt_.PhotonCand_pt[j], evt_.PhotonCand_eta[j], evt_.PhotonCand_phi[j], evt_.PhotonCand_e[j] );
    evt_.PairGamma_pair.emplace_back( npair );
    evt_.PairGamma_pho.emplace_back( j );
    evt_.PairGamma_mass.emplace_back( ( l1+l2+pho ).M() );
  }

  nCandidates_++;

  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void
GammaGammaLL::beginJob()
{
  // Filling the ntuple
  evt_.attach( tree_, leptonsType_, runOnMC_ );
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
  if ( !hltPrescale_.init( iRun, iSetup, hltMenuLabel_, changed ) ) {
    throw cms::Exception( "GammaGammaLL" ) << " prescales extraction failure with process name " << hltMenuLabel_;
  }
  // Initialise HLTConfigProvider
  hltConfig_ = hltPrescale_.hltConfigProvider();
  if ( !hltConfig_.init( iRun, iSetup, hltMenuLabel_, changed ) ) {
    throw cms::Exception( "GammaGammaLL" ) << " config extraction failure with process name " << hltMenuLabel_;
  }
  else if ( hltConfig_.size() == 0 ) {
    edm::LogError( "GammaGammaLL" ) << "HLT config size error";
  }
}

// ------------ method called when ending the processing of a run  ------------
void
GammaGammaLL::endRun( const edm::Run&, const edm::EventSetup& )
{
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
