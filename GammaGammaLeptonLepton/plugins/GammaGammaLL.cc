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

//
// constructors and destructor
//
GammaGammaLL::GammaGammaLL(const edm::ParameterSet& iConfig) :
  tree_(0),
  fetchMuons_(false), fetchElectrons_(false),
  fetchProtons_       (iConfig.getParameter<bool>("fetchProtons")),
  hltMenuLabel_       (iConfig.getParameter<std::string>("HLTMenuTag")),
  triggersList_       (iConfig.getParameter<std::vector<std::string> >("triggersList")),
  triggerResultsToken_(consumes<edm::TriggerResults>                   (iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  pileupToken_        (consumes< edm::View<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupInfo"))),
  recoVertexToken_    (consumes< edm::View<reco::Vertex> >             (iConfig.getParameter<edm::InputTag>("vertexTag"))),
  recoTrackToken_     (consumes< edm::View<reco::Track> >              (iConfig.getParameter<edm::InputTag>("trackTag"))),
  muonToken_          (consumes< edm::View<pat::Muon> >                (iConfig.getParameter<edm::InputTag>("muonTag"))),
  eleToken_           (consumes< edm::View<pat::Electron> >            (iConfig.getParameter<edm::InputTag>("electronTag"))),
  eleLooseIdMapToken_ (consumes< edm::ValueMap<bool> >                 (iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes< edm::ValueMap<bool> >                 (iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_ (consumes< edm::ValueMap<bool> >                 (iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  jetToken_           (consumes< edm::View<pat::Jet> >                 (iConfig.getParameter<edm::InputTag>("jetTag"))),
  metToken_           (consumes< edm::View<pat::MET> >                 (iConfig.getParameter<edm::InputTag>("metTag"))),
  totemRPHitToken_    (consumes< edm::DetSetVector<TotemRPLocalTrack> >(iConfig.getParameter<edm::InputTag>("totemRPLocalTrackTag"))),
  photonToken_        (consumes< edm::View<pat::Photon> >              (iConfig.getParameter<edm::InputTag>("photonTag"))),
  runOnMC_            (iConfig.getUntrackedParameter<bool>("runOnMC", false)),
  printCandidates_    (iConfig.getUntrackedParameter<bool>("printCandidates", false)),
  sqrts_              (iConfig.getParameter<double>("sqrtS")),
  useLegacyVertexing_ (iConfig.getParameter<bool>("useLegacyVertexing")),
  maxExTrkVtx_        (iConfig.getUntrackedParameter<unsigned int>("maxExtraTracks", MAX_ET)),
  hltPrescale_        (iConfig, consumesCollector(), *this),
  // Pileup input tags
  mcPileupFile_       (iConfig.getParameter<std::string>("mcpufile")),
  dataPileupFile_     (iConfig.getParameter<std::string>("datapufile")),
  mcPileupPath_       (iConfig.getParameter<std::string>("mcpupath")),
  dataPileupPath_     (iConfig.getParameter<std::string>("datapupath"))
{
  //now do what ever initialization is needed

  hlts_ = HLTMatcher(triggersList_);
  nHLT = triggersList_.size();	
  
  // Generator level
  if (runOnMC_) {
    genToken_ = consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"));
    minPtMC_ = iConfig.getUntrackedParameter<double>("MCAcceptPtCut", 10.);
    minEtaMC_ = iConfig.getUntrackedParameter<double>("MCAcceptEtaCut", 2.5);
  }

  // Leptons input tags
  const std::string ltype = iConfig.getParameter<std::string>("leptonsType");
  if (ltype=="ElectronMuon") {
    leptonsType_ = GammaGammaLL::ElectronMuon;
    fetchElectrons_ = true;
    fetchMuons_ = true;
  }
  else if (ltype=="Muon") {
    leptonsType_ = GammaGammaLL::Dimuon;
    fetchMuons_ = true;
  }
  else if (ltype=="Electron") {
    leptonsType_ = GammaGammaLL::Dielectron;
    fetchElectrons_ = true;
  }
  else throw cms::Exception("GammaGammaLL") << "'LeptonsType' parameter should either be:\n"
                                            << "   * 'ElectronMuon'       (for mixed leptons pair)\n"
                                            << "   * 'Electron' or 'Muon' (for same-flavour leptons)";

  // Pileup reweighting utilities
  if (runOnMC_) {
    lumiWeights_ = edm::LumiReWeighting(mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_);
  }

  // Book the output tree
  usesResource( "TFileService" );
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntp1", "ntp1");
}

GammaGammaLL::~GammaGammaLL()
{}

//
// member functions
//

void
GammaGammaLL::lookAtTriggers(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get the trigger information from the event
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken(triggerResultsToken_, hltResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);

  std::ostringstream os;
  os << "Trigger names: " << std::endl;
  for (unsigned int i=0; i<trigNames.size(); i++) {
    os << "--> " << trigNames.triggerNames().at(i) << std::endl;

    const int trigNum = hlts_.TriggerNum(trigNames.triggerNames().at(i));
    if (trigNum<0) continue; // Trigger didn't match the interesting ones

    HLT_Accept[trigNum] = hltResults->accept(i);

    // extract prescale value for this path
    if (runOnMC_) { HLT_Prescl[trigNum] = 1.; continue; } //FIXME
    int prescale_set = hltPrescale_.prescaleSet(iEvent, iSetup);
    HLT_Prescl[trigNum] = (prescale_set<0) ? 0. : hltConfig_.prescaleValue(prescale_set, trigNames.triggerNames().at(i)); //FIXME
  }
  LogDebug("GammaGammaLL") << os.str();
}

// ------------ method called for each event  ------------
void
GammaGammaLL::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Kalman filtering
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", KalVtx_);

  // First initialization of the variables
  clearTree();

  Weight = 1.;
  
  // Run and BX information
  BX = iEvent.bunchCrossing();
  Run = iEvent.id().run();
  LumiSection = iEvent.luminosityBlock();
  EventNum = iEvent.id().event();
  
  // High level trigger information retrieval  
  lookAtTriggers(iEvent, iSetup);
  
  LogDebug("GammaGammaLL") << "Passed trigger filtering stage";

  // Generator level information
  if (runOnMC_) {
    analyzeMCEventContent(iEvent);

    // Pileup information
    //const edm::EventBase* iEventB = static_cast<const edm::EventBase*>(&iEvent);
    //Weight = lumiWeights_.weight(*iEventB);

    edm::Handle< edm::View<PileupSummaryInfo> > pu_info;
    iEvent.getByToken(pileupToken_, pu_info);

    int npv0true = 0;
    for (unsigned int i=0; i<pu_info->size(); i++) {
      edm::Ptr<PileupSummaryInfo> PVI = pu_info->ptrAt(i);

      const int beamXing = PVI->getBunchCrossing(),
                npvtrue = PVI->getTrueNumInteractions();

      if(beamXing == 0) npv0true += npvtrue;
    }

    Weight = lumiWeights_.weight( npv0true );
    LogDebug("GammaGammaLL") << "Passed Pileup retrieval stage";
  } // run on MC?

  muonsMomenta_.clear();
  electronsMomenta_.clear();
  muonTransientTracks_.clear();
  eleTransientTracks_.clear();

  if (fetchMuons_) fetchMuons(iEvent);
  if (fetchElectrons_) fetchElectrons(iEvent);

  fetchPhotons(iEvent);

  foundPairInEvent_ = false;

  newVertexInfoRetrieval(iEvent);
  //legacyVertexInfoRetrieval(iEvent);

  /////
  if (!foundPairInEvent_) return; // avoid to unpack RP/jet/MET if no dilepton candidate has been found
  /////
  
  if (fetchProtons_) fetchProtons(iEvent);


  // Get the Jet collection from the event
  // PAT

  edm::Handle<edm::View<pat::Jet> > jetColl;
  iEvent.getByToken(jetToken_, jetColl);

  double totalJetEnergy = 0.,
         HEJet_pt = 0., HEJet_eta = 0., HEJet_phi = 0., HEJet_e = 0.;

  for (unsigned int i=0; i<jetColl->size(); i++) {
    edm::Ptr<pat::Jet> jet = jetColl->ptrAt(i);

    JetCand_e[nJetCand] = jet->energy();
    JetCand_pt[nJetCand] = jet->pt();
    JetCand_eta[nJetCand] = jet->eta();
    JetCand_phi[nJetCand] = jet->phi();
    totalJetEnergy += JetCand_e[nJetCand];
    // Find kinematics quantities associated to the highest energy jet
    if (JetCand_e[nJetCand]>HEJet_e) {
      HEJet_e = JetCand_e[nJetCand];
      HEJet_pt = JetCand_pt[nJetCand];
      HEJet_eta = JetCand_eta[nJetCand];
      HEJet_phi = JetCand_phi[nJetCand];
    }
    nJetCand++;
  }
  HighestJet_pt = HEJet_pt;
  HighestJet_eta = HEJet_eta;
  HighestJet_phi = HEJet_phi;
  HighestJet_e = HEJet_e;
  SumJet_e = totalJetEnergy;

  LogDebug("GammaGammaLL") << "Passed Loop on jets";

  // Missing ET
  edm::Handle< edm::View<pat::MET> > MET;
  iEvent.getByToken(metToken_, MET); 
  const edm::View<pat::MET>* metColl = MET.product(); 
  edm::View<pat::MET>::const_iterator met = metColl->begin();
  
  Etmiss = met->et();
  Etmiss_phi = met->phi();
  Etmiss_significance = met->significance();

  LogDebug("GammaGammaLL") << "Passed MET retrieval stage";

  if (printCandidates_) {
    std::cout << "Event " << Run << ":" << EventNum << " has " << nCandidatesInEvent << " leptons pair(s) candidate(s) (vertex mult. : " << nPrimVertexCand << ")" << std::endl;
  }

  tree_->Fill();
}

void
GammaGammaLL::analyzeMCEventContent(const edm::Event& iEvent)
{
  edm::Handle< edm::View<reco::GenParticle> > genPartColl;
  iEvent.getByToken(genToken_, genPartColl);
 
  for (reco::GenParticleCollection::const_iterator genPart=genPartColl->begin(); genPart!=genPartColl->end(); genPart++) {

    if (!genPart->isPromptFinalState()) continue;

    // check the particles out of acceptance
    if (genPart->pt()<minPtMC_ || (minEtaMC_!=-1. && fabs(genPart->eta())>minEtaMC_)) {
      if (fabs(genPart->pdgId())==13) nGenMuonCandOutOfAccept++;
      if (fabs(genPart->pdgId())==11) nGenEleCandOutOfAccept++;
      if (fabs(genPart->pdgId())==22) nGenPhotCandOutOfAccept++;
      if (genPart->pdgId()!=2212) continue; // we keep the forward protons
    }

    // generated outgoing protons
    if (genPart->pdgId()==2212 && nGenProtCand<MAX_GENPRO) {
      GenProtCand_pt[nGenProtCand] = genPart->pt();
      GenProtCand_eta[nGenProtCand] = genPart->eta();
      GenProtCand_phi[nGenProtCand] = genPart->phi();
      GenProtCand_e[nGenProtCand] = genPart->energy();
      GenProtCand_status[nGenProtCand] = genPart->status();
      nGenProtCand++;
    }

    // generated central dimuon system
    if (fabs(genPart->pdgId())==13 && nGenMuonCand<MAX_GENMU) {
      GenMuonCand_pt[nGenMuonCand] = genPart->pt();
      GenMuonCand_eta[nGenMuonCand] = genPart->eta();
      GenMuonCand_phi[nGenMuonCand] = genPart->phi();
      GenMuonCand_e[nGenMuonCand] = genPart->energy();
      nGenMuonCand++;
    }

    // generated central dimuon system
    if (fabs(genPart->pdgId())==11 && nGenEleCand<MAX_GENELE) {
      GenEleCand_pt[nGenEleCand] = genPart->pt();
      GenEleCand_eta[nGenEleCand] = genPart->eta();
      GenEleCand_phi[nGenEleCand] = genPart->phi(); 
      GenEleCand_e[nGenEleCand] = genPart->energy();
      nGenEleCand++;
    }

    // generated inner photon line
    if (genPart->pdgId()==22 && nGenPhotCand<MAX_GENPHO) {
      GenPhotCand_e[nGenPhotCand] = genPart->energy();
      GenPhotCand_pt[nGenPhotCand] = genPart->pt();
      GenPhotCand_eta[nGenPhotCand] = genPart->eta();
      GenPhotCand_phi[nGenPhotCand] = genPart->phi();
      nGenPhotCand++;
    }
    if (genPart->pdgId()==2212 && fabs(genPart->pz())>3000.) {
      // Kinematic quantities computation
      // xi = fractional momentum loss
      if (genPart->pz()>0.) xi = 1.-genPart->pz()/sqrts_;
      else xi = 1.+genPart->pz()/sqrts_;
      t = -(std::pow(genPart->pt(), 2)+std::pow(genPart->mass()*xi, 2))/(1.-xi);
    }
  }

  bool foundGenCandPairInEvent = false;

  TLorentzVector l1, l2;
  // electron+muon
  if (leptonsType_==ElectronMuon and nGenMuonCand==1 and nGenEleCand==1) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE(GenMuonCand_pt[0], GenMuonCand_eta[0], GenMuonCand_phi[0], GenMuonCand_e[0]);
    l2.SetPtEtaPhiE(GenEleCand_pt[0], GenEleCand_eta[0], GenEleCand_phi[0], GenEleCand_e[0]);
    foundGenCandPairInEvent = true;
  }
  // dielectron
  else if (leptonsType_==Dielectron and nGenEleCand==2) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE(GenEleCand_pt[0], GenEleCand_eta[0], GenEleCand_phi[0], GenEleCand_e[0]);
    l2.SetPtEtaPhiE(GenEleCand_pt[1], GenEleCand_eta[1], GenEleCand_phi[1], GenEleCand_e[1]);
    foundGenCandPairInEvent = true;
  }
  // dimuon
  else if (leptonsType_==Dimuon and nGenMuonCand==2) { // FIXME maybe a bit tight according to the newer PU conditions?
    l1.SetPtEtaPhiE(GenMuonCand_pt[0], GenMuonCand_eta[0], GenMuonCand_phi[0], GenMuonCand_e[0]);
    l2.SetPtEtaPhiE(GenMuonCand_pt[1], GenMuonCand_eta[1], GenMuonCand_phi[1], GenMuonCand_e[1]);
    foundGenCandPairInEvent = true;      	
  }
  if (foundGenCandPairInEvent) {
    const TLorentzVector pair = l1+l2;
    GenPair_mass = pair.M();
    GenPair_pt = pair.Pt();
    GenPair_eta = pair.Eta();
    GenPair_phi = pair.Phi();

    double dphi = fabs(l1.Phi()-l2.Phi());
    // dphi lies in [-pi, pi]
    while (dphi<-pi) { dphi += 2.*pi; }
    while (dphi> pi) { dphi -= 2.*pi; }
    GenPair_dphi = dphi;

    GenPair_dpt = fabs(l1.Pt()-l2.Pt());
    GenPair_3Dangle = (l1.Angle(l2.Vect()))/pi;
  }
}

void
GammaGammaLL::fetchMuons(const edm::Event& iEvent)
{
  // Get the muons collection from the event
  // PAT muons
  edm::Handle<edm::View<pat::Muon> > muonColl;
  edm::View<pat::Muon>::const_iterator muon;

  iEvent.getByToken(muonToken_, muonColl);

  for (unsigned int i=0; i<muonColl->size() && nMuonCand<MAX_MUONS; i++) {
    const edm::Ptr<pat::Muon> muon = muonColl->ptrAt(i);

    MuonCand_pt[nMuonCand] = muon->pt();
    MuonCand_eta[nMuonCand] = muon->eta();
    MuonCand_phi[nMuonCand] = muon->phi();
    MuonCand_e[nMuonCand] = muon->energy();
    MuonCand_charge[nMuonCand] = muon->charge();
    MuonCand_dxy[nMuonCand] = muon->dB();
    MuonCand_nstatseg[nMuonCand] = muon->numberOfMatchedStations();

    MuonCand_vtxx[nMuonCand] = muon->vertex().x();
    MuonCand_vtxy[nMuonCand] = muon->vertex().y();
    MuonCand_vtxz[nMuonCand] = muon->vertex().z();

    MuonCand_isglobal[nMuonCand] = muon->isGlobalMuon();
    MuonCand_istracker[nMuonCand] = muon->isTrackerMuon();
    MuonCand_isstandalone[nMuonCand] = muon->isStandAloneMuon();
    MuonCand_ispfmuon[nMuonCand] = muon->isPFMuon();

    TLorentzVector leptonptmp;  
    leptonptmp.SetPtEtaPhiM(muon->pt(), muon->eta(), muon->phi(), muon->mass());

    if (MuonCand_istracker[nMuonCand]) {
      MuonCand_npxlhits[nMuonCand] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();
      MuonCand_ntrklayers[nMuonCand] = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      leptonptmp.SetPtEtaPhiM(muon->innerTrack()->pt(), muon->innerTrack()->eta(), muon->innerTrack()->phi(), muon->mass());
    }
    muonsMomenta_.insert(std::pair<int,TLorentzVector>(nMuonCand, leptonptmp));

    if (MuonCand_isglobal[nMuonCand] && MuonCand_istracker[nMuonCand]) {
      MuonCand_innerTrackPt[nMuonCand] = muon->innerTrack()->pt();
      MuonCand_innerTrackEta[nMuonCand] = muon->innerTrack()->eta();
      MuonCand_innerTrackPhi[nMuonCand] = muon->innerTrack()->phi();
      MuonCand_innerTrackVtxz[nMuonCand] = muon->innerTrack()->vertex().z();
      MuonCandTrack_nmuchits[nMuonCand] = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
      MuonCandTrack_chisq[nMuonCand] = muon->globalTrack()->normalizedChi2();
      const bool istight = ( MuonCand_ispfmuon[nMuonCand]
                        and (MuonCandTrack_chisq[nMuonCand]<10.)
                        and (MuonCandTrack_nmuchits[nMuonCand]>=1)
                        and (MuonCand_nstatseg[nMuonCand]>=2)
                        and (MuonCand_dxy[nMuonCand]<.2)
                        and (MuonCand_npxlhits[nMuonCand]>0)
                        and (MuonCand_ntrklayers[nMuonCand]>5));
      MuonCand_istight[nMuonCand] = istight;

      muonTransientTracks_[nMuonCand] = KalVtx_->build(*muon->innerTrack());
    }

    nMuonCand++;
  }
  nLeptonCand += nMuonCand;
  LogDebug("GammaGammaLL") << "Passed Muon retrieval stage. Got " << nMuonCand << " muon(s)";
}

void
GammaGammaLL::fetchElectrons(const edm::Event& iEvent)
{
  // Get the electrons collection from the event
  edm::Handle<edm::View<pat::Electron> > eleColl;
  iEvent.getByToken(eleToken_, eleColl);

  edm::Handle< edm::ValueMap<bool> > loose_id_decisions, medium_id_decisions, tight_id_decisions; 
  iEvent.getByToken(eleLooseIdMapToken_, loose_id_decisions);
  iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_, tight_id_decisions);

  for (unsigned int j=0; j<eleColl->size(); j++) {
    const edm::Ptr<pat::Electron> electron = eleColl->ptrAt(j);

    EleCand_et[nEleCand] = electron->et();
    EleCand_eta[nEleCand] = electron->eta();
    EleCand_phi[nEleCand] = electron->phi();
    EleCand_e[nEleCand] = electron->energy();
    EleCand_charge[nEleCand] = electron->charge();

    EleCand_vtxx[nEleCand] = electron->vertex().x();
    EleCand_vtxy[nEleCand] = electron->vertex().y();
    EleCand_vtxz[nEleCand] = electron->vertex().z();

    TLorentzVector leptonptmp;
    leptonptmp.SetPtEtaPhiM(electron->et(), electron->eta(), electron->phi(), electron->mass());

    if (electron->closestCtfTrackRef().isNonnull()) { // Only for pat::Electron
      EleCand_innerTrackPt[nEleCand] = electron->closestCtfTrackRef()->pt();
      EleCand_innerTrackEta[nEleCand] = electron->closestCtfTrackRef()->eta();
      EleCand_innerTrackPhi[nEleCand] = electron->closestCtfTrackRef()->phi();
      EleCand_innerTrackVtxz[nEleCand] = electron->closestCtfTrackRef()->vertex().z();
      leptonptmp.SetPtEtaPhiM(EleCand_innerTrackPt[nEleCand], EleCand_innerTrackEta[nEleCand], EleCand_innerTrackPhi[nEleCand], electron->mass());
      //eleTransientTrack[nEleCand] = KalVtx_->build(*electron->closestCtfTrackRef());
    }
    eleTransientTracks_[nEleCand] = KalVtx_->build(*electron->gsfTrack());

    electronsMomenta_.insert(std::pair<int,TLorentzVector>(nEleCand, leptonptmp));

    EleCand_deltaPhi[nEleCand] = electron->deltaPhiSuperClusterTrackAtVtx();
    EleCand_deltaEta[nEleCand] = electron->deltaEtaSuperClusterTrackAtVtx();
    EleCand_HoverE[nEleCand] = electron->hcalOverEcal();
    EleCand_trackiso[nEleCand] = electron->dr03TkSumPt() / electron->et();
    EleCand_ecaliso[nEleCand] = electron->dr03EcalRecHitSumEt() / electron->et();
    EleCand_hcaliso[nEleCand] = electron->dr03HcalTowerSumEt() / electron->et();
    EleCand_sigmaIetaIeta[nEleCand] = electron->sigmaIetaIeta();
    EleCand_convDist[nEleCand] = fabs(electron->convDist()); 
    EleCand_convDcot[nEleCand] = fabs(electron->convDcot()); 
    EleCand_ecalDriven[nEleCand] = electron->ecalDrivenSeed();

    //reco::GsfElectron* electron 
    EleCand_looseID[nEleCand] = (*loose_id_decisions)[electron];
    EleCand_mediumID[nEleCand] = (*medium_id_decisions)[electron];
    EleCand_tightID[nEleCand] = (*tight_id_decisions)[electron];

    nEleCand++;
  }
  nLeptonCand += nEleCand;
  LogDebug("GammaGammaLL") << "Passed Electron retrieval stage. Got " << nEleCand << " electron(s)";
}

void
GammaGammaLL::fetchPhotons(const edm::Event& iEvent)
{
  // Get the photons collection from the event
  edm::Handle< edm::View<pat::Photon> > photonColl;
  iEvent.getByToken(photonToken_, photonColl);

  for (unsigned int i=0; i<photonColl->size(); i++) {
    const edm::Ptr<pat::Photon> photon = photonColl->ptrAt(i);

    PhotonCand_pt[nPhotonCand] = photon->pt();
    PhotonCand_eta[nPhotonCand] = photon->eta();
    PhotonCand_phi[nPhotonCand] = photon->phi();
    PhotonCand_e[nPhotonCand] = photon->energy();
    PhotonCand_r9[nPhotonCand] = photon->r9();

    PhotonCand_drtrue[nPhotonCand] = -999.;
    PhotonCand_detatrue[nPhotonCand] = -999.;
    PhotonCand_dphitrue[nPhotonCand] = -999.;
    if (runOnMC_) {
      double photdr = 999., photdeta = 999., photdphi = 999.;
      double endphotdr = 999., endphotdeta = 999., endphotdphi = 999.;
      for (int j=0; j<nGenPhotCand; j++) { // matching with the 'true' photon object from MC
        photdeta = (PhotonCand_eta[nPhotonCand]-GenPhotCand_eta[j]);
        photdphi = (PhotonCand_phi[nPhotonCand]-GenPhotCand_phi[j]);
        photdr = sqrt(std::pow(photdeta, 2)+std::pow(photdphi, 2));
        if (photdr<endphotdr) {
          endphotdr = photdr;
          endphotdeta = photdeta;
          endphotdphi = photdphi;
        }
      }
      PhotonCand_detatrue[nPhotonCand] = endphotdeta;
      PhotonCand_dphitrue[nPhotonCand] = endphotdphi;
      PhotonCand_drtrue[nPhotonCand] = endphotdr;
    }
    nPhotonCand++;
    LogDebug("GammaGammaLL") << "Passed photons retrieval stage. Got " << nPhotonCand << " photon(s)";
  }
}

void
GammaGammaLL::fetchProtons(const edm::Event& iEvent)
{
  // Forward proton tracks
  edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rplocaltracks; 
  iEvent.getByToken(totemRPHitToken_, rplocaltracks);

  nLocalProtCand = 0;
  for (edm::DetSetVector<TotemRPLocalTrack>::const_iterator rplocaltrack=rplocaltracks->begin(); rplocaltrack!=rplocaltracks->end(); rplocaltrack++) {
    const unsigned int det_id = rplocaltrack->detId();
    const unsigned short arm = (det_id%100==2), // 0->F, 1->N (3/103->F, 2/102->N)
                         side = (det_id/100); // 0->L, 1->R (2/3->L, 102/103->R)
    for (edm::DetSet<TotemRPLocalTrack>::const_iterator proton=rplocaltrack->begin(); proton!=rplocaltrack->end(); proton++) {
      if (!proton->isValid()) continue;
      LocalProtCand_x[nLocalProtCand] = (proton->getX0())/1.e3;
      LocalProtCand_y[nLocalProtCand] = (proton->getY0())/1.e3;
      LocalProtCand_z[nLocalProtCand] = (proton->getZ0())/1.e3;
      LocalProtCand_xSigma[nLocalProtCand] = (proton->getX0Sigma())/1.e3;
      LocalProtCand_ySigma[nLocalProtCand] = (proton->getY0Sigma())/1.e3;
      LocalProtCand_Tx[nLocalProtCand] = proton->getTx();
      LocalProtCand_Ty[nLocalProtCand] = proton->getTy();
      LocalProtCand_TxSigma[nLocalProtCand] = proton->getTxSigma();
      LocalProtCand_TySigma[nLocalProtCand] = proton->getTySigma();
      LocalProtCand_arm[nLocalProtCand] = arm;
      LocalProtCand_side[nLocalProtCand] = side;
      nLocalProtCand++;
      LogDebug("GammaGammaLL") << "Proton track candidate with origin: (" << proton->getX0() << ", " << proton->getY0() << ", " << proton->getZ0() << ") extracted!";
    }
  }
  LogDebug("GammaGammaLL") << "Passed TOTEM RP info retrieval stage. Got " << nLocalProtCand << " local track(s)";
}

void
GammaGammaLL::newVertexInfoRetrieval(const edm::Event& iEvent)
{
  iEvent.getByToken(recoTrackToken_, trackColl_);

  if (leptonsType_==ElectronMuon) {
    for (int i=0; i<nMuonCand; i++) {
      for (int j=0; j<nEleCand; j++) {
        if (MuonCand_charge[i]*EleCand_charge[j]>0) continue;

        foundPairInEvent_ = newTracksInfoRetrieval(i, j);
      }
    }
  }
  else if (leptonsType_==Dielectron) {
    for (int i=0; i<nEleCand; i++) {
      for (int j=i+1; j<nEleCand; j++) {
        if (EleCand_charge[i]*EleCand_charge[j]>0) continue;

        foundPairInEvent_ = newTracksInfoRetrieval(i, j);
      }
    }
  }
  else if (leptonsType_==Dimuon) {
    for (int i=0; i<nMuonCand; i++) {
      for (int j=i+1; j<nMuonCand; j++) {
        if (MuonCand_charge[i]*MuonCand_charge[j]>0) continue;

        foundPairInEvent_ = newTracksInfoRetrieval(i, j);
      }
    }
  }

}

bool
GammaGammaLL::newTracksInfoRetrieval(int l1id, int l2id)
{
  double l1cand_pt, l2cand_pt;
  std::vector<reco::TransientTrack> translepttrks;

  if (leptonsType_==ElectronMuon) {
    translepttrks.push_back(muonTransientTracks_[l1id]);
    translepttrks.push_back(eleTransientTracks_[l2id]);
    l1cand_pt = MuonCand_innerTrackPt[l1id];
    l2cand_pt = EleCand_innerTrackPt[l2id];
  }
  else if (leptonsType_==Dielectron) {
    translepttrks.push_back(eleTransientTracks_[l1id]);
    translepttrks.push_back(eleTransientTracks_[l2id]);
    l1cand_pt = EleCand_innerTrackPt[l1id];
    l2cand_pt = EleCand_innerTrackPt[l2id];
  }
  else if (leptonsType_==Dimuon) {
    translepttrks.push_back(muonTransientTracks_[l1id]);
    translepttrks.push_back(muonTransientTracks_[l2id]);
    l1cand_pt = MuonCand_innerTrackPt[l1id];
    l2cand_pt = MuonCand_innerTrackPt[l2id];
  }
  else throw cms::Exception("GammaGammaLL") << "Invalid leptons type: " << leptonsType_;

  if (translepttrks.size()<2) return false; // just in case...

  const bool use_smoothing = true;
  KalmanVertexFitter fitter(use_smoothing); 
  TransientVertex dileptonVertex;
  try {
    dileptonVertex = fitter.vertex(translepttrks);
  } catch (...) { return false; }

  if (!dileptonVertex.isValid()) return false; // only keep the pairs with valid vertex

  KalmanVertexCand_x[nPair] = dileptonVertex.position().x(); 
  KalmanVertexCand_y[nPair] = dileptonVertex.position().y(); 
  KalmanVertexCand_z[nPair] = dileptonVertex.position().z(); 

  // Count nearby tracks
  int closesttrkid = -1, closesthighpuritytrkid = -1;
  double closesttrkdxyz = 999., closesthighpuritytrkdxyz = 999.;

  for (unsigned int i=0; i<trackColl_->size() && nExtraTracks<MAX_ET; i++) {
    edm::Ptr<reco::Track> track = trackColl_->ptrAt(i);

    if (track->pt()==l1cand_pt or track->pt()==l2cand_pt) continue;

    const double vtxdst = sqrt(std::pow((track->vertex().x()-KalmanVertexCand_x[nPair]),2)
                             + std::pow((track->vertex().y()-KalmanVertexCand_y[nPair]),2)
                             + std::pow((track->vertex().z()-KalmanVertexCand_z[nPair]),2));
    if (vtxdst<0.1) Pair_extratracks1mm[nPair]++;
    if (vtxdst<0.2) Pair_extratracks2mm[nPair]++;
    if (vtxdst<0.3) Pair_extratracks3mm[nPair]++;
    if (vtxdst<0.4) Pair_extratracks4mm[nPair]++;
    if (vtxdst<0.5) Pair_extratracks5mm[nPair]++;
    if (vtxdst<1.0) Pair_extratracks1cm[nPair]++;
    if (vtxdst<2.0) Pair_extratracks2cm[nPair]++;
    if (vtxdst<3.0) Pair_extratracks3cm[nPair]++;
    if (vtxdst<4.0) Pair_extratracks4cm[nPair]++;
    if (vtxdst<5.0) Pair_extratracks5cm[nPair]++;
    if (vtxdst<10.) Pair_extratracks10cm[nPair]++;
    if (vtxdst<closesttrkdxyz) {
      closesttrkid = i;
      closesttrkdxyz = vtxdst;
    }
    if (track->quality(reco::TrackBase::highPurity)==1 and vtxdst<closesthighpuritytrkdxyz) {
      closesthighpuritytrkid = i;
      closesthighpuritytrkdxyz = vtxdst;
    }

    // Save track properties if within 5mm
    if (vtxdst<0.5) {
      ExtraTrack_pair[nExtraTracks] = nPair;
      ExtraTrack_purity[nExtraTracks] = track->quality(reco::TrackBase::highPurity);
      ExtraTrack_nhits[nExtraTracks] = track->numberOfValidHits();

      ExtraTrack_px[nExtraTracks] = track->px();
      ExtraTrack_py[nExtraTracks] = track->py();
      ExtraTrack_pz[nExtraTracks] = track->pz();
      ExtraTrack_charge[nExtraTracks] = track->charge();
      ExtraTrack_chi2[nExtraTracks] = track->chi2();
      ExtraTrack_ndof[nExtraTracks] = track->ndof();
      ExtraTrack_vtxdxyz[nExtraTracks] = vtxdst;
      ExtraTrack_vtxT[nExtraTracks] = sqrt(std::pow(track->vertex().x()-KalmanVertexCand_x[nPair], 2)
                                        +  std::pow(track->vertex().y()-KalmanVertexCand_y[nPair], 2));
      ExtraTrack_vtxZ[nExtraTracks] = fabs(track->vertex().z()-KalmanVertexCand_z[nPair]);
      ExtraTrack_x[nExtraTracks] = track->vertex().x();  
      ExtraTrack_y[nExtraTracks] = track->vertex().y();  
      ExtraTrack_z[nExtraTracks] = track->vertex().z();  

      nExtraTracks++;
    }
  }
  ClosestExtraTrack_vtxdxyz[nPair] = closesttrkdxyz;
  ClosestExtraTrack_id[nPair] = closesttrkid;
  ClosestHighPurityExtraTrack_vtxdxyz[nPair] = closesthighpuritytrkdxyz;
  ClosestHighPurityExtraTrack_id[nPair] = closesthighpuritytrkid;

  TLorentzVector l1, l2;
  if (leptonsType_==ElectronMuon) {
    l1.SetPtEtaPhiE(MuonCand_pt[l1id], MuonCand_eta[l1id], MuonCand_phi[l1id], MuonCand_e[l1id]);
    l2.SetPtEtaPhiE(EleCand_et[l2id], EleCand_eta[l2id], EleCand_phi[l2id], EleCand_e[l2id]);
  }
  else if (leptonsType_==Dielectron) {
    l1.SetPtEtaPhiE(EleCand_et[l1id], EleCand_eta[l1id], EleCand_phi[l1id], EleCand_e[l1id]);
    l2.SetPtEtaPhiE(EleCand_et[l2id], EleCand_eta[l2id], EleCand_phi[l2id], EleCand_e[l2id]);
  }
  else if (leptonsType_==Dimuon) {
    l1.SetPtEtaPhiE(MuonCand_pt[l1id], MuonCand_eta[l1id], MuonCand_phi[l1id], MuonCand_e[l1id]);
    l2.SetPtEtaPhiE(MuonCand_pt[l2id], MuonCand_eta[l2id], MuonCand_phi[l2id], MuonCand_e[l2id]);
  }
  else throw cms::Exception("GammaGammaLL") << "Invalid leptons type: " << leptonsType_;

  const TLorentzVector pair(l1+l2);

  Pair_p[nPair] = pair.P();
  Pair_pt[nPair] = pair.Pt();
  Pair_mass[nPair] = pair.M();
  Pair_phi[nPair] = pair.Phi();
  Pair_eta[nPair] = pair.Eta();

  double dphi = fabs(l1.Phi()-l2.Phi());
  // dphi lies in [-pi, pi]
  while (dphi<-pi) { dphi += 2.*pi; }
  while (dphi> pi) { dphi -= 2.*pi; }
  Pair_dphi[nPair] = dphi;

  Pair_dpt[nPair] = fabs(l1.Pt()-l2.Pt());
  Pair_3Dangle[nPair] = (l1.Angle(l2.Vect()))/pi;

  for (int j=0; j<nPhotonCand; j++) {
    TLorentzVector pho; pho.SetPtEtaPhiE(PhotonCand_pt[j], PhotonCand_eta[j], PhotonCand_phi[j], PhotonCand_e[j]);
    PairGamma_mass[nPair][j] = (l1+l2+pho).M();
    //std::cout << "Photon " << j << " added to give a mass = " << PairGamma_mass[nPair][j] << std::endl;
  }

  nPair++;
  return true;
}

void
GammaGammaLL::legacyVertexInfoRetrieval(const edm::Event& iEvent)
{
  // Get the vertex collection from the event
  edm::Handle< edm::View<reco::Vertex> > recoVertexColl;
  iEvent.getByToken(recoVertexToken_, recoVertexColl);

  int vtxind = 0; // Primary vertex index (used in loop over vertices)

  // Vertex tracks multiplicity
  nPrimVertexCand = recoVertexColl->size();

  nExtraTracks = 0;

  if (nLeptonCand>=2) return;

  // Enough leptons candidates to go deeper and analyze the primary vertices 
  int etind = 0; // Extra tracks on vertex index (used in loop over tracks)

  PrimaryVertexSelector vtx(muonsMomenta_, electronsMomenta_);

  for (int i=0; i<nPrimVertexCand && vtxind<MAX_VTX; i++) {
    const edm::Ptr<reco::Vertex> vertex = recoVertexColl->ptrAt( i );

    PrimVertexCand_tracks[vtxind] = vertex->tracksSize();
    if (static_cast<unsigned int>(PrimVertexCand_tracks[vtxind]-2)>maxExTrkVtx_) continue; // cut on the upper number of extra tracks

    nExtraTracks = 0;
    nQualityExtraTrack = 0;

    // Find association between tracks and lepton kinematics
    vtx.feedTracks(vertex->tracks_begin(), vertex->tracks_end());

    const PrimaryVertexSelector::MatchedLeptonsMap mu_match = vtx.matchedMuons(),
                                                   ele_match = vtx.matchedElectrons();

    PrimVertexCand_id[vtxind] = vtxind;
    PrimVertexCand_x[vtxind] = vertex->x();
    PrimVertexCand_y[vtxind] = vertex->y();
    PrimVertexCand_z[vtxind] = vertex->z();
    PrimVertexCand_chi2[vtxind] = vertex->chi2();
    PrimVertexCand_ndof[vtxind] = vertex->ndof();

    PrimVertexCand_matchedtracks[vtxind] = vtx.matchedMuons().size()+vtx.matchedElectrons().size();

    Pair_lepton1[vtxind] = Pair_lepton2[vtxind] = -1;

    TLorentzVector l1, l2;
    bool foundPairOnVertex = false;
    double minDist = 999.;

    if (leptonsType_==ElectronMuon) {

      if (ele_match.size()==0 or mu_match.size()==0) { // not enough muons or electrons candidates on the vertex
        edm::LogWarning("GammaGammaLL") << "Not enough electrons (" << ele_match.size() << ") or muons (" << mu_match.size() << ") arising from the primary vertex !";
        continue;
      }

      for (unsigned int i=0; i<mu_match.size(); i++) {
        const int lep1 = mu_match.at(i).first;

        for (unsigned int j=0; j<ele_match.size(); j++) {
          const int lep2 = ele_match.at(j).first;
          if (MuonCand_charge[lep1]*EleCand_charge[lep2]>0) continue;
          foundPairOnVertex = true;
          const double leptonsDist = sqrt(pow(MuonCand_vtxx[lep1]-EleCand_vtxx[lep2],2)+
                                          pow(MuonCand_vtxy[lep1]-EleCand_vtxy[lep2],2)+
                                          pow(MuonCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
          if (leptonsDist<minDist) {
            minDist = leptonsDist;
            Pair_lepton1[vtxind] = lep1;
            Pair_lepton2[vtxind] = lep2;
          }
        }
      }
      l1.SetPtEtaPhiE(MuonCand_pt[Pair_lepton1[vtxind]], MuonCand_eta[Pair_lepton1[vtxind]], MuonCand_phi[Pair_lepton1[vtxind]], MuonCand_e[Pair_lepton1[vtxind]]);
      l2.SetPtEtaPhiE(EleCand_et[Pair_lepton2[vtxind]], EleCand_eta[Pair_lepton2[vtxind]], EleCand_phi[Pair_lepton2[vtxind]], EleCand_e[Pair_lepton2[vtxind]]);

      if (muonTransientTracks_.find(Pair_lepton1[vtxind])!=muonTransientTracks_.end()
        && eleTransientTracks_.find(Pair_lepton2[vtxind])!=eleTransientTracks_.end())
      {
        std::vector<reco::TransientTrack> trans_trks;
        trans_trks.push_back(muonTransientTracks_[Pair_lepton1[vtxind]]);
        trans_trks.push_back(eleTransientTracks_[Pair_lepton2[vtxind]]);

        KalmanVertexFitter fitter(true);
        TransientVertex tr_vtx = fitter.vertex(trans_trks);
        if (tr_vtx.isValid()) {
          KalmanVertexCand_x[vtxind] = tr_vtx.position().x();
          KalmanVertexCand_y[vtxind] = tr_vtx.position().y();
          KalmanVertexCand_z[vtxind] = tr_vtx.position().z();
        }
      }
    }
    else if (leptonsType_==Dielectron) {

      if (ele_match.size()<2) continue; // not enough electrons candidates on the vertex

      for (unsigned int i=0; i<ele_match.size(); i++) {
        const int lep1 = ele_match.at(i).first;
        for (unsigned int j=i+1; j<ele_match.size(); j++) {
          const int lep2 = ele_match.at(j).first;
          if (EleCand_charge[lep1]*EleCand_charge[lep2]>0) continue;
          foundPairOnVertex = true;
          const double leptonsDist = sqrt(pow(EleCand_vtxx[lep1]-EleCand_vtxx[lep2],2)+
                                          pow(EleCand_vtxy[lep1]-EleCand_vtxy[lep2],2)+
                                          pow(EleCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
          if (leptonsDist<minDist) {
            minDist = leptonsDist;
            Pair_lepton1[vtxind] = lep1;
            Pair_lepton2[vtxind] = lep2;
          }
        }
      }
      l1.SetPtEtaPhiE(EleCand_et[Pair_lepton1[vtxind]], EleCand_eta[Pair_lepton1[vtxind]], EleCand_phi[Pair_lepton1[vtxind]], EleCand_e[Pair_lepton1[vtxind]]);
      l2.SetPtEtaPhiE(EleCand_et[Pair_lepton2[vtxind]], EleCand_eta[Pair_lepton2[vtxind]], EleCand_phi[Pair_lepton2[vtxind]], EleCand_e[Pair_lepton2[vtxind]]);

      if (eleTransientTracks_.find(Pair_lepton1[vtxind])!=eleTransientTracks_.end()
       && eleTransientTracks_.find(Pair_lepton2[vtxind])!=eleTransientTracks_.end())
      {
        std::vector<reco::TransientTrack> trans_trks;
        trans_trks.push_back(eleTransientTracks_[Pair_lepton1[vtxind]]);
        trans_trks.push_back(eleTransientTracks_[Pair_lepton2[vtxind]]);

        KalmanVertexFitter fitter(true);
        TransientVertex tr_vtx = fitter.vertex(trans_trks);
        if (tr_vtx.isValid()) {
          KalmanVertexCand_x[vtxind] = tr_vtx.position().x();
          KalmanVertexCand_y[vtxind] = tr_vtx.position().y();
          KalmanVertexCand_z[vtxind] = tr_vtx.position().z();
        }
      }
    }
    else if (leptonsType_==Dimuon) { // Looks at dimuons
      if (vtx.matchedMuons().size()<2) continue; // not enough muons candidates on the vertex

      for (unsigned int i=0; i<mu_match.size(); i++) {
        const int lep1 = mu_match.at(i).first;

        for (unsigned int j=i+1; j<mu_match.size(); j++) {
          const int lep2 = mu_match.at(j).first;

          if (MuonCand_charge[lep1]*MuonCand_charge[lep2]>0) continue;
          foundPairOnVertex = true;
          const double leptonsDist = sqrt(pow(MuonCand_vtxx[lep1]-MuonCand_vtxx[lep2],2)+
                                          pow(MuonCand_vtxy[lep1]-MuonCand_vtxy[lep2],2)+
                                          pow(MuonCand_vtxz[lep1]-MuonCand_vtxz[lep2],2));
          if (leptonsDist<minDist) {
            minDist = leptonsDist;
            Pair_lepton1[vtxind] = lep1;
            Pair_lepton2[vtxind] = lep2;
          }
        }
      }
      l1.SetPtEtaPhiE(MuonCand_pt[Pair_lepton1[vtxind]], MuonCand_eta[Pair_lepton1[vtxind]], MuonCand_phi[Pair_lepton1[vtxind]], MuonCand_e[Pair_lepton1[vtxind]]);
      l2.SetPtEtaPhiE(MuonCand_pt[Pair_lepton2[vtxind]], MuonCand_eta[Pair_lepton2[vtxind]], MuonCand_phi[Pair_lepton2[vtxind]], MuonCand_e[Pair_lepton2[vtxind]]);

      if (muonTransientTracks_.find(Pair_lepton1[vtxind])!=muonTransientTracks_.end()
       && muonTransientTracks_.find(Pair_lepton2[vtxind])!=muonTransientTracks_.end())
      {
        std::vector<reco::TransientTrack> trans_trks;
        trans_trks.push_back(muonTransientTracks_[Pair_lepton1[vtxind]]);
        trans_trks.push_back(muonTransientTracks_[Pair_lepton2[vtxind]]);

        KalmanVertexFitter fitter(true);
        TransientVertex tr_vtx = fitter.vertex(trans_trks);
        if (tr_vtx.isValid()) {
          KalmanVertexCand_x[vtxind] = tr_vtx.position().x();
          KalmanVertexCand_y[vtxind] = tr_vtx.position().y();
          KalmanVertexCand_z[vtxind] = tr_vtx.position().z();
        }
      }
    }

    double closesttrkdxyz = 999., closesthighpuritytrkdxyz = 999., closesttrkdxyz_kalman = 999.;
    int closesttrkid = -1;
    // Loop on all the tracks matched with this vertex
    for (reco::Vertex::trackRef_iterator track=vertex->tracks_begin();
        track!=vertex->tracks_end() && etind<MAX_ET && etind<(int)maxExTrkVtx_; track++) {

      const reco::TrackRef trk = track->castTo<reco::TrackRef>();
      const int match_mu = vtx.matchedMuon(trk),
                match_ele = vtx.matchedElectron(trk);
      // check if track was matched to any of the leptons in the collection
      if (leptonsType_==GammaGammaLL::ElectronMuon    and (match_mu ==Pair_lepton1[vtxind] or match_ele==Pair_lepton2[vtxind])) continue;
      else if (leptonsType_==GammaGammaLL::Dimuon     and (match_mu ==Pair_lepton1[vtxind] or match_mu ==Pair_lepton2[vtxind])) continue;
      else if (leptonsType_==GammaGammaLL::Dielectron and (match_ele==Pair_lepton1[vtxind] or match_ele==Pair_lepton2[vtxind])) continue;

      ExtraTrack_pair[etind] = vtxind;
      const double vtxt = sqrt(std::pow(((*track)->vertex().x()-vertex->x()),2)
                             + std::pow(((*track)->vertex().y()-vertex->y()),2)),
                   vtxdst = sqrt(vtxt*vtxt
                               + std::pow(((*track)->vertex().z()-vertex->z()),2));

      ExtraTrack_purity[etind] = (*track)->quality(reco::TrackBase::highPurity);
      ExtraTrack_nhits[etind] = (*track)->numberOfValidHits();

      ExtraTrack_px[etind] = (*track)->px();
      ExtraTrack_py[etind] = (*track)->py();
      ExtraTrack_pz[etind] = (*track)->pz();
      ExtraTrack_charge[etind] = (*track)->charge();
      ExtraTrack_chi2[etind] = (*track)->chi2();
      ExtraTrack_ndof[etind] = (*track)->ndof();
      ExtraTrack_vtxdxyz[etind] = vtxdst;
      ExtraTrack_vtxT[etind] = vtxt;
      ExtraTrack_vtxZ[etind] = fabs((*track)->vertex().z()-vertex->z());
      ExtraTrack_x[etind] = (*track)->vertex().x();
      ExtraTrack_y[etind] = (*track)->vertex().y();
      ExtraTrack_z[etind] = (*track)->vertex().z();

      if (vtxdst<0.1) Pair_extratracks1mm[vtxind]++;
      if (vtxdst<0.2) Pair_extratracks2mm[vtxind]++;
      if (vtxdst<0.3) Pair_extratracks3mm[vtxind]++;
      if (vtxdst<0.4) Pair_extratracks4mm[vtxind]++;
      if (vtxdst<0.5) Pair_extratracks5mm[vtxind]++;
      if (vtxdst<1.0) Pair_extratracks1cm[vtxind]++;
      if (vtxdst<2.0) Pair_extratracks2cm[vtxind]++;
      if (vtxdst<3.0) Pair_extratracks3cm[vtxind]++;
      if (vtxdst<4.0) Pair_extratracks4cm[vtxind]++;
      if (vtxdst<5.0) Pair_extratracks5cm[vtxind]++;
      if (vtxdst<10.) Pair_extratracks10cm[vtxind]++;

      // minimum distance track to "simple" vertex
      if (vtxdst<closesttrkdxyz) {
        closesttrkdxyz = vtxdst;
        closesttrkid = etind;
      }

      // minimum distance track to associated Kalman vertex
      if (KalmanVertexCand_x[vtxind]>-999. and KalmanVertexCand_y[vtxind]>-999. and KalmanVertexCand_z[vtxind]>-999.) {
        const double vtxdst_kalman = sqrt(std::pow(((*track)->vertex().x()-KalmanVertexCand_x[vtxind]),2)+
                                          std::pow(((*track)->vertex().y()-KalmanVertexCand_y[vtxind]),2)+
                                          std::pow(((*track)->vertex().z()-KalmanVertexCand_z[vtxind]),2));
        if (vtxdst_kalman<closesttrkdxyz_kalman) {
          closesttrkdxyz_kalman = vtxdst_kalman;
          //closesttrkid_kalman = etind;
        }
      }

      // minimum distance high purity track to "simple" vertex
      if (ExtraTrack_purity[etind]==1 && ExtraTrack_nhits[etind]>=3) {
        nQualityExtraTrack++;
        if (vtxdst<closesthighpuritytrkdxyz) {
          closesthighpuritytrkdxyz = vtxdst;
          //closesthighpuritytrkid = etind;
        }
      }
      etind++;
    }
    ClosestExtraTrack_vtxdxyz[vtxind] = closesttrkdxyz;
    ClosestExtraTrack_id[vtxind] = closesttrkid;
    ClosestHighPurityExtraTrack_vtxdxyz[vtxind] = closesthighpuritytrkdxyz;
    ClosestHighPurityExtraTrack_id[vtxind] = closesttrkid;
    ClosestExtraTrackKalman_vtxdxyz[vtxind] = closesttrkdxyz_kalman;

    if (foundPairOnVertex) {
      Pair_mindist[vtxind] = minDist;

      std::ostringstream os; os << "---> Found a lepton pair on vertex!" << std::endl;
      os << "Matched muons : " << std::endl;
      for (PrimaryVertexSelector::MatchedLeptonsMap::const_iterator it=vtx.matchedMuons().begin(); it!=vtx.matchedMuons().end(); ++it) {
        os << "-> " << it->first << std::endl;
      }
      os << "Matched electrons : " << std::endl;
      for (PrimaryVertexSelector::MatchedLeptonsMap::const_iterator it=vtx.matchedElectrons().begin(); it!=vtx.matchedElectrons().end(); ++it) {
        os << "-> " << it->first << std::endl;
      }
      LogDebug("GammaGammaLL") << os.str();

      const TLorentzVector pair = l1+l2;
      Pair_p[vtxind] = pair.P();
      Pair_pt[vtxind] = pair.Pt();
      Pair_mass[vtxind] = pair.M();
      Pair_phi[vtxind] = pair.Phi();
      Pair_eta[vtxind] = pair.Eta();

      double dphi = fabs(l1.Phi()-l2.Phi());
      // dphi lies in [-pi, pi]
      while (dphi<-pi) { dphi += 2.*pi; }
      while (dphi> pi) { dphi -= 2.*pi; }
      Pair_dphi[vtxind] = dphi;

      Pair_dpt[vtxind] = fabs(l1.Pt()-l2.Pt());
      Pair_3Dangle[vtxind] = (l1.Angle(l2.Vect()))/pi;

      for (int j=0; j<nPhotonCand; j++) {
        TLorentzVector pho; pho.SetPtEtaPhiE(PhotonCand_pt[j], PhotonCand_eta[j], PhotonCand_phi[j], PhotonCand_e[j]);
        PairGamma_mass[vtxind][j] = (l1+l2+pho).M();
        //std::cout << "Photon " << j << " added to give a mass = " << PairGamma_mass[vtxind][j] << std::endl;
      }
      PrimVertexCand_hasdil[vtxind] = 1;

      nCandidatesInEvent++;
      nCandidates++;
      foundPairInEvent_ = true;
    }

    vtxind++;
  }
  nFilteredPrimVertexCand = vtxind;
  LogDebug("GammaGammaLL") << "Passed the loop on vertices : vtxind = " << vtxind;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaLL::beginJob()
{
  // Filling the ntuple

  tree_->Branch("Run", &Run, "Run/I");
  tree_->Branch("LumiSection", &LumiSection, "LumiSection/I");
  tree_->Branch("BX", &BX, "BX/I");
  tree_->Branch("EventNum", &EventNum, "EventNum/I");

  tree_->Branch("nHLT", &nHLT, "nHLT/I");
  tree_->Branch("HLT_Accept", HLT_Accept, "HLT_Accept[nHLT]/I");
  tree_->Branch("HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I");
  tree_->Branch("HLT_Name", &HLT_Name); *HLT_Name = triggersList_;
  
  if (fetchMuons_) {
    tree_->Branch("nMuonCand", &nMuonCand, "nMuonCand/I");
    tree_->Branch("MuonCand_pt", MuonCand_pt, "MuonCand_pt[nMuonCand]/D");
    tree_->Branch("MuonCand_eta", MuonCand_eta, "MuonCand_eta[nMuonCand]/D");
    tree_->Branch("MuonCand_phi", MuonCand_phi, "MuonCand_phi[nMuonCand]/D");
    tree_->Branch("MuonCand_e", MuonCand_e, "MuonCand_e[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackPt", MuonCand_innerTrackPt, "MuonCand_innerTrackPt[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackEta", MuonCand_innerTrackEta, "MuonCand_innerTrackEta[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackPhi", MuonCand_innerTrackPhi, "MuonCand_innerTrackPhi[nMuonCand]/D");
    tree_->Branch("MuonCand_charge", MuonCand_charge, "MuonCand_charge[nMuonCand]/I");
    tree_->Branch("MuonCand_vtxx", MuonCand_vtxx, "MuonCand_vtxx[nMuonCand]/D");
    tree_->Branch("MuonCand_vtxy", MuonCand_vtxy, "MuonCand_vtxy[nMuonCand]/D");
    tree_->Branch("MuonCand_vtxz", MuonCand_vtxz, "MuonCand_vtxz[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackVtxz", MuonCand_innerTrackVtxz, "MuonCand_innerTrackVtxz[nMuonCand]/D");
    tree_->Branch("MuonCand_dxy", MuonCand_dxy, "MuonCand_dxy[nMuonCand]/D");
    tree_->Branch("MuonCand_nstatseg", MuonCand_nstatseg, "MuonCand_nstatseg[nMuonCand]/I");
    tree_->Branch("MuonCand_ntrklayers", MuonCand_ntrklayers, "MuonCand_ntrklayers[nMuonCand]/I");
    tree_->Branch("MuonCand_npxlhits", MuonCand_npxlhits, "MuonCand_npxlhits[nMuonCand]/I");
    tree_->Branch("MuonCand_isglobal", MuonCand_isglobal, "MuonCand_isglobal[nMuonCand]/I");
    tree_->Branch("MuonCand_istracker", MuonCand_istracker, "MuonCand_istracker[nMuonCand]/I");
    tree_->Branch("MuonCand_isstandalone", MuonCand_isstandalone, "MuonCand_isstandalone[nMuonCand]/I");
    tree_->Branch("MuonCand_ispfmuon", MuonCand_ispfmuon, "MuonCand_ispfmuon[nMuonCand]/I");
    tree_->Branch("MuonCand_istight", MuonCand_istight, "MuonCand_istight[nMuonCand]/I");
    tree_->Branch("MuonCandTrack_nmuchits", MuonCandTrack_nmuchits, "MuonCandTrack_nmuchits[nMuonCand]/I");
    tree_->Branch("MuonCandTrack_chisq", MuonCandTrack_chisq, "MuonCandTrack_chisq[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackPt", MuonCand_innerTrackPt, "MuonCand_innerTrackPt[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackEta", MuonCand_innerTrackEta, "MuonCand_innerTrackEta[nMuonCand]/D");
    tree_->Branch("MuonCand_innerTrackPhi", MuonCand_innerTrackPhi, "MuonCand_innerTrackPhi[nMuonCand]/D");
    if (runOnMC_) {
      tree_->Branch("nGenMuonCand", &nGenMuonCand, "nGenMuonCand/I");
      tree_->Branch("nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept, "nGenMuonCandOutOfAccept/I");    
      tree_->Branch("GenMuonCand_pt", GenMuonCand_pt, "GenMuonCand_pt[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_eta", GenMuonCand_eta, "GenMuonCand_eta[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_phi", GenMuonCand_phi, "GenMuonCand_phi[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_e", GenMuonCand_e, "GenMuonCand_e[nGenMuonCand]/D");
    }
  }
  
  if (fetchElectrons_) {
    tree_->Branch("nEleCand", &nEleCand, "nEleCand/I");
    tree_->Branch("EleCand_et", EleCand_et, "EleCand_et[nEleCand]/D");
    tree_->Branch("EleCand_eta", EleCand_eta, "EleCand_eta[nEleCand]/D");
    tree_->Branch("EleCand_phi", EleCand_phi, "EleCand_phi[nEleCand]/D");
    tree_->Branch("EleCand_e", EleCand_e, "EleCand_e[nEleCand]/D");
    tree_->Branch("EleCand_innerTrackPt", EleCand_innerTrackPt, "EleCand_innerTrackPt[nEleCand]/D");
    tree_->Branch("EleCand_innerTrackEta", EleCand_innerTrackEta, "EleCand_innerTrackEta[nEleCand]/D");
    tree_->Branch("EleCand_innerTrackPhi", EleCand_innerTrackPhi, "EleCand_innerTrackPhi[nEleCand]/D");
    tree_->Branch("EleCand_charge", EleCand_charge, "EleCand_charge[nEleCand]/I");
    tree_->Branch("EleCand_vtxx", EleCand_vtxx, "EleCand_vtxx[nEleCand]/D");
    tree_->Branch("EleCand_vtxy", EleCand_vtxy, "EleCand_vtxy[nEleCand]/D");
    tree_->Branch("EleCand_vtxz", EleCand_vtxz, "EleCand_vtxz[nEleCand]/D");
    tree_->Branch("EleCand_innerTrackVtxz", EleCand_innerTrackVtxz, "EleCand_innerTrackVtxz[nEleCand]/D");
    tree_->Branch("EleCand_deltaPhi", EleCand_deltaPhi, "EleCand_deltaPhi[nEleCand]/D"); 
    tree_->Branch("EleCand_deltaEta", EleCand_deltaEta, "EleCand_deltaEta[nEleCand]/D"); 
    tree_->Branch("EleCand_HoverE", EleCand_HoverE, "EleCand_HoverE[nEleCand]/D"); 
    tree_->Branch("EleCand_trackiso", EleCand_trackiso, "EleCand_trackiso[nEleCand]/D"); 
    tree_->Branch("EleCand_ecaliso", EleCand_ecaliso," EleCand_ecaliso[nEleCand]/D"); 
    tree_->Branch("EleCand_hcaliso", EleCand_hcaliso," EleCand_hcaliso[nEleCand]/D"); 
    tree_->Branch("EleCand_sigmaIetaIeta", EleCand_sigmaIetaIeta, "EleCand_sigmaIetaIeta[nEleCand]/D"); 
    tree_->Branch("EleCand_convDist", EleCand_convDist, "EleCand_convDist[nEleCand]/D"); 
    tree_->Branch("EleCand_convDcot", EleCand_convDcot, "EleCand_convDcot[nEleCand]/D"); 
    tree_->Branch("EleCand_ecalDriven", EleCand_ecalDriven, "EleCand_ecalDriven[nEleCand]/D");  
    tree_->Branch("EleCand_mediumID", EleCand_mediumID, "EleCand_mediumID[nEleCand]/I");    
    tree_->Branch("EleCand_looseID", EleCand_looseID, "EleCand_looseID[nEleCand]/I");   
    tree_->Branch("EleCand_tightID", EleCand_tightID, "EleCand_tightID[nEleCand]/I");   
    tree_->Branch("EleCand_innerTrackPt", EleCand_innerTrackPt, "EleCand_innerTrackPt[nEleCand]/D");
    tree_->Branch("EleCand_innerTrackEta", EleCand_innerTrackEta, "EleCand_innerTrackEta[nEleCand]/D");
    tree_->Branch("EleCand_innerTrackPhi", EleCand_innerTrackPhi, "EleCand_innerTrackPhi[nEleCand]/D");
    if (runOnMC_) {
      tree_->Branch("nGenEleCand", &nGenEleCand, "nGenEleCand/I");
      tree_->Branch("nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept, "nGenEleCandOutOfAccept/I");    
      tree_->Branch("GenEleCand_pt", GenEleCand_pt, "GenEleCand_pt[nGenEleCand]/D");
      tree_->Branch("GenEleCand_eta", GenEleCand_eta, "GenEleCand_eta[nGenEleCand]/D");
      tree_->Branch("GenEleCand_phi", GenEleCand_phi, "GenEleCand_phi[nGenEleCand]/D");
      tree_->Branch("GenEleCand_e", GenEleCand_e, "GenEleCand_e[nGenEleCand]/D");
    }
  }
  tree_->Branch("nPhotonCand", &nPhotonCand, "nPhotonCand/I");
  tree_->Branch("PhotonCand_pt", PhotonCand_pt, "PhotonCand_pt[nPhotonCand]/D");
  tree_->Branch("PhotonCand_eta", PhotonCand_eta, "PhotonCand_eta[nPhotonCand]/D");
  tree_->Branch("PhotonCand_phi", PhotonCand_phi, "PhotonCand_phi[nPhotonCand]/D");
  tree_->Branch("PhotonCand_e", PhotonCand_e, "PhotonCand_e[nPhotonCand]/D");
  tree_->Branch("PhotonCand_r9", PhotonCand_r9, "PhotonCand_r9[nPhotonCand]/D");
  tree_->Branch("PhotonCand_drtrue", PhotonCand_drtrue, "PhotonCand_drtrue[nPhotonCand]/D"); 
  tree_->Branch("PhotonCand_detatrue", PhotonCand_detatrue, "PhotonCand_detatrue[nPhotonCand]/D"); 
  tree_->Branch("PhotonCand_dphitrue", PhotonCand_dphitrue, "PhotonCand_dphitrue[nPhotonCand]/D"); 
  if (runOnMC_) {
    tree_->Branch("nGenPhotCand", &nGenPhotCand, "nGenPhotCand/I");    
    tree_->Branch("nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept, "nGenPhotCandOutOfAccept/I");    
    tree_->Branch("GenPhotCand_pt", GenPhotCand_pt, "GenPhotCand_pt[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_eta", GenPhotCand_eta, "GenPhotCand_eta[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_phi", GenPhotCand_phi, "GenPhotCand_phi[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_e", GenPhotCand_e, "GenPhotCand_e[nGenPhotCand]/D");
    tree_->Branch("nGenProtCand", &nGenProtCand, "nGenProtCand/I");    
    tree_->Branch("GenProtCand_pt", GenProtCand_pt, "GenProtCand_pt[nGenProtCand]/D");
    tree_->Branch("GenProtCand_eta", GenProtCand_eta, "GenProtCand_eta[nGenProtCand]/D");
    tree_->Branch("GenProtCand_phi", GenProtCand_phi, "GenProtCand_phi[nGenProtCand]/D");
    tree_->Branch("GenProtCand_e", GenProtCand_e, "GenProtCand_e[nGenProtCand]/D");
    tree_->Branch("GenProtCand_status", GenProtCand_status, "GenProtCand_status[nGenProtCand]/I");
  }
  
  // Primary vertices' information
  tree_->Branch("nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/I");
  tree_->Branch("nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nFilteredPrimVertexCand/I");
  tree_->Branch("PrimVertexCand_id", PrimVertexCand_id, "PrimVertexCand_id[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/I");

  // Lepton pairs' information
  if (useLegacyVertexing_) {
    tree_->Branch("PrimVertexCand_matchedtracks", PrimVertexCand_matchedtracks, "PrimVertexCand_matchedtracks[nPrimVertexCand]/I");
    tree_->Branch("PrimVertexCand_unmatchedtracks", PrimVertexCand_unmatchedtracks, "PrimVertexCand_unmatchedtracks[nPrimVertexCand]/I");
    tree_->Branch("PrimVertexCand_hasdil", PrimVertexCand_hasdil, "PrimVertexCand_hasdil[nPrimVertexCand]/I");
    tree_->Branch("Pair_lepton1", Pair_lepton1, "Pair_lepton1[nPrimVertexCand]/I");
    tree_->Branch("Pair_lepton2", Pair_lepton2, "Pair_lepton2[nPrimVertexCand]/I");
    tree_->Branch("Pair_mindist", Pair_mindist, "Pair_mindist[nPrimVertexCand]/D");
    tree_->Branch("Pair_mass", Pair_mass, "Pair_mass[nPrimVertexCand]/D");
    tree_->Branch("Pair_pt", Pair_pt, "Pair_pt[nPrimVertexCand]/D");
    tree_->Branch("Pair_eta", Pair_eta, "Pair_eta[nPrimVertexCand]/D");
    tree_->Branch("Pair_phi", Pair_phi, "Pair_phi[nPrimVertexCand]/D");
    tree_->Branch("Pair_p", Pair_p, "Pair_p[nPrimVertexCand]/D");
    tree_->Branch("Pair_dpt", Pair_dpt, "Pair_dpt[nPrimVertexCand]/D");
    tree_->Branch("Pair_dphi", Pair_dphi, "Pair_dphi[nPrimVertexCand]/D");
    tree_->Branch("Pair_3Dangle", Pair_3Dangle, "Pair_3Dangle[nPrimVertexCand]/D");
    tree_->Branch("Pair_extratracks1mm", Pair_extratracks1mm, "Pair_extratracks1mm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks2mm", Pair_extratracks2mm, "Pair_extratracks2mm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks3mm", Pair_extratracks3mm, "Pair_extratracks3mm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks4mm", Pair_extratracks4mm, "Pair_extratracks4mm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks5mm", Pair_extratracks5mm, "Pair_extratracks5mm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks1cm", Pair_extratracks1cm, "Pair_extratracks1cm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks2cm", Pair_extratracks2cm, "Pair_extratracks2cm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks3cm", Pair_extratracks3cm, "Pair_extratracks3cm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks4cm", Pair_extratracks4cm, "Pair_extratracks4cm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks5cm", Pair_extratracks5cm, "Pair_extratracks5cm[nPrimVertexCand]/I");
    tree_->Branch("Pair_extratracks10cm", Pair_extratracks10cm, "Pair_extratracks10cm[nPrimVertexCand]/I");
    tree_->Branch("PairGamma_mass", PairGamma_mass, "PairGamma_mass[nPrimVertexCand][nPhotonCand]/D");
    // Kalman dilepton vertex information
    tree_->Branch("KalmanVertexCand_x", KalmanVertexCand_x, "KalmanVertexCand_x[nPrimVertexCand]/D");
    tree_->Branch("KalmanVertexCand_y", KalmanVertexCand_y, "KalmanVertexCand_y[nPrimVertexCand]/D");
    tree_->Branch("KalmanVertexCand_z", KalmanVertexCand_z, "KalmanVertexCand_z[nPrimVertexCand]/D");
  }
  else {
    tree_->Branch("nPair", &nPair, "nPair/I");
    tree_->Branch("Pair_lepton1", Pair_lepton1, "Pair_lepton1[nPair]/I");
    tree_->Branch("Pair_lepton2", Pair_lepton2, "Pair_lepton2[nPair]/I");
    tree_->Branch("Pair_mindist", Pair_mindist, "Pair_mindist[nPair]/D");
    tree_->Branch("Pair_mass", Pair_mass, "Pair_mass[nPair]/D");
    tree_->Branch("Pair_pt", Pair_pt, "Pair_pt[nPair]/D");
    tree_->Branch("Pair_eta", Pair_eta, "Pair_eta[nPair]/D");
    tree_->Branch("Pair_phi", Pair_phi, "Pair_phi[nPair]/D");
    tree_->Branch("Pair_p", Pair_p, "Pair_p[nPair]/D");
    tree_->Branch("Pair_dpt", Pair_dpt, "Pair_dpt[nPair]/D");
    tree_->Branch("Pair_dphi", Pair_dphi, "Pair_dphi[nPair]/D");
    tree_->Branch("Pair_3Dangle", Pair_3Dangle, "Pair_3Dangle[nPair]/D");
    tree_->Branch("Pair_extratracks1mm", Pair_extratracks1mm, "Pair_extratracks1mm[nPair]/I");
    tree_->Branch("Pair_extratracks2mm", Pair_extratracks2mm, "Pair_extratracks2mm[nPair]/I");
    tree_->Branch("Pair_extratracks3mm", Pair_extratracks3mm, "Pair_extratracks3mm[nPair]/I");
    tree_->Branch("Pair_extratracks4mm", Pair_extratracks4mm, "Pair_extratracks4mm[nPair]/I");
    tree_->Branch("Pair_extratracks5mm", Pair_extratracks5mm, "Pair_extratracks5mm[nPair]/I");
    tree_->Branch("Pair_extratracks1cm", Pair_extratracks1cm, "Pair_extratracks1cm[nPair]/I");
    tree_->Branch("Pair_extratracks2cm", Pair_extratracks2cm, "Pair_extratracks2cm[nPair]/I");
    tree_->Branch("Pair_extratracks3cm", Pair_extratracks3cm, "Pair_extratracks3cm[nPair]/I");
    tree_->Branch("Pair_extratracks4cm", Pair_extratracks4cm, "Pair_extratracks4cm[nPair]/I");
    tree_->Branch("Pair_extratracks5cm", Pair_extratracks5cm, "Pair_extratracks5cm[nPair]/I");
    tree_->Branch("Pair_extratracks10cm", Pair_extratracks10cm, "Pair_extratracks10cm[nPair]/I");
    tree_->Branch("PairGamma_mass", PairGamma_mass, "PairGamma_mass[nPair][nPhotonCand]/D");
    // Kalman dilepton vertex information
    tree_->Branch("KalmanVertexCand_x", KalmanVertexCand_x, "KalmanVertexCand_x[nPair]/D");
    tree_->Branch("KalmanVertexCand_y", KalmanVertexCand_y, "KalmanVertexCand_y[nPair]/D");
    tree_->Branch("KalmanVertexCand_z", KalmanVertexCand_z, "KalmanVertexCand_z[nPair]/D");
  }
  if (runOnMC_) {
    tree_->Branch("GenPair_mass", &GenPair_mass, "GenPair_mass/D");
    tree_->Branch("GenPair_pt", &GenPair_pt, "GenPair_pt/D");
    tree_->Branch("GenPair_eta", &GenPair_eta, "GenPair_eta/D");
    tree_->Branch("GenPair_phi", &GenPair_phi, "GenPair_phi/D");
    tree_->Branch("GenPair_dpt", &GenPair_dpt, "GenPair_dpt/D");
    tree_->Branch("GenPair_dphi", &GenPair_dphi, "GenPair_dphi/D");
    tree_->Branch("GenPair_3Dangle", &GenPair_3Dangle, "GenPair_3Dangle/D");
  }

  if (fetchProtons_) {
    tree_->Branch("nLocalProtCand", &nLocalProtCand, "nLocalProtCand/I");
    tree_->Branch("LocalProtCand_x", LocalProtCand_x, "LocalProtCand_x[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_y", LocalProtCand_y, "LocalProtCand_y[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_z", LocalProtCand_z, "LocalProtCand_z[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_xSigma", LocalProtCand_xSigma, "LocalProtCand_xSigma[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_ySigma", LocalProtCand_ySigma, "LocalProtCand_ySigma[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_arm", LocalProtCand_arm, "LocalProtCand_arm[nLocalProtCand]/I");
    tree_->Branch("LocalProtCand_side", LocalProtCand_side, "LocalProtCand_side[nLocalProtCand]/I");
    tree_->Branch("LocalProtCand_Tx", LocalProtCand_Tx, "LocalProtCand_Tx[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_Ty", LocalProtCand_Ty, "LocalProtCand_Ty[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_TxSigma", LocalProtCand_TxSigma, "LocalProtCand_TxSigma[nLocalProtCand]/D");
    tree_->Branch("LocalProtCand_TySigma", LocalProtCand_TySigma, "LocalProtCand_TySigma[nLocalProtCand]/D");
  }

  // Extra tracks on vertex's information
  tree_->Branch("nExtraTracks", &nExtraTracks, "nExtraTracks/I");
  tree_->Branch("ExtraTrack_pair", ExtraTrack_pair, "ExtraTrack_pair[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_px", ExtraTrack_px, "ExtraTrack_px[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_py", ExtraTrack_py, "ExtraTrack_py[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_pz", ExtraTrack_pz, "ExtraTrack_pz[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_chi2", ExtraTrack_chi2, "ExtraTrack_chi2[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz, "ExtraTrack_vtxdxyz[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_vtxT", ExtraTrack_vtxT, "ExtraTrack_vtxT[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_vtxZ", ExtraTrack_vtxZ, "ExtraTrack_vtxZ[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_x", ExtraTrack_x, "ExtraTrack_x[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_y", ExtraTrack_y, "ExtraTrack_y[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_z", ExtraTrack_z, "ExtraTrack_z[nExtraTracks]/D");
  tree_->Branch("nQualityExtraTrack", &nQualityExtraTrack, "nQualityExtraTrack/I");
  if (useLegacyVertexing_) {
    tree_->Branch("ClosestExtraTrack_vtxdxyz",ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz[nPrimVertexCand]/D");
    tree_->Branch("ClosestExtraTrack_id",ClosestExtraTrack_id,"ClosestExtraTrack_id[nPrimVertexCand]/I");
    tree_->Branch("ClosestHighPurityExtraTrack_vtxdxyz",ClosestHighPurityExtraTrack_vtxdxyz,"ClosestHighPurityExtraTrack_vtxdxyz[nPrimVertexCand]/D");
    tree_->Branch("ClosestHighPurityExtraTrack_id",ClosestHighPurityExtraTrack_id,"ClosestHighPurityExtraTrack_id[nPrimVertexCand]/I");
  }
  else {
    tree_->Branch("ClosestExtraTrack_vtxdxyz",ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz[nPair]/D");
    tree_->Branch("ClosestExtraTrack_id",ClosestExtraTrack_id,"ClosestExtraTrack_id[nPair]/I");
    tree_->Branch("ClosestHighPurityExtraTrack_vtxdxyz",ClosestHighPurityExtraTrack_vtxdxyz,"ClosestHighPurityExtraTrack_vtxdxyz[nPair]/D");
    tree_->Branch("ClosestHighPurityExtraTrack_id",ClosestHighPurityExtraTrack_id,"ClosestHighPurityExtraTrack_id[nPair]/I");
  }
  
  // Jets/MET information
  tree_->Branch("nJetCand", &nJetCand, "nJetCand/I");
  tree_->Branch("JetCand_pt", JetCand_pt, "JetCand_pt[nJetCand]/D");
  tree_->Branch("JetCand_eta", JetCand_eta, "JetCand_eta[nJetCand]/D");
  tree_->Branch("JetCand_phi", JetCand_phi, "JetCand_phi[nJetCand]/D");
  tree_->Branch("JetCand_e", JetCand_e, "JetCand_e[nJetCand]/D");
  tree_->Branch("HighestJet_pt", &HighestJet_pt, "HighestJet_pt/D");
  tree_->Branch("HighestJet_eta", &HighestJet_eta, "HighestJet_eta/D");
  tree_->Branch("HighestJet_phi", &HighestJet_phi, "HighestJet_phi/D");
  tree_->Branch("HighestJet_e", &HighestJet_e, "HighestJet_e/D");
  tree_->Branch("SumJet_e", &SumJet_e, "SumJet_e/D");
  tree_->Branch("Etmiss", &Etmiss, "Etmiss/D");
  tree_->Branch("Etmiss_phi", &Etmiss_phi, "Etmiss_phi/D");
  tree_->Branch("Etmiss_significance", &Etmiss_significance, "Etmiss_significance/D");

  // Pileup reweighting
  tree_->Branch("Weight", &Weight, "Weight/D");
  tree_->Branch("PUWeightTrue", &PUWeightTrue, "PUWeightTrue/D");

  nCandidates = 0;
}

void
GammaGammaLL::clearTree()
{
  nCandidatesInEvent = 0;
  nPair = 0;
  nPrimVertexCand = nFilteredPrimVertexCand = -1;
  nMuonCand = nEleCand = nLeptonCand = 0;
  nExtraTracks = nQualityExtraTrack = 0;
  nJetCand = 0;
  nGenMuonCand = nGenMuonCandOutOfAccept = 0;
  nGenEleCand = nGenEleCandOutOfAccept = 0;
  nGenPhotCand = nGenPhotCandOutOfAccept = 0;
  nGenProtCand = 0;
  nPhotonCand = 0;

  //LHCFillNum = LHCBeamMode = -1;

  HPS_acc420b1 = HPS_acc220b1 = HPS_acc420and220b1 = HPS_acc420or220b1 = -1;
  HPS_acc420b2 = HPS_acc220b2 = HPS_acc420and220b2 = HPS_acc420or220b2 = -1;
  GenPair_pt = GenPair_mass = GenPair_phi = GenPair_eta = -999.;
  GenPair_dphi = GenPair_dpt = GenPair_3Dangle = 0.;

  for (unsigned int i=0; i<MAX_LL; i++) {
    MuonCand_pt[i] = MuonCand_eta[i] = MuonCand_phi[i] = MuonCand_e[i] = -999.;
    MuonCand_charge[i] = -999;
    MuonCand_vtxx[i] = MuonCand_vtxy[i] = MuonCand_vtxz[i] = -999.;
    MuonCand_innerTrackPt[i] = MuonCand_innerTrackEta[i] = MuonCand_innerTrackPhi[i] = -999.;
    MuonCand_innerTrackVtxz[i] = -999.;
    MuonCand_npxlhits[i] = MuonCand_nstatseg[i] = MuonCand_ntrklayers[i] = -999;
    MuonCand_dxy[i] = -999.;
    MuonCand_isglobal[i] = MuonCand_istracker[i] = MuonCand_isstandalone[i] = MuonCand_ispfmuon[i] = -999;
    MuonCand_istight[i] = -999;
    MuonCandTrack_nmuchits[i] = -999;
    MuonCandTrack_chisq[i] = -999.;
    EleCand_e[i] = EleCand_et[i] = EleCand_phi[i] = EleCand_eta[i] = -999.;
    EleCand_charge[i] = -999;
    EleCand_vtxx[i] = EleCand_vtxy[i] = EleCand_vtxz[i] = -999.;
    EleCand_innerTrackPt[i] = EleCand_innerTrackEta[i] = EleCand_innerTrackPhi[i] = -999.;
    EleCand_innerTrackVtxz[i] = -999.;
    EleCand_deltaPhi[i] = EleCand_deltaEta[i] = EleCand_HoverE[i] = -999.;
    EleCand_trackiso[i] = EleCand_ecaliso[i] = EleCand_hcaliso[i] = EleCand_sigmaIetaIeta[i] = -999.;
    EleCand_convDist[i] = EleCand_convDcot[i] = EleCand_ecalDriven[i] = -999.;
    EleCand_tightID[i] = EleCand_mediumID[i] = EleCand_looseID[i] = -1;
  }
  for (unsigned int i=0; i<MAX_ET && i<maxExTrkVtx_; i++) {
    ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
    ExtraTrack_charge[i] = ExtraTrack_ndof[i] = -999;
    ExtraTrack_chi2[i] = ExtraTrack_vtxdxyz[i] = -999.;
    ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
    ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
    ExtraTrack_pair[i] = -1;
  }
  for (unsigned int i=0; i<MAX_PAIRS; i++) {
    Pair_lepton1[i] = Pair_lepton2[i] = -1;
    Pair_mindist[i] = Pair_p[i] = Pair_pt[i] = Pair_mass[i] = Pair_phi[i] = Pair_eta[i] = -999.;
    Pair_dphi[i] = Pair_dpt[i] = Pair_3Dangle[i] = -999.;
    Pair_extratracks1mm[i] = Pair_extratracks2mm[i] = Pair_extratracks3mm[i] = 0;
    Pair_extratracks4mm[i] = Pair_extratracks5mm[i] = Pair_extratracks1cm[i] = 0;
    Pair_extratracks2cm[i] = Pair_extratracks3cm[i] = Pair_extratracks4cm[i] = 0;
    Pair_extratracks5cm[i] = Pair_extratracks10cm[i] = 0;    
  }
  for (unsigned int i=0; i<MAX_VTX; i++) {
    PrimVertexCand_id[i] = -1;
    PrimVertexCand_tracks[i] = PrimVertexCand_matchedtracks[i] = PrimVertexCand_unmatchedtracks[i] = PrimVertexCand_hasdil[i] = 0;
    PrimVertexCand_x[i] = PrimVertexCand_y[i] = PrimVertexCand_z[i] = -999.;
    PrimVertexCand_chi2[i] = PrimVertexCand_ndof[i] = -999.;
    KalmanVertexCand_x[i] = KalmanVertexCand_y[i] = KalmanVertexCand_z[i] = -999.;
    ClosestExtraTrackKalman_vtxdxyz[i] = 999.;
  }
  for (unsigned int i=0; i<MAX_PHO; i++) { 
    PhotonCand_e[i] = PhotonCand_pt[i] = PhotonCand_eta[i] = PhotonCand_phi[i] = PhotonCand_r9[i] = -999.;
    PhotonCand_detatrue[i] = PhotonCand_dphitrue[i] = PhotonCand_drtrue[i] = -999.;
  }
  for (unsigned int i=0; i<MAX_JETS; i++) {
    JetCand_pt[i] = JetCand_phi[i] = JetCand_eta[i] = JetCand_e[i] = -999;
  }
  HighestJet_pt = HighestJet_eta = HighestJet_phi = HighestJet_e = -999.;
  SumJet_e = 0.;
  Etmiss = Etmiss_phi = Etmiss_significance = -999.;
  for (unsigned int i=0; i<MAX_LOCALPCAND; i++) {
    LocalProtCand_x[i] = LocalProtCand_y[i] = LocalProtCand_z[i] = -999.;
    LocalProtCand_xSigma[i] = LocalProtCand_ySigma[i] = -999.;
    LocalProtCand_Tx[i] = LocalProtCand_Ty[i] = -999.;
    LocalProtCand_TxSigma[i] = LocalProtCand_TySigma[i] = -999.;
    LocalProtCand_arm[i] = LocalProtCand_side[i] = -1;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaLL::endJob() 
{
  std::cout << "==> Number of candidates in the dataset : " << nCandidates << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
GammaGammaLL::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed = true;
  if (!hltPrescale_.init(iRun, iSetup, hltMenuLabel_, changed)) {
    throw cms::Exception("GammaGammaLL") << " prescales extraction failure with process name " << hltMenuLabel_;
  }
  // Initialise HLTConfigProvider
  hltConfig_ = hltPrescale_.hltConfigProvider();
  if (!hltConfig_.init(iRun, iSetup, hltMenuLabel_, changed)) {
    throw cms::Exception("GammaGammaLL") << " config extraction failure with process name " << hltMenuLabel_;
  }
  else if (hltConfig_.size()<=0) {
    edm::LogError("GammaGammaLL") << "HLT config size error";
  }
}

// ------------ method called when ending the processing of a run  ------------
void 
GammaGammaLL::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GammaGammaLL::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GammaGammaLL::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GammaGammaLL::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaGammaLL);
