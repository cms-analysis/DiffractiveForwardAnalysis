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
  fetchMuons_(false), fetchElectrons_(false),
  fetchProtons_       (iConfig.getParameter<bool>("fetchProtons")),
  verb_               (iConfig.getUntrackedParameter<unsigned int>("verbosity", 1)),
  outputFile_         (iConfig.getUntrackedParameter<std::string>("outfilename", "output.root")),
  leptonsType_        (iConfig.getParameter< std::vector<std::string> >("LeptonsType")),
  hltMenuLabel_       (iConfig.getParameter<std::string>("HLTMenuLabel")),
  triggersList_       (iConfig.getParameter<std::vector<std::string> >("TriggersList")),
  triggerResultsToken_(consumes<edm::TriggerResults>              (iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  //beamSpotToken_      (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotLabel"))),
  recoVertexToken_    (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("RecoVertexLabel"))),
  muonToken_          (consumes< edm::View<pat::Muon> >           (iConfig.getUntrackedParameter<edm::InputTag>("GlobalMuonCollectionLabel", std::string("muons")))),
  eleToken_           (consumes< edm::View<pat::Electron> >       (iConfig.getUntrackedParameter<edm::InputTag>("GlobalEleCollectionLabel", std::string("gsfElectrons")))),
  eleLooseIdMapToken_ (consumes< edm::ValueMap<bool> >            (iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes< edm::ValueMap<bool> >            (iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_ (consumes< edm::ValueMap<bool> >            (iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  //conversionsToken_   (consumes<reco::ConversionCollection>       (iConfig.getParameter<edm::InputTag>("conversionsLabel"))),
  pileupToken_        (consumes< std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileupInfo"))),
  jetToken_           (consumes< edm::View<pat::Jet> >            (iConfig.getParameter<edm::InputTag>("JetCollectionLabel"))),
  metToken_           (consumes< edm::View<pat::MET> >            (iConfig.getParameter<edm::InputTag>("MetLabel"))),
  totemRPHitToken_    (consumes< edm::DetSetVector<TotemRPLocalTrack> >(iConfig.getParameter<edm::InputTag>("totemRPLocalTrackLabel"))),
  runOnMC_            (iConfig.getUntrackedParameter<bool>("RunOnMC", false)),
  printCandidates_    (iConfig.getUntrackedParameter<bool>("PrintCandidates", false)),
  sqrts_              (iConfig.getParameter<double>("SqrtS")),
  maxExTrkVtx_        (iConfig.getUntrackedParameter<unsigned int>("maxExtraTracks", 1000)),
  hltPrescale_        (iConfig, consumesCollector(), *this),
  lumiWeights_(0),
  photonToken_        (consumes< edm::View<pat::Photon> >         (iConfig.getParameter<edm::InputTag>("photonLabel"))),
  pflowToken_         (consumes< edm::View<reco::PFCandidate> >   (iConfig.getUntrackedParameter<edm::InputTag>("PFLabel", std::string("particleFlow")))) //FIXME will not work for miniAOD
  //pflowToken_         (consumes< edm::View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("PFLabel", std::string("particleFlow"))))
{
  //now do what ever initialization is needed
  
  hlts_ = new HLTMatcher(triggersList_);
  nHLT = triggersList_.size();	
  
  // Generator level
  if (runOnMC_) {
    genToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticlesCollectionLabel"));
    minPtMC_ = iConfig.getUntrackedParameter<double>("MCAcceptPtCut", 10.);
    minEtaMC_ = iConfig.getUntrackedParameter<double>("MCAcceptEtaCut", 2.5);

    // Pileup input tags
    mcPileupFile_ = iConfig.getUntrackedParameter<std::string>("mcpufile", "PUHistos.root");
    mcPileupPath_ = iConfig.getUntrackedParameter<std::string>("mcpupath", "pileup");
    dataPileupFile_ = iConfig.getUntrackedParameter<std::string>("datapufile", "PUHistos_duplicated.root");
    dataPileupPath_ = iConfig.getUntrackedParameter<std::string>("datapupath", "pileup");  
  }

  // Leptons input tags
  for (std::vector<std::string>::const_iterator it=leptonsType_.begin(); it!=leptonsType_.end(); it++) {
    if (*it=="Muon") fetchMuons_ = true;
    else if (*it=="Electron") fetchElectrons_ = true;
    else {
      throw cms::Exception("GammaGammaLL") << "'LeptonsType' parameter should either be:\n"
                                           << "   * ['Electron'], or ['Muon']  (for same-flavour leptons)\n"
                                           << "   * ['Electron', 'Muon']       (for electron-muon pairs)";
    }
  }

  file_ = new TFile(outputFile_.c_str(), "recreate");
  file_->cd();
  // tree definition
  tree_ = new TTree("ntp1", "ntp1");

  // HPS acceptance file readout definition
  if (runOnMC_) {
    // edm::FileInPath myDataFile("FastSimulation/ProtonTaggers/data/acceptance_420_220.root");  
    /*myDataFile = new edm::FileInPath("FastSimulation/ForwardDetectors/data/acceptance_420_220.root");
    fullAcceptancePath = myDataFile->fullPath();
    std::cout << "Opening " << fullAcceptancePath << std::endl;
    f = new TFile(fullAcceptancePath.c_str());
    if (f->Get("description")!=NULL) {
      std::cout << "Description found: " << f->Get("description")->GetTitle() << std::endl;
      std::cout << "Reading acceptance tables " << std::endl;
    }
    helper420beam1.Init(*f, "a420");
    helper420beam2.Init(*f, "a420_b2");
    helper220beam1.Init(*f, "a220");
    helper220beam2.Init(*f, "a220_b2");
    helper420a220beam1.Init(*f, "a420a220");
    helper420a220beam2.Init(*f, "a420a220_b2");*/

    lumiWeights_ = new edm::LumiReWeighting(mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_);
    //lumiWeights_->setPileupSummaryInfoInputTag(iConfig.getParameter<edm::InputTag>("pileupInfo"));
  }
}

GammaGammaLL::~GammaGammaLL()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  /*const edm::ParameterSet& pset = edm::getProcessParameterSet();
  TList *list = tree_->GetUserInfo();
  list->Add(new TObjString(pset.dump().c_str()));*/
  file_->Write();
  file_->Close();

  if (hlts_) delete hlts_;
  if (tree_) delete tree_;
  if (lumiWeights_) delete lumiWeights_;
}

//
// member functions
//

void
GammaGammaLL::lookAtTriggers(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  int trigNum;

  // Get the trigger information from the event
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken(triggerResultsToken_, hltResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);

  std::ostringstream os;
  if (verb_>2) { os << "Trigger names: " << std::endl; }
  for (unsigned int i=0; i<trigNames.size(); i++) {
    if (verb_>2) os << "--> " << trigNames.triggerNames().at(i) << std::endl;
    if (!hlts_) continue;
    trigNum = hlts_->TriggerNum(trigNames.triggerNames().at(i));
    if (trigNum<0) continue; // Trigger didn't match the interesting ones
    HLT_Accept[trigNum] = hltResults->accept(i);
    //if (runOnMC_) { HLT_Prescl[trigNum] = 1.; continue; } //FIXME
    int prescale_set = hltPrescale_.prescaleSet(iEvent, iSetup);
    HLT_Prescl[trigNum] = (prescale_set<0) ? 0. : hltConfig_.prescaleValue(prescale_set, trigNames.triggerNames().at(i)); //FIXME
  }
  if (verb_>2) edm::LogWarning("GammaGammaLL") << os.str();
}

// ------------ method called for each event  ------------
void
GammaGammaLL::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

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
  
  if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed trigger filtering stage";

  // beam spot information
  /*iEvent.getByToken(beamSpotToken_, beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());*/

  // Get the vertex collection from the event
  iEvent.getByToken(recoVertexToken_, recoVertexColl_);
  const reco::VertexCollection* vertices = recoVertexColl_.product();

  // Generator level information
  if (runOnMC_) {
    iEvent.getByToken(genToken_, genPartColl_);
    
    for (genPart=genPartColl_->begin(); genPart!=genPartColl_->end(); genPart++) {
      if (!genPart->isPromptFinalState()) continue;
      if (genPart->pt()<minPtMC_ || (minEtaMC_!=-1. && fabs(genPart->eta())>minEtaMC_)) {
        if (fabs(genPart->pdgId())==13) nGenMuonCandOutOfAccept++;
        if (fabs(genPart->pdgId())==11) nGenEleCandOutOfAccept++;
        if (fabs(genPart->pdgId())==22) nGenPhotCandOutOfAccept++;
        if (genPart->pdgId()!=2212) continue; // we keep the forward protons
      }
      if (genPart->pdgId()==2212 && nGenProtCand<MAX_GENPRO) {
        GenProtCand_p[nGenProtCand] = genPart->p();
        GenProtCand_px[nGenProtCand] = genPart->px();
        GenProtCand_py[nGenProtCand] = genPart->py();
        GenProtCand_pz[nGenProtCand] = genPart->pz();
        GenProtCand_pt[nGenProtCand] = genPart->pt();
        GenProtCand_eta[nGenProtCand] = genPart->eta();
        GenProtCand_phi[nGenProtCand] = genPart->phi();
	GenProtCand_status[nGenProtCand] = genPart->status();
        nGenProtCand++;
      }
      if (fabs(genPart->pdgId())==13 && nGenMuonCand<MAX_GENMU) {
        GenMuonCand_p[nGenMuonCand] = genPart->p();
        GenMuonCand_px[nGenMuonCand] = genPart->px();
        GenMuonCand_py[nGenMuonCand] = genPart->py();
        GenMuonCand_pz[nGenMuonCand] = genPart->pz();
        GenMuonCand_pt[nGenMuonCand] = genPart->pt();
        GenMuonCand_eta[nGenMuonCand] = genPart->eta();
        GenMuonCand_phi[nGenMuonCand] = genPart->phi();
        nGenMuonCand++;
      }
      if (fabs(genPart->pdgId())==11 && nGenEleCand<MAX_GENELE) {
        GenEleCand_p[nGenEleCand] = genPart->p();
        GenEleCand_px[nGenEleCand] = genPart->px();
        GenEleCand_py[nGenEleCand] = genPart->py();
        GenEleCand_pz[nGenEleCand] = genPart->pz();
        GenEleCand_pt[nGenEleCand] = genPart->pt();
        GenEleCand_eta[nGenEleCand] = genPart->eta();
        GenEleCand_phi[nGenEleCand] = genPart->phi(); 
        nGenEleCand++;
      }
      if (genPart->pdgId()==22 && nGenPhotCand<MAX_GENPHO) {
        GenPhotCand_e[nGenPhotCand] = genPart->energy();
        GenPhotCand_p[nGenPhotCand] = genPart->p();
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
        t = -(std::pow(genPart->pt(), 2)+std::pow(MASS_P*xi, 2))/(1.-xi);
        
        // FIXME to be updated with the new CT-PPS simulation (high-PU)
        // HPS acceptance computation
        /*if (genPart->pz()>0.) {
          HPS_acc420b1 = helper420beam1.GetAcceptance(t, xi, genPart->phi());
          HPS_acc220b1 = helper220beam1.GetAcceptance(t, xi, genPart->phi());
          HPS_acc420and220b1 = helper420a220beam1.GetAcceptance(t, xi, genPart->phi());
          HPS_acc420or220b1 = HPS_acc420b1 + HPS_acc220b1 - HPS_acc420and220b1;
        }
        else {
          HPS_acc420b2 = helper420beam2.GetAcceptance(t, xi, genPart->phi());
          HPS_acc220b2 = helper220beam2.GetAcceptance(t, xi, genPart->phi());
          HPS_acc420and220b2 = helper420a220beam2.GetAcceptance(t, xi, genPart->phi());
          HPS_acc420or220b2 = HPS_acc420b2 + HPS_acc220b2 - HPS_acc420and220b2;
        }*/
      }
    }

    bool foundGenCandPairInEvent = false;
  
    TLorentzVector l1, l2;
    // electron+muon
    if ((fetchElectrons_ && fetchMuons_) and (nGenMuonCand==1 && nGenEleCand==1)) { // FIXME maybe a bit tight according to the newer PU conditions?
      l1.SetXYZM(GenMuonCand_px[0], GenMuonCand_py[0], GenMuonCand_pz[0], MASS_MU);
      l2.SetXYZM(GenEleCand_px[0], GenEleCand_py[0], GenEleCand_pz[0], MASS_E);
      foundGenCandPairInEvent = true;
    }
    // dielectron
    else if (fetchElectrons_ and nGenEleCand==2) { // FIXME maybe a bit tight according to the newer PU conditions?
      l1.SetXYZM(GenEleCand_px[0], GenEleCand_py[0], GenEleCand_pz[0], MASS_E);
      l2.SetXYZM(GenEleCand_px[1], GenEleCand_py[1], GenEleCand_pz[1], MASS_E);
      foundGenCandPairInEvent = true;
    }
    // dimuon
    else if (fetchMuons_ and nGenMuonCand==2) { // FIXME maybe a bit tight according to the newer PU conditions?
      l1.SetXYZM(GenMuonCand_px[0], GenMuonCand_py[0], GenMuonCand_pz[0], MASS_MU);
      l2.SetXYZM(GenMuonCand_px[1], GenMuonCand_py[1], GenMuonCand_pz[1], MASS_MU);
      foundGenCandPairInEvent = true;      	
    }
    if (foundGenCandPairInEvent) {
      pair = l1+l2;
      GenPair_p = pair.P();
      GenPair_pt = pair.Pt();
      GenPair_mass = pair.M();
      GenPair_phi = pair.Phi();
      GenPair_eta = pair.Eta();
      dphi = fabs(l1.Phi()-l2.Phi());
      GenPair_dphi = (dphi<pi) ? dphi : 2.*pi-dphi; // dphi lies in [-pi, pi]
      GenPair_dpt = fabs(l1.Pt()-l2.Pt());
      GenPair_3Dangle = (l1.Angle(l2.Vect()))/pi;
    }

    // Pileup information
    const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
    if (lumiWeights_) Weight = lumiWeights_->weight(*iEventB);

    if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed Pileup retrieval stage";
  } // run on MC?
  
  std::map<int, TLorentzVector> muonsMomenta, electronsMomenta;
  muonsMomenta.clear();
  electronsMomenta.clear();
  
  TLorentzVector leptonptmp;
  
  // Get the muons collection from the event
  if (fetchMuons_) {
    iEvent.getByToken(muonToken_, muonColl_);
    for (muon=muonColl_->begin(); muon!=muonColl_->end() && nMuonCand<MAX_MUONS; muon++) {
      MuonCand_p[nMuonCand] = muon->p();
      MuonCand_px[nMuonCand] = muon->px();
      MuonCand_py[nMuonCand] = muon->py();
      MuonCand_pz[nMuonCand] = muon->pz();
      MuonCand_pt[nMuonCand] = muon->pt();
      MuonCand_eta[nMuonCand] = muon->eta();
      MuonCand_phi[nMuonCand] = muon->phi();
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

      if (MuonCand_istracker[nMuonCand]) {
	MuonCand_npxlhits[nMuonCand] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();
	MuonCand_ntrklayers[nMuonCand] = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
	leptonptmp.SetXYZM(muon->innerTrack()->px(), muon->innerTrack()->py(), muon->innerTrack()->pz(), muon->mass());
      }
      else {
	leptonptmp.SetXYZM(muon->px(), muon->py(), muon->pz(), muon->mass());
      }
      muonsMomenta.insert(std::pair<int,TLorentzVector>(nMuonCand, leptonptmp));

      if (MuonCand_isglobal[nMuonCand] && MuonCand_istracker[nMuonCand]) {
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
      } 

      nMuonCand++;
    }
    nLeptonCand += nMuonCand;
    if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed Muon retrieval stage. Got " << nMuonCand << " muon(s)";
  } // fetch muons?

  // Get the electrons collection from the event
  if (fetchElectrons_) {
    iEvent.getByToken(eleToken_, eleColl_);
    // New 2012 electron ID variables conversions
    //iEvent.getByToken(conversionsToken_, conversions_h);

    edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    iEvent.getByToken(eleLooseIdMapToken_, loose_id_decisions);
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
    iEvent.getByToken(eleTightIdMapToken_, tight_id_decisions);

    for (unsigned int j=0; j<eleColl_->size(); j++) {
      const auto electron = eleColl_->ptrAt(j);

      EleCand_e[nEleCand] = electron->energy();
      EleCand_et[nEleCand] = electron->et();
      EleCand_px[nEleCand] = electron->px();
      EleCand_py[nEleCand] = electron->py();
      EleCand_pz[nEleCand] = electron->pz();
      EleCand_p[nEleCand] = electron->p();
      EleCand_phi[nEleCand] = electron->phi();
      EleCand_eta[nEleCand] = electron->eta();
      EleCand_charge[nEleCand] = electron->charge();
  	  
      EleCand_vtxx[nEleCand] = electron->vertex().x();
      EleCand_vtxy[nEleCand] = electron->vertex().y();
      EleCand_vtxz[nEleCand] = electron->vertex().z();
  	  
      leptonptmp.SetXYZM(electron->px(), electron->py(), electron->pz(), electron->mass());

      if (electron->closestCtfTrackRef().isNonnull()) { // Only for pat::Electron
        EleCandTrack_p[nEleCand] = electron->closestCtfTrackRef()->p(); 
        EleCandTrack_pt[nEleCand] = electron->closestCtfTrackRef()->pt();  
        EleCandTrack_eta[nEleCand] = electron->closestCtfTrackRef()->eta();  
        EleCandTrack_phi[nEleCand] = electron->closestCtfTrackRef()->phi();  
        EleCandTrack_vtxz[nEleCand] = electron->closestCtfTrackRef()->vertex().z();
	leptonptmp.SetPtEtaPhiM(EleCandTrack_pt[nEleCand], EleCandTrack_eta[nEleCand], EleCandTrack_phi[nEleCand], electron->mass());
      }

      electronsMomenta.insert(std::pair<int,TLorentzVector>(nEleCand, leptonptmp));

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
    if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed Electron retrieval stage. Got " << nEleCand << " electron(s)";
  } // fetch electrons?  
  
  // Get the PFlow collection from the event
  //iEvent.getByToken(pflowToken_, pflowColl_);

  // Get the photons collection from the event
  iEvent.getByToken(photonToken_, photonColl_);
  for (photon=photonColl_->begin(); photon!=photonColl_->end(); photon++) { 
    PhotonCand_p[nPhotonCand] = photon->p();
    PhotonCand_px[nPhotonCand] = photon->px();
    PhotonCand_py[nPhotonCand] = photon->py();
    PhotonCand_pz[nPhotonCand] = photon->pz();
    PhotonCand_pt[nPhotonCand] = photon->pt();
    PhotonCand_eta[nPhotonCand] = photon->eta();
    PhotonCand_phi[nPhotonCand] = photon->phi();
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
    if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed photons retrieval stage. Got " << nPhotonCand << " photon(s)";
  }

  // Vertex tracks multiplicity
  vtxind = 0;
  nPrimVertexCand = vertices->size();

  bool foundPairInEvent = false;
  
  if (nLeptonCand>=2) {
    // Enough leptons candidates to go deeper and analyze the primary vertices 
    etind = 0;

    double closesttrkdxyz = 999., closesthighpuritytrkdxyz = 999.;
    for (vertex=vertices->begin(); vertex!=vertices->end() && vtxind<MAX_VTX; ++vertex) {
      PrimaryVertexSelector vtx(leptonsType_, muonsMomenta, electronsMomenta);

      nLeptonsInPrimVertex = 0;
      nExtraTracks = 0;
      nQualityExtraTrack = 0;
      
      vtx.SetPosition(vertex->x(), vertex->y(), vertex->z());

      PrimVertexCand_id[vtxind] = vtxind;
      PrimVertexCand_x[vtxind] = vtx.Position.X();
      PrimVertexCand_y[vtxind] = vtx.Position.Y();
      PrimVertexCand_z[vtxind] = vtx.Position.Z();
      PrimVertexCand_chi2[vtxind] = vertex->chi2();
      PrimVertexCand_ndof[vtxind] = vertex->ndof();
      
      double closesttrkid = -1.;
      // Loop on all the tracks matched with this vertex
      for (reco::Vertex::trackRef_iterator track=vertex->tracks_begin();
          track!=vertex->tracks_end() && etind<MAX_ET && etind<(int)maxExTrkVtx_; track++) {
	const int leptonId_ = vtx.AddTrack((*track).castTo<reco::TrackRef>());
	if (leptonId_==-1) { // Track was not matched to any of the leptons in the collection
          ExtraTrack_vtxId[etind] = vtxind;
	  const double vtxdst = sqrt(std::pow(((*track)->vertex().x()-vtx.Position.X()),2)+
                                std::pow(((*track)->vertex().y()-vtx.Position.Y()),2)+
                                std::pow(((*track)->vertex().z()-vtx.Position.Z()),2));
	  
	  ExtraTrack_purity[etind] = (*track)->quality(reco::TrackBase::highPurity);
	  ExtraTrack_nhits[etind] = (*track)->numberOfValidHits();
	  
	  ExtraTrack_p[etind] = (*track)->p();
	  ExtraTrack_px[etind] = (*track)->px();
	  ExtraTrack_py[etind] = (*track)->py();
	  ExtraTrack_pz[etind] = (*track)->pz();
	  ExtraTrack_pt[etind] = (*track)->pt();
	  ExtraTrack_eta[etind] = (*track)->eta();
	  ExtraTrack_phi[etind] = (*track)->phi();
	  ExtraTrack_charge[etind] = (*track)->charge();
	  ExtraTrack_chi2[etind] = (*track)->chi2();
	  ExtraTrack_ndof[etind] = (*track)->ndof();
          ExtraTrack_vtxdxyz[etind] = vtxdst;
          ExtraTrack_vtxT[etind] = sqrt(std::pow((*track)->vertex().x()-vtx.Position.X(),2)+
					std::pow((*track)->vertex().y()-vtx.Position.Y(),2));
          ExtraTrack_vtxZ[etind] = fabs((*track)->vertex().z()-vtx.Position.Z());
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
          if (vtxdst<closesttrkdxyz) {
            closesttrkdxyz = vtxdst;
            closesttrkid = etind;
          }
          
	  if (ExtraTrack_purity[etind]==1 && ExtraTrack_nhits[etind]>=3) {
	    nQualityExtraTrack++;
	    if (vtxdst<closesthighpuritytrkdxyz) {
	      closesthighpuritytrkdxyz = vtxdst;
              //closesthighpuritytrkid = etind;
	    }
	  }
          etind++;
        }
        else {
          if (verb_>1) { edm::LogWarning("GammaGammaLL") << "-- Lepton match on vertex " << vtxind << " --> " << leptonId_; }
          nLeptonsInPrimVertex++;
        }
      }
      ClosestExtraTrack_vtxdxyz[vtxind] = closesttrkdxyz;
      ClosestExtraTrack_id[vtxind] = closesttrkid;
      ClosestHighPurityExtraTrack_vtxdxyz[vtxind] = closesthighpuritytrkdxyz;
      ClosestHighPurityExtraTrack_id[vtxind] = closesttrkid;
      
      PrimVertexCand_tracks[vtxind] = vtx.nTracks;
      PrimVertexCand_matchedtracks[vtxind] = vtx.nMatchedTracks;
      PrimVertexCand_unmatchedtracks[vtxind] = vtx.nUnmatchedTracks;

      if (nLeptonsInPrimVertex<2) continue;
      
      // At this stage we have at least two matched leptons track on the vertex
      Pair_candidates[vtxind][0] = Pair_lepton1[vtxind] = -1;
      Pair_candidates[vtxind][1] = Pair_lepton2[vtxind] = -1;
      
      if (PrimVertexCand_unmatchedtracks[vtxind]>(int)maxExTrkVtx_) continue; // cut on the upper number of extra tracks

      TLorentzVector l1, l2;
      bool foundPairOnVertex = false;
      double minDist = 999.;

      if (fetchElectrons_ && fetchMuons_) { // Looks at electron+muon
	// Not enough muons or electrons candidates on the vertex
	if (vtx.Electrons()==0 or vtx.Muons()==0) {
          if (verb_>2) { edm::LogWarning("GammaGammaLL") << "Not enough electrons (" << vtx.Electrons() << ") or muons (" << vtx.Muons() << ") arising from the primary vertex !"; }
	  continue;
	}
	for (unsigned int i=0; i<vtx.MatchedMuons.size(); i++) {
	  const int lep1 = vtx.MatchedMuons[i];

	  TVector3 vtxlep1(MuonCand_vtxx[lep1], MuonCand_vtxy[lep1], MuonCand_vtxz[lep1]);
	  MuonCand_dz[lep1] = vtx.dZ(vtxlep1, lep1);
	  MuonCand_istight[lep1]&= (MuonCand_dz[lep1]<.5);

	  for (unsigned int j=0; j<vtx.MatchedElectrons.size(); j++) {
            const int lep2 = vtx.MatchedElectrons[j];
            if (MuonCand_charge[lep1]*EleCand_charge[lep2]>0) continue;
            foundPairOnVertex = true;
            const double leptonsDist = sqrt(pow(MuonCand_vtxx[lep1]-EleCand_vtxx[lep2],2)+
			                    pow(MuonCand_vtxy[lep1]-EleCand_vtxy[lep2],2)+
                                            pow(MuonCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
            if (leptonsDist<minDist) {
              minDist = leptonsDist;
              Pair_candidates[vtxind][0] = Pair_lepton1[vtxind] = lep1;
              Pair_candidates[vtxind][1] = Pair_lepton2[vtxind] = lep2;
            }
          }
        }
        l1.SetXYZM(MuonCand_px[Pair_lepton1[vtxind]],
                   MuonCand_py[Pair_lepton1[vtxind]],
                   MuonCand_pz[Pair_lepton1[vtxind]],
                   MASS_MU);
        l2.SetXYZM(EleCand_px[Pair_lepton2[vtxind]],
                   EleCand_py[Pair_lepton2[vtxind]],
                   EleCand_pz[Pair_lepton2[vtxind]],
                   MASS_E);
      }
      else if (fetchElectrons_) { // Looks at dielectrons
	// Not enough electrons candidates on the vertex
	if (vtx.MatchedElectrons.size()<2) continue;
        for (unsigned int i=0; i<vtx.MatchedElectrons.size(); i++) {
          const int lep1 = vtx.MatchedElectrons[i];
          for (unsigned int j=i+1; j<vtx.MatchedElectrons.size(); j++) {
            const int lep2 = vtx.MatchedElectrons[j];
            if (EleCand_charge[lep1]*EleCand_charge[lep2]>0) continue;
            foundPairOnVertex = true;
            const double leptonsDist = sqrt(pow(EleCand_vtxx[lep1]-EleCand_vtxx[lep2],2)+
                                            pow(EleCand_vtxy[lep1]-EleCand_vtxy[lep2],2)+
                                            pow(EleCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
            if (leptonsDist<minDist) {
              minDist = leptonsDist;
              Pair_candidates[vtxind][0] = Pair_lepton1[vtxind] = lep1;
              Pair_candidates[vtxind][1] = Pair_lepton2[vtxind] = lep2;
            }
          }
        }
        l1.SetXYZM(EleCand_px[Pair_lepton1[vtxind]],
                   EleCand_py[Pair_lepton1[vtxind]],
                   EleCand_pz[Pair_lepton1[vtxind]],
                   MASS_E);
        l2.SetXYZM(EleCand_px[Pair_lepton2[vtxind]],
                   EleCand_py[Pair_lepton2[vtxind]],
                   EleCand_pz[Pair_lepton2[vtxind]],
                   MASS_E);
      }
      else if (fetchMuons_) { // Looks at dimuons
	// Not enough muons candidates on the vertex
	if (vtx.MatchedMuons.size()<2) continue;
        for (unsigned int i=0; i<vtx.MatchedMuons.size(); i++) {
          const int lep1 = vtx.MatchedMuons[i];

	  TVector3 vtxlep1(MuonCand_vtxx[lep1], MuonCand_vtxy[lep1], MuonCand_vtxz[lep1]);
	  MuonCand_dz[lep1] = vtx.dZ(vtxlep1, lep1);
	  MuonCand_istight[lep1]&= (MuonCand_dz[lep1]<.5);

          for (unsigned int j=i+1; j<vtx.MatchedMuons.size(); j++) {
            const int lep2 = vtx.MatchedMuons[j];

	    TVector3 vtxlep2(MuonCand_vtxx[lep2], MuonCand_vtxy[lep2], MuonCand_vtxz[lep2]);
	    MuonCand_dz[lep2] = vtx.dZ(vtxlep2, lep2);
	    MuonCand_istight[lep2]&= (MuonCand_dz[lep2]<.5);

            if (MuonCand_charge[lep1]*MuonCand_charge[lep2]>0) continue;
            foundPairOnVertex = true;
            const double leptonsDist = sqrt(pow(MuonCand_vtxx[lep1]-MuonCand_vtxx[lep2],2)+
                                            pow(MuonCand_vtxy[lep1]-MuonCand_vtxy[lep2],2)+
                                            pow(MuonCand_vtxz[lep1]-MuonCand_vtxz[lep2],2));
            if (leptonsDist<minDist) {
              minDist = leptonsDist;
              Pair_candidates[vtxind][0] = Pair_lepton1[vtxind] = lep1;
              Pair_candidates[vtxind][1] = Pair_lepton2[vtxind] = lep2;
            }
          }
        }
        l1.SetXYZM(MuonCand_px[Pair_lepton1[vtxind]],
                   MuonCand_py[Pair_lepton1[vtxind]],
                   MuonCand_pz[Pair_lepton1[vtxind]],
                   MASS_MU);
        l2.SetXYZM(MuonCand_px[Pair_lepton2[vtxind]],
                   MuonCand_py[Pair_lepton2[vtxind]],
                   MuonCand_pz[Pair_lepton2[vtxind]],
                   MASS_MU);
      }

      if (foundPairOnVertex) {
	Pair_mindist[vtxind] = minDist;
        if (verb_>1) {
      	  std::ostringstream os; os << "---> Found a lepton pair on vertex!" << std::endl;
          os << "Matched muons : " << std::endl;
          for (unsigned int i=0; i<vtx.MatchedMuons.size(); i++) { os << "-> " << vtx.MatchedMuons[i] << std::endl; }
          os << "Matched electrons : " << std::endl;
          for (unsigned int i=0; i<vtx.MatchedElectrons.size(); i++) { os << "-> " << vtx.MatchedElectrons[i] << std::endl; }
          edm::LogInfo("GammaGammaLL") << os.str();
        }
        pair = l1+l2;
        Pair_p[vtxind] = pair.P();
        Pair_pt[vtxind] = pair.Pt();
        Pair_mass[vtxind] = pair.M();
        Pair_phi[vtxind] = pair.Phi();
        Pair_eta[vtxind] = pair.Eta();
        dphi = fabs(l1.Phi()-l2.Phi());
        Pair_dphi[vtxind] = (dphi<pi) ? dphi : 2.*pi-dphi; // dphi lies in [-pi, pi]
        Pair_dpt[vtxind] = fabs(l1.Pt()-l2.Pt());
        Pair_3Dangle[vtxind] = (l1.Angle(l2.Vect()))/pi;
        
        for (int j=0; j<nPhotonCand; j++) {
          PairGamma_mass[vtxind][j] = sqrt(
            std::pow(l1.P() +l2.P() +PhotonCand_p[j], 2)
           -std::pow(l1.Px()+l2.Px()+PhotonCand_px[j], 2)
           -std::pow(l1.Py()+l2.Py()+PhotonCand_py[j], 2)
           -std::pow(l1.Pz()+l2.Pz()+PhotonCand_pz[j], 2));
          //std::cout << "Photon " << j << " added to give a mass = " << PairGamma_mass[vtxind][j] << std::endl;
        }
        PrimVertexCand_hasdil[vtxind] = 1;
        
	nCandidatesInEvent++;
	nCandidates++;
        foundPairInEvent = true;
      }

      vtxind++;
    }
    nFilteredPrimVertexCand = vtxind;
    if (verb_>2) edm::LogInfo("GammaGammaLL") << "end of the loop on the vertices : vtxind = " << vtxind;
  }
  nExtraTracks = etind;

  if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed Loop on vertices";
  
  /////
  if (!foundPairInEvent) return; // avoid to unpack RP/jet/MET if no dilepton candidate has been found
  /////
  
  if (fetchProtons_) {
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
        if (verb_>2) edm::LogInfo("GammaGammaLL") << "Proton track candidate with origin: (" << proton->getX0() << ", " << proton->getY0() << ", " << proton->getZ0() << ") extracted!";
      }
    }
    if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed TOTEM RP info retrieval stage. Got " << nLocalProtCand << " local track(s)";
  }

  // Get the Jet collection from the event
  // PAT
  iEvent.getByToken(jetToken_, jetColl_);
  for (jet=jetColl_->begin(); jet!=jetColl_->end() && nJetCand<MAX_JETS; jet++) {
    JetCand_e[nJetCand] = jet->energy();
    JetCand_px[nJetCand] = jet->px();
    JetCand_py[nJetCand] = jet->py();
    JetCand_pz[nJetCand] = jet->pz();
    JetCand_phi[nJetCand] = jet->phi();
    JetCand_eta[nJetCand] = jet->eta();
    totalJetEnergy += JetCand_e[nJetCand];
    // Find kinematics quantities associated to the highest energy jet
    if (JetCand_e[nJetCand]>HEJet_e) {
      HEJet_e = JetCand_e[nJetCand];
      HEJet_eta = JetCand_eta[nJetCand];
      HEJet_phi = JetCand_phi[nJetCand];
    }
    nJetCand++;
  }
  HighestJet_e = HEJet_e;
  HighestJet_eta = HEJet_eta;
  HighestJet_phi = HEJet_phi;
  SumJet_e = totalJetEnergy;

  if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed Loop on jets";

  // Missing ET
  iEvent.getByToken(metToken_, MET); 
  const edm::View<pat::MET>* metColl = MET.product(); 
  met = metColl->begin();
  
  Etmiss = met->et();
  Etmiss_phi = met->phi();
  Etmiss_x = met->px();
  Etmiss_y = met->py();
  Etmiss_z = met->pz();
  Etmiss_significance = met->significance();

  if (verb_>1) edm::LogInfo("GammaGammaLL") << "Passed MET retrieval stage";

  if (verb_>2 or printCandidates_) {
    std::cout << "Event " << Run << ":" << EventNum << " has " << nCandidatesInEvent << " leptons pair(s) candidate(s) (vertex mult. : " << nPrimVertexCand << ")" << std::endl;
  }
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaLL::beginJob()
{
  // Booking the ntuple

  tree_->Branch("Run", &Run, "Run/I");
  tree_->Branch("LumiSection", &LumiSection, "LumiSection/I");
  tree_->Branch("BX", &BX, "BX/I");
  tree_->Branch("EventNum", &EventNum, "EventNum/I");

  tree_->Branch("nHLT", &nHLT, "nHLT/I");
  tree_->Branch("HLT_Accept", HLT_Accept, "HLT_Accept[nHLT]/I");
  tree_->Branch("HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I");
  tree_->Branch("HLT_Name", &HLT_Name);
  *HLT_Name = triggersList_;
  
  if (fetchMuons_) {
    tree_->Branch("nMuonCand", &nMuonCand, "nMuonCand/I");
    tree_->Branch("MuonCand_px", MuonCand_px, "MuonCand_px[nMuonCand]/D");
    tree_->Branch("MuonCand_py", MuonCand_py, "MuonCand_py[nMuonCand]/D");
    tree_->Branch("MuonCand_pz", MuonCand_pz, "MuonCand_pz[nMuonCand]/D");
    tree_->Branch("MuonCand_p", MuonCand_p, "MuonCand_p[nMuonCand]/D");
    tree_->Branch("MuonCand_pt", MuonCand_pt, "MuonCand_pt[nMuonCand]/D");
    tree_->Branch("MuonCand_eta", MuonCand_eta, "MuonCand_eta[nMuonCand]/D");
    tree_->Branch("MuonCand_phi", MuonCand_phi, "MuonCand_phi[nMuonCand]/D");
    tree_->Branch("MuonCand_charge", MuonCand_charge, "MuonCand_charge[nMuonCand]/I");
    tree_->Branch("MuonCand_vtxx", MuonCand_vtxx, "MuonCand_vtxx[nMuonCand]/D");
    tree_->Branch("MuonCand_vtxy", MuonCand_vtxy, "MuonCand_vtxy[nMuonCand]/D");
    tree_->Branch("MuonCand_vtxz", MuonCand_vtxz, "MuonCand_vtxz[nMuonCand]/D");
    tree_->Branch("MuonCand_dxy", MuonCand_dxy, "MuonCand_dxy[nMuonCand]/D");
    tree_->Branch("MuonCand_dz", MuonCand_dz, "MuonCand_dz[nMuonCand]/D");
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
    if (runOnMC_) {
      tree_->Branch("nGenMuonCand", &nGenMuonCand, "nGenMuonCand/I");
      tree_->Branch("nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept, "nGenMuonCandOutOfAccept/I");    
      tree_->Branch("GenMuonCand_p", GenMuonCand_p, "GenMuonCand_p[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_px", GenMuonCand_px, "GenMuonCand_px[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_py", GenMuonCand_py, "GenMuonCand_py[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_pz", GenMuonCand_pz, "GenMuonCand_pz[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_pt", GenMuonCand_pt, "GenMuonCand_pt[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_eta", GenMuonCand_eta, "GenMuonCand_eta[nGenMuonCand]/D");
      tree_->Branch("GenMuonCand_phi", GenMuonCand_phi, "GenMuonCand_phi[nGenMuonCand]/D");
    }
  }
  
  if (fetchElectrons_) {
    tree_->Branch("nEleCand", &nEleCand, "nEleCand/I");
    tree_->Branch("EleCand_px", EleCand_px, "EleCand_px[nEleCand]/D");
    tree_->Branch("EleCand_py", EleCand_py, "EleCand_py[nEleCand]/D");
    tree_->Branch("EleCand_pz", EleCand_pz, "EleCand_pz[nEleCand]/D");
    tree_->Branch("EleCand_p", EleCand_p, "EleCand_p[nEleCand]/D");
    tree_->Branch("EleCand_e", EleCand_e, "EleCand_e[nEleCand]/D");
    tree_->Branch("EleCand_et", EleCand_et, "EleCand_et[nEleCand]/D");
    tree_->Branch("EleCand_eta", EleCand_eta, "EleCand_eta[nEleCand]/D");
    tree_->Branch("EleCand_phi", EleCand_phi, "EleCand_phi[nEleCand]/D");
    tree_->Branch("EleCand_charge", EleCand_charge, "EleCand_charge[nEleCand]/I");
    tree_->Branch("EleCand_vtxx", EleCand_vtxx, "EleCand_vtxx[nEleCand]/D");
    tree_->Branch("EleCand_vtxy", EleCand_vtxy, "EleCand_vtxy[nEleCand]/D");
    tree_->Branch("EleCand_vtxz", EleCand_vtxz, "EleCand_vtxz[nEleCand]/D");
    tree_->Branch("EleCandTrack_p", EleCandTrack_p, "EleCandTrack_p[nEleCand]/D");
    tree_->Branch("EleCandTrack_pt", EleCandTrack_pt, "EleCandTrack_pt[nEleCand]/D");
    tree_->Branch("EleCandTrack_eta", EleCandTrack_eta, "EleCandTrack_eta[nEleCand]/D");
    tree_->Branch("EleCandTrack_phi", EleCandTrack_phi, "EleCandTrack_phi[nEleCand]/D");
    tree_->Branch("EleCandTrack_vtxz", EleCandTrack_vtxz, "EleCandTrack_vtxz[nEleCand]/D");
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
    if (runOnMC_) {
      tree_->Branch("nGenEleCand", &nGenEleCand, "nGenEleCand/I");
      tree_->Branch("nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept, "nGenEleCandOutOfAccept/I");    
      tree_->Branch("GenEleCand_p", GenEleCand_p, "GenEleCand_p[nGenEleCand]/D");
      tree_->Branch("GenEleCand_px", GenEleCand_px, "GenEleCand_px[nGenEleCand]/D");
      tree_->Branch("GenEleCand_py", GenEleCand_py, "GenEleCand_py[nGenEleCand]/D");
      tree_->Branch("GenEleCand_pz", GenEleCand_pz, "GenEleCand_pz[nGenEleCand]/D");
      tree_->Branch("GenEleCand_pt", GenEleCand_pt, "GenEleCand_pt[nGenEleCand]/D");
      tree_->Branch("GenEleCand_eta", GenEleCand_eta, "GenEleCand_eta[nGenEleCand]/D");
      tree_->Branch("GenEleCand_phi", GenEleCand_phi, "GenEleCand_phi[nGenEleCand]/D");
    }
  }
  tree_->Branch("nPhotonCand", &nPhotonCand, "nPhotonCand/I");
  tree_->Branch("PhotonCand_pt", PhotonCand_pt, "PhotonCand_pt[nPhotonCand]/D");
  tree_->Branch("PhotonCand_eta", PhotonCand_eta, "PhotonCand_eta[nPhotonCand]/D");
  tree_->Branch("PhotonCand_phi", PhotonCand_phi, "PhotonCand_phi[nPhotonCand]/D");
  tree_->Branch("PhotonCand_r9", PhotonCand_r9, "PhotonCand_r9[nPhotonCand]/D");
  tree_->Branch("PhotonCand_drtrue", PhotonCand_drtrue, "PhotonCand_drtrue[nPhotonCand]/D"); 
  tree_->Branch("PhotonCand_detatrue", PhotonCand_detatrue, "PhotonCand_detatrue[nPhotonCand]/D"); 
  tree_->Branch("PhotonCand_dphitrue", PhotonCand_dphitrue, "PhotonCand_dphitrue[nPhotonCand]/D"); 
  if (runOnMC_) {
    tree_->Branch("nGenPhotCand", &nGenPhotCand, "nGenPhotCand/I");    
    tree_->Branch("nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept, "nGenPhotCandOutOfAccept/I");    
    tree_->Branch("GenPhotCand_e", GenPhotCand_e, "GenPhotCand_e[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_p", GenPhotCand_p, "GenPhotCand_p[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_pt", GenPhotCand_pt, "GenPhotCand_pt[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_eta", GenPhotCand_eta, "GenPhotCand_eta[nGenPhotCand]/D");
    tree_->Branch("GenPhotCand_phi", GenPhotCand_phi, "GenPhotCand_phi[nGenPhotCand]/D");
    tree_->Branch("nGenProtCand", &nGenProtCand, "nGenProtCand/I");    
    tree_->Branch("GenProtCand_p", GenProtCand_p, "GenProtCand_p[nGenProtCand]/D");
    tree_->Branch("GenProtCand_px", GenProtCand_px, "GenProtCand_px[nGenProtCand]/D");
    tree_->Branch("GenProtCand_py", GenProtCand_py, "GenProtCand_py[nGenProtCand]/D");
    tree_->Branch("GenProtCand_pz", GenProtCand_pz, "GenProtCand_pz[nGenProtCand]/D");
    tree_->Branch("GenProtCand_pt", GenProtCand_pt, "GenProtCand_pt[nGenProtCand]/D");
    tree_->Branch("GenProtCand_eta", GenProtCand_eta, "GenProtCand_eta[nGenProtCand]/D");
    tree_->Branch("GenProtCand_phi", GenProtCand_phi, "GenProtCand_phi[nGenProtCand]/D");
    tree_->Branch("GenProtCand_status", GenProtCand_status, "GenProtCand_status[nGenProtCand]/I");
  }
  
  // Primary vertices' information
  tree_->Branch("nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/I");
  tree_->Branch("nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nFilteredPrimVertexCand/I");
  tree_->Branch("PrimVertexCand_id", PrimVertexCand_id, "PrimVertexCand_id[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_hasdil", PrimVertexCand_hasdil, "PrimVertexCand_hasdil[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_matchedtracks", PrimVertexCand_matchedtracks, "PrimVertexCand_matchedtracks[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_unmatchedtracks", PrimVertexCand_unmatchedtracks, "PrimVertexCand_unmatchedtracks[nPrimVertexCand]/I");
  tree_->Branch("PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D");
  tree_->Branch("PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/I");

  // Lepton pairs' information
  tree_->Branch("Pair_candidates", Pair_candidates, "Pair_candidates[nPrimVertexCand][2]/I");
  tree_->Branch("Pair_lepton1", Pair_lepton1, "Pair_lepton1[nPrimVertexCand]/I");
  tree_->Branch("Pair_lepton2", Pair_lepton2, "Pair_lepton2[nPrimVertexCand]/I");
  tree_->Branch("Pair_mindist", Pair_mindist, "Pair_mindist[nPrimVertexCand]/D");
  tree_->Branch("Pair_p", Pair_p, "Pair_p[nPrimVertexCand]/D");
  tree_->Branch("Pair_pt", Pair_pt, "Pair_pt[nPrimVertexCand]/D");
  tree_->Branch("Pair_dpt", Pair_dpt, "Pair_dpt[nPrimVertexCand]/D");
  tree_->Branch("Pair_mass", Pair_mass, "Pair_mass[nPrimVertexCand]/D");
  tree_->Branch("Pair_eta", Pair_eta, "Pair_eta[nPrimVertexCand]/D");
  tree_->Branch("Pair_phi", Pair_phi, "Pair_phi[nPrimVertexCand]/D");
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

  if (!runOnMC_) {
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

  if (runOnMC_) {
    tree_->Branch("GenPair_p", &GenPair_p, "GenPair_p/D");
    tree_->Branch("GenPair_pt", &GenPair_pt, "GenPair_pt/D");
    tree_->Branch("GenPair_dpt", &GenPair_dpt, "GenPair_dpt/D");
    tree_->Branch("GenPair_mass", &GenPair_mass, "GenPair_mass/D");
    tree_->Branch("GenPair_eta", &GenPair_eta, "GenPair_eta/D");
    tree_->Branch("GenPair_phi", &GenPair_phi, "GenPair_phi/D");
    tree_->Branch("GenPair_dphi", &GenPair_dphi, "GenPair_dphi/D");
    tree_->Branch("GenPair_3Dangle", &GenPair_3Dangle, "GenPair_3Dangle[nPrimVertexCand]/D");
    /*tree_->Branch("HPS_acc420b1", &HPS_acc420b1, "HPS_acc420b1/D");
    tree_->Branch("HPS_acc220b1", &HPS_acc220b1, "HPS_acc220b1/D");
    tree_->Branch("HPS_acc420and220b1", &HPS_acc420and220b1, "HPS_acc420and220b1/D");
    tree_->Branch("HPS_acc420or220b1", &HPS_acc420or220b1, "HPS_acc420or220b1/D");
    tree_->Branch("HPS_acc420b2", &HPS_acc420b2, "HPS_acc420b2/D");
    tree_->Branch("HPS_acc220b2", &HPS_acc220b2, "HPS_acc220b2/D");
    tree_->Branch("HPS_acc420and220b2", &HPS_acc420and220b2, "HPS_acc420and220b2/D");
    tree_->Branch("HPS_acc420or220b2", &HPS_acc420or220b2, "HPS_acc420or220b2/D");*/
  }
  // Extra tracks on vertex's information
  tree_->Branch("nExtraTracks", &nExtraTracks, "nExtraTracks/I");
  tree_->Branch("ExtraTrack_vtxId", ExtraTrack_vtxId, "ExtraTrack_vtxId[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/I");
  tree_->Branch("ExtraTrack_p", ExtraTrack_p, "ExtraTrack_p[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_pt", ExtraTrack_pt, "ExtraTrack_pt[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_px", ExtraTrack_px, "ExtraTrack_px[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_py", ExtraTrack_py, "ExtraTrack_py[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_pz", ExtraTrack_pz, "ExtraTrack_pz[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_eta", ExtraTrack_eta, "ExtraTrack_eta[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_phi", ExtraTrack_phi, "ExtraTrack_phi[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_chi2", ExtraTrack_chi2, "ExtraTrack_chi2[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz, "ExtraTrack_vtxdxyz[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_vtxT", ExtraTrack_vtxT, "ExtraTrack_vtxT[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_vtxZ", ExtraTrack_vtxZ, "ExtraTrack_vtxZ[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_x", ExtraTrack_x, "ExtraTrack_x[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_y", ExtraTrack_y, "ExtraTrack_y[nExtraTracks]/D");
  tree_->Branch("ExtraTrack_z", ExtraTrack_z, "ExtraTrack_z[nExtraTracks]/D");
  tree_->Branch("nQualityExtraTrack", &nQualityExtraTrack, "nQualityExtraTrack/I");
  tree_->Branch("ClosestExtraTrack_vtxdxyz",ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz[nPrimVertexCand]/D");
  tree_->Branch("ClosestExtraTrack_id",ClosestExtraTrack_id,"ClosestExtraTrack_id[nPrimVertexCand]/I");
  tree_->Branch("ClosestHighPurityExtraTrack_vtxdxyz",ClosestHighPurityExtraTrack_vtxdxyz,"ClosestHighPurityExtraTrack_vtxdxyz[nPrimVertexCand]/D");
  tree_->Branch("ClosestHighPurityExtraTrack_id",ClosestHighPurityExtraTrack_id,"ClosestHighPurityExtraTrack_id[nPrimVertexCand]/I");
  
  // Jets/MET information
  tree_->Branch("nJetCand", &nJetCand, "nJetCand/I");
  tree_->Branch("JetCand_px", JetCand_px, "JetCand_px[nJetCand]/D");
  tree_->Branch("JetCand_py", JetCand_py, "JetCand_py[nJetCand]/D");
  tree_->Branch("JetCand_pz", JetCand_pz, "JetCand_pz[nJetCand]/D");
  tree_->Branch("JetCand_e", JetCand_e, "JetCand_e[nJetCand]/D");
  tree_->Branch("JetCand_eta", JetCand_eta, "JetCand_eta[nJetCand]/D");
  tree_->Branch("JetCand_phi", JetCand_phi, "JetCand_phi[nJetCand]/D");
  tree_->Branch("HighestJet_e", &HighestJet_e, "HighestJet_e/D");
  tree_->Branch("HighestJet_eta", &HighestJet_eta, "HighestJet_eta/D");
  tree_->Branch("HighestJet_phi", &HighestJet_phi, "HighestJet_phi/D");
  tree_->Branch("SumJet_e", &SumJet_e, "SumJet_e/D");
  tree_->Branch("Etmiss", &Etmiss, "Etmiss/D");
  tree_->Branch("Etmiss_phi", &Etmiss_phi, "Etmiss_phi/D");
  tree_->Branch("Etmiss_x", &Etmiss_x, "Etmiss_x/D");
  tree_->Branch("Etmiss_y", &Etmiss_y, "Etmiss_y/D");
  tree_->Branch("Etmiss_z", &Etmiss_z, "Etmiss_z/D"); 
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
  GenPair_p = GenPair_pt = GenPair_mass = GenPair_phi = GenPair_eta = -999.;
  GenPair_dphi = GenPair_dpt = GenPair_3Dangle = 0.;

  for (unsigned int i=0; i<MAX_LL; i++) {
    MuonCand_p[i] = MuonCand_px[i] = MuonCand_py[i] = MuonCand_pz[i] = -999.;
    MuonCand_pt[i] = MuonCand_eta[i] = MuonCand_phi[i] = -999.;
    MuonCand_charge[i] = -999;
    MuonCand_vtxx[i] = MuonCand_vtxy[i] = MuonCand_vtxz[i] = -999.;
    MuonCand_npxlhits[i] = MuonCand_nstatseg[i] = MuonCand_ntrklayers[i] = -999;
    MuonCand_dxy[i] = MuonCand_dz[i] = -999.;
    MuonCand_isglobal[i] = MuonCand_istracker[i] = MuonCand_isstandalone[i] = MuonCand_ispfmuon[i] = -999;
    MuonCand_istight[i] = -999;
    MuonCandTrack_nmuchits[i] = -999;
    MuonCandTrack_chisq[i] = -999.;
    EleCand_e[i] = EleCand_et[i] = EleCand_px[i] = EleCand_py[i] = EleCand_pz[i] = -999.;
    EleCand_p[i] = EleCand_phi[i] = EleCand_eta[i] = -999.;
    EleCand_charge[i] = -999;
    EleCand_vtxx[i] = EleCand_vtxy[i] = EleCand_vtxz[i] = -999.;
    EleCandTrack_p[i] = EleCandTrack_pt[i] = EleCandTrack_eta[i] = EleCandTrack_phi[i] = -999.; 
    EleCandTrack_vtxz[i] = -999.;
    EleCand_deltaPhi[i] = EleCand_deltaEta[i] = EleCand_HoverE[i] = -999.;
    EleCand_trackiso[i] = EleCand_ecaliso[i] = EleCand_hcaliso[i] = EleCand_sigmaIetaIeta[i] = -999.;
    EleCand_convDist[i] = EleCand_convDcot[i] = EleCand_ecalDriven[i] = -999.;
    EleCand_tightID[i] = EleCand_mediumID[i] = EleCand_looseID[i] = -1;
  }
  for (unsigned int i=0; i<MAX_ET && i<maxExTrkVtx_; i++) {
    ExtraTrack_p[i] = ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
    ExtraTrack_pt[i] = ExtraTrack_eta[i] = ExtraTrack_phi[i] = -999.;
    ExtraTrack_charge[i] = ExtraTrack_ndof[i] = -999;
    ExtraTrack_chi2[i] = ExtraTrack_vtxdxyz[i] = -999.;
    ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
    ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
    ExtraTrack_vtxId[i] = -1;
  }
  for (unsigned int i=0; i<MAX_PAIRS; i++) {
    Pair_candidates[i][0] = Pair_candidates[i][1] = -1;
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
  }
  for (unsigned int i=0; i<MAX_PHO; i++) { 
    PhotonCand_p[i] = PhotonCand_pt[i] = -999.;
    PhotonCand_px[i] = PhotonCand_py[i] = PhotonCand_pz[i] = -999.;
    PhotonCand_eta[i] = PhotonCand_phi[i] = PhotonCand_r9[i] = -999.;
    PhotonCand_detatrue[i] = PhotonCand_dphitrue[i] = PhotonCand_drtrue[i] = -999.;
  }
  for (unsigned int i=0; i<MAX_JETS; i++) {
    JetCand_px[i] = JetCand_py[i] = JetCand_pz[i] = -999.;
    JetCand_e[i] = JetCand_phi[i] = JetCand_eta[i] = -999;
  }
  HighestJet_e = HighestJet_eta = HighestJet_phi = -999.;
  HEJet_e = SumJet_e = totalJetEnergy = 0.;
  Etmiss = Etmiss_x = Etmiss_y = Etmiss_z = Etmiss_significance = -999.;
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
