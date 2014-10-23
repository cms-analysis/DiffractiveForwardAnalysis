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

// WP80 cuts
const float MAX_MissingHits      = 0.0;
const float MIN_Dist             = 0.02;
const float MIN_Dcot             = 0.02;
const float cut_EB_trackRel03    = 0.09;
const float cut_EB_ecalRel03     = 0.07;
const float cut_EB_hcalRel03     = 0.10;
const float cut_EB_sigmaIetaIeta = 0.01;
const float cut_EB_deltaPhi      = 0.06;
const float cut_EB_deltaEta      = 0.004;
const float cut_EB_HoverE        = 0.04;
const float cut_EE_trackRel03    = 0.04;
const float cut_EE_ecalRel03     = 0.05;
const float cut_EE_hcalRel03     = 0.025;
const float cut_EE_sigmaIetaIeta = 0.03;
const float cut_EE_deltaPhi      = 0.03;
const float cut_EE_deltaEta      = 0.007; 
const float cut_EE_HoverE        = 0.025;

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaLL.h"
//
// constructors and destructor
//
GammaGammaLL::GammaGammaLL(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  _fetchMuons = false;
  _fetchElectrons = false;
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outfilename", "output.root");
  
  hltMenuLabel_ = iConfig.getParameter<std::string>("HLTMenuLabel");
  triggersList_ = iConfig.getParameter<std::vector<std::string> >("TriggersList");
  _hlts = new HLTmatches(triggersList_);
  nHLT = triggersList_.size();
	
  recoVertexLabel_ = iConfig.getParameter<edm::InputTag>("RecoVertexLabel");

  maxExTrkVtx_ = iConfig.getUntrackedParameter<UInt_t>("maxExtraTracks", 1000);
  
  // Generator level
  sqrts_ = iConfig.getParameter<Double_t>("SqrtS");
  runOnMC_ = iConfig.getUntrackedParameter<bool>("RunOnMC", true);
  genLabel_ = iConfig.getParameter<edm::InputTag>("GenParticlesCollectionLabel");
  minPtMC_ = iConfig.getUntrackedParameter<Double_t>("MCAcceptPtCut", 20.);
  minEtaMC_ = iConfig.getUntrackedParameter<Double_t>("MCAcceptEtaCut", 2.5);
  
  // Pileup input tags
  pileupLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("pileupInfo", std::string("addPileupInfo"));
  mcPileupFile_ = iConfig.getUntrackedParameter<std::string>("mcpufile", "PUHistos.root");
  mcPileupPath_ = iConfig.getUntrackedParameter<std::string>("mcpupath", "pileup");
  dataPileupFile_ = iConfig.getUntrackedParameter<std::string>("datapufile", "PUHistos_duplicated.root");
  dataPileupPath_ = iConfig.getUntrackedParameter<std::string>("datapupath", "pileup");
  
  // Leptons input tags
  leptonsType_ = iConfig.getParameter< std::vector<std::string> >("LeptonsType");
  for (i=0; i<leptonsType_.size(); i++) {
    if (leptonsType_[i]=="Muon") _fetchMuons = true;
    else if (leptonsType_[i]=="Electron") _fetchElectrons = true;
  }
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalMuonCollectionLabel", std::string("muons"));
  eleLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalEleCollectionLabel", std::string("gsfElectrons"));
  isoValLabel_ = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");
  beamSpotLabel_ = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  conversionsLabel_ = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
  rhoIsoLabel_ = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
  printCandidates_ = iConfig.getUntrackedParameter<bool>("PrintCandidates", false);
  
  // Particle flow input tag
  pflowLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("PFLabel", std::string("particleFlow"));
  
  // Jets/MET input tags
  jetLabel_ = iConfig.getParameter<edm::InputTag>("JetCollectionLabel");
  metLabel_ = iConfig.getParameter<edm::InputTag>("MetLabel");
  
  file_ = new TFile(outputFile_.c_str(), "recreate");
  file_->cd();
  // tree definition
  tree_ = new TTree("ntp1", "ntp1");

  //log_hist = new TH1D("log", "", 500, -25., 25.);
  logfile = new std::ofstream("log_file.out");
      
  // HPS acceptance file readout definition
  if (runOnMC_) {
    // edm::FileInPath myDataFile("FastSimulation/ProtonTaggers/data/acceptance_420_220.root");  
    myDataFile = new edm::FileInPath("FastSimulation/ForwardDetectors/data/acceptance_420_220.root");
    fullAcceptancePath = myDataFile->fullPath();
    std::cout << "Opening " << fullAcceptancePath << std::endl;
    f = new TFile(fullAcceptancePath.c_str());
    if (f->Get("description") != NULL) {
      std::cout << "Description found: " << f->Get("description")->GetTitle() << std::endl;
      std::cout << "Reading acceptance tables " << std::endl;
    }
    helper420beam1.Init(*f, "a420");
    helper420beam2.Init(*f, "a420_b2");
    helper220beam1.Init(*f, "a220");
    helper220beam2.Init(*f, "a220_b2");
    helper420a220beam1.Init(*f, "a420a220");
    helper420a220beam2.Init(*f, "a420a220_b2");

    LumiWeights = new edm::LumiReWeighting(mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_);
  }

}


GammaGammaLL::~GammaGammaLL()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  const edm::ParameterSet& pset = edm::getProcessParameterSet();
  TList *list = tree_->GetUserInfo();
  list->Add(new TObjString(pset.dump().c_str()));
  file_->Write();
  file_->Close();

  logfile->close();

  delete _hlts;
  delete tree_;

}


//
// member functions
//

void
GammaGammaLL::LookAtTriggers(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Int_t trigNum;
	
  // Get the trigger information from the event
  iEvent.getByLabel(edm::InputTag("TriggerResults","",hltMenuLabel_),hltResults_) ; 
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltResults_);

  for (unsigned int i=0; i<trigNames.size(); i++) {
    //std::cout << "--> " << trigNames.triggerNames().at(i) << std::endl;
    trigNum = _hlts->TriggerNum(trigNames.triggerNames().at(i));
    if (trigNum==-1) continue; // Trigger didn't match the interesting ones
    HLT_Accept[trigNum] = hltResults_->accept(i) ? 1 : 0;
    HLT_Prescl[trigNum] = hltConfig_.prescaleValue(iEvent, iSetup, trigNames.triggerNames().at(i));
    
    //LF FIXME need to think about that implementation...
    /*if (trigNames.triggerNames().at(i).find("CaloIdL")) {} // Leading lepton
      else if (trigNames.triggerNames().at(i).find("CaloIdT")) {} // Trailing lepton*/
    //std::cout << "*-------> " << trigNames.triggerNames().at(i).substr(0, trigNames.triggerNames().at(i).find_last_of("_"));
    //HLT_LeadingLepton_Prescl[] = hltConfig_.prescaleValue(event, iSetup, "HLT_Mu8Ele17L");  
  }
}

// ------------ method called for each event  ------------
void
GammaGammaLL::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //std::cout << "Beginning First init" << std::endl;

  // First initialization of the variables
  nCandidatesInEvent = 0;
  nPrimVertexCand = nFilteredPrimVertexCand = -1;
  nMuonCand = nEleCand = nLeptonCand = 0;
  nExtraTracks = nQualityExtraTrack = 0;
  nJetCand = 0;
  nGenMuonCand = nGenMuonCandOutOfAccept = 0;
  nGenEleCand = nGenEleCandOutOfAccept = 0;
  nGenPhotCand = nGenPhotCandOutOfAccept = 0;
  nGenProtCand = 0;
  nPFPhotonCand = 0;

  HPS_acc420b1 = HPS_acc220b1 = HPS_acc420and220b1 = HPS_acc420or220b1 = -1;
  HPS_acc420b2 = HPS_acc220b2 = HPS_acc420and220b2 = HPS_acc420or220b2 = -1;
  GenPair_p = GenPair_pt = GenPair_mass = GenPair_phi = GenPair_eta = -999.;
  GenPair_dphi = GenPair_dpt = GenPair_3Dangle = 0.;

  closesttrkdxyz = closesthighpuritytrkdxyz = 999.;
  for (i=0; i<MAX_LL; i++) {
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
    EleCand_wp80[i] = EleCand_mediumID[i] = EleCand_looseID[i] = -1;
  }
  for (i=0; i<MAX_ET && i<(UInt_t)maxExTrkVtx_; i++) {
    ExtraTrack_p[i] = ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
    ExtraTrack_pt[i] = ExtraTrack_eta[i] = ExtraTrack_phi[i] = -999.;
    ExtraTrack_charge[i] = ExtraTrack_ndof[i] = -999;
    ExtraTrack_chi2[i] = ExtraTrack_vtxdxyz[i] = -999.;
    ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
    ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
    ExtraTrack_vtxId[i] = -1;
  }
  for (i=0; i<MAX_PAIRS; i++) {
    Pair_candidates[i][0] = Pair_candidates[i][1] = -1;
    Pair_mindist[i] = Pair_p[i] = Pair_pt[i] = Pair_mass[i] = Pair_phi[i] = Pair_eta[i] = -999.;
    Pair_dphi[i] = Pair_dpt[i] = Pair_3Dangle[i] = -999.;
    Pair_extratracks1mm[i] = Pair_extratracks2mm[i] = Pair_extratracks3mm[i] = 0;
    Pair_extratracks4mm[i] = Pair_extratracks5mm[i] = Pair_extratracks1cm[i] = 0;
    Pair_extratracks2cm[i] = Pair_extratracks3cm[i] = Pair_extratracks4cm[i] = 0;
    Pair_extratracks5cm[i] = Pair_extratracks10cm[i] = 0;    
  }
  for (i=0; i<MAX_VTX; i++) {
    PrimVertexCand_id[i] = -1;
    PrimVertexCand_tracks[i] = PrimVertexCand_matchedtracks[i] = PrimVertexCand_unmatchedtracks[i] = PrimVertexCand_hasdil[i] = 0;
    PrimVertexCand_x[i] = PrimVertexCand_y[i] = PrimVertexCand_z[i] = -999.;
    PrimVertexCand_chi2[i] = PrimVertexCand_ndof[i] = -999.;
  }
  for (i=0; i<MAX_PHO; i++) { 
    PFPhotonCand_p[i] = PFPhotonCand_pt[i] = -999.;
    PFPhotonCand_px[i] = PFPhotonCand_py[i] = PFPhotonCand_pz[i] = -999.;
    PFPhotonCand_eta[i] = PFPhotonCand_phi[i] = -999.;
    PFPhotonCand_detatrue[i] = PFPhotonCand_dphitrue[i] = PFPhotonCand_drtrue[i] = -999.;
  }
  for (i=0; i<MAX_JETS; i++) {
    JetCand_px[i] = JetCand_py[i] = JetCand_pz[i] = -999.;
    JetCand_e[i] = JetCand_phi[i] = JetCand_eta[i] = -999;
  }
  HighestJet_e = HighestJet_eta = HighestJet_phi = -999.;
  HEJet_e = SumJet_e = totalJetEnergy = 0.;
  Etmiss = Etmiss_x = Etmiss_y = Etmiss_z = Etmiss_significance = -999.;
  
  Weight = 1.;
  
  foundPairInEvent = false;
  foundGenCandPairInEvent = false;
  
  muonsMomenta.clear();
  electronsMomenta.clear();
  
  _leptonptmp = new TLorentzVector();
  
  //std::cout << "Passed First init of the variables" << std::endl;

  // Run and BX information
  BX = iEvent.bunchCrossing();
  Run = iEvent.id().run();
  LumiSection = iEvent.luminosityBlock();
  EventNum = iEvent.id().event();
  
  // High level trigger information retrieval  
  LookAtTriggers(iEvent, iSetup);
  
  // beam spot information
  iEvent.getByLabel(beamSpotLabel_, beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // Get the vertex collection from the event
  iEvent.getByLabel(recoVertexLabel_, recoVertexColl);
  const reco::VertexCollection* vertices = recoVertexColl.product();

  // rho for isolation 
  iEvent.getByLabel(rhoIsoLabel_, rhoIso_h); 
  rhoIso = *(rhoIso_h.product()); 
  
  //std::cout << "Passed Isolation" << std::endl;

  // Generator level information
  if (runOnMC_) {
    iEvent.getByLabel(genLabel_, genPartColl);
    
    for (genPart=genPartColl->begin(); genPart!=genPartColl->end(); genPart++) {
      if (genPart->pt()<minPtMC_ || (minEtaMC_!=-1. && fabs(genPart->eta())>minEtaMC_)) {
        if (fabs(genPart->pdgId())==13) nGenMuonCandOutOfAccept++;
        if (fabs(genPart->pdgId())==11) nGenEleCandOutOfAccept++;
        if (fabs(genPart->pdgId())==22) nGenPhotCandOutOfAccept++;
        continue;
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
      foundGenCandPairInEvent = false;
      if (_fetchElectrons && _fetchMuons) { // Looks at electron+muon
        if(nGenMuonCand!=1 || nGenEleCand!=1) continue; // FIXME maybe a bit tight according to the newer PU conditions?
        l1.SetXYZM(GenMuonCand_px[0], GenMuonCand_py[0], GenMuonCand_pz[0], MASS_MU);
        l2.SetXYZM(GenEleCand_px[0], GenEleCand_py[0], GenEleCand_pz[0], MASS_E);
        foundGenCandPairInEvent = true;
      }
      else if (_fetchElectrons) { // Looks at dielectrons
        if(nGenEleCand!=2) continue; // FIXME maybe a bit tight according to the newer PU conditions?
        l1.SetXYZM(GenEleCand_px[0], GenEleCand_py[0], GenEleCand_pz[0], MASS_E);      	
        l2.SetXYZM(GenEleCand_px[1], GenEleCand_py[1], GenEleCand_pz[1], MASS_E);      	
        foundGenCandPairInEvent = true;      	
      }
      else if (_fetchMuons) { // Looks at dimuons
        if(nGenMuonCand!=2) continue; // FIXME maybe a bit tight according to the newer PU conditions?
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
      if(genPart->pdgId()==2212 && fabs(genPart->pz())>3000.) {
        // Kinematic quantities computation
        // xi = fractional momentum loss
        if (genPart->pz()>0.) xi = 1.-genPart->pz()/sqrts_;
        else xi = 1.+genPart->pz()/sqrts_;
        t = -(std::pow(genPart->pt(), 2)+std::pow(MASS_P*xi, 2))/(1.-xi);
        
        // HPS acceptance computation
        if(genPart->pz()>0.) {
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
        }
      }
    }
  }

  // Pileup information
  if (runOnMC_) {
    iEvent.getByLabel(pileupLabel_, pileupInfo);

    // This part is optional if the distributions are already generated
    sum_nvtx = 0;
    npv = npvtrue = npvm1true = npvp1true = npv0true = npv0 = 0;
    for(PVI=pileupInfo->begin(); PVI!=pileupInfo->end(); PVI++) {
      beamXing = PVI->getBunchCrossing();

      npv = PVI->getPU_NumInteractions();
      npvtrue = PVI->getTrueNumInteractions();
      sum_nvtx += npvtrue;
      
      if(beamXing == -1) npvm1true+=npvtrue;
      if(beamXing == 0) {
        npv0 += npv;
        npv0true += npvtrue;
      }
      if(beamXing == 1) npvp1true+=npvtrue;
    }
    nTruePUforPUWeightBX0 = npv0true;
    //
    
    const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
    Weight = LumiWeights->weight(*iEventB);
  }
  
  //std::cout << "Passed Pileup" << std::endl;

  // Get the muons collection from the event
  if (_fetchMuons) {
    iEvent.getByLabel(muonLabel_, muonColl);
    for (muon=muonColl->begin(); muon!=muonColl->end() && nMuonCand<MAX_MUONS; muon++) {
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
	_leptonptmp->SetXYZM(muon->innerTrack()->px(), muon->innerTrack()->py(), muon->innerTrack()->pz(), muon->mass());
      }
      else {
	_leptonptmp->SetXYZM(muon->px(), muon->py(), muon->pz(), muon->mass());
      }
      muonsMomenta.insert(std::pair<Int_t,TLorentzVector>(nMuonCand, *_leptonptmp));

      if (MuonCand_isglobal[nMuonCand] && MuonCand_istracker[nMuonCand]) {
	MuonCandTrack_nmuchits[nMuonCand] = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
	MuonCandTrack_chisq[nMuonCand] = muon->globalTrack()->normalizedChi2();
	istight = true;
	istight&= MuonCand_ispfmuon[nMuonCand];
	istight&= (MuonCandTrack_chisq[nMuonCand]<10.);
	istight&= (MuonCandTrack_nmuchits[nMuonCand]>=1);
	istight&= (MuonCand_nstatseg[nMuonCand]>=2);
	istight&= (MuonCand_dxy[nMuonCand]<.2);
	istight&= (MuonCand_npxlhits[nMuonCand]>0);
	istight&= (MuonCand_ntrklayers[nMuonCand]>5);
	MuonCand_istight[nMuonCand] = istight;
      } 

      nMuonCand++;
    }
  }
  
  //std::cout << "Passed Muon" << std::endl;

  // Get the electrons collection from the event
  if (_fetchElectrons) {
    iEvent.getByLabel(eleLabel_, eleColl);
    // RECO electrons
    iEvent.getByLabel(InputTag("gsfElectrons"), eleCollRECO);
    // iso deposits
    IsoDepositVals isoVals(isoValLabel_.size());
    // New 2012 electron ID variables conversions
    iEvent.getByLabel(conversionsLabel_, conversions_h);

    for (size_t j = 0; j<isoValLabel_.size(); j++) { 
      iEvent.getByLabel(isoValLabel_[j], isoVals[j]); 
    }
        
    for (electron=eleColl->begin(); electron!=eleColl->end() && nEleCand<MAX_ELE; electron++) {
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
  	  
      _leptonptmp->SetXYZM(electron->px(), electron->py(), electron->pz(), electron->mass());

      if(electron->closestCtfTrackRef().isNonnull()) { // Only for PAT::Electron
        EleCandTrack_p[nEleCand] = electron->closestCtfTrackRef()->p(); 
        EleCandTrack_pt[nEleCand] = electron->closestCtfTrackRef()->pt();  
        EleCandTrack_eta[nEleCand] = electron->closestCtfTrackRef()->eta();  
        EleCandTrack_phi[nEleCand] = electron->closestCtfTrackRef()->phi();  
        EleCandTrack_vtxz[nEleCand] = electron->closestCtfTrackRef()->vertex().z();
	_leptonptmp->SetPtEtaPhiM(EleCandTrack_pt[nEleCand], EleCandTrack_eta[nEleCand], EleCandTrack_phi[nEleCand], electron->mass());
      }

      electronsMomenta.insert(std::pair<Int_t,TLorentzVector>(nEleCand, *_leptonptmp));

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

      EleCand_wp80[nEleCand] = 0;
      EleCand_mediumID[nEleCand] = 0;
      EleCand_looseID[nEleCand] = 0;
      bool select = false;
      bool selectmedium = false;
      bool selectloose = false;
      if (electron->ecalDrivenSeed()==1)  select = true;
      if (EleCand_convDist[nEleCand]>=MIN_Dist or EleCand_convDcot[nEleCand]>=MIN_Dcot) select = true;
      if (electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=MAX_MissingHits) select = true;
      if(select) {
        if(electron->isEB()) {
          if (fabs(EleCand_deltaPhi[nEleCand])>cut_EB_deltaPhi) select = false;
          if (fabs(EleCand_deltaEta[nEleCand])>cut_EB_deltaEta) select = false;
          if (EleCand_HoverE[nEleCand]>cut_EB_HoverE) select = false;
          if (EleCand_trackiso[nEleCand]>cut_EB_trackRel03) select = false;
          if (EleCand_ecaliso[nEleCand]>cut_EB_ecalRel03) select = false;
          if (EleCand_hcaliso[nEleCand]>cut_EB_hcalRel03) select = false;
          if (EleCand_sigmaIetaIeta[nEleCand] > cut_EB_sigmaIetaIeta) select = false;
        }
        else if(electron->isEE()) {
          if (fabs(EleCand_deltaPhi[nEleCand])>cut_EE_deltaPhi) select = false;
          if (fabs(EleCand_deltaEta[nEleCand])>cut_EE_deltaEta) select = false;
          if (EleCand_HoverE[nEleCand]>cut_EE_HoverE) select = false;
          if (EleCand_trackiso[nEleCand]>cut_EE_trackRel03) select = false;
          if (EleCand_ecaliso[nEleCand]>cut_EE_ecalRel03) select = false;
          if (EleCand_hcaliso[nEleCand]>cut_EE_hcalRel03) select = false;
          if (EleCand_sigmaIetaIeta[nEleCand] > cut_EE_sigmaIetaIeta) select = false;
        }
      }
      if(select) EleCand_wp80[nEleCand] = 1;
      
      // get reference to RECO electron
      reco::GsfElectronRef ele(eleCollRECO, nEleCand);
      
      iso_ch =  (*(isoVals)[0])[ele];
      iso_em = (*(isoVals)[1])[ele];
      iso_nh = (*(isoVals)[2])[ele];

      if(PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, ele) == true) {
        selectmedium = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, recoVertexColl, iso_ch, iso_em, iso_nh, rhoIso, ElectronEffectiveArea::kEleEAData2012);
        selectloose = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele, conversions_h, beamSpot, recoVertexColl, iso_ch, iso_em, iso_nh, rhoIso, ElectronEffectiveArea::kEleEAData2012);
      }
      if(selectmedium) EleCand_mediumID[nEleCand] = 1;
      if(selectloose) EleCand_looseID[nEleCand] = 1;
      
      nEleCand++;
    }
  }

  delete _leptonptmp;
  
  nLeptonCand += (_fetchMuons) ? nMuonCand : 0;
  nLeptonCand += (_fetchElectrons) ? nEleCand : 0;
  
  //std::cout << "Passed Electron" << std::endl;
  
  // Get the PFlow collection from the event
  iEvent.getByLabel(pflowLabel_, pflowColl);
  for(pflow=pflowColl->begin(); pflow!=pflowColl->end(); pflow++) { 
    parttype = reco::PFCandidate::ParticleType(pflow->particleId()); 
    if(parttype==4 && nPFPhotonCand<MAX_PHO) { 
      if(nPFPhotonCand==0) { 
        leadingphotpx = pflow->px(); 
        leadingphotpy = pflow->py();  
        leadingphotpz = pflow->pz();  
        leadingphotp = pflow->p();  
      } 
      PFPhotonCand_p[nPFPhotonCand] = pflow->p();
      PFPhotonCand_px[nPFPhotonCand] = pflow->px();
      PFPhotonCand_py[nPFPhotonCand] = pflow->py();
      PFPhotonCand_pz[nPFPhotonCand] = pflow->pz();
      PFPhotonCand_pt[nPFPhotonCand] = pflow->pt();
      PFPhotonCand_eta[nPFPhotonCand] = pflow->eta();
      PFPhotonCand_phi[nPFPhotonCand] = pflow->phi();
      PFPhotonCand_drtrue[nPFPhotonCand] = -999.;
      PFPhotonCand_detatrue[nPFPhotonCand] = -999.;
      PFPhotonCand_dphitrue[nPFPhotonCand] = -999.;
      photdr = 999.;
      endphotdr = endphotdeta = endphotdphi = 999.;
      if (runOnMC_) {
	for (Int_t j=0; j<nGenPhotCand; j++) { // matching with the 'true' photon object from MC
	  photdeta = (PFPhotonCand_eta[nPFPhotonCand]-GenPhotCand_eta[j]);
	  photdphi = (PFPhotonCand_phi[nPFPhotonCand]-GenPhotCand_phi[j]);
	  photdr = sqrt(std::pow(photdeta, 2)+std::pow(photdphi, 2));
	  if (photdr<endphotdr) {
	    endphotdr = photdr;
	    endphotdeta = photdeta;
	    endphotdphi = photdphi;
	  }
	}
	PFPhotonCand_detatrue[nPFPhotonCand] = endphotdeta;
	PFPhotonCand_dphitrue[nPFPhotonCand] = endphotdphi;
	PFPhotonCand_drtrue[nPFPhotonCand] = endphotdr;
      }
      nPFPhotonCand++;
    } 
  }
  //std::cout << "Passed PF photons" << std::endl;

  vtxind = 0;
  if (nLeptonCand>=2) {
    // Enough leptons candidates to go deeper and analyze the primary vertices
    
    _leptonType = new TString();
    
    etind = 0;
    nPrimVertexCand = vertices->size();

    for (vertex=vertices->begin(); vertex!=vertices->end() && vtxind<MAX_VTX; ++vertex, ++vtxind) {
      PrimaryVertex *_vtx;
      Int_t leptonId_;

      //(*logfile) << vertex->z() << std::endl;

      nLeptonsInPrimVertex = 0;
      nExtraTracks = 0;
      nQualityExtraTrack = 0;
      foundPairOnVertex = false;
      
      _vtx = new PrimaryVertex(leptonsType_, muonsMomenta, electronsMomenta);
      _vtx->SetPosition(vertex->x(), vertex->y(), vertex->z());

      PrimVertexCand_id[vtxind] = vtxind;
      PrimVertexCand_x[vtxind] = _vtx->Position.X();
      PrimVertexCand_y[vtxind] = _vtx->Position.Y();
      PrimVertexCand_z[vtxind] = _vtx->Position.Z();
      PrimVertexCand_chi2[vtxind] = vertex->chi2();
      PrimVertexCand_ndof[vtxind] = vertex->ndof();
      
      closesttrkid = closesthighpuritytrkid = -1;
      // Loop on all the tracks matched with this vertex
      reco::Vertex::trackRef_iterator _track;
      for (_track=vertex->tracks_begin(); _track!=vertex->tracks_end() && etind<MAX_ET; _track++) {
	leptonId_ = _vtx->AddTrack((*_track).castTo<reco::TrackRef>(), *_leptonType);
	if (leptonId_==-1) { // Track was not matched to any of the leptons in the collection
          ExtraTrack_vtxId[etind] = vtxind;
	  vtxdst = sqrt(std::pow(((*_track)->vertex().x()-_vtx->Position.X()),2)+
			std::pow(((*_track)->vertex().y()-_vtx->Position.Y()),2)+
			std::pow(((*_track)->vertex().z()-_vtx->Position.Z()),2));
	  
	  ExtraTrack_purity[etind] = (*_track)->quality(reco::TrackBase::highPurity);
	  ExtraTrack_nhits[etind] = (*_track)->numberOfValidHits();
	  
	  ExtraTrack_p[etind] = (*_track)->p();
	  ExtraTrack_px[etind] = (*_track)->px();
	  ExtraTrack_py[etind] = (*_track)->py();
	  ExtraTrack_pz[etind] = (*_track)->pz();
	  ExtraTrack_pt[etind] = (*_track)->pt();
	  ExtraTrack_eta[etind] = (*_track)->eta();
	  ExtraTrack_phi[etind] = (*_track)->phi();
	  ExtraTrack_charge[etind] = (*_track)->charge();
	  ExtraTrack_chi2[etind] = (*_track)->chi2();
	  ExtraTrack_ndof[etind] = (*_track)->ndof();
          ExtraTrack_vtxdxyz[etind] = vtxdst;
          ExtraTrack_vtxT[etind] = sqrt(std::pow((*_track)->vertex().x()-_vtx->Position.X(),2)+
					std::pow((*_track)->vertex().y()-_vtx->Position.Y(),2));
          ExtraTrack_vtxZ[etind] = fabs((*_track)->vertex().z()-_vtx->Position.Z());
          ExtraTrack_x[etind] = (*_track)->vertex().x();
          ExtraTrack_y[etind] = (*_track)->vertex().y();
          ExtraTrack_z[etind] = (*_track)->vertex().z();
          
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
              closesthighpuritytrkid = etind;
	    }
	  }
          etind++;
        }
        else {
#ifdef DEBUG
          std::cout << "-- Lepton match on vertex " << vtxind << " --> " << leptonId_ << std::endl;
#endif
          nLeptonsInPrimVertex++;
        }
      }
      ClosestExtraTrack_vtxdxyz[vtxind] = closesttrkdxyz;
      ClosestExtraTrack_id[vtxind] = closesttrkid;
      ClosestHighPurityExtraTrack_vtxdxyz[vtxind] = closesthighpuritytrkdxyz;
      ClosestHighPurityExtraTrack_id[vtxind] = closesttrkid;
      
      PrimVertexCand_tracks[vtxind] = _vtx->nTracks;
      PrimVertexCand_matchedtracks[vtxind] = _vtx->nMatchedTracks;
      PrimVertexCand_unmatchedtracks[vtxind] = _vtx->nUnmatchedTracks;

      if (nLeptonsInPrimVertex<2) continue;
      
      // At this stage we have at least two matched leptons track on the vertex
      Pair_candidates[vtxind][0] = -1;
      Pair_candidates[vtxind][1] = -1;
      
      if (_fetchElectrons && _fetchMuons) { // Looks at electron+muon
	// Not enough muons or electrons candidates on the vertex
	if (_vtx->Electrons()==0 or _vtx->Muons()==0) {
#ifdef DEBUG
          std::cout << "Not enough electrons (" << _vtx->Electrons() << ") or muons (" << _vtx->Muons() << ") arising from the primary vertex !" << std::endl;
#endif
	  continue;
	}
	minDist = 999.;
	for(UInt_t i=0; i<_vtx->MatchedMuons.size(); i++) {
	  lep1 = _vtx->MatchedMuons[i];

	  TVector3 vtxlep1(MuonCand_vtxx[lep1], MuonCand_vtxy[lep1], MuonCand_vtxz[lep1]);
	  MuonCand_dz[lep1] = _vtx->dZ(vtxlep1, lep1);
	  MuonCand_istight[lep1]&= (MuonCand_dz[lep1]<.5);

	  for(UInt_t j=0; j<_vtx->MatchedElectrons.size(); j++) {
            lep2 = _vtx->MatchedElectrons[j];
            if (MuonCand_charge[lep1]*EleCand_charge[lep2]>0) {
              continue;
            }
            foundPairOnVertex = true;
            leptonsDist=sqrt(pow(MuonCand_vtxx[lep1]-EleCand_vtxx[lep2],2)+
			     pow(MuonCand_vtxy[lep1]-EleCand_vtxy[lep2],2)+
			     pow(MuonCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
            if (leptonsDist<minDist) {
              minDist = leptonsDist;
              Pair_candidates[vtxind][0] = lep1;
              Pair_candidates[vtxind][1] = lep2;
            }
          }
        }
        if (Pair_candidates[vtxind][0]==-1 || Pair_candidates[vtxind][1]==-1) {
          foundPairOnVertex = false;
        }
        l1.SetXYZM(MuonCand_px[Pair_candidates[vtxind][0]],
                   MuonCand_py[Pair_candidates[vtxind][0]],
                   MuonCand_pz[Pair_candidates[vtxind][0]],
                   MASS_MU);
        l2.SetXYZM(EleCand_px[Pair_candidates[vtxind][1]],
                   EleCand_py[Pair_candidates[vtxind][1]],
                   EleCand_pz[Pair_candidates[vtxind][1]],
                   MASS_E);
      }
      else if (_fetchElectrons) { // Looks at dielectrons
	// Not enough electrons candidates on the vertex
	if (_vtx->MatchedElectrons.size()<2) continue;
        minDist = 999.;
        for(UInt_t i=0; i<_vtx->MatchedElectrons.size(); i++) {
          lep1 = _vtx->MatchedElectrons[i];
          for(UInt_t j=i+1; j<_vtx->MatchedElectrons.size(); j++) {
            lep2 = _vtx->MatchedElectrons[j];
            if (EleCand_charge[lep1]*EleCand_charge[lep2]>0) {
              continue;
            }
            foundPairOnVertex = true;
            leptonsDist=sqrt(pow(EleCand_vtxx[lep1]-EleCand_vtxx[lep2],2)+
			     pow(EleCand_vtxy[lep1]-EleCand_vtxy[lep2],2)+
			     pow(EleCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
            if (leptonsDist<minDist) {
              minDist = leptonsDist;
              Pair_candidates[vtxind][0] = lep1;
              Pair_candidates[vtxind][1] = lep2;
            }
          }
        }
        if (Pair_candidates[vtxind][0]==-1 || Pair_candidates[vtxind][1]==-1) {
          foundPairOnVertex = false;
        }
        l1.SetXYZM(EleCand_px[Pair_candidates[vtxind][0]],
                   EleCand_py[Pair_candidates[vtxind][0]],
                   EleCand_pz[Pair_candidates[vtxind][0]],
                   MASS_E);
        l2.SetXYZM(EleCand_px[Pair_candidates[vtxind][1]],
                   EleCand_py[Pair_candidates[vtxind][1]],
                   EleCand_pz[Pair_candidates[vtxind][1]],
                   MASS_E);
      }
      else if (_fetchMuons) { // Looks at dimuons
	// Not enough muons candidates on the vertex
	if (_vtx->MatchedMuons.size()<2) continue;
        minDist = 999.;
        for(UInt_t i=0; i<_vtx->MatchedMuons.size(); i++) {
          lep1 = _vtx->MatchedMuons[i];

	  TVector3 vtxlep1(MuonCand_vtxx[lep1], MuonCand_vtxy[lep1], MuonCand_vtxz[lep1]);
	  MuonCand_dz[lep1] = _vtx->dZ(vtxlep1, lep1);
	  MuonCand_istight[lep1]&= (MuonCand_dz[lep1]<.5);

          for(UInt_t j=i+1; j<_vtx->MatchedMuons.size(); j++) {
            lep2 = _vtx->MatchedMuons[j];

	    TVector3 vtxlep2(MuonCand_vtxx[lep2], MuonCand_vtxy[lep2], MuonCand_vtxz[lep2]);
	    MuonCand_dz[lep2] = _vtx->dZ(vtxlep2, lep2);
	    MuonCand_istight[lep2]&= (MuonCand_dz[lep2]<.5);

            if (MuonCand_charge[lep1]*MuonCand_charge[lep2]>0) {
              continue;
            }
            foundPairOnVertex = true;
            leptonsDist=sqrt(pow(MuonCand_vtxx[lep1]-MuonCand_vtxx[lep2],2)+
			     pow(MuonCand_vtxy[lep1]-MuonCand_vtxy[lep2],2)+
			     pow(MuonCand_vtxz[lep1]-MuonCand_vtxz[lep2],2));
            if (leptonsDist<minDist) {
              minDist = leptonsDist;
              Pair_candidates[vtxind][0] = lep1;
              Pair_candidates[vtxind][1] = lep2;
            }
          }
        }
        if (Pair_candidates[vtxind][0]==-1 || Pair_candidates[vtxind][1]==-1) {
          foundPairOnVertex = false;
        }
        l1.SetXYZM(MuonCand_px[Pair_candidates[vtxind][0]],
                   MuonCand_py[Pair_candidates[vtxind][0]],
                   MuonCand_pz[Pair_candidates[vtxind][0]],
                   MASS_MU);
        l2.SetXYZM(MuonCand_px[Pair_candidates[vtxind][1]],
                   MuonCand_py[Pair_candidates[vtxind][1]],
                   MuonCand_pz[Pair_candidates[vtxind][1]],
                   MASS_MU);
      }

#ifdef DEBUG
      std::cout << "=====> did we find a pair on this vertex ? " << foundPairOnVertex << std::endl;
#endif
      if (foundPairOnVertex) {
	Pair_mindist[vtxind] = minDist;
#ifdef DEBUG
      	std::cout << "Matched muons : " << std::endl;
      	for (i=0; i<_vtx->MatchedMuons.size(); i++) {
      	  std::cout << "-> " << _vtx->MatchedMuons[i] << std::endl;
      	}
      	std::cout << "Matched electrons : " << std::endl;
      	for (i=0; i<_vtx->MatchedElectrons.size(); i++) {
      	  std::cout << "-> " << _vtx->MatchedElectrons[i] << std::endl;
        }
#endif
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
        
        for (Int_t j=0; j<nPFPhotonCand; j++) {
          pairgmass = sqrt(std::pow(l1.P() +l2.P() +PFPhotonCand_p[j], 2)-
			   std::pow(l1.Px()+l2.Px()+PFPhotonCand_px[j], 2)-
			   std::pow(l1.Py()+l2.Py()+PFPhotonCand_py[j], 2)-
			   std::pow(l1.Pz()+l2.Pz()+PFPhotonCand_pz[j], 2));
          //std::cout << "Photon " << j << " added to give a mass = " << pairgmass << std::endl;
          PairGamma_mass[vtxind][j] = pairgmass;
        }
        PrimVertexCand_hasdil[vtxind] = 1;
        
	nCandidatesInEvent++;
	nCandidates++;
        foundPairInEvent = true;
      }
      
      delete _vtx;
    }
    nFilteredPrimVertexCand = vtxind;
    //std::cout << "end of the loop on the vertices : vtxind = " << vtxind << std::endl;
    delete _leptonType;
  }
  nExtraTracks = etind;
  
  //std::cout << "Passed Loop on vertices" << std::endl;
  
  // Get the Jet collection from the event
  // PAT
  iEvent.getByLabel(jetLabel_, jetColl);
  for (jet=jetColl->begin(); jet!=jetColl->end() && nJetCand<MAX_JETS; jet++) {
    JetCand_e[nJetCand] = jet->energy();
    JetCand_px[nJetCand] = jet->px();
    JetCand_py[nJetCand] = jet->py();
    JetCand_pz[nJetCand] = jet->pz();
    JetCand_phi[nJetCand] = jet->phi();
    JetCand_eta[nJetCand] = jet->eta();
    totalJetEnergy += JetCand_e[nJetCand];
    // Find kinematics quantities associated to the highest energy jet
    if(JetCand_e[nJetCand]>HEJet_e) {
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

  //std::cout << "Passed Loop on jets" << std::endl;

  // Missing ET
  iEvent.getByLabel(metLabel_, MET); 
  const reco::PFMETCollection* metColl = MET.product(); 
  met = metColl->begin();
  
  Etmiss = met->et();
  Etmiss_phi = met->phi();
  Etmiss_x = met->px();
  Etmiss_y = met->py();
  Etmiss_z = met->pz();
  Etmiss_significance = met->significance();

  //std::cout << "Passed MET" << std::endl;

  if (foundPairInEvent) {
    if (printCandidates_) {
      std::cout << "Event " << Run << ":" << EventNum << " has " << nCandidatesInEvent << " leptons pair(s) candidate(s) (vertex mult. : " << nPrimVertexCand << ")" << std::endl;
    }
    tree_->Fill();
  }
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
  /*tree_->Branch("AvgInstDelLumi", &AvgInstDelLumi, "AvgInstDelLumi/D");
  tree_->Branch("BunchInstLumi", &BunchInstLumi, "BunchInstLumi[3]/D");*/

  tree_->Branch("nHLT", &nHLT, "nHLT/I");
  tree_->Branch("HLT_Accept", HLT_Accept, "HLT_Prescl[nHLT]/I");
  tree_->Branch("HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I");
  
  if (_fetchMuons) {
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
  
  if (_fetchElectrons) {
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
    tree_->Branch("EleCand_wp80", EleCand_wp80, "EleCand_wp80[nEleCand]/I");   
    tree_->Branch("EleCand_mediumID", EleCand_mediumID, "EleCand_mediumID[nEleCand]/I");    
    tree_->Branch("EleCand_looseID", EleCand_looseID, "EleCand_looseID[nEleCand]/I");   
    /*tree_->Branch("EleCand_looseid", EleCand_looseid,"EleCand_looseid[nEleCand]/I");
    tree_->Branch("EleCand_likelihoodid", EleCand_likelihoodid,"EleCand_likelihoodid[nEleCand]/D");
    tree_->Branch("EleCand_robustid", EleCand_robustid,"EleCand_robustid[nEleCand]/I");*/
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
  tree_->Branch("nPFPhotonCand", &nPFPhotonCand, "nPFPhotonCand/I");
  tree_->Branch("PFPhotonCand_pt", PFPhotonCand_pt, "PFPhotonCand_pt[nPFPhotonCand]/D");
  tree_->Branch("PFPhotonCand_eta", PFPhotonCand_eta, "PFPhotonCand_eta[nPFPhotonCand]/D");
  tree_->Branch("PFPhotonCand_phi", PFPhotonCand_phi, "PFPhotonCand_phi[nPFPhotonCand]/D");
  tree_->Branch("PFPhotonCand_drtrue", PFPhotonCand_drtrue, "PFPhotonCand_drtrue[nPFPhotonCand]/D"); 
  tree_->Branch("PFPhotonCand_detatrue", PFPhotonCand_detatrue, "PFPhotonCand_detatrue[nPFPhotonCand]/D"); 
  tree_->Branch("PFPhotonCand_dphitrue", PFPhotonCand_dphitrue, "PFPhotonCand_dphitrue[nPFPhotonCand]/D"); 
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
  tree_->Branch("nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nPrimVertexCand/I");
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
  tree_->Branch("PairGamma_mass", PairGamma_mass, "PairGamma_mass[nPrimVertexCand][nPFPhotonCand]/D");
  if (runOnMC_) {
    tree_->Branch("GenPair_p", &GenPair_p, "GenPair_p/D");
    tree_->Branch("GenPair_pt", &GenPair_pt, "GenPair_pt/D");
    tree_->Branch("GenPair_dpt", &GenPair_dpt, "GenPair_dpt/D");
    tree_->Branch("GenPair_mass", &GenPair_mass, "GenPair_mass/D");
    tree_->Branch("GenPair_eta", &GenPair_eta, "GenPair_eta/D");
    tree_->Branch("GenPair_phi", &GenPair_phi, "GenPair_phi/D");
    tree_->Branch("GenPair_dphi", &GenPair_dphi, "GenPair_dphi/D");
    tree_->Branch("GenPair_3Dangle", &GenPair_3Dangle, "GenPair_3Dangle[nPrimVertexCand]/D");
    tree_->Branch("HPS_acc420b1", &HPS_acc420b1, "HPS_acc420b1/D");
    tree_->Branch("HPS_acc220b1", &HPS_acc220b1, "HPS_acc220b1/D");
    tree_->Branch("HPS_acc420and220b1", &HPS_acc420and220b1, "HPS_acc420and220b1/D");
    tree_->Branch("HPS_acc420or220b1", &HPS_acc420or220b1, "HPS_acc420or220b1/D");
    tree_->Branch("HPS_acc420b2", &HPS_acc420b2, "HPS_acc420b2/D");
    tree_->Branch("HPS_acc220b2", &HPS_acc220b2, "HPS_acc220b2/D");
    tree_->Branch("HPS_acc420and220b2", &HPS_acc420and220b2, "HPS_acc420and220b2/D");
    tree_->Branch("HPS_acc420or220b2", &HPS_acc420or220b2, "HPS_acc420or220b2/D");
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
  tree_->Branch("nTruePUforPUWeight",&nTruePUforPUWeight,"nTruePUforPUWeight/I");
  tree_->Branch("nTruePUafterPUWeight",&nTruePUafterPUWeight,"nTruePUafterPUWeight/D");
  tree_->Branch("nTruePUforPUWeightBXM1", &nTruePUforPUWeightBXM1, "nTruePUforPUWeightBXM1/I");
  tree_->Branch("nTruePUafterPUWeightBXM1", &nTruePUafterPUWeightBXM1, "nTruePUafterPUWeightBXM1/D");
  tree_->Branch("nTruePUforPUWeightBXP1", &nTruePUforPUWeightBXP1, "nTruePUforPUWeightBXP1/I"); 
  tree_->Branch("nTruePUafterPUWeightBXP1", &nTruePUafterPUWeightBXP1, "nTruePUafterPUWeightBXP1/D"); 
  tree_->Branch("nTruePUforPUWeightBX0", &nTruePUforPUWeightBX0, "nTruePUforPUWeightBX0/I");
  tree_->Branch("nTruePUafterPUWeightBX0", &nTruePUafterPUWeightBX0, "nTruePUafterPUWeightBX0/D");
  tree_->Branch("Weight", &Weight, "Weight/D");
  tree_->Branch("PUWeightTrue", &PUWeightTrue, "PUWeightTrue/D");

  nCandidates = 0;
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
  bool changed;
  std::string triggerName_;
  
  changed = true;
  if (hltConfig_.init(iRun, iSetup, hltMenuLabel_, changed)) {
    if (changed) {
      // check if trigger name in (new) config
      triggerName_ = "HLT_DoubleMu0";
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
        const UInt_t n = hltConfig_.size();
        const UInt_t triggerIndex = hltConfig_.triggerIndex(triggerName_);
        if (triggerIndex>=n) {
          std::cout << "GammaGammaMuMu::analyze:"
                    << " TriggerName " << triggerName_ 
                    << " not available in (new) config!"
		    << std::endl;
          //std::cout << "Available TriggerNames are: " << std::endl;
          //hltConfig_.dump("Triggers");
        }
      }
      //hltConfig_.dump("Streams");
      //hltConfig_.dump("Datasets");
      //hltConfig_.dump("PrescaleTable");
      //hltConfig_.dump("ProcessPSet");
    }
  }
  else {
    std::cout << "GammaGammaMuMu::beginRun:"
              << " config extraction failure with process name "
              << hltMenuLabel_
	      << std::endl;
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

/// Class PrimaryVertex

//
// constructors and destructor
//
PrimaryVertex::PrimaryVertex(std::vector<std::string>& _leptonsType, std::map<Int_t,TLorentzVector>& _muonsMomenta, std::map<Int_t,TLorentzVector>& _electronsMomenta) :
  nTracks(0),
  nMatchedTracks(0),
  nUnmatchedTracks(0),
  nMatchedMuons(0),
  nMatchedElectrons(0),
  FetchMuons(false),
  FetchElectrons(false)
{
  //LeptonsType = _leptonsType;
  for (i=0; i<_leptonsType.size(); i++) {
    if (_leptonsType[i]=="Muon") {
      FetchMuons = true;
    }
    else if (_leptonsType[i]=="Electron") {
      FetchElectrons = true;
    }
  }
  MuonMomenta = _muonsMomenta;
  ElectronMomenta = _electronsMomenta;
}

PrimaryVertex::~PrimaryVertex() {
}

void
PrimaryVertex::SetPosition(Double_t _x, Double_t _y, Double_t _z)
{
  Position.SetXYZ(_x, _y, _z);
#ifdef DEBUG
  std::cout << "[PrimaryVertex] SetPosition :: Vertex located at (" << Position.x() << ", " << Position.y() << ", " << Position.z() << ")" << std::endl;
#endif
}

/**
 * \brief Matches a track arising from the vertex with a lepton track from the
 *  internal collection
 */
Int_t
PrimaryVertex::AddTrack(const reco::TrackRef& _track, TString& _leptonType)
{
  nTracks++; // total number of tracks matched with the vertex
  std::map<Int_t,TLorentzVector>::iterator lep;
  for (lep=MuonMomenta.begin(); lep!=MuonMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "muon";
    MatchedMuons.push_back(lep->first);
    nMatchedMuons++;
    nMatchedTracks++;
    return lep->first;
  }
  for (lep=ElectronMomenta.begin(); lep!=ElectronMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "electron";
    MatchedElectrons.push_back(lep->first);
    nMatchedElectrons++;
    nMatchedTracks++;
    return lep->first;
  }
  nUnmatchedTracks++;
  return -1;
}

Double_t
PrimaryVertex::dZ(TVector3 _vmu, Int_t _muind)
{
  TLorentzVector m(MuonMomenta[_muind]);
  return (_vmu.Z()-Position.Z())-((_vmu.X()-Position.X())*m.Px()+(_vmu.Y()-Position.Y())*m.Py())/m.Pt()*m.Pz()/m.Pt();
}

HLTmatches::HLTmatches(std::vector<std::string> _HLTlist)
{
  for (i=0; i<_HLTlist.size(); i++) {
    HLTnames.push_back(_HLTlist[i].substr(0, _HLTlist[i].find_first_of("*")));
  }
#ifdef DEBUG
  for (i=0; i<HLTnames.size(); i++) {
    std::cout << i << " ==> " << HLTnames[i] << std::endl;
  }
#endif
}

HLTmatches::~HLTmatches()
{
}

Int_t
HLTmatches::TriggerNum(std::string _trigName)
{
  for (i=0; i<HLTnames.size(); i++) {
    if (_trigName.find(HLTnames[i])!=std::string::npos) return i;
  }
  return -1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaGammaLL);
