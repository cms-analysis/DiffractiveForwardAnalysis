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
// $Id: GammaGammaLL.cc,v 1.1 2013/03/20 17:56:13 lforthom Exp $
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
  outPileupFile_ = iConfig.getUntrackedParameter<std::string>("outpufilename", "test.root");
  
  // Leptons input tags
  leptonsType_ = iConfig.getParameter< std::vector<std::string> >("LeptonsType");
  for (i=0; i<leptonsType_.size(); i++) {
    if (leptonsType_[i]=="Muon") _fetchMuons = true;
    else if (leptonsType_[i]=="Electron") _fetchElectrons = true;
  }
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalMuonCollectionLabel", std::string("muons"));
  eleLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalEleCollectionLabel", std::string("gsfElectrons"));
  isoValInputTag_ = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags"); 
  beamSpotInputTag_ = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  conversionsInputTag_ = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
  rhoIsoInputTag_ = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag"); 
  
  file = new TFile(outputFile_.c_str(), "recreate");
  file->cd();
  // tree definition
  tree = new TTree("ntp1", "ntp1");

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
    
    LumiWeights = new edm::Lumi3DReWeighting(mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_, outPileupFile_);
    LumiWeights->weight3D_init(1.0); 
  }

}


GammaGammaLL::~GammaGammaLL()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  const edm::ParameterSet& pset = edm::getProcessParameterSet();
  TList *list = tree->GetUserInfo();
  list->Add(new TObjString(pset.dump().c_str()));
  file->Write();
  file->Close();

	delete _hlts;
  delete tree;

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

  // First initialization of the variables
  nPrimVertexCand = nFilteredPrimVertexCand = 0;
  nMuonCand = nEleCand = nLeptonCand = 0;
  nExtraTracks = nQualityExtraTrack = 0;
  nGenMuonCand = nGenMuonCandOutOfAccept = 0;
  nGenEleCand = nGenEleCandOutOfAccept = 0;
  
  HPS_acc420b1 = HPS_acc220b1 = HPS_acc420and220b1 = HPS_acc420or220b1 = -1;
  HPS_acc420b2 = HPS_acc220b2 = HPS_acc420and220b2 = HPS_acc420or220b2 = -1;
  GenPair_p = GenPair_pt = GenPair_mass = GenPair_phi = GenPair_eta = -999.;
  GenPair_dphi = GenPair_dpt = GenPair_3Dangle = 0.;
  for (i=0; i<MAX_MUONS; i++) {
    MuonCand_p[i] = MuonCand_px[i] = MuonCand_py[i] = MuonCand_pz[i] = -999.;
    MuonCand_pt[i] = MuonCand_eta[i] = MuonCand_phi[i] = -999.;
    MuonCand_charge[i] = -999;
    MuonCand_vtxx[i] = MuonCand_vtxy[i] = MuonCand_vtxz[i] = -999.;
    MuonCand_isglobal[i] = MuonCand_istracker[i] = MuonCand_isstandalone[i] = -999;
  }
  for (i=0; i<MAX_ELE; i++) {
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
  for (i=0; i<MAX_ET; i++) {
  	ExtraTrack_p[i] = ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
  	ExtraTrack_pt[i] = ExtraTrack_eta[i] = ExtraTrack_phi[i] = -999.;
  	ExtraTrack_charge[etind] = -999;
  	ExtraTrack_chi2[i] = ExtraTrack_ndof[i] = ExtraTrack_vtxdxyz[i] = -999.;
    ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
    ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
  }
  for (i=0; i<MAX_VTX; i++) {
    Pair_candidates[i][0] = Pair_candidates[vtxind][1] = PrimVertexCand_id[i] = -1;
    PrimVertexCand_tracks[i] = PrimVertexCand_matchedtracks[i] = PrimVertexCand_unmatchedtracks[i] = 0;
    Pair_mindist[i] = Pair_p[i] = Pair_pt[i] = Pair_mass[i] = Pair_phi[i] = Pair_eta[i] = -999.;
    Pair_dphi[i] = Pair_dpt[i] = Pair_3Dangle[i] = -999.;
    PrimVertexCand_x[i] = PrimVertexCand_y[vtxind] = PrimVertexCand_z[vtxind] = -999.;
    PrimVertexCand_chi2[i] = PrimVertexCand_ndof[i] = -999.;
  }
  
  Weight3D = 1.;
  
  foundPairInEvent = false;
  foundGenCandPairInEvent = false;
  
  muonsMomenta.clear();
  electronsMomenta.clear();
  
  _leptonptmp = new TLorentzVector();
  
  // Run and BX information
  BX = iEvent.bunchCrossing();
  Run = iEvent.id().run();
  LumiSection = iEvent.luminosityBlock();
  EventNum = iEvent.id().event();
  
  // High level trigger information retrieval  
  LookAtTriggers(iEvent, iSetup);
  
  // beam spot information
  iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // Get the vertex collection from the event
  iEvent.getByLabel(recoVertexLabel_, recoVertexColl);
  const reco::VertexCollection* vertices = recoVertexColl.product();

  // rho for isolation 
  iEvent.getByLabel(rhoIsoInputTag_, rhoIso_h); 
  rhoIso = *(rhoIso_h.product()); 

  // Generator level information
  if (runOnMC_) {
    iEvent.getByLabel(genLabel_, genPartColl);
    
    for (genPart=genPartColl->begin(); genPart!=genPartColl->end(); genPart++) {
      if (genPart->pt()<minPtMC_ || (minEtaMC_!=-1. && fabs(genPart->eta())>minEtaMC_)) {
        if (fabs(genPart->pdgId())==13) nGenMuonCandOutOfAccept += 1;
        if (fabs(genPart->pdgId())==11) nGenEleCandOutOfAccept += 1;
        continue;
      }
      if (genPart->pdgId()!=2212 && genPart->status()!=1) continue;
      if (fabs(genPart->pdgId())==13 && nGenMuonCand<MAX_GENMU) {
        GenMuonCand_p[nGenMuonCand] = genPart->p();
        GenMuonCand_px[nGenMuonCand] = genPart->px();
        GenMuonCand_py[nGenMuonCand] = genPart->py();
        GenMuonCand_pz[nGenMuonCand] = genPart->pz();
        GenMuonCand_pt[nGenMuonCand] = genPart->pt();
        GenMuonCand_eta[nGenMuonCand] = genPart->eta();
        GenMuonCand_phi[nGenMuonCand] = genPart->phi();
        
        nGenMuonCand += 1;
      }
      if (fabs(genPart->pdgId())==11 && nGenEleCand<MAX_GENELE) {
        GenEleCand_p[nGenEleCand] = genPart->p();
        GenEleCand_px[nGenEleCand] = genPart->px();
        GenEleCand_py[nGenEleCand] = genPart->py();
        GenEleCand_pz[nGenEleCand] = genPart->pz();
        GenEleCand_pt[nGenEleCand] = genPart->pt();
        GenEleCand_eta[nGenEleCand] = genPart->eta();
        GenEleCand_phi[nGenEleCand] = genPart->phi();
        
        nGenEleCand += 1;
      }
      foundGenCandPairInEvent = true;
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
    Weight3D = LumiWeights->weight3D(*iEventB);
  }

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
      
      _leptonptmp->SetXYZM(muon->px(), muon->py(), muon->pz(), muon->mass());
      muonsMomenta.insert(std::pair<Int_t,TLorentzVector>(nMuonCand, *_leptonptmp));
      
      MuonCand_vtxx[nMuonCand] = muon->vertex().x();
      MuonCand_vtxy[nMuonCand] = muon->vertex().y();
      MuonCand_vtxz[nMuonCand] = muon->vertex().z();
      
      MuonCand_isglobal[nMuonCand] = muon->isGlobalMuon();
      MuonCand_istracker[nMuonCand] = muon->isTrackerMuon();
      MuonCand_isstandalone[nMuonCand] = muon->isStandAloneMuon();

      nMuonCand += 1;
    }
  }
  
  // Get the electrons collection from the event
  if (_fetchElectrons) {
    iEvent.getByLabel(eleLabel_, eleColl);
    // RECO electrons
    iEvent.getByLabel(InputTag("gsfElectrons"), eleCollRECO);
    // iso deposits
    IsoDepositVals isoVals(isoValInputTag_.size());
    // New 2012 electron ID variables conversions
    iEvent.getByLabel(conversionsInputTag_, conversions_h);

    for (size_t j = 0; j<isoValInputTag_.size(); j++) { 
      iEvent.getByLabel(isoValInputTag_[j], isoVals[j]); 
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
  	  
      _leptonptmp->SetXYZM(electron->px(), electron->py(), electron->pz(), electron->mass());
      electronsMomenta.insert(std::pair<Int_t,TLorentzVector>(nEleCand, *_leptonptmp));

      EleCand_vtxx[nEleCand] = electron->vertex().x();
      EleCand_vtxy[nEleCand] = electron->vertex().y();
      EleCand_vtxz[nEleCand] = electron->vertex().z();
  	  
      if(electron->closestCtfTrackRef().isNonnull()) { // Only for PAT::Electron
        EleCandTrack_p[nEleCand] = electron->closestCtfTrackRef()->p(); 
        EleCandTrack_pt[nEleCand] = electron->closestCtfTrackRef()->pt();  
        EleCandTrack_eta[nEleCand] = electron->closestCtfTrackRef()->eta();  
        EleCandTrack_phi[nEleCand] = electron->closestCtfTrackRef()->phi();  
        EleCandTrack_vtxz[nEleCand] = electron->closestCtfTrackRef()->vertex().z();
      }
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
        selectmedium = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, recoVertexColl, iso_ch, iso_em, iso_nh, rhoIso);
        selectloose = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele, conversions_h, beamSpot, recoVertexColl, iso_ch, iso_em, iso_nh, rhoIso);
      }
      if(selectmedium) EleCand_mediumID[nEleCand] = 1;
      if(selectloose) EleCand_looseID[nEleCand] = 1;
      
  	  nEleCand++;
    }
  }
  
  delete _leptonptmp;
  
  nLeptonCand += (_fetchMuons) ? nMuonCand : 0;
  nLeptonCand += (_fetchElectrons) ? nEleCand : 0;
	vtxind = 0;

  if (nLeptonCand>=2) {
    // Enough leptons candidates to go deeper and analyze the primary vertices
        
	  _leptonType = new TString();

		nPrimVertexCand = vertices->size();
		
    for (vertex=vertices->begin(); vertex!=vertices->end() && vtxind<MAX_VTX; vertex++) {
      nLeptonsInPrimVertex = 0;
      nExtraTracks = 0;
      nQualityExtraTrack = 0;
      foundPairOnVertex = false;
      etind = 0;
      
      PrimaryVertex* _vtx = new PrimaryVertex(leptonsType_, muonsMomenta, electronsMomenta);
      _vtx->SetPosition(vertex->x(), vertex->y(), vertex->z());

      // Loop on all the tracks matched with this vertex    
      reco::Vertex::trackRef_iterator _track;
    	for (_track=vertex->tracks_begin(); _track!=vertex->tracks_end() && nExtraTracks<MAX_ET; _track++) {
    	  Int_t leptonId_ = _vtx->AddTrack((*_track).castTo<reco::TrackRef>(), *_leptonType);
        if (leptonId_==-1) { // Track was not matched to any of the leptons in the collection
        	ExtraTrack_purity[etind] = (*_track)->quality(reco::TrackBase::highPurity);
        	ExtraTrack_nhits[etind] = (*_track)->numberOfValidHits();
        	if (ExtraTrack_purity[etind]==1 && ExtraTrack_nhits[etind]>=3) {
        		nQualityExtraTrack += 1;
        	}
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
          ExtraTrack_vtxdxyz[etind] = sqrt(std::pow(((*_track)->vertex().x()-_vtx->Position.X()),2)
                                         + std::pow(((*_track)->vertex().y()-_vtx->Position.Y()),2)
                                         + std::pow(((*_track)->vertex().z()-_vtx->Position.Z()),2));
          ExtraTrack_vtxT[etind] = sqrt(std::pow((*_track)->vertex().x()-_vtx->Position.X(),2)
                                      + std::pow((*_track)->vertex().y()-_vtx->Position.Y(),2));
          ExtraTrack_vtxZ[etind] = fabs((*_track)->vertex().z()-_vtx->Position.Z());
          ExtraTrack_x[etind] = (*_track)->vertex().x();
          ExtraTrack_y[etind] = (*_track)->vertex().y();
          ExtraTrack_z[etind] = (*_track)->vertex().z();
        	etind += 1;
        }
        nExtraTracks = etind;
        nLeptonsInPrimVertex += 1;
    	}
    	if (nLeptonsInPrimVertex<2) continue;
    	
    	// At this stage we have at least two matched leptons track on the vertex
      Pair_candidates[vtxind][0] = -1;
      Pair_candidates[vtxind][1] = -1;

    	if (_fetchElectrons && _fetchMuons) { // Looks at electron+muon
    	  // Not enough muons or electrons candidates on the vertex
    	  if (_vtx->MatchedElectrons.size()==0 or _vtx->MatchedMuons.size()==0) {
    	    continue;
    	  }
        minDist = 999.;
        for(UInt_t i=0; i<_vtx->MatchedMuons.size(); i++) {
          lep1 = _vtx->MatchedMuons[i];
          for(UInt_t j=0; j<_vtx->MatchedElectrons.size(); j++) {
            lep2 = _vtx->MatchedElectrons[j];
            if (MuonCand_charge[lep1]*EleCand_charge[lep2]>0) {
              continue;
            }
            foundPairOnVertex = true;
            leptonsDist=sqrt(pow(MuonCand_vtxx[lep1]-EleCand_vtxx[lep2],2)
                            +pow(MuonCand_vtxy[lep1]-EleCand_vtxy[lep2],2)
                            +pow(MuonCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
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
            leptonsDist=sqrt(pow(EleCand_vtxx[lep1]-EleCand_vtxx[lep2],2)
                            +pow(EleCand_vtxy[lep1]-EleCand_vtxy[lep2],2)
                            +pow(EleCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
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
          for(UInt_t j=i+1; j<_vtx->MatchedMuons.size(); j++) {
            lep2 = _vtx->MatchedMuons[j];
            if (MuonCand_charge[lep1]*MuonCand_charge[lep2]>0) {
              continue;
            }
            foundPairOnVertex = true;
            leptonsDist=sqrt(pow(MuonCand_vtxx[lep1]-MuonCand_vtxx[lep2],2)
                            +pow(MuonCand_vtxy[lep1]-MuonCand_vtxy[lep2],2)
                            +pow(MuonCand_vtxz[lep1]-MuonCand_vtxz[lep2],2));
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

      if (foundPairOnVertex) {
        std::cout << "Pair : ("
                  << Pair_candidates[vtxind][0] << ", "
                  << Pair_candidates[vtxind][1]
                  << ") at distance " << minDist
                  << std::endl;
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
        
        PrimVertexCand_id[vtxind] = vtxind;
        PrimVertexCand_x[vtxind] = _vtx->Position.X();
        PrimVertexCand_y[vtxind] = _vtx->Position.Y();
        PrimVertexCand_z[vtxind] = _vtx->Position.Z();
				PrimVertexCand_tracks[vtxind] = _vtx->nTracks;
        PrimVertexCand_matchedtracks[vtxind] = _vtx->nMatchedTracks;
        PrimVertexCand_unmatchedtracks[vtxind] = _vtx->nUnmatchedTracks;
        PrimVertexCand_chi2[vtxind] = vertex->chi2();
        PrimVertexCand_ndof[vtxind] = vertex->ndof();
        
		    nCandidates += 1;
        foundPairInEvent = true;
      }
      
      vtxind += 1;

      delete _vtx;
    }
    nFilteredPrimVertexCand = vtxind;
    delete _leptonType;
  }
  
  if (foundPairInEvent) {
    tree->Fill();
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaLL::beginJob()
{
  // Booking the ntuple

  tree->Branch("Run", &Run, "Run/I");
  tree->Branch("LumiSection", &LumiSection, "LumiSection/I");
  tree->Branch("BX", &BX, "BX/I");
  tree->Branch("EventNum", &EventNum, "EventNum/I");
  /*tree->Branch("AvgInstDelLumi", &AvgInstDelLumi, "AvgInstDelLumi/D");
  tree->Branch("BunchInstLumi", &BunchInstLumi, "BunchInstLumi[3]/D");*/

  tree->Branch("nHLT", &nHLT, "nHLT/I");
  tree->Branch("HLT_Accept", HLT_Accept, "HLT_Prescl[nHLT]/I");
  tree->Branch("HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I");
  
  if (_fetchMuons) {
    tree->Branch("nMuonCand", &nMuonCand, "nMuonCand/I");
    tree->Branch("MuonCand_px", MuonCand_px, "MuonCand_px[nMuonCand]/D");
    tree->Branch("MuonCand_py", MuonCand_py, "MuonCand_py[nMuonCand]/D");
    tree->Branch("MuonCand_pz", MuonCand_pz, "MuonCand_pz[nMuonCand]/D");
    tree->Branch("MuonCand_p", MuonCand_p, "MuonCand_p[nMuonCand]/D");
    tree->Branch("MuonCand_pt", MuonCand_pt, "MuonCand_pt[nMuonCand]/D");
    tree->Branch("MuonCand_eta", MuonCand_eta, "MuonCand_eta[nMuonCand]/D");
    tree->Branch("MuonCand_phi", MuonCand_phi, "MuonCand_phi[nMuonCand]/D");
    tree->Branch("MuonCand_charge", MuonCand_charge, "MuonCand_charge[nMuonCand]/I");
    tree->Branch("MuonCand_vtxx", MuonCand_vtxx, "MuonCand_vtxx[nMuonCand]/D");
    tree->Branch("MuonCand_vtxy", MuonCand_vtxy, "MuonCand_vtxy[nMuonCand]/D");
    tree->Branch("MuonCand_vtxz", MuonCand_vtxz, "MuonCand_vtxz[nMuonCand]/D");
    tree->Branch("nGenMuonCand", &nGenMuonCand, "nGenMuonCand/I");
    tree->Branch("nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept, "nGenMuonCandOutOfAccept/I");    
    tree->Branch("GenMuonCand_p", GenMuonCand_p, "GenMuonCand_p[nGenMuonCand]/D");
    tree->Branch("GenMuonCand_px", GenMuonCand_px, "GenMuonCand_px[nGenMuonCand]/D");
    tree->Branch("GenMuonCand_py", GenMuonCand_py, "GenMuonCand_py[nGenMuonCand]/D");
    tree->Branch("GenMuonCand_pz", GenMuonCand_pz, "GenMuonCand_pz[nGenMuonCand]/D");
    tree->Branch("GenMuonCand_pt", GenMuonCand_pt, "GenMuonCand_pt[nGenMuonCand]/D");
    tree->Branch("GenMuonCand_eta", GenMuonCand_eta, "GenMuonCand_eta[nGenMuonCand]/D");
    tree->Branch("GenMuonCand_phi", GenMuonCand_phi, "GenMuonCand_phi[nGenMuonCand]/D");
  }
  
  if (_fetchElectrons) {
    tree->Branch("nEleCand", &nEleCand, "nEleCand/I");
    tree->Branch("EleCand_px", EleCand_px, "EleCand_px[nEleCand]/D");
    tree->Branch("EleCand_py", EleCand_py, "EleCand_py[nEleCand]/D");
    tree->Branch("EleCand_pz", EleCand_pz, "EleCand_pz[nEleCand]/D");
    tree->Branch("EleCand_p", EleCand_p, "EleCand_p[nEleCand]/D");
    tree->Branch("EleCand_e", EleCand_e, "EleCand_e[nEleCand]/D");
    tree->Branch("EleCand_et", EleCand_et, "EleCand_et[nEleCand]/D");
    tree->Branch("EleCand_eta", EleCand_eta, "EleCand_eta[nEleCand]/D");
    tree->Branch("EleCand_phi", EleCand_phi, "EleCand_phi[nEleCand]/D");
    tree->Branch("EleCand_charge", EleCand_charge, "EleCand_charge[nEleCand]/I");
    tree->Branch("EleCand_vtxx", EleCand_vtxx, "EleCand_vtxx[nEleCand]/D");
    tree->Branch("EleCand_vtxy", EleCand_vtxy, "EleCand_vtxy[nEleCand]/D");
    tree->Branch("EleCand_vtxz", EleCand_vtxz, "EleCand_vtxz[nEleCand]/D");
    tree->Branch("EleCandTrack_p", EleCandTrack_p, "EleCandTrack_p[nEleCand]/D");
    tree->Branch("EleCandTrack_pt", EleCandTrack_pt, "EleCandTrack_pt[nEleCand]/D");
    tree->Branch("EleCandTrack_eta", EleCandTrack_eta, "EleCandTrack_eta[nEleCand]/D");
    tree->Branch("EleCandTrack_phi", EleCandTrack_phi, "EleCandTrack_phi[nEleCand]/D");
    tree->Branch("EleCandTrack_vtxz", EleCandTrack_vtxz, "EleCandTrack_vtxz[nEleCand]/D");
    tree->Branch("EleCand_deltaPhi", EleCand_deltaPhi, "EleCand_deltaPhi[nEleCand]/D"); 
    tree->Branch("EleCand_deltaEta", EleCand_deltaEta, "EleCand_deltaEta[nEleCand]/D"); 
    tree->Branch("EleCand_HoverE", EleCand_HoverE, "EleCand_HoverE[nEleCand]/D"); 
    tree->Branch("EleCand_trackiso", EleCand_trackiso, "EleCand_trackiso[nEleCand]/D"); 
    tree->Branch("EleCand_ecaliso", EleCand_ecaliso," EleCand_ecaliso[nEleCand]/D"); 
    tree->Branch("EleCand_hcaliso", EleCand_hcaliso," EleCand_hcaliso[nEleCand]/D"); 
    tree->Branch("EleCand_sigmaIetaIeta", EleCand_sigmaIetaIeta, "EleCand_sigmaIetaIeta[nEleCand]/D"); 
    tree->Branch("EleCand_convDist", EleCand_convDist, "EleCand_convDist[nEleCand]/D"); 
    tree->Branch("EleCand_convDcot", EleCand_convDcot, "EleCand_convDcot[nEleCand]/D"); 
    tree->Branch("EleCand_ecalDriven", EleCand_ecalDriven, "EleCand_ecalDriven[nEleCand]/D");  
    tree->Branch("EleCand_wp80", EleCand_wp80, "EleCand_wp80[nEleCand]/I");   
    tree->Branch("EleCand_mediumID", EleCand_mediumID, "EleCand_mediumID[nEleCand]/I");    
    tree->Branch("EleCand_looseID", EleCand_looseID, "EleCand_looseID[nEleCand]/I");   
    /*tree->Branch("EleCand_looseid", EleCand_looseid,"EleCand_looseid[nEleCand]/I");
    tree->Branch("EleCand_likelihoodid", EleCand_likelihoodid,"EleCand_likelihoodid[nEleCand]/D");
    tree->Branch("EleCand_robustid", EleCand_robustid,"EleCand_robustid[nEleCand]/I");*/
    tree->Branch("nGenEleCand", &nGenEleCand, "nGenEleCand/I");
    tree->Branch("nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept, "nGenEleCandOutOfAccept/I");    
    tree->Branch("GenEleCand_p", GenEleCand_p, "GenEleCand_p[nGenEleCand]/D");
    tree->Branch("GenEleCand_px", GenEleCand_px, "GenEleCand_px[nGenEleCand]/D");
    tree->Branch("GenEleCand_py", GenEleCand_py, "GenEleCand_py[nGenEleCand]/D");
    tree->Branch("GenEleCand_pz", GenEleCand_pz, "GenEleCand_pz[nGenEleCand]/D");
    tree->Branch("GenEleCand_pt", GenEleCand_pt, "GenEleCand_pt[nGenEleCand]/D");
    tree->Branch("GenEleCand_eta", GenEleCand_eta, "GenEleCand_eta[nGenEleCand]/D");
    tree->Branch("GenEleCand_phi", GenEleCand_phi, "GenEleCand_phi[nGenEleCand]/D");
  }
  
	// Primary vertices' information
  tree->Branch("nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/I");
  tree->Branch("PrimVertexCand_id", PrimVertexCand_id, "PrimVertexCand_id[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_matchedtracks", PrimVertexCand_matchedtracks, "PrimVertexCand_matchedtracks[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_unmatchedtracks", PrimVertexCand_unmatchedtracks, "PrimVertexCand_unmatchedtracks[nPrimVertexCand]/I");
  tree->Branch("PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D");
  tree->Branch("PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/I");
  tree->Branch("nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nFilteredPrimVertexCand/I");

  // Lepton pairs' information
  tree->Branch("Pair_candidates", Pair_candidates, "Pair_candidates[nPrimVertexCand][2]/I");
  tree->Branch("Pair_mindist", Pair_mindist, "Pair_mindist[nPrimVertexCand]/D");
  tree->Branch("Pair_p", Pair_p, "Pair_p[nPrimVertexCand]/D");
  tree->Branch("Pair_pt", Pair_pt, "Pair_pt[nPrimVertexCand]/D");
  tree->Branch("Pair_dpt", Pair_dpt, "Pair_dpt[nPrimVertexCand]/D");
  tree->Branch("Pair_mass", Pair_mass, "Pair_mass[nPrimVertexCand]/D");
  tree->Branch("Pair_eta", Pair_eta, "Pair_eta[nPrimVertexCand]/D");
  tree->Branch("Pair_phi", Pair_phi, "Pair_phi[nPrimVertexCand]/D");
  tree->Branch("Pair_dphi", Pair_dphi, "Pair_dphi[nPrimVertexCand]/D");
  tree->Branch("Pair_3Dangle", Pair_3Dangle, "Pair_3Dangle[nPrimVertexCand]/D");
  tree->Branch("GenPair_p", &GenPair_p, "GenPair_p/D");
  tree->Branch("GenPair_pt", &GenPair_pt, "GenPair_pt/D");
  tree->Branch("GenPair_dpt", &GenPair_dpt, "GenPair_dpt/D");
  tree->Branch("GenPair_mass", &GenPair_mass, "GenPair_mass/D");
  tree->Branch("GenPair_eta", &GenPair_eta, "GenPair_eta/D");
  tree->Branch("GenPair_phi", &GenPair_phi, "GenPair_phi/D");
  tree->Branch("GenPair_dphi", &GenPair_dphi, "GenPair_dphi/D");
  tree->Branch("GenPair_3Dangle", &GenPair_3Dangle, "GenPair_3Dangle[nPrimVertexCand]/D");
  tree->Branch("HPS_acc420b1", &HPS_acc420b1, "HPS_acc420b1/D");
  tree->Branch("HPS_acc220b1", &HPS_acc220b1, "HPS_acc220b1/D");
  tree->Branch("HPS_acc420and220b1", &HPS_acc420and220b1, "HPS_acc420and220b1/D");
  tree->Branch("HPS_acc420or220b1", &HPS_acc420or220b1, "HPS_acc420or220b1/D");
  tree->Branch("HPS_acc420b2", &HPS_acc420b2, "HPS_acc420b2/D");
  tree->Branch("HPS_acc220b2", &HPS_acc220b2, "HPS_acc220b2/D");
  tree->Branch("HPS_acc420and220b2", &HPS_acc420and220b2, "HPS_acc420and220b2/D");
  tree->Branch("HPS_acc420or220b2", &HPS_acc420or220b2, "HPS_acc420or220b2/D");

  // Extra tracks on vertex's information
  tree->Branch("nExtraTracks", &nExtraTracks, "nExtraTracks/I");
  tree->Branch("ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I");
  tree->Branch("ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/I");
  tree->Branch("ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I");
  tree->Branch("ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/I");
  tree->Branch("ExtraTrack_p", ExtraTrack_p, "ExtraTrack_p[nExtraTracks]/D");
  tree->Branch("ExtraTrack_pt", ExtraTrack_pt, "ExtraTrack_pt[nExtraTracks]/D");
  tree->Branch("ExtraTrack_px", ExtraTrack_px, "ExtraTrack_px[nExtraTracks]/D");
  tree->Branch("ExtraTrack_py", ExtraTrack_py, "ExtraTrack_py[nExtraTracks]/D");
  tree->Branch("ExtraTrack_pz", ExtraTrack_pz, "ExtraTrack_pz[nExtraTracks]/D");
  tree->Branch("ExtraTrack_eta", ExtraTrack_eta, "ExtraTrack_eta[nExtraTracks]/D");
  tree->Branch("ExtraTrack_phi", ExtraTrack_phi, "ExtraTrack_phi[nExtraTracks]/D");
  tree->Branch("ExtraTrack_chi2", ExtraTrack_chi2, "ExtraTrack_chi2[nExtraTracks]/D");
  tree->Branch("ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz, "ExtraTrack_vtxdxyz[nExtraTracks]/D");
  tree->Branch("ExtraTrack_vtxT", ExtraTrack_vtxT, "ExtraTrack_vtxT[nExtraTracks]/D");
  tree->Branch("ExtraTrack_vtxZ", ExtraTrack_vtxZ, "ExtraTrack_vtxZ[nExtraTracks]/D");
  tree->Branch("ExtraTrack_x", ExtraTrack_x, "ExtraTrack_x[nExtraTracks]/D");
  tree->Branch("ExtraTrack_y", ExtraTrack_y, "ExtraTrack_y[nExtraTracks]/D");
  tree->Branch("ExtraTrack_z", ExtraTrack_z, "ExtraTrack_z[nExtraTracks]/D");
  tree->Branch("nQualityExtraTrack", &nQualityExtraTrack, "nQualityExtraTrack/I");

  // Pileup reweighting
  tree->Branch("nTruePUforPUWeight",&nTruePUforPUWeight,"nTruePUforPUWeight/I");
  tree->Branch("nTruePUafterPUWeight",&nTruePUafterPUWeight,"nTruePUafterPUWeight/D");
  tree->Branch("nTruePUforPUWeightBXM1", &nTruePUforPUWeightBXM1, "nTruePUforPUWeightBXM1/I");
  tree->Branch("nTruePUafterPUWeightBXM1", &nTruePUafterPUWeightBXM1, "nTruePUafterPUWeightBXM1/D");
  tree->Branch("nTruePUforPUWeightBXP1", &nTruePUforPUWeightBXP1, "nTruePUforPUWeightBXP1/I"); 
  tree->Branch("nTruePUafterPUWeightBXP1", &nTruePUafterPUWeightBXP1, "nTruePUafterPUWeightBXP1/D"); 
  tree->Branch("nTruePUforPUWeightBX0", &nTruePUforPUWeightBX0, "nTruePUforPUWeightBX0/I");
  tree->Branch("nTruePUafterPUWeightBX0", &nTruePUafterPUWeightBX0, "nTruePUafterPUWeightBX0/D");
  tree->Branch("Weight3D", &Weight3D, "Weight3D/D");
  tree->Branch("PUWeightTrue", &PUWeightTrue, "PUWeightTrue/D");

  nCandidates = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaLL::endJob() 
{
	std::cout << "==> number of candidates in the dataset : " << nCandidates << std::endl;
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
                    << " not available in (new) config!" << std::endl;
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
              << hltMenuLabel_ << std::endl;
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
PrimaryVertex::SetPosition(Double_t _x, Double_t _y, Double_t _z) {
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
PrimaryVertex::AddTrack(const reco::TrackRef& _track, TString& _leptonType) {
	nTracks += 1; // total number of tracks matched with the vertex
  std::map<Int_t,TLorentzVector>::iterator lep;
  for (lep=MuonMomenta.begin(); lep!=MuonMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "muon";
    MatchedMuons.push_back(lep->first);
    nMatchedMuons += 1;
    nMatchedTracks += 1;
    return lep->first;
  }
  for (lep=ElectronMomenta.begin(); lep!=ElectronMomenta.end(); lep++) {
    if (fabs(_track->p()-lep->second.P())>.01) continue;
    if (fabs(_track->pt()-lep->second.Pt())>.01) continue;
    if (fabs(_track->eta()-lep->second.Eta())>.01) continue;
    if (fabs(_track->phi()-lep->second.Phi())>.01) continue;
    _leptonType = "electron";
    MatchedElectrons.push_back(lep->first);
    nMatchedElectrons += 1;
    nMatchedTracks += 1;
    return lep->first;
  }
  nUnmatchedTracks += 1;
  return -1;
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
