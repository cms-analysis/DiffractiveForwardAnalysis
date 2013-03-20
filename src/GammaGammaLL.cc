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
// $Id$
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
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outfilename", "output.root");
  leptonsType_ = iConfig.getParameter< std::vector<std::string> >("LeptonsType");

  _fetchMuons = false;
  _fetchElectrons = false;
  std::cout << "LeptonsType : " << std::endl;
  for (UInt_t i=0; i<leptonsType_.size(); i++) {
    if (leptonsType_[i]=="Muon") _fetchMuons = true;
    else if (leptonsType_[i]=="Electron") _fetchElectrons = true;
    std::cout << " -> " << leptonsType_[i] << std::endl;
  }

  recoVertexLabel_ = iConfig.getParameter<edm::InputTag>("RecoVertexLabel");
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalMuonCollectionLabel", std::string("muons"));
  
  // Electron input tags
  eleLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("GlobalEleCollectionLabel", std::string("gsfElectrons"));
  isoValInputTag_ = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags"); 
  beamSpotInputTag_ = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  conversionsInputTag_ = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
  rhoIsoInputTag_ = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag"); 
  
  file = new TFile(outputFile_.c_str(), "recreate");
  file->cd();
  // tree definition
  tree = new TTree("ntp1", "ntp1");

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

  delete tree;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GammaGammaLL::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // First initialization of the variables
  nPrimVertexCand = 0;
  nUnfilteredPrimVertexCand = 0;
  nMuonCand = 0;
  nEleCand = 0;
  nLeptonCand = 0;
  foundPair = false;
  muonsMomenta.clear();
  electronsMomenta.clear();
  
  _leptonptmp = new TLorentzVector();
  
  // beam spot information
  iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // Get the vertex collection from the event
  iEvent.getByLabel(recoVertexLabel_, recoVertexColl);
  const reco::VertexCollection* vertices = recoVertexColl.product();

  // rho for isolation 
  iEvent.getByLabel(rhoIsoInputTag_, rhoIso_h); 
  double rhoIso = *(rhoIso_h.product()); 

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
  	  
      EleCandTrack_p[nEleCand] = -999.;
      EleCandTrack_pt[nEleCand] = -999.; 
      EleCandTrack_eta[nEleCand] = -999.; 
      EleCandTrack_phi[nEleCand] = -999.; 
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
      
      double iso_ch =  (*(isoVals)[0])[ele];
      double iso_em = (*(isoVals)[1])[ele];
      double iso_nh = (*(isoVals)[2])[ele];

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
  //nLeptonCand = nMuonCand + nEleCand;
  //std::cout << "Number of leptons candidates : " << nLeptonCand << std::endl;

  if (nLeptonCand>=2) {
    // Enough leptons candidates to go deeper and analyze the primary vertices
        
	  _leptonType = new TString();

		nUnfilteredPrimVertexCand = vertices->size();
		
    for (vertex=vertices->begin(); vertex!=vertices->end(); vertex++) {
      nLeptonsInPrimVertex = 0;
      PrimaryVertex* _vtx = new PrimaryVertex(leptonsType_, muonsMomenta, electronsMomenta);
      _vtx->SetPosition(vertex->x(), vertex->y(), vertex->z());

      // Loop on all the tracks matched with this vertex    
      reco::Vertex::trackRef_iterator _track;
    	for (_track=vertex->tracks_begin(); _track!=vertex->tracks_end(); _track++) {
    	  Int_t leptonId_ = _vtx->AddTrack((*_track).castTo<reco::TrackRef>(), *_leptonType);
        if (leptonId_==-1) { // Track was not matched to any of the leptons in the collection
          continue;
        }
        std::cout << "--> Lepton " << leptonId_ << " (type : " << *_leptonType << ")" << std::endl;
        nLeptonsInPrimVertex += 1;
    	}
    	if (nLeptonsInPrimVertex<2) continue;
    	
    	// At this stage we have at least two matched leptons track on the vertex
      Int_t lep1, lep2;
      Pair_candidates[nPrimVertexCand][0] = -1;
      Pair_candidates[nPrimVertexCand][1] = -1;

    	if (_fetchElectrons && _fetchMuons) { // Looks at electron+muon
    	  if (_vtx->MatchedElectrons.size()==0 or _vtx->MatchedMuons.size()==0) {
    	    continue;
    	  }
    	  // Enough muons candidates on the vertex
    	  if (_vtx->MatchedElectrons.size()>=1 or _vtx->MatchedMuons.size()>=1) {
          minDist = 999.;
          foundPair = false;
          for(UInt_t i=0; i<_vtx->MatchedMuons.size(); i++) {
            lep1 = _vtx->MatchedMuons[i];
            for(UInt_t j=0; j<_vtx->MatchedElectrons.size(); j++) {
              lep2 = _vtx->MatchedElectrons[j];
              if (MuonCand_charge[lep1]*EleCand_charge[lep2]>0) {
                continue;
              }
              foundPair=true;
              leptonsDist=sqrt(pow(MuonCand_vtxx[lep1]-EleCand_vtxx[lep2],2)
                              +pow(MuonCand_vtxy[lep1]-EleCand_vtxy[lep2],2)
                              +pow(MuonCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
              if (leptonsDist<minDist) {
                minDist = leptonsDist;
                Pair_candidates[nPrimVertexCand][0] = lep1;
                Pair_candidates[nPrimVertexCand][1] = lep2;
              }
            }
          }
          if (Pair_candidates[nPrimVertexCand][0]==-1
           || Pair_candidates[nPrimVertexCand][1]==-1) foundPair = false;
          l1.SetXYZM(MuonCand_px[Pair_candidates[nPrimVertexCand][0]],
                     MuonCand_py[Pair_candidates[nPrimVertexCand][0]],
                     MuonCand_pz[Pair_candidates[nPrimVertexCand][0]],
                     MASS_MU);
          l2.SetXYZM(EleCand_px[Pair_candidates[nPrimVertexCand][1]],
                     EleCand_py[Pair_candidates[nPrimVertexCand][1]],
                     EleCand_pz[Pair_candidates[nPrimVertexCand][1]],
                     MASS_E);
        }    	  
    	}
    	else if (_fetchElectrons) { // Looks at dielectrons
    	  if (_vtx->MatchedElectrons.size()<2) {
      	  // Not enough electrons candidates on the vertex
    	    continue;
    	  }
    	  // Enough electrons candidates on the vertex
    	  if (_vtx->MatchedElectrons.size()>=2) {
          minDist = 999.;
          foundPair = false;
          for(UInt_t i=0; i<_vtx->MatchedElectrons.size(); i++) {
            lep1 = _vtx->MatchedElectrons[i];
            for(UInt_t j=i+1; j<_vtx->MatchedElectrons.size(); j++) {
              lep2 = _vtx->MatchedElectrons[j];
              if (EleCand_charge[lep1]*EleCand_charge[lep2]>0) {
                continue;
              }
              foundPair=true;
              leptonsDist=sqrt(pow(EleCand_vtxx[lep1]-EleCand_vtxx[lep2],2)
                              +pow(EleCand_vtxy[lep1]-EleCand_vtxy[lep2],2)
                              +pow(EleCand_vtxz[lep1]-EleCand_vtxz[lep2],2));
              if (leptonsDist<minDist) {
                minDist = leptonsDist;
                Pair_candidates[nPrimVertexCand][0] = lep1;
                Pair_candidates[nPrimVertexCand][1] = lep2;
              }
            }
          }
          if (Pair_candidates[nPrimVertexCand][0]==-1
           || Pair_candidates[nPrimVertexCand][1]==-1) foundPair = false;
          l1.SetXYZM(EleCand_px[Pair_candidates[nPrimVertexCand][0]],
                     EleCand_py[Pair_candidates[nPrimVertexCand][0]],
                     EleCand_pz[Pair_candidates[nPrimVertexCand][0]],
                     MASS_E);
          l2.SetXYZM(EleCand_px[Pair_candidates[nPrimVertexCand][1]],
                     EleCand_py[Pair_candidates[nPrimVertexCand][1]],
                     EleCand_pz[Pair_candidates[nPrimVertexCand][1]],
                     MASS_E);
        }
    	}
    	else if (_fetchMuons) { // Looks at dimuons
    	  // Enough muons candidates on the vertex
    	  if (_vtx->MatchedMuons.size()>=2) {
          minDist = 999.;
          foundPair = false;
          for(UInt_t i=0; i<_vtx->MatchedMuons.size(); i++) {
            lep1 = _vtx->MatchedMuons[i];
            for(UInt_t j=i+1; j<_vtx->MatchedMuons.size(); j++) {
              lep2 = _vtx->MatchedMuons[j];
              if (MuonCand_charge[lep1]*MuonCand_charge[lep2]>0) {
                continue;
              }
              foundPair=true;
              leptonsDist=sqrt(pow(MuonCand_vtxx[lep1]-MuonCand_vtxx[lep2],2)
                              +pow(MuonCand_vtxy[lep1]-MuonCand_vtxy[lep2],2)
                              +pow(MuonCand_vtxz[lep1]-MuonCand_vtxz[lep2],2));
              if (leptonsDist<minDist) {
                minDist = leptonsDist;
                Pair_candidates[nPrimVertexCand][0] = lep1;
                Pair_candidates[nPrimVertexCand][1] = lep2;
              }
            }
          }
          if (Pair_candidates[nPrimVertexCand][0]==-1
           || Pair_candidates[nPrimVertexCand][1]==-1) foundPair = false;
          l1.SetXYZM(MuonCand_px[Pair_candidates[nPrimVertexCand][0]],
                     MuonCand_py[Pair_candidates[nPrimVertexCand][0]],
                     MuonCand_pz[Pair_candidates[nPrimVertexCand][0]],
                     MASS_MU);
          l2.SetXYZM(MuonCand_px[Pair_candidates[nPrimVertexCand][1]],
                     MuonCand_py[Pair_candidates[nPrimVertexCand][1]],
                     MuonCand_pz[Pair_candidates[nPrimVertexCand][1]],
                     MASS_MU);
        }
    	}

      if (foundPair) {
        std::cout << "Pair : (" << Pair_candidates[nPrimVertexCand][0] << ", " << Pair_candidates[nPrimVertexCand][1] << ") at distance " << minDist << std::endl;

      	std::cout << "Matched muons : " << std::endl;
      	for (UInt_t i=0; i<_vtx->MatchedMuons.size(); i++) {
      	  std::cout << "-> " << _vtx->MatchedMuons[i] << std::endl;
      	}
      	std::cout << "Matched electrons : " << std::endl;
      	for (UInt_t i=0; i<_vtx->MatchedElectrons.size(); i++) {
      	  std::cout << "-> " << _vtx->MatchedElectrons[i] << std::endl;
      	}
      	/*if (_vtx->Flatten()) {
      	  foundDiLeptonPair = true;
      	}*/
        
        Pair_p[nPrimVertexCand] = (l1+l2).P();
        Pair_pt[nPrimVertexCand] = (l1+l2).Pt();
        Pair_mass[nPrimVertexCand] = (l1+l2).M();
        Pair_phi[nPrimVertexCand] = (l1+l2).Phi();
        Pair_eta[nPrimVertexCand] = (l1+l2).Eta();
        double dphi = fabs(l1.Phi()-l2.Phi());
        Pair_dphi[nPrimVertexCand] = (dphi<pi) ? dphi : (2.*pi)-dphi;
        Pair_dpt[nPrimVertexCand] = fabs(l1.Pt()-l2.Pt());
        Pair_3Dangle[nPrimVertexCand] = (l1.Angle(l2.Vect()))/pi;
      }
      
      nPrimVertexCand += 1;

      delete _vtx;
    }
    /*std::cout << "Number of primary vertices in the event : " << nPrimVertexCand << std::endl;
    std::cout << "Number of muons in the event : " << nMuonCand << std::endl;
    std::cout << "Number of electrons in the event : " << nEleCand << std::endl;*/
    delete _leptonType;
  }
  
  if (foundPair) {
    tree->Fill();
    nCandidates += 1;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaLL::beginJob()
{
  // Booking the ntuple
  tree->Branch("nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/I");
  tree->Branch("nUnfilteredPrimVertexCand", &nUnfilteredPrimVertexCand, "nUnfilteredPrimVertexCand/I");
  if (_fetchMuons) {
    tree->Branch("nMuonCand", &nMuonCand, "nMuonCand/I");
    tree->Branch("MuonCand_px", MuonCand_px,"MuonCand_px[nMuonCand]/D");
    tree->Branch("MuonCand_py", MuonCand_py,"MuonCand_py[nMuonCand]/D");
    tree->Branch("MuonCand_pz", MuonCand_pz,"MuonCand_pz[nMuonCand]/D");
    tree->Branch("MuonCand_p", MuonCand_p,"MuonCand_p[nMuonCand]/D");
    tree->Branch("MuonCand_pt", MuonCand_pt,"MuonCand_pt[nMuonCand]/D");
    tree->Branch("MuonCand_eta", MuonCand_eta,"MuonCand_eta[nMuonCand]/D");
    tree->Branch("MuonCand_phi", MuonCand_phi,"MuonCand_phi[nMuonCand]/D");
    tree->Branch("MuonCand_charge", MuonCand_charge,"MuonCand_charge[nMuonCand]/I");
    tree->Branch("MuonCand_vtxx", MuonCand_vtxx,"MuonCand_vtxx[nMuonCand]/D");
    tree->Branch("MuonCand_vtxy", MuonCand_vtxy,"MuonCand_vtxy[nMuonCand]/D");
    tree->Branch("MuonCand_vtxz", MuonCand_vtxz,"MuonCand_vtxz[nMuonCand]/D");
  }
  if (_fetchElectrons) {
    tree->Branch("nEleCand", &nEleCand,"nEleCand/I");
    tree->Branch("EleCand_px", EleCand_px,"EleCand_px[nEleCand]/D");
    tree->Branch("EleCand_py", EleCand_py,"EleCand_py[nEleCand]/D");
    tree->Branch("EleCand_pz", EleCand_pz,"EleCand_pz[nEleCand]/D");
    tree->Branch("EleCand_p", EleCand_p,"EleCand_p[nEleCand]/D");
    tree->Branch("EleCand_e", EleCand_e,"EleCand_e[nEleCand]/D");
    tree->Branch("EleCand_et", EleCand_et,"EleCand_et[nEleCand]/D");
    tree->Branch("EleCand_eta", EleCand_eta,"EleCand_eta[nEleCand]/D");
    tree->Branch("EleCand_phi", EleCand_phi,"EleCand_phi[nEleCand]/D");
    tree->Branch("EleCand_charge", EleCand_charge,"EleCand_charge[nEleCand]/I");
    tree->Branch("EleCand_vtxx",EleCand_vtxx,"EleCand_vtxx[nEleCand]/D");
    tree->Branch("EleCand_vtxy",EleCand_vtxy,"EleCand_vtxy[nEleCand]/D");
    tree->Branch("EleCand_vtxz",EleCand_vtxz,"EleCand_vtxz[nEleCand]/D");
    tree->Branch("EleCandTrack_p",EleCandTrack_p,"EleCandTrack_p[nEleCand]/D");
    tree->Branch("EleCandTrack_pt",EleCandTrack_pt,"EleCandTrack_pt[nEleCand]/D");
    tree->Branch("EleCandTrack_eta",EleCandTrack_eta,"EleCandTrack_eta[nEleCand]/D");
    tree->Branch("EleCandTrack_phi",EleCandTrack_phi,"EleCandTrack_phi[nEleCand]/D");
    tree->Branch("EleCandTrack_vtxz",EleCandTrack_vtxz,"EleCandTrack_vtxz[nEleCand]/D");
    tree->Branch("EleCand_deltaPhi",EleCand_deltaPhi,"EleCand_deltaPhi[nEleCand]/D"); 
    tree->Branch("EleCand_deltaEta",EleCand_deltaEta,"EleCand_deltaEta[nEleCand]/D"); 
    tree->Branch("EleCand_HoverE",EleCand_HoverE,"EleCand_HoverE[nEleCand]/D"); 
    tree->Branch("EleCand_trackiso",EleCand_trackiso,"EleCand_trackiso[nEleCand]/D"); 
    tree->Branch("EleCand_ecaliso",EleCand_ecaliso,"EleCand_ecaliso[nEleCand]/D"); 
    tree->Branch("EleCand_hcaliso",EleCand_hcaliso,"EleCand_hcaliso[nEleCand]/D"); 
    tree->Branch("EleCand_sigmaIetaIeta",EleCand_sigmaIetaIeta,"EleCand_sigmaIetaIeta[nEleCand]/D"); 
    tree->Branch("EleCand_convDist",EleCand_convDist,"EleCand_convDist[nEleCand]/D"); 
    tree->Branch("EleCand_convDcot",EleCand_convDcot,"EleCand_convDcot[nEleCand]/D"); 
    tree->Branch("EleCand_ecalDriven",EleCand_ecalDriven,"EleCand_ecalDriven[nEleCand]/D");  
    tree->Branch("EleCand_wp80", EleCand_wp80, "EleCand_wp80[nEleCand]/I");   
    tree->Branch("EleCand_mediumID", EleCand_mediumID, "EleCand_mediumID[nEleCand]/I");    
    tree->Branch("EleCand_looseID", EleCand_looseID, "EleCand_looseID[nEleCand]/I");   
    /*tree->Branch("EleCand_looseid", EleCand_looseid,"EleCand_looseid[nEleCand]/I");
    tree->Branch("EleCand_likelihoodid", EleCand_likelihoodid,"EleCand_likelihoodid[nEleCand]/D");
    tree->Branch("EleCand_robustid", EleCand_robustid,"EleCand_robustid[nEleCand]/I");*/
  }
  tree->Branch("Pair_candidates", &Pair_candidates, "Pair_candidates[nPrimVertexCand][2]/I");
  tree->Branch("Pair_p",Pair_p,"Pair_p[nPrimVertexCand]/D");
  tree->Branch("Pair_pt",Pair_pt,"Pair_pt[nPrimVertexCand]/D");
  tree->Branch("Pair_dpt",Pair_dpt,"Pair_dpt[nPrimVertexCand]/D");
  tree->Branch("Pair_mass",Pair_mass,"Pair_mass[nPrimVertexCand]/D");
  tree->Branch("Pair_eta",Pair_eta,"Pair_eta[nPrimVertexCand]/D");
  tree->Branch("Pair_phi",Pair_phi,"Pair_phi[nPrimVertexCand]/D");
  tree->Branch("Pair_dphi",Pair_dphi,"Pair_dphi[nPrimVertexCand]/D");
  tree->Branch("Pair_3Dangle",Pair_3Dangle,"Pair_3Dangle[nPrimVertexCand]/D");

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
GammaGammaLL::beginRun(edm::Run const&, edm::EventSetup const&)
{
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
  nMatchedMuons(0),
  nMatchedElectrons(0),
  nMatchedTracks(0),
  FetchMuons(false),
  FetchElectrons(false)
{
  //LeptonsType = _leptonsType;
  for (UInt_t i=0; i<_leptonsType.size(); i++) {
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
  return -1;
}

/*bool
PrimaryVertex::ComputeKinematicQuantities() {
  if (nMatchedTracks!=2) {
#ifdef DEBUG
    std::cout << "[PrimaryVertex] ComputeKinematicQuantities :: " << nMatchedTracks << " tracks were matched to the vertex!" << std::endl;
#endif
    return false;
  }
  if (FetchMuons && FetchElectrons) {
    if (nMatchedMuons==1 && nMatchedElectrons==1) {
      Double_t invariantMass = (MatchedMuons[0]+MatchedElectrons[0]).M();
      Double_t pTpair = (MatchedMuons[0]+MatchedElectrons[0]).Pt();
      std::cout << "[PrimaryVertex] ComputeKinematicQuantities :: Invariant mass of the dilepton pair : " << invariantMass << std::endl;
      std::cout << "[PrimaryVertex] ComputeKinematicQuantities :: Transverse momentum of the dilepton pair : " << pTpair << std::endl;
    }
  }
  else if (FetchMuons) {
    
  }
  else if (FetchElectrons) {
    
  }
  return true;
}*/

bool
PrimaryVertex::Flatten() {
  /*if (!ComputeKinematicQuantities()) {
#ifdef DEBUG
    std::cout << "[PrimaryVertex] Flatten :: Impossible to compute the kinematic quantities on the vertex!" << std::endl;
#endif
    return false;
  }*/
  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaGammaLL);
