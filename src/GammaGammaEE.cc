// -*- C++ -*-
//
// Package:    GammaGammaEE
// Class:      GammaGammaEE
// 
/**\class GammaGammaEE GammaGammaEE.cc GammaGammaLeptonLepton/GammaGammaEE/src/GammaGammaEE.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaEE.cc,v 1.1 2006/09/21 10:27:18 anonymous Exp $
//
//


// system include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaEE.h"

// user include files
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
// Muons:
#include <DataFormats/TrackReco/interface/Track.h>
// Electrons
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

// C++
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

using namespace std;
using namespace edm;
using namespace reco;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GammaGammaEE::GammaGammaEE(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  recTrackLabel      = pset.getParameter<edm::InputTag>("RecoTrackLabel");
  theGLBMuonLabel    = pset.getParameter<edm::InputTag>("GlobalMuonCollectionLabel");
  thePixelGsfELabel  = pset.getParameter<edm::InputTag>("ElectronCollectionLabel");
  theJetLabel        = pset.getParameter<edm::InputTag>("JetCollectionLabel");
  theMetLabel        = pset.getParameter<edm::InputTag>("MetLabel");
  thePhotonLabel     = pset.getParameter<edm::InputTag>("PhotonCollectionLabel");


  eldetmax           = pset.getParameter<double>("DielectronMaxdEt");
  eldphimin          = pset.getParameter<double>("DielectronMindphi");
  njetsmax           = pset.getParameter<int>("JetMultMax");
  highestjetemax     = pset.getParameter<double>("HighestJetEMax");
  sumjetemax         = pset.getParameter<double>("SumJetEMax");

  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");

  //  nEvt=0;
  ELEMAX=10;
  JETMAX=30;
  
  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

  thetree->Branch("nEleCand",&nEleCand,"nEleCand/I");
  thetree->Branch("EleCand_px",EleCand_px,"EleCand_px[nEleCand]/D");
  thetree->Branch("EleCand_py",EleCand_py,"EleCand_py[nEleCand]/D");
  thetree->Branch("EleCand_pz",EleCand_pz,"EleCand_pz[nEleCand]/D");
  thetree->Branch("EleCand_p",EleCand_p,"EleCand_p[nEleCand]/D");
  thetree->Branch("EleCand_e",EleCand_e,"EleCand_e[nEleCand]/D");
  thetree->Branch("EleCand_et",EleCand_et,"EleCand_et[nEleCand]/D");
  thetree->Branch("EleCand_eta",EleCand_eta,"EleCand_eta[nEleCand]/D");
  thetree->Branch("EleCand_phi",EleCand_phi,"EleCand_phi[nEleCand]/D");
  thetree->Branch("EleCand_charge",EleCand_charge,"EleCand_charge[nEleCand]/D");

  thetree->Branch("nJetCand",&nJetCand,"nJetCand/I");
  thetree->Branch("JetCand_px",JetCand_px,"JetCand_px[nJetCand]/D");
  thetree->Branch("JetCand_py",JetCand_py,"JetCand_py[nJetCand]/D");
  thetree->Branch("JetCand_pz",JetCand_pz,"JetCand_pz[nJetCand]/D");
  thetree->Branch("JetCand_e",JetCand_e,"JetCand_e[nJetCand]/D");
  thetree->Branch("JetCand_eta",JetCand_eta,"JetCand_eta[nJetCand]/D");
  thetree->Branch("JetCand_phi",JetCand_phi,"JetCand_phi[nJetCand]/D");

  thetree->Branch("ElEl_mass",&ElEl_mass,"ElEl_mass/D");
  thetree->Branch("ElEl_dphi",&ElEl_dphi,"ElEl_dphi/D");
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");
}


GammaGammaEE::~GammaGammaEE()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaEE::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nEleCand=0;
  nJetCand=0;

  ElEl_mass = -1;
  ElEl_dphi = -1;

  bool passed = true;

 //using namespace edm;
  using reco::TrackCollection;
  
  // Get the electron track collection from the event
  edm::Handle<reco::PixelMatchGsfElectronCollection> pTracks;
  event.getByLabel(thePixelGsfELabel,pTracks);
  const reco::PixelMatchGsfElectronCollection* electrons = pTracks.product();
  reco::PixelMatchGsfElectronCollection::const_iterator electron;

  if(electrons->size() == 2)
    {
      for ( electron = electrons->begin(); electron != electrons->end() && nEleCand<ELEMAX; ++electron )
	{
	  EleCand_e[nEleCand]=electron->energy();
	  EleCand_et[nEleCand]=electron->et();
	  EleCand_px[nEleCand]=electron->px();
	  EleCand_py[nEleCand]=electron->py();
	  EleCand_pz[nEleCand]=electron->pz();
	  EleCand_p[nEleCand]=electron->p();
	  EleCand_phi[nEleCand]=electron->phi();
	  EleCand_eta[nEleCand]=electron->eta();
	  EleCand_charge[nEleCand]=electron->charge();
	  nEleCand++;
	}

      // Calculate invariant mass and delta-phi
      if(EleCand_charge[0]*EleCand_charge[1]<0)
	{
	  double mass = pow(EleCand_p[0]+EleCand_p[1],2);
	  mass-=pow(EleCand_px[0]+EleCand_px[1],2);
	  mass-=pow(EleCand_py[0]+EleCand_py[1],2);
	  mass-=pow(EleCand_pz[0]+EleCand_pz[1],2);
	  ElEl_mass = sqrt(mass);

	  double dphi = fabs(EleCand_phi[0]-EleCand_phi[1]);
	  if(dphi < 3.14159)
	    ElEl_dphi = dphi;
	  else
	    ElEl_dphi = (2.0*3.14159)-dphi;
	}
    }

  // Get the Jet collection from the event
  edm::Handle<reco::CaloJetCollection> pJets;
  event.getByLabel(theJetLabel,pJets);
  const reco::CaloJetCollection* jets = pJets.product();
  reco::CaloJetCollection::const_iterator jet;

  // Get the MET collection from the event
  edm::Handle<reco::CaloMETCollection> pMET;
  event.getByLabel(theMetLabel,pMET);
  const reco::CaloMETCollection* mets = pMET.product();
  reco::CaloMETCollection::const_iterator met;

  double highestejet = -1.0;
  double totalejet = -1.0;

  // If this event contains a di-mu/e/gamma candidate, look at Jets & MET
  if(nEleCand == 2)
    {
      for ( jet = jets->begin(); jet != jets->end() && nJetCand<JETMAX; ++jet )
	{
	  JetCand_e[nJetCand]=jet->energy();
	  JetCand_px[nJetCand]=jet->px();
	  JetCand_py[nJetCand]=jet->py();
	  JetCand_pz[nJetCand]=jet->pz();
	  JetCand_phi[nJetCand]=jet->phi();
	  JetCand_eta[nJetCand]=jet->eta();

	  totalejet = totalejet + JetCand_e[nJetCand];
	  if(JetCand_e[nJetCand] > highestejet)
	    highestejet = JetCand_e[nJetCand];

	  nJetCand++;
	}

      met = mets->begin();
      float e_met = met->energy();
      Etmiss = e_met;
    }

  // Check for di-objects
  if(nEleCand != 2)
    passed = false;
  else
    {
     if(ElEl_dphi < eldphimin)
       passed = false;
     if(fabs(EleCand_et[0]-EleCand_et[1]) > eldetmax)
       passed = false;
    }

  // "Exclusivity" cuts
  if(nJetCand > njetsmax || 
     highestejet > highestjetemax || 
     totalejet > sumjetemax)
    passed = false;

  if(passed == true)
    thetree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaEE::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaEE::endJob() {
  thefile->Write();
  thefile->Close();
}
  
