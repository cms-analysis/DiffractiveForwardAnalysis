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
// $Id: GammaGammaEE.cc,v 1.20 2008/06/30 12:38:44 jjhollar Exp $
//
//


// system include files
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h" 
#include "DataFormats/PatCandidates/interface/Electron.h" 
#include "DataFormats/PatCandidates/interface/Tau.h" 
#include "DataFormats/PatCandidates/interface/Photon.h" 
#include "DataFormats/PatCandidates/interface/MET.h" 

#include "DataFormats/Common/interface/TriggerResults.h"   
#include "FWCore/Framework/interface/TriggerNames.h"   

#include "FWCore/Framework/interface/ESHandle.h" 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/CaloJet.h" 
#include "DataFormats/EgammaCandidates/interface/Electron.h" 
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h" 
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"   
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h" 
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h" 
#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/MuonReco/interface/Muon.h" 
#include "DataFormats/MuonReco/interface/MuonFwd.h"  
#include "DataFormats/METReco/interface/CaloMET.h" 
#include "DataFormats/METReco/interface/CaloMETFwd.h"  
#include "DataFormats/METReco/interface/CaloMETCollection.h" 
#include "DataFormats/EgammaCandidates/interface/Photon.h" 
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" 
#include "DataFormats/CaloTowers/interface/CaloTower.h" 
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"  
#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/Candidate/interface/CandidateFwd.h"  
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" 


#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaEE.h"

// user include files
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
// Muons:
#include <DataFormats/TrackReco/interface/Track.h>
// Electrons
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

// Vertexing 
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrack.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h" 
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h" 
#include "SimTracker/Records/interface/TrackAssociatorRecord.h" 
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" 

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
  theCaloTowLabel    = pset.getParameter<edm::InputTag>("CaloTowerLabel");

  eldetmax           = pset.getParameter<double>("DielectronMaxdEt");
  eldphimin          = pset.getParameter<double>("DielectronMindphi");
  drisocalo          = pset.getParameter<double>("CaloTowerdR");

  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");

  //  nEvt=0;
  ELEMAX=10;
  JETMAX=30;
  TRACKMAX=100;

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
  thetree->Branch("HighestJet_e",&HighestJet_e,"HighestJet_e/D");
  thetree->Branch("HighestJet_eta",&HighestJet_eta,"HighestJet_eta/D"); 
  thetree->Branch("HighestJet_phi",&HighestJet_phi,"HighestJet_phi/D"); 
  thetree->Branch("SumJet_e",&SumJet_e,"SumJet_e/D");

  thetree->Branch("nCaloCand",&nCaloCand,"nCaloCand/I");
  thetree->Branch("CaloTower_e",CaloTower_e,"CaloTower_e[nCaloCand]/D");
  thetree->Branch("CaloTower_et",CaloTower_et,"CaloTower_et[nCaloCand]/D");
  thetree->Branch("CaloTower_eta",CaloTower_eta,"CaloTower_eta[nCaloCand]/D"); 
  thetree->Branch("CaloTower_phi",CaloTower_phi,"CaloTower_phi[nCaloCand]/D"); 
  thetree->Branch("CaloTower_dr",CaloTower_dr,"CaloTower_dr[nCaloCand]/D");
  thetree->Branch("HighestCaloTower_e",&HighestCaloTower_e,"HighestCaloTower_e/D");
  thetree->Branch("HighestCaloTower_eta",&HighestCaloTower_eta,"HighestCaloTower_eta/D");
  thetree->Branch("HighestCaloTower_phi",&HighestCaloTower_phi,"HighestCaloTower_phi/D"); 
  thetree->Branch("HighestCaloTower_dr",&HighestCaloTower_dr,"HighestCaloTower_dr/D");
  thetree->Branch("HighestEtCaloTower_et",&HighestEtCaloTower_et,"HighestEtCaloTower_et/D");
  thetree->Branch("HighestEtCaloTower_eta",&HighestEtCaloTower_eta,"HighestEtCaloTower_eta/D");
  thetree->Branch("HighestEtCaloTower_phi",&HighestEtCaloTower_phi,"HighestEtCaloTower_phi/D"); 
  thetree->Branch("HighestEtCaloTower_dr",&HighestEtCaloTower_dr,"HighestEtCaloTower_dr/D");
  thetree->Branch("SumCalo_e",&SumCalo_e,"SumCalo_e/D");

  thetree->Branch("nExtraCaloTowersE1",&nExtraCaloTowersE1,"nExtraCaloTowersE1/I");
  thetree->Branch("nExtraCaloTowersE2",&nExtraCaloTowersE2,"nExtraCaloTowersE2/I");
  thetree->Branch("nExtraCaloTowersE3",&nExtraCaloTowersE3,"nExtraCaloTowersE3/I"); 
  thetree->Branch("nExtraCaloTowersE4",&nExtraCaloTowersE4,"nExtraCaloTowersE4/I"); 
  thetree->Branch("nExtraCaloTowersE5",&nExtraCaloTowersE5,"nExtraCaloTowersE5/I"); 
  thetree->Branch("nExtraCaloTowersE6",&nExtraCaloTowersE6,"nExtraCaloTowersE6/I");    
  thetree->Branch("nExtraCaloTowersE7",&nExtraCaloTowersE7,"nExtraCaloTowersE7/I");    
  thetree->Branch("nExtraCaloTowersE8",&nExtraCaloTowersE8,"nExtraCaloTowersE8/I");    
  thetree->Branch("nExtraCaloTowersE9",&nExtraCaloTowersE9,"nExtraCaloTowersE9/I");    

  thetree->Branch("nExtraCaloTowersEt0pt1",&nExtraCaloTowersEt0pt1,"nExtraCaloTowersEt0pt1/I");  
  thetree->Branch("nExtraCaloTowersEt0pt2",&nExtraCaloTowersEt0pt2,"nExtraCaloTowersEt0pt2/I");  
  thetree->Branch("nExtraCaloTowersEt0pt5",&nExtraCaloTowersEt0pt5,"nExtraCaloTowersEt0pt5/I");  
  thetree->Branch("nExtraCaloTowersEt1",&nExtraCaloTowersEt1,"nExtraCaloTowersEt1/I");  
  thetree->Branch("nExtraCaloTowersEt2",&nExtraCaloTowersEt2,"nExtraCaloTowersEt2/I");  

  thetree->Branch("nTrackCand",&nTrackCand,"nTrackCand/I");
  thetree->Branch("TrackCand_px",TrackCand_px,"TrackCand_px[nTrackCand]/D");
  thetree->Branch("TrackCand_py",TrackCand_py,"TrackCand_py[nTrackCand]/D");
  thetree->Branch("TrackCand_pz",TrackCand_pz,"TrackCand_pz[nTrackCand]/D");
  thetree->Branch("TrackCand_p",TrackCand_p,"TrackCand_p[nTrackCand]/D");
  thetree->Branch("TrackCand_pt",TrackCand_pt,"TrackCand_pt[nTrackCand]/D");
  thetree->Branch("TrackCand_eta",TrackCand_eta,"TrackCand_eta[nTrackCand]/D");
  thetree->Branch("TrackCand_phi",TrackCand_phi,"TrackCand_phi[nTrackCand]/D");
  thetree->Branch("TrackCand_vtxdxyz",TrackCand_vtxdxyz,"TrackCand_vtxdxyz[nTrackCand]/D"); 
  thetree->Branch("ClosestExtraTrack_vtxdxyz",&ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz/D"); 

  thetree->Branch("ElEl_mass",&ElEl_mass,"ElEl_mass/D");
  thetree->Branch("ElEl_dphi",&ElEl_dphi,"ElEl_dphi/D");
  thetree->Branch("ElEl_vtxx",&ElEl_vtxx,"ElEl_vtxx/D");
  thetree->Branch("ElEl_vtxy",&ElEl_vtxy,"ElEl_vtxy/D"); 
  thetree->Branch("ElEl_vtxz",&ElEl_vtxz,"ElEl_vtxz/D"); 
  thetree->Branch("ElEl_vtxchi2dof",&ElEl_vtxchi2dof,"ElEl_vtxchi2dof/D");
  thetree->Branch("ElEl_vtxisvalid",&ElEl_vtxisvalid,"ElEl_vtxisvalid/I");
  thetree->Branch("ElEl_extratracks2cm",&ElEl_extratracks2cm,"ElEl_extratracks2cm/I"); 
  thetree->Branch("ElEl_extratracks5cm",&ElEl_extratracks5cm,"ElEl_extratracks5cm/I");  
  thetree->Branch("ElEl_extratracks10cm",&ElEl_extratracks10cm,"ElEl_extratracks10cm/I");  

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");

  thetree->Branch("HLT2Electron5_L1R_NI",&HLT2Electron5_L1R_NI,"HLT2Electron5_L1R_NI/I"); 
  thetree->Branch("HLT2ElectronExclusive",&HLT2ElectronExclusive,"HLT2ElectronExclusive/I");  

  //  thetree->Branch("evweight",&evweight,"evweight/D");
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
  nCaloCand=0;
  nTrackCand=0;

  nExtraCaloTowersE1=0; 
  nExtraCaloTowersE2=0; 
  nExtraCaloTowersE3=0;  
  nExtraCaloTowersE4=0;  
  nExtraCaloTowersE5=0;
  nExtraCaloTowersE6=0;  
  nExtraCaloTowersE7=0;   
  nExtraCaloTowersE8=0;   
  nExtraCaloTowersE9=0;    
  nExtraCaloTowersEt0pt1=0;  
  nExtraCaloTowersEt0pt2=0; 
  nExtraCaloTowersEt0pt5=0;   
  nExtraCaloTowersEt1=0;  
  nExtraCaloTowersEt2=0;   
  HitInZDC=0;
  HitInCastor=0;

  ElEl_mass = -1;
  ElEl_dphi = -1;

  bool passed = true;

 //using namespace edm;
  using reco::TrackCollection;
  
  //  Handle< double> weightHandle; 
  //  event.getByLabel ("weight", weightHandle); 
  //  evweight = * weightHandle; 

  // Get the trigger information from the event 
  edm::Handle<edm::TriggerResults> hltResults ;  
  event.getByLabel(InputTag("TriggerResults::HLT"),hltResults) ;  
  trigNames.init(*hltResults) ;  
  for (unsigned int i=0; i<trigNames.size(); i++)   
    {  
      if ( trigNames.triggerNames().at(i) == "HLT2Electron5_L1R_NI" )        
        {   
          if ( hltResults->accept(i) )   
            HLT2Electron5_L1R_NI = 1; 
          else 
            HLT2Electron5_L1R_NI = 0; 
        }   
      if ( trigNames.triggerNames().at(i) == "HLT2ElectronExclusive" )  
        {    
          if ( hltResults->accept(i) )   
            HLT2ElectronExclusive = 1; 
          else 
            HLT2ElectronExclusive = 0; 
        }    
    } 

  // Get the electron collection from the event

  // PAT
  edm::Handle<edm::View<pat::Electron> > electrons; 
  event.getByLabel(thePixelGsfELabel,electrons); 
  edm::View<pat::Electron>::const_iterator electron;

  // AOD
  //  edm::Handle<reco::PixelMatchGsfElectronCollection> pTracks;
  //  event.getByLabel(thePixelGsfELabel,pTracks);
  //  const reco::PixelMatchGsfElectronCollection* electrons = pTracks.product();
  //  reco::PixelMatchGsfElectronCollection::const_iterator electron;

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

	  EleCandTrack_p[nEleCand] = electron->gsfTrack()->p(); 

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

  // PAT 
  edm::Handle<edm::View<pat::Jet> > jets;  
  event.getByLabel(theJetLabel,jets);  
  edm::View<pat::Jet>::const_iterator jet; 
 
  // AOD 
  //  edm::Handle<reco::CaloJetCollection> pJets; 
  //  event.getByLabel(theJetLabel,pJets); 
  //  const reco::CaloJetCollection* jets = pJets.product(); 
  //  reco::CaloJetCollection::const_iterator jet; 

  // Get the MET collection from the event

  // PAT
  edm::Handle<edm::View<pat::MET> > mets;  
  event.getByLabel(theMetLabel,mets);  
  edm::View<pat::MET>::const_iterator met; 
 
  // AOD 
  //  edm::Handle<reco::CaloMETCollection> pMET; 
  //  event.getByLabel(theMetLabel,pMET); 
  //  const reco::CaloMETCollection* mets = pMET.product(); 
  //  reco::CaloMETCollection::const_iterator met; 

  // Get the CaloTower collection from the event 
  edm::Handle<CaloTowerCollection> caloTowers;  
  event.getByLabel(theCaloTowLabel,caloTowers);  
  const CaloTowerCollection* towers = caloTowers.product();  
  CaloTowerCollection::const_iterator calo;  

  // Get the track collection from the event
  edm::Handle<reco::TrackCollection> recoTracks;
  event.getByLabel(recTrackLabel, recoTracks);
  const TrackCollection* tracks = recoTracks.product();
  TrackCollection::const_iterator track;

  double highestejet = -1.0;
  double highestejeteta = -999.0;
  double highestejetphi = -999.0;
  double totalejet = -1.0;
  double highestetower = -1.0; 
  double highestetowerdr = -999.0;
  double highestetowereta = -999.0;
  double highestetowerphi = -999.0;
  double highestettower = -1.0; 
  double highestettowerdr = -999.0;
  double highestettowereta = -999.0;
  double highestettowerphi = -999.0;
  double totalecalo = -1.0; 
  double closesttrkdxyz = 999.0; 

  // If this event contains a di-mu/e/gamma candidate, look at Jets & MET & CaloTowers & Tracks
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
	    {
	      highestejet = JetCand_e[nJetCand];
	      highestejeteta = JetCand_eta[nJetCand];
	      highestejetphi = JetCand_phi[nJetCand];
	    }
	  nJetCand++;
	}

      HighestJet_e = highestejet;
      HighestJet_eta = highestejeteta;
      HighestJet_phi = highestejetphi;
      SumJet_e = totalejet;

      met = mets->begin();
      float e_met = met->energy();
      Etmiss = e_met;

      for (calo = towers->begin(); calo != towers->end(); ++calo )
	{
	  CaloTower_e[nCaloCand]=calo->energy(); 
	  CaloTower_et[nCaloCand]=calo->et();
	  CaloTower_phi[nCaloCand]=calo->phi(); 
	  CaloTower_eta[nCaloCand]=calo->eta(); 
	  

	  float calodr1 = sqrt(((CaloTower_eta[nCaloCand]-EleCand_eta[0])*(CaloTower_eta[nCaloCand]-EleCand_eta[0])) + 
			       ((CaloTower_phi[nCaloCand]-EleCand_phi[0])*(CaloTower_phi[nCaloCand]-EleCand_phi[0])));
	  float calodr2 = sqrt(((CaloTower_eta[nCaloCand]-EleCand_eta[1])*(CaloTower_eta[nCaloCand]-EleCand_eta[1])) + 
			       ((CaloTower_phi[nCaloCand]-EleCand_phi[1])*(CaloTower_phi[nCaloCand]-EleCand_phi[1])));
	  
	  if(calodr1 < calodr2)
	    CaloTower_dr[nCaloCand] = calodr1;
	  else
	    CaloTower_dr[nCaloCand] = calodr2;

	  totalecalo = totalecalo + CaloTower_e[nCaloCand]; 
	  if(CaloTower_e[nCaloCand] > highestetower) 
	    {
	      highestetower = CaloTower_e[nCaloCand]; 
	      highestetowereta = CaloTower_eta[nCaloCand];
	      highestetowerphi = CaloTower_phi[nCaloCand];
	      highestetowerdr = CaloTower_dr[nCaloCand];
	    }

	  if(CaloTower_et[nCaloCand] > highestettower) 
	    {
	      highestettower = CaloTower_et[nCaloCand]; 
	      highestettowereta = CaloTower_eta[nCaloCand];
	      highestettowerphi = CaloTower_phi[nCaloCand];
	      highestettowerdr = CaloTower_dr[nCaloCand];
	    }

	  if(CaloTower_dr[nCaloCand] > drisocalo)
	    {
	      if(CaloTower_e[nCaloCand] > 1.0)
		nExtraCaloTowersE1++;
              if(CaloTower_e[nCaloCand] > 2.0) 
                nExtraCaloTowersE2++; 
              if(CaloTower_e[nCaloCand] > 3.0) 
                nExtraCaloTowersE3++; 
              if(CaloTower_e[nCaloCand] > 4.0) 
                nExtraCaloTowersE4++; 
              if(CaloTower_e[nCaloCand] > 5.0) 
                nExtraCaloTowersE5++; 
              if(CaloTower_e[nCaloCand] > 6.0)    
                nExtraCaloTowersE6++;    
              if(CaloTower_e[nCaloCand] > 7.0)    
                nExtraCaloTowersE7++;    
              if(CaloTower_e[nCaloCand] > 8.0)    
                nExtraCaloTowersE8++;    
              if(CaloTower_e[nCaloCand] > 9.0)    
                nExtraCaloTowersE9++;    

	      if(CaloTower_et[nCaloCand] > 0.1)
		nExtraCaloTowersEt0pt1++;
              if(CaloTower_et[nCaloCand] > 0.2) 
                nExtraCaloTowersEt0pt2++; 
              if(CaloTower_et[nCaloCand] > 0.5) 
                nExtraCaloTowersEt0pt5++; 
              if(CaloTower_et[nCaloCand] > 1.0) 
                nExtraCaloTowersEt1++; 
              if(CaloTower_et[nCaloCand] > 2.0) 
                nExtraCaloTowersEt2++; 

	    }
	  
	  nCaloCand++;

	}      

      SumCalo_e = totalecalo;
      HighestCaloTower_e = highestetower;
      HighestCaloTower_eta = highestetowereta;
      HighestCaloTower_phi = highestetowerphi;
      HighestCaloTower_dr = highestetowerdr;
      HighestEtCaloTower_et = highestettower;
      HighestEtCaloTower_eta = highestettowereta;
      HighestEtCaloTower_phi = highestettowerphi;
      HighestEtCaloTower_dr = highestettowerdr;
      
    }

  // Check for particles in ZDC/Castor acceptance. 
  // Use MC truth for now, replace with real RECO when available
  double MCPar_px,MCPar_py,MCPar_pz,MCPar_e,MCPar_eta,MCPar_mass;
  int MCPar_pdgid;

  Handle<GenParticleCollection> genParticles; 
  event.getByLabel( "genParticles", genParticles ); 
  for ( size_t i = 0; i < genParticles->size(); ++ i ) 
    {
      const Candidate & p = (*genParticles)[ i ];
      MCPar_pdgid=p.pdgId();
      MCPar_eta=p.eta();
      MCPar_px=p.px();
      MCPar_py=p.py();
      MCPar_pz=p.pz();
      MCPar_mass=p.mass();
      MCPar_e = sqrt(MCPar_mass*MCPar_mass + (MCPar_px*MCPar_px + MCPar_py*MCPar_py + MCPar_pz*MCPar_pz));

      if(MCPar_pdgid == 22 && abs(MCPar_eta) > 8.6 && MCPar_e > 20.0) 
	HitInZDC++;
      if(MCPar_pdgid == 2112 && abs(MCPar_eta) > 8.6 && MCPar_e > 50.0)
	HitInZDC++;
      if((MCPar_pdgid != 22 && MCPar_pdgid != 2112) && (abs(MCPar_eta) > 5.2 && abs(MCPar_eta) < 6.6))
	HitInCastor++;
    }

  // Now do vertexing and track counting 
  edm::ESHandle<TransientTrackBuilder> theVtx; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theVtx); 
  //  vector < reco::TransientTrack > etrks; 
  vector<TransientTrack> transetrks;  
  reco::TrackCollection * etrks = new reco::TrackCollection; 

  // First get "electron" tracks 
  bool isElectron = false; 
  for( track = tracks->begin(); track != tracks->end(); ++ track )  
    {  
      isElectron = false; 
      for(int j = 0;j < nEleCand; j++) 
        { 
          if(EleCandTrack_p[j] == track->p()) 
            { 
              isElectron = true; 
              etrks->push_back( *track ); 
              TransientTrack tmptrk = (*theVtx).build( *track ); 
              transetrks.push_back( tmptrk ); 
            } 
        } 
    } 
 
  // If 2 electrons, make a vertex 
  if(transetrks.size() == 2)  
    {  
      KalmanVertexFitter fitter(true);  
      TransientVertex elelVertex = fitter.vertex(transetrks);  
      if(elelVertex.isValid()) 
        { 
          ElEl_vtxx = elelVertex.position().x();  
          ElEl_vtxy = elelVertex.position().y();  
          ElEl_vtxz = elelVertex.position().z();  
          ElEl_vtxchi2dof = elelVertex.normalisedChiSquared(); 
          ElEl_vtxisvalid = 1; 
        } 
      else 
        { 
          ElEl_vtxx = 0;   
          ElEl_vtxy = 0;   
          ElEl_vtxz = 0;   
          ElEl_vtxchi2dof = 0; 
          ElEl_vtxisvalid = 0; 
        } 
      // OK, now go back and count "extra" tracks on the dielectron vertex 
      for(track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track)  
        {  
          if(track->p() == EleCandTrack_p[0] || track->p() == EleCandTrack_p[1]) 
            continue; 
           
          TrackCand_p[nTrackCand]=track->p();  
          TrackCand_px[nTrackCand]=track->px();  
          TrackCand_py[nTrackCand]=track->py();  
          TrackCand_pz[nTrackCand]=track->pz();  
          TrackCand_pt[nTrackCand]=track->pt();  
          TrackCand_eta[nTrackCand]=track->eta();  
          TrackCand_phi[nTrackCand]=track->phi();  
          TrackCand_charge[nTrackCand]=track->charge();  
          TrackCand_vtxdxyz[nTrackCand] = sqrt(((track->vertex().x() - ElEl_vtxx)*(track->vertex().x() - ElEl_vtxx)) +  
					       ((track->vertex().y() - ElEl_vtxy)*(track->vertex().y() - ElEl_vtxy)) + 
					       ((track->vertex().z() - ElEl_vtxz)*(track->vertex().z() - ElEl_vtxz))); 
           
          if(TrackCand_vtxdxyz[nTrackCand] < 2) 
            ElEl_extratracks2cm++; 
          if(TrackCand_vtxdxyz[nTrackCand] < 5)  
            ElEl_extratracks5cm++;  
          if(TrackCand_vtxdxyz[nTrackCand] < 10)  
            ElEl_extratracks10cm++;  
          if(TrackCand_vtxdxyz[nTrackCand] < closesttrkdxyz) 
            closesttrkdxyz = TrackCand_vtxdxyz[nTrackCand]; 
 
          nTrackCand++;   
        }  
      ClosestExtraTrack_vtxdxyz = closesttrkdxyz; 
    }  
  else  
    {  
      ElEl_vtxx = 0;  
      ElEl_vtxy = 0;  
      ElEl_vtxz = 0;  
      ElEl_vtxchi2dof = 0;  
      ElEl_vtxisvalid = 0;  
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
  
