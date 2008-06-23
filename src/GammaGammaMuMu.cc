// -*- C++ -*-
//
// Package:    GammaGammaMuMu
// Class:      GammaGammaMuMu
// 
/**\class GammaGammaMuMu GammaGammaMuMu.cc GammaGammaLeptonLepton/GammaGammaMuMu/src/GammaGammaMuMu.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaMuMu.cc,v 1.17 2008/03/18 14:43:32 jjhollar Exp $
//
//


// system include files
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h" 
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"  
#include "DataFormats/Common/interface/Ref.h"  
 
#include "DataFormats/Common/interface/TriggerResults.h"  
#include "FWCore/Framework/interface/TriggerNames.h"  
  
 
#include "FWCore/Framework/interface/ESHandle.h" 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/EgammaCandidates/interface/Electron.h"  
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"   
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

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMu.h"

// user include files
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
// Muons:
#include <DataFormats/TrackReco/interface/Track.h>
// Electrons
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"  


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
GammaGammaMuMu::GammaGammaMuMu(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  recTrackLabel      = pset.getParameter<edm::InputTag>("RecoTrackLabel");
  theGLBMuonLabel    = pset.getParameter<edm::InputTag>("GlobalMuonCollectionLabel");
  thePixelGsfELabel  = pset.getParameter<edm::InputTag>("ElectronCollectionLabel");
  theJetLabel        = pset.getParameter<edm::InputTag>("JetCollectionLabel");
  theMetLabel        = pset.getParameter<edm::InputTag>("MetLabel");
  thePhotonLabel     = pset.getParameter<edm::InputTag>("PhotonCollectionLabel");
  theCaloTowLabel    = pset.getParameter<edm::InputTag>("CaloTowerLabel");

  mudptmax           = pset.getParameter<double>("DimuonMaxdpt");
  mudphimin          = pset.getParameter<double>("DimuonMindphi");
  drisocalo          = pset.getParameter<double>("CaloTowerdR");

  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");

  edm::FileInPath myDataFile("FastSimulation/ProtonTaggers/data/acceptance_420_220.root");  
  std::string fullPath = myDataFile.fullPath();  
  std::cout << "Opening " << fullPath << std::endl;  
  TFile f(fullPath.c_str());  
  if (f.Get("description") != NULL)  
    std::cout << "Description found: " << f.Get("description")->GetTitle() << std::endl;  
    
  std::cout << "Reading acceptance tables " << std::endl;  
 
  helper420beam1.Init(f, "a420");  
  helper420beam2.Init(f, "a420_b2");  
  helper220beam1.Init(f, "a220");  
  helper220beam2.Init(f, "a220_b2");  
  helper420a220beam1.Init(f, "a420a220");  
  helper420a220beam2.Init(f, "a420a220_b2");  


  //  nEvt=0;
  MUONMAX=10;
  JETMAX=30;
  TRACKMAX=100;

  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

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

  thetree->Branch("nMuonCand",&nMuonCand,"nMuonCand/I");
  thetree->Branch("MuonCand_px",MuonCand_px,"MuonCand_px[nMuonCand]/D");
  thetree->Branch("MuonCand_py",MuonCand_py,"MuonCand_py[nMuonCand]/D");
  thetree->Branch("MuonCand_pz",MuonCand_pz,"MuonCand_pz[nMuonCand]/D");
  thetree->Branch("MuonCand_p",MuonCand_p,"MuonCand_p[nMuonCand]/D");
  thetree->Branch("MuonCand_pt",MuonCand_pt,"MuonCand_pt[nMuonCand]/D");
  thetree->Branch("MuonCand_eta",MuonCand_eta,"MuonCand_eta[nMuonCand]/D");
  thetree->Branch("MuonCand_phi",MuonCand_phi,"MuonCand_phi[nMuonCand]/D");

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

  thetree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  thetree->Branch("MuMu_dphi",&MuMu_dphi,"MuMu_dphi/D");

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");

  //  thetree->Branch("evweight",&evweight,"evweight/D"); 
}


GammaGammaMuMu::~GammaGammaMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaMuMu::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nMuonCand=0;
  nJetCand=0;
  nCaloCand=0;
  nTrackCand=0;

  nExtraCaloTowersE1=0;
  nExtraCaloTowersE2=0;
  nExtraCaloTowersE3=0; 
  nExtraCaloTowersE4=0; 
  nExtraCaloTowersE5=0; 
  nExtraCaloTowersEt0pt1=0; 
  nExtraCaloTowersEt0pt2=0;
  nExtraCaloTowersEt0pt5=0;  
  nExtraCaloTowersEt1=0; 
  nExtraCaloTowersEt2=0;  
  HitInZDC=0;
  HitInCastor=0;

  MuMu_mass = -1;
  MuMu_dphi = -1;

  bool passed = true;

 //using namespace edm;
  using reco::TrackCollection;

  //  Handle< double> weightHandle;
  //  event.getByLabel ("weight", weightHandle);
  //  evweight = * weightHandle;

  // Get the muon track collection from the event
  Handle<reco::MuonCollection> muons;
  event.getByLabel(theGLBMuonLabel, muons);
  reco::MuonCollection::const_iterator muon;

  if(muons->size() == 2)
    {
      for (muon = muons->begin(); muon != muons->end() && nMuonCand<MUONMAX; ++muon)
	{
	  MuonCand_p[nMuonCand]=muon->p();
	  MuonCand_px[nMuonCand]=muon->px();
	  MuonCand_py[nMuonCand]=muon->py();
	  MuonCand_pz[nMuonCand]=muon->pz();
	  MuonCand_pt[nMuonCand]=muon->pt();
	  MuonCand_eta[nMuonCand]=muon->eta();
	  MuonCand_phi[nMuonCand]=muon->phi();
	  MuonCand_charge[nMuonCand]=muon->charge();
	  nMuonCand++;
	}  

      // Calculate invariant mass and delta-phi
      if(MuonCand_charge[0]*MuonCand_charge[1]<0)
	{
	  double mass = pow(MuonCand_p[0]+MuonCand_p[1],2);
	  mass-=pow(MuonCand_px[0]+MuonCand_px[1],2);
	  mass-=pow(MuonCand_py[0]+MuonCand_py[1],2);
	  mass-=pow(MuonCand_pz[0]+MuonCand_pz[1],2);
	  MuMu_mass = sqrt(mass);

	  double dphi = fabs(MuonCand_phi[0]-MuonCand_phi[1]);
	  if(dphi < 3.14159)
	    MuMu_dphi = dphi;
	  else
	    MuMu_dphi = (2.0*3.14159)-dphi;
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

  // If this event contains a di-mu/e/gamma candidate, look at Jets & MET & CaloTowers & Tracks
  if(nMuonCand == 2)
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
	  
	  float calodr1 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[0])*(CaloTower_eta[nCaloCand]-MuonCand_eta[0])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[0])*(CaloTower_phi[nCaloCand]-MuonCand_phi[0])));
	  float calodr2 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[1])*(CaloTower_eta[nCaloCand]-MuonCand_eta[1])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[1])*(CaloTower_phi[nCaloCand]-MuonCand_phi[1])));
	  
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

      for(track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track)
	{
	  TrackCand_p[nTrackCand]=track->p();
	  TrackCand_px[nTrackCand]=track->px();
	  TrackCand_py[nTrackCand]=track->py();
	  TrackCand_pz[nTrackCand]=track->pz();
	  TrackCand_pt[nTrackCand]=track->pt();
	  TrackCand_eta[nTrackCand]=track->eta();
	  TrackCand_phi[nTrackCand]=track->phi();
	  TrackCand_charge[nTrackCand]=track->charge();
	  nTrackCand++; 
	}
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

      if(MCPar_pdgid == 2212 && MCPar_pz > 3000.0) 
        { 
          double MCPar_pt = sqrt(MCPar_px*MCPar_px + MCPar_py*MCPar_py); 
          double phi = p.phi(); 
          double mp = 0.938272029; 
          // ... compute kinimatical variable  
  
          float xi  = 1.0;    // fractional momentum loss  
          if (MCPar_pz>0)  
            xi -= MCPar_pz/7000.0;  
          else  
            xi += MCPar_pz/7000.0;  
  
          double t   = (-MCPar_pt*MCPar_pt - mp*mp*xi*xi) / (1-xi); // "t"  
 
          float acc420b1, acc220b1, acc420and220b1, acc420or220b1; // beam 1 (clockwise)  
          float acc420b2, acc220b2, acc420and220b2, acc420or220b2; // beam 2 (anti-clockwise)  
  
          acc420b1 = acc220b1 = acc420and220b1 = acc420or220b1 = 0;  
          acc420b2 = acc220b2 = acc420and220b2 = acc420or220b2 = 0;  
 
          if(MCPar_pz > 0) 
            { 
              acc420b1       = helper420beam1.GetAcceptance(t, xi, phi);  
              acc220b1       = helper220beam1.GetAcceptance(t, xi, phi);  
              acc420and220b1 = helper420a220beam1.GetAcceptance(t, xi, phi);  
              acc420or220b1  = acc420b1 + acc220b1 - acc420and220b1;  
            } 
          else 
            { 
              acc420b2       = helper420beam2.GetAcceptance(t, xi, phi);  
              acc220b2       = helper220beam2.GetAcceptance(t, xi, phi);  
              acc420and220b2 = helper420a220beam2.GetAcceptance(t, xi, phi);  
              acc420or220b2  = acc420b2 + acc220b2 - acc420and220b2;  
            } 
 
        } 

    }

  // Check for di-objects
  if(nMuonCand != 2)
    passed = false;
  else
    {
      if(MuMu_dphi < mudphimin) 
	passed = false;
      if(fabs(MuonCand_pt[0]-MuonCand_pt[1]) > mudptmax)
	passed = false;      
    }

  // "Exclusivity" cuts

  if(passed == true)
    thetree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuMu::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuMu::endJob() {
  thefile->Write();
  thefile->Close();
}
  
