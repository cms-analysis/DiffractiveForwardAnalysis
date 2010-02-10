 // -*- C++ -*-
//
// Package:    CollisionsMuMu
// Class:      CollisionsMuMu
// 
/**\class CollisionsMuMu CollisionsMuMu.cc GammaGammaLeptonLepton/CollisionsMuMu/src/CollisionsMuMu.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: CollisionsMuMu.cc,v 1.8 2010/02/09 15:45:24 jjhollar Exp $
//
//


// system include files
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h" 
#include "DataFormats/PatCandidates/interface/Electron.h" 
#include "DataFormats/PatCandidates/interface/Tau.h" 
#include "DataFormats/PatCandidates/interface/Photon.h" 
#include "DataFormats/PatCandidates/interface/MET.h" 

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h" 
 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "FWCore/ParameterSet/interface/ParameterDescriptionNode.h" 
#include "DataFormats/Common/interface/Ref.h"  
 
#include "DataFormats/Common/interface/TriggerResults.h"  
#include "FWCore/Framework/interface/TriggerNames.h"  
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
   
#include "FWCore/Framework/interface/ESHandle.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/EgammaCandidates/interface/Electron.h"  
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"   
#include "DataFormats/TrackReco/interface/Track.h"  
#include "DataFormats/TrackReco/interface/TrackFwd.h"  
#include "DataFormats/MuonReco/interface/Muon.h"  
#include "DataFormats/MuonReco/interface/MuonFwd.h"   
#include "DataFormats/MuonReco/interface/MuonSelectors.h"    
#include "DataFormats/METReco/interface/CaloMET.h"  
#include "DataFormats/METReco/interface/CaloMETFwd.h"   
#include "DataFormats/METReco/interface/CaloMETCollection.h"  
#include "DataFormats/EgammaCandidates/interface/Photon.h"  
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"  
#include "DataFormats/CaloTowers/interface/CaloTower.h"  
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"   
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Candidate/interface/Candidate.h"  
#include "DataFormats/Candidate/interface/CandidateFwd.h"   
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/CollisionsMuMu.h"

// user include files
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
// Muons:
#include <DataFormats/TrackReco/interface/Track.h>
// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

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

// for HF application
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/HcalTowerAlgo/src/HcalHardcodeGeometryData.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"  


// C++
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TMath.h>

using namespace HepMC;
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
CollisionsMuMu::CollisionsMuMu(const edm::ParameterSet& pset)
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

  //  nEvt=0;
  MUONMAX=10;
  JETMAX=30;
  TRACKMAX=100;
  CALOMAX=5000;

  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

  thetree->Branch("Run",&Run,"Run/I");
  thetree->Branch("LumiSection",&LumiSection,"LumiSection/I");
  thetree->Branch("GenProcessId",&GenProcessId,"GenProcessId/I");

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
  thetree->Branch("MuonCand_charge",MuonCand_charge,"MuonCand_charge[nMuonCand]/I");
  thetree->Branch("MuonCand_tmlsloosemuonid",MuonCand_tmlsloosemuonid,"MuonCand_tmlsloosemuonid[nMuonCand]/I");
  thetree->Branch("MuonCand_tm2dloosemuid",MuonCand_tm2dloosemuid,"MuonCand_tm2dloosemuid[nMuonCand]/I");
  thetree->Branch("MuonCand_arbmuid",MuonCand_arbmuid,"MuonCand_arbmuid[nMuonCand]/I");
  thetree->Branch("MuonCand_isglobal",MuonCand_isglobal,"MuonCand_isglobal[nMuonCand]/I");
  thetree->Branch("MuonCand_istracker",MuonCand_istracker,"MuonCand_istracker[nMuonCand]/I"); 
  thetree->Branch("MuonCand_isstandalone",MuonCand_isstandalone,"MuonCand_isstandalone[nMuonCand]/I"); 
  thetree->Branch("MuonCand_hcalisor3",MuonCand_hcalisor3,"MuonCand_hcalisor3[nMuonCand]/D"); 
  thetree->Branch("MuonCand_ecalisor3",MuonCand_ecalisor3,"MuonCand_ecalisor3[nMuonCand]/D");  
  thetree->Branch("MuonCand_trkisor3",MuonCand_trkisor3,"MuonCand_trkisor3[nMuonCand]/D");  
  thetree->Branch("MuonCand_hcalisor5",MuonCand_hcalisor5,"MuonCand_hcalisor5[nMuonCand]/D");  
  thetree->Branch("MuonCand_ecalisor5",MuonCand_ecalisor5,"MuonCand_ecalisor5[nMuonCand]/D");  
  thetree->Branch("MuonCand_trkisor5",MuonCand_trkisor5,"MuonCand_trkisor5[nMuonCand]/D"); 
  thetree->Branch("MuonCand_vtxx",MuonCand_vtxx,"MuonCand_vtxx[nMuonCand]/D"); 
  thetree->Branch("MuonCand_vtxy",MuonCand_vtxy,"MuonCand_vtxy[nMuonCand]/D");  
  thetree->Branch("MuonCand_vtxz",MuonCand_vtxz,"MuonCand_vtxz[nMuonCand]/D");  

  thetree->Branch("MuonCand_trkpt", MuonCand_trkpt, "MuonCand_trkpt[nMuonCand]/D");
  thetree->Branch("MuonCand_trketa", MuonCand_trketa, "MuonCand_trketa[nMuonCand]/D"); 
  thetree->Branch("MuonCand_trkphi", MuonCand_trkphi, "MuonCand_trkphi[nMuonCand]/D");   
  thetree->Branch("MuonCand_samuonpt", MuonCand_samuonpt, "MuonCand_samuonpt[nMuonCand]/D");
  thetree->Branch("MuonCand_samuoneta", MuonCand_samuoneta, "MuonCand_samuoneta[nMuonCand]/D");  
  thetree->Branch("MuonCand_samuonphi", MuonCand_samuonphi, "MuonCand_samuonphi[nMuonCand]/D"); 
  thetree->Branch("MuonCand_timein", MuonCand_timein, "MuonCand_timein[nMuonCand]/D");  
  thetree->Branch("MuonCand_timeout", MuonCand_timeout, "MuonCand_timeout[nMuonCand]/D");  
  thetree->Branch("MuonCand_timeinerr", MuonCand_timeinerr, "MuonCand_timeinerr[nMuonCand]/D");   
  thetree->Branch("MuonCand_timeouterr", MuonCand_timeouterr, "MuonCand_timeouterr[nMuonCand]/D");   

  thetree->Branch("nCaloCand",&nCaloCand,"nCaloCand/I");
  thetree->Branch("CaloTower_e",CaloTower_e,"CaloTower_e[nCaloCand]/D");
  thetree->Branch("CaloTower_emE",CaloTower_emE,"CaloTower_emE[nCaloCand]/D");
  thetree->Branch("CaloTower_hadE",CaloTower_hadE,"CaloTower_hadE[nCaloCand]/D");
  thetree->Branch("CaloTower_outE",CaloTower_outE,"CaloTower_outE[nCaloCand]/D");
  thetree->Branch("CaloTower_et",CaloTower_et,"CaloTower_et[nCaloCand]/D");
  thetree->Branch("CaloTower_eta",CaloTower_eta,"CaloTower_eta[nCaloCand]/D"); 
  thetree->Branch("CaloTower_phi",CaloTower_phi,"CaloTower_phi[nCaloCand]/D"); 
  thetree->Branch("CaloTower_dr",CaloTower_dr,"CaloTower_dr[nCaloCand]/D");

  thetree->Branch("CaloTower_ID",CaloTower_ID,"CaloTower_ID[nCaloCand]/I");
  thetree->Branch("CaloTower_x",CaloTower_x,"CaloTower_x[nCaloCand]/D");
  thetree->Branch("CaloTower_y",CaloTower_y,"CaloTower_y[nCaloCand]/D");
  thetree->Branch("CaloTower_z",CaloTower_z,"CaloTower_z[nCaloCand]/D");
  thetree->Branch("CaloTower_t",CaloTower_t,"CaloTower_t[nCaloCand]/D");

  thetree->Branch("CaloTower_badhcalcells",CaloTower_badhcalcells,"CaloTower_badhcalcells[nCaloCand]/I"); 
  thetree->Branch("CaloTower_problemhcalcells",CaloTower_problemhcalcells,"CaloTower_problemhcalcells[nCaloCand]/I"); 
  thetree->Branch("CaloTower_badecalcells",CaloTower_badecalcells,"CaloTower_badecalcells[nCaloCand]/I");  
  thetree->Branch("CaloTower_problemecalcells",CaloTower_problemecalcells,"CaloTower_problemecalcells[nCaloCand]/I");  

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

  thetree->Branch("nExtraCaloTowersE0hf", &nExtraCaloTowersE0hf, "nExtraCaloTowersE0hf/I");
  thetree->Branch("nExtraCaloTowersE1hf", &nExtraCaloTowersE1hf, "nExtraCaloTowersE1hf/I"); 
  thetree->Branch("nExtraCaloTowersE2hf", &nExtraCaloTowersE2hf, "nExtraCaloTowersE12hf/I"); 
  thetree->Branch("nExtraCaloTowersE3hf", &nExtraCaloTowersE3hf, "nExtraCaloTowersE13hf/I");  
  thetree->Branch("nExtraCaloTowersE4hf", &nExtraCaloTowersE4hf, "nExtraCaloTowersE14hf/I");  
  thetree->Branch("nExtraCaloTowersE5hf", &nExtraCaloTowersE5hf, "nExtraCaloTowersE15hf/I");  

  thetree->Branch("nExtraCaloTowersE1he", &nExtraCaloTowersE1he, "nExtraCaloTowersE1he/I"); 
  thetree->Branch("nExtraCaloTowersE2he", &nExtraCaloTowersE2he, "nExtraCaloTowersE2he/I");  
  thetree->Branch("nExtraCaloTowersE3he", &nExtraCaloTowersE3he, "nExtraCaloTowersE3he/I");  
  thetree->Branch("nExtraCaloTowersE4he", &nExtraCaloTowersE4he, "nExtraCaloTowersE4he/I");   
  thetree->Branch("nExtraCaloTowersE5he", &nExtraCaloTowersE5he, "nExtraCaloTowersE5he/I");   

  thetree->Branch("nExtraCaloTowersE2hb", &nExtraCaloTowersE2hb, "nExtraCaloTowersE2hb/I"); 
  thetree->Branch("nExtraCaloTowersE3hb", &nExtraCaloTowersE3hb, "nExtraCaloTowersE3hb/I");  
  thetree->Branch("nExtraCaloTowersE4hb", &nExtraCaloTowersE4hb, "nExtraCaloTowersE4hb/I");  
  thetree->Branch("nExtraCaloTowersE5hb", &nExtraCaloTowersE5hb, "nExtraCaloTowersE5hb/I");   

  thetree->Branch("nExtraCaloTowersEt0pt1",&nExtraCaloTowersEt0pt1,"nExtraCaloTowersEt0pt1/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt2",&nExtraCaloTowersEt0pt2,"nExtraCaloTowersEt0pt2/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt5",&nExtraCaloTowersEt0pt5,"nExtraCaloTowersEt0pt5/I"); 
  thetree->Branch("nExtraCaloTowersEt1",&nExtraCaloTowersEt1,"nExtraCaloTowersEt1/I"); 
  thetree->Branch("nExtraCaloTowersEt2",&nExtraCaloTowersEt2,"nExtraCaloTowersEt2/I"); 
  thetree->Branch("nExtraCaloTowersEt3",&nExtraCaloTowersEt3,"nExtraCaloTowersEt3/I");  
  thetree->Branch("nExtraCaloTowersEt4",&nExtraCaloTowersEt4,"nExtraCaloTowersEt4/I");  

  thetree->Branch("HF_plus_energy",&HF_plus_energy,"HF_plus_energy/D");
  thetree->Branch("HF_minus_energy",&HF_minus_energy,"HF_minus_energy/D");
  thetree->Branch("HF_plus_time",&HF_plus_time,"HF_plus_time/D");
  thetree->Branch("HF_minus_time",&HF_minus_time,"HF_minus_time/D");

  thetree->Branch("nTrackCand",&nTrackCand,"nTrackCand/I");
  thetree->Branch("nExtraTrackCand",&nExtraTrackCand,"nExtraTrackCand/I"); 
  thetree->Branch("TrackCand_px",TrackCand_px,"TrackCand_px[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_py",TrackCand_py,"TrackCand_py[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_pz",TrackCand_pz,"TrackCand_pz[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_p",TrackCand_p,"TrackCand_p[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_pt",TrackCand_pt,"TrackCand_pt[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_eta",TrackCand_eta,"TrackCand_eta[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_phi",TrackCand_phi,"TrackCand_phi[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_vtxdxyz",TrackCand_vtxdxyz,"TrackCand_vtxdxyz[nExtraTrackCand]/D");
  thetree->Branch("ClosestExtraTrack_vtxdxyz",&ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz/D");
  thetree->Branch("nTrackCandPassQCDCuts",&nTrackCandPassQCDCuts,"nTrackCandPassQCDCuts/I");

  thetree->Branch("nVertexCand",&nVertexCand,"nVertexCand/I");
  thetree->Branch("VertexCand_x",&VertexCand_x,"VertexCand_x[nVertexCand]/D");
  thetree->Branch("VertexCand_y",&VertexCand_y,"VertexCand_y[nVertexCand]/D");
  thetree->Branch("VertexCand_z",&VertexCand_z,"VertexCand_z[nVertexCand]/D");
  thetree->Branch("VertexCand_tracks",&VertexCand_tracks,"VertexCand_tracks[nVertexCand]/I");
  thetree->Branch("VertexCand_chi2",&VertexCand_chi2,"VertexCand_chi2[nVertexCand]/D");
  thetree->Branch("VertexCand_ndof",&VertexCand_ndof,"VertexCand_ndof[nVertexCand]/I");
  
  thetree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  thetree->Branch("MuMu_dphi",&MuMu_dphi,"MuMu_dphi/D");
  thetree->Branch("MuMu_vtxx",&MuMu_vtxx,"MuMu_vtxx/D");
  thetree->Branch("MuMu_vtxy",&MuMu_vtxy,"MuMu_vtxy/D"); 
  thetree->Branch("MuMu_vtxz",&MuMu_vtxz,"MuMu_vtxz/D"); 
  thetree->Branch("MuMu_vtxchi2dof",&MuMu_vtxchi2dof,"MuMu_vtxchi2dof/D");
  thetree->Branch("MuMu_vtxisvalid",&MuMu_vtxisvalid,"MuMu_vtxisvalid/I");

  thetree->Branch("MuMu_extratracks1mm",&MuMu_extratracks1mm,"MuMu_extratracks1mm/I");
  thetree->Branch("MuMu_extratracks2mm",&MuMu_extratracks2mm,"MuMu_extratracks2mm/I");
  thetree->Branch("MuMu_extratracks3mm",&MuMu_extratracks3mm,"MuMu_extratracks3mm/I");
  thetree->Branch("MuMu_extratracks4mm",&MuMu_extratracks4mm,"MuMu_extratracks4mm/I"); 
  thetree->Branch("MuMu_extratracks5mm",&MuMu_extratracks5mm,"MuMu_extratracks5mm/I");  
  thetree->Branch("MuMu_extratracks1cm",&MuMu_extratracks1cm,"MuMu_extratracks1cm/I"); 
  thetree->Branch("MuMu_extratracks2cm",&MuMu_extratracks2cm,"MuMu_extratracks2cm/I");
  thetree->Branch("MuMu_extratracks5cm",&MuMu_extratracks5cm,"MuMu_extratracks5cm/I"); 
  thetree->Branch("MuMu_extratracks10cm",&MuMu_extratracks10cm,"MuMu_extratracks10cm/I"); 

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");

  thetree->Branch("HLTMinBiasBSC",&HLTMinBiasBSC,"HLTMinBiasBSC/I");
  thetree->Branch("HLTMinBiasBSCOR",&HLTMinBiasBSCOR,"HLTMinBiasBSCOR/I");
  thetree->Branch("HLTMinBiasPixelSingleTrack",&HLTMinBiasPixelSingleTrack,"HLTMinBiasPixelSingleTrack/I");

  thetree->Branch("HLT2MuonNonIso",&HLT2MuonNonIso,"HLT2MuonNonIso/I");
  thetree->Branch("HLT1MuonPrescalePt3",&HLT1MuonPrescalePt3,"HLT1MuonPrescalePt3/I"); 
  thetree->Branch("L1TechnicalTriggers",L1TechnicalTriggers,"L1TechnicalTriggers[128]/I");
  thetree->Branch("BX",&BX,"BX/I");

  //  thetree->Branch("evweight",&evweight,"evweight/D"); 
}


CollisionsMuMu::~CollisionsMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CollisionsMuMu::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nMuonCand=0;
  nJetCand=0;
  nCaloCand=0;
  nTrackCand=0;
  nTrackCandPassQCDCuts = 0;
  nExtraTrackCand=0;
  nVertexCand=0;

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
  nExtraCaloTowersEt3=0;  
  nExtraCaloTowersEt4=0;   
  nExtraCaloTowersE0hf=0;
  nExtraCaloTowersE1hf=0;
  nExtraCaloTowersE2hf=0;  
  nExtraCaloTowersE3hf=0; 
  nExtraCaloTowersE4hf=0; 
  nExtraCaloTowersE5hf=0;   
  nExtraCaloTowersE1he=0; 
  nExtraCaloTowersE2he=0; 
  nExtraCaloTowersE3he=0;  
  nExtraCaloTowersE4he=0;  
  nExtraCaloTowersE5he=0;   
  nExtraCaloTowersE2hb=0; 
  nExtraCaloTowersE3hb=0; 
  nExtraCaloTowersE4hb=0;  
  nExtraCaloTowersE5hb=0;   

  HitInZDC=0;
  HitInCastor=0;

  MuMu_mass = -1;
  MuMu_dphi = -1;

  MuMu_extratracks1mm = 0;
  MuMu_extratracks2mm = 0;
  MuMu_extratracks3mm = 0;
  MuMu_extratracks4mm = 0;
  MuMu_extratracks5mm = 0;
  MuMu_extratracks1cm = 0; 
  MuMu_extratracks2cm = 0;
  MuMu_extratracks5cm = 0;
  MuMu_extratracks10cm = 0;

 //using namespace edm;
  using reco::TrackCollection;

  //  Handle< double> weightHandle;
  //  event.getByLabel ("weight", weightHandle);
  //  evweight = * weightHandle;

  BX = event.bunchCrossing();
  Run = event.id().run();
  LumiSection = event.luminosityBlock();
  
  GenProcessId = -1;
  Handle<HepMCProduct> mcevt;
  event.getByLabel(InputTag("generator"), mcevt);
  if(mcevt.isValid())
    {
      const HepMC::GenEvent* genEvent = mcevt->GetEvent() ;
      GenProcessId = genEvent->signal_process_id(); 
    }

  // L1 technical triggers
  edm::Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  edm::Handle<L1GlobalTriggerObjectMapRecord> L1GTOMRec;
  event.getByLabel(InputTag("gtDigis::RECO"), L1GTRR);
  event.getByLabel(InputTag("hltL1GtObjectMap::HLT"), L1GTOMRec);
  if (L1GTRR.isValid()) {
    cout << "L1GTRR is valid" << endl;
    DecisionWord gtDecisionWord = L1GTRR->decisionWord();
    const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = L1GTRR->technicalTriggerWord();
    const unsigned int numberTechnicalTriggerBits(technicalTriggerWordBeforeMask.size());
    for (unsigned int iBit = 0; iBit < numberTechnicalTriggerBits; ++iBit) {
      int techTrigger = (int) technicalTriggerWordBeforeMask.at(iBit);
      L1TechnicalTriggers[iBit] = techTrigger;
    }
  }

  edm::Handle<edm::TriggerResults> hltResults ; 
  event.getByLabel(InputTag("TriggerResults::HLT"),hltResults) ; 
  trigNames.init(*hltResults) ; 
  for (unsigned int i=0; i<trigNames.size(); i++)  
    //    {if(hltResults->accept(i)==1)} cout<<"bit "<<i<<" = \t"<<trigNames.triggerNames().at(i)<<" accepted "<<endl;}
    { 
      if ( trigNames.triggerNames().at(i) == "HLT_MinBiasBSC" )
        {
          if ( hltResults->accept(i) )
            HLTMinBiasBSC = 1;
          else
            HLTMinBiasBSC = 0;
        }
      if ( trigNames.triggerNames().at(i) == "HLT_MinBiasBSC_OR" )
        {
          if ( hltResults->accept(i) )
            HLTMinBiasBSCOR = 1;
          else
            HLTMinBiasBSCOR = 0;
        }
      if ( trigNames.triggerNames().at(i) == "HLT_MinBiasPixel_SingleTrack" )
        {
          if ( hltResults->accept(i) )
            HLTMinBiasPixelSingleTrack = 1;
          else
            HLTMinBiasPixelSingleTrack = 0;
        }
    }

  // Get the electron collection from the event
  edm::Handle<reco::GsfElectronCollection> electrons;
  event.getByLabel(thePixelGsfELabel,electrons);
//  const reco::PixelMatchGsfElectronCollection* electrons = pTracks.product();
  reco::GsfElectronCollection::const_iterator electron;
  //  cout << "Found " << electrons->size() << " electrons" << endl; 

  /*
  // PAT
  edm::Handle<edm::View<pat::Muon> > muons; 
  event.getByLabel(theGLBMuonLabel,muons); 
  edm::View<pat::Muon>::const_iterator muon;
  */
  // AOD
  Handle<reco::MuonCollection> muons;
  event.getByLabel(theGLBMuonLabel, muons);
  reco::MuonCollection::const_iterator muon;
  //  cout << "Found " << muons->size() << " muons" << endl; 


  // Get the electron collection from the event

  //  if(muons->size() == 2)
  if(muons->size() > 0)
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
	  MuonCand_vtxx[nMuonCand]=muon->vertex().x();
	  MuonCand_vtxy[nMuonCand]=muon->vertex().y(); 
          MuonCand_vtxz[nMuonCand]=muon->vertex().z(); 

	  // Muon ID
// 	  MuonCand_tmlsloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationLoose);
// 	  MuonCand_tm2dloosemuid[nMuonCand]=muon::isGoodMuon(*muon, muon::TM2DCompatibilityLoose);
// 	  MuonCand_arbmuid[nMuonCand]=muon::isGoodMuon(*muon, muon::AllArbitrated);
	  MuonCand_isglobal[nMuonCand]=muon->isGlobalMuon();
	  MuonCand_istracker[nMuonCand]=muon->isTrackerMuon(); 
	  MuonCand_isstandalone[nMuonCand]=muon->isStandAloneMuon(); 

	  // Isolation
	  MuonCand_hcalisor3[nMuonCand]=muon->isolationR03().hadEt; 
	  MuonCand_ecalisor3[nMuonCand]=muon->isolationR03().emEt;  
	  MuonCand_trkisor3[nMuonCand]=muon->isolationR03().nTracks;  
	  MuonCand_hcalisor5[nMuonCand]=muon->isolationR05().hadEt;  
	  MuonCand_ecalisor5[nMuonCand]=muon->isolationR05().emEt;   
	  MuonCand_trkisor5[nMuonCand]=muon->isolationR05().nTracks;   

	  if(MuonCand_isglobal[nMuonCand] || MuonCand_istracker[nMuonCand])
	    MuonCandTrack_p[nMuonCand] = muon->track()->p(); 
	  else
	    MuonCandTrack_p[nMuonCand] = muon->p();

	  if(muon->innerTrack().isNonnull())
	    {
	      MuonCand_trkpt[nMuonCand]=muon->innerTrack()->pt(); 
	      MuonCand_trketa[nMuonCand]=muon->innerTrack()->eta();  
	      MuonCand_trkphi[nMuonCand]=muon->innerTrack()->phi();  
	    }
	  if(muon->outerTrack().isNonnull())
	    {
	      MuonCand_samuonpt[nMuonCand]=muon->outerTrack()->pt();  
	      MuonCand_samuoneta[nMuonCand]=muon->outerTrack()->eta();  
	      MuonCand_samuonphi[nMuonCand]=muon->outerTrack()->phi();  
	    }

	  MuonCand_timeout[nMuonCand]=muon->time().timeAtIpInOut;
	  MuonCand_timein[nMuonCand]=muon->time().timeAtIpOutIn;
	  MuonCand_timeouterr[nMuonCand]=muon->time().timeAtIpInOutErr; 
          MuonCand_timeinerr[nMuonCand]=muon->time().timeAtIpOutInErr; 

	  nMuonCand++;
	}  

      // Calculate invariant mass and delta-phi
      if(nMuonCand == 2)
	{
	  //	  if(MuonCand_charge[0]*MuonCand_charge[1]<0)
	  //	    {
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
	  
	  //	    }
	}
    }


  // Get the Jet collection from the event
  /*
  // PAT
  edm::Handle<edm::View<pat::Jet> > jets; 
  event.getByLabel(theJetLabel,jets); 
  edm::View<pat::Jet>::const_iterator jet;
  */

  // AOD
  edm::Handle<reco::CaloJetCollection> pJets;
  event.getByLabel(theJetLabel,pJets);
  const reco::CaloJetCollection* jets = pJets.product();
  reco::CaloJetCollection::const_iterator jet;
  //  cout << "Found " << jets->size() << " jets" << endl;

  // Get the MET collection from the event

  /*
  // PAT
  edm::Handle<edm::View<pat::MET> > mets; 
  event.getByLabel(theMetLabel,mets); 
  edm::View<pat::MET>::const_iterator met;
  */

  // AOD
  edm::Handle<reco::CaloMETCollection> pMET;
  event.getByLabel(theMetLabel,pMET);
  const reco::CaloMETCollection* mets = pMET.product();
  reco::CaloMETCollection::const_iterator met;


  // Get the CaloTower collection from the event
  edm::Handle<CaloTowerCollection> caloTowers; 
  event.getByLabel(theCaloTowLabel,caloTowers); 
  const CaloTowerCollection* towers = caloTowers.product(); 
  CaloTowerCollection::const_iterator calo; 
  //  cout << "Found " << towers->size() << " towers" << endl;

  // Get the reconstructed hits from the frame.
  edm::Handle<HFRecHitCollection> hf;
  if( !event.getByLabel("hfreco",hf) ){
    std::cout << "Could not get rec hits! Tried with label: hfreco" << std::endl;
  }

  ESHandle<CaloGeometry> geometry ;
  iSetup.get<CaloGeometryRecord>().get(geometry);


  // Get the vertex collection from the event
  edm::Handle<reco::VertexCollection> recoVertexs;
//  event.getByLabel(recVertexLabel, recoVertexs);
  event.getByLabel(InputTag("offlinePrimaryVertices"), recoVertexs);
  const VertexCollection* vertexs = recoVertexs.product();
  VertexCollection::const_iterator vertex_i;
  //  cout << "Found " << vertexs->size() << " vertices" << endl;

  for (vertex_i = vertexs->begin(); vertex_i != vertexs->end(); vertex_i++){
	  VertexCand_x[nVertexCand] = vertex_i->x();
          VertexCand_y[nVertexCand] = vertex_i->y();
          VertexCand_z[nVertexCand] = vertex_i->z();
	  VertexCand_tracks[nVertexCand] = vertex_i->tracksSize();
	  VertexCand_chi2[nVertexCand] = vertex_i->chi2();
	  VertexCand_ndof[nVertexCand] = vertex_i->ndof();
	  nVertexCand++;	
  }


  // Get the track collection from the event
    edm::Handle<reco::TrackCollection> recoTracks;
//    event.getByLabel(recTrackLabel, recoTracks);
    event.getByLabel(InputTag("generalTracks"), recoTracks);
    const TrackCollection* tracks = recoTracks.product();
    TrackCollection::const_iterator track;
    //    cout << "Found " << tracks->size() << " tracks" << endl;

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
  if(nMuonCand > -1)
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

      for (calo = towers->begin(); calo != towers->end() && nCaloCand < CALOMAX; ++calo )
	{
	  CaloTower_e[nCaloCand]=calo->energy();
	  CaloTower_et[nCaloCand]=calo->et();
	  CaloTower_phi[nCaloCand]=calo->phi(); 
	  CaloTower_eta[nCaloCand]=calo->eta(); 
	  CaloTower_emE[nCaloCand]=calo->emEnergy();  //ecal
	  CaloTower_hadE[nCaloCand]=calo->hadEnergy();  //hcal
	  CaloTower_outE[nCaloCand]=calo->outerEnergy(); //ho

	  GlobalPoint emPosition=calo->emPosition();
          GlobalPoint hadPosition=calo->hadPosition();
	  //	  cout<<"em. = "<<calo->emEnergy()<<"\t had. = "<<calo->hadEnergy()<<"\t eta = "<<calo->eta()<<endl;
	  if(CaloTower_emE[nCaloCand]>0){ //if ECAL
		CaloTower_x[nCaloCand]=emPosition.x();
                CaloTower_y[nCaloCand]=emPosition.y();
                CaloTower_z[nCaloCand]=emPosition.z();
                CaloTower_t[nCaloCand]=calo->ecalTime();
		if(fabs(CaloTower_eta[nCaloCand])>=0   && fabs(CaloTower_eta[nCaloCand])<1.4) {CaloTower_ID[nCaloCand]=1;} //cout<<" -> EB"<<endl;}//EB
		if(fabs(CaloTower_eta[nCaloCand])>=1.4 && fabs(CaloTower_eta[nCaloCand])<2.95){CaloTower_ID[nCaloCand]=2;} //cout<<" -> EE"<<endl;}//EE
		if(fabs(CaloTower_eta[nCaloCand])>=2.95 && fabs(CaloTower_eta[nCaloCand])<5.2){CaloTower_ID[nCaloCand]=3;} //cout<<" -> EF"<<endl;}//EF ??
	  }

	  else if(CaloTower_hadE[nCaloCand]>0){
                CaloTower_x[nCaloCand]=hadPosition.x();
                CaloTower_y[nCaloCand]=hadPosition.y();
                CaloTower_z[nCaloCand]=hadPosition.z();
                CaloTower_t[nCaloCand]=calo->hcalTime();
		if(fabs(CaloTower_eta[nCaloCand])>=0   && fabs(CaloTower_eta[nCaloCand])<1.4) {CaloTower_ID[nCaloCand]=4;} // cout<<" -> HB"<<endl;}//HB
		if(fabs(CaloTower_eta[nCaloCand])>=1.4 && fabs(CaloTower_eta[nCaloCand])<2.95) {CaloTower_ID[nCaloCand]=5;} // cout<<" -> HE"<<endl;}//HE
		if(fabs(CaloTower_eta[nCaloCand])>=2.95 && fabs(CaloTower_eta[nCaloCand])<5.2) {CaloTower_ID[nCaloCand]=6;} // cout<<" -> HF"<<endl;}//HF ??
	  }
	  //	  else cout <<"PBLM ?????"<<endl;

	  CaloTower_badhcalcells[nCaloCand]=calo->numBadHcalCells();
	  CaloTower_problemhcalcells[nCaloCand]=calo->numProblematicHcalCells(); 
          CaloTower_badecalcells[nCaloCand]=calo->numBadEcalCells(); 
          CaloTower_problemecalcells[nCaloCand]=calo->numProblematicEcalCells();  

	  float calodr1 = 999.;
	  float calodr2 = 999.;
	  if(nMuonCand > 0)
	    calodr1 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[0])*(CaloTower_eta[nCaloCand]-MuonCand_eta[0])) + 
			   ((CaloTower_phi[nCaloCand]-MuonCand_phi[0])*(CaloTower_phi[nCaloCand]-MuonCand_phi[0])));
	  if(nMuonCand > 1)
	    calodr2 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[1])*(CaloTower_eta[nCaloCand]-MuonCand_eta[1])) + 
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
              if(CaloTower_et[nCaloCand] > 3.0)   
                nExtraCaloTowersEt3++;   
              if(CaloTower_et[nCaloCand] > 4.0)   
                nExtraCaloTowersEt4++;   

              if(CaloTower_e[nCaloCand] > 0.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)  
		nExtraCaloTowersE0hf++; 
              if(CaloTower_e[nCaloCand] > 1.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)   
		nExtraCaloTowersE1hf++; 
              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)   
		nExtraCaloTowersE2hf++;
	      if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)   
                nExtraCaloTowersE3hf++;  
              if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)    
                nExtraCaloTowersE4hf++;  
              if(CaloTower_e[nCaloCand] > 5.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)    
                nExtraCaloTowersE5hf++; 

              if(CaloTower_e[nCaloCand] > 1.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)      
		nExtraCaloTowersE1he++;  
              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)   
		nExtraCaloTowersE2he++;  
	      if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)   
		nExtraCaloTowersE3he++;   
	      if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)    
                nExtraCaloTowersE4he++;   
              if(CaloTower_e[nCaloCand] > 5.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)    
                nExtraCaloTowersE5he++;    

              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)   
		nExtraCaloTowersE2hb++;  
              if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)   
		nExtraCaloTowersE3hb++;  
              if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)   
		nExtraCaloTowersE4hb++;   
              if(CaloTower_e[nCaloCand] > 5.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)    
                nExtraCaloTowersE5hb++;    
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

      //HF treatment
      double sumE_plus(0), sumE_minus(0), time_plus(0), time_minus(0);
      for( unsigned int i=0; i<hf->size(); i++ )
	{
	  double energy = (*hf)[i].energy();
	  if(energy < 1.)continue;
	  double time = (*hf)[i].time();
	  HcalDetId cell( (*hf)[i].id());

         const CaloCellGeometry* cellGeometry =
	   geometry->getSubdetectorGeometry (cell)->getGeometry (cell);
         if( cellGeometry == 0 ) std::cout << "No cell geometry " << cell.rawId() << std::endl;
         double fEta = cellGeometry->getPosition().eta();
         //double fPhi = cellGeometry->getPosition().phi();

	 if(fEta > 0){sumE_plus+=energy; time_plus+=energy*time;}
	 if(fEta < 0){sumE_minus+=energy; time_minus+=energy*time;}
	}
      HF_plus_energy =sumE_plus;
      HF_minus_energy=sumE_minus;
      HF_plus_time = (sumE_plus > 0) ? time_plus/sumE_plus : -99;
      HF_minus_time= (sumE_minus > 0) ? time_minus/sumE_minus : -99;
    }

      for(track = tracks->begin(); track != tracks->end() && nExtraTrackCand<TRACKMAX; ++ track) 
        { 
//	  if(track->p() == MuonCandTrack_p[0] || track->p() == MuonCandTrack_p[1])
//	    continue;
	  
          TrackCand_p[nExtraTrackCand]=track->p(); 
          TrackCand_px[nExtraTrackCand]=track->px(); 
          TrackCand_py[nExtraTrackCand]=track->py(); 
          TrackCand_pz[nExtraTrackCand]=track->pz(); 
          TrackCand_pt[nExtraTrackCand]=track->pt(); 
          TrackCand_eta[nExtraTrackCand]=track->eta(); 
          TrackCand_phi[nExtraTrackCand]=track->phi(); 
          TrackCand_charge[nExtraTrackCand]=track->charge(); 

	  if((TrackCand_pt[nExtraTrackCand] > 0.5) && (fabs(TrackCand_eta[nExtraTrackCand]) < 2.0) && (fabs(track->vertex().z() - VertexCand_z[0]) < 1))
	    nTrackCandPassQCDCuts++;

          nExtraTrackCand++;  
        } 

      //    cout << "-----------" << endl;
    thetree->Fill();

}

void 
CollisionsMuMu::fillDescriptions(ConfigurationDescriptions & descriptions) { 
   
  descriptions.setComment("Exclusive dimuon cosmic ray EDAnalyzer."); 
   
  edm::ParameterSetDescription iDesc;   
 
  iDesc.add<edm::InputTag>("GlobalMuonCollectionLabel", edm::InputTag("muons"))->setComment("input muon collection"); 
  iDesc.add<edm::InputTag>("CaloTowerLabel", edm::InputTag("towerMaker"))->setComment("input calo tower collection");  
  iDesc.add<edm::InputTag>("RecoTrackLabel", edm::InputTag("cosmictrackfinderP5"))->setComment("input track collection");  
  iDesc.add<edm::InputTag>("JetCollectionLabel", edm::InputTag("sisCone5CaloJets"))->setComment("input jet collection");  
  iDesc.add<edm::InputTag>("ElectronCollectionLabel", edm::InputTag("pixelMatchGsfElectrons"))->setComment("input electron collection");  
  iDesc.add<edm::InputTag>("PhotonCollectionLabel", edm::InputTag("correctedPhotons"))->setComment("input photon collection");  
  iDesc.add<edm::InputTag>("MetLabel", edm::InputTag("met"))->setComment("input MET collection");    
  iDesc.add<double>("CaloTowerdR", 0.3)->setComment("Minimum delta-R to use for finding extra towers");   
  iDesc.add<double>("DimuonMindphi", 0.0)->setComment("Minimum delta-phi of dimuon pair");   
  iDesc.add<double>("DimuonMaxdpt", 2000.0)->setComment("Maximum delta-pT of dimuon pair");   
  iDesc.addOptionalUntracked<std::string>("outfilename", ("mumu.recolevel.root"))->setComment("output flat ntuple file name");   
 
  descriptions.add("ParameterDescriptionsForCollisionsMuMu", iDesc); 
} 


// ------------ method called once each job just before starting event loop  ------------
void 
CollisionsMuMu::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CollisionsMuMu::endJob() {
  thefile->Write();
  thefile->Close();
}
  
