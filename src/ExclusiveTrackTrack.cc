// -*- C++ -*-
//
// Package:    ExclusiveTrackTrack
// Class:      ExclusiveTrackTrack
// 
/**\class ExclusiveTrackTrack ExclusiveTrackTrack.cc GammaGammaLeptonLepton/ExclusiveTrackTrack/src/ExclusiveTrackTrack.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: ExclusiveTrackTrack.cc,v 1.14 2010/09/12 13:55:03 jjhollar Exp $
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
#include "DataFormats/HLTReco/interface/TriggerEvent.h" 
#include "FWCore/Common/interface/TriggerNames.h"   

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/RecoCandidate/interface/CaloRecHitCandidate.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

#include "DataFormats/CastorReco/interface/CastorTower.h"  
#include "DataFormats/HcalRecHit/interface/CastorRecHit.h"

#include "FWCore/Framework/interface/ESHandle.h" 
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h" 
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/CaloJet.h" 
#include "DataFormats/EgammaCandidates/interface/Electron.h" 
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h" 
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"   
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h" 
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" 
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"


#include "SimDataFormats/CrossingFrame/interface/MixCollection.h" // for PU

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ExclusiveTrackTrack.h"

// user include files
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
// Muons:
#include <DataFormats/TrackReco/interface/Track.h>
// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// dEdX
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

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
using namespace trigger;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ExclusiveTrackTrack::ExclusiveTrackTrack(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  recTrackLabel      = pset.getParameter<edm::InputTag>("RecoTrackLabel");
  theCaloTowLabel    = pset.getParameter<edm::InputTag>("CaloTowerLabel");
  recCastorTowerLabel = pset.getParameter<edm::InputTag>("CastorTowerLabel");
  recCastorRecHitsLabel = pset.getParameter<edm::InputTag>("CastorRecHitsLabel");
  recZDCRecHitsLabel = pset.getParameter<edm::InputTag>("ZDCRecHitsLabel");
  drisocalo          = pset.getParameter<double>("CaloTowerdR");
  fillallmc          = pset.getParameter<bool>("FillAllMCParticles"); 
  requiretwotracks   = pset.getParameter<bool>("RequireTwoTracks");

  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");

  //  nEvt=0;
  TRACKMAX=100;
  MCPARMAX=500;

  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

  thetree->Branch("GenProcessId",&GenProcessId,"GenProcessId/I");
  thetree->Branch("GenHasRho",&GenHasRho,"GenHasRho/I");

  thetree->Branch("Run",&Run,"Run/I");
  thetree->Branch("LumiSection",&LumiSection,"LumiSection/I");
  thetree->Branch("BX",&BX,"BX/I");

  thetree->Branch("L1TechnicalTriggers",L1TechnicalTriggers,"L1TechnicalTriggers[128]/I");
  thetree->Branch("HLTZeroBiasPixelSingleTrack",&HLTZeroBiasPixelSingleTrack,"HLTZeroBiasPixelSingleTrack/I");
  thetree->Branch("HLTZeroBias",&HLTZeroBias,"HLTZeroBias/I");
  thetree->Branch("HLT_L1_BscMinBiasOR_BptxPlusORMinus",&HLT_L1_BscMinBiasOR_BptxPlusORMinus,"HLT_L1_BscMinBiasOR_BptxPlusORMinus/I");
  thetree->Branch("HLTPhysicsDeclared",&HLTPhysicsDeclared,"HLTPhysicsDeclared/I");


  thetree->Branch("nTrackCand",&nTrackCand,"nTrackCand/I");
  thetree->Branch("TrackCand_px",TrackCand_px,"TrackCand_px[nTrackCand]/D");
  thetree->Branch("TrackCand_py",TrackCand_py,"TrackCand_py[nTrackCand]/D");
  thetree->Branch("TrackCand_pz",TrackCand_pz,"TrackCand_pz[nTrackCand]/D");
  thetree->Branch("TrackCand_p",TrackCand_p,"TrackCand_p[nTrackCand]/D");
  thetree->Branch("TrackCand_e",TrackCand_e,"TrackCand_e[nTrackCand]/D");
  thetree->Branch("TrackCand_pt",TrackCand_pt,"TrackCand_pt[nTrackCand]/D");
  thetree->Branch("TrackCand_eta",TrackCand_eta,"TrackCand_eta[nTrackCand]/D");
  thetree->Branch("TrackCand_phi",TrackCand_phi,"TrackCand_phi[nTrackCand]/D");
  thetree->Branch("TrackCand_charge",TrackCand_charge,"TrackCand_charge[nTrackCand]/I");
  thetree->Branch("TrackCand_chi2",TrackCand_chi2,"TrackCand_chi2[nTrackCand]/D");
  thetree->Branch("TrackCand_ndof",TrackCand_ndof,"TrackCand_ndof[nTrackCand]/D");
  thetree->Branch("TrackCand_purity",TrackCand_purity,"TrackCand_purity[nTrackCand]/I"); 
  thetree->Branch("TrackCand_dEdx",TrackCand_dEdx,"TrackCand_dEdx[nTrackCand]/D");
  thetree->Branch("TrackCand_nSaturated",TrackCand_nSaturated,"TrackCand_nSaturated[nTrackCand]/I");
  thetree->Branch("TrackCand_nMeasurements",TrackCand_nMeasurements,"TrackCand_nMeasurements[nTrackCand]/I");

  thetree->Branch("nCaloCand",&nCaloCand,"nCaloCand/I");
  thetree->Branch("CaloTower_e",CaloTower_e,"CaloTower_e[nCaloCand]/D");
  thetree->Branch("CaloTower_et",CaloTower_et,"CaloTower_et[nCaloCand]/D");
  thetree->Branch("CaloTower_eta",CaloTower_eta,"CaloTower_eta[nCaloCand]/D"); 
  thetree->Branch("CaloTower_phi",CaloTower_phi,"CaloTower_phi[nCaloCand]/D"); 
  thetree->Branch("CaloTower_dr",CaloTower_dr,"CaloTower_dr[nCaloCand]/D");
  thetree->Branch("CaloTower_emE",CaloTower_emE,"CaloTower_emE[nCaloCand]/D"); 
  thetree->Branch("CaloTower_hadE",CaloTower_hadE,"CaloTower_hadE[nCaloCand]/D"); 
  thetree->Branch("CaloTower_outE",CaloTower_outE,"CaloTower_outE[nCaloCand]/D"); 
  thetree->Branch("CaloTower_ID",CaloTower_ID,"CaloTower_ID[nCaloCand]/I"); 
  thetree->Branch("HighestCaloTower_e",&HighestCaloTower_e,"HighestCaloTower_e/D");
  thetree->Branch("HighestCaloTower_eta",&HighestCaloTower_eta,"HighestCaloTower_eta/D");
  thetree->Branch("HighestCaloTower_phi",&HighestCaloTower_phi,"HighestCaloTower_phi/D"); 
  thetree->Branch("HighestCaloTower_dr",&HighestCaloTower_dr,"HighestCaloTower_dr/D");
  thetree->Branch("HighestEtCaloTower_et",&HighestEtCaloTower_et,"HighestEtCaloTower_et/D");
  thetree->Branch("HighestEtCaloTower_eta",&HighestEtCaloTower_eta,"HighestEtCaloTower_eta/D");
  thetree->Branch("HighestEtCaloTower_phi",&HighestEtCaloTower_phi,"HighestEtCaloTower_phi/D"); 
  thetree->Branch("HighestEtCaloTower_dr",&HighestEtCaloTower_dr,"HighestEtCaloTower_dr/D");
  thetree->Branch("SumCalo_e",&SumCalo_e,"SumCalo_e/D");
  thetree->Branch("SumHFPlus_e",&SumHFPlus_e,"SumHFPlus_e/D");
  thetree->Branch("SumHFMinus_e",&SumHFMinus_e,"SumHFMinus_e/D");


  thetree->Branch("nCastorTowerCand",&nCastorTowerCand,"nCastorTowerCand/I");
  thetree->Branch("CastorTower_e",CastorTower_e,"CastorTower_e[nCastorTowerCand]/D");
  thetree->Branch("CastorTower_eta",CastorTower_eta,"CastorTower_eta[nCastorTowerCand]/D");
  thetree->Branch("CastorTower_phi",CastorTower_phi,"CastorTower_phi[nCastorTowerCand]/D");
  thetree->Branch("CastorTower_emratio",CastorTower_emratio,"CastorTower_emratio[nCastorTowerCand]/D");
  thetree->Branch("HighestCastorTowerFwd_e",&HighestCastorTowerFwd_e,"HighestCastorTowerFwd_e/D");
  thetree->Branch("HighestCastorTowerBwd_e",&HighestCastorTowerBwd_e,"HighestCastorTowerBwd_e/D");
  thetree->Branch("SumCastorFwd_e",&SumCastorFwd_e,"SumCastorFwd_e/D");
  thetree->Branch("SumCastorBwd_e",&SumCastorBwd_e,"SumCastorBwd_e/D");

  thetree->Branch("nZDChitCand", &nZDChitCand, "nZDChitCand/I");
  thetree->Branch("ZDChit_section", ZDChit_section, "ZDChit_section[nZDChitCand]/I");
  thetree->Branch("ZDChit_energy", ZDChit_energy, "ZDChit_energy[nZDChitCand]/D");
  thetree->Branch("ZDChit_time", ZDChit_time, "ZDChit_time[nZDChitCand]/D");
  thetree->Branch("ZDChit_side", ZDChit_side, "ZDChit_side[nZDChitCand]/I");
  thetree->Branch("ZDCsumEMplus", &ZDCsumEMplus, "ZDCsumEMplus/D");
  thetree->Branch("ZDCsumHADplus", &ZDCsumHADplus, "ZDCsumHADplus/D");
  thetree->Branch("ZDCsumEMminus", &ZDCsumEMminus, "ZDCsumEMminus/D");
  thetree->Branch("ZDCsumHADminus", &ZDCsumHADminus, "ZDCsumHADminus/D");
  thetree->Branch("CASTORsumRecHitsE", &CASTORsumRecHitsE, "CASTORsumRecHitsE/D");

  thetree->Branch("nExtraCaloTowersE1",&nExtraCaloTowersE1,"nExtraCaloTowersE1/I");
  thetree->Branch("nExtraCaloTowersE2",&nExtraCaloTowersE2,"nExtraCaloTowersE2/I");
  thetree->Branch("nExtraCaloTowersE3",&nExtraCaloTowersE3,"nExtraCaloTowersE3/I"); 
  thetree->Branch("nExtraCaloTowersE4",&nExtraCaloTowersE4,"nExtraCaloTowersE4/I"); 
  thetree->Branch("nExtraCaloTowersE5",&nExtraCaloTowersE5,"nExtraCaloTowersE5/I"); 
  thetree->Branch("nExtraCaloTowersE6",&nExtraCaloTowersE6,"nExtraCaloTowersE6/I");    
  thetree->Branch("nExtraCaloTowersE7",&nExtraCaloTowersE7,"nExtraCaloTowersE7/I");    
  thetree->Branch("nExtraCaloTowersE8",&nExtraCaloTowersE8,"nExtraCaloTowersE8/I");    
  thetree->Branch("nExtraCaloTowersE9",&nExtraCaloTowersE9,"nExtraCaloTowersE9/I");    

  thetree->Branch("nExtraCaloTowersE0pt6eb",&nExtraCaloTowersE0pt6eb, "nExtraCaloTowersE0pt6eb/I"); 
  thetree->Branch("nExtraCaloTowersE2pt45ee", &nExtraCaloTowersE2pt45ee, "nExtraCaloTowersE2pt45ee/I"); 
  thetree->Branch("nExtraCaloTowersE1pt25hb", &nExtraCaloTowersE1pt25hb, "nExtraCaloTowersE1pt25hb/I"); 
  thetree->Branch("nExtraCaloTowersE1pt9he", &nExtraCaloTowersE1pt9he, "nExtraCaloTowersE1pt9he/I"); 
  thetree->Branch("nExtraCaloTowersE4pt5hfp", &nExtraCaloTowersE4pt5hfp, "nExtraCaloTowersE4pt5hfp/I"); 
  thetree->Branch("nExtraCaloTowersE4pt0hfm", &nExtraCaloTowersE4pt0hfm, "nExtraCaloTowersE4pt0hfm/I");  

  thetree->Branch("nExtraCaloTowersEt0pt1",&nExtraCaloTowersEt0pt1,"nExtraCaloTowersEt0pt1/I");  
  thetree->Branch("nExtraCaloTowersEt0pt2",&nExtraCaloTowersEt0pt2,"nExtraCaloTowersEt0pt2/I");  
  thetree->Branch("nExtraCaloTowersEt0pt5",&nExtraCaloTowersEt0pt5,"nExtraCaloTowersEt0pt5/I");  
  thetree->Branch("nExtraCaloTowersEt1",&nExtraCaloTowersEt1,"nExtraCaloTowersEt1/I");  
  thetree->Branch("nExtraCaloTowersEt2",&nExtraCaloTowersEt2,"nExtraCaloTowersEt2/I");  
  thetree->Branch("nExtraCaloTowersEt3",&nExtraCaloTowersEt3,"nExtraCaloTowersEt3/I");  
  thetree->Branch("nExtraCaloTowersEt4",&nExtraCaloTowersEt4,"nExtraCaloTowersEt4/I");  

  thetree->Branch("TrTr_mass",&TrTr_mass,"TrTr_mass/D");
  thetree->Branch("TrTr_dphi",&TrTr_dphi,"TrTr_dphi/D");
  thetree->Branch("TrTr_dpt",&TrTr_dpt,"TrTr_dpt/D");
  thetree->Branch("TrTr_pt",&TrTr_pt,"TrTr_pt/D");
  thetree->Branch("TrTr_eta",&TrTr_eta,"TrTr_eta/D");
  thetree->Branch("TrTr_Kalmanvtxx",&TrTr_Kalmanvtxx,"TrTr_Kalmanvtxx/D");  
  thetree->Branch("TrTr_Kalmanvtxy",&TrTr_Kalmanvtxy,"TrTr_Kalmanvtxy/D");  
  thetree->Branch("TrTr_Kalmanvtxz",&TrTr_Kalmanvtxz,"TrTr_Kalmanvtxz/D");  
  thetree->Branch("TrTr_Kalmanvtxchi2dof",&TrTr_Kalmanvtxchi2dof,"TrTr_Kalmanvtxchi2dof/D");  
  thetree->Branch("TrTr_Kalmanvtxisvalid",&TrTr_Kalmanvtxisvalid,"TrTr_Kalmanvtxisvalid/I");  

  thetree->Branch("TrTr_mass_khyp",&TrTr_mass_khyp,"TrTr_mass_khyp/D");
  thetree->Branch("TrTr_pt_khyp",&TrTr_pt_khyp,"TrTr_pt_khyp/D");
  thetree->Branch("TrTr_eta_khyp",&TrTr_eta_khyp,"TrTr_eta_khyp/D");

  thetree->Branch("nVertexCand",&nVertexCand,"nVertexCand/I");
  thetree->Branch("VertexCand_x",&VertexCand_x,"VertexCand_x[nVertexCand]/D");
  thetree->Branch("VertexCand_y",&VertexCand_y,"VertexCand_y[nVertexCand]/D");
  thetree->Branch("VertexCand_z",&VertexCand_z,"VertexCand_z[nVertexCand]/D");
  thetree->Branch("VertexCand_tracks",&VertexCand_tracks,"VertexCand_tracks[nVertexCand]/I");
  thetree->Branch("VertexCand_chi2",&VertexCand_chi2,"VertexCand_chi2[nVertexCand]/D");
  thetree->Branch("VertexCand_ndof",&VertexCand_ndof,"VertexCand_ndof[nVertexCand]/D");

  if(fillallmc == true) 
    { 
      thetree->Branch("nMCPar",&nMCPar,"nMCPar/I");  
      thetree->Branch("MCPar_px",MCPar_px,"MCPar_px[nMCPar]/D");  
      thetree->Branch("MCPar_py",MCPar_py,"MCPar_py[nMCPar]/D");  
      thetree->Branch("MCPar_pz",MCPar_pz,"MCPar_pz[nMCPar]/D");  
      thetree->Branch("MCPar_eta",MCPar_eta,"MCPar_eta[nMCPar]/D");  
      thetree->Branch("MCPar_phi",MCPar_phi,"MCPar_phi[nMCPar]/D");  
      thetree->Branch("MCPar_pdgid",MCPar_pdgid,"MCPar_pdgid[nMCPar]/I"); 
      thetree->Branch("MCPar_status",MCPar_status,"MCPar_status[nMCPar]/I"); 
    } 

}


ExclusiveTrackTrack::~ExclusiveTrackTrack()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ExclusiveTrackTrack::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nCaloCand=0;
  nTrackCand=0;
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
  nExtraCaloTowersE0pt6eb=0; 
  nExtraCaloTowersE2pt45ee=0; 
  nExtraCaloTowersE1pt25hb=0;  
  nExtraCaloTowersE1pt9he=0;  
  nExtraCaloTowersE4pt5hfp=0; 
  nExtraCaloTowersE4pt0hfm=0;  
  nCastorTowerCand=0;
  nZDChitCand=0;
  ZDCsumHADminus=0;
  ZDCsumEMminus=0;
  ZDCsumHADplus=0;
  ZDCsumEMplus=0;
  CASTORsumRecHitsE=0;

  nCastorTowerCand=0;

  double highestcastortowerfwd = -999.0;
  double highestcastortowerbwd = -999.0;
  double totalecastorfwd = 0.0;
  double totalecastorbwd = 0.0;

  SumHFPlus_e=0.0;
  SumHFMinus_e=0.0;

  GenProcessId = -1;
  Handle<HepMCProduct> mcevt;
  event.getByLabel(InputTag("generator"), mcevt);
  if(mcevt.isValid())
    {
      const HepMC::GenEvent* genEvent = mcevt->GetEvent() ;
      GenProcessId = genEvent->signal_process_id();
    }

  GenHasRho = 0;
  Handle<GenParticleCollection> genParticles;
  event.getByLabel( InputTag("genParticles"), genParticles );
  if(genParticles.isValid())
    { 
      nMCPar=0;

      for ( size_t i = 0; i < genParticles->size(); ++ i )
	{
	  const Candidate & p = (*genParticles)[ i ];
	  int MCPar_rhocandid=p.pdgId();
	  if(MCPar_rhocandid == 113)
	    GenHasRho = 1;

	  if(fillallmc == true)
	    {
	      MCPar_pdgid[i]=p.pdgId(); 
	      MCPar_eta[i]=p.eta(); 
	      MCPar_phi[i]=p.phi();
	      MCPar_px[i]=p.px(); 
	      MCPar_py[i]=p.py(); 
	      MCPar_pz[i]=p.pz(); 
	      MCPar_mass[i]=p.mass(); 
	      MCPar_status[i]=p.status();
	      nMCPar++;
	    } 
	}
    }

  BX = event.bunchCrossing();
  Run = event.id().run();
  LumiSection = event.luminosityBlock();

  // L1 technical triggers
  edm::Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  edm::Handle<L1GlobalTriggerObjectMapRecord> L1GTOMRec;
  event.getByLabel(InputTag("gtDigis::RECO"), L1GTRR);
  event.getByLabel(InputTag("hltL1GtObjectMap::HLT"), L1GTOMRec);
  if (L1GTRR.isValid()) {
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
  //  trigNames.init(*hltResults) ;
  const edm::TriggerNames & trigNames = event.triggerNames(*hltResults); 
  for (unsigned int i=0; i<trigNames.size(); i++)
    //    {if(hltResults->accept(i)==1)} cout<<"bit "<<i<<" = \t"<<trigNames.triggerNames().at(i)<<" accepted "<<endl;}
    {
      if ( trigNames.triggerNames().at(i) == "HLT_ZeroBias" )
        {
          if ( hltResults->accept(i) )
            HLTZeroBias = 1;
          else
            HLTZeroBias = 0;
        }
      if ( trigNames.triggerNames().at(i) == "HLT_ZeroBiasPixel_SingleTrack" )
        {
          if ( hltResults->accept(i) )
            HLTZeroBiasPixelSingleTrack = 1;
          else
            HLTZeroBiasPixelSingleTrack = 0;
        }
      if ( trigNames.triggerNames().at(i) == "HLT_L1_BscMinBiasOR_BptxPlusORMinus" )
	{
	  if ( hltResults->accept(i) ) 
	    HLT_L1_BscMinBiasOR_BptxPlusORMinus = 1;
          else 
	    HLT_L1_BscMinBiasOR_BptxPlusORMinus = 0;
	}
      if ( trigNames.triggerNames().at(i) == "HLT_PhysicsDeclared" )
        {
          if ( hltResults->accept(i) )
            HLTPhysicsDeclared = 1;
          else
            HLTPhysicsDeclared = 0;
        }
    }

  edm::Handle<reco::VertexCollection> recoVertexs;
  event.getByLabel(InputTag("offlinePrimaryVertices"), recoVertexs);
  const VertexCollection* vertexs = recoVertexs.product();
  VertexCollection::const_iterator vertex_i;

  for (vertex_i = vertexs->begin(); vertex_i != vertexs->end(); vertex_i++){
    VertexCand_x[nVertexCand] = vertex_i->x();
    VertexCand_y[nVertexCand] = vertex_i->y();
    VertexCand_z[nVertexCand] = vertex_i->z();
    VertexCand_tracks[nVertexCand] = vertex_i->tracksSize();
    VertexCand_chi2[nVertexCand] = vertex_i->chi2();
    VertexCand_ndof[nVertexCand] = vertex_i->ndof();
    nVertexCand++;
  }


  TrTr_mass = -1;
  TrTr_dphi = -1;
  TrTr_dpt = -999.;
  TrTr_pt = -999.;

  bool passed = true;

 //using namespace edm;
  using reco::TrackCollection;

  // Now do vertexing and track counting
  edm::ESHandle<TransientTrackBuilder> theVtx;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theVtx);
  vector<TransientTrack> transtrks; 
  reco::TrackCollection * pitrks = new reco::TrackCollection;
  
  // Get the track collection from the event
  edm::Handle<reco::TrackCollection> recoTracks;
  event.getByLabel(recTrackLabel, recoTracks);
  const TrackCollection* tracks = recoTracks.product();
  TrackCollection::const_iterator track;

  // Get the dE/dX map from the event
  edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle;
  event.getByLabel("dedxHarmonic2", dEdxTrackHandle);
  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxTrackHandle.product();

  if((tracks->size() == 2) || (requiretwotracks == false))
    {
      for ( track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++track )
	{
	  reco::TrackRef trackref  = reco::TrackRef( recoTracks, nTrackCand );
	  //You can access dE/dx Estimation of your track with:
	  TrackCand_dEdx[nTrackCand] = dEdxTrack[trackref].dEdx();
	  TrackCand_nSaturated[nTrackCand] = dEdxTrack[trackref].numberOfSaturatedMeasurements();
	  TrackCand_nMeasurements[nTrackCand] = dEdxTrack[trackref].numberOfMeasurements();

	  TrackCand_pt[nTrackCand]=track->pt();
	  TrackCand_px[nTrackCand]=track->px();
	  TrackCand_py[nTrackCand]=track->py();
	  TrackCand_pz[nTrackCand]=track->pz();
	  TrackCand_p[nTrackCand]=track->p();
	  TrackCand_phi[nTrackCand]=track->phi();
	  TrackCand_eta[nTrackCand]=track->eta();
	  TrackCand_charge[nTrackCand]=track->charge(); 
	  TrackCand_chi2[nTrackCand]=track->chi2();
	  TrackCand_ndof[nTrackCand]=track->ndof();
          TrackCand_purity[nTrackCand]=track->quality(TrackBase::highPurity); 

	  pitrks->push_back( *track );
	  TransientTrack tmptrk = (*theVtx).build( *track );
	  transtrks.push_back( tmptrk );

	  nTrackCand++;
	}

      // Calculate invariant mass and delta-phi
      //      if(TrackCand_charge[0]*TrackCand_charge[1]<0)
      if(nTrackCand == 2)
	{
	  TLorentzVector pi1, pi2, rho;
	  pi1.SetXYZM(TrackCand_px[0], TrackCand_py[0], TrackCand_pz[0], 0.1396);
	  pi2.SetXYZM(TrackCand_px[1], TrackCand_py[1], TrackCand_pz[1], 0.1396);
	  rho = pi1 + pi2;

	  TrTr_mass = rho.M();
	  TrTr_pt = rho.Pt();
	  TrTr_eta = rho.Eta();

          TLorentzVector k1, k2, phi;
          k1.SetXYZM(TrackCand_px[0], TrackCand_py[0], TrackCand_pz[0], 0.492);
          k2.SetXYZM(TrackCand_px[1], TrackCand_py[1], TrackCand_pz[1], 0.492);
          phi = k1 + k2;

          TrTr_mass_khyp = phi.M();
          TrTr_pt_khyp = phi.Pt();
          TrTr_eta_khyp = phi.Eta();
	  

	  TrTr_dpt = TrackCand_pt[0]-TrackCand_pt[1];
	  double dphi = fabs(TrackCand_phi[0]-TrackCand_phi[1]);
	  if(dphi < 3.14159)
	    TrTr_dphi = dphi;
	  else
	    TrTr_dphi = (2.0*3.14159)-dphi;

          KalmanVertexFitter fitter(true);
          TransientVertex pipiVertex = fitter.vertex(transtrks);
          if(pipiVertex.isValid())
            {
              TrTr_Kalmanvtxx = pipiVertex.position().x();
              TrTr_Kalmanvtxy = pipiVertex.position().y();
              TrTr_Kalmanvtxz = pipiVertex.position().z();
              TrTr_Kalmanvtxchi2dof = pipiVertex.normalisedChiSquared();
              TrTr_Kalmanvtxisvalid = 1;
            }
          else
            {
              TrTr_Kalmanvtxx = 99;
              TrTr_Kalmanvtxy = 99;
              TrTr_Kalmanvtxz = 99;
              TrTr_Kalmanvtxchi2dof = 9999;
	      TrTr_Kalmanvtxisvalid = 0;
            }
	}
    }

  // Get the CaloTower collection from the event 
  edm::Handle<CaloTowerCollection> caloTowers;  
  event.getByLabel(theCaloTowLabel,caloTowers);  
  const CaloTowerCollection* towers = caloTowers.product();  
  CaloTowerCollection::const_iterator calo;  

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
     if((nTrackCand == 2) || (requiretwotracks == false))
    {
      for (calo = towers->begin(); calo != towers->end(); ++calo )
	{
	  CaloTower_e[nCaloCand]=calo->energy(); 
	  CaloTower_et[nCaloCand]=calo->et();
	  CaloTower_phi[nCaloCand]=calo->phi(); 
	  CaloTower_eta[nCaloCand]=calo->eta(); 
          CaloTower_emE[nCaloCand]=calo->emEnergy(); 
          CaloTower_hadE[nCaloCand]=calo->hadEnergy();  //hcal 
          CaloTower_outE[nCaloCand]=calo->outerEnergy(); //ho 
          GlobalPoint emPosition=calo->emPosition(); 
          GlobalPoint hadPosition=calo->hadPosition(); 
          CaloTowerDetId idTower=calo->id(); 
           
          size_t numRecHits = calo->constituentsSize(); 
          bool isEB(false),isEE(false),isHB(false),isHE(false),isHF(false),isHO(false); 
          for (size_t j = 0; j < numRecHits; j++) { 
            DetId RecHitDetID=calo->constituent(j); 
	    DetId::Detector DetNum=RecHitDetID.det(); 
            if( DetNum == DetId::Hcal){ 
              HcalDetId HcalID = RecHitDetID; 
              int HcalNum =  HcalID.subdetId(); 
              if(HcalNum == HcalForward ){isHF=true;} 
              if(HcalNum == HcalBarrel ) {isHB=true;} 
              if(HcalNum == HcalEndcap ) {isHE=true;} 
              if(HcalNum == HcalOuter ) {isHO=true;} 
            } 
            if( DetNum == DetId::Ecal){ 
              int EcalNum = RecHitDetID.subdetId(); 
              if(EcalNum==1){isEB=true;} 
              if(EcalNum==2){isEE=true;} 
            } 
          } 
          if(isHF&&!isHB&&!isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=1;} //HF 
          else if(!isHF&&!isHB&&!isHE&&isEB&&!isEE){CaloTower_ID[nCaloCand]=2;} //EB 
          else if(!isHF&&!isHB&&!isHE&&!isEB&&isEE){CaloTower_ID[nCaloCand]=3;} //EE 
          else if(!isHF&&isHB&&!isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=4;} //HB 
          else if(!isHF&&!isHB&&isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=5;} //HE 
          else if(!isHF&&isHB&&!isHE&&isEB&&!isEE){CaloTower_ID[nCaloCand]=6;} // HB+EB 
          else if(!isHF&&!isHB&&isHE&&isEB&&!isEE){CaloTower_ID[nCaloCand]=7;} // HE+EB 
          else if(!isHF&&!isHB&&isHE&&!isEB&&isEE){CaloTower_ID[nCaloCand]=8;} // HE+EE 
          else if(!isHF&&isHB&&isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=9;} // HB+HE 

	  float calodr1 = sqrt(((CaloTower_eta[nCaloCand]-TrackCand_eta[0])*(CaloTower_eta[nCaloCand]-TrackCand_eta[0])) + 
			       ((CaloTower_phi[nCaloCand]-TrackCand_phi[0])*(CaloTower_phi[nCaloCand]-TrackCand_phi[0])));
	  float calodr2 = sqrt(((CaloTower_eta[nCaloCand]-TrackCand_eta[1])*(CaloTower_eta[nCaloCand]-TrackCand_eta[1])) + 
			       ((CaloTower_phi[nCaloCand]-TrackCand_phi[1])*(CaloTower_phi[nCaloCand]-TrackCand_phi[1])));
	  
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

              // Thresholds tuned on ZeroBias data 
              if(CaloTower_emE[nCaloCand] > 0.6 && CaloTower_ID[nCaloCand] == 1)  
                nExtraCaloTowersE0pt6eb=0;  
              if(CaloTower_emE[nCaloCand] > 2.45 && CaloTower_ID[nCaloCand] == 2) 
                nExtraCaloTowersE2pt45ee=0;  
              if(CaloTower_hadE[nCaloCand] > 1.25 && CaloTower_ID[nCaloCand] == 4) 
                nExtraCaloTowersE1pt25hb=0;   
              if(CaloTower_hadE[nCaloCand] > 1.9 && CaloTower_ID[nCaloCand] == 5)  
                nExtraCaloTowersE1pt9he=0;   
              if(CaloTower_e[nCaloCand] > 4.5 && CaloTower_eta[nCaloCand] > 2.95) 
		nExtraCaloTowersE4pt5hfp=0;  
              if(CaloTower_e[nCaloCand] > 4.0 && CaloTower_eta[nCaloCand] < -2.95)       
		nExtraCaloTowersE4pt0hfm=0;   

	      if(CaloTower_eta[nCaloCand] > 3.0)
		SumHFPlus_e += CaloTower_e[nCaloCand];
	      if(CaloTower_eta[nCaloCand] < -3.0)
		SumHFMinus_e += CaloTower_e[nCaloCand];
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

  // Get the ZDC rechits collection from the event
  edm::Handle<ZDCRecHitCollection> recoZDChits;
  event.getByLabel(recZDCRecHitsLabel, recoZDChits);
  const ZDCRecHitCollection* zdchits = recoZDChits.product();
  ZDCRecHitCollection::const_iterator zdchit;

  for ( zdchit = zdchits->begin(); zdchit != zdchits->end(); ++zdchit )
    {
      HcalZDCDetId id(zdchit->id());
      int Side      = (zdchit->id()).zside();
      int Section   = (zdchit->id()).section();

      ZDChit_section[nZDChitCand] = Section;
      ZDChit_energy[nZDChitCand] = zdchit->energy();
      ZDChit_time[nZDChitCand] = zdchit->time();
      ZDChit_side[nZDChitCand] = Side;

      if((Section == 1) && (Side == 1))
	{
	  ZDCsumEMplus = ZDCsumEMplus + ZDChit_energy[nZDChitCand];
	}
      if((Section == 1) && (Side == -1))
	{
	  ZDCsumEMminus = ZDCsumEMminus + ZDChit_energy[nZDChitCand];
	}
      if((Section == 2) && (Side == 1))
	{
	  ZDCsumHADplus = ZDCsumHADplus + ZDChit_energy[nZDChitCand];
	}
      if((Section == 2) && (Side == -1))
	{
	  ZDCsumHADminus = ZDCsumHADminus + ZDChit_energy[nZDChitCand];
	}
      nZDChitCand++;
    }

  // Now CASTOR rechits
  edm::Handle<CastorRecHitCollection> recoCASTORhits;
  event.getByLabel(recCastorRecHitsLabel, recoCASTORhits);
  const CastorRecHitCollection* castorhits = recoCASTORhits.product();
  CastorRecHitCollection::const_iterator castorhit;
  for ( castorhit = castorhits->begin(); castorhit != castorhits->end(); ++castorhit )
    {
      CASTORsumRecHitsE += castorhit->energy();
    }

  // Now CASTOR towers
  // Get the CASTOR towers collection from the event
  edm::Handle<reco::CastorTowerCollection> recoCastorTowers;
  event.getByLabel(recCastorTowerLabel, recoCastorTowers);
  
  if(recoCastorTowers.isValid())
    {
      const CastorTowerCollection* castortowers = recoCastorTowers.product();
      CastorTowerCollection::const_iterator castortower;
      
      for ( castortower = castortowers->begin(); castortower != castortowers->end(); ++castortower )
	{
	  CastorTower_e[nCastorTowerCand] = castortower->energy();
	  CastorTower_eta[nCastorTowerCand] = castortower->eta();
	  CastorTower_phi[nCastorTowerCand] = castortower->phi();
	  CastorTower_emratio[nCastorTowerCand] = castortower->fem();
	  
	  if(CastorTower_eta[nCastorTowerCand] > 0)
	    {
	      totalecastorfwd+=CastorTower_e[nCastorTowerCand];
	      if(CastorTower_e[nCastorTowerCand] > highestcastortowerfwd)
		highestcastortowerfwd = CastorTower_e[nCastorTowerCand];
	    }
	  if(CastorTower_eta[nCastorTowerCand] < 0)
	    {
	      totalecastorbwd+=CastorTower_e[nCastorTowerCand];
	      if(CastorTower_e[nCastorTowerCand] > highestcastortowerbwd)
		highestcastortowerbwd = CastorTower_e[nCastorTowerCand];
	    }
	  
	  nCastorTowerCand++;
	}
      
      HighestCastorTowerFwd_e = highestcastortowerfwd;
      HighestCastorTowerBwd_e = highestcastortowerbwd;
      SumCastorFwd_e = totalecastorfwd;
      SumCastorBwd_e = totalecastorbwd;
    }
     
     // Check for di-objects
  if((nTrackCand != 2) && (requiretwotracks == true))
    passed = false;
  else
    {
    }

  // "Exclusivity" cuts
  if(passed == true)
    thetree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ExclusiveTrackTrack::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExclusiveTrackTrack::endJob() {
  const edm::ParameterSet &thepset = edm::getProcessParameterSet();
  TList *list = thetree->GetUserInfo();
  list->Add(new TObjString(thepset.dump().c_str()));
  thefile->Write();
  thefile->Close();
}
  
