#ifndef DiffractiveForwardAnalysis_GammaGammaLL8TeV_h
#define DiffractiveForwardAnalysis_GammaGammaLL8TeV_h

// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

// "Includes sandbox" : please delete me or merge me once the tests are done
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
//#include "FWCore/Framework/interface/ESHandle.h"

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Electrons collection
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"  
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"   
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#define MAX_LL    50 // Maximum number of leptons per event
#define MAX_MUONS 25 // Maximum number of muons per event
#define MAX_ELE   25 // Maximum number of electrons per event
#define MAX_PAIRS 25 // Maximum number of leptons pairs per event

#define MASS_MU 0.1057
#define MASS_E  0.000511
#define pi 3.14159265359

typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals; 

//
// class declaration
//

class GammaGammaLL : public edm::EDAnalyzer {
   public:
      explicit GammaGammaLL(const edm::ParameterSet&);
      ~GammaGammaLL();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::Handle<reco::BeamSpot> beamspot_h;
      edm::InputTag beamSpotInputTag_;

      // Electron Id
      std::vector<edm::InputTag> isoValInputTag_; 
      edm::InputTag conversionsInputTag_;
      edm::InputTag rhoIsoInputTag_;
      
      std::vector<std::string> leptonsType_;
      std::string outputFile_;
      edm::InputTag recoVertexLabel_, muonLabel_, eleLabel_;
      bool _fetchMuons, _fetchElectrons;
      
      bool foundPair;
  	  Double_t leptonsDist, minDist; 
      Int_t nMuonCand, nEleCand;
      Int_t nLeptonCand, nPrimVertexCand, nLeptonsInPrimVertex, nUnfilteredPrimVertexCand, nCandidates;
      TLorentzVector* _leptonptmp;
      TLorentzVector l1, l2;
      std::map<Int_t, TLorentzVector> muonsMomenta, electronsMomenta;
      
      // Tree contents
      
      // Muon quantities
      Double_t MuonCand_px[MAX_LL], MuonCand_py[MAX_LL], MuonCand_pz[MAX_LL];
      Double_t MuonCand_p[MAX_LL], MuonCand_pt[MAX_LL];
      Double_t MuonCand_eta[MAX_LL], MuonCand_phi[MAX_LL];
      Double_t MuonCand_vtxx[MAX_LL], MuonCand_vtxy[MAX_LL], MuonCand_vtxz[MAX_LL];
      Double_t MuonCand_charge[MAX_LL];
      Int_t MuonCand_isglobal[MAX_LL], MuonCand_istracker[MAX_LL], MuonCand_isstandalone[MAX_LL];

      // Electron quantities
      Double_t EleCand_px[MAX_LL], EleCand_py[MAX_LL], EleCand_pz[MAX_LL];
      Double_t EleCand_p[MAX_LL], EleCand_e[MAX_LL], EleCand_et[MAX_LL];
      Double_t EleCand_eta[MAX_LL], EleCand_phi[MAX_LL];
      Double_t EleCand_vtxx[MAX_LL], EleCand_vtxy[MAX_LL], EleCand_vtxz[MAX_LL];
      Double_t EleCand_charge[MAX_LL];
		  Double_t EleCandTrack_p[MAX_LL], EleCandTrack_pt[MAX_LL];
		  Double_t EleCandTrack_eta[MAX_LL], EleCandTrack_phi[MAX_LL];
      Double_t EleCandTrack_vtxz[MAX_LL]; 
      Double_t EleCand_deltaPhi[MAX_LL], EleCand_deltaEta[MAX_LL];
      Double_t EleCand_HoverE[MAX_LL];
      Double_t EleCand_trackiso[MAX_LL], EleCand_ecaliso[MAX_LL], EleCand_hcaliso[MAX_LL];
      Double_t EleCand_sigmaIetaIeta[MAX_LL];
      Double_t EleCand_convDist[MAX_LL], EleCand_convDcot[MAX_LL];
      Int_t EleCand_ecalDriven[MAX_LL]; 
      Int_t EleCand_wp80[MAX_LL];
      Int_t EleCand_mediumID[MAX_LL], EleCand_looseID[MAX_LL];
      
      // Pair quantities
			Int_t Pair_candidates[MAX_PAIRS][2];
			Double_t Pair_p[MAX_PAIRS], Pair_pt[MAX_PAIRS], Pair_dpt[MAX_PAIRS];
			Double_t Pair_mass[MAX_PAIRS], Pair_dphi[MAX_PAIRS];
			Double_t Pair_eta[MAX_PAIRS], Pair_phi[MAX_PAIRS], Pair_3Dangle[MAX_PAIRS];

      edm::Handle<reco::VertexCollection> recoVertexColl;
      reco::VertexCollection::const_iterator vertex;
      TString* _leptonType;
      
      // PAT muons
      edm::Handle<edm::View<pat::Muon> > muonColl;
      edm::View<pat::Muon>::const_iterator muon;
      // AOD muons
      /*edm::Handle<reco::MuonCollection> muonColl;
      reco::MuonCollection::const_iterator muon;*/
      
      // PAT electrons
      edm::Handle<edm::View<pat::Electron> > eleColl;
      edm::View<pat::Electron>::const_iterator electron;
      // RECO electrons
      edm::Handle<reco::GsfElectronCollection> eleCollRECO;
      //reco::GsfElectronCollection::const_iterator reco_electron;
      edm::Handle<reco::ConversionCollection> conversions_h;
      edm::Handle<double> rhoIso_h; 
      
};

class PrimaryVertex : public reco::Vertex {
  public:
    explicit PrimaryVertex(std::vector<std::string>&, std::map<Int_t,TLorentzVector>&, std::map<Int_t,TLorentzVector>&);
    ~PrimaryVertex();
    void SetPosition(Double_t, Double_t, Double_t);
    Int_t AddTrack(const reco::TrackRef&, TString&);
    //bool ComputeKinematicQuantities();
    bool Flatten();

    std::vector<Int_t> MatchedMuons, MatchedElectrons;
    
  private:
    //std::vector<std::string> LeptonsType;
    Int_t nMatchedMuons, nMatchedElectrons;
    Int_t nMatchedTracks;
    bool FetchMuons, FetchElectrons;
    std::map<Int_t,TLorentzVector> MuonMomenta;
    std::map<Int_t,TLorentzVector> ElectronMomenta;
    TVector3 Position;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
TFile* file;
TTree* tree;

#endif
