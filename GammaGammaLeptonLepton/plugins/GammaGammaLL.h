#ifndef DiffractiveForwardAnalysis_GammaGammaLL_h
#define DiffractiveForwardAnalysis_GammaGammaLL_h

// system include files
#include <fstream>
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

// HLT information
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// Generator level collection
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// Pileup
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "DiffractiveForwardAnalysis/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Electrons collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Photons collection
#include "DataFormats/PatCandidates/interface/Photon.h"

// Particle flow collection
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// Jets/MET collection
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

// CT-PPS objects
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

// HPS acceptance
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/PrimaryVertexSelector.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/HLTMatcher.h"

// LHC fill information
//#include "DataFormats/Common/interface/ConditionsInEdm.h" // L1 method
//#include "CondFormats/RunInfo/interface/FillInfo.h"
//#include "CondFormats/DataRecord/interface/FillInfoRcd.h" // db method

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>

#define MAX_HLT    10   // Maximum number of HLT to check
#define MAX_LL     50   // Maximum number of leptons per event
#define MAX_MUONS  25   // Maximum number of muons per event
#define MAX_ELE    25   // Maximum number of electrons per event
#define MAX_PHO    50   // Maximum number of photons per event
#define MAX_PAIRS  25   // Maximum number of leptons pairs per event
#define MAX_VTX    1000 // Maximum number of primary vertices per event
#define MAX_ET     10000// Maximum number of extra tracks per event
#define MAX_GENMU  25   // Maximum number of generator level muons per event
#define MAX_GENELE 25   // Maximum number of generator level electrons per event
#define MAX_GENPHO 10   // Maximum number of generator level photons per event
#define MAX_GENPRO 8    // Maximum number of generator level protons per event
#define MAX_JETS   30   // Maximum number of jets per event
#define MAX_LOCALPCAND 10 // Maximum number of reconstructed local tracks in RPs
#define MAX_LOCALPPAIRCAND 5 // Maximum number of reconstructed local tracks pairs in RPs

#define MASS_MU 0.1057
#define MASS_E  0.000511
#define MASS_P  0.938272029
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
      
      virtual void lookAtTriggers(const edm::Event&, const edm::EventSetup&);
      void clearTree();

      // ----------member data ---------------------------

      bool fetchMuons_, fetchElectrons_, fetchProtons_;
      
      unsigned int verb_;

      std::ofstream *logfile_;
      
      // Input tags
      std::string outputFile_;
      std::vector<std::string> leptonsType_;
      std::string hltMenuLabel_;
      std::vector<std::string> triggersList_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<reco::VertexCollection> recoVertexToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
      edm::EDGetTokenT< edm::View<pat::Muon> > muonToken_;
      edm::EDGetTokenT< edm::View<pat::Electron> > eleToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      // rhoIsoLabel_;
      edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT< edm::View<pat::Jet> > jetToken_;
      edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
      edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > totemRPHitToken_;
      std::vector<edm::InputTag> isoValLabel_; 
      bool runOnMC_, printCandidates_;
      double minPtMC_, minEtaMC_;
      double sqrts_;
      unsigned int maxExTrkVtx_;

      // Extended LHC run information
      /*edm::EDGetTokenT<edm::ConditionsInRunBlock> condInRunToken_;
      edm::Handle<edm::ConditionsInRunBlock> condInRunBlock_;*/

      // Beam spot
      edm::Handle<reco::BeamSpot> beamspot_h;

      // Generator level information
      edm::Handle<reco::GenParticleCollection> genPartColl_;
      reco::GenParticleCollection::const_iterator genPart;
      std::string fullAcceptancePath;
      edm::FileInPath *myDataFile;
      // HPS acceptance tables
      TFile* f;
      AcceptanceTableHelper helper420beam1, helper220beam1, helper420a220beam1;
      AcceptanceTableHelper helper420beam2, helper220beam2, helper420a220beam2;

      // Trigger information
      HLTMatcher* hlts_;
      HLTConfigProvider hltConfig_;
      HLTPrescaleProvider hltPrescale_;

      // Pileup information
      edm::LumiReWeighting *lumiWeights_; 
      edm::Handle< std::vector<PileupSummaryInfo> >  pileupInfo;
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      int sum_nvtx, beamXing;
      int npv, npvtrue, npvm1true, npvp1true, npv0true, npv0;
      std::string mcPileupFile_, mcPileupPath_, dataPileupFile_, dataPileupPath_;
      
      TLorentzVector leptonptmp_;
      
      // Isolation
      double rhoIso;
      double iso_ch, iso_em, iso_nh; // Electron isolation quantities
      int vtxind; // Primary vertex index (used in loop over vertices)
      int etind; // Extra tracks on vertex index (used in loop over tracks)

      // Vertices
      edm::Handle<reco::VertexCollection> recoVertexColl_;
      reco::VertexCollection::const_iterator vertex;
      
      // PAT muons
      edm::Handle<edm::View<pat::Muon> > muonColl_;
      edm::View<pat::Muon>::const_iterator muon;
      // AOD muons
      /*edm::Handle<reco::MuonCollection> muonColl_;
      reco::MuonCollection::const_iterator muon;*/
      
      // PAT electrons
      edm::Handle<edm::View<pat::Electron> > eleColl_;
      edm::View<pat::Electron>::const_iterator electron;
      // RECO electrons
      edm::Handle<reco::ConversionCollection> conversions_h;
      edm::Handle<double> rhoIso_h; 
    	
      TLorentzVector pair;
      double dphi;
     
      // PAT photons
      edm::EDGetTokenT< edm::View<pat::Photon> > photonToken_;
      edm::Handle< edm::View<pat::Photon> > photonColl_;
      edm::View<pat::Photon>::const_iterator photon;
 
      // Particle Flow
      edm::EDGetTokenT< edm::View<reco::PFCandidate> > pflowToken_;
      edm::Handle< edm::View<reco::PFCandidate> > pflowColl_;
      edm::View<reco::PFCandidate>::const_iterator pflow; //FIXME will not work for miniAOD
      /*edm::EDGetTokenT< edm::View<pat::PackedCandidate> > pflowToken_;
      edm::Handle<edm::View<pat::PackedCandidate> > pflowColl_;
      edm::View<pat::PackedCandidate>::const_iterator pflow;*/
      
      // Jets/MET
      edm::Handle<edm::View<pat::Jet> > jetColl_;
      edm::View<pat::Jet>::const_iterator jet;
      edm::Handle< edm::View<pat::MET> > MET; 
      edm::View<pat::MET>::const_iterator met; 
      double HEJet_e, HEJet_eta, HEJet_phi;
      double totalJetEnergy;

      ////// Tree contents //////
      
      // Run/event quantities
      int BX, Run, LumiSection, EventNum;
      //int LHCFillNum, LHCBeamMode;
      //double AvgInstDelLumi, BunchInstLumi[3]; 
      
      // HLT quantities
      int nHLT;
      int HLT_Accept[MAX_HLT], HLT_Prescl[MAX_HLT];
      //char* HLT_Name[MAX_HLT];
      std::vector<std::string>* HLT_Name;
      //std::map<int,std::string>* HLT_Name;
      /*int nHLTLeptonCand[MAX_HLT];
      double HLTLeptonCand_pt[2][MAX_HLT];
      double HLTLeptonCand_eta[2][MAX_HLT];
      double HLTLeptonCand_phi[2][MAX_HLT];
      int HLTLeptonCand_charge[2][MAX_HLT];
      int HLT_LeadingLepton[MAX_HLT], HLT_TrailingLepton[MAX_HLT];
      int HLT_LeadingLepton_Prescl[MAX_HLT], HLT_TrailingLepton_Prescl[MAX_HLT];*/
      
      // Generator level quantities
      int nGenMuonCand, nGenMuonCandOutOfAccept;
      double GenMuonCand_px[MAX_GENMU], GenMuonCand_py[MAX_GENMU], GenMuonCand_pz[MAX_GENMU];
      double GenMuonCand_p[MAX_GENMU], GenMuonCand_pt[MAX_GENMU];
      double GenMuonCand_eta[MAX_GENMU], GenMuonCand_phi[MAX_GENMU];
      int nGenEleCand, nGenEleCandOutOfAccept;
      double GenEleCand_px[MAX_GENELE], GenEleCand_py[MAX_GENELE], GenEleCand_pz[MAX_GENELE];
      double GenEleCand_p[MAX_GENELE], GenEleCand_pt[MAX_GENELE];
      double GenEleCand_eta[MAX_GENELE], GenEleCand_phi[MAX_GENELE];
      double GenPair_p, GenPair_pt, GenPair_mass;
      double GenPair_phi, GenPair_eta;
      double GenPair_dphi, GenPair_dpt, GenPair_3Dangle;
      int nGenPhotCand, nGenPhotCandOutOfAccept;
      double GenPhotCand_p[MAX_GENPHO], GenPhotCand_e[MAX_GENPHO];
      double GenPhotCand_pt[MAX_GENPHO], GenPhotCand_eta[MAX_GENPHO], GenPhotCand_phi[MAX_GENPHO];
      int nGenProtCand;
      double GenProtCand_p[MAX_GENPRO], GenProtCand_px[MAX_GENPRO], GenProtCand_py[MAX_GENPRO], GenProtCand_pz[MAX_GENPRO];
      double GenProtCand_pt[MAX_GENPRO], GenProtCand_eta[MAX_GENPRO], GenProtCand_phi[MAX_GENPRO];
      int GenProtCand_status[MAX_GENPRO];

      // HPS quantities
      double xi, t;
      double HPS_acc420b1, HPS_acc220b1, HPS_acc420and220b1, HPS_acc420or220b1; // beam 1 (clockwise)  
      double HPS_acc420b2, HPS_acc220b2, HPS_acc420and220b2, HPS_acc420or220b2; // beam 2 (anti-clockwise)  

      int nLeptonCand, nLeptonsInPrimVertex, nCandidates, nCandidatesInEvent;

      // Pileup reweighting quantities
      double nTruePUafterPUWeight;
      double nTruePUafterPUWeightBXM1, nTruePUafterPUWeightBXP1, nTruePUafterPUWeightBX0;
      double PUWeightTrue;
      int nTruePUforPUWeight;
      int nTruePUforPUWeightBXM1, nTruePUforPUWeightBXP1, nTruePUforPUWeightBX0;
      double Weight;

      // Muon quantities
      int nMuonCand;
      double MuonCand_px[MAX_LL], MuonCand_py[MAX_LL], MuonCand_pz[MAX_LL];
      double MuonCand_p[MAX_LL], MuonCand_pt[MAX_LL];
      double MuonCand_eta[MAX_LL], MuonCand_phi[MAX_LL];
      double MuonCand_vtxx[MAX_LL], MuonCand_vtxy[MAX_LL], MuonCand_vtxz[MAX_LL];
      int MuonCand_charge[MAX_LL];
      double MuonCand_dxy[MAX_LL], MuonCand_dz[MAX_LL];
      int MuonCand_nstatseg[MAX_LL], MuonCand_npxlhits[MAX_LL], MuonCand_ntrklayers[MAX_LL];
      double MuonCand_[MAX_LL];
      int MuonCandTrack_nmuchits[MAX_LL];
      double MuonCandTrack_chisq[MAX_LL];
      int MuonCand_isglobal[MAX_LL], MuonCand_istracker[MAX_LL], MuonCand_isstandalone[MAX_LL], MuonCand_ispfmuon[MAX_LL];
      int MuonCand_istight[MAX_LL];

      // Electron quantities
      int nEleCand;
      double EleCand_px[MAX_LL], EleCand_py[MAX_LL], EleCand_pz[MAX_LL];
      double EleCand_p[MAX_LL], EleCand_e[MAX_LL], EleCand_et[MAX_LL];
      double EleCand_eta[MAX_LL], EleCand_phi[MAX_LL];
      double EleCand_vtxx[MAX_LL], EleCand_vtxy[MAX_LL], EleCand_vtxz[MAX_LL];
      int EleCand_charge[MAX_LL];
      double EleCandTrack_p[MAX_LL], EleCandTrack_pt[MAX_LL];
      double EleCandTrack_eta[MAX_LL], EleCandTrack_phi[MAX_LL];
      double EleCandTrack_vtxz[MAX_LL]; 
      double EleCand_deltaPhi[MAX_LL], EleCand_deltaEta[MAX_LL];
      double EleCand_HoverE[MAX_LL];
      double EleCand_trackiso[MAX_LL], EleCand_ecaliso[MAX_LL], EleCand_hcaliso[MAX_LL];
      double EleCand_sigmaIetaIeta[MAX_LL];
      double EleCand_convDist[MAX_LL], EleCand_convDcot[MAX_LL];
      int EleCand_ecalDriven[MAX_LL]; 
      int EleCand_tightID[MAX_LL], EleCand_mediumID[MAX_LL], EleCand_looseID[MAX_LL];
      
      // Photon quantities
      int nPhotonCand;
      double PhotonCand_px[MAX_PHO], PhotonCand_py[MAX_PHO], PhotonCand_pz[MAX_PHO];
      double PhotonCand_p[MAX_PHO], PhotonCand_pt[MAX_PHO];
      double PhotonCand_eta[MAX_PHO], PhotonCand_phi[MAX_PHO], PhotonCand_r9[MAX_PHO];
      double PhotonCand_drtrue[MAX_PHO], PhotonCand_detatrue[MAX_PHO], PhotonCand_dphitrue[MAX_PHO];
      
      // Pair quantities
      int Pair_candidates[MAX_PAIRS][2], Pair_lepton1[MAX_PAIRS], Pair_lepton2[MAX_PAIRS];
      double Pair_mindist[MAX_PAIRS];
      double Pair_p[MAX_PAIRS], Pair_pt[MAX_PAIRS], Pair_dpt[MAX_PAIRS];
      double Pair_mass[MAX_PAIRS], Pair_dphi[MAX_PAIRS];
      double Pair_eta[MAX_PAIRS], Pair_phi[MAX_PAIRS], Pair_3Dangle[MAX_PAIRS];
      double PairGamma_mass[MAX_PAIRS][MAX_PHO];
      // Extra tracks
      int Pair_extratracks1mm[MAX_PAIRS], Pair_extratracks2mm[MAX_PAIRS];
      int Pair_extratracks3mm[MAX_PAIRS], Pair_extratracks4mm[MAX_PAIRS];
      int Pair_extratracks5mm[MAX_PAIRS], Pair_extratracks1cm[MAX_PAIRS];
      int Pair_extratracks2cm[MAX_PAIRS], Pair_extratracks3cm[MAX_PAIRS];
      int Pair_extratracks4cm[MAX_PAIRS], Pair_extratracks5cm[MAX_PAIRS];
      int Pair_extratracks10cm[MAX_PAIRS];
      
      // Vertex quantities
      int nPrimVertexCand;
      int PrimVertexCand_id[MAX_VTX], PrimVertexCand_hasdil[MAX_VTX];
      double PrimVertexCand_x[MAX_VTX], PrimVertexCand_y[MAX_VTX], PrimVertexCand_z[MAX_VTX];
      int PrimVertexCand_tracks[MAX_VTX], PrimVertexCand_matchedtracks[MAX_VTX], PrimVertexCand_unmatchedtracks[MAX_VTX];
      double PrimVertexCand_chi2[MAX_VTX];
      int PrimVertexCand_ndof[MAX_VTX];
      int nFilteredPrimVertexCand;
      
      // Extra tracks on vertex quantities
      int nExtraTracks;
      int ExtraTrack_purity[MAX_ET], ExtraTrack_nhits[MAX_ET];
      int ExtraTrack_charge[MAX_ET], ExtraTrack_ndof[MAX_ET];
      int ExtraTrack_vtxId[MAX_ET];
      double ExtraTrack_p[MAX_ET], ExtraTrack_pt[MAX_ET];
      double ExtraTrack_px[MAX_ET], ExtraTrack_py[MAX_ET], ExtraTrack_pz[MAX_ET];
      double ExtraTrack_eta[MAX_ET], ExtraTrack_phi[MAX_ET];
      double ExtraTrack_chi2[MAX_ET];
      double ExtraTrack_vtxdxyz[MAX_ET];
      double ExtraTrack_vtxT[MAX_ET], ExtraTrack_vtxZ[MAX_ET];
      double ExtraTrack_x[MAX_ET], ExtraTrack_y[MAX_ET], ExtraTrack_z[MAX_ET];
      double ClosestExtraTrack_vtxdxyz[MAX_VTX], ClosestHighPurityExtraTrack_vtxdxyz[MAX_VTX];
      int ClosestExtraTrack_id[MAX_VTX], ClosestHighPurityExtraTrack_id[MAX_VTX];
      int nQualityExtraTrack;

      // Jets/MET quantities
      int nJetCand;
      double JetCand_px[MAX_JETS], JetCand_py[MAX_JETS], JetCand_pz[MAX_JETS];
      double JetCand_e[MAX_JETS], JetCand_eta[MAX_JETS], JetCand_phi[MAX_JETS];
      double HighestJet_e, HighestJet_eta, HighestJet_phi;
      double SumJet_e;
      double Etmiss, Etmiss_phi, Etmiss_x, Etmiss_y, Etmiss_z, Etmiss_significance;

      // CTPPS quantities
      int nLocalProtCand;
      double LocalProtCand_x[MAX_LOCALPCAND], LocalProtCand_y[MAX_LOCALPCAND], LocalProtCand_z[MAX_LOCALPCAND];
      double LocalProtCand_xSigma[MAX_LOCALPCAND], LocalProtCand_ySigma[MAX_LOCALPCAND];
      double LocalProtCand_xi[MAX_LOCALPCAND];
      double LocalProtCand_Tx[MAX_LOCALPCAND], LocalProtCand_Ty[MAX_LOCALPCAND];
      double LocalProtCand_TxSigma[MAX_LOCALPCAND], LocalProtCand_TySigma[MAX_LOCALPCAND];
      int LocalProtCand_arm[MAX_LOCALPCAND], LocalProtCand_side[MAX_LOCALPCAND];

      int nLocalProtPairCand;
      double LocalProtPairCand_mass[MAX_LOCALPPAIRCAND], LocalProtPairCand_pt[MAX_LOCALPPAIRCAND], LocalProtPairCand_y[MAX_LOCALPPAIRCAND];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
TFile* file_;
TTree* tree_;

#endif
