#ifndef DiffractiveForwardAnalysis_GammaGammaLeptonLepton_AnalysisEvent_h
#define DiffractiveForwardAnalysis_GammaGammaLeptonLepton_AnalysisEvent_h

#include <vector>
#include <string>

#include "TTree.h"

namespace ggll
{
  enum TreeType {
    invalidTree = -1,
    DiMuon,
    DiElectron,
    ElectronMuon
  };

  class AnalysisEvent
  {
    public:
      AnalysisEvent() :
        HLT_Name( new std::vector<std::string>() ) {
        clear();
      }
      ~AnalysisEvent() {
        if ( HLT_Name ) delete HLT_Name;
      }

      ////// Tree contents //////

      // Run/event quantities
      unsigned int BX, Run, LumiSection, EventNum;
      //int LHCFillNum, LHCBeamMode;
      //double AvgInstDelLumi, BunchInstLumi[3];

      // HLT quantities
      std::vector<int> HLT_Accept, HLT_Prescl;
      std::vector<std::string>* HLT_Name;
      /*std::vector<int> nHLTLeptonCand;
      std::vector<double> HLTLeptonCand_pt[2];
      std::vector<double> HLTLeptonCand_eta[2];
      std::vector<double> HLTLeptonCand_phi[2];
      std::vector<int> HLTLeptonCand_charge[2];
      std::vector<int> HLT_LeadingLepton, HLT_TrailingLepton;
      std::vector<int> HLT_LeadingLepton_Prescl, HLT_TrailingLepton_Prescl;*/

      // Generator level quantities
      unsigned int nGenMuonCandOutOfAccept;
      std::vector<double> GenMuonCand_pt, GenMuonCand_eta, GenMuonCand_phi, GenMuonCand_e;
      unsigned int nGenEleCandOutOfAccept;
      std::vector<double> GenEleCand_pt, GenEleCand_eta, GenEleCand_phi, GenEleCand_e;
      double GenPair_pt, GenPair_eta, GenPair_phi, GenPair_mass;
      double GenPair_dphi, GenPair_dpt, GenPair_3Dangle;
      unsigned int nGenPhotCandOutOfAccept;
      std::vector<double> GenPhotCand_pt, GenPhotCand_eta, GenPhotCand_phi, GenPhotCand_e;
      std::vector<double> GenProtCand_pt, GenProtCand_eta, GenProtCand_phi, GenProtCand_e;
      std::vector<int> GenProtCand_status;

      // Pileup reweighting quantities
      double PUWeightTrue, Weight;

      // Muon quantities
      std::vector<double> MuonCand_pt, MuonCand_eta, MuonCand_phi, MuonCand_e;
      std::vector<double> MuonCand_innerTrackPt, MuonCand_innerTrackEta, MuonCand_innerTrackPhi;
      std::vector<double> MuonCand_innerTrackVtxz;
      std::vector<double> MuonCand_vtxx, MuonCand_vtxy, MuonCand_vtxz;
      std::vector<int> MuonCand_charge;
      std::vector<double> MuonCand_dxy;
      std::vector<int> MuonCand_nstatseg, MuonCand_npxlhits, MuonCand_ntrklayers;
      std::vector<int> MuonCandTrack_nmuchits;
      std::vector<double> MuonCandTrack_chisq;
      std::vector<int> MuonCand_isglobal, MuonCand_istracker, MuonCand_isstandalone, MuonCand_ispfmuon;
      std::vector<int> MuonCand_istight;

      // Electron quantities
      unsigned int nEleCand;
      std::vector<double> EleCand_et, EleCand_eta, EleCand_phi, EleCand_e;
      std::vector<double> EleCand_vtxx, EleCand_vtxy, EleCand_vtxz;
      std::vector<double> EleCand_innerTrackPt, EleCand_innerTrackEta, EleCand_innerTrackPhi;
      std::vector<double> EleCand_innerTrackVtxz;
      std::vector<int> EleCand_charge;
      std::vector<double> EleCand_deltaPhi, EleCand_deltaEta;
      std::vector<double> EleCand_HoverE;
      std::vector<double> EleCand_trackiso, EleCand_ecaliso, EleCand_hcaliso;
      std::vector<double> EleCand_sigmaIetaIeta;
      std::vector<double> EleCand_convDist, EleCand_convDcot;
      std::vector<int> EleCand_ecalDriven;
      std::vector<int> EleCand_tightID, EleCand_mediumID;

      // Photon quantities
      unsigned int nPhotonCand;
      std::vector<double> PhotonCand_pt, PhotonCand_eta, PhotonCand_phi, PhotonCand_e;
      std::vector<double> PhotonCand_r9;
      std::vector<double> PhotonCand_drtrue, PhotonCand_detatrue, PhotonCand_dphitrue;
      std::vector<int> PhotonCand_tightID, PhotonCand_mediumID;

      // Pair quantities
      unsigned int nPair;
      std::vector<int> Pair_lepton1, Pair_lepton2;
      std::vector<double> Pair_pt, Pair_eta, Pair_phi, Pair_mass;
      std::vector<double> Pair_dpt, Pair_dphi, Pair_3Dangle;
      //std::vector<double> Pair_mindist;

      unsigned int nPairGamma;
      std::vector<int> PairGamma_pair, PairGamma_pho;
      std::vector<double> PairGamma_mass;

      // Extra tracks
      std::vector<unsigned int> Pair_extratracks0p5mm;
      std::vector<unsigned int> Pair_extratracks1mm, Pair_extratracks2mm;
      std::vector<unsigned int> Pair_extratracks3mm, Pair_extratracks4mm;
      std::vector<unsigned int> Pair_extratracks5mm, Pair_extratracks1cm;
      std::vector<unsigned int> Pair_extratracks2cm, Pair_extratracks3cm;
      std::vector<unsigned int> Pair_extratracks4cm, Pair_extratracks5cm;
      std::vector<unsigned int> Pair_extratracks10cm;

      // Vertex quantities
      unsigned int nPrimVertexCand;
      std::vector<int> PrimVertexCand_id, PrimVertexCand_hasdil;
      std::vector<double> PrimVertexCand_x, PrimVertexCand_y, PrimVertexCand_z;
      std::vector<unsigned int> PrimVertexCand_tracks;
      std::vector<double> PrimVertexCand_chi2;
      std::vector<unsigned int> PrimVertexCand_ndof;
      std::vector<double> KalmanVertexCand_x, KalmanVertexCand_y, KalmanVertexCand_z;
      unsigned int nFilteredPrimVertexCand;

      // Extra tracks on vertex quantities
      unsigned int nExtraTracks;
      std::vector<int> ExtraTrack_pair;
      std::vector<int> ExtraTrack_purity;
      std::vector<unsigned int> ExtraTrack_nhits;
      std::vector<int> ExtraTrack_charge;
      std::vector<unsigned int> ExtraTrack_ndof;
      std::vector<double> ExtraTrack_px, ExtraTrack_py, ExtraTrack_pz;
      std::vector<double> ExtraTrack_chi2;
      std::vector<double> ExtraTrack_vtxdxyz;
      std::vector<double> ExtraTrack_vtxT, ExtraTrack_vtxZ;
      std::vector<double> ExtraTrack_x, ExtraTrack_y, ExtraTrack_z;

      std::vector<double> ClosestExtraTrack_vtxdxyz, ClosestHighPurityExtraTrack_vtxdxyz;
      std::vector<int> ClosestExtraTrack_id, ClosestHighPurityExtraTrack_id;
      unsigned int nQualityExtraTrack;

      // Jets/MET quantities
      std::vector<double> JetCand_pt, JetCand_eta, JetCand_phi, JetCand_e;
      double HighestJet_pt, HighestJet_eta, HighestJet_phi, HighestJet_e;
      double SumJet_e;
      double Etmiss, Etmiss_phi, Etmiss_significance;

      // CTPPS quantities
      unsigned int nLocalProtCand;
      std::vector<double> LocalProtCand_x, LocalProtCand_y, LocalProtCand_z;
      std::vector<double> LocalProtCand_xSigma, LocalProtCand_ySigma;
      std::vector<double> LocalProtCand_Tx, LocalProtCand_Ty;
      std::vector<double> LocalProtCand_TxSigma, LocalProtCand_TySigma;
      std::vector<int> LocalProtCand_arm, LocalProtCand_pot;

      void clear() {
        // event-level branches
        BX = Run = LumiSection = EventNum = 0;

        // high-level trigger
        HLT_Accept.clear(); HLT_Prescl.clear();
        HLT_Name->clear();

        // gen-level information
        nGenMuonCandOutOfAccept = 0;
        GenMuonCand_pt.clear(); GenMuonCand_eta.clear(); GenMuonCand_phi.clear(); GenMuonCand_e.clear();
        nGenEleCandOutOfAccept = 0;
        GenEleCand_pt.clear(); GenEleCand_eta.clear(); GenEleCand_phi.clear(); GenEleCand_e.clear();
        GenPair_pt = GenPair_eta = GenPair_phi = GenPair_mass = -999.;
        GenPair_dphi = GenPair_dpt = GenPair_3Dangle = -999.;
        nGenPhotCandOutOfAccept = 0;
        GenPhotCand_pt.clear(); GenPhotCand_eta.clear(); GenPhotCand_phi.clear(); GenPhotCand_e.clear();
        GenProtCand_pt.clear(); GenProtCand_eta.clear(); GenProtCand_phi.clear(); GenProtCand_e.clear();
        GenProtCand_status.clear();

        PUWeightTrue = Weight = 0.;

        //LHCFillNum = LHCBeamMode = -1;

        // single lepton candidates
        MuonCand_pt.clear(); MuonCand_eta.clear(); MuonCand_phi.clear(); MuonCand_e.clear();
        MuonCand_innerTrackPt.clear(); MuonCand_innerTrackEta.clear(); MuonCand_innerTrackPhi.clear();
        MuonCand_innerTrackVtxz.clear();
        MuonCand_vtxx.clear(); MuonCand_vtxy.clear(); MuonCand_vtxz.clear();
        MuonCand_charge.clear();
        MuonCand_dxy.clear();
        MuonCand_nstatseg.clear(); MuonCand_npxlhits.clear(); MuonCand_ntrklayers.clear();
        MuonCandTrack_nmuchits.clear();
        MuonCandTrack_chisq.clear();
        MuonCand_isglobal.clear(); MuonCand_istracker.clear(); MuonCand_isstandalone.clear(); MuonCand_ispfmuon.clear();
        MuonCand_istight.clear();
        EleCand_et.clear(); EleCand_eta.clear(); EleCand_phi.clear(); EleCand_e.clear();
        EleCand_vtxx.clear(); EleCand_vtxy.clear(); EleCand_vtxz.clear();
        EleCand_innerTrackPt.clear(); EleCand_innerTrackEta.clear(); EleCand_innerTrackPhi.clear();
        EleCand_innerTrackVtxz.clear();
        EleCand_charge.clear();
        EleCand_deltaPhi.clear(); EleCand_deltaEta.clear();
        EleCand_HoverE.clear();
        EleCand_trackiso.clear(); EleCand_ecaliso.clear(); EleCand_hcaliso.clear();
        EleCand_sigmaIetaIeta.clear();
        EleCand_convDist.clear(); EleCand_convDcot.clear();
        EleCand_ecalDriven.clear();
        EleCand_tightID.clear(); EleCand_mediumID.clear();

        // single photon candidates
        PhotonCand_pt.clear(); PhotonCand_eta.clear(); PhotonCand_phi.clear(); PhotonCand_e.clear();
        PhotonCand_r9.clear();
        PhotonCand_drtrue.clear(); PhotonCand_detatrue.clear(); PhotonCand_dphitrue.clear();
        PhotonCand_tightID.clear(); PhotonCand_mediumID.clear();

        // dilepton pair candidates
        Pair_lepton1.clear(); Pair_lepton2.clear();
        Pair_pt.clear(); Pair_mass.clear(); Pair_phi.clear(); Pair_eta.clear();
        Pair_dpt.clear(); Pair_dphi.clear(); Pair_3Dangle.clear();
        //Pair_mindist.clear();

        // dilepton pair + associated photon candidates
        PairGamma_pair.clear();
        PairGamma_pho.clear();
        PairGamma_mass.clear();

        // extra tracks associated to the central system vertex
        Pair_extratracks0p5mm.clear();
        Pair_extratracks1mm.clear(); Pair_extratracks2mm.clear(); Pair_extratracks3mm.clear();
        Pair_extratracks4mm.clear(); Pair_extratracks5mm.clear(); Pair_extratracks1cm.clear();
        Pair_extratracks2cm.clear(); Pair_extratracks3cm.clear(); Pair_extratracks4cm.clear();
        Pair_extratracks5cm.clear(); Pair_extratracks10cm.clear();

        // offline primary vertices
        PrimVertexCand_id.clear(); PrimVertexCand_hasdil.clear();
        PrimVertexCand_x.clear(); PrimVertexCand_y.clear(); PrimVertexCand_z.clear();
        PrimVertexCand_tracks.clear();
        PrimVertexCand_chi2.clear();
        PrimVertexCand_ndof.clear();
        KalmanVertexCand_x.clear(); KalmanVertexCand_y.clear(); KalmanVertexCand_z.clear();
        nFilteredPrimVertexCand = 0;

        // extra tracks associated to the central system
        ExtraTrack_pair.clear();
        ExtraTrack_purity.clear(); ExtraTrack_nhits.clear();
        ExtraTrack_charge.clear();
        ExtraTrack_ndof.clear();
        ExtraTrack_px.clear(); ExtraTrack_py.clear(); ExtraTrack_pz.clear();
        ExtraTrack_chi2.clear(); ExtraTrack_vtxdxyz.clear();
        ExtraTrack_vtxT.clear(); ExtraTrack_vtxZ.clear();
        ExtraTrack_x.clear(); ExtraTrack_y.clear(); ExtraTrack_z.clear();

        nQualityExtraTrack = 0;

        ClosestExtraTrack_vtxdxyz.clear(); ClosestHighPurityExtraTrack_vtxdxyz.clear();
        ClosestExtraTrack_id.clear(); ClosestHighPurityExtraTrack_id.clear();

        // jets collection
        JetCand_pt.clear(); JetCand_eta.clear(); JetCand_phi.clear(); JetCand_e.clear();
        HighestJet_pt = HighestJet_eta = HighestJet_phi = HighestJet_e = -999.;
        SumJet_e = 0.;

        // missing ET
        Etmiss = Etmiss_phi = Etmiss_significance = -999.;

        // CTPPS strips leaves
        LocalProtCand_x.clear(); LocalProtCand_y.clear(); LocalProtCand_z.clear();
        LocalProtCand_xSigma.clear(); LocalProtCand_ySigma.clear();
        LocalProtCand_Tx.clear(); LocalProtCand_Ty.clear();
        LocalProtCand_TxSigma.clear(); LocalProtCand_TySigma.clear();
        LocalProtCand_arm.clear(); LocalProtCand_pot.clear();
      }
      void attach( TTree* tree, TreeType tt, bool mc ) {
        if ( !tree ) return;

        tree->Branch( "Run", &Run );
        tree->Branch( "LumiSection", &LumiSection );
        tree->Branch( "BX", &BX );
        tree->Branch( "EventNum", &EventNum );

        tree->Branch( "HLT_Accept", &HLT_Accept );
        tree->Branch( "HLT_Prescl", &HLT_Prescl );
        tree->Branch( "HLT_Name", &HLT_Name);

        if ( tt == ElectronMuon || tt == DiMuon ) {
          tree->Branch( "MuonCand_pt", &MuonCand_pt );
          tree->Branch( "MuonCand_eta", &MuonCand_eta );
          tree->Branch( "MuonCand_phi", &MuonCand_phi );
          tree->Branch( "MuonCand_e", &MuonCand_e );
          tree->Branch( "MuonCand_charge", &MuonCand_charge );
          tree->Branch( "MuonCand_vtxx", &MuonCand_vtxx );
          tree->Branch( "MuonCand_vtxy", &MuonCand_vtxy );
          tree->Branch( "MuonCand_vtxz", &MuonCand_vtxz );
          tree->Branch( "MuonCand_dxy", &MuonCand_dxy );
          tree->Branch( "MuonCand_nstatseg", &MuonCand_nstatseg );
          tree->Branch( "MuonCand_ntrklayers", &MuonCand_ntrklayers );
          tree->Branch( "MuonCand_npxlhits", &MuonCand_npxlhits );
          tree->Branch( "MuonCand_isglobal", &MuonCand_isglobal );
          tree->Branch( "MuonCand_istracker", &MuonCand_istracker );
          tree->Branch( "MuonCand_isstandalone", &MuonCand_isstandalone );
          tree->Branch( "MuonCand_ispfmuon", &MuonCand_ispfmuon );
          tree->Branch( "MuonCand_istight", &MuonCand_istight );
          tree->Branch( "MuonCandTrack_nmuchits", &MuonCandTrack_nmuchits );
          tree->Branch( "MuonCandTrack_chisq", &MuonCandTrack_chisq );
          tree->Branch( "MuonCand_innerTrackPt", &MuonCand_innerTrackPt );
          tree->Branch( "MuonCand_innerTrackEta", &MuonCand_innerTrackEta );
          tree->Branch( "MuonCand_innerTrackPhi", &MuonCand_innerTrackPhi );
          tree->Branch( "MuonCand_innerTrackVtxz", &MuonCand_innerTrackVtxz );
          if ( mc ) {
            tree->Branch( "nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept );
            tree->Branch( "GenMuonCand_pt", &GenMuonCand_pt );
            tree->Branch( "GenMuonCand_eta", &GenMuonCand_eta );
            tree->Branch( "GenMuonCand_phi", &GenMuonCand_phi );
            tree->Branch( "GenMuonCand_e", &GenMuonCand_e );
          }
        }

        if ( tt == ElectronMuon || tt == DiElectron ) {
          tree->Branch( "EleCand_et", &EleCand_et );
          tree->Branch( "EleCand_eta", &EleCand_eta );
          tree->Branch( "EleCand_phi", &EleCand_phi );
          tree->Branch( "EleCand_e", &EleCand_e );
          tree->Branch( "EleCand_charge", &EleCand_charge );
          tree->Branch( "EleCand_vtxx", &EleCand_vtxx );
          tree->Branch( "EleCand_vtxy", &EleCand_vtxy );
          tree->Branch( "EleCand_vtxz", &EleCand_vtxz );
          tree->Branch( "EleCand_deltaPhi", &EleCand_deltaPhi );
          tree->Branch( "EleCand_deltaEta", &EleCand_deltaEta );
          tree->Branch( "EleCand_HoverE", &EleCand_HoverE );
          tree->Branch( "EleCand_trackiso", &EleCand_trackiso );
          tree->Branch( "EleCand_ecaliso", &EleCand_ecaliso );
          tree->Branch( "EleCand_hcaliso", &EleCand_hcaliso );
          tree->Branch( "EleCand_sigmaIetaIeta", &EleCand_sigmaIetaIeta );
          tree->Branch( "EleCand_convDist", &EleCand_convDist );
          tree->Branch( "EleCand_convDcot", &EleCand_convDcot );
          tree->Branch( "EleCand_ecalDriven", &EleCand_ecalDriven );
          tree->Branch( "EleCand_mediumID", &EleCand_mediumID );
          tree->Branch( "EleCand_tightID", &EleCand_tightID );
          tree->Branch( "EleCand_innerTrackPt", &EleCand_innerTrackPt );
          tree->Branch( "EleCand_innerTrackEta", &EleCand_innerTrackEta );
          tree->Branch( "EleCand_innerTrackPhi", &EleCand_innerTrackPhi );
          tree->Branch( "EleCand_innerTrackVtxz", &EleCand_innerTrackVtxz );
          if ( mc ) {
            tree->Branch( "nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept );
            tree->Branch( "GenEleCand_pt", &GenEleCand_pt );
            tree->Branch( "GenEleCand_eta", &GenEleCand_eta );
            tree->Branch( "GenEleCand_phi", &GenEleCand_phi );
            tree->Branch( "GenEleCand_e", &GenEleCand_e );
          }
        }
        tree->Branch( "PhotonCand_pt", &PhotonCand_pt );
        tree->Branch( "PhotonCand_eta", &PhotonCand_eta );
        tree->Branch( "PhotonCand_phi", &PhotonCand_phi );
        tree->Branch( "PhotonCand_e", &PhotonCand_e );
        tree->Branch( "PhotonCand_r9", &PhotonCand_r9 );
        tree->Branch( "PhotonCand_drtrue", &PhotonCand_drtrue );
        tree->Branch( "PhotonCand_detatrue", &PhotonCand_detatrue );
        tree->Branch( "PhotonCand_dphitrue", &PhotonCand_dphitrue );
        tree->Branch( "PhotonCand_mediumID", &PhotonCand_mediumID );
        tree->Branch( "PhotonCand_tightID", &PhotonCand_tightID );
        if ( mc ) {
          tree->Branch( "nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept );
          tree->Branch( "GenPhotCand_pt", &GenPhotCand_pt );
          tree->Branch( "GenPhotCand_eta", &GenPhotCand_eta );
          tree->Branch( "GenPhotCand_phi", &GenPhotCand_phi );
          tree->Branch( "GenPhotCand_e", &GenPhotCand_e );

          tree->Branch( "GenProtCand_pt", &GenProtCand_pt );
          tree->Branch( "GenProtCand_eta", &GenProtCand_eta );
          tree->Branch( "GenProtCand_phi", &GenProtCand_phi );
          tree->Branch( "GenProtCand_e", &GenProtCand_e );
          tree->Branch( "GenProtCand_status", &GenProtCand_status );
        }

        // Primary vertices' information
        tree->Branch( "nFilteredPrimVertexCand", &nFilteredPrimVertexCand );
        tree->Branch( "PrimVertexCand_id", &PrimVertexCand_id );
        tree->Branch( "PrimVertexCand_x", &PrimVertexCand_x );
        tree->Branch( "PrimVertexCand_y", &PrimVertexCand_y );
        tree->Branch( "PrimVertexCand_z", &PrimVertexCand_z );
        tree->Branch( "PrimVertexCand_chi2", &PrimVertexCand_chi2 );
        tree->Branch( "PrimVertexCand_ndof", &PrimVertexCand_ndof );
        tree->Branch( "PrimVertexCand_tracks", &PrimVertexCand_tracks );

        // Lepton pairs' information
        tree->Branch( "Pair_lepton1", &Pair_lepton1 );
        tree->Branch( "Pair_lepton2", &Pair_lepton2 );
        //tree->Branch( "Pair_mindist", &Pair_mindist );
        tree->Branch( "Pair_mass", &Pair_mass );
        tree->Branch( "Pair_pt", &Pair_pt );
        tree->Branch( "Pair_eta", &Pair_eta );
        tree->Branch( "Pair_phi", &Pair_phi );
        tree->Branch( "Pair_dpt", &Pair_dpt );
        tree->Branch( "Pair_dphi", &Pair_dphi );
        tree->Branch( "Pair_3Dangle", &Pair_3Dangle );
        tree->Branch( "Pair_extratracks0p5mm", &Pair_extratracks0p5mm );
        tree->Branch( "Pair_extratracks1mm", &Pair_extratracks1mm );
        tree->Branch( "Pair_extratracks2mm", &Pair_extratracks2mm );
        tree->Branch( "Pair_extratracks3mm", &Pair_extratracks3mm );
        tree->Branch( "Pair_extratracks4mm", &Pair_extratracks4mm );
        tree->Branch( "Pair_extratracks5mm", &Pair_extratracks5mm );
        tree->Branch( "Pair_extratracks1cm", &Pair_extratracks1cm );
        tree->Branch( "Pair_extratracks2cm", &Pair_extratracks2cm );
        tree->Branch( "Pair_extratracks3cm", &Pair_extratracks3cm );
        tree->Branch( "Pair_extratracks4cm", &Pair_extratracks4cm );
        tree->Branch( "Pair_extratracks5cm", &Pair_extratracks5cm );
        tree->Branch( "Pair_extratracks10cm", &Pair_extratracks10cm );
        // Kalman dilepton vertex information
        tree->Branch( "KalmanVertexCand_x", &KalmanVertexCand_x );
        tree->Branch( "KalmanVertexCand_y", &KalmanVertexCand_y );
        tree->Branch( "KalmanVertexCand_z", &KalmanVertexCand_z );

        tree->Branch( "PairGamma_pair", &PairGamma_pair );
        tree->Branch( "PairGamma_pho", &PairGamma_pho );
        tree->Branch( "PairGamma_mass", &PairGamma_mass );
        if ( mc ) {
          tree->Branch( "GenPair_mass", &GenPair_mass );
          tree->Branch( "GenPair_pt", &GenPair_pt );
          tree->Branch( "GenPair_eta", &GenPair_eta );
          tree->Branch( "GenPair_phi", &GenPair_phi );
          tree->Branch( "GenPair_dpt", &GenPair_dpt );
          tree->Branch( "GenPair_dphi", &GenPair_dphi );
          tree->Branch( "GenPair_3Dangle", &GenPair_3Dangle );
        }

        if ( !mc ) {
          tree->Branch( "LocalProtCand_x", &LocalProtCand_x );
          tree->Branch( "LocalProtCand_y", &LocalProtCand_y );
          tree->Branch( "LocalProtCand_z", &LocalProtCand_z );
          tree->Branch( "LocalProtCand_xSigma", &LocalProtCand_xSigma );
          tree->Branch( "LocalProtCand_ySigma", &LocalProtCand_ySigma );
          tree->Branch( "LocalProtCand_arm", &LocalProtCand_arm );
          tree->Branch( "LocalProtCand_pot", &LocalProtCand_pot );
          tree->Branch( "LocalProtCand_Tx", &LocalProtCand_Tx );
          tree->Branch( "LocalProtCand_Ty", &LocalProtCand_Ty );
          tree->Branch( "LocalProtCand_TxSigma", &LocalProtCand_TxSigma );
          tree->Branch( "LocalProtCand_TySigma", &LocalProtCand_TySigma );
        }

        // Extra tracks on vertex's information
        tree->Branch( "ExtraTrack_pair", &ExtraTrack_pair );
        tree->Branch( "ExtraTrack_purity", &ExtraTrack_purity );
        tree->Branch( "ExtraTrack_nhits", &ExtraTrack_nhits );
        tree->Branch( "ExtraTrack_charge", &ExtraTrack_charge );
        tree->Branch( "ExtraTrack_ndof", &ExtraTrack_ndof );
        tree->Branch( "ExtraTrack_px", &ExtraTrack_px );
        tree->Branch( "ExtraTrack_py", &ExtraTrack_py );
        tree->Branch( "ExtraTrack_pz", &ExtraTrack_pz );
        tree->Branch( "ExtraTrack_chi2", &ExtraTrack_chi2 );
        tree->Branch( "ExtraTrack_vtxdxyz", &ExtraTrack_vtxdxyz );
        tree->Branch( "ExtraTrack_vtxT", &ExtraTrack_vtxT );
        tree->Branch( "ExtraTrack_vtxZ", &ExtraTrack_vtxZ );
        tree->Branch( "ExtraTrack_x", &ExtraTrack_x );
        tree->Branch( "ExtraTrack_y", &ExtraTrack_y );
        tree->Branch( "ExtraTrack_z", &ExtraTrack_z );
        tree->Branch( "nQualityExtraTrack", &nQualityExtraTrack );

        // indexed by pair id
        tree->Branch( "ClosestExtraTrack_vtxdxyz", &ClosestExtraTrack_vtxdxyz );
        tree->Branch( "ClosestExtraTrack_id", &ClosestExtraTrack_id );
        tree->Branch( "ClosestHighPurityExtraTrack_vtxdxyz", &ClosestHighPurityExtraTrack_vtxdxyz );
        tree->Branch( "ClosestHighPurityExtraTrack_id", &ClosestHighPurityExtraTrack_id );

        // Jets/MET information
        tree->Branch( "JetCand_pt", &JetCand_pt );
        tree->Branch( "JetCand_eta", &JetCand_eta );
        tree->Branch( "JetCand_phi", &JetCand_phi );
        tree->Branch( "JetCand_e", &JetCand_e );
        tree->Branch( "HighestJet_pt", &HighestJet_pt );
        tree->Branch( "HighestJet_eta", &HighestJet_eta );
        tree->Branch( "HighestJet_phi", &HighestJet_phi );
        tree->Branch( "HighestJet_e", &HighestJet_e );
        tree->Branch( "SumJet_e", &SumJet_e );
        tree->Branch( "Etmiss", &Etmiss );
        tree->Branch( "Etmiss_phi", &Etmiss_phi );
        tree->Branch( "Etmiss_significance", &Etmiss_significance );

        // Pileup reweighting
        tree->Branch( "Weight", &Weight );
        tree->Branch( "PUWeightTrue", &PUWeightTrue );
      }
      void load( TTree* tree, TreeType tt, bool mc ) {
        if ( !tree ) return;

        tree->SetBranchAddress( "Run", &Run );
        tree->SetBranchAddress( "LumiSection", &LumiSection );
        tree->SetBranchAddress( "BX", &BX );
        tree->SetBranchAddress( "EventNum", &EventNum );

        tree->SetBranchAddress( "HLT_Accept", &HLT_Accept );
        tree->SetBranchAddress( "HLT_Prescl", &HLT_Prescl );
        tree->SetBranchAddress( "HLT_Name", &HLT_Name );

        if ( tt == ElectronMuon || tt == DiMuon ) {
          tree->SetBranchAddress( "MuonCand_pt", &MuonCand_pt );
          tree->SetBranchAddress( "MuonCand_eta", &MuonCand_eta );
          tree->SetBranchAddress( "MuonCand_phi", &MuonCand_phi );
          tree->SetBranchAddress( "MuonCand_e", &MuonCand_e );
          tree->SetBranchAddress( "MuonCand_charge", &MuonCand_charge );
          tree->SetBranchAddress( "MuonCand_vtxx", &MuonCand_vtxx );
          tree->SetBranchAddress( "MuonCand_vtxy", &MuonCand_vtxy );
          tree->SetBranchAddress( "MuonCand_vtxz", &MuonCand_vtxz );
          tree->SetBranchAddress( "MuonCand_dxy", &MuonCand_dxy );
          tree->SetBranchAddress( "MuonCand_nstatseg", &MuonCand_nstatseg );
          tree->SetBranchAddress( "MuonCand_ntrklayers", &MuonCand_ntrklayers );
          tree->SetBranchAddress( "MuonCand_npxlhits", &MuonCand_npxlhits );
          tree->SetBranchAddress( "MuonCand_isglobal", &MuonCand_isglobal );
          tree->SetBranchAddress( "MuonCand_istracker", &MuonCand_istracker );
          tree->SetBranchAddress( "MuonCand_isstandalone", &MuonCand_isstandalone );
          tree->SetBranchAddress( "MuonCand_ispfmuon", &MuonCand_ispfmuon );
          tree->SetBranchAddress( "MuonCand_istight", &MuonCand_istight );
          tree->SetBranchAddress( "MuonCandTrack_nmuchits", &MuonCandTrack_nmuchits );
          tree->SetBranchAddress( "MuonCandTrack_chisq", &MuonCandTrack_chisq );
          tree->SetBranchAddress( "MuonCand_innerTrackPt", &MuonCand_innerTrackPt );
          tree->SetBranchAddress( "MuonCand_innerTrackEta", &MuonCand_innerTrackEta );
          tree->SetBranchAddress( "MuonCand_innerTrackPhi", &MuonCand_innerTrackPhi );
          tree->SetBranchAddress( "MuonCand_innerTrackVtxz", &MuonCand_innerTrackVtxz );
          if ( mc ) {
            tree->SetBranchAddress( "nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept );
            tree->SetBranchAddress( "GenMuonCand_pt", &GenMuonCand_pt );
            tree->SetBranchAddress( "GenMuonCand_eta", &GenMuonCand_eta );
            tree->SetBranchAddress( "GenMuonCand_phi", &GenMuonCand_phi );
            tree->SetBranchAddress( "GenMuonCand_e", &GenMuonCand_e );
          }
        }

        if ( tt == ElectronMuon || tt == DiElectron ) {
          tree->SetBranchAddress( "EleCand_et", &EleCand_et );
          tree->SetBranchAddress( "EleCand_eta", &EleCand_eta );
          tree->SetBranchAddress( "EleCand_phi", &EleCand_phi );
          tree->SetBranchAddress( "EleCand_e", &EleCand_e );
          tree->SetBranchAddress( "EleCand_charge", &EleCand_charge );
          tree->SetBranchAddress( "EleCand_vtxx", &EleCand_vtxx );
          tree->SetBranchAddress( "EleCand_vtxy", &EleCand_vtxy );
          tree->SetBranchAddress( "EleCand_vtxz", &EleCand_vtxz );
          tree->SetBranchAddress( "EleCand_deltaPhi", &EleCand_deltaPhi );
          tree->SetBranchAddress( "EleCand_deltaEta", &EleCand_deltaEta );
          tree->SetBranchAddress( "EleCand_HoverE", &EleCand_HoverE );
          tree->SetBranchAddress( "EleCand_trackiso", &EleCand_trackiso );
          tree->SetBranchAddress( "EleCand_ecaliso", &EleCand_ecaliso );
          tree->SetBranchAddress( "EleCand_hcaliso", &EleCand_hcaliso );
          tree->SetBranchAddress( "EleCand_sigmaIetaIeta", &EleCand_sigmaIetaIeta );
          tree->SetBranchAddress( "EleCand_convDist", &EleCand_convDist );
          tree->SetBranchAddress( "EleCand_convDcot", &EleCand_convDcot );
          tree->SetBranchAddress( "EleCand_ecalDriven", &EleCand_ecalDriven );
          tree->SetBranchAddress( "EleCand_mediumID", &EleCand_mediumID );
          tree->SetBranchAddress( "EleCand_tightID", &EleCand_tightID );
          tree->SetBranchAddress( "EleCand_innerTrackPt", &EleCand_innerTrackPt );
          tree->SetBranchAddress( "EleCand_innerTrackEta", &EleCand_innerTrackEta );
          tree->SetBranchAddress( "EleCand_innerTrackPhi", &EleCand_innerTrackPhi );
          tree->SetBranchAddress( "EleCand_innerTrackVtxz", &EleCand_innerTrackVtxz );
          if ( mc ) {
            tree->SetBranchAddress( "nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept );
            tree->SetBranchAddress( "GenEleCand_pt", &GenEleCand_pt );
            tree->SetBranchAddress( "GenEleCand_eta", &GenEleCand_eta );
            tree->SetBranchAddress( "GenEleCand_phi", &GenEleCand_phi );
            tree->SetBranchAddress( "GenEleCand_e", &GenEleCand_e );
          }
        }
        tree->SetBranchAddress( "PhotonCand_pt", &PhotonCand_pt );
        tree->SetBranchAddress( "PhotonCand_eta", &PhotonCand_eta );
        tree->SetBranchAddress( "PhotonCand_phi", &PhotonCand_phi );
        tree->SetBranchAddress( "PhotonCand_e", &PhotonCand_e );
        tree->SetBranchAddress( "PhotonCand_r9", &PhotonCand_r9 );
        tree->SetBranchAddress( "PhotonCand_drtrue", &PhotonCand_drtrue );
        tree->SetBranchAddress( "PhotonCand_detatrue", &PhotonCand_detatrue );
        tree->SetBranchAddress( "PhotonCand_dphitrue", &PhotonCand_dphitrue );
        tree->SetBranchAddress( "PhotonCand_mediumID", &PhotonCand_mediumID );
        tree->SetBranchAddress( "PhotonCand_tightID", &PhotonCand_tightID );
        if ( mc ) {
          tree->SetBranchAddress( "nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept );
          tree->SetBranchAddress( "GenPhotCand_pt", &GenPhotCand_pt );
          tree->SetBranchAddress( "GenPhotCand_eta", &GenPhotCand_eta );
          tree->SetBranchAddress( "GenPhotCand_phi", &GenPhotCand_phi );
          tree->SetBranchAddress( "GenPhotCand_e", &GenPhotCand_e );

          tree->SetBranchAddress( "GenProtCand_pt", &GenProtCand_pt );
          tree->SetBranchAddress( "GenProtCand_eta", &GenProtCand_eta );
          tree->SetBranchAddress( "GenProtCand_phi", &GenProtCand_phi );
          tree->SetBranchAddress( "GenProtCand_e", &GenProtCand_e );
          tree->SetBranchAddress( "GenProtCand_status", &GenProtCand_status );
        }

        // Primary vertices' information
        tree->SetBranchAddress( "nFilteredPrimVertexCand", &nFilteredPrimVertexCand );
        tree->SetBranchAddress( "PrimVertexCand_id", &PrimVertexCand_id );
        tree->SetBranchAddress( "PrimVertexCand_x", &PrimVertexCand_x );
        tree->SetBranchAddress( "PrimVertexCand_y", &PrimVertexCand_y );
        tree->SetBranchAddress( "PrimVertexCand_z", &PrimVertexCand_z );
        tree->SetBranchAddress( "PrimVertexCand_chi2", &PrimVertexCand_chi2 );
        tree->SetBranchAddress( "PrimVertexCand_ndof", &PrimVertexCand_ndof );
        tree->SetBranchAddress( "PrimVertexCand_tracks", &PrimVertexCand_tracks );

        // Lepton pairs' information
        tree->SetBranchAddress( "Pair_lepton1", &Pair_lepton1 );
        tree->SetBranchAddress( "Pair_lepton2", &Pair_lepton2 );
        //tree->SetBranchAddress( "Pair_mindist", &Pair_mindist );
        tree->SetBranchAddress( "Pair_mass", &Pair_mass );
        tree->SetBranchAddress( "Pair_pt", &Pair_pt );
        tree->SetBranchAddress( "Pair_eta", &Pair_eta );
        tree->SetBranchAddress( "Pair_phi", &Pair_phi );
        tree->SetBranchAddress( "Pair_dpt", &Pair_dpt );
        tree->SetBranchAddress( "Pair_dphi", &Pair_dphi );
        tree->SetBranchAddress( "Pair_3Dangle", &Pair_3Dangle );
        tree->SetBranchAddress( "Pair_extratracks0p5mm", &Pair_extratracks0p5mm );
        tree->SetBranchAddress( "Pair_extratracks1mm", &Pair_extratracks1mm );
        tree->SetBranchAddress( "Pair_extratracks2mm", &Pair_extratracks2mm );
        tree->SetBranchAddress( "Pair_extratracks3mm", &Pair_extratracks3mm );
        tree->SetBranchAddress( "Pair_extratracks4mm", &Pair_extratracks4mm );
        tree->SetBranchAddress( "Pair_extratracks5mm", &Pair_extratracks5mm );
        tree->SetBranchAddress( "Pair_extratracks1cm", &Pair_extratracks1cm );
        tree->SetBranchAddress( "Pair_extratracks2cm", &Pair_extratracks2cm );
        tree->SetBranchAddress( "Pair_extratracks3cm", &Pair_extratracks3cm );
        tree->SetBranchAddress( "Pair_extratracks4cm", &Pair_extratracks4cm );
        tree->SetBranchAddress( "Pair_extratracks5cm", &Pair_extratracks5cm );
        tree->SetBranchAddress( "Pair_extratracks10cm", &Pair_extratracks10cm );
        // Kalman dilepton vertex information
        tree->SetBranchAddress( "KalmanVertexCand_x", &KalmanVertexCand_x );
        tree->SetBranchAddress( "KalmanVertexCand_y", &KalmanVertexCand_y );
        tree->SetBranchAddress( "KalmanVertexCand_z", &KalmanVertexCand_z );

        tree->SetBranchAddress( "PairGamma_pair", &PairGamma_pair );
        tree->SetBranchAddress( "PairGamma_pho", &PairGamma_pho );
        tree->SetBranchAddress( "PairGamma_mass", &PairGamma_mass );
        if ( mc ) {
          tree->SetBranchAddress( "GenPair_mass", &GenPair_mass );
          tree->SetBranchAddress( "GenPair_pt", &GenPair_pt );
          tree->SetBranchAddress( "GenPair_eta", &GenPair_eta );
          tree->SetBranchAddress( "GenPair_phi", &GenPair_phi );
          tree->SetBranchAddress( "GenPair_dpt", &GenPair_dpt );
          tree->SetBranchAddress( "GenPair_dphi", &GenPair_dphi );
          tree->SetBranchAddress( "GenPair_3Dangle", &GenPair_3Dangle );
        }

        if ( !mc ) {
          tree->SetBranchAddress( "LocalProtCand_x", &LocalProtCand_x );
          tree->SetBranchAddress( "LocalProtCand_y", &LocalProtCand_y );
          tree->SetBranchAddress( "LocalProtCand_z", &LocalProtCand_z );
          tree->SetBranchAddress( "LocalProtCand_xSigma", &LocalProtCand_xSigma );
          tree->SetBranchAddress( "LocalProtCand_ySigma", &LocalProtCand_ySigma );
          tree->SetBranchAddress( "LocalProtCand_arm", &LocalProtCand_arm );
          tree->SetBranchAddress( "LocalProtCand_pot", &LocalProtCand_pot );
          tree->SetBranchAddress( "LocalProtCand_Tx", &LocalProtCand_Tx );
          tree->SetBranchAddress( "LocalProtCand_Ty", &LocalProtCand_Ty );
          tree->SetBranchAddress( "LocalProtCand_TxSigma", &LocalProtCand_TxSigma );
          tree->SetBranchAddress( "LocalProtCand_TySigma", &LocalProtCand_TySigma );
        }

        // Extra tracks on vertex's information
        tree->SetBranchAddress( "ExtraTrack_pair", &ExtraTrack_pair );
        tree->SetBranchAddress( "ExtraTrack_purity", &ExtraTrack_purity );
        tree->SetBranchAddress( "ExtraTrack_nhits", &ExtraTrack_nhits );
        tree->SetBranchAddress( "ExtraTrack_charge", &ExtraTrack_charge );
        tree->SetBranchAddress( "ExtraTrack_ndof", &ExtraTrack_ndof );
        tree->SetBranchAddress( "ExtraTrack_px", &ExtraTrack_px );
        tree->SetBranchAddress( "ExtraTrack_py", &ExtraTrack_py );
        tree->SetBranchAddress( "ExtraTrack_pz", &ExtraTrack_pz );
        tree->SetBranchAddress( "ExtraTrack_chi2", &ExtraTrack_chi2 );
        tree->SetBranchAddress( "ExtraTrack_vtxdxyz", &ExtraTrack_vtxdxyz );
        tree->SetBranchAddress( "ExtraTrack_vtxT", &ExtraTrack_vtxT );
        tree->SetBranchAddress( "ExtraTrack_vtxZ", &ExtraTrack_vtxZ );
        tree->SetBranchAddress( "ExtraTrack_x", &ExtraTrack_x );
        tree->SetBranchAddress( "ExtraTrack_y", &ExtraTrack_y );
        tree->SetBranchAddress( "ExtraTrack_z", &ExtraTrack_z );
        tree->SetBranchAddress( "nQualityExtraTrack", &nQualityExtraTrack );
        tree->SetBranchAddress( "ClosestExtraTrack_vtxdxyz", &ClosestExtraTrack_vtxdxyz );
        tree->SetBranchAddress( "ClosestExtraTrack_id", &ClosestExtraTrack_id );
        tree->SetBranchAddress( "ClosestHighPurityExtraTrack_vtxdxyz", &ClosestHighPurityExtraTrack_vtxdxyz );
        tree->SetBranchAddress( "ClosestHighPurityExtraTrack_id", &ClosestHighPurityExtraTrack_id );

        // Jets/MET information
        tree->SetBranchAddress( "JetCand_pt", &JetCand_pt );
        tree->SetBranchAddress( "JetCand_eta", &JetCand_eta );
        tree->SetBranchAddress( "JetCand_phi", &JetCand_phi );
        tree->SetBranchAddress( "JetCand_e", &JetCand_e );
        tree->SetBranchAddress( "HighestJet_pt", &HighestJet_pt );
        tree->SetBranchAddress( "HighestJet_eta", &HighestJet_eta );
        tree->SetBranchAddress( "HighestJet_phi", &HighestJet_phi );
        tree->SetBranchAddress( "HighestJet_e", &HighestJet_e );
        tree->SetBranchAddress( "SumJet_e", &SumJet_e );
        tree->SetBranchAddress( "Etmiss", &Etmiss );
        tree->SetBranchAddress( "Etmiss_phi", &Etmiss_phi );
        tree->SetBranchAddress( "Etmiss_significance", &Etmiss_significance );

        // Pileup reweighting
        tree->SetBranchAddress( "Weight", &Weight );
        tree->SetBranchAddress( "PUWeightTrue", &PUWeightTrue );
      }
  };
}

#endif

