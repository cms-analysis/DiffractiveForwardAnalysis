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

      //////  Leaves size  //////

      /// Maximum number of HLT to check
      static constexpr unsigned int MAX_HLT = 10;
      /// Maximum number of leptons per event
      static constexpr unsigned int MAX_LL = 50;
      /// Maximum number of muons per event
      static constexpr unsigned int MAX_MUONS = 25;
      /// Maximum number of electrons per event
      static constexpr unsigned int MAX_ELE = 25;
      /// Maximum number of photons per event
      static constexpr unsigned int MAX_PHO = 50;
      /// Maximum number of leptons pairs per event
      static constexpr unsigned int MAX_PAIRS = 25;
      static constexpr unsigned int MAX_PAIRPHO = 25;
      /// Maximum number of primary vertices per event
      static constexpr unsigned int MAX_VTX = 150;
      /// Maximum number of extra tracks per event
      static constexpr unsigned int MAX_ET = 1000;
      /// Maximum number of generator level muons per event
      static constexpr unsigned int MAX_GENMU = 25;
      /// Maximum number of generator level electrons per event
      static constexpr unsigned int MAX_GENELE = 25;
      /// Maximum number of generator level photons per event
      static constexpr unsigned int MAX_GENPHO = 10;
      /// Maximum number of generator level protons per event
      static constexpr unsigned int MAX_GENPRO = 8;
      /// Maximum number of jets per event
      static constexpr unsigned int MAX_JETS = 40;
      /// Maximum number of reconstructed local tracks in RPs
      static constexpr unsigned int MAX_LOCALPCAND = 25;
      /// Maximum number of reconstructed local tracks pairs in RPs
      static constexpr unsigned int MAX_LOCALPPAIRCAND = 10;

      ////// Tree contents //////

      // Run/event quantities
      unsigned int BX, Run, LumiSection, EventNum;
      //int LHCFillNum, LHCBeamMode;
      //double AvgInstDelLumi, BunchInstLumi[3];

      // HLT quantities
      unsigned int nHLT;
      int HLT_Accept[MAX_HLT], HLT_Prescl[MAX_HLT];
      std::vector<std::string>* HLT_Name;
      /*int nHLTLeptonCand[MAX_HLT];
      double HLTLeptonCand_pt[2][MAX_HLT];
      double HLTLeptonCand_eta[2][MAX_HLT];
      double HLTLeptonCand_phi[2][MAX_HLT];
      int HLTLeptonCand_charge[2][MAX_HLT];
      int HLT_LeadingLepton[MAX_HLT], HLT_TrailingLepton[MAX_HLT];
      int HLT_LeadingLepton_Prescl[MAX_HLT], HLT_TrailingLepton_Prescl[MAX_HLT];*/

      // Generator level quantities
      unsigned int nGenMuonCand, nGenMuonCandOutOfAccept;
      double GenMuonCand_pt[MAX_GENMU], GenMuonCand_eta[MAX_GENMU], GenMuonCand_phi[MAX_GENMU], GenMuonCand_e[MAX_GENMU];
      unsigned int nGenEleCand, nGenEleCandOutOfAccept;
      double GenEleCand_pt[MAX_GENELE], GenEleCand_eta[MAX_GENELE], GenEleCand_phi[MAX_GENELE], GenEleCand_e[MAX_GENELE];
      double GenPair_pt, GenPair_eta, GenPair_phi, GenPair_mass;
      double GenPair_dphi, GenPair_dpt, GenPair_3Dangle;
      unsigned int nGenPhotCand, nGenPhotCandOutOfAccept;
      double GenPhotCand_pt[MAX_GENPHO], GenPhotCand_eta[MAX_GENPHO], GenPhotCand_phi[MAX_GENPHO], GenPhotCand_e[MAX_GENPHO];
      unsigned int nGenProtCand;
      double GenProtCand_pt[MAX_GENPRO], GenProtCand_eta[MAX_GENPRO], GenProtCand_phi[MAX_GENPRO], GenProtCand_e[MAX_GENPHO];
      int GenProtCand_status[MAX_GENPRO];

      // Pileup reweighting quantities
      double PUWeightTrue, Weight;

      // Muon quantities
      unsigned int nMuonCand;
      double MuonCand_pt[MAX_LL], MuonCand_eta[MAX_LL], MuonCand_phi[MAX_LL], MuonCand_e[MAX_LL];
      double MuonCand_innerTrackPt[MAX_LL], MuonCand_innerTrackEta[MAX_LL], MuonCand_innerTrackPhi[MAX_LL];
      double MuonCand_innerTrackVtxz[MAX_LL];
      double MuonCand_vtxx[MAX_LL], MuonCand_vtxy[MAX_LL], MuonCand_vtxz[MAX_LL];
      int MuonCand_charge[MAX_LL];
      double MuonCand_dxy[MAX_LL];
      int MuonCand_nstatseg[MAX_LL], MuonCand_npxlhits[MAX_LL], MuonCand_ntrklayers[MAX_LL];
      int MuonCandTrack_nmuchits[MAX_LL];
      double MuonCandTrack_chisq[MAX_LL];
      int MuonCand_isglobal[MAX_LL], MuonCand_istracker[MAX_LL], MuonCand_isstandalone[MAX_LL], MuonCand_ispfmuon[MAX_LL];
      int MuonCand_istight[MAX_LL];

      // Electron quantities
      unsigned int nEleCand;
      double EleCand_et[MAX_LL], EleCand_eta[MAX_LL], EleCand_phi[MAX_LL], EleCand_e[MAX_LL];
      double EleCand_vtxx[MAX_LL], EleCand_vtxy[MAX_LL], EleCand_vtxz[MAX_LL];
      double EleCand_innerTrackPt[MAX_LL], EleCand_innerTrackEta[MAX_LL], EleCand_innerTrackPhi[MAX_LL];
      double EleCand_innerTrackVtxz[MAX_LL];
      int EleCand_charge[MAX_LL];
      double EleCand_deltaPhi[MAX_LL], EleCand_deltaEta[MAX_LL];
      double EleCand_HoverE[MAX_LL];
      double EleCand_trackiso[MAX_LL], EleCand_ecaliso[MAX_LL], EleCand_hcaliso[MAX_LL];
      double EleCand_sigmaIetaIeta[MAX_LL];
      double EleCand_convDist[MAX_LL], EleCand_convDcot[MAX_LL];
      int EleCand_ecalDriven[MAX_LL];
      int EleCand_tightID[MAX_LL], EleCand_mediumID[MAX_LL];

      // Photon quantities
      unsigned int nPhotonCand;
      double PhotonCand_pt[MAX_PHO], PhotonCand_eta[MAX_PHO], PhotonCand_phi[MAX_PHO], PhotonCand_e[MAX_PHO];
      double PhotonCand_r9[MAX_PHO];
      double PhotonCand_drtrue[MAX_PHO], PhotonCand_detatrue[MAX_PHO], PhotonCand_dphitrue[MAX_PHO];
      int PhotonCand_tightID[MAX_PHO], PhotonCand_mediumID[MAX_PHO];

      // Pair quantities
      unsigned int nPair;
      int Pair_lepton1[MAX_PAIRS], Pair_lepton2[MAX_PAIRS];
      double Pair_pt[MAX_PAIRS], Pair_eta[MAX_PAIRS], Pair_phi[MAX_PAIRS], Pair_mass[MAX_PAIRS];
      double Pair_dpt[MAX_PAIRS], Pair_dphi[MAX_PAIRS], Pair_3Dangle[MAX_PAIRS];
      //double Pair_mindist[MAX_PAIRS];

      unsigned int nPairGamma;
      int PairGamma_pair[MAX_PHO];
      double PairGamma_mass[MAX_PHO];

      // Extra tracks
      unsigned int Pair_extratracks0p5mm[MAX_PAIRS];
      unsigned int Pair_extratracks1mm[MAX_PAIRS], Pair_extratracks2mm[MAX_PAIRS];
      unsigned int Pair_extratracks3mm[MAX_PAIRS], Pair_extratracks4mm[MAX_PAIRS];
      unsigned int Pair_extratracks5mm[MAX_PAIRS], Pair_extratracks1cm[MAX_PAIRS];
      unsigned int Pair_extratracks2cm[MAX_PAIRS], Pair_extratracks3cm[MAX_PAIRS];
      unsigned int Pair_extratracks4cm[MAX_PAIRS], Pair_extratracks5cm[MAX_PAIRS];
      unsigned int Pair_extratracks10cm[MAX_PAIRS];

      // Vertex quantities
      unsigned int nPrimVertexCand;
      int PrimVertexCand_id[MAX_VTX], PrimVertexCand_hasdil[MAX_VTX];
      double PrimVertexCand_x[MAX_VTX], PrimVertexCand_y[MAX_VTX], PrimVertexCand_z[MAX_VTX];
      unsigned int PrimVertexCand_tracks[MAX_VTX];
      double PrimVertexCand_chi2[MAX_VTX];
      unsigned int PrimVertexCand_ndof[MAX_VTX];
      double KalmanVertexCand_x[MAX_VTX], KalmanVertexCand_y[MAX_VTX], KalmanVertexCand_z[MAX_VTX];
      unsigned int nFilteredPrimVertexCand;

      // Extra tracks on vertex quantities
      unsigned int nExtraTracks;
      int ExtraTrack_pair[MAX_ET];
      int ExtraTrack_purity[MAX_ET];
      unsigned int ExtraTrack_nhits[MAX_ET];
      int ExtraTrack_charge[MAX_ET];
      unsigned int ExtraTrack_ndof[MAX_ET];
      double ExtraTrack_px[MAX_ET], ExtraTrack_py[MAX_ET], ExtraTrack_pz[MAX_ET];
      double ExtraTrack_chi2[MAX_ET];
      double ExtraTrack_vtxdxyz[MAX_ET];
      double ExtraTrack_vtxT[MAX_ET], ExtraTrack_vtxZ[MAX_ET];
      double ExtraTrack_x[MAX_ET], ExtraTrack_y[MAX_ET], ExtraTrack_z[MAX_ET];
      double ClosestExtraTrack_vtxdxyz[MAX_PAIRS], ClosestHighPurityExtraTrack_vtxdxyz[MAX_PAIRS];
      int ClosestExtraTrack_id[MAX_PAIRS], ClosestHighPurityExtraTrack_id[MAX_PAIRS];
      unsigned int nQualityExtraTrack;

      // Jets/MET quantities
      unsigned int nJetCand;
      double JetCand_pt[MAX_JETS], JetCand_eta[MAX_JETS], JetCand_phi[MAX_JETS], JetCand_e[MAX_JETS];
      double HighestJet_pt, HighestJet_eta, HighestJet_phi, HighestJet_e;
      double SumJet_e;
      double Etmiss, Etmiss_phi, Etmiss_significance;

      // CTPPS quantities
      unsigned int nLocalProtCand;
      double LocalProtCand_x[MAX_LOCALPCAND], LocalProtCand_y[MAX_LOCALPCAND], LocalProtCand_z[MAX_LOCALPCAND];
      double LocalProtCand_xSigma[MAX_LOCALPCAND], LocalProtCand_ySigma[MAX_LOCALPCAND];
      double LocalProtCand_Tx[MAX_LOCALPCAND], LocalProtCand_Ty[MAX_LOCALPCAND];
      double LocalProtCand_TxSigma[MAX_LOCALPCAND], LocalProtCand_TySigma[MAX_LOCALPCAND];
      int LocalProtCand_arm[MAX_LOCALPCAND], LocalProtCand_pot[MAX_LOCALPCAND];

      void clear() {
        // event-level branches
        BX = Run = LumiSection = EventNum = 0;

        // high-level trigger
        nHLT = 0;
        for ( unsigned int i = 0; i < MAX_HLT; ++i ) {
          HLT_Accept[i] = HLT_Prescl[i] = -1;
        }

        // gen-level information
        nGenMuonCand = nGenMuonCandOutOfAccept = 0;
        for ( unsigned int i = 0; i < MAX_GENMU; ++i ) {
          GenMuonCand_pt[i] = GenMuonCand_eta[i] = GenMuonCand_phi[i] = GenMuonCand_e[i] = -999.;
        }
        nGenEleCand = nGenEleCandOutOfAccept = 0;
        for ( unsigned int i = 0; i < MAX_GENELE; ++i ) {
          GenEleCand_pt[i] = GenEleCand_eta[i] = GenEleCand_phi[i] = GenEleCand_e[i] = -999.;
        }
        GenPair_pt = GenPair_eta = GenPair_phi = GenPair_mass = -999.;
        GenPair_dphi = GenPair_dpt = GenPair_3Dangle = -999.;
        nGenPhotCand = nGenPhotCandOutOfAccept = 0;
        for ( unsigned int i = 0; i < MAX_GENPHO; ++i ) {
          GenPhotCand_pt[i] = GenPhotCand_eta[i] = GenPhotCand_phi[i] = GenPhotCand_e[i] = -999.;
        }
        nGenProtCand = 0;
        for ( unsigned int i = 0; i < MAX_GENPRO; ++i ) {
          GenProtCand_pt[i] = GenProtCand_eta[i] = GenProtCand_phi[i] = GenProtCand_e[i] = -999.;
          GenProtCand_status[i] = -1;
        }

        PUWeightTrue = Weight = 0.;

        //LHCFillNum = LHCBeamMode = -1;

        // single lepton candidates
        nMuonCand = nEleCand = 0;
        for ( unsigned int i = 0; i < MAX_LL; ++i ) {
          MuonCand_pt[i] = MuonCand_eta[i] = MuonCand_phi[i] = MuonCand_e[i] = -999.;
          MuonCand_innerTrackPt[i] = MuonCand_innerTrackEta[i] = MuonCand_innerTrackPhi[i] = -999.;
          MuonCand_innerTrackVtxz[i] = -999.;
          MuonCand_vtxx[i] = MuonCand_vtxy[i] = MuonCand_vtxz[i] = -999.;
          MuonCand_charge[i] = 0;
          MuonCand_dxy[i] = -999.;
          MuonCand_nstatseg[i] = MuonCand_npxlhits[i] = MuonCand_ntrklayers[i] = -999;
          MuonCandTrack_nmuchits[i] = -999;
          MuonCandTrack_chisq[i] = -999.;
          MuonCand_isglobal[i] = MuonCand_istracker[i] = MuonCand_isstandalone[i] = MuonCand_ispfmuon[i] = -999;
          MuonCand_istight[i] = -999;
          EleCand_et[i] = EleCand_eta[i] = EleCand_phi[i] = EleCand_e[i] = -999.;
          EleCand_vtxx[i] = EleCand_vtxy[i] = EleCand_vtxz[i] = -999.;
          EleCand_innerTrackPt[i] = EleCand_innerTrackEta[i] = EleCand_innerTrackPhi[i] = -999.;
          EleCand_innerTrackVtxz[i] = -999.;
          EleCand_charge[i] = 0;
          EleCand_deltaPhi[i] = EleCand_deltaEta[i] = -999.;
          EleCand_HoverE[i] = -999.;
          EleCand_trackiso[i] = EleCand_ecaliso[i] = EleCand_hcaliso[i] = -999.;
          EleCand_sigmaIetaIeta[i] = -999.;
          EleCand_convDist[i] = EleCand_convDcot[i] = -999.;
          EleCand_ecalDriven[i] = -999;
          EleCand_tightID[i] = EleCand_mediumID[i] = -1;
        }

        // single photon candidates
        nPhotonCand = 0;
        for ( unsigned int i = 0; i < MAX_PHO; ++i ) {
          PhotonCand_pt[i] = PhotonCand_eta[i] = PhotonCand_phi[i] = PhotonCand_e[i] = -999.;
          PhotonCand_r9[i] = -999.;
          PhotonCand_drtrue[i] = PhotonCand_detatrue[i] = PhotonCand_dphitrue[i] = -999.;
          PhotonCand_tightID[i] = PhotonCand_mediumID[i] = -1;
        }

        // dilepton pair candidates
        nPair = 0;
        for ( unsigned int i = 0; i < MAX_PAIRS; ++i ) {
          Pair_lepton1[i] = Pair_lepton2[i] = -1;
          Pair_pt[i] = Pair_mass[i] = Pair_phi[i] = Pair_eta[i] = -999.;
          Pair_dpt[i] = Pair_dphi[i] = Pair_3Dangle[i] = -999.;
          //Pair_mindist[i] = -999.;
        }

        // dilepton pair + associated photon candidates
        nPairGamma = 0;
        for ( unsigned int i = 0; i < MAX_PHO; ++i ) {
          PairGamma_pair[i] = -1;
          PairGamma_mass[i] = -999.;
        }

        // extra tracks associated to the central system vertex
        for ( unsigned int i = 0; i < MAX_PAIRS; ++i ) {
          Pair_extratracks0p5mm[i] = 0;
          Pair_extratracks1mm[i] = Pair_extratracks2mm[i] = Pair_extratracks3mm[i] = 0;
          Pair_extratracks4mm[i] = Pair_extratracks5mm[i] = Pair_extratracks1cm[i] = 0;
          Pair_extratracks2cm[i] = Pair_extratracks3cm[i] = Pair_extratracks4cm[i] = 0;
          Pair_extratracks5cm[i] = Pair_extratracks10cm[i] = 0;
        }

        // offline primary vertices
        nPrimVertexCand = 0;
        for ( unsigned int i = 0; i < MAX_VTX; ++i ) {
          PrimVertexCand_id[i] = PrimVertexCand_hasdil[i] = -1;
          PrimVertexCand_x[i] = PrimVertexCand_y[i] = PrimVertexCand_z[i] = -999.;
          PrimVertexCand_tracks[i] = 0;
          PrimVertexCand_chi2[i] = -999.;
          PrimVertexCand_ndof[i] = 0;
          KalmanVertexCand_x[i] = KalmanVertexCand_y[i] = KalmanVertexCand_z[i] = -999.;
        }
        nFilteredPrimVertexCand = 0;

        // extra tracks associated to the central system
        nExtraTracks = 0;
        for ( unsigned int i = 0; i < MAX_ET; ++i ) {
          ExtraTrack_pair[i] = -1;
          ExtraTrack_purity[i] = ExtraTrack_nhits[i] = -1;
          ExtraTrack_charge[i] = -999;
          ExtraTrack_ndof[i] = 0;
          ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
          ExtraTrack_chi2[i] = ExtraTrack_vtxdxyz[i] = -999.;
          ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
          ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
        }
        for ( unsigned int i = 0; i < MAX_PAIRS; ++i ) {
          ClosestExtraTrack_vtxdxyz[i] = ClosestHighPurityExtraTrack_vtxdxyz[i] = -999.;
          ClosestExtraTrack_id[i] = ClosestHighPurityExtraTrack_id[i] = -1;
        }
        nQualityExtraTrack = 0;

        // jets collection
        nJetCand = 0;
        for ( unsigned int i = 0; i < MAX_JETS; ++i ) {
          JetCand_pt[i] = JetCand_eta[i] = JetCand_phi[i] = JetCand_e[i] = -999;
        }
        HighestJet_pt = HighestJet_eta = HighestJet_phi = HighestJet_e = -999.;
        SumJet_e = 0.;

        // missing ET
        Etmiss = Etmiss_phi = Etmiss_significance = -999.;

        // CTPPS strips leaves
        nLocalProtCand = 0;
        for ( unsigned int i = 0; i < MAX_LOCALPCAND; ++i ) {
          LocalProtCand_x[i] = LocalProtCand_y[i] = LocalProtCand_z[i] = -999.;
          LocalProtCand_xSigma[i] = LocalProtCand_ySigma[i] = -999.;
          LocalProtCand_Tx[i] = LocalProtCand_Ty[i] = -999.;
          LocalProtCand_TxSigma[i] = LocalProtCand_TySigma[i] = -999.;
          LocalProtCand_arm[i] = LocalProtCand_pot[i] = -1;
        }
      }
      void attach( TTree* tree, TreeType tt, bool mc ) {
        if ( !tree ) return;

        tree->Branch( "Run", &Run, "Run/i" );
        tree->Branch( "LumiSection", &LumiSection, "LumiSection/i" );
        tree->Branch( "BX", &BX, "BX/i" );
        tree->Branch( "EventNum", &EventNum, "EventNum/i" );

        tree->Branch( "nHLT", &nHLT, "nHLT/i" );
        tree->Branch( "HLT_Accept", HLT_Accept, "HLT_Accept[nHLT]/I" );
        tree->Branch( "HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I" );
        tree->Branch( "HLT_Name", &HLT_Name);

        if ( tt == ElectronMuon || tt == DiMuon ) {
          tree->Branch( "nMuonCand", &nMuonCand, "nMuonCand/i" );
          tree->Branch( "MuonCand_pt", MuonCand_pt, "MuonCand_pt[nMuonCand]/D" );
          tree->Branch( "MuonCand_eta", MuonCand_eta, "MuonCand_eta[nMuonCand]/D" );
          tree->Branch( "MuonCand_phi", MuonCand_phi, "MuonCand_phi[nMuonCand]/D" );
          tree->Branch( "MuonCand_e", MuonCand_e, "MuonCand_e[nMuonCand]/D" );
          tree->Branch( "MuonCand_charge", MuonCand_charge, "MuonCand_charge[nMuonCand]/I" );
          tree->Branch( "MuonCand_vtxx", MuonCand_vtxx, "MuonCand_vtxx[nMuonCand]/D" );
          tree->Branch( "MuonCand_vtxy", MuonCand_vtxy, "MuonCand_vtxy[nMuonCand]/D" );
          tree->Branch( "MuonCand_vtxz", MuonCand_vtxz, "MuonCand_vtxz[nMuonCand]/D" );
          tree->Branch( "MuonCand_dxy", MuonCand_dxy, "MuonCand_dxy[nMuonCand]/D" );
          tree->Branch( "MuonCand_nstatseg", MuonCand_nstatseg, "MuonCand_nstatseg[nMuonCand]/I" );
          tree->Branch( "MuonCand_ntrklayers", MuonCand_ntrklayers, "MuonCand_ntrklayers[nMuonCand]/I" );
          tree->Branch( "MuonCand_npxlhits", MuonCand_npxlhits, "MuonCand_npxlhits[nMuonCand]/I" );
          tree->Branch( "MuonCand_isglobal", MuonCand_isglobal, "MuonCand_isglobal[nMuonCand]/I" );
          tree->Branch( "MuonCand_istracker", MuonCand_istracker, "MuonCand_istracker[nMuonCand]/I" );
          tree->Branch( "MuonCand_isstandalone", MuonCand_isstandalone, "MuonCand_isstandalone[nMuonCand]/I" );
          tree->Branch( "MuonCand_ispfmuon", MuonCand_ispfmuon, "MuonCand_ispfmuon[nMuonCand]/I" );
          tree->Branch( "MuonCand_istight", MuonCand_istight, "MuonCand_istight[nMuonCand]/I" );
          tree->Branch( "MuonCandTrack_nmuchits", MuonCandTrack_nmuchits, "MuonCandTrack_nmuchits[nMuonCand]/I" );
          tree->Branch( "MuonCandTrack_chisq", MuonCandTrack_chisq, "MuonCandTrack_chisq[nMuonCand]/D" );
          tree->Branch( "MuonCand_innerTrackPt", MuonCand_innerTrackPt, "MuonCand_innerTrackPt[nMuonCand]/D" );
          tree->Branch( "MuonCand_innerTrackEta", MuonCand_innerTrackEta, "MuonCand_innerTrackEta[nMuonCand]/D" );
          tree->Branch( "MuonCand_innerTrackPhi", MuonCand_innerTrackPhi, "MuonCand_innerTrackPhi[nMuonCand]/D" );
          tree->Branch( "MuonCand_innerTrackVtxz", MuonCand_innerTrackVtxz, "MuonCand_innerTrackVtxz[nMuonCand]/D" );
          if ( mc ) {
            tree->Branch( "nGenMuonCand", &nGenMuonCand, "nGenMuonCand/i" );
            tree->Branch( "nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept, "nGenMuonCandOutOfAccept/i" );
            tree->Branch( "GenMuonCand_pt", GenMuonCand_pt, "GenMuonCand_pt[nGenMuonCand]/D" );
            tree->Branch( "GenMuonCand_eta", GenMuonCand_eta, "GenMuonCand_eta[nGenMuonCand]/D" );
            tree->Branch( "GenMuonCand_phi", GenMuonCand_phi, "GenMuonCand_phi[nGenMuonCand]/D" );
            tree->Branch( "GenMuonCand_e", GenMuonCand_e, "GenMuonCand_e[nGenMuonCand]/D" );
          }
        }

        if ( tt == ElectronMuon || tt == DiElectron ) {
          tree->Branch( "nEleCand", &nEleCand, "nEleCand/i" );
          tree->Branch( "EleCand_et", EleCand_et, "EleCand_et[nEleCand]/D" );
          tree->Branch( "EleCand_eta", EleCand_eta, "EleCand_eta[nEleCand]/D" );
          tree->Branch( "EleCand_phi", EleCand_phi, "EleCand_phi[nEleCand]/D" );
          tree->Branch( "EleCand_e", EleCand_e, "EleCand_e[nEleCand]/D" );
          tree->Branch( "EleCand_charge", EleCand_charge, "EleCand_charge[nEleCand]/I" );
          tree->Branch( "EleCand_vtxx", EleCand_vtxx, "EleCand_vtxx[nEleCand]/D" );
          tree->Branch( "EleCand_vtxy", EleCand_vtxy, "EleCand_vtxy[nEleCand]/D" );
          tree->Branch( "EleCand_vtxz", EleCand_vtxz, "EleCand_vtxz[nEleCand]/D" );
          tree->Branch( "EleCand_deltaPhi", EleCand_deltaPhi, "EleCand_deltaPhi[nEleCand]/D" );
          tree->Branch( "EleCand_deltaEta", EleCand_deltaEta, "EleCand_deltaEta[nEleCand]/D" );
          tree->Branch( "EleCand_HoverE", EleCand_HoverE, "EleCand_HoverE[nEleCand]/D" );
          tree->Branch( "EleCand_trackiso", EleCand_trackiso, "EleCand_trackiso[nEleCand]/D" );
          tree->Branch( "EleCand_ecaliso", EleCand_ecaliso," EleCand_ecaliso[nEleCand]/D" );
          tree->Branch( "EleCand_hcaliso", EleCand_hcaliso," EleCand_hcaliso[nEleCand]/D" );
          tree->Branch( "EleCand_sigmaIetaIeta", EleCand_sigmaIetaIeta, "EleCand_sigmaIetaIeta[nEleCand]/D" );
          tree->Branch( "EleCand_convDist", EleCand_convDist, "EleCand_convDist[nEleCand]/D" );
          tree->Branch( "EleCand_convDcot", EleCand_convDcot, "EleCand_convDcot[nEleCand]/D" );
          tree->Branch( "EleCand_ecalDriven", EleCand_ecalDriven, "EleCand_ecalDriven[nEleCand]/D" );
          tree->Branch( "EleCand_mediumID", EleCand_mediumID, "EleCand_mediumID[nEleCand]/I" );
          tree->Branch( "EleCand_tightID", EleCand_tightID, "EleCand_tightID[nEleCand]/I" );
          tree->Branch( "EleCand_innerTrackPt", EleCand_innerTrackPt, "EleCand_innerTrackPt[nEleCand]/D" );
          tree->Branch( "EleCand_innerTrackEta", EleCand_innerTrackEta, "EleCand_innerTrackEta[nEleCand]/D" );
          tree->Branch( "EleCand_innerTrackPhi", EleCand_innerTrackPhi, "EleCand_innerTrackPhi[nEleCand]/D" );
          tree->Branch( "EleCand_innerTrackVtxz", EleCand_innerTrackVtxz, "EleCand_innerTrackVtxz[nEleCand]/D" );
          if ( mc ) {
            tree->Branch( "nGenEleCand", &nGenEleCand, "nGenEleCand/i" );
            tree->Branch( "nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept, "nGenEleCandOutOfAccept/I" );
            tree->Branch( "GenEleCand_pt", GenEleCand_pt, "GenEleCand_pt[nGenEleCand]/D" );
            tree->Branch( "GenEleCand_eta", GenEleCand_eta, "GenEleCand_eta[nGenEleCand]/D" );
            tree->Branch( "GenEleCand_phi", GenEleCand_phi, "GenEleCand_phi[nGenEleCand]/D" );
            tree->Branch( "GenEleCand_e", GenEleCand_e, "GenEleCand_e[nGenEleCand]/D" );
          }
        }
        tree->Branch( "nPhotonCand", &nPhotonCand, "nPhotonCand/i" );
        tree->Branch( "PhotonCand_pt", PhotonCand_pt, "PhotonCand_pt[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_eta", PhotonCand_eta, "PhotonCand_eta[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_phi", PhotonCand_phi, "PhotonCand_phi[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_e", PhotonCand_e, "PhotonCand_e[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_r9", PhotonCand_r9, "PhotonCand_r9[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_drtrue", PhotonCand_drtrue, "PhotonCand_drtrue[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_detatrue", PhotonCand_detatrue, "PhotonCand_detatrue[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_dphitrue", PhotonCand_dphitrue, "PhotonCand_dphitrue[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_mediumID", PhotonCand_mediumID, "PhotonCand_mediumID[nPhotonCand]/I" );
        tree->Branch( "PhotonCand_tightID", PhotonCand_tightID, "PhotonCand_tightID[nPhotonCand]/I" );
        if ( mc ) {
          tree->Branch( "nGenPhotCand", &nGenPhotCand, "nGenPhotCand/i" );
          tree->Branch( "nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept, "nGenPhotCandOutOfAccept/I" );
          tree->Branch( "GenPhotCand_pt", GenPhotCand_pt, "GenPhotCand_pt[nGenPhotCand]/D" );
          tree->Branch( "GenPhotCand_eta", GenPhotCand_eta, "GenPhotCand_eta[nGenPhotCand]/D" );
          tree->Branch( "GenPhotCand_phi", GenPhotCand_phi, "GenPhotCand_phi[nGenPhotCand]/D" );
          tree->Branch( "GenPhotCand_e", GenPhotCand_e, "GenPhotCand_e[nGenPhotCand]/D" );

          tree->Branch( "nGenProtCand", &nGenProtCand, "nGenProtCand/i" );
          tree->Branch( "GenProtCand_pt", GenProtCand_pt, "GenProtCand_pt[nGenProtCand]/D" );
          tree->Branch( "GenProtCand_eta", GenProtCand_eta, "GenProtCand_eta[nGenProtCand]/D" );
          tree->Branch( "GenProtCand_phi", GenProtCand_phi, "GenProtCand_phi[nGenProtCand]/D" );
          tree->Branch( "GenProtCand_e", GenProtCand_e, "GenProtCand_e[nGenProtCand]/D" );
          tree->Branch( "GenProtCand_status", GenProtCand_status, "GenProtCand_status[nGenProtCand]/I" );
        }

        // Primary vertices' information
        tree->Branch( "nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/i" );
        tree->Branch( "nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nFilteredPrimVertexCand/i" );
        tree->Branch( "PrimVertexCand_id", PrimVertexCand_id, "PrimVertexCand_id[nPrimVertexCand]/I" );
        tree->Branch( "PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/i" );
        tree->Branch( "PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/i" );

        // Lepton pairs' information
        tree->Branch( "nPair", &nPair, "nPair/i" );
        tree->Branch( "Pair_lepton1", Pair_lepton1, "Pair_lepton1[nPair]/I" );
        tree->Branch( "Pair_lepton2", Pair_lepton2, "Pair_lepton2[nPair]/I" );
        //tree->Branch( "Pair_mindist", Pair_mindist, "Pair_mindist[nPair]/D" );
        tree->Branch( "Pair_mass", Pair_mass, "Pair_mass[nPair]/D" );
        tree->Branch( "Pair_pt", Pair_pt, "Pair_pt[nPair]/D" );
        tree->Branch( "Pair_eta", Pair_eta, "Pair_eta[nPair]/D" );
        tree->Branch( "Pair_phi", Pair_phi, "Pair_phi[nPair]/D" );
        tree->Branch( "Pair_dpt", Pair_dpt, "Pair_dpt[nPair]/D" );
        tree->Branch( "Pair_dphi", Pair_dphi, "Pair_dphi[nPair]/D" );
        tree->Branch( "Pair_3Dangle", Pair_3Dangle, "Pair_3Dangle[nPair]/D" );
        tree->Branch( "Pair_extratracks0p5mm", Pair_extratracks0p5mm, "Pair_extratracks0p5mm[nPair]/i" );
        tree->Branch( "Pair_extratracks1mm", Pair_extratracks1mm, "Pair_extratracks1mm[nPair]/i" );
        tree->Branch( "Pair_extratracks2mm", Pair_extratracks2mm, "Pair_extratracks2mm[nPair]/i" );
        tree->Branch( "Pair_extratracks3mm", Pair_extratracks3mm, "Pair_extratracks3mm[nPair]/i" );
        tree->Branch( "Pair_extratracks4mm", Pair_extratracks4mm, "Pair_extratracks4mm[nPair]/i" );
        tree->Branch( "Pair_extratracks5mm", Pair_extratracks5mm, "Pair_extratracks5mm[nPair]/i" );
        tree->Branch( "Pair_extratracks1cm", Pair_extratracks1cm, "Pair_extratracks1cm[nPair]/i" );
        tree->Branch( "Pair_extratracks2cm", Pair_extratracks2cm, "Pair_extratracks2cm[nPair]/i" );
        tree->Branch( "Pair_extratracks3cm", Pair_extratracks3cm, "Pair_extratracks3cm[nPair]/i" );
        tree->Branch( "Pair_extratracks4cm", Pair_extratracks4cm, "Pair_extratracks4cm[nPair]/i" );
        tree->Branch( "Pair_extratracks5cm", Pair_extratracks5cm, "Pair_extratracks5cm[nPair]/i" );
        tree->Branch( "Pair_extratracks10cm", Pair_extratracks10cm, "Pair_extratracks10cm[nPair]/i" );
        // Kalman dilepton vertex information
        tree->Branch( "KalmanVertexCand_x", KalmanVertexCand_x, "KalmanVertexCand_x[nPair]/D" );
        tree->Branch( "KalmanVertexCand_y", KalmanVertexCand_y, "KalmanVertexCand_y[nPair]/D" );
        tree->Branch( "KalmanVertexCand_z", KalmanVertexCand_z, "KalmanVertexCand_z[nPair]/D" );

        tree->Branch( "nPairGamma", &nPairGamma, "nPairGamma/i" );
        tree->Branch( "PairGamma_pair", PairGamma_pair, "PairGamma_pair[nPairGamma]/I" );
        tree->Branch( "PairGamma_mass", PairGamma_mass, "PairGamma_mass[nPairGamma]/D" );
        if ( mc ) {
          tree->Branch( "GenPair_mass", &GenPair_mass, "GenPair_mass/D" );
          tree->Branch( "GenPair_pt", &GenPair_pt, "GenPair_pt/D" );
          tree->Branch( "GenPair_eta", &GenPair_eta, "GenPair_eta/D" );
          tree->Branch( "GenPair_phi", &GenPair_phi, "GenPair_phi/D" );
          tree->Branch( "GenPair_dpt", &GenPair_dpt, "GenPair_dpt/D" );
          tree->Branch( "GenPair_dphi", &GenPair_dphi, "GenPair_dphi/D" );
          tree->Branch( "GenPair_3Dangle", &GenPair_3Dangle, "GenPair_3Dangle/D" );
        }

        if ( !mc ) {
          tree->Branch( "nLocalProtCand", &nLocalProtCand, "nLocalProtCand/i" );
          tree->Branch( "LocalProtCand_x", LocalProtCand_x, "LocalProtCand_x[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_y", LocalProtCand_y, "LocalProtCand_y[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_z", LocalProtCand_z, "LocalProtCand_z[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_xSigma", LocalProtCand_xSigma, "LocalProtCand_xSigma[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_ySigma", LocalProtCand_ySigma, "LocalProtCand_ySigma[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_arm", LocalProtCand_arm, "LocalProtCand_arm[nLocalProtCand]/I" );
          tree->Branch( "LocalProtCand_pot", LocalProtCand_pot, "LocalProtCand_pot[nLocalProtCand]/I" );
          tree->Branch( "LocalProtCand_Tx", LocalProtCand_Tx, "LocalProtCand_Tx[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_Ty", LocalProtCand_Ty, "LocalProtCand_Ty[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_TxSigma", LocalProtCand_TxSigma, "LocalProtCand_TxSigma[nLocalProtCand]/D" );
          tree->Branch( "LocalProtCand_TySigma", LocalProtCand_TySigma, "LocalProtCand_TySigma[nLocalProtCand]/D" );
        }

        // Extra tracks on vertex's information
        tree->Branch( "nExtraTracks", &nExtraTracks, "nExtraTracks/i" );
        tree->Branch( "ExtraTrack_pair", ExtraTrack_pair, "ExtraTrack_pair[nExtraTracks]/I" );
        tree->Branch( "ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I" );
        tree->Branch( "ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/i" );
        tree->Branch( "ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I" );
        tree->Branch( "ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/i" );
        tree->Branch( "ExtraTrack_px", ExtraTrack_px, "ExtraTrack_px[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_py", ExtraTrack_py, "ExtraTrack_py[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_pz", ExtraTrack_pz, "ExtraTrack_pz[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_chi2", ExtraTrack_chi2, "ExtraTrack_chi2[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz, "ExtraTrack_vtxdxyz[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_vtxT", ExtraTrack_vtxT, "ExtraTrack_vtxT[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_vtxZ", ExtraTrack_vtxZ, "ExtraTrack_vtxZ[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_x", ExtraTrack_x, "ExtraTrack_x[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_y", ExtraTrack_y, "ExtraTrack_y[nExtraTracks]/D" );
        tree->Branch( "ExtraTrack_z", ExtraTrack_z, "ExtraTrack_z[nExtraTracks]/D" );
        tree->Branch( "nQualityExtraTrack", &nQualityExtraTrack, "nQualityExtraTrack/i" );
        tree->Branch( "ClosestExtraTrack_vtxdxyz", ClosestExtraTrack_vtxdxyz, "ClosestExtraTrack_vtxdxyz[nPair]/D" );
        tree->Branch( "ClosestExtraTrack_id", ClosestExtraTrack_id, "ClosestExtraTrack_id[nPair]/I" );
        tree->Branch( "ClosestHighPurityExtraTrack_vtxdxyz", ClosestHighPurityExtraTrack_vtxdxyz, "ClosestHighPurityExtraTrack_vtxdxyz[nPair]/D" );
        tree->Branch( "ClosestHighPurityExtraTrack_id", ClosestHighPurityExtraTrack_id, "ClosestHighPurityExtraTrack_id[nPair]/I" );

        // Jets/MET information
        tree->Branch( "nJetCand", &nJetCand, "nJetCand/i" );
        tree->Branch( "JetCand_pt", JetCand_pt, "JetCand_pt[nJetCand]/D" );
        tree->Branch( "JetCand_eta", JetCand_eta, "JetCand_eta[nJetCand]/D" );
        tree->Branch( "JetCand_phi", JetCand_phi, "JetCand_phi[nJetCand]/D" );
        tree->Branch( "JetCand_e", JetCand_e, "JetCand_e[nJetCand]/D" );
        tree->Branch( "HighestJet_pt", &HighestJet_pt, "HighestJet_pt/D" );
        tree->Branch( "HighestJet_eta", &HighestJet_eta, "HighestJet_eta/D" );
        tree->Branch( "HighestJet_phi", &HighestJet_phi, "HighestJet_phi/D" );
        tree->Branch( "HighestJet_e", &HighestJet_e, "HighestJet_e/D" );
        tree->Branch( "SumJet_e", &SumJet_e, "SumJet_e/D" );
        tree->Branch( "Etmiss", &Etmiss, "Etmiss/D" );
        tree->Branch( "Etmiss_phi", &Etmiss_phi, "Etmiss_phi/D" );
        tree->Branch( "Etmiss_significance", &Etmiss_significance, "Etmiss_significance/D" );

        // Pileup reweighting
        tree->Branch( "Weight", &Weight, "Weight/D" );
        tree->Branch( "PUWeightTrue", &PUWeightTrue, "PUWeightTrue/D" );
      }
      void load( TTree* tree, TreeType tt, bool mc ) {
        if ( !tree ) return;

        tree->SetBranchAddress( "Run", &Run );
        tree->SetBranchAddress( "LumiSection", &LumiSection );
        tree->SetBranchAddress( "BX", &BX );
        tree->SetBranchAddress( "EventNum", &EventNum );

        tree->SetBranchAddress( "nHLT", &nHLT );
        tree->SetBranchAddress( "HLT_Accept", HLT_Accept );
        tree->SetBranchAddress( "HLT_Prescl", HLT_Prescl );
        tree->SetBranchAddress( "HLT_Name", &HLT_Name );

        if ( tt == ElectronMuon || tt == DiMuon ) {
          tree->SetBranchAddress( "nMuonCand", &nMuonCand );
          tree->SetBranchAddress( "MuonCand_pt", MuonCand_pt );
          tree->SetBranchAddress( "MuonCand_eta", MuonCand_eta );
          tree->SetBranchAddress( "MuonCand_phi", MuonCand_phi );
          tree->SetBranchAddress( "MuonCand_e", MuonCand_e );
          tree->SetBranchAddress( "MuonCand_charge", MuonCand_charge );
          tree->SetBranchAddress( "MuonCand_vtxx", MuonCand_vtxx );
          tree->SetBranchAddress( "MuonCand_vtxy", MuonCand_vtxy );
          tree->SetBranchAddress( "MuonCand_vtxz", MuonCand_vtxz );
          tree->SetBranchAddress( "MuonCand_dxy", MuonCand_dxy );
          tree->SetBranchAddress( "MuonCand_nstatseg", MuonCand_nstatseg );
          tree->SetBranchAddress( "MuonCand_ntrklayers", MuonCand_ntrklayers );
          tree->SetBranchAddress( "MuonCand_npxlhits", MuonCand_npxlhits );
          tree->SetBranchAddress( "MuonCand_isglobal", MuonCand_isglobal );
          tree->SetBranchAddress( "MuonCand_istracker", MuonCand_istracker );
          tree->SetBranchAddress( "MuonCand_isstandalone", MuonCand_isstandalone );
          tree->SetBranchAddress( "MuonCand_ispfmuon", MuonCand_ispfmuon );
          tree->SetBranchAddress( "MuonCand_istight", MuonCand_istight );
          tree->SetBranchAddress( "MuonCandTrack_nmuchits", MuonCandTrack_nmuchits );
          tree->SetBranchAddress( "MuonCandTrack_chisq", MuonCandTrack_chisq );
          tree->SetBranchAddress( "MuonCand_innerTrackPt", MuonCand_innerTrackPt );
          tree->SetBranchAddress( "MuonCand_innerTrackEta", MuonCand_innerTrackEta );
          tree->SetBranchAddress( "MuonCand_innerTrackPhi", MuonCand_innerTrackPhi );
          tree->SetBranchAddress( "MuonCand_innerTrackVtxz", MuonCand_innerTrackVtxz );
          if ( mc ) {
            tree->SetBranchAddress( "nGenMuonCand", &nGenMuonCand );
            tree->SetBranchAddress( "nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept );
            tree->SetBranchAddress( "GenMuonCand_pt", GenMuonCand_pt );
            tree->SetBranchAddress( "GenMuonCand_eta", GenMuonCand_eta );
            tree->SetBranchAddress( "GenMuonCand_phi", GenMuonCand_phi );
            tree->SetBranchAddress( "GenMuonCand_e", GenMuonCand_e );
          }
        }

        if ( tt == ElectronMuon || tt == DiElectron ) {
          tree->SetBranchAddress( "nEleCand", &nEleCand );
          tree->SetBranchAddress( "EleCand_et", EleCand_et );
          tree->SetBranchAddress( "EleCand_eta", EleCand_eta );
          tree->SetBranchAddress( "EleCand_phi", EleCand_phi );
          tree->SetBranchAddress( "EleCand_e", EleCand_e );
          tree->SetBranchAddress( "EleCand_charge", EleCand_charge );
          tree->SetBranchAddress( "EleCand_vtxx", EleCand_vtxx );
          tree->SetBranchAddress( "EleCand_vtxy", EleCand_vtxy );
          tree->SetBranchAddress( "EleCand_vtxz", EleCand_vtxz );
          tree->SetBranchAddress( "EleCand_deltaPhi", EleCand_deltaPhi );
          tree->SetBranchAddress( "EleCand_deltaEta", EleCand_deltaEta );
          tree->SetBranchAddress( "EleCand_HoverE", EleCand_HoverE );
          tree->SetBranchAddress( "EleCand_trackiso", EleCand_trackiso );
          tree->SetBranchAddress( "EleCand_ecaliso", EleCand_ecaliso );
          tree->SetBranchAddress( "EleCand_hcaliso", EleCand_hcaliso );
          tree->SetBranchAddress( "EleCand_sigmaIetaIeta", EleCand_sigmaIetaIeta );
          tree->SetBranchAddress( "EleCand_convDist", EleCand_convDist );
          tree->SetBranchAddress( "EleCand_convDcot", EleCand_convDcot );
          tree->SetBranchAddress( "EleCand_ecalDriven", EleCand_ecalDriven );
          tree->SetBranchAddress( "EleCand_mediumID", EleCand_mediumID );
          tree->SetBranchAddress( "EleCand_tightID", EleCand_tightID );
          tree->SetBranchAddress( "EleCand_innerTrackPt", EleCand_innerTrackPt );
          tree->SetBranchAddress( "EleCand_innerTrackEta", EleCand_innerTrackEta );
          tree->SetBranchAddress( "EleCand_innerTrackPhi", EleCand_innerTrackPhi );
          tree->SetBranchAddress( "EleCand_innerTrackVtxz", EleCand_innerTrackVtxz );
          if ( mc ) {
            tree->SetBranchAddress( "nGenEleCand", &nGenEleCand );
            tree->SetBranchAddress( "nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept );
            tree->SetBranchAddress( "GenEleCand_pt", GenEleCand_pt );
            tree->SetBranchAddress( "GenEleCand_eta", GenEleCand_eta );
            tree->SetBranchAddress( "GenEleCand_phi", GenEleCand_phi );
            tree->SetBranchAddress( "GenEleCand_e", GenEleCand_e );
          }
        }
        tree->SetBranchAddress( "nPhotonCand", &nPhotonCand );
        tree->SetBranchAddress( "PhotonCand_pt", PhotonCand_pt );
        tree->SetBranchAddress( "PhotonCand_eta", PhotonCand_eta );
        tree->SetBranchAddress( "PhotonCand_phi", PhotonCand_phi );
        tree->SetBranchAddress( "PhotonCand_e", PhotonCand_e );
        tree->SetBranchAddress( "PhotonCand_r9", PhotonCand_r9 );
        tree->SetBranchAddress( "PhotonCand_drtrue", PhotonCand_drtrue );
        tree->SetBranchAddress( "PhotonCand_detatrue", PhotonCand_detatrue );
        tree->SetBranchAddress( "PhotonCand_dphitrue", PhotonCand_dphitrue );
        tree->SetBranchAddress( "PhotonCand_mediumID", PhotonCand_mediumID );
        tree->SetBranchAddress( "PhotonCand_tightID", PhotonCand_tightID );
        if ( mc ) {
          tree->SetBranchAddress( "nGenPhotCand", &nGenPhotCand );
          tree->SetBranchAddress( "nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept );
          tree->SetBranchAddress( "GenPhotCand_pt", GenPhotCand_pt );
          tree->SetBranchAddress( "GenPhotCand_eta", GenPhotCand_eta );
          tree->SetBranchAddress( "GenPhotCand_phi", GenPhotCand_phi );
          tree->SetBranchAddress( "GenPhotCand_e", GenPhotCand_e );

          tree->SetBranchAddress( "nGenProtCand", &nGenProtCand );
          tree->SetBranchAddress( "GenProtCand_pt", GenProtCand_pt );
          tree->SetBranchAddress( "GenProtCand_eta", GenProtCand_eta );
          tree->SetBranchAddress( "GenProtCand_phi", GenProtCand_phi );
          tree->SetBranchAddress( "GenProtCand_e", GenProtCand_e );
          tree->SetBranchAddress( "GenProtCand_status", GenProtCand_status );
        }

        // Primary vertices' information
        tree->SetBranchAddress( "nPrimVertexCand", &nPrimVertexCand );
        tree->SetBranchAddress( "nFilteredPrimVertexCand", &nFilteredPrimVertexCand );
        tree->SetBranchAddress( "PrimVertexCand_id", PrimVertexCand_id );
        tree->SetBranchAddress( "PrimVertexCand_x", PrimVertexCand_x );
        tree->SetBranchAddress( "PrimVertexCand_y", PrimVertexCand_y );
        tree->SetBranchAddress( "PrimVertexCand_z", PrimVertexCand_z );
        tree->SetBranchAddress( "PrimVertexCand_chi2", PrimVertexCand_chi2 );
        tree->SetBranchAddress( "PrimVertexCand_ndof", PrimVertexCand_ndof );
        tree->SetBranchAddress( "PrimVertexCand_tracks", PrimVertexCand_tracks );

        // Lepton pairs' information
        tree->SetBranchAddress( "nPair", &nPair );
        tree->SetBranchAddress( "Pair_lepton1", Pair_lepton1 );
        tree->SetBranchAddress( "Pair_lepton2", Pair_lepton2 );
        //tree->SetBranchAddress( "Pair_mindist", Pair_mindist );
        tree->SetBranchAddress( "Pair_mass", Pair_mass );
        tree->SetBranchAddress( "Pair_pt", Pair_pt );
        tree->SetBranchAddress( "Pair_eta", Pair_eta );
        tree->SetBranchAddress( "Pair_phi", Pair_phi );
        tree->SetBranchAddress( "Pair_dpt", Pair_dpt );
        tree->SetBranchAddress( "Pair_dphi", Pair_dphi );
        tree->SetBranchAddress( "Pair_3Dangle", Pair_3Dangle );
        tree->SetBranchAddress( "Pair_extratracks0p5mm", Pair_extratracks0p5mm );
        tree->SetBranchAddress( "Pair_extratracks1mm", Pair_extratracks1mm );
        tree->SetBranchAddress( "Pair_extratracks2mm", Pair_extratracks2mm );
        tree->SetBranchAddress( "Pair_extratracks3mm", Pair_extratracks3mm );
        tree->SetBranchAddress( "Pair_extratracks4mm", Pair_extratracks4mm );
        tree->SetBranchAddress( "Pair_extratracks5mm", Pair_extratracks5mm );
        tree->SetBranchAddress( "Pair_extratracks1cm", Pair_extratracks1cm );
        tree->SetBranchAddress( "Pair_extratracks2cm", Pair_extratracks2cm );
        tree->SetBranchAddress( "Pair_extratracks3cm", Pair_extratracks3cm );
        tree->SetBranchAddress( "Pair_extratracks4cm", Pair_extratracks4cm );
        tree->SetBranchAddress( "Pair_extratracks5cm", Pair_extratracks5cm );
        tree->SetBranchAddress( "Pair_extratracks10cm", Pair_extratracks10cm );
        // Kalman dilepton vertex information
        tree->SetBranchAddress( "KalmanVertexCand_x", KalmanVertexCand_x );
        tree->SetBranchAddress( "KalmanVertexCand_y", KalmanVertexCand_y );
        tree->SetBranchAddress( "KalmanVertexCand_z", KalmanVertexCand_z );

        tree->SetBranchAddress( "nPairGamma", &nPairGamma );
        tree->SetBranchAddress( "PairGamma_pair", PairGamma_pair );
        tree->SetBranchAddress( "PairGamma_mass", PairGamma_mass );
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
          tree->SetBranchAddress( "nLocalProtCand", &nLocalProtCand );
          tree->SetBranchAddress( "LocalProtCand_x", LocalProtCand_x );
          tree->SetBranchAddress( "LocalProtCand_y", LocalProtCand_y );
          tree->SetBranchAddress( "LocalProtCand_z", LocalProtCand_z );
          tree->SetBranchAddress( "LocalProtCand_xSigma", LocalProtCand_xSigma );
          tree->SetBranchAddress( "LocalProtCand_ySigma", LocalProtCand_ySigma );
          tree->SetBranchAddress( "LocalProtCand_arm", LocalProtCand_arm );
          tree->SetBranchAddress( "LocalProtCand_pot", LocalProtCand_pot );
          tree->SetBranchAddress( "LocalProtCand_Tx", LocalProtCand_Tx );
          tree->SetBranchAddress( "LocalProtCand_Ty", LocalProtCand_Ty );
          tree->SetBranchAddress( "LocalProtCand_TxSigma", LocalProtCand_TxSigma );
          tree->SetBranchAddress( "LocalProtCand_TySigma", LocalProtCand_TySigma );
        }

        // Extra tracks on vertex's information
        tree->SetBranchAddress( "nExtraTracks", &nExtraTracks );
        tree->SetBranchAddress( "ExtraTrack_pair", ExtraTrack_pair );
        tree->SetBranchAddress( "ExtraTrack_purity", ExtraTrack_purity );
        tree->SetBranchAddress( "ExtraTrack_nhits", ExtraTrack_nhits );
        tree->SetBranchAddress( "ExtraTrack_charge", ExtraTrack_charge );
        tree->SetBranchAddress( "ExtraTrack_ndof", ExtraTrack_ndof );
        tree->SetBranchAddress( "ExtraTrack_px", ExtraTrack_px );
        tree->SetBranchAddress( "ExtraTrack_py", ExtraTrack_py );
        tree->SetBranchAddress( "ExtraTrack_pz", ExtraTrack_pz );
        tree->SetBranchAddress( "ExtraTrack_chi2", ExtraTrack_chi2 );
        tree->SetBranchAddress( "ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz );
        tree->SetBranchAddress( "ExtraTrack_vtxT", ExtraTrack_vtxT );
        tree->SetBranchAddress( "ExtraTrack_vtxZ", ExtraTrack_vtxZ );
        tree->SetBranchAddress( "ExtraTrack_x", ExtraTrack_x );
        tree->SetBranchAddress( "ExtraTrack_y", ExtraTrack_y );
        tree->SetBranchAddress( "ExtraTrack_z", ExtraTrack_z );
        tree->SetBranchAddress( "nQualityExtraTrack", &nQualityExtraTrack );
        tree->SetBranchAddress( "ClosestExtraTrack_vtxdxyz", ClosestExtraTrack_vtxdxyz );
        tree->SetBranchAddress( "ClosestExtraTrack_id", ClosestExtraTrack_id );
        tree->SetBranchAddress( "ClosestHighPurityExtraTrack_vtxdxyz", ClosestHighPurityExtraTrack_vtxdxyz );
        tree->SetBranchAddress( "ClosestHighPurityExtraTrack_id", ClosestHighPurityExtraTrack_id );

        // Jets/MET information
        tree->SetBranchAddress( "nJetCand", &nJetCand );
        tree->SetBranchAddress( "JetCand_pt", JetCand_pt );
        tree->SetBranchAddress( "JetCand_eta", JetCand_eta );
        tree->SetBranchAddress( "JetCand_phi", JetCand_phi );
        tree->SetBranchAddress( "JetCand_e", JetCand_e );
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

