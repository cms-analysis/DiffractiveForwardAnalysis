#ifndef DiffractiveForwardAnalysis_GammaGammaLeptonLepton_AnalysisEvent_h
#define DiffractiveForwardAnalysis_GammaGammaLeptonLepton_AnalysisEvent_h

#include <vector>
#include <string>

class TTree;

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
      AnalysisEvent();
      void clear();
      void attach( TTree*, TreeType, bool );
      void retrieve( TTree*, TreeType, bool );

      //////  Leaves size  //////

      /// Maximum number of HLT to check
      static constexpr unsigned short MAX_HLT = 10;
      /// Maximum number of leptons per event
      static constexpr unsigned short MAX_LL = 50;
      /// Maximum number of muons per event
      static constexpr unsigned short MAX_MUONS = 25;
      /// Maximum number of electrons per event
      static constexpr unsigned short MAX_ELE = 25;
      /// Maximum number of photons per event
      static constexpr unsigned short MAX_PHO = 50;
      /// Maximum number of leptons pairs per event
      static constexpr unsigned short MAX_PAIRS = 25;
      static constexpr unsigned short MAX_PAIRPHO = 25;
      /// Maximum number of primary vertices per event
      static constexpr unsigned short MAX_VTX = 150;
      /// Maximum number of extra tracks per event
      static constexpr unsigned short MAX_ET = 1000;
      /// Maximum number of generator level muons per event
      static constexpr unsigned short MAX_GENMU = 25;
      /// Maximum number of generator level electrons per event
      static constexpr unsigned short MAX_GENELE = 25;
      /// Maximum number of generator level photons per event
      static constexpr unsigned short MAX_GENPHO = 10;
      /// Maximum number of generator level protons per event
      static constexpr unsigned short MAX_GENPRO = 8;
      /// Maximum number of jets per event
      static constexpr unsigned short MAX_JETS = 40;
      /// Maximum number of reconstructed local tracks in RPs
      static constexpr unsigned short MAX_LOCALPCAND = 25;
      /// Maximum number of reconstructed local tracks pairs in RPs
      static constexpr unsigned short MAX_LOCALPPAIRCAND = 10;

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
      double GenMuonCand_pt[MAX_GENMU], GenMuonCand_eta[MAX_GENMU], GenMuonCand_phi[MAX_GENMU], GenMuonCand_e[MAX_GENMU];
      int nGenEleCand, nGenEleCandOutOfAccept;
      double GenEleCand_pt[MAX_GENELE], GenEleCand_eta[MAX_GENELE], GenEleCand_phi[MAX_GENELE], GenEleCand_e[MAX_GENELE];
      double GenPair_pt, GenPair_eta, GenPair_phi, GenPair_mass;
      double GenPair_dphi, GenPair_dpt, GenPair_3Dangle;
      int nGenPhotCand, nGenPhotCandOutOfAccept;
      double GenPhotCand_pt[MAX_GENPHO], GenPhotCand_eta[MAX_GENPHO], GenPhotCand_phi[MAX_GENPHO], GenPhotCand_e[MAX_GENPHO];
      int nGenProtCand;
      double GenProtCand_pt[MAX_GENPRO], GenProtCand_eta[MAX_GENPRO], GenProtCand_phi[MAX_GENPRO], GenProtCand_e[MAX_GENPHO];
      int GenProtCand_status[MAX_GENPRO];

      int nLeptonCand, nCandidates, nCandidatesInEvent;

      // Pileup reweighting quantities
      double nTruePUafterPUWeight;
      double nTruePUafterPUWeightBXM1, nTruePUafterPUWeightBXP1, nTruePUafterPUWeightBX0;
      double PUWeightTrue;
      int nTruePUforPUWeight;
      int nTruePUforPUWeightBXM1, nTruePUforPUWeightBXP1, nTruePUforPUWeightBX0;
      double Weight;

      // Muon quantities
      int nMuonCand;
      double MuonCand_pt[MAX_LL], MuonCand_eta[MAX_LL], MuonCand_phi[MAX_LL], MuonCand_e[MAX_LL];
      double MuonCand_innerTrackPt[MAX_LL], MuonCand_innerTrackEta[MAX_LL], MuonCand_innerTrackPhi[MAX_LL];
      double MuonCand_innerTrackVtxz[MAX_LL];
      double MuonCand_vtxx[MAX_LL], MuonCand_vtxy[MAX_LL], MuonCand_vtxz[MAX_LL];
      int MuonCand_charge[MAX_LL];
      double MuonCand_dxy[MAX_LL];
      int MuonCand_nstatseg[MAX_LL], MuonCand_npxlhits[MAX_LL], MuonCand_ntrklayers[MAX_LL];
      double MuonCand_[MAX_LL];
      int MuonCandTrack_nmuchits[MAX_LL];
      double MuonCandTrack_chisq[MAX_LL];
      int MuonCand_isglobal[MAX_LL], MuonCand_istracker[MAX_LL], MuonCand_isstandalone[MAX_LL], MuonCand_ispfmuon[MAX_LL];
      int MuonCand_istight[MAX_LL];

      // Electron quantities
      int nEleCand;
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
      int EleCand_vetoID[MAX_LL], EleCand_tightID[MAX_LL], EleCand_mediumID[MAX_LL], EleCand_looseID[MAX_LL];

      // Photon quantities
      int nPhotonCand;
      double PhotonCand_pt[MAX_PHO], PhotonCand_eta[MAX_PHO], PhotonCand_phi[MAX_PHO], PhotonCand_e[MAX_PHO];
      double PhotonCand_r9[MAX_PHO];
      double PhotonCand_drtrue[MAX_PHO], PhotonCand_detatrue[MAX_PHO], PhotonCand_dphitrue[MAX_PHO];

      // Pair quantities
      int nPair;
      int Pair_lepton1[MAX_PAIRS], Pair_lepton2[MAX_PAIRS];
      double Pair_pt[MAX_PAIRS], Pair_eta[MAX_PAIRS], Pair_phi[MAX_PAIRS], Pair_mass[MAX_PAIRS];
      double Pair_dpt[MAX_PAIRS], Pair_dphi[MAX_PAIRS], Pair_3Dangle[MAX_PAIRS];
      double Pair_mindist[MAX_PAIRS];

      int nPairGamma;
      int PairGamma_pair[MAX_PHO];
      double PairGamma_mass[MAX_PHO];

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

      double KalmanVertexCand_x[MAX_VTX], KalmanVertexCand_y[MAX_VTX], KalmanVertexCand_z[MAX_VTX];
      double ClosestExtraTrackKalman_vtxdxyz[MAX_VTX];

      // Extra tracks on vertex quantities
      int nExtraTracks;
      int ExtraTrack_pair[MAX_ET];
      int ExtraTrack_purity[MAX_ET], ExtraTrack_nhits[MAX_ET];
      int ExtraTrack_charge[MAX_ET], ExtraTrack_ndof[MAX_ET];
      double ExtraTrack_px[MAX_ET], ExtraTrack_py[MAX_ET], ExtraTrack_pz[MAX_ET];
      double ExtraTrack_chi2[MAX_ET];
      double ExtraTrack_vtxdxyz[MAX_ET];
      double ExtraTrack_vtxT[MAX_ET], ExtraTrack_vtxZ[MAX_ET];
      double ExtraTrack_x[MAX_ET], ExtraTrack_y[MAX_ET], ExtraTrack_z[MAX_ET];
      double ClosestExtraTrack_vtxdxyz[MAX_VTX], ClosestHighPurityExtraTrack_vtxdxyz[MAX_VTX];
      int ClosestExtraTrack_id[MAX_VTX], ClosestHighPurityExtraTrack_id[MAX_VTX];
      int nQualityExtraTrack;

      // Jets/MET quantities
      int nJetCand;
      double JetCand_e[MAX_JETS], JetCand_pt[MAX_JETS], JetCand_eta[MAX_JETS], JetCand_phi[MAX_JETS];
      double HighestJet_pt, HighestJet_eta, HighestJet_phi, HighestJet_e;
      double SumJet_e;
      double Etmiss, Etmiss_phi, Etmiss_significance;

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
}

#endif
