#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "TTree.h"

namespace ggll
{
  AnalysisEvent::AnalysisEvent()
  {
    clear();
  }

  void
  AnalysisEvent::clear()
  {
    nCandidates = 0;
    nPair = 0;
    nPrimVertexCand = nFilteredPrimVertexCand = -1;
    nMuonCand = nEleCand = nLeptonCand = 0;
    nExtraTracks = nQualityExtraTrack = 0;
    nJetCand = 0;
    nGenMuonCand = nGenMuonCandOutOfAccept = 0;
    nGenEleCand = nGenEleCandOutOfAccept = 0;
    nGenPhotCand = nGenPhotCandOutOfAccept = 0;
    nGenProtCand = 0;
    nPhotonCand = nPairGamma = 0;

    //LHCFillNum = LHCBeamMode = -1;

    GenPair_pt = GenPair_mass = GenPair_phi = GenPair_eta = -999.;
    GenPair_dphi = GenPair_dpt = GenPair_3Dangle = 0.;

    for ( unsigned int i = 0; i < MAX_LL; ++i ) {
      MuonCand_pt[i] = MuonCand_eta[i] = MuonCand_phi[i] = MuonCand_e[i] = -999.;
      MuonCand_charge[i] = -999;
      MuonCand_vtxx[i] = MuonCand_vtxy[i] = MuonCand_vtxz[i] = -999.;
      MuonCand_innerTrackPt[i] = MuonCand_innerTrackEta[i] = MuonCand_innerTrackPhi[i] = -999.;
      MuonCand_innerTrackVtxz[i] = -999.;
      MuonCand_npxlhits[i] = MuonCand_nstatseg[i] = MuonCand_ntrklayers[i] = -999;
      MuonCand_dxy[i] = -999.;
      MuonCand_isglobal[i] = MuonCand_istracker[i] = MuonCand_isstandalone[i] = MuonCand_ispfmuon[i] = -999;
      MuonCand_istight[i] = -999;
      MuonCandTrack_nmuchits[i] = -999;
      MuonCandTrack_chisq[i] = -999.;
      EleCand_e[i] = EleCand_et[i] = EleCand_phi[i] = EleCand_eta[i] = -999.;
      EleCand_charge[i] = -999;
      EleCand_vtxx[i] = EleCand_vtxy[i] = EleCand_vtxz[i] = -999.;
      EleCand_innerTrackPt[i] = EleCand_innerTrackEta[i] = EleCand_innerTrackPhi[i] = -999.;
      EleCand_innerTrackVtxz[i] = -999.;
      EleCand_deltaPhi[i] = EleCand_deltaEta[i] = EleCand_HoverE[i] = -999.;
      EleCand_trackiso[i] = EleCand_ecaliso[i] = EleCand_hcaliso[i] = EleCand_sigmaIetaIeta[i] = -999.;
      EleCand_convDist[i] = EleCand_convDcot[i] = EleCand_ecalDriven[i] = -999.;
      EleCand_vetoID[i] = EleCand_tightID[i] = EleCand_mediumID[i] = EleCand_looseID[i] = -1;
    }
    for ( unsigned int i = 0; i < MAX_ET; ++i ) {
      ExtraTrack_px[i] = ExtraTrack_py[i] = ExtraTrack_pz[i] = -999.;
      ExtraTrack_charge[i] = ExtraTrack_ndof[i] = -999;
      ExtraTrack_chi2[i] = ExtraTrack_vtxdxyz[i] = -999.;
      ExtraTrack_vtxT[i] = ExtraTrack_vtxZ[i] = -999.;
      ExtraTrack_x[i] = ExtraTrack_y[i] = ExtraTrack_z[i] = -999.;
      ExtraTrack_pair[i] = -1;
    }
    for ( unsigned int i = 0; i < MAX_PAIRS; ++i ) {
      Pair_lepton1[i] = Pair_lepton2[i] = -1;
      Pair_mindist[i] = Pair_pt[i] = Pair_mass[i] = Pair_phi[i] = Pair_eta[i] = -999.;
      Pair_dphi[i] = Pair_dpt[i] = Pair_3Dangle[i] = -999.;
      Pair_extratracks1mm[i] = Pair_extratracks2mm[i] = Pair_extratracks3mm[i] = 0;
      Pair_extratracks4mm[i] = Pair_extratracks5mm[i] = Pair_extratracks1cm[i] = 0;
      Pair_extratracks2cm[i] = Pair_extratracks3cm[i] = Pair_extratracks4cm[i] = 0;
      Pair_extratracks5cm[i] = Pair_extratracks10cm[i] = 0;    
    }
    for ( unsigned int i = 0; i < MAX_VTX; ++i ) {
      PrimVertexCand_id[i] = -1;
      PrimVertexCand_tracks[i] = PrimVertexCand_matchedtracks[i] = PrimVertexCand_unmatchedtracks[i] = PrimVertexCand_hasdil[i] = 0;
      PrimVertexCand_x[i] = PrimVertexCand_y[i] = PrimVertexCand_z[i] = -999.;
      PrimVertexCand_chi2[i] = PrimVertexCand_ndof[i] = -999.;
      KalmanVertexCand_x[i] = KalmanVertexCand_y[i] = KalmanVertexCand_z[i] = -999.;
      ClosestExtraTrackKalman_vtxdxyz[i] = 999.;
    }
    for ( unsigned int i = 0; i < MAX_PHO; ++i ) { 
      PhotonCand_e[i] = PhotonCand_pt[i] = PhotonCand_eta[i] = PhotonCand_phi[i] = PhotonCand_r9[i] = -999.;
      PhotonCand_detatrue[i] = PhotonCand_dphitrue[i] = PhotonCand_drtrue[i] = -999.;
    }
    for ( unsigned int i = 0; i < MAX_PAIRPHO; ++i ) {
      PairGamma_pair[i] = -1;
      PairGamma_mass[i] = -999.;
    }
    for ( unsigned int i = 0; i < MAX_JETS; ++i ) {
      JetCand_pt[i] = JetCand_phi[i] = JetCand_eta[i] = JetCand_e[i] = -999;
    }
    HighestJet_pt = HighestJet_eta = HighestJet_phi = HighestJet_e = -999.;
    SumJet_e = 0.;
    Etmiss = Etmiss_phi = Etmiss_significance = -999.;
    for ( unsigned int i = 0; i < MAX_LOCALPCAND; ++i ) {
      LocalProtCand_x[i] = LocalProtCand_y[i] = LocalProtCand_z[i] = -999.;
      LocalProtCand_xSigma[i] = LocalProtCand_ySigma[i] = -999.;
      LocalProtCand_Tx[i] = LocalProtCand_Ty[i] = -999.;
      LocalProtCand_TxSigma[i] = LocalProtCand_TySigma[i] = -999.;
      LocalProtCand_arm[i] = LocalProtCand_side[i] = -1;
    }
  }

  void
  AnalysisEvent::attach( TTree* tree, TreeType tt, bool mc )
  {
    tree->Branch( "Run", &Run, "Run/I" );
    tree->Branch( "LumiSection", &LumiSection, "LumiSection/I" );
    tree->Branch( "BX", &BX, "BX/I" );
    tree->Branch( "EventNum", &EventNum, "EventNum/I" );

    tree->Branch( "nHLT", &nHLT, "nHLT/I" );
    tree->Branch( "HLT_Accept", HLT_Accept, "HLT_Accept[nHLT]/I" );
    tree->Branch( "HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I" );
    tree->Branch( "HLT_Name", &HLT_Name);

    if ( tt == ElectronMuon || tt == DiMuon ) {
      tree->Branch( "nMuonCand", &nMuonCand, "nMuonCand/I" );
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
        tree->Branch( "nGenMuonCand", &nGenMuonCand, "nGenMuonCand/I" );
        tree->Branch( "nGenMuonCandOutOfAccept", &nGenMuonCandOutOfAccept, "nGenMuonCandOutOfAccept/I" );
        tree->Branch( "GenMuonCand_pt", GenMuonCand_pt, "GenMuonCand_pt[nGenMuonCand]/D" );
        tree->Branch( "GenMuonCand_eta", GenMuonCand_eta, "GenMuonCand_eta[nGenMuonCand]/D" );
        tree->Branch( "GenMuonCand_phi", GenMuonCand_phi, "GenMuonCand_phi[nGenMuonCand]/D" );
        tree->Branch( "GenMuonCand_e", GenMuonCand_e, "GenMuonCand_e[nGenMuonCand]/D" );
      }
    }

    if ( tt == ElectronMuon || tt == DiElectron ) {
      tree->Branch( "nEleCand", &nEleCand, "nEleCand/I" );
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
      tree->Branch( "EleCand_looseID", EleCand_looseID, "EleCand_looseID[nEleCand]/I" );
      tree->Branch( "EleCand_tightID", EleCand_tightID, "EleCand_tightID[nEleCand]/I" );
      tree->Branch( "EleCand_vetoID", EleCand_vetoID, "EleCand_vetoID[nEleCand]/I" );
      tree->Branch( "EleCand_innerTrackPt", EleCand_innerTrackPt, "EleCand_innerTrackPt[nEleCand]/D" );
      tree->Branch( "EleCand_innerTrackEta", EleCand_innerTrackEta, "EleCand_innerTrackEta[nEleCand]/D" );
      tree->Branch( "EleCand_innerTrackPhi", EleCand_innerTrackPhi, "EleCand_innerTrackPhi[nEleCand]/D" );
      tree->Branch( "EleCand_innerTrackVtxz", EleCand_innerTrackVtxz, "EleCand_innerTrackVtxz[nEleCand]/D" );
      if ( mc ) {
        tree->Branch( "nGenEleCand", &nGenEleCand, "nGenEleCand/I" );
        tree->Branch( "nGenEleCandOutOfAccept", &nGenEleCandOutOfAccept, "nGenEleCandOutOfAccept/I" );
        tree->Branch( "GenEleCand_pt", GenEleCand_pt, "GenEleCand_pt[nGenEleCand]/D" );
        tree->Branch( "GenEleCand_eta", GenEleCand_eta, "GenEleCand_eta[nGenEleCand]/D" );
        tree->Branch( "GenEleCand_phi", GenEleCand_phi, "GenEleCand_phi[nGenEleCand]/D" );
        tree->Branch( "GenEleCand_e", GenEleCand_e, "GenEleCand_e[nGenEleCand]/D" );
      }
    }
    tree->Branch( "nPhotonCand", &nPhotonCand, "nPhotonCand/I" );
    tree->Branch( "PhotonCand_pt", PhotonCand_pt, "PhotonCand_pt[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_eta", PhotonCand_eta, "PhotonCand_eta[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_phi", PhotonCand_phi, "PhotonCand_phi[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_e", PhotonCand_e, "PhotonCand_e[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_r9", PhotonCand_r9, "PhotonCand_r9[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_drtrue", PhotonCand_drtrue, "PhotonCand_drtrue[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_detatrue", PhotonCand_detatrue, "PhotonCand_detatrue[nPhotonCand]/D" );
    tree->Branch( "PhotonCand_dphitrue", PhotonCand_dphitrue, "PhotonCand_dphitrue[nPhotonCand]/D" );
    if ( mc ) {
      tree->Branch( "nGenPhotCand", &nGenPhotCand, "nGenPhotCand/I" );
      tree->Branch( "nGenPhotCandOutOfAccept", &nGenPhotCandOutOfAccept, "nGenPhotCandOutOfAccept/I" );
      tree->Branch( "GenPhotCand_pt", GenPhotCand_pt, "GenPhotCand_pt[nGenPhotCand]/D" );
      tree->Branch( "GenPhotCand_eta", GenPhotCand_eta, "GenPhotCand_eta[nGenPhotCand]/D" );
      tree->Branch( "GenPhotCand_phi", GenPhotCand_phi, "GenPhotCand_phi[nGenPhotCand]/D" );
      tree->Branch( "GenPhotCand_e", GenPhotCand_e, "GenPhotCand_e[nGenPhotCand]/D" );
      tree->Branch( "nGenProtCand", &nGenProtCand, "nGenProtCand/I" );
      tree->Branch( "GenProtCand_pt", GenProtCand_pt, "GenProtCand_pt[nGenProtCand]/D" );
      tree->Branch( "GenProtCand_eta", GenProtCand_eta, "GenProtCand_eta[nGenProtCand]/D" );
      tree->Branch( "GenProtCand_phi", GenProtCand_phi, "GenProtCand_phi[nGenProtCand]/D" );
      tree->Branch( "GenProtCand_e", GenProtCand_e, "GenProtCand_e[nGenProtCand]/D" );
      tree->Branch( "GenProtCand_status", GenProtCand_status, "GenProtCand_status[nGenProtCand]/I" );
    }

    // Primary vertices' information
    tree->Branch( "nPrimVertexCand", &nPrimVertexCand, "nPrimVertexCand/I" );
    tree->Branch( "nFilteredPrimVertexCand", &nFilteredPrimVertexCand, "nFilteredPrimVertexCand/I" );
    tree->Branch( "PrimVertexCand_id", PrimVertexCand_id, "PrimVertexCand_id[nPrimVertexCand]/I" );
    tree->Branch( "PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D" );
    tree->Branch( "PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D" );
    tree->Branch( "PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D" );
    tree->Branch( "PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D" );
    tree->Branch( "PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/I" );
    tree->Branch( "PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/I" );

    // Lepton pairs' information
    tree->Branch( "nPair", &nPair, "nPair/I" );
    tree->Branch( "Pair_lepton1", Pair_lepton1, "Pair_lepton1[nPair]/I" );
    tree->Branch( "Pair_lepton2", Pair_lepton2, "Pair_lepton2[nPair]/I" );
    tree->Branch( "Pair_mindist", Pair_mindist, "Pair_mindist[nPair]/D" );
    tree->Branch( "Pair_mass", Pair_mass, "Pair_mass[nPair]/D" );
    tree->Branch( "Pair_pt", Pair_pt, "Pair_pt[nPair]/D" );
    tree->Branch( "Pair_eta", Pair_eta, "Pair_eta[nPair]/D" );
    tree->Branch( "Pair_phi", Pair_phi, "Pair_phi[nPair]/D" );
    tree->Branch( "Pair_dpt", Pair_dpt, "Pair_dpt[nPair]/D" );
    tree->Branch( "Pair_dphi", Pair_dphi, "Pair_dphi[nPair]/D" );
    tree->Branch( "Pair_3Dangle", Pair_3Dangle, "Pair_3Dangle[nPair]/D" );
    tree->Branch( "Pair_extratracks1mm", Pair_extratracks1mm, "Pair_extratracks1mm[nPair]/I" );
    tree->Branch( "Pair_extratracks2mm", Pair_extratracks2mm, "Pair_extratracks2mm[nPair]/I" );
    tree->Branch( "Pair_extratracks3mm", Pair_extratracks3mm, "Pair_extratracks3mm[nPair]/I" );
    tree->Branch( "Pair_extratracks4mm", Pair_extratracks4mm, "Pair_extratracks4mm[nPair]/I" );
    tree->Branch( "Pair_extratracks5mm", Pair_extratracks5mm, "Pair_extratracks5mm[nPair]/I" );
    tree->Branch( "Pair_extratracks1cm", Pair_extratracks1cm, "Pair_extratracks1cm[nPair]/I" );
    tree->Branch( "Pair_extratracks2cm", Pair_extratracks2cm, "Pair_extratracks2cm[nPair]/I" );
    tree->Branch( "Pair_extratracks3cm", Pair_extratracks3cm, "Pair_extratracks3cm[nPair]/I" );
    tree->Branch( "Pair_extratracks4cm", Pair_extratracks4cm, "Pair_extratracks4cm[nPair]/I" );
    tree->Branch( "Pair_extratracks5cm", Pair_extratracks5cm, "Pair_extratracks5cm[nPair]/I" );
    tree->Branch( "Pair_extratracks10cm", Pair_extratracks10cm, "Pair_extratracks10cm[nPair]/I" );
    // Kalman dilepton vertex information
    tree->Branch( "KalmanVertexCand_x", KalmanVertexCand_x, "KalmanVertexCand_x[nPair]/D" );
    tree->Branch( "KalmanVertexCand_y", KalmanVertexCand_y, "KalmanVertexCand_y[nPair]/D" );
    tree->Branch( "KalmanVertexCand_z", KalmanVertexCand_z, "KalmanVertexCand_z[nPair]/D" );

    tree->Branch( "nPairGamma", &nPairGamma, "nPairGamma/I" );
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
      tree->Branch( "nLocalProtCand", &nLocalProtCand, "nLocalProtCand/I" );
      tree->Branch( "LocalProtCand_x", LocalProtCand_x, "LocalProtCand_x[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_y", LocalProtCand_y, "LocalProtCand_y[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_z", LocalProtCand_z, "LocalProtCand_z[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_xSigma", LocalProtCand_xSigma, "LocalProtCand_xSigma[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_ySigma", LocalProtCand_ySigma, "LocalProtCand_ySigma[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_arm", LocalProtCand_arm, "LocalProtCand_arm[nLocalProtCand]/I" );
      tree->Branch( "LocalProtCand_side", LocalProtCand_side, "LocalProtCand_side[nLocalProtCand]/I" );
      tree->Branch( "LocalProtCand_Tx", LocalProtCand_Tx, "LocalProtCand_Tx[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_Ty", LocalProtCand_Ty, "LocalProtCand_Ty[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_TxSigma", LocalProtCand_TxSigma, "LocalProtCand_TxSigma[nLocalProtCand]/D" );
      tree->Branch( "LocalProtCand_TySigma", LocalProtCand_TySigma, "LocalProtCand_TySigma[nLocalProtCand]/D" );
    }

    // Extra tracks on vertex's information
    tree->Branch( "nExtraTracks", &nExtraTracks, "nExtraTracks/I" );
    tree->Branch( "ExtraTrack_pair", ExtraTrack_pair, "ExtraTrack_pair[nExtraTracks]/I" );
    tree->Branch( "ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I" );
    tree->Branch( "ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/I" );
    tree->Branch( "ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I" );
    tree->Branch( "ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/I" );
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
    tree->Branch( "nQualityExtraTrack", &nQualityExtraTrack, "nQualityExtraTrack/I" );
    tree->Branch( "ClosestExtraTrack_vtxdxyz",ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz[nPair]/D" );
    tree->Branch( "ClosestExtraTrack_id",ClosestExtraTrack_id,"ClosestExtraTrack_id[nPair]/I" );
    tree->Branch( "ClosestHighPurityExtraTrack_vtxdxyz",ClosestHighPurityExtraTrack_vtxdxyz,"ClosestHighPurityExtraTrack_vtxdxyz[nPair]/D" );
    tree->Branch( "ClosestHighPurityExtraTrack_id",ClosestHighPurityExtraTrack_id,"ClosestHighPurityExtraTrack_id[nPair]/I" );

    // Jets/MET information
    tree->Branch( "nJetCand", &nJetCand, "nJetCand/I" );
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
}
