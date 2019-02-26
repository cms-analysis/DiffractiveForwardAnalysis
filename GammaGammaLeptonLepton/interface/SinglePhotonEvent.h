#ifndef DiffractiveForwardAnalysis_GammaGammaLeptonLepton_SinglePhotonEvent_h
#define DiffractiveForwardAnalysis_GammaGammaLeptonLepton_SinglePhotonEvent_h

#include <vector>
#include <string>

#include "TTree.h"

namespace gggx
{
  class SinglePhotonEvent
  {
    public:
      SinglePhotonEvent() : pHLT_Name( nullptr ) {
        clear();
      }
      ~SinglePhotonEvent() {}

      //////  Leaves size  //////

      /// Maximum number of HLT to check
      static constexpr unsigned int MAX_HLT = 10;
      /// Maximum number of photons per event
      static constexpr unsigned int MAX_PHO = 50;
      /// Maximum number of primary vertices per event
      static constexpr unsigned int MAX_VTX = 150;
      /// Maximum number of generator level photons per event
      static constexpr unsigned int MAX_GENPHO = 10;
      /// Maximum number of generator level protons per event
      static constexpr unsigned int MAX_GENPRO = 8;
      /// Maximum number of jets per event
      static constexpr unsigned int MAX_JETS = 40;
      /// Maximum number of reconstructed local tracks in RPs
      static constexpr unsigned int MAX_FWDTRKCAND = 25;

      ////// Tree contents //////

      // Run/event quantities
      unsigned int BX, Run, LumiSection, EventNum;
      //int LHCFillNum, LHCBeamMode;
      //double AvgInstDelLumi, BunchInstLumi[3];

      // HLT quantities
      unsigned int nHLT;
      int HLT_Accept[MAX_HLT], HLT_Prescl[MAX_HLT];
      std::vector<std::string> HLT_Name, *pHLT_Name;

      // Generator level quantities
      unsigned int nGenPhotCand, nGenPhotCandOutOfAccept;
      double GenPhotCand_pt[MAX_GENPHO], GenPhotCand_eta[MAX_GENPHO], GenPhotCand_phi[MAX_GENPHO], GenPhotCand_e[MAX_GENPHO];
      unsigned int nGenProtCand;
      double GenProtCand_pt[MAX_GENPRO], GenProtCand_eta[MAX_GENPRO], GenProtCand_phi[MAX_GENPRO], GenProtCand_e[MAX_GENPHO];
      int GenProtCand_status[MAX_GENPRO];

      // Pileup reweighting quantities
      double PUWeightTrue, Weight;

      // Photon quantities
      unsigned int nPhotonCand;
      double PhotonCand_pt[MAX_PHO], PhotonCand_eta[MAX_PHO], PhotonCand_phi[MAX_PHO], PhotonCand_e[MAX_PHO];
      double PhotonCand_r9[MAX_PHO];
      double PhotonCand_drtrue[MAX_PHO], PhotonCand_detatrue[MAX_PHO], PhotonCand_dphitrue[MAX_PHO];
      int PhotonCand_wp80id[MAX_PHO], PhotonCand_wp90id[MAX_PHO];
      int PhotonCand_electronveto[MAX_PHO], PhotonCand_pixelseed[MAX_PHO];

      // Vertex quantities
      unsigned int nPrimVertexCand;
      double PrimVertexCand_x[MAX_VTX], PrimVertexCand_y[MAX_VTX], PrimVertexCand_z[MAX_VTX];
      unsigned int PrimVertexCand_tracks[MAX_VTX];
      double PrimVertexCand_chi2[MAX_VTX];
      unsigned int PrimVertexCand_ndof[MAX_VTX];
      unsigned int nFilteredPrimVertexCand;

      // Jets/MET quantities
      unsigned int nJetCand;
      double JetCand_pt[MAX_JETS], JetCand_eta[MAX_JETS], JetCand_phi[MAX_JETS], JetCand_e[MAX_JETS];
      double HighestJet_pt, HighestJet_eta, HighestJet_phi, HighestJet_e;
      double SumJet_e;
      double Etmiss, Etmiss_phi, Etmiss_significance;

      // PPS quantities
      unsigned int nFwdTrkCand;
      double FwdTrkCand_x[MAX_FWDTRKCAND], FwdTrkCand_y[MAX_FWDTRKCAND];
      double FwdTrkCand_xSigma[MAX_FWDTRKCAND], FwdTrkCand_ySigma[MAX_FWDTRKCAND];
      int FwdTrkCand_arm[MAX_FWDTRKCAND], FwdTrkCand_station[MAX_FWDTRKCAND], FwdTrkCand_pot[MAX_FWDTRKCAND];

      void clear() {
        // event-level branches
        BX = Run = LumiSection = EventNum = 0;

        // high-level trigger
        nHLT = 0;
        HLT_Name.clear();
        for ( unsigned int i = 0; i < MAX_HLT; ++i )
          HLT_Accept[i] = HLT_Prescl[i] = -1;

        // gen-level information
        nGenPhotCand = nGenPhotCandOutOfAccept = 0;
        for ( unsigned int i = 0; i < MAX_GENPHO; ++i )
          GenPhotCand_pt[i] = GenPhotCand_eta[i] = GenPhotCand_phi[i] = GenPhotCand_e[i] = -999.;
        nGenProtCand = 0;
        for ( unsigned int i = 0; i < MAX_GENPRO; ++i ) {
          GenProtCand_pt[i] = GenProtCand_eta[i] = GenProtCand_phi[i] = GenProtCand_e[i] = -999.;
          GenProtCand_status[i] = -1;
        }

        PUWeightTrue = Weight = 0.;

        //LHCFillNum = LHCBeamMode = -1;

        // single photon candidates
        nPhotonCand = 0;
        for ( unsigned int i = 0; i < MAX_PHO; ++i ) {
          PhotonCand_pt[i] = PhotonCand_eta[i] = PhotonCand_phi[i] = PhotonCand_e[i] = -999.;
          PhotonCand_r9[i] = -999.;
          PhotonCand_drtrue[i] = PhotonCand_detatrue[i] = PhotonCand_dphitrue[i] = -999.;
          PhotonCand_wp80id[i] = PhotonCand_wp90id[i] = -1;
          PhotonCand_electronveto[i] = PhotonCand_pixelseed[i] = -1;
        }

        // offline primary vertices
        nPrimVertexCand = 0;
        for ( unsigned int i = 0; i < MAX_VTX; ++i ) {
          PrimVertexCand_x[i] = PrimVertexCand_y[i] = PrimVertexCand_z[i] = -999.;
          PrimVertexCand_tracks[i] = 0;
          PrimVertexCand_chi2[i] = -999.;
          PrimVertexCand_ndof[i] = 0;
        }
        nFilteredPrimVertexCand = 0;

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
        nFwdTrkCand = 0;
        for ( unsigned int i = 0; i < MAX_FWDTRKCAND; ++i ) {
          FwdTrkCand_x[i] = FwdTrkCand_y[i] = FwdTrkCand_xSigma[i] = FwdTrkCand_ySigma[i] = -999.;
          FwdTrkCand_arm[i] = FwdTrkCand_station[i] = FwdTrkCand_pot[i] = -1;
        }
      }
      void attach( TTree* tree, bool mc ) {
        if ( !tree )
          return;

        tree->Branch( "Run", &Run, "Run/i" );
        tree->Branch( "LumiSection", &LumiSection, "LumiSection/i" );
        tree->Branch( "BX", &BX, "BX/i" );
        tree->Branch( "EventNum", &EventNum, "EventNum/i" );

        tree->Branch( "nHLT", &nHLT, "nHLT/i" );
        tree->Branch( "HLT_Accept", HLT_Accept, "HLT_Accept[nHLT]/I" );
        tree->Branch( "HLT_Prescl", HLT_Prescl, "HLT_Prescl[nHLT]/I" );
        pHLT_Name = &HLT_Name;
        tree->Branch( "HLT_Name", "std::vector<std::string>", &pHLT_Name );

        tree->Branch( "nPhotonCand", &nPhotonCand, "nPhotonCand/i" );
        tree->Branch( "PhotonCand_pt", PhotonCand_pt, "PhotonCand_pt[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_eta", PhotonCand_eta, "PhotonCand_eta[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_phi", PhotonCand_phi, "PhotonCand_phi[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_e", PhotonCand_e, "PhotonCand_e[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_r9", PhotonCand_r9, "PhotonCand_r9[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_drtrue", PhotonCand_drtrue, "PhotonCand_drtrue[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_detatrue", PhotonCand_detatrue, "PhotonCand_detatrue[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_dphitrue", PhotonCand_dphitrue, "PhotonCand_dphitrue[nPhotonCand]/D" );
        tree->Branch( "PhotonCand_wp80id", PhotonCand_wp80id, "PhotonCand_wp80id[nPhotonCand]/I" );
        tree->Branch( "PhotonCand_wp90id", PhotonCand_wp90id, "PhotonCand_wp90id[nPhotonCand]/I" );
        tree->Branch( "PhotonCand_electronveto", PhotonCand_electronveto, "PhotonCand_electronveto[nPhotonCand]/I" );
        tree->Branch( "PhotonCand_pixelseed", PhotonCand_pixelseed, "PhotonCand_pixelseed[nPhotonCand]/I" );

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
        tree->Branch( "PrimVertexCand_x", PrimVertexCand_x, "PrimVertexCand_x[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_y", PrimVertexCand_y, "PrimVertexCand_y[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_z", PrimVertexCand_z, "PrimVertexCand_z[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_chi2", PrimVertexCand_chi2, "PrimVertexCand_chi2[nPrimVertexCand]/D" );
        tree->Branch( "PrimVertexCand_ndof", PrimVertexCand_ndof, "PrimVertexCand_ndof[nPrimVertexCand]/i" );
        tree->Branch( "PrimVertexCand_tracks", PrimVertexCand_tracks, "PrimVertexCand_tracks[nPrimVertexCand]/i" );

        if ( !mc ) {
          tree->Branch( "nFwdTrkCand", &nFwdTrkCand, "nFwdTrkCand/i" );
          tree->Branch( "FwdTrkCand_x", FwdTrkCand_x, "FwdTrkCand_x[nFwdTrkCand]/D" );
          tree->Branch( "FwdTrkCand_y", FwdTrkCand_y, "FwdTrkCand_y[nFwdTrkCand]/D" );
          tree->Branch( "FwdTrkCand_xSigma", FwdTrkCand_xSigma, "FwdTrkCand_xSigma[nFwdTrkCand]/D" );
          tree->Branch( "FwdTrkCand_ySigma", FwdTrkCand_ySigma, "FwdTrkCand_ySigma[nFwdTrkCand]/D" );
          tree->Branch( "FwdTrkCand_arm", FwdTrkCand_arm, "FwdTrkCand_arm[nFwdTrkCand]/I" );
          tree->Branch( "FwdTrkCand_station", FwdTrkCand_station, "FwdTrkCand_station[nFwdTrkCand]/I" );
          tree->Branch( "FwdTrkCand_pot", FwdTrkCand_pot, "FwdTrkCand_pot[nFwdTrkCand]/I" );
        }

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
      void load( TTree* tree, bool mc ) {
        if ( !tree ) return;

        tree->SetBranchAddress( "Run", &Run );
        tree->SetBranchAddress( "LumiSection", &LumiSection );
        tree->SetBranchAddress( "BX", &BX );
        tree->SetBranchAddress( "EventNum", &EventNum );

        tree->SetBranchAddress( "nHLT", &nHLT );
        tree->SetBranchAddress( "HLT_Accept", HLT_Accept );
        tree->SetBranchAddress( "HLT_Prescl", HLT_Prescl );
        tree->SetBranchAddress( "HLT_Name", &pHLT_Name );

        tree->SetBranchAddress( "nPhotonCand", &nPhotonCand );
        tree->SetBranchAddress( "PhotonCand_pt", PhotonCand_pt );
        tree->SetBranchAddress( "PhotonCand_eta", PhotonCand_eta );
        tree->SetBranchAddress( "PhotonCand_phi", PhotonCand_phi );
        tree->SetBranchAddress( "PhotonCand_e", PhotonCand_e );
        tree->SetBranchAddress( "PhotonCand_r9", PhotonCand_r9 );
        tree->SetBranchAddress( "PhotonCand_drtrue", PhotonCand_drtrue );
        tree->SetBranchAddress( "PhotonCand_detatrue", PhotonCand_detatrue );
        tree->SetBranchAddress( "PhotonCand_dphitrue", PhotonCand_dphitrue );
        tree->SetBranchAddress( "PhotonCand_wp80id", PhotonCand_wp80id );
        tree->SetBranchAddress( "PhotonCand_wp90id", PhotonCand_wp90id );
        tree->SetBranchAddress( "PhotonCand_electronveto", PhotonCand_electronveto );
        tree->SetBranchAddress( "PhotonCand_pixelseed", PhotonCand_pixelseed );

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
        tree->SetBranchAddress( "PrimVertexCand_x", PrimVertexCand_x );
        tree->SetBranchAddress( "PrimVertexCand_y", PrimVertexCand_y );
        tree->SetBranchAddress( "PrimVertexCand_z", PrimVertexCand_z );
        tree->SetBranchAddress( "PrimVertexCand_chi2", PrimVertexCand_chi2 );
        tree->SetBranchAddress( "PrimVertexCand_ndof", PrimVertexCand_ndof );
        tree->SetBranchAddress( "PrimVertexCand_tracks", PrimVertexCand_tracks );

        if ( !mc ) {
          tree->SetBranchAddress( "nFwdTrkCand", &nFwdTrkCand );
          tree->SetBranchAddress( "FwdTrkCand_x", FwdTrkCand_x );
          tree->SetBranchAddress( "FwdTrkCand_y", FwdTrkCand_y );
          tree->SetBranchAddress( "FwdTrkCand_xSigma", FwdTrkCand_xSigma );
          tree->SetBranchAddress( "FwdTrkCand_ySigma", FwdTrkCand_ySigma );
          tree->SetBranchAddress( "FwdTrkCand_arm", FwdTrkCand_arm );
          tree->SetBranchAddress( "FwdTrkCand_station", FwdTrkCand_station );
          tree->SetBranchAddress( "FwdTrkCand_pot", FwdTrkCand_pot );
        }

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

