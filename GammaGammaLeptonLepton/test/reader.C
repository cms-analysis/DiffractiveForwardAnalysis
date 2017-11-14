#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"

#include <iostream>

using namespace std;

void reader( const char* filename = "output.root" )
{
  auto file = TFile::Open( filename );
  auto tree = dynamic_cast<TTree*>( file->Get( "ggll_aod/ntp1" ) );
  ggll::AnalysisEvent evt;
  evt.load( tree, ggll::DiMuon, false );

  TH1D h_pair_mass( "h_pair_mass", "m(#mu^{+}#mu^{-})\\Events\\GeV?.2f", 200, 100., 2100. );

  for ( unsigned long long i = 0; i < tree->GetEntriesFast(); ++i ) {
    tree->GetEntry( i );
    //cout << ">>> event " << i << endl;
    /*cout << "- fired triggers:" << endl;
    for ( const auto& hlt : *evt.HLT_Name ) {
      cout << "  *) " << hlt << endl;
    }*/

    //cout << "- dilepton pairs:" << endl;
    for ( unsigned int j = 0; j < evt.nPair; ++j ) {
      //cout << "  *) pair " << j << " has invariant mass = " << evt.Pair_mass[j] << endl;

      if ( evt.Pair_extratracks0p5mm[j] != 0 ) continue;
      if ( fabs( evt.KalmanVertexCand_z[j] ) > 15. ) continue;
      if ( 1.-fabs( evt.Pair_dphi[j] )/M_PI > 0.009 ) continue;
      if ( evt.Pair_mass[j] < 110. ) continue;

      cout << "CANDIDATE!!!" << endl;

      const unsigned int l1 = evt.Pair_lepton1[j], l2 = evt.Pair_lepton2[j];
      const double xip = ( evt.MuonCand_pt[l1]*exp( +evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( +evt.MuonCand_eta[l2] ) ) / 13.e3;
      const double xim = ( evt.MuonCand_pt[l1]*exp( -evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( -evt.MuonCand_eta[l2] ) ) / 13.e3;
      cout << "central system xip/xim = " << xip << " / " << xim << endl;

      h_pair_mass.Fill( evt.Pair_mass[j] );

      //cout << evt.Pair_extratracks2mm[j] << endl;

      /*auto kvc = TVector3( evt.KalmanVertexCand_x[j], evt.KalmanVertexCand_y[j], evt.KalmanVertexCand_z[j] );
      double min_dist = 999.;
      int min_dist_vtx = -1;
      for ( unsigned int k = 0; k < evt.nPrimVertexCand; ++k ) {
        auto pvc = TVector3( evt.PrimVertexCand_x[k], evt.PrimVertexCand_y[k], evt.PrimVertexCand_z[k] );
        const double dist = ( kvc-pvc ).Mag();
        if ( dist < min_dist ) {
          min_dist = dist;
          min_dist_vtx = k;
        }
      }*/
    }
  }

  {
    Canvas c( "cand_pair_mass", "CMS+TOTEM Preliminary 2017, #sqrt{s} = 13 TeV" );
    h_pair_mass.Draw();
    c.Prettify( &h_pair_mass );
    c.Save( "pdf" );
  }
}
