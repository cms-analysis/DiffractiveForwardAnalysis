#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include <iostream>

using namespace std;

void reader( const char* filename = "output.root" )
{
  auto file = TFile::Open( filename );
  auto tree = dynamic_cast<TTree*>( file->Get( "ggll_aod/ntp1" ) );
  ggll::AnalysisEvent evt;
  evt.load( tree, ggll::DiMuon, false );

  for ( unsigned long long i = 0; i < tree->GetEntriesFast(); ++i ) {
    tree->GetEntry( i );
    for ( unsigned int j = 0; j < evt.nPair; ++j ) {
      cout << "---> pair " << j << " has invariant mass = " << evt.Pair_mass[j] << endl;
    }
  }
}
