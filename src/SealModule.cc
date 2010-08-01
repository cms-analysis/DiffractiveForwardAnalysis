#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMu.h"
//#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/CosmicsMuMu.h"
//#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/CollisionsMuMu.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ExclusiveTrackTrack.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ExclusiveDDbar.h" 
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ExclusiveDDbarSemiLept.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ZeroBiasAnalyzer.h"

//DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(GammaGammaMuMu);
//DEFINE_FWK_MODULE(CosmicsMuMu);
//DEFINE_FWK_MODULE(CollisionsMuMu);
DEFINE_FWK_MODULE(ExclusiveTrackTrack);
DEFINE_FWK_MODULE(ExclusiveDDbar); 
DEFINE_FWK_MODULE(ExclusiveDDbarSemiLept);
DEFINE_FWK_MODULE(ZeroBiasAnalyzer);
