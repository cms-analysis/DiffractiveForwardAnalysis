#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMu.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaEE.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMuMC.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaEEMC.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/CosmicsMuMu.h"

DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(GammaGammaMuMu);
DEFINE_ANOTHER_FWK_MODULE(GammaGammaEE);
DEFINE_FWK_MODULE(GammaGammaMuMuMC);
DEFINE_ANOTHER_FWK_MODULE(GammaGammaEEMC);
DEFINE_ANOTHER_FWK_MODULE(CosmicsMuMu);
