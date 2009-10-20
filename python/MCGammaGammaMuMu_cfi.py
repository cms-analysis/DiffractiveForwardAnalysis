import FWCore.ParameterSet.Config as cms

mcgamgammumuanalysis = cms.EDFilter("GammaGammaMuMuMC",
    outfilename = cms.untracked.string('mumu.genlevel.root'),
    FillAllMCParticles = cms.bool(True)
)



