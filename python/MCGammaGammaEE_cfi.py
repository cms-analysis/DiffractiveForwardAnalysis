import FWCore.ParameterSet.Config as cms

mcgamgameeanalysis = cms.EDFilter("GammaGammaEEMC",
    outfilename = cms.untracked.string('ee.genlevel.root'),
    FillAllMCParticles = cms.bool(False)
)



