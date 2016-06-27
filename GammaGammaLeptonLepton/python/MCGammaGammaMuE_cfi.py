import FWCore.ParameterSet.Config as cms

mcgamgammueanalysis = cms.EDAnalyzer("GammaGammaMuEMC",
    outfilename = cms.untracked.string('mue.genlevel.root'),
    FillAllMCParticles = cms.bool(True)
)



