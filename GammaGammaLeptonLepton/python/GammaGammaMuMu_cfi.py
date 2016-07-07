import FWCore.ParameterSet.Config as cms

ggmumu = cms.EDAnalyzer(
    'GammaGammaMuMu',
    SqrtS = cms.double(13000.),
    HLTMenuLabel = cms.string("HLT"),
    LeptonsType = cms.InputTag('muon'),
    RunOnMC = cms.untracked.bool(True),
    RunOnProtons = cms.untracked.bool(True),
    MCAcceptPtCut = cms.untracked.double(0.),
    MCAcceptEtaCut = cms.untracked.double(-1.),
    GenParticlesCollectionLabel = cms.InputTag('genParticles'),
    PrintCandidates = cms.untracked.bool(False),
)
