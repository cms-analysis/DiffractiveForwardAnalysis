import FWCore.ParameterSet.Config as cms

singlePhotonTreeProducer = cms.EDAnalyzer('SinglePhotonTreeProducer',
    hltMenuTag = cms.string('HLT'),
    runOnMC = cms.bool(True),
    triggersList = cms.vstring(),
    photonsTag = cms.InputTag('slimmedPhotons'),
    ppsTracksTag = cms.InputTag('ctppsLocalTrackLiteProducer'),
)
