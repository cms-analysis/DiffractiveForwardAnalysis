import FWCore.ParameterSet.Config as cms

singlePhotonTreeProducer = cms.EDAnalyzer('SinglePhotonTreeProducer',
    #--- general steering parameters
    runOnMC = cms.bool(True),
    minPhotonMult = cms.uint32(1), # minimum number of photons in the event
    minFwdTrks = cms.uint32(1), # minimum number of tracks in PPS stations
    #--- input collections
    hltMenuTag = cms.string('HLT'),
    triggersList = cms.vstring(),
    photonsTag = cms.InputTag('slimmedPhotons'),
    ppsTracksTag = cms.InputTag('ctppsLocalTrackLiteProducer'),
    verticesTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jetsTag = cms.InputTag('slimmedJetsAK8'),
    metsTag = cms.InputTag('slimmedMETs'),
)
