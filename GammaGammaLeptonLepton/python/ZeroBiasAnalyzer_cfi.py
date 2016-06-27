import FWCore.ParameterSet.Config as cms

zerobiasanalysis = cms.EDAnalyzer("ZeroBiasAnalyzer",
    outfilename = cms.untracked.string('zerobias_zdcthresh.root'),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    CastorTowerLabel = cms.InputTag("CastorTowerReco"),
    CastorRecHitsLabel = cms.InputTag("castorreco"),
    ZDCRecHitsLabel = cms.InputTag("zdcreco"),                                  
    RecoTrackLabel = cms.InputTag("generalTracks"),
    HLTMenuLabel = cms.string("HLT")
)



