import FWCore.ParameterSet.Config as cms

gamgammueanalysis = cms.EDAnalyzer("GammaGammaMuE",
    ElectronCollectionLabel = cms.InputTag("selectedPatElectrons"),
    outfilename = cms.untracked.string('mue.pat.root'),
    JetCollectionLabel = cms.InputTag("selectedPatJets"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("selectedPatMuons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    CastorTowerLabel = cms.InputTag("CastorTowerReco"),
    ZDCRecHitsLabel = cms.InputTag("zdcreco"),                             
    CastorRecHitsLabel = cms.InputTag("castorreco"),
    CaloTowerdR = cms.double(0.3),
##    MetLabel = cms.InputTag("met"),
    MetLabel = cms.InputTag("pfMet"),                                  
    MinMuEVertexSeparation = cms.double(0.1),
    KeepSameSign = cms.bool(False),
    HLTMenuLabel = cms.string("HLT"),
##    MCPileupDist = cms.untracked.string("PUHistos.root"),
##    DataPileupDist = cms.untracked.string("PUHistos.root"),
)



