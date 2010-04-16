import FWCore.ParameterSet.Config as cms

gamgammumuanalysis = cms.EDFilter("GammaGammaMuMu",
    ElectronCollectionLabel = cms.InputTag("selectedPatElectrons"),
    outfilename = cms.untracked.string('mumu.pat.root'),
    JetCollectionLabel = cms.InputTag("selectedPatJets"),
    PhotonCollectionLabel = cms.InputTag("selectedPatPhotons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("selectedPatMuons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    CastorTowerLabel = cms.InputTag("CastorTowerReco"),
    ZDCRecHitsLabel = cms.InputTag("zdcreco"),                             
    CastorRecHitsLabel = cms.InputTag("castorreco"),
    CaloTowerdR = cms.double(0.3),
    DimuonMindphi = cms.double(0.0),
    MetLabel = cms.InputTag("met"),
    DimuonMaxdpt = cms.double(2000.0),
    KeepSameSignDimuons = cms.bool(False),
    AlgoNames = cms.vstring('TriggerMuonFromGlobalMuonZ', 'GlobalMuonFromTrackerTrackZ','TrackerTrackFromStandaloneMuonZ',
                            'TriggerMuonFromGlobalMuonJpsi', 'GlobalMuonFromTrackerTrackJpsi', 'TrackerTrackFromStandaloneMuonJpsi'),
    HLTMenuLabel = cms.string("HLT")                                  
)



