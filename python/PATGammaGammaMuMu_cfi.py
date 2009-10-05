import FWCore.ParameterSet.Config as cms

gamgammumuanalysis = cms.EDFilter("GammaGammaMuMu",
    ElectronCollectionLabel = cms.InputTag("selectedLayer1Electrons"),
    outfilename = cms.untracked.string('mumu.pat.root'),
    JetCollectionLabel = cms.InputTag("selectedLayer1Jets"),
    PhotonCollectionLabel = cms.InputTag("selectedLayer1Photons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("selectedLayer1Muons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    CastorTowerLabel = cms.InputTag("CastorFastTowerReco"),
    CaloTowerdR = cms.double(0.3),
    DimuonMindphi = cms.double(0.0),
    MetLabel = cms.InputTag("met"),
    DimuonMaxdpt = cms.double(2000.0),
    KeepSameSignDimuons = cms.bool(False),
    AlgoNames = cms.vstring('TriggerMuonFromGlobalMuonJpsi', 'GlobalMuonFromTrackerTrackJpsi', 'TrackerTrackFromStandaloneMuonJpsi'),
    HLTMenuLabel = cms.string("HLT8E29")                                  
)



