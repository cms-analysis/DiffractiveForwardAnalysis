import FWCore.ParameterSet.Config as cms

gamgameeanalysis = cms.EDFilter("GammaGammaEE",
    ElectronCollectionLabel = cms.InputTag("selectedLayer1Electrons"),
    outfilename = cms.untracked.string('ee.pat.root'),
    DielectronMindphi = cms.double(0.0),
    JetCollectionLabel = cms.InputTag("selectedLayer1Jets"),
    PhotonCollectionLabel = cms.InputTag("selectedLayer1Photons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("selectedLayer1Muons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    CaloTowerdR = cms.double(0.3),
    DielectronMaxdEt = cms.double(2000.0),
    MetLabel = cms.InputTag("selectedLayer1METs")
)



