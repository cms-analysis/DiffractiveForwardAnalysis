import FWCore.ParameterSet.Config as cms

gamgameeanalysis = cms.EDFilter("GammaGammaEE",
    ElectronCollectionLabel = cms.InputTag("pixelMatchGsfElectrons"),
    outfilename = cms.untracked.string('ee.recolevel.root'),
    DielectronMindphi = cms.double(0.0),
    JetCollectionLabel = cms.InputTag("sisCone5CaloJets"),
    PhotonCollectionLabel = cms.InputTag("correctedPhotons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("muons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    CaloTowerdR = cms.double(0.3),
    DielectronMaxdEt = cms.double(2000.0),
    MetLabel = cms.InputTag("met")
)



