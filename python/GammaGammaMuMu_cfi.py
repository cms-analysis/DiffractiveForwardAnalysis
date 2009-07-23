import FWCore.ParameterSet.Config as cms

gamgammumuanalysis = cms.EDFilter("GammaGammaMuMu",
    ElectronCollectionLabel = cms.InputTag("pixelMatchGsfElectrons"),
    outfilename = cms.untracked.string('mumu.recolevel.root'),
    JetCollectionLabel = cms.InputTag("sisCone5CaloJets"),
    PhotonCollectionLabel = cms.InputTag("correctedPhotons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("muons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    CastorTowerLabel = cms.InputTag("CastorFastTowerReco"),
    CaloTowerdR = cms.double(0.3),
    DimuonMindphi = cms.double(0.0),
    MetLabel = cms.InputTag("met"),
    DimuonMaxdpt = cms.double(2000.0),
    HLTMenuLabel = cms.string("HLT8E29")                                  
)



