import FWCore.ParameterSet.Config as cms

gamgammumuanalysis = cms.EDAnalyzer("GammaGammaMuMu",
    ElectronCollectionLabel = cms.InputTag("selectedPatElectrons"),
    #ElectronCollectionLabel = cms.InputTag("selectedPatElectronsPFlow"),
    outfilename = cms.untracked.string('mumu.pat.root'),
    JetCollectionLabel = cms.InputTag("selectedPatJets"),
    #JetCollectionLabel = cms.InputTag("selectedPatJetsPFlow"),
    PhotonCollectionLabel = cms.InputTag("selectedPatPhotons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("selectedPatMuons"),
    #GlobalMuonCollectionLabel = cms.InputTag("selectedPatMuonsPFlow"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    CastorTowerLabel = cms.InputTag("CastorTowerReco"),
    ZDCRecHitsLabel = cms.InputTag("zdcreco"),                             
    CastorRecHitsLabel = cms.InputTag("castorreco"),
    CaloTowerdR = cms.double(0.3),
    DimuonMindphi = cms.double(0.0),
    MetLabel = cms.InputTag("pfMet"),
                                    
    DimuonMaxdpt = cms.double(2000.0),
    MinMuMuVertexSeparation = cms.double(0.1),
    KeepSameSignDimuons = cms.bool(False),
    ReadMCPileup = cms.bool(False),
    ReadMCEffCorrections = cms.bool(False),
    ReadMCEffCorrectionsByCharge = cms.bool(False),
    ReadmcEffCorrectionsBySignedEta = cms.bool(False),
    AlgoNames = cms.vstring('Excl_Data_TrackProbeP_JPsiZ',
                            'Excl_Data_TrackProbeM_JPsiZ',
                            'HLT_DoubleMu3_Data_TrackProbeP_JPsiZ',
                            'HLT_DoubleMu3_Data_TrackProbeM_JPsiZ',
                            'Excl_MC_TrackProbeP_JPsiZ',
                            'Excl_MC_TrackProbeM_JPsiZ',
                            'HLT_DoubleMu3_MC_TrackProbeP_JPsiZ',
                            'HLT_DoubleMu3_MC_TrackProbeM_JPsiZ'),
    HLTMenuLabel = cms.string("HLT")                                  
)



