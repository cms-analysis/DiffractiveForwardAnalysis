import FWCore.ParameterSet.Config as cms

ggll_aod = cms.EDAnalyzer('GammaGammaLL',
    # General parameters
    leptonsType = cms.InputTag('ElectronMuon'),
    #maxExtraTracks = cms.untracked.uint32(10000),
    sqrtS = cms.double(13.e3), # in GeV
    fetchProtons = cms.bool(False), # retrieve the TOTEM/PPS info from the files (data only!)
    printCandidates = cms.bool(False),

    # MC tweaks
    runOnMC = cms.bool(True),
    MCAcceptPtCut = cms.untracked.double(0.),
    MCAcceptEtaCut = cms.untracked.double(-1.),

    # HLT selection
    triggerEvent = cms.InputTag('patTriggerEvent'),
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # Input collections
    #muonTag = cms.InputTag("muons"), # RECO ; needs recompilation!
    muonTag = cms.InputTag("patMuons"), # PAT
    #electronTag = cms.InputTag("gsfElectrons"), # RECO ; needs recompilation!
    electronTag = cms.InputTag("patElectrons"), # PAT
    vertexTag = cms.InputTag('offlinePrimaryVertices'),
    trackTag = cms.InputTag('generalTracks'),
    jetTag = cms.InputTag('patJets'),
    metTag = cms.InputTag('patMETs'),
    photonTag = cms.InputTag('selectedPatPhotons'),
    ppsLocalTrackTag = cms.InputTag('ctppsLocalTrackLiteProducer'),
    ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP"),
    ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP"),

    genParticleTag = cms.InputTag('genParticles'),

    # Pileup reweighting
    pileupInfo = cms.InputTag('addPileupInfo'),
    mcpufile = cms.string('PUHistos_mc.root'),
    mcpupath = cms.string('pileup'),
    datapufile = cms.string('PUHistos_data.root'),
    datapupath = cms.string('pileup'),
    fixedGridRhoFastjetAllLabel = cms.InputTag('fixedGridRhoFastjetAll'),
)

ggll = ggll_aod.clone() ## for backward-compatibility

ggll_aod_pflow = ggll_aod.clone(
    muonTag = cms.InputTag("selectedPatMuonsPFlow"),
    electronTag = cms.InputTag("selectedPatElectronsPFlow"),
    jetTag = cms.InputTag('selectedPatJetsPFlow'),
    metTag = cms.InputTag('pfMet'),
)

ggll_miniaod = ggll_aod.clone(
    vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    muonTag = cms.InputTag("slimmedMuons"), # PAT
    electronTag = cms.InputTag("slimmedElectrons"), # PAT
    jetTag = cms.InputTag('slimmedJetsAK8'),
    photonTag = cms.InputTag('slimmedPhotons'),
    metTag = cms.InputTag('slimmedMETs'), # PAT
    genParticleTag = cms.InputTag('prunedGenParticles'),
    #pfTag = cms.InputTag('packedPFCandidates'),
    pileupInfo = cms.InputTag('slimmedAddPileupInfo'),
)
