import FWCore.ParameterSet.Config as cms

ggll_aod = cms.EDAnalyzer(
    'GammaGammaLL',

    # General parameters
    leptonsType = cms.InputTag('electronmuon'),
    #maxExtraTracks = cms.untracked.uint32(10000),
    sqrtS = cms.double(13.e3), # in GeV
    fetchProtons = cms.bool(False), # retrieve the TOTEM/PPS info from the files (data only!)
    printCandidates = cms.bool(False),
    useLegacyVertexing = cms.bool(False),

    # MC tweaks
    runOnMC = cms.bool(True),
    MCAcceptPtCut = cms.untracked.double(0.),
    MCAcceptEtaCut = cms.untracked.double(-1.),

    # HLT selection
    HLTMenuTag = cms.string('HLT'),
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
    totemRPLocalTrackTag = cms.InputTag('totemRPLocalTrackFitter'),
    genParticleTag = cms.InputTag('genParticles'),

    # Pileup reweighting
    pileupInfo = cms.InputTag('addPileupInfo'),
    mcpufile = cms.string('PUHistos_mc.root'),
    mcpupath = cms.string('pileup'),
    datapufile = cms.string('PUHistos_data.root'),
    datapupath = cms.string('pileup'),

    # Electron ID
    fixedGridRhoFastjetAllLabel = cms.InputTag('fixedGridRhoFastjetAll'),
    eleIdLabels = cms.PSet(
       looseLabel = cms.InputTag('cutBasedElectronID-Spring15-25ns-V1-standalone-loose'),
       mediumLabel = cms.InputTag('cutBasedElectronID-Spring15-25ns-V1-standalone-medium'),
       tightLabel = cms.InputTag('cutBasedElectronID-Spring15-25ns-V1-standalone-tight'),
       vetoLabel = cms.InputTag('cutBasedElectronID-Spring15-25ns-V1-standalone-veto'),
    ),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    #eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    #eleLooseMVAIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    #eleTightMVAIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
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
