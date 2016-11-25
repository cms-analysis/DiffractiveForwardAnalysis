import FWCore.ParameterSet.Config as cms

ggll_aod = cms.EDAnalyzer(
    'GammaGammaLL',

    # General parameters
    leptonsType = cms.InputTag('electronmuon'),
    maxExtraTracks = cms.untracked.uint32(10000),
    sqrtS = cms.double(13.e3), # in GeV
    fetchProtons = cms.bool(False), # retrieve the TOTEM/PPS info from the files (data only!)
    printCandidates = cms.untracked.bool(False),
    useLegacyVertexing = cms.bool(False),

    # MC tweaks
    runOnMC = cms.untracked.bool(True),
    MCAcceptPtCut = cms.untracked.double(0.),
    MCAcceptEtaCut = cms.untracked.double(-1.),

    # HLT selection
    HLTMenuTag = cms.string('HLT'),
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # Input collections
    #muonTag = cms.InputTag("muons"), # RECO ; needs recompilation!
    muonTag = cms.InputTag("selectedPatMuons"), # PAT
    #electronTag = cms.InputTag("gsfElectrons"), # RECO ; needs recompilation!
    electronTag = cms.InputTag("selectedPatElectrons"), # PAT
    vertexTag = cms.InputTag('offlinePrimaryVertices'),
    trackTag = cms.InputTag('generalTracks'),
    jetTag = cms.InputTag('selectedPatJets'),
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

    # Lepton ID
    isoValInputTags = cms.VInputTag(
        cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
        cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
        cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
    ),
    #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
)
ggll = ggll_aod.clone() ## for backward-compatibility

ggll_aod_pflow = ggll_aod.clone(
    muonTag = cms.untracked.InputTag("selectedPatMuonsPFlow"),
    electronTag = cms.untracked.InputTag("selectedPatElectronsPFlow"),
    jetTag = cms.InputTag('selectedPatJetsPFlow'),
    metTag = cms.InputTag('pfMet'),
)

ggll_miniaod = ggll_aod.clone(
    vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    muonTag = cms.untracked.InputTag("slimmedMuons"), # PAT
    electronTag = cms.untracked.InputTag("slimmedElectrons"), # PAT
    jetTag = cms.InputTag('slimmedJetsAK8'),
    photonTag = cms.InputTag('slimmedPhotons'),
    metTag = cms.InputTag('slimmedMETs'), # PAT
    genParticleTag = cms.InputTag('prunedGenParticles'),
    pfTag = cms.untracked.InputTag('packedPFCandidates'),
    pileupInfo = cms.InputTag('slimmedAddPileupInfo'),
)
