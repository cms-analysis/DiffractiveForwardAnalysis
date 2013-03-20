import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

#########################
#    General options    #
#########################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")
process.options   = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/data/Run2012A/MuEG/AOD/20Nov2012-v2/00000/2ED48DC8-CA34-E211-86FF-003048FFCBFC.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/00C95B60-E470-E211-B576-00261894382D.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/020ABE37-0371-E211-8D8B-003048679168.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/02189C21-F170-E211-AADA-00259059391E.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/02ABD60A-F970-E211-874A-0026189438BD.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/02F0D309-F970-E211-AEF9-0026189438F8.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/0458153E-FB70-E211-B707-00261894389C.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/061ED436-0371-E211-AEC7-003048FFD71A.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/02B908B0-0471-E211-89C0-0025905964C2.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/02DE7336-0271-E211-8F3D-002618943935.root',
#'/store/mc/Summer12_DR53X/Wbb_FullyHadronic_8TeV_madgraph/AODSIM/PU_S10_START53_V7C-v1/30000/0440159C-0571-E211-B578-003048678E92.root',
#      'file:/afs/cern.ch/work/l/lforthom/private/18F94A2F-D502-E111-AF68-003048678BAC.root',
#      'rfio:/castor/cern.ch/user/j/jjhollar/ExclWW/GamGamWW_SM_WithFF_START44_Fall11PU_RECO.root',
    )
)

#########################
#      Preskimming      #
#########################
##process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
## need to include Global Tag stuff!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START52_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(15),
    maxd0 = cms.double(2)
)

process.muonFilter=cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("muons"),
    minNumber = cms.uint32(1)
)

#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from Configuration.EventContent.EventContent_cff import *
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patExtraAodEventContent
from PhysicsTools.PatAlgos.tools.coreTools import *

#
# rho value for isolation
#
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
#
# particle flow isolation
#
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offlinePrimaryVertices*_*_*',
        'keep *_*Muons*_*_*',
        'keep *_*Electrons*_*_*',
        #'keep *_*Photons*_*_*',
        #'keep *_*Jets*_*_*',
        #'keep recoPFCandidates_particleFlow_*_*'
    ),
)

#removeCleaning(process)
removeMCMatching(process, ['All'])

# JH - testing
#from PhysicsTools.PatAlgos.tools.pfTools import *
#postfix = "PFlow"
#usePF2PAT(process,
#    runPF2PAT = True,
#    jetAlgo = "AK5",
#    runOnMC = False,
#    postfix = postfix
#)
# end JH
#useGsfElectrons(process, postfix)

#########################
#       Analysis        #
#########################
#process.ggll = cms.EDAnalyzer('GammaGammaLL',
#process.ggll.LeptonsType = cms.vstring('Muon')
process.ggll.LeptonsType = cms.vstring('Electron', 'Muon')
#process.ggll.LeptonsType = cms.vstring('Electron')
process.ggll.RecoVertexLabel = cms.InputTag("offlinePrimaryVertices")
#process.ggll.GlobalMuonCollectionLabel = cms.untracked.InputTag("muons") # RECO
#process.ggll.GlobalMuonCollectionLabel = cms.untracked.InputTag("selectedPatMuonsPFlow"), # PAT (particle flow)
process.ggll.GlobalMuonCollectionLabel = cms.untracked.InputTag("selectedPatMuons") # PAT
#process.ggll.GlobalEleCollectionLabel = cms.untracked.InputTag("gsfElectrons") # RECO
#process.ggll.GlobalEleCollectionLabel = cms.untracked.InputTag("selectedPatElectronsPFlow") # PAT (particle flow)
process.ggll.GlobalEleCollectionLabel = cms.untracked.InputTag("selectedPatElectrons") # PAT

process.p = cms.Path(
    #process.scrapingVeto+
    process.kt6PFJetsForIsolation+
    process.pfiso+
    process.primaryVertexFilter+
    process.patDefaultSequence+
    #getattr(process,"patPF2PATSequence"+postfix)+
    process.ggll
)
