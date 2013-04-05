import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

#########################
#    General options    #
#########################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/data/Run2012A/MuEG/AOD/20Nov2012-v2/00000/2ED48DC8-CA34-E211-86FF-003048FFCBFC.root',
    )
)

#########################
#        Triggers       #
#########################
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_Mu10_Ele10_CaloIdL_*', 'HLT_Mu8_Ele17_*', 'HLT_Mu17_Ele8_*']

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START52_V9::All')
process.GlobalTag.globaltag = cms.string('FT_R_53_V6::All')
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

process.muonFilter = cms.EDFilter("CandViewCountFilter",
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
# Particle flow isolation
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
    ),
)

#removeCleaning(process)
removeMCMatching(process, ['All'])

#########################
#       Analysis        #
#########################
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")
process.ggll.TriggersList = process.hltFilter.HLTPaths
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
process.ggll.RunOnMC = cms.untracked.bool(False)
process.ggll.outfilename = cms.untracked.string('output.root')

process.p = cms.Path(
    #process.scrapingVeto+
    process.kt6PFJetsForIsolation+
    process.pfiso+
    process.primaryVertexFilter+
    process.patDefaultSequence+
    process.ggll
)
