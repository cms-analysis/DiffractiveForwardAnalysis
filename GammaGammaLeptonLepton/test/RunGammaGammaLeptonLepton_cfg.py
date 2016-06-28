import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

runOnMC = False
useAOD = True # AOD or MiniAOD?

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(250) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/20000/002F5062-346F-E211-BF00-1CC1DE04DF20.root',
      #'/store/mc/Run2015D/DoubleEG/AOD/04Dec2015-v1/10000/04D11E1B-BB9E-E511-AC1A-047D7B881D62.root',
      #'/store/data/Run2016B/DoubleMuon/AOD/PromptReco-v2/000/273/150/00000/3E460221-D919-E611-AE4F-02163E014142.root',
      #'/store/data/Run2016B/DoubleMuon/AOD/PromptReco-v2/000/273/158/00000/0EF12F56-EC19-E611-AF65-02163E01385D.root',
      #'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/080636C6-723B-E611-92AB-02163E01429D.root',
      '/store/data/Run2016B/DoubleEG/AOD/PromptReco-v2/000/273/158/00000/006772B7-E019-E611-AEBE-02163E014583.root',
    ),
    firstEvent = cms.untracked.uint32(540)
)

#########################
#        Triggers       #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ['HLT_Mu17_Mu8_*']
#process.hltFilter.HLTPaths = ['HLT_Mu10_Ele10_CaloIdL_*', 'HLT_Mu8_Ele17_*', 'HLT_Mu17_Ele8_*']
process.hltFilter.HLTPaths = ['HLT_Ele17_Ele12_*', 'HLT_Ele23_Ele12_*']

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START52_V9::All')
#process.GlobalTag.globaltag = cms.string('FT_R_53_V6::All')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#process.scrapingVeto = cms.EDFilter("FilterOutScraping",
#    applyfilter = cms.untracked.bool(True),
#    debugOn = cms.untracked.bool(False),
#    numtrack = cms.untracked.uint32(10),
#    thresh = cms.untracked.double(0.2)
#)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
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
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patExtraAodEventContent
#from PhysicsTools.PatAlgos.tools.coreTools import *

#
# rho value for isolation
#
#from RecoJets.JetProducers.kt4PFJets_cfi import *
#process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
#
# Particle flow isolation
#
#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
#process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offline*PrimaryVertices*_*_*',
        'keep *_*Muons*_*_*',
        'keep *_*Electrons*_*_*',
        'keep *_*Photons*_*_*',
        'keep *_*Jets*_*_*',
        'keep *_*MET*_*_*',
        'keep recoPFCandidates_particleFlow_*_*',
        #*patEventContentNoCleaning
    ),
)

# Particle flow
#from PhysicsTools.PatAlgos.tools.pfTools import *
#postfix = "PFlow"
##jetAlgo="AK5" ##FIXME!! 
#jetAlgo="AK4"
##usePFBRECO(process,runPFBRECO=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix) 
#usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix)
#useGsfElectrons(process,postfix)
#removeCleaning(process)
#removeMCMatching(process, ['All'])

#from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
from PhysicsTools.PatAlgos.tools.coreTools import runOnData, removeMCMatching

if not runOnMC:
    names = ['Photons', 'Electrons', 'Muons', 'Jets', 'METs']
    #runOnData(process, names=names)
    removeMCMatching(process, names=names)

#process.tauMatchBoostedPFlow = cms.PSet()
#process.patTausBoostedPFlow = cms.PSet()
#usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix)
#useGsfElectrons(process,postfix)

#########################
#      Electron ID      #
#########################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, setupVIDElectronSelection, setupAllVIDIdsInModule, DataFormat

switchOnVIDElectronIdProducer(process, DataFormat.AOD)
setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff', setupVIDElectronSelection)

#########################
#       Analysis        #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")

process.ggll.TriggersList = process.hltFilter.HLTPaths
#process.ggll.LeptonsType = cms.vstring('Muon')
#process.ggll.LeptonsType = cms.vstring('Electron', 'Muon')
process.ggll.LeptonsType = cms.vstring('Electron')
process.ggll.RecoVertexLabel = cms.InputTag('offlinePrimaryVertices')
#process.ggll.GlobalMuonCollectionLabel = cms.untracked.InputTag("muons") # RECO
#process.ggll.GlobalMuonCollectionLabel = cms.untracked.InputTag("selectedPatMuonsPFlow"), # PAT (particle flow)
process.ggll.GlobalMuonCollectionLabel = cms.untracked.InputTag("selectedPatMuons") # PAT
#process.ggll.GlobalEleCollectionLabel = cms.untracked.InputTag("gsfElectrons") # RECO
#process.ggll.GlobalEleCollectionLabel = cms.untracked.InputTag("selectedPatElectronsPFlow") # PAT (particle flow)
process.ggll.GlobalEleCollectionLabel = cms.untracked.InputTag("selectedPatElectrons") # PAT
process.ggll.RunOnMC = cms.untracked.bool(runOnMC)
process.ggll.outfilename = cms.untracked.string('output.root')

process.p = cms.Path(
    #process.scrapingVeto*
    #process.kt6PFJetsForIsolation*
    #process.pfiso*
    process.primaryVertexFilter*
    process.patDefaultSequence*
    process.egmGsfElectronIDSequence*
    #getattr(process,"patPF2PATSequence"+postfix)*
    process.ggll
)
#getattr(process,"pfNoElectron"+postfix).enable = True
