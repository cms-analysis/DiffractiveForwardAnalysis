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

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
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
      #'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/666D4E0A-6F3B-E611-B6C0-02163E01383E.root',
      #'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/080636C6-723B-E611-92AB-02163E01429D.root',
      #'/store/data/Run2016B/DoubleEG/AOD/PromptReco-v2/000/273/158/00000/006772B7-E019-E611-AEBE-02163E014583.root',
      #'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/3A85B622-713B-E611-8BA0-02163E011A7C.root',
      '/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/778/00000/10EEA4DB-C53C-E611-9676-02163E01374C.root',
      #'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/777/00000/3EFB3C3E-903C-E611-BB5B-02163E01378A.root',
      #'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/276/283/00000/EC53DC17-E544-E611-8625-02163E011E99.root',
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ['HLT_Mu17_Mu8_*']
process.hltFilter.HLTPaths = cms.vstring(
    'HLT_DoubleMu33NoFiltersNoVtx_v*', 'HLT_DoubleMu38NoFiltersNoVtx_v*',
    'HLT_Mu17_Mu8_v*', #'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*',
    'HLT_Mu20_Mu10_v*'
)
#process.hltFilter.HLTPaths = ['HLT_Mu10_Ele10_CaloIdL_*', 'HLT_Mu8_Ele17_*', 'HLT_Mu17_Ele8_*']
#process.hltFilter.HLTPaths = ['HLT_Ele17_Ele12_*', 'HLT_Ele23_Ele12_*']

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

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

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offline*PrimaryVertices*_*_*',
        'keep *_selectedPatMuons*_*_*',
        'keep *_*lectron*_*_*',
        'keep *_selectedPatElectrons*_*_*',
        'keep *_selectedPatPhotons*_*_*',
        'keep *_selectedPatJets*_*_*',
        'keep *_*MET*_*_*',
        'keep *_*particleFlow*_*_*',
        #*patEventContentNoCleaning
    ),
)

from DiffractiveForwardAnalysis.GammaGammaLeptonLepton.RemovePATMCMatching_cfi import removePATMCMatching

if not runOnMC:
    #names = ['Photons', 'Electrons', 'Muons', 'Jets', 'METs']
    removePATMCMatching(process)

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

process.ggll_aod.TriggersList = process.hltFilter.HLTPaths
process.ggll_aod.LeptonsType = cms.vstring('Muon')
#process.ggll_aod.LeptonsType = cms.vstring('Electron', 'Muon')
#process.ggll_aod.LeptonsType = cms.vstring('Electron')
process.ggll_aod.RunOnMC = cms.untracked.bool(runOnMC)
process.ggll_aod.outfilename = cms.untracked.string('output.root')
process.ggll_aod.fetchProtons = cms.bool(True)

process.p = cms.Path(
    process.hltFilter*
    #process.scrapingVeto*
    #process.kt6PFJetsForIsolation*
    #process.pfiso*
    #process.primaryVertexFilter*
    process.egmGsfElectronIDSequence*
    process.patDefaultSequence*
    #getattr(process,"patPF2PATSequence"+postfix)*
    process.ggll_aod
)
