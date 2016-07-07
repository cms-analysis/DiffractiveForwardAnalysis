import FWCore.ParameterSet.Config as cms

process = cms.Process("ggmumu")

#runOnMC = True
runOnMC = False
runOnProtons = True

#########################
#    General options    #
#########################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/080636C6-723B-E611-92AB-02163E01429D.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/16858380-6F3B-E611-AE0A-02163E0118AD.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/3A85B622-713B-E611-8BA0-02163E011A7C.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/4ADDCF1D-6D3B-E611-8674-02163E0138F1.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/582C24A7-733B-E611-88D2-02163E0133C5.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/666D4E0A-6F3B-E611-B6C0-02163E01383E.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/7E013031-6F3B-E611-9FCF-02163E0137D8.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/A433FE41-743B-E611-B494-02163E0133AD.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/BE0F4FD7-793B-E611-B6DB-02163E0140FB.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/BE7447BE-803B-E611-BA73-02163E014417.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/C8AF2067-693B-E611-AF0D-02163E011E25.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/CA0F96A1-6B3B-E611-B25C-02163E014417.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/D2ACF96A-713B-E611-9133-02163E013839.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/DCD3A6F9-6D3B-E611-BD12-02163E013802.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/E0A72A4A-763B-E611-B27D-02163E01439B.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/E2C84243-793B-E611-B492-02163E014655.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/E84C6A2A-7C3B-E611-BC92-02163E011942.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/F41A11F9-773B-E611-A30D-02163E013802.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/F66395A8-713B-E611-8E3F-02163E011DD7.root',
'/store/data/Run2016C/DoubleMuon/AOD/PromptReco-v2/000/275/657/00000/FCF024C6-773B-E611-9FDB-02163E01452A.root'

#                                '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/275/774/00000/BEDEAA86-7B3C-E611-B7A4-02163E013787.root'

#'/store/data/Run2016B/DoubleMuon/AOD/PromptReco-v2/000/273/730/00000/7EFAF1BB-9321-E611-97D5-02163E0142D8.root'
#'/store/user/jjhollar/Test2016/GammaGammaMuMu-Elastic-RunIIWinter15R-0001.root'
#                                '/store/data/Run2016B/DoubleMuon/RECO/PromptReco-v2/000/274/969/00000/FCA084F5-C832-E611-A6F2-02163E011F89.root'
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/000CB19C-6CF2-E511-8F52-001517F7F510.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/0072D43A-6DF2-E511-87F9-001E67A3EF48.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/0078F4B5-9EF2-E511-9A52-0025905A6126.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/009D2C6E-84F2-E511-B93D-0CC47A4D7644.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/0446C76C-84F2-E511-9543-0025905A6136.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/06429E88-6BF2-E511-A9C1-001E67DFF609.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/065C16D4-83F2-E511-8A1A-0CC47A78A30E.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/086E5071-6BF2-E511-BA36-001E67A3EA7A.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/08950D54-72F2-E511-AC52-0CC47A4D762A.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/08C38B56-6EF2-E511-8DEC-001E67DDD0B9.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/0E02EF6A-84F2-E511-85A3-0CC47A4C8E3C.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/10614CF9-4BF2-E511-8AA6-001E67DDC88A.root',
#'/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/20000/12002B5C-6EF2-E511-8932-001E67A3EF48.root'

    ),
#    firstEvent = cms.untracked.uint32(1)
)

#########################
#        Triggers       #
#########################
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ['HLT_Mu8_Ele17_*', 'HLT_Mu17_Ele8_*']

process.hltFilter.HLTPaths = ['HLT_DoubleMu38NoFiltersNoVtx_*']
#process.hltFilter.HLTPaths = ['HLT_Mu50_v*'] 

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START52_V9::All')
process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v8')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

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
#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
#process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offlinePrimaryVertices*_*_*',
        'keep *_*Muons*_*_*',
        'keep *_*Electrons*_*_*',
        #'keep *_*Photons*_*_*',
        'keep *_*Jets*_*_*',
        'keep *_*MET*_*_*',
        'keep recoPFCandidates_particleFlow_*_*',
        #*patEventContentNoCleaning
    ),
)

# Particle flow
from PhysicsTools.PatAlgos.tools.pfTools import *
#postfix = "PFlow"
#jetAlgo="AK4" 
#usePFBRECO(process,runPFBRECO=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix) 
#usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix)
#useGsfElectrons(process,postfix)
#removeCleaning(process)
#removeMCMatching(process, ['All'])

#########################
#       Analysis        #
#########################
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaMuMu_cfi")
process.ggmumu.TriggersList = process.hltFilter.HLTPaths
process.ggmumu.LeptonsType = cms.vstring('Muon')
#process.ggmumu.LeptonsType = cms.vstring('Electron', 'Muon')
process.ggmumu.RecoVertexLabel = cms.InputTag("offlinePrimaryVertices")
process.ggmumu.GlobalMuonCollectionLabel = cms.untracked.InputTag("muons") # RECO
#process.ggmumu.GlobalMuonCollectionLabel = cms.untracked.InputTag("selectedPatMuonsPFlow"), # PAT (particle flow)
#process.ggmumu.GlobalMuonCollectionLabel = cms.untracked.InputTag("selectedPatMuons") # PAT
process.ggmumu.GlobalEleCollectionLabel = cms.untracked.InputTag("selectedPatElectrons") # PAT
process.ggmumu.RunOnMC = cms.untracked.bool(runOnMC)
process.ggmumu.RunOnProtons = cms.untracked.bool(runOnProtons)

#process.ggmumu.outfilename = cms.untracked.string('output_mumu_dymumumc_test.root')
#process.ggmumu.outfilename = cms.untracked.string('output_mumu_data_test.root')
#process.ggmumu.outfilename = cms.untracked.string('output_mumu_mc_test.root')
process.ggmumu.outfilename = cms.untracked.string('output_doublemudata_275657_rp_test.root')


#########################
#       Protons         #
#########################
# RP geometry
#process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
#process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")


process.p = cms.Path(
 #   process.kt6PFJetsForIsolation+
#    process.pfiso+
#    process.patDefaultSequence+
#    getattr(process,"patPF2PATSequence"+postfix)+
    process.hltFilter +
#    process.recoCTPPS + 
    process.ggmumu
)
#getattr(process,"pfNoElectron"+postfix).enable = True

#print process.dumpPython()
