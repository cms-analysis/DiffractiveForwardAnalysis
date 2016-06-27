import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('PATLayer0Summary')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# DB efficiency stuff - uncomment if running on MC
##process.load("CondCore.DBCommon.CondDBCommon_cfi")
##process.CondDBCommon.connect = 'sqlite_file:ExclJPsiZMuonChargePhysicsPerformance7TeVRun2010B.db'
##process.load("MuonAnalysis.TagAndProbe.ExclMuonPerformaceMergeESSource_cfi")
##process.load("MuonAnalysis.TagAndProbe.ExclMuonPerformaceMergeESProducer_cfi")
# End of DB stuff

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
#                                '/store/relval/CMSSW_4_4_0/RelValZMM/GEN-SIM-RECO/START44_V5-v2/0051/0407CEA7-ECE6-E011-B360-002618943826.root',
#                                '/store/relval/CMSSW_4_4_0/RelValZMM/GEN-SIM-RECO/START44_V5-v2/0045/E095BF22-11E6-E011-8B1F-002618943821.root',
#                                '/store/relval/CMSSW_4_4_0/RelValZMM/GEN-SIM-RECO/START44_V5-v2/0045/82EECF40-10E6-E011-A6D9-0018F3D09624.root',
#                                '/store/relval/CMSSW_4_4_0/RelValZMM/GEN-SIM-RECO/START44_V5-v2/0045/7C73A5A1-0AE6-E011-992F-00261894380A.root'

##                                'rfio:/castor/cern.ch/user/j/jjhollar/ExclWW/GamGamWW_Anomalous1_Lambda500GeV_A0W0point0002_NoFF_START44_Fall11PU_RAW.root'
    'rfio:/castor/cern.ch/user/j/jjhollar/ExclWW/GamGamWW_SM_WithFF_START44_Fall11PU_RECO.root'
##                                'file:/tmp/jjhollar/CalcHEP_GamGamWW_ElEl_SM_RECO_10kevents.root'
                                                                )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_44_V12::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load CASTOR FastSim
#process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaMuE_cfi")

# require scraping filter
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
                                src =cms.InputTag("muons"), minNumber = cms.uint32(1))


# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_Mu10_Ele10_CaloIdL_*', 'HLT_Mu8_Ele17_*', 'HLT_Mu17_Ele8_*']

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring("drop *")
                               )

# PAT Layer 0+1
from Configuration.EventContent.EventContent_cff import *
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("PhysicsTools.PatAlgos.patSequences_cff")
removeCleaning(process)
removeMCMatching(process, ['All'])
#restrictInputToAOD(process, ['All'])

# JH - testing
#from PhysicsTools.PatAlgos.tools.pfTools import *
#postfix = "PFlow"
#jetAlgo="AK5"
#usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix)
# end JH


# Set to True if running on MC  
##process.gamgammueanalysis.ReadMCEffCorrections = True
##process.gamgammueanalysis.outfilename = "SignalMuEAnalyzer_Anomalous1_Lambda500GeV_A0W0point0002_NoFF_START44_Fall11PU.root"
process.gamgammueanalysis.outfilename = "/tmp/jjhollar/MuEAnalyzer_nonPF.root"

# Put it all together
process.p = cms.Path(
#    process.hltFilter 
    process.patDefaultSequence  
#    getattr(process,"patPF2PATSequence"+postfix) 
    + process.gamgammueanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

#print process.dumpPython()
