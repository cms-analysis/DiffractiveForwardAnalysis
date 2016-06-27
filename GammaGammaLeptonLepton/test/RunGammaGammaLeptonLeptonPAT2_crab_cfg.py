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
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
                                       SkipEvent = cms.untracked.vstring('ProductNotFound')
                                       )
# DB efficiency stuff - uncomment if running on MC
##process.load("CondCore.DBCommon.CondDBCommon_cfi")
##process.CondDBCommon.connect = 'sqlite_file:ExclJPsiZMuonChargePhysicsPerformance7TeVRun2010B.db'
##process.load("MuonAnalysis.TagAndProbe.ExclMuonPerformaceMergeESSource_cfi")
##process.load("MuonAnalysis.TagAndProbe.ExclMuonPerformaceMergeESProducer_cfi")
# End of DB stuff

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
#    				'data.root'
#    'rfio:/castor/cern.ch/user/j/jjhollar/ExclWW/GamGamWW_SM_WithFF_START44_Fall11PU_RECO.root',
#    'file:/tmp/lforthom/MuMuAnalyzer_withPU_ElElMu-15GeV_0010.root',
    'file:/tmp/lforthom/step2_RAW2DIGI_L1Reco_RECO_InelElmumu15GeV_0087.root',
                               ),
			    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_44_V13::All')
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load CASTOR FastSim
#process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaMuMu_cfi")

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
#process.hltFilter.HLTPaths = ['HLT_DoubleMu7_Acoplanarity03*']
#process.hltFilter.HLTPaths = ['HLT_DoubleMu*_Acoplanarity*','HLT_Mu13_Mu8_*','HLT_Mu17_Mu8_*']
process.hltFilter.HLTPaths = ['HLT_Mu13_Mu8_*','HLT_Mu17_Mu8_*']
process.hltFilter.throw = cms.bool(False)

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

#offline vertices with deterministic annealing. Should become the default as of 4_2_0_pre7. Requires > V01-04-04      RecoVertex/PrimaryVertexProducer
process.load("RecoVertex.Configuration.RecoVertex_cff")
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.offlinePrimaryVerticesDA = process.offlinePrimaryVertices.clone()


#process.output.outputCommands.extend(AODEventContent.outputCommands)

# Set to True if running on MC  
process.gamgammumuanalysis.ReadMCEffCorrections = False
process.gamgammumuanalysis.ReadMCEffCorrectionsByCharge = False
process.gamgammumuanalysis.ReadmcEffCorrectionsBySignedEta = False
process.gamgammumuanalysis.outfilename = "MuMuAnalyzer.root" 
#FIXME =============================================================
process.gamgammumuanalysis.maxExtraTracks = cms.untracked.int32(15)
process.gamgammumuanalysis.runningOnData = cms.untracked.bool(False)
#===================================================================

if (process.gamgammumuanalysis.runningOnData == False):
    print "Including the computation of the PU weights"
    #process.gamgammumuanalysis.mcpufile = cms.string("PileupHistos_genPU-Fall11_0110.root")
    #process.gamgammumuanalysis.mcpupath = cms.untracked.string("genPUdist/TNPUTrue")
    #process.gamgammumuanalysis.mcpupath = cms.untracked.string("genPUdist/TPU")
    #process.gamgammumuanalysis.mcpupath = cms.untracked.string("genPUdist/TNVTX")
    #process.gamgammumuanalysis.mcpufile = cms.string("PU_MC-Fall11_1110.root")
    #process.gamgammumuanalysis.mcpupath = cms.untracked.string("pu")
    #process.gamgammumuanalysis.datapufile = cms.string("PileupHistos_data-2011.root")
    #process.gamgammumuanalysis.datapufileA = cms.string("PileupHistos_data-2011A.root")
    #process.gamgammumuanalysis.datapufileB = cms.string("PileupHistos_data-2011B.root")
    #process.gamgammumuanalysis.datapufile = cms.string("PileupHistos_data-2011_HNTruth.root")
    #process.gamgammumuanalysis.datapufileA = cms.string("PileupHistos_data-2011A_HNTruth.root")
    #process.gamgammumuanalysis.datapufileB = cms.string("PileupHistos_data-2011B_HNTruth.root")
    #process.gamgammumuanalysis.datapufile = cms.string("PileupHistos_data-2011_HNnonTruth.root")
    #process.gamgammumuanalysis.datapufileA = cms.string("PileupHistos_data-2011A_HNnonTruth.root")
    #process.gamgammumuanalysis.datapufileB = cms.string("PileupHistos_data-2011B_HNnonTruth.root")
    #process.gamgammumuanalysis.datapufile = cms.string("MyDataPileupHistogram.root")
    #process.gamgammumuanalysis.datapufileA = cms.string("MyDataPileupHistogramA.root")
    #process.gamgammumuanalysis.datapufileB = cms.string("MyDataPileupHistogramB.root")
    process.gamgammumuanalysis.mcpufile = cms.string("PU_MC-Fall11_1710.root")
    process.gamgammumuanalysis.mcpupath = cms.untracked.string("pileup")
    process.gamgammumuanalysis.datapufile = cms.string("Cert_160404-178078_7TeV_Collisions11_JSON_v3.pileupTruth_v2_finebin.root")
    process.gamgammumuanalysis.datapufileA = cms.string("Cert_160404-173692_7TeV_Collisions11_JSON_v3.pileupTruth_v2_finebin.root")
    process.gamgammumuanalysis.datapufileB = cms.string("Cert_175832-178078_7TeV_Collisions11_JSON_v3.pileupTruth_v2_finebin.root")
    process.gamgammumuanalysis.datapupath = cms.untracked.string("pileup")

else:
    process.gamgammumuanalysis.mcpufile = cms.string("")
    process.gamgammumuanalysis.mcpupath = cms.untracked.string("")
    process.gamgammumuanalysis.datapufile = cms.string("")
    process.gamgammumuanalysis.datapufileA = cms.string("")
    process.gamgammumuanalysis.datapufileB = cms.string("")
    process.gamgammumuanalysis.datapupath = cms.untracked.string("")
    


																						
# Put it all together
process.p = cms.Path(
#    process.mcgamgammumuanalysis
#    process.CastorFastReco
    process.hltFilter 
    + process.offlinePrimaryVerticesDA
    + process.patDefaultSequence  
    + process.gamgammumuanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

#print process.dumpPython()
