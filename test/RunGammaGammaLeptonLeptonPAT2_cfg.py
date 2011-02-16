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
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/EEE64447-9F9B-DF11-A573-001617C3B69C.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/E04E9705-8A9B-DF11-8EDB-003048D37560.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/CA804E38-8E9B-DF11-A825-001D09F252DA.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/C4A1CB31-9D9B-DF11-805E-003048F1BF66.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/BA38B11B-C39B-DF11-B826-001D09F25208.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/8AD5A5DC-969B-DF11-A088-001D09F2441B.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/60D7391C-9B9B-DF11-A1E6-001D09F231C9.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/445A9EF8-989B-DF11-BC84-0030487CD7C0.root',
                             '/store/data/Run2010A/MuOnia/RECO/v4/000/141/956/00BE9355-909B-DF11-A111-001617C3B69C.root'
                                )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V9::All')
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
process.hltFilter.HLTPaths = ['HLT_DoubleMu3']

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

# Output definition - this is for testing "PATtuples"
#process.output = cms.OutputModule("PoolOutputModule",
#                                  outputCommands = cms.untracked.vstring("drop *"),    
#                                  outputCommands = cms.untracked.vstring('drop *',
#                                                                         'keep *_selectedPatPhotons_*_*',
#                                                                         'keep *_selectedPatElectrons_*_*',
#                                                                         'keep *_selectedPatMuons_*_*',
#                                                                         'keep *_selectedPatTaus_*_*',
#                                                                         'keep *_selectedPatJets_*_*',
#                                                                         'keep *_layer1METs_*_*',
#                                                                         'keep *_selectedPatPFParticles_*_*',
#                                                                         'keep *_CastorFastTowerReco_*_*',
#                                                                         *patExtraAodEventContent ),
#                                  fileName = cms.untracked.string('file:/tmp/jjhollar/mumu.pattuple.root'),
#                                  dataset = cms.untracked.PSet(
#    dataTier = cms.untracked.string('PAT'),
#    filterName = cms.untracked.string('')
#    ),
#)

#process.output.outputCommands.extend(AODEventContent.outputCommands)

# Set to True if running on MC  
##process.gamgammumuanalysis.ReadMCEffCorrections = True
process.gamgammumuanalysis.outfilename = "DimuonAnalyzer.root"

# Put it all together
process.p = cms.Path(
#    process.mcgamgammumuanalysis
#    process.CastorFastReco
    process.hltFilter 
    + process.patDefaultSequence  
    + process.gamgammumuanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

#print process.dumpPython()
