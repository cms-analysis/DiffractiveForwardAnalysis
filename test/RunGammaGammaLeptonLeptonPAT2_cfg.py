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

# DB efficiency stuff
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance.db'
process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMUJPSI_WP'),
    label = cms.untracked.string('GLBMUJPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMUJPSI_TABLE'),
    label = cms.untracked.string('GLBMUJPSI_TABLE')
            )))

process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
                                       process.CondDBCommon,
                                       toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRKEFFMUJPSI_WP'),
    label = cms.untracked.string('TRKEFFMUJPSI_WP')
             ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKEFFMUJPSI_TABLE'),
    label = cms.untracked.string('TRKEFFMUJPSI_TABLE')
             )))

process.PoolDBESSource3 = cms.ESSource("PoolDBESSource",
                                       process.CondDBCommon,
                                       toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRGMUJPSI_WP'),
              label = cms.untracked.string('TRGMUJPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRGMUJPSI_TABLE'),
    label = cms.untracked.string('TRGMUJPSI_TABLE')
              )))

process.MuonPerformanceESProducer_GlobalMuon.PayloadName = "GLBMUJPSI_TABLE"
process.MuonPerformanceESProducer_GlobalMuon.WorkingPointName = "GLBMUJPSI_WP"
process.MuonPerformanceESProducer_GlobalMuon.ComponentName = "GlobalMuonFromTrackerTrackJpsi"
process.MuonPerformanceESProducer_TrackerTrackMuon.PayloadName = "TRKEFFMUJPSI_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon.WorkingPointName = "TRKEFFMUJPSI_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon.ComponentName = "TrackerTrackFromStandaloneMuonJpsi"
process.MuonPerformanceESProducer_TriggerMuon.PayloadName = "TRGMUJPSI_TABLE"
process.MuonPerformanceESProducer_TriggerMuon.WorkingPointName = "TRGMUJPSI_WP"
process.MuonPerformanceESProducer_TriggerMuon.ComponentName = "TriggerMuonFromGlobalMuonJpsi"
# End of DB stuff


# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/DiffChiB_10tev_RECO_8E29.root'
#    'file:/tmp/jjhollar/DiffChiB_10tev_RECO_8E29.root'
   'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamGamMuMu_LPAIRelastic_10tev_RECO_8E29.root'
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamP_Upsilon2Smumu_STARLIGHT_10tev_RECO_8E29.root',
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamP_Upsilon1Smumu_STARLIGHT_10tev_RECO_8E29.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP31X_V2::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load CASTOR FastSim
process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaMuMu_cfi")
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.MCGammaGammaMuMu_cfi") 
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaEE_cfi")

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring("drop *")
                               )

# PAT Layer 0+1
from Configuration.EventContent.EventContent_cff import *
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoLayer1Cleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
process.load("PhysicsTools.PatAlgos.patSequences_cff")
# replacements currently needed to make the jets work
process.allLayer1Jets.addDiscriminators    = False
process.allLayer1Jets.discriminatorSources = []
removeCleaning(process)

# Output definition - this is for testing "PATtuples"
#process.output = cms.OutputModule("PoolOutputModule",
#                                  outputCommands = cms.untracked.vstring("drop *"),    
#                                  outputCommands = cms.untracked.vstring('drop *',
#                                                                         'keep *_selectedLayer1Photons_*_*',
#                                                                         'keep *_selectedLayer1Electrons_*_*',
#                                                                         'keep *_selectedLayer1Muons_*_*',
#                                                                         'keep *_selectedLayer1Taus_*_*',
#                                                                         'keep *_selectedLayer1Jets_*_*',
#                                                                         'keep *_layer1METs_*_*',
#                                                                         'keep *_selectedLayer1PFParticles_*_*',
#                                                                         'keep *_CastorFastTowerReco_*_*',
#                                                                         *patExtraAodEventContent ),
#                                  fileName = cms.untracked.string('file:/tmp/jjhollar/mumu.pattuple.root'),
#                                  dataset = cms.untracked.PSet(
#    dataTier = cms.untracked.string('PAT'),
#    filterName = cms.untracked.string('')
#    ),
#)

#process.output.outputCommands.extend(AODEventContent.outputCommands)

process.gamgammumuanalysis.outfilename = "GamGamMuMu_LPAIRelastic_10tev_ANAL_8E29.root"

# Put it all together
process.p = cms.Path(
      process.CastorFastReco
    + process.hltFilter
    + process.patDefaultSequence  
    + process.gamgammumuanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

print process.dumpPython()
