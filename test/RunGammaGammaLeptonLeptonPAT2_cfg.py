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
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance.db'
#process.CondDBCommon.connect = 'sqlite_file:../../../MuonAnalysis/TagAndProbe/test/performanceDB/MuonPhysicsPerformance.db'
#process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")

##process.PoolDBESSource = cms.ESSource("PoolDBESSource",
##                                      process.CondDBCommon,
##                                      toGet = cms.VPSet(
##    cms.PSet(
##    record = cms.string('PerformanceWPRecord'),
##    tag = cms.string('GLBMUJPSI_OCTXTEST_WP'),
##    label = cms.untracked.string('GLBMUJPSI_OCTXTEST_WP')
##    ),
##    cms.PSet(
##    record = cms.string('PerformancePayloadRecord'),
##    tag = cms.string('GLBMUJPSI_OCTXTEST_TABLE'),
##    label = cms.untracked.string('GLBMUJPSI_OCTXTEST_TABLE')
##    ),
##    cms.PSet(
##    record = cms.string('PerformanceWPRecord'),
##    tag = cms.string('GLBMUZ_OCTXTEST_WP'),
##    label = cms.untracked.string('GLBMUZ_OCTXTEST_WP')
##    ),
##    cms.PSet(
##    record = cms.string('PerformancePayloadRecord'),
##    tag = cms.string('GLBMUZ_OCTXTEST_TABLE'),
##    label = cms.untracked.string('GLBMUZ_OCTXTEST_TABLE')
##    )))

##process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
##                                       process.CondDBCommon,
##                                       toGet = cms.VPSet(
##    cms.PSet(
##    record = cms.string('PerformanceWPRecord'),
##    tag = cms.string('TRKEFFMUJPSI_OCTXTEST_WP'),
##    label = cms.untracked.string('TRKEFFMUJPSI_OCTXTEST_WP')
##    ),
##    cms.PSet(
##    record = cms.string('PerformancePayloadRecord'),
##    tag = cms.string('TRKEFFMUJPSI_OCTXTEST_TABLE'),
##    label = cms.untracked.string('TRKEFFMUJPSI_OCTXTEST_TABLE')
##    ),
##    cms.PSet(
##    record = cms.string('PerformanceWPRecord'),
##    tag = cms.string('TRKEFFMUZ_OCTXTEST_WP'),
##    label = cms.untracked.string('TRKEFFMUZ_OCTXTEST_WP')
##    ),
##    cms.PSet(
##    record = cms.string('PerformancePayloadRecord'),
##            tag = cms.string('TRKEFFMUZ_OCTXTEST_TABLE'),
##    label = cms.untracked.string('TRKEFFMUZ_OCTXTEST_TABLE')
##    )))

##process.PoolDBESSource3 = cms.ESSource("PoolDBESSource",
##                                       process.CondDBCommon,
##                                       toGet = cms.VPSet(
##    cms.PSet(
##    record = cms.string('PerformanceWPRecord'),
##    tag = cms.string('TRGMUJPSI_OCTXTEST_WP'),
##    label = cms.untracked.string('TRGMUJPSI_OCTXTEST_WP')
##    ),
##    cms.PSet(
##    record = cms.string('PerformancePayloadRecord'),
##    tag = cms.string('TRGMUJPSI_OCTXTEST_TABLE'),
##    label = cms.untracked.string('TRGMUJPSI_OCTXTEST_TABLE')
##    ),
##    cms.PSet(
##    record = cms.string('PerformanceWPRecord'),
##    tag = cms.string('TRGMUZ_OCTXTEST_WP'),
##    label = cms.untracked.string('TRGMUZ_OCTXTEST_WP')
##    ),
##    cms.PSet(
##    record = cms.string('PerformancePayloadRecord'),
##    tag = cms.string('TRGMUZ_OCTXTEST_TABLE'),
##    label = cms.untracked.string('TRGMUZ_OCTXTEST_TABLE')
##    )))

#
# change inside the source
#
##process.MuonPerformanceESProducer_GlobalMuon1.PayloadName = "GLBMUZ_OCTXTEST_TABLE"
##process.MuonPerformanceESProducer_GlobalMuon1.WorkingPointName = "GLBMUZ_OCTXTEST_WP"
##process.MuonPerformanceESProducer_GlobalMuon1.ComponentName = "GlobalMuonFromTrackerTrackZ"
##process.MuonPerformanceESProducer_TrackerTrackMuon1.PayloadName = "TRKEFFMUZ_OCTXTEST_TABLE"
##process.MuonPerformanceESProducer_TrackerTrackMuon1.WorkingPointName = "TRKEFFMUZ_OCTXTEST_WP"
##process.MuonPerformanceESProducer_TrackerTrackMuon1.ComponentName = "TrackerTrackFromStandaloneMuonZ"
##process.MuonPerformanceESProducer_TriggerMuon1.PayloadName = "TRGMUZ_OCTXTEST_TABLE"
##process.MuonPerformanceESProducer_TriggerMuon1.WorkingPointName = "TRGMUZ_OCTXTEST_WP"
##process.MuonPerformanceESProducer_TriggerMuon1.ComponentName = "TriggerMuonFromGlobalMuonZ"

##process.MuonPerformanceESProducer_GlobalMuon2.PayloadName = "GLBMUJPSI_OCTXTEST_TABLE"
##process.MuonPerformanceESProducer_GlobalMuon2.WorkingPointName = "GLBMUJPSI_OCTXTEST_WP"
##process.MuonPerformanceESProducer_GlobalMuon2.ComponentName = "GlobalMuonFromTrackerTrackJpsi"
##process.MuonPerformanceESProducer_TrackerTrackMuon2.PayloadName = "TRKEFFMUJPSI_OCTXTEST_TABLE"
##process.MuonPerformanceESProducer_TrackerTrackMuon2.WorkingPointName = "TRKEFFMUJPSI_OCTXTEST_WP"
##process.MuonPerformanceESProducer_TrackerTrackMuon2.ComponentName = "TrackerTrackFromStandaloneMuonJpsi"
##process.MuonPerformanceESProducer_TriggerMuon2.PayloadName = "TRGMUJPSI_OCTXTEST_TABLE"
##process.MuonPerformanceESProducer_TriggerMuon2.WorkingPointName = "TRGMUJPSI_OCTXTEST_WP"
##process.MuonPerformanceESProducer_TriggerMuon2.ComponentName = "TriggerMuonFromGlobalMuonJpsi"
 
# End of DB stuff


# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
#    'file:/tmp/schul/ups1S_mumu_SIM-RECO_7TeV_8E29_startup-mc.root'
    'file:/tmp/schul/chib_SIM-RECO_7TeV_8E29_startup-mc.root'
#        '/store/mc/Summer09/ppMuMuX/GEN-SIM-RECO/MC_31X_V3_SD_DoubleMu-v1/0003/2056F3A0-79AA-DE11-988C-003048D462C8.root'        
#    'file:/tmp/schul/ups1S_mumu_SIM-RECO_7TeV_8E29_startup-mc.root'
#        '/store/mc/Summer09/ppMuMuX/GEN-SIM-RECO/MC_31X_V3-v1/0007/DA87304D-2A8E-DE11-99F3-001AA009FB99.root'        
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/DiffChiB_10tev_RECO_8E29.root'
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamP_Upsilon1Smumu_STARLIGHT_10tev_RECO_8E29.root'
#    'file:/tmp/jjhollar/inelastictest_GEN.root'
#    'file:/tmp/schul/inelgamma_mumu_SIM-RECO_0PU_8E29_startup-mc.root'
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/DiffChiB_10tev_RECO_8E29.root'
#    'file:/tmp/jjhollar/DiffChiB_10tev_RECO_8E29.root'
#   'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamGamMuMu_LPAIRelastic_10tev_RECO_8E29.root'
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamP_Upsilon2Smumu_STARLIGHT_10tev_RECO_8E29.root',
#    'rfio:/castor/cern.ch/user/j/jjhollar/312signalMC/GamP_Upsilon1Smumu_STARLIGHT_10tev_RECO_8E29.root'
#'rfio:/castor/cern.ch/user/r/roland/FAMC10TeV/SmallSample/PYTHIA6_DYmumu_M5_20_10TeV/PYTHIA6_DYmumu_M5_20_10TeV_cff_py_GEN_FASTSIM.root'
    
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
#from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")
# replacements currently needed to make the jets work
process.allLayer1Jets.addDiscriminators    = False
process.allLayer1Jets.discriminatorSources = []
removeCleaning(process)
#restrictInputToAOD(process, ['All'])

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

process.gamgammumuanalysis.outfilename = "Chib_7TeV_8E29_withmumug.root"

# Put it all together
process.p = cms.Path(
#    process.mcgamgammumuanalysis
      process.CastorFastReco
    + process.hltFilter
    + process.patDefaultSequence  
    + process.gamgammumuanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

#print process.dumpPython()
