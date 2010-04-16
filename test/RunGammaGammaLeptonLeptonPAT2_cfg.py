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
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F4C92A98-163C-DF11-9788-0030487C7392.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F427D642-173C-DF11-A909-0030487C60AE.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/E27821C3-0C3C-DF11-9BD9-0030487CD718.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/D87D5469-2E3C-DF11-A470-000423D99896.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/B647CAD9-0E3C-DF11-886F-0030487CD716.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/A860D55E-193C-DF11-BE29-0030487C60AE.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/9884BC11-0C3C-DF11-8F9C-000423D986C4.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/92684831-233C-DF11-ABA0-0030487CD16E.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/90269E76-0D3C-DF11-A1A0-0030487CD840.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/8CAE3014-133C-DF11-A05D-000423D174FE.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/8C51BAC6-1A3C-DF11-A0EE-000423D94A04.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/8C042B04-2D3C-DF11-939F-0030487CD178.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/80471A6B-0E3C-DF11-8DCD-0030487C6A66.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/762824C3-0C3C-DF11-A4FD-0030487CD6D2.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/6A3533F5-103C-DF11-B3AA-00304879BAB2.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/4C8979D2-073C-DF11-B97B-000423D6AF24.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/26C8DED9-0E3C-DF11-9D83-0030487CD7B4.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/181C44F7-093C-DF11-A9CB-001D09F24FEC.root',
                                '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/0AA7C390-0F3C-DF11-BD65-000423D998BA.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START3X_V21::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load CASTOR FastSim
#process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaMuMu_cfi")
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.MCGammaGammaMuMu_cfi") 
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaEE_cfi")

# Physics declared and scraping removal
# require physics declared
process.physDecl = cms.EDFilter("PhysDecl",
                                applyfilter = cms.untracked.bool(True)
                                )

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

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring("drop *")
                               )

# PAT Layer 0+1
from Configuration.EventContent.EventContent_cff import *
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
#from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
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

process.gamgammumuanalysis.outfilename = "run132440torun132606_MinimumBiasPromptReco.root"

# Put it all together
process.p = cms.Path(
#    process.mcgamgammumuanalysis
#      process.CastorFastReco
    process.hltFilter*
    process.scrapingVeto*
    process.physDecl
#    process.muonFilter
    + process.patDefaultSequence  
    + process.gamgammumuanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

#print process.dumpPython()
