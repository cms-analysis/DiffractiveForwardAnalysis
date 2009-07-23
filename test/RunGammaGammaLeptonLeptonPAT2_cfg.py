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

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
    'file:/tmp/jjhollar/GamGamMuMu_LPAIRelastic_10tev_RECO_8E29.root'
#        '/store/relval/CMSSW_3_1_1/RelValZMM/GEN-SIM-RECO/STARTUP31X_V1-v2/0002/C8CEE598-CB6B-DE11-871F-001D09F2905B.root',
#        '/store/relval/CMSSW_3_1_1/RelValZMM/GEN-SIM-RECO/STARTUP31X_V1-v2/0002/642F8176-C96B-DE11-9D10-000423D98BE8.root',
#        '/store/relval/CMSSW_3_1_1/RelValZMM/GEN-SIM-RECO/STARTUP31X_V1-v2/0002/46AA6A11-D46B-DE11-A614-001D09F25438.root',
#        '/store/relval/CMSSW_3_1_1/RelValZMM/GEN-SIM-RECO/STARTUP31X_V1-v2/0002/443BC1DD-CC6B-DE11-804C-000423D98EC4.root'        
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
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaEE_cfi")

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")

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
#                                  )

#process.output.outputCommands.extend(AODEventContent.outputCommands)

# Put it all together
process.p = cms.Path(
      process.CastorFastReco
    + process.hltFilter
    + process.patDefaultSequence  
    + process.gamgammumuanalysis
#   The output module here is only needed for making 'PATtuples'/skims
#   + process.output
    )

