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

# End of DB stuff


# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
    '/store/data/BeamCommissioning09/MinimumBias/RECO/Jan29ReReco-v2/0020/E8593279-0A0E-DF11-A36D-001A9281171E.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP31X_V2::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load CASTOR FastSim
process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATExclusiveTrackTrack_cfi")
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.MCGammaGammaVV_cfi") 
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaEE_cfi")

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")

process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ['HLT_MinBiasPixel_SingleTrack']
process.hltFilter.HLTPaths = ['HLT_PhysicsDeclared']

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

#process.output.outputCommands.extend(AODEventContent.outputCommands)

process.excltrktrkanalysis.outfilename = "ExclTrackTrack_900GeVJan29rerecoMC.root"

# Put it all together
process.p = cms.Path(
    process.hltFilter
    + process.CastorFastReco
    #    process.patDefaultSequence  
    + process.excltrktrkanalysis
    )

#print process.dumpPython()
