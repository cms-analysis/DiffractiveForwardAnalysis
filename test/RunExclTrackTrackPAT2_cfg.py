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
    #    '/store/relval/CMSSW_3_5_2/RelValMinBias/GEN-SIM-RECO/START3X_V21-v1/0016/1880B254-D91E-DF11-B071-003048678B1C.root'
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/BEF77CC5-C71D-DF11-A707-00237DA15C7C.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/B6B3BCD7-C61D-DF11-9201-00237DA1CD92.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/A81017F4-C91D-DF11-BFBF-00237DA12CEC.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/A40AE4C3-D21D-DF11-B94F-001B78E1096E.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/68DC12EB-CD1D-DF11-94F8-00237DA1DDE4.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/44E91AE1-C81D-DF11-83D5-00237DA13C2E.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/4254A372-C81D-DF11-B090-00237DA12CA0.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/2239603C-C31D-DF11-A5D7-0022640D5E2E.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/18thFebPreProd_351p1-v3/0000/10B3EBA3-CC1D-DF11-83DD-00237DA13C16.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR09_R_35X_V2::All')
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

process.excltrktrkanalysis.outfilename = "/tmp/jjhollar/ExclTrackTrack_900GeVJan29rerecoMC.root"

# Put it all together
process.p = cms.Path(
    process.hltFilter
    #    + process.CastorFastReco
    #    process.patDefaultSequence  
    + process.excltrktrkanalysis
    )

#print process.dumpPython()
