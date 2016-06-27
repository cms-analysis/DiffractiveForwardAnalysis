import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("ANALYSIS")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
                                '/store/data/Run2010B/MinimumBias/RECO/PromptReco-v2/000/148/952/FC01E3BC-8AE1-DF11-8874-0030487CD7CA.root'
                                )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load CASTOR FastSim
##process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.ZeroBiasAnalyzer_cfi")

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")

process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_ZeroBias']

#process.out = cms.OutputModule("PoolOutputModule",
#                               outputCommands = cms.untracked.vstring("drop *")
#                               )

# PAT Layer 0+1
##process.excltrktrkanalysis.outfilename = "/tmp/jjhollar/ExclTrackTrack_900GeVJan29rerecoMC.root"

# Put it all together
process.p = cms.Path(
    process.hltFilter +
    process.zerobiasanalysis
    )

