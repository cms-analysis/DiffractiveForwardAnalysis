import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
'rfio:/castor/cern.ch/user/j/jjhollar/22XsignalMC/gamgammumu_IDEALV12_pileup2_fastsim.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V12::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# Load CASTOR FastSim
process.load("FastSimulation.ForwardDetectors.CastorTowerProducer_cfi")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaMuMu_cfi")
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATGammaGammaEE_cfi")

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

# Put it all together
process.p = cms.Path(
    process.CastorTowerReco
    + process.hltFilter
    + process.patLayer0  
    + process.patLayer1
    + process.gamgammumuanalysis
#    + process.gamgameeanalysis
    )


