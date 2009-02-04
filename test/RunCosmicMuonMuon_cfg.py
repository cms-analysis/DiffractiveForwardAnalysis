import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
'/store/data/Commissioning08/MinimumBias/RECO/v1/000/067/141/00518942-BDA0-DD11-A40A-000423D95030.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.Cosmics_cfi")

process.gamgammumuanalysis.outfilename = "craftrun67818.superpointing.exclmumu.root"

# Put it all together
process.p = cms.Path(
    process.gamgammumuanalysis
    )


