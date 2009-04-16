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
        '/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V9_TrackingPointing_225-v3/0006/1A2DA3FE-7C00-DE11-B52C-001A92971B78.root',
                '/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V9_TrackingPointing_225-v3/0006/E0B7EC56-7B00-DE11-84B8-00304867903E.root',
                '/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V9_TrackingPointing_225-v3/0008/8022B2FE-A300-DE11-AD49-001731AF66AF.root',
                '/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V9_TrackingPointing_225-v3/0008/96385CED-8E00-DE11-A509-001A92971B92.root',
                '/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V9_TrackingPointing_225-v3/0015/720EDC48-1402-DE11-8925-001A9281170E.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.Cosmics_cfi")

process.gamgammumuanalysis.outfilename = "craftrun68949.superpointing.exclmumu.root"

# Put it all together
process.p = cms.Path(
    process.gamgammumuanalysis
    )


