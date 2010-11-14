import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
                                'rfio:/castor/cern.ch/user/j/jjhollar/ExclPiPi/STARLIGHT_ExclusiveRhoToPiPi_RAW2DIGI_L1Reco_RECO_VALIDATION.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATExclusiveTrackTrack_cfi")

# If running on MC, Llad CASTOR FastSim
#process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")
#process.excltrktrkanalysis.CastorTowerLabel = "CastorFastTowerReco"

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
#process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_ZeroBiasPixel_SingleTrack','HLT_ZeroBias']

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring("drop *")
                               )

process.excltrktrkanalysis.outfilename = "TestExclusivePiPiAnalyzer.root"

# Put it all together
process.p = cms.Path(
    #    process.hltFilter
    #   process.CastorFastReco
    process.excltrktrkanalysis
    )

#print process.dumpPython()
