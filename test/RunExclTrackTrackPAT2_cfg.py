import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
            '/store/data/Run2010A/ZeroBias/RECO/May27thReReco_v1/0174/FE0369EF-956A-DF11-ADD6-00E0817917F5.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR09_R_35X_V2::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# If running on MC, Llad CASTOR FastSim
#process.load("FastSimulation.ForwardDetectors.CastorFastReco_cff")
#process.excltrktrkanalysis.CastorTowerLabel = "CastorFastTowerReco"

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.PATExclusiveTrackTrack_cfi")

# Trigger
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
#process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_ZeroBiasPixel_SingleTrack','HLT_ZeroBias']

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring("drop *")
                               )

process.excltrktrkanalysis.outfilename = "ExclTrackTrack_June9_Commissioning10.root"

# Put it all together
process.p = cms.Path(
    process.hltFilter
    #    process.CastorFastReco
    + process.excltrktrkanalysis
    )

#print process.dumpPython()
