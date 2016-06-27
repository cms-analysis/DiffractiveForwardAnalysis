import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
# Field-on 900 GeV MinBias
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8D_900GeV-v1/0005/E4B3A7BE-3AD7-DE11-9230-002618943939.root'

# 11 golden events (BSC coincidence) , bad tracking
#    'file:skimI/bit40or41skim.root'

# 68 events from 122314 (BSC coincidence), reReco
#    'file:skimI/reRecoOutput_122314.root'

# 142 events from run 122314 (BSc trigger + tracker activity) , bad tracking
#        'file:skimII/BSC_activity_58.root',
#        'file:skimII/BSC_activity_61.root',
#        'file:skimII/BSC_activity_62.root',
#        'file:skimII/BSC_activity_63.root',
#        'file:skimII/BSC_activity_64.root'

# 142 events from run 122314 (BSc trigger + tracker activity) , GOOD tracking (RETRACK & REVERTEX)
#	'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/skimII/Run122314_BSCSkim_MinBiasPD_ReTracking.root'

# 142 events from run 122314 (should be same)
#        'file:skimII_v6/BSC_activity_57.root',
#        'file:skimII_v6/BSC_activity_59.root',
#        'file:skimII_v6/BSC_activity_60.root',
#        'file:skimII_v6/BSC_activity_62.root',
#        'file:skimII_v6/BSC_activity_63.root',
#        'file:skimII_v6/BSC_activity_64.root'

# 142 events ReReco by me
#	'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/skimII_v6/Skim_PersoReco.root'
#       'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/skimII_v6/Skim_PersoReco2.root'

# 33810 events from run 122294 (BSC trigger + tracker activity),bad tracking
#        'file:skimII/BSC_activity_1.root',
#        'file:skimII/BSC_activity_2.root',
#        'file:skimII/BSC_activity_10.root',
#        'file:skimII/BSC_activity_11.root',
#        'file:skimII/BSC_activity_12.root',
#        'file:skimII/BSC_activity_13.root',
#        'file:skimII/BSC_activity_34.root',
#        'file:skimII/BSC_activity_35.root',
#        'file:skimII/BSC_activity_36.root',
#        'file:skimII/BSC_activity_50.root',
#        'file:skimII/BSC_activity_51.root',
#        'file:skimII/BSC_activity_52.root'

# MONTECARLO MINBIAS, 0T, TOB-track only
	#'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_1.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_2.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_3.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_4.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_5.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_6.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_7.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_8.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_9.root',
         #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_10.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_11.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_12.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_13.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_14.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_15.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_16.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_17.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_18.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_19.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_20.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_21.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_22.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_23.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_24.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_25.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_26.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_27.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_28.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_29.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_30.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_31.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_32.root',
       #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_33.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_34.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_35.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_36.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_37.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_38.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_39.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_40.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_41.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_42.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_43.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_44.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_45.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_46.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_47.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_48.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_49.root',
        #'file:/home/fynu/schul/scratch/data_analyses/CMSSW_3_3_3/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/mc/MinBias900GeV_NoField_v2_TOBONLY_50.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25000))


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('FIRSTCOLL::All')
process.load("Configuration.StandardSequences.MagneticField_0T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.CollisionsMuMu_cfi")

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring("drop *")
                               )

#process.load('DPGAnalysis/SiStripTools/largesistripclusterevents_cfi')
#process.largeSiStripClusterEvents.absoluteThreshold = 50
#process.largeSiStripClusterEvents.moduleThreshold = 1000
#process.clusters = cms.Path(process.largeSiStripClusterEvents)
#process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")

process.gamgammumuanalysis.outfilename = "minbias.900gev.startupv8d.root"

# Put it all together
process.p = cms.Path(
#    process.hltFilter +
#    process.largeSiStripClusterEvents +
    process.gamgammumuanalysis
    )
#print process.dumpPython()
