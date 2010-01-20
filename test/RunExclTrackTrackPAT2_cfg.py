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
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/FE69944A-D2EC-DE11-83AD-00163EC20501.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/FACDCBE9-C3EC-DE11-96F6-002481DE4C92.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F898E834-A7EC-DE11-9D06-001EC9ED4FAA.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F88F177A-C8EC-DE11-A974-00187186E871.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F6A0ED38-CAEC-DE11-984D-002481DE47D0.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F4EA1131-BBEC-DE11-B04F-00163E051101.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F4E4C00D-B8EC-DE11-AFFE-0025B3E025B6.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F4CE101C-CAEC-DE11-8CC3-00187186EC17.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F49488C8-C1EC-DE11-B444-001EC9E14A75.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F417F642-C5EC-DE11-9BD0-002481CFEA5C.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F2C66841-C5EC-DE11-AC9C-002481DE4B5A.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F2B7E6D5-C1EC-DE11-AF93-002481DE4CC2.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F2A8AC05-D6EC-DE11-AE47-002481CFEA0E.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F209BC2D-A7EC-DE11-B9EE-001EC9E1671C.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F0F8DDD6-CBEC-DE11-A268-001A6478AB0C.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/F0DA1677-ACEC-DE11-9B23-00163EC20401.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/EE55122B-B1EC-DE11-8A54-001EC9E12EDB.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/ECF3AF84-D2EC-DE11-8650-00187186E6F6.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/EC96F440-C5EC-DE11-8825-002481CFE648.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E860E2ED-C3EC-DE11-B544-002481CFE864.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E64F7C62-D2EC-DE11-89F4-00163E041101.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E63AC380-AFEC-DE11-B210-001EC9ED4FB2.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E4885230-ACEC-DE11-993F-002481DE4B5A.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E43DF281-B5EC-DE11-A903-00163EC10301.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E0C562E5-A8EC-DE11-BFA2-001EC9ED88D4.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E06481F7-D3EC-DE11-A146-00163E031401.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E035AF9C-D5EC-DE11-AF12-00163E050901.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E0143EC5-C1EC-DE11-BFD5-001EC9ED4F72.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/E004D30D-B8EC-DE11-A5E1-002481CFE888.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/DEBBC482-B0EC-DE11-8034-00163E060901.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/DC16C566-D2EC-DE11-8757-00163EC20501.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/DAAF3943-CCEC-DE11-BDB5-002481CFE834.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/DA4D95E8-C3EC-DE11-B70C-0025B3268576.root',
    '/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8K_900GeV-v1/0019/D6442764-AAEC-DE11-96F1-001EC9ED7E42.root'
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
process.hltFilter.HLTPaths = ['HLT_MinBiasPixel_SingleTrack']

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

process.excltrktrkanalysis.outfilename = "/tmp/jjhollar/ExclTrackTrack_900GeVMC.root"

# Put it all together
process.p = cms.Path(
    #    process.hltFilter
    process.patDefaultSequence  
    + process.excltrktrkanalysis
    )

#print process.dumpPython()
