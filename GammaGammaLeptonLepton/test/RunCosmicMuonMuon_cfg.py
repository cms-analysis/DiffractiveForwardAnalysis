import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("gamgam2leplepanalysis")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0017/969F073D-DB7E-DE11-994A-001EC9D2887E.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0015/CC3B681A-FE7D-DE11-B550-001A4BA910A0.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0015/882CC080-067E-DE11-9952-0030487CDA68.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0014/EA9ED8D5-E77D-DE11-9DD5-001A4BA566A6.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0014/6C5904DE-F97D-DE11-BFA6-0030487CDB2C.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0014/500CDC39-EE7D-DE11-AEB5-001A4BA82FF8.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0014/4AA8C161-E17D-DE11-B478-001EC9D8B179.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0013/7EAA79D2-D77D-DE11-938C-001A4BA82FF8.root',
                '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0013/0C5E722F-DC7D-DE11-AC0D-00093D121C2E.root'
        
#        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0017/F8F72E92-DB7E-DE11-86BC-0019BB3FD4CC.root',
#        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0014/1A1FB738-EC7D-DE11-BB8E-001EC9D82BC7.root'
        
#        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing_322-v1/0017/583EB59C-DB7E-DE11-BF0C-001EC9D8D987.root'
        
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# Load configuration stuff
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Load analysis modules
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.Cosmics_cfi")

process.gamgammumuanalysis.outfilename = "/tmp/jjhollar/craftrun109049.superpointing.exclmumu.root"

# Put it all together
process.p = cms.Path(
    process.gamgammumuanalysis
    )


