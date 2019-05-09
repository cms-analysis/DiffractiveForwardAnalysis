import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

runOnMC = False
useAOD = True # AOD or MiniAOD?

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True),
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/data/Run2016G/DoubleEG/AOD/23Sep2016-v1/100000/0042DBD3-BA8E-E611-919E-002481ACDAA8.root',
#'/store/data/Run2017C/DoubleMuon/AOD/12Sep2017-v1/10000/029F251F-B1A2-E711-AAC3-001E67792890.root',
#'/store/data/Run2018B/DoubleMuon/MINIAOD/17Sep2018-v1/00000/8E1342C6-AA35-9049-B101-B5B595EAAEE2.root'
#'/store/data/Run2017C/DoubleMuon/AOD/17Nov2017-v1/30001/30DDB6DA-CBD8-E711-AF2A-A4BF0112BCB4.root'
#'file:/tmp/jjhollar/8816F63B-C0D5-E711-B32B-002590D9D9F0.root'

#'/store/data/Run2017D/DoubleMuon/AOD/17Nov2017-v1/30000/6482D69E-48D6-E711-9AF7-008CFAF71FB4.root'
#'file:/tmp/jjhollar/6482D69E-48D6-E711-9AF7-008CFAF71FB4.root'
#'file:/tmp/jjhollar/F2C53386-ABFA-4A4F-B851-0065858DB53C.root'
#'file:/tmp/jjhollar/14D52021-5DDE-E711-8A48-02163E01453B.root'
#
# at CERN eos
#'/store/data/Run2017C/DoubleMuon/AOD/PromptReco-v3/000/301/998/00000/8E0BF23F-8F8E-E711-9774-02163E019B27.root'
#
#'/store/data/Run2017C/DoubleMuon/AOD/17Nov2017-v1/30000/90A083BD-CBD8-E711-A440-A4BF0108B5F2.root'
#'file:/tmp/jjhollar/90A083BD-CBD8-E711-A440-A4BF0108B5F2.root'
'file:pickevents.root'
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = cms.vstring(
    'HLT_DoubleMu43NoFiltersNoVtx_*',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_*',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_*'
#    'HLT_DoubleEle33_CaloIdL_MW_v*',
#    'HLT_Ele27_HighEta_Ele20_Mass55_v*',
#    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*',
)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATuple.root'),
    outputCommands = cms.untracked.vstring(
#        'drop *',
#        'keep *_offline*PrimaryVertices*_*_*',
#        'keep *_selectedPatMuons*_*_*',
#        'keep *_*lectron*_*_*',
#        'keep *_selectedPatElectrons*_*_*',
#        'keep *_selectedPat*Photons*_*_*',
#        'keep *_selectedPatJets*_*_*',
#        'keep *_*MET*_*_*',
#        'keep *_*particleFlow*_*_*',
        'keep *_*_*_*'
    ),
)
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
patAlgosToolsTask.add(process.patCandidatesTask)

process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
patAlgosToolsTask.add(process.selectedPatCandidatesTask)

from PhysicsTools.PatAlgos.tools.coreTools import runOnData
if not runOnMC:
    runOnData( process )

#########################
#      Electron ID      #
#########################

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, setupVIDElectronSelection, setupAllVIDIdsInModule, DataFormat
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDElectronIdProducer(process, DataFormat.AOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff', setupVIDElectronSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff', setupVIDElectronSelection)

#########################
#       Photon ID       #
#########################

switchOnVIDPhotonIdProducer(process, DataFormat.AOD)
setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff', setupVIDPhotonSelection)

#########################
#     Proton RECO       #
#########################
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

# get optics from Wagners's DB tag
from CondCore.CondDB.CondDB_cfi import *
process.CondDBOptics = CondDB.clone( connect = 'frontier://FrontierProd/CMS_CONDITIONS' )
process.PoolDBESSourceOptics = cms.ESSource("PoolDBESSource",
    process.CondDBOptics,
    DumpStat = cms.untracked.bool(False),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('CTPPSOpticsRcd'),
        tag = cms.string("PPSOpticalFunctions_offline_v1")
    )),
)

# get alignment from Clemencia's file
from CondCore.CondDB.CondDB_cfi import *
process.CondDBAlignment = CondDB.clone( connect = 'sqlite_file:/afs/cern.ch/user/c/cmora/public/CTPPSDB/AlignmentSQlite/CTPPSRPRealAlignment_table_v26Apr.db' )
process.PoolDBESSourceAlignment = cms.ESSource("PoolDBESSource",
    process.CondDBAlignment,
    #timetype = cms.untracked.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('RPRealAlignmentRecord'),
        tag = cms.string('CTPPSRPAlignment_real_table_v26A19')
    ))
)

# get LHCInfo from DB tag
#from CondCore.CondDB.CondDB_cfi import *
#CondDB.connect = 'frontier://FrontierProd/CMS_CONDITIONS'
#process.PoolDBESSource = cms.ESSource("PoolDBESSource",
#    CondDB,
#    DumpStat = cms.untracked.bool(False),
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('LHCInfoRcd'),
#        tag = cms.string("LHCInfoEndFill_prompt_v2")
#    ))
#)


#JH - ESPrefer to get optical functions from CTPPSOpticalFunctionsESSource instead of global tag for now                                                                           
#process.es_prefer_ppsOptics = cms.ESPrefer("CTPPSOpticalFunctionsESSource","ctppsOpticalFunctionsESSource")
process.es_prefer_ppsOptics = cms.ESPrefer("PoolDBESSource","PoolDBESSourceOptics")

# For testing on old prompt reco data - need to rerun pixel tracking + LiteTracks first
# Should not be needed for final version running on re-RECO data
process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","ggll")
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(True)
process.ctppsProtons.doSingleRPReconstruction = cms.bool(True)
process.ctppsProtons.doMultiRPReconstruction = cms.bool(True)

#########################
#       Analysis        #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")

process.ggll_aod.triggersList = process.hltFilter.HLTPaths
process.ggll_aod.leptonsType = cms.string('Muon')
#process.ggll_aod.leptonsType = cms.string('ElectronMuon')
#process.ggll_aod.leptonsType = cms.string('Electron')
process.ggll_aod.runOnMC = cms.bool(runOnMC)
process.ggll_aod.fetchProtons = cms.bool(True)
process.ggll_aod.saveExtraTracks = cms.bool(False)
process.ggll_aod.year = cms.string('2017')


# E/gamma identification
process.ggll_aod.eleIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
    tightLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
)
process.ggll_aod.phoIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp90'),
    tightLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp80'),
)
#process.ggll_aod.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")
#process.ggll_aod.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80")
#process.ggll_aod.phoMediumIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90")
#process.ggll_aod.phoTightIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp80")

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.hltFilter*
    process.egmPhotonIDSequence*
    process.egmGsfElectronIDSequence*
    # Rerun lots of things!
    process.totemRPLocalReconstruction *
    process.ctppsDiamondLocalTracks*
    process.ctppsPixelLocalReconstruction *
    process.ctppsLocalTrackLiteProducer *
    process.ctppsProtons *
    process.ggll_aod
)

#process.outpath = cms.EndPath(process.out, patAlgosToolsTask)
process.outpath = cms.EndPath(patAlgosToolsTask)

print process.dumpPython()
