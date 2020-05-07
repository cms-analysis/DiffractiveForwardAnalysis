import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

runOnMC = True
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
process.MessageLogger.cerr.FwkReport.reportEvery = 1

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'file:/tmp/jjhollar/ctppsSimElEl_xangle120_divergence30.root'
#'file:/tmp/jjhollar/ctppsSimElEl_xangle150_divergence0.root'
#/store/data/Run2018A/DoubleMuon/AOD/17Sep2018-v2/60000/32D7FDE9-748B-ED43-A45D-587019CC5D92.root'
 #       'file:/tmp/jjhollar/32D7FDE9-748B-ED43-A45D-587019CC5D92.root'
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
# declare global tag
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one                                     
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, "106X_dataRun2_testPPS_v1")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc");

#JH - for tests of 2018 alignment candidate
#from CalibPPS.ESProducers.ctppsRPAlignmentCorrectionsDataESSourceXML_cfi import *
process.ctppsRPAlignmentCorrectionsDataESSourceXML = cms.ESSource("CTPPSRPAlignmentCorrectionsDataESSourceXML",
    verbosity = cms.untracked.uint32(0),

    MeasuredFiles = cms.vstring(),
    RealFiles = cms.vstring(),
    MisalignedFiles = cms.vstring()
)

process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = cms.vstring(
  # internal pixel alignment
  "CalibPPS/ESProducers/data/alignment/RPixGeometryCorrections-2017-2018.xml",

  # global alignment of tracking RPs (strips and pixels)
  "CalibPPS/ESProducers/data/alignment/alignment_export_2019_05_09.1.xml",

  # alignment of timing RPs (diamonds and UFSDs)
  "CalibPPS/ESProducers/data/alignment/timing_RP_2017_preTS2_v2.xml",
  "CalibPPS/ESProducers/data/alignment/timing_RP_2017_postTS2_v2.xml",
  "CalibPPS/ESProducers/data/alignment/timing_RP_2018_v2.xml"
)

process.esPreferLocalAlignment = cms.ESPrefer("CTPPSRPAlignmentCorrectionsDataESSourceXML", "ctppsRPAlignmentCorrectionsDataESSourceXML")





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

#JH - ESPrefer to get optical functions from CTPPSOpticalFunctionsESSource instead of global tag for now

# For testing on old prompt reco data - need to rerun pixel tracking + LiteTracks first
# Should not be needed for final version running on re-RECO data
process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","CTPPSFastSimulation")
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(False)

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
process.ggll_aod.year = cms.string('2017MC')


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
    fileName = cms.string('output_xangle120.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.hltFilter*
    process.egmPhotonIDSequence*
    process.egmGsfElectronIDSequence*
    # For testing on Prompt Reco. Rerun pixel+diamond tracking & LiteTracks
    #    process.totemRPLocalReconstruction*
    #    process.ctppsPixelLocalReconstruction*
    #    process.ctppsDiamondLocalTracks*
    #    process.ctppsLocalTrackLiteProducer*
    # Only run high-level proton reco
    #    process.ctppsProtons *                                                                                            
    # Rerun lots of things!
#    process.totemRPUVPatternFinder * 
#    process.totemRPLocalTrackFitter * 
#    process.ctppsDiamondRecHits * 
#    process.ctppsDiamondLocalTracks * 
#    process.ctppsPixelLocalTracks * 
#    process.ctppsLocalTrackLiteProducer * 
#    process.ctppsProtons *
    process.ggll_aod
)

#process.outpath = cms.EndPath(process.out, patAlgosToolsTask)
process.outpath = cms.EndPath(patAlgosToolsTask)

print process.dumpPython()
