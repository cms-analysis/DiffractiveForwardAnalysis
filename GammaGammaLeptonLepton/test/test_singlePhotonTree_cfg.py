import FWCore.ParameterSet.Config as cms

process = cms.Process('analysis')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/data/Run2016G/DoubleEG/AOD/23Sep2016-v1/100000/0042DBD3-BA8E-E611-919E-002481ACDAA8.root',
#'/store/data/Run2017C/DoubleMuon/AOD/12Sep2017-v1/10000/029F251F-B1A2-E711-AAC3-001E67792890.root',
#'/store/data/Run2017B/SinglePhoton/MINIAOD/17Nov2017-v1/30000/10F66FEF-B7D5-E711-9003-008CFAC93F3C.root',
'/store/data/Run2017B/SinglePhoton/MINIAOD/12Sep2017-v1/100000/F0CB68CC-0EA3-E711-A877-0025905B85BE.root',
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
#process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltHighLevel.HLTPaths = cms.vstring(
    'HLT_Photon33_v*',
    'HLT_Photon50_v*',
    'HLT_Photon75_v*',
    'HLT_Photon90_v*',
)
#process.hltHighLevel.throw = cms.bool(False)

#########################
#      Preskimming      #
#########################

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#########################
#       Photon ID       #
#########################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff', setupVIDPhotonSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff', setupVIDPhotonSelection)

#########################
#    Tree production    #
#########################

process.load('DiffractiveForwardAnalysis.GammaGammaLeptonLepton.singlePhotonTreeProducer_cfi')
process.singlePhotonTreeProducer.runOnMC = cms.bool(False)
process.singlePhotonTreeProducer.triggerResults = process.hltHighLevel.TriggerResultsTag
process.singlePhotonTreeProducer.triggersList = process.hltHighLevel.HLTPaths
process.singlePhotonTreeProducer.phoWP80IdMap = cms.InputTag('egmPhotonIDs:mvaPhoID-RunIIFall17-v1-wp80')
process.singlePhotonTreeProducer.phoWP90IdMap = cms.InputTag('egmPhotonIDs:mvaPhoID-RunIIFall17-v1-wp90')

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.hltHighLevel*
    process.egmPhotonIDSequence*
    process.singlePhotonTreeProducer
)
