import FWCore.ParameterSet.Config as cms

ggll = cms.EDAnalyzer(
    'GammaGammaLL',
    SqrtS = cms.double(8000.),
    HLTMenuLabel = cms.string("HLT"),
    LeptonsType = cms.InputTag('electron', 'muon'),
    isoValInputTags = cms.VInputTag(
        cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
        cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
        cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
    ),
    beamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    conversionsInputTag = cms.InputTag('allConversions'),
    rhoIsoInputTag = cms.InputTag('kt6PFJetsForIsolation', 'rho'),
    RunOnMC = cms.untracked.bool(True),
    MCAcceptPtCut = cms.untracked.double(0.),
    MCAcceptEtaCut = cms.untracked.double(-1.),
    GenParticlesCollectionLabel = cms.InputTag('genParticles'),
)
