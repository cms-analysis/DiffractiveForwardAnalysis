import FWCore.ParameterSet.Config as cms

ggll = cms.EDAnalyzer(
    'GammaGammaLL',
    LeptonsType = cms.InputTag('electron', 'muon'),
    isoValInputTags = cms.VInputTag(
        cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
        cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
        cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
    ),
    beamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    conversionsInputTag = cms.InputTag('allConversions'),
    rhoIsoInputTag = cms.InputTag('kt6PFJetsForIsolation', 'rho'),
)
