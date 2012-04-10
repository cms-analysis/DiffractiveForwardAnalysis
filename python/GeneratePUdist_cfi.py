import FWCore.ParameterSet.Config as cms

genPUdist = cms.EDAnalyzer('GeneratePUdist',
                           RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
##                           PileupFile = cms.string("PUHistos.root"),
##                           DataPileupFile = cms.string("PUHistos.root"),
)
