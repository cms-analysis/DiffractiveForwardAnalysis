import FWCore.ParameterSet.Config as cms

genPUdist = cms.EDAnalyzer('GeneratePUdist',
                           RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
                           outfilename = cms.untracked.string('puhistos.root'),
##                           PileupFile = cms.string("PUHistos.root"),
##                           DataPileupFile = cms.string("PUHistos.root"),
)
