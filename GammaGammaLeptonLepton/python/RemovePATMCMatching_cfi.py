def removePATMCMatching(process, names=['All']):
    from PhysicsTools.PatAlgos.tools.coreTools import runOnData, removeMCMatching
    removeMCMatching(process, names)

    from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import electronMatch
    getattr(process, 'makePatElectrons').remove(electronMatch)

    from PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi import muonMatch
    getattr(process, 'makePatMuons').remove(muonMatch)

    from PhysicsTools.PatAlgos.mcMatchLayer0.tauMatch_cfi import tauMatch
    from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff import tauGenJets, tauGenJetsSelectorAllHadrons, tauGenJetMatch
    taus = getattr(process, 'makePatTaus')
    taus.remove(tauMatch)
    taus.remove(tauGenJets)
    taus.remove(tauGenJetsSelectorAllHadrons)
    taus.remove(tauGenJetMatch)

    from PhysicsTools.PatAlgos.mcMatchLayer0.photonMatch_cfi import photonMatch
    getattr(process, 'makePatPhotons').remove(photonMatch)

    from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetPartonMatch, patJetGenJetMatch
    from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import patJetFlavourIdLegacy, patJetFlavourId
    jets = getattr(process, 'makePatJets')
    jets.remove(patJetPartonMatch)
    jets.remove(patJetGenJetMatch)
    jets.remove(patJetFlavourIdLegacy)
    jets.remove(patJetFlavourId)

