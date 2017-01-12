### (Cut-based) electron ID
#### Summer16 recipe
To be introduced in `8_0_X` (`X`&ge;10)

Following the [Egamma POG recipe](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_8_0):
```sh
git cms-merge-topic ikrav:egm_id_80X_v3
scram b -j 10
```
And replace the ```ele*IdMap``` input tags in ```DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi``` by the new electron ID tags:
- ```eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),```
- ```eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")```
- ```eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")```
- ```eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")```
