ROOT ntuple producer for the search of central exclusive two-photon induced events at CMS

***IMPORTANT NOTE***

This branch requires a version of CMSSW from CMSSW_9_4_0 on.

### Electron/photon ID

The `Egamma` POG recommendations for multivariate electron and photon identification 2017+ data analysis is, as listed [here](https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_Recipe_for_regul):

```sh
cmsrel CMSSW_9_4_0
cd CMSSW_9_4_0/src
cmsenv
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3 # electron identification
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP # photon identification
scram b -j 8
# Add the area containing the MVA weights (from cms-data, to appear in `external`).
# Note: the `external` area appears after `scram build` is run at least once, as above
cd $CMSSW_BASE/external
# below, you may have a different architecture, this is just one example from lxplus
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
# Go back to the src/
cd $CMSSW_BASE/src
```

### Analysis package

Once the Egamma MVA ID packages are retrieved, clone this repository in your `$CMSSW_BASE/src` environment:

```sh
git clone git@github.com:cms-analysis/DiffractiveForwardAnalysis.git
```

If the authentication fails in any of this or previous stage, it might be that you did not map your `lxplus`' public key to your GitHub account.
In that case, follow [this procedure](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/).
