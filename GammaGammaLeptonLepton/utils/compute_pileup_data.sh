#!/bin/sh
pileupCalc.py -i ~/work/photon-analysis/CMSSW_8_0_8_patch1/src/Cert_274388-275125_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 71300 --maxPileupBin 50 --numPileupBins 50 Cert_274388-275125_13TeV_PromptReco_Collisions16_Pileup.root
