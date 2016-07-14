#!/bin/sh
#JSON_FILE = ~/work/photon-analysis/CMSSW_8_0_8_patch1/src/Cert_274388-275125_13TeV_PromptReco_Collisions16_JSON.txt
JSON_FILE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt
PILEUP_FILE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
pileupCalc.py -i ${JSON_FILE} --inputLumiJSON ${PILEUP_FILE} --calcMode true --minBiasXsec 71300 --maxPileupBin 50 --numPileupBins 50 Cert_274388-275125_13TeV_PromptReco_Collisions16_Pileup.root
