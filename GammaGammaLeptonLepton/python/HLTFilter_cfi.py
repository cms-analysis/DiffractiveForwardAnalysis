import FWCore.ParameterSet.Config as cms
import copy

# Trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
hltFilter = copy.deepcopy(hltHighLevel)
hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
hltFilter.HLTPaths = ['HLT_L1DoubleMuOpen','HLT_L1DoubleMuOpen_Tight','HLT_Mu3','HLT_DoubleMu0']
hltFilter.throw = cms.bool(False)
