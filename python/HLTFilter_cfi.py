import FWCore.ParameterSet.Config as cms
import copy

# Trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
hltFilter = copy.deepcopy(hltHighLevel)
hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT8E29")
hltFilter.HLTPaths = ['HLT_DoubleMu3','HLT_Mu3','HLT_DoubleMu0']
