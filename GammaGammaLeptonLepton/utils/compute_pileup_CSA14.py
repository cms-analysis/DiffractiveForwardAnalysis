#!/usr/bin/python

from ROOT import TH1D
from SimGeneral.MixingModule.mix_CSA14_50ns_PoissonOOTPU_cfi import mix

prob = TH1D('pileup', '', 50, 0., 50.)

i = 0
for value in mix.input.nbPileupEvents.probValue:
    prob.Fill(i, value)
    i += 1

prob.SaveAs('pileup_CSA14_50ns_PoissonOOTPU.root')
