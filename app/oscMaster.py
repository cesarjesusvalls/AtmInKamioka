import numpy as np
import os
from wrapper_Prob3 import BargerPropagator

class oscMaster:
    def __init__(self,  cfg:dict):
        PREM_file = cfg['data']['filepath']
        
        self.bp_true = BargerPropagator.PyBargerPropagator(PREM_file)

        self.t12   = cfg['osc3']['oscpar'].get('t12', 0.308)
        self.t13   = 0.022
        self.t23   = 0.58
        self.dm21  = cfg['osc3']['oscpar'].get('dm21', 7.54E-5)
        self.mAtm  = 2.4e-3
        self.delta = 4.53

        self.ks = cfg['osc3'].get('ks', 1)
        self.nutype = cfg['osc3'].get('nutype', 1)
        self.prod_height = cfg['osc3'].get('prod_height', 20.0)

        self.nu_flavor = {1: '#nu_{e}', 2: '#nu_{#mu}', 3: '#nu_{#tau}', -1: '#bar{#nu}_{e}', -2: '#bar{#nu}_{#mu}', -3: '#bar{#nu}_{#tau}'}
        self.MO = {1: 'NO', -1: 'IO'}

    def get_osc_prob(self, energy, cosZ, nutype, flav_ini, flav_end):
        self.bp_true.definePath(cosZ, self.prod_height, True)
        self.bp_true.setMNS(self.t12, self.t13, self.t23, self.dm21, self.mAtm, self.delta, energy, self.ks, nutype)
        self.bp_true.propagate(nutype)
        return self.bp_true.getProb(flav_ini, flav_end)