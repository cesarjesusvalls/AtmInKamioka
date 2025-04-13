from ROOT import TCanvas
from ROOT import TH1D
import ROOT
import numpy as np
import os
from core import newEarth
from wrapper_Prob3 import BargerPropagator

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

class oscProb:
    def __init__(self, cfg:dict):
        self._cfg = cfg
        self.kappa = cfg['dial'].get('kappa', 0.2)
        self.deltaR = cfg['dial']['deltaR']
        self._earth = newEarth(cfg, self.kappa, self.deltaR)
        self._outdir = cfg['output']['osc-outdir']
        os.makedirs(self._outdir, exist_ok=True)

        self._FULL_PREM = self._earth.raw_PREM
        self._dump_prem()

        self.bp_true = None
        self.hist_true = None

        # Define the oscillation parameters
        self.t12 = self._cfg['osc3']['oscpar'].get('t12', 0.308)
        self.t13 = self._cfg['osc3']['oscpar'].get('t13', 0.0215)
        self.t23 = self._cfg['osc3']['oscpar'].get('t23', 0.5)
        self.dm21 = self._cfg['osc3']['oscpar'].get('dm21', 7.54E-5)
        self.mAtm = self._cfg['osc3']['oscpar'].get('mAtm', 2.50E-3)
        self.delta = self._cfg['osc3']['oscpar'].get('delta', 0.0)

        self.cosz = self._cfg['osc3']['zenith'].get('cosz', -1.0)
        self.ks = self._cfg['osc3'].get('ks', 1)
        self.nutype = self._cfg['osc3'].get('nutype', 1)
        self.prod_height = self._cfg['osc3'].get('prod_height', 20.0)

        self.nu_origin = self._cfg['osc3'].get('nu_origin', 2)
        self.nu_end = self._cfg['osc3'].get('nu_end', 1)

        self.nu_flavor = {1: '#nu_{e}', 2: '#nu_{#mu}', 3: '#nu_{#tau}', -1: '#bar{#nu}_{e}', -2: '#bar{#nu}_{#mu}', -3: '#bar{#nu}_{#tau}'}
        self.MO = {1: 'NO', -1: 'IO'}

        self.Enu_min = self._cfg['osc3'].get('Enu_min', 1)
        self.Enu_max = self._cfg['osc3'].get('Enu_max', 21)

        self.L = ROOT.TLegend(0.1, 0.75, 0.3, 0.9)
        self._prepare_hist()

    def _dump_prem(self):
        np.savetxt(f'{self._outdir}/FULL_PREM.dat', self._FULL_PREM, fmt='%.3f', delimiter='\t')
        
    def _load_barger_prop(self):
        # Create an instance of BargerPropagator
        self.bp_true = BargerPropagator.PyBargerPropagator(f'{self._outdir}/FULL_PREM.dat')
        self.bp_true.definePath(self.cosz, self.prod_height, True)

        self._compute_osc_prob()
    def _prepare_hist(self):
        self.hist_true = TH1D('hist_true',
                              self.nu_flavor[self.nutype*self.nu_origin] + f"#rightarrow" + self.nu_flavor[self.nutype*self.nu_end]
                              + f"cosZ={self.cosz} sin23={self.t23}" + self.MO[self.mAtm/abs(self.mAtm)], 10000, self.Enu_min, self.Enu_max)
        self.hist_mod = TH1D('hist_mod',
                              self.nu_flavor[self.nutype * self.nu_origin] + f"#rightarrow" + self.nu_flavor[self.nutype * self.nu_end]
                              + f"cosZ={self.cosz} sin23={self.t23}" + self.MO[self.mAtm / abs(self.mAtm)], 10000,
                              self.Enu_min, self.Enu_max)
        self.hist_diff = TH1D('hist_diff', '', 10000, self.Enu_min, self.Enu_max)

        self.hist_true.SetLineColor(ROOT.kOrange+2)
        self.hist_mod.SetLineColor(ROOT.kAzure-1)
        self.hist_diff.SetLineColor(ROOT.kBlack)

        self.L.AddEntry(self.hist_true, 'True Earth', 'l')
        self.L.AddEntry(self.hist_mod, 'Modified Earth', 'l')

    # Compute some oscillation probabilities
    def _compute_osc_prob(self):
        for i in range(1, self.hist_true.GetNbinsX()):
            energy = self.hist_true.GetBinCenter(i)
            self.bp_true.setMNS(self.t12, self.t13, self.t23, self.dm21, self.mAtm, self.delta, energy, self.ks, self.nutype)
            self.bp_true.propagate(self.nutype)
            self.hist_true.SetBinContent(i, self.bp_true.getProb(self.nu_origin, self.nu_end))

    def get_osc_prob(self, energy, cosZ, nutype, flav_ini, flav_end):
        self.bp_true.definePath(cosZ, self.prod_height, True)
        self.bp_true.setMNS(self.t12, self.t13, self.t23, self.dm21, self.mAtm, self.delta, energy, self.ks, nutype)
        self.bp_true.propagate(nutype)
        return self.bp_true.getProb(flav_ini, flav_end)

    def set_Earth_model(self):
        self._earth = newEarth(self._cfg, self._earth.kappa, self._earth.deltaR)
        np.savetxt(f'{self._outdir}/TEMP_PREM.dat', self._earth.mod_PREM, fmt='%.3f', delimiter='\t')
        self.bp_true = BargerPropagator.PyBargerPropagator(f'{self._outdir}/TEMP_PREM.dat')
