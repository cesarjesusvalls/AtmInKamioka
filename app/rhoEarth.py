import numpy as np
import matplotlib.pyplot as plt
import os
from core import newEarth

class showDensityComparison:
    def __init__(self, cfg:dict):
        self._cfg = cfg
        self._outdir = cfg['output']['outdir']
        os.makedirs(self._outdir, exist_ok=True)
        self.kappa = cfg['dial'].get('kappa', 0.2)
        self.deltaR = cfg['dial']['deltaR']
        self._earth = newEarth(cfg, self.kappa, self.deltaR)
        self._FULL_PREM = self._earth.raw_PREM
        # if self._earth.kappa == 0:
        #     Alter_PREM = np.loadtxt('/Users/seanxia/Documents/Research/nu_Tomography/HKTomography/config/PREM.dat')
        #     Alter_PREM_R = Alter_PREM[:, 0]
        #     Alter_PREM_R_fine = np.linspace(0, 6371, 6370)
        #     Alter_PREM_RHO = Alter_PREM[:, 1]

        #     Alter_PREM_RHO_fine = np.zeros(len(Alter_PREM_R_fine))
        #     index_coarse = 0
        #     for i in range(len(Alter_PREM_R)):
        #         Alter_PREM_RHO_fine[i] = Alter_PREM_RHO[index_coarse]
        #         if Alter_PREM5_R[i] > Alter_PREM_R[index_coarse]:
        #             index_coarse += 1

        #     self._ALTER_PREM = np.column_stack((Alter_PREM_R_fine, Alter_PREM_RHO_fine))

        # else:
        self._ALTER_PREM = self._earth.mod_PREM

    def plotDensity(self):
        plt.plot(self._ALTER_PREM[:,0], self._ALTER_PREM[:, 1], color='orange', label='ALTER_PREM')
        plt.plot(self._FULL_PREM[:, 0], self._FULL_PREM[:, 1], color='blue', label='FULL_PREM')
        plt.xlabel('Radius (km)')
        plt.ylabel('Earth Density (g/cm$^3$)')
        plt.legend()
        plt.savefig(f'{self._outdir}/ALTER_PREM_vs_FULL_PREM_kappa_{self._earth.kappa}_dRIC_{self._earth.deltaR[0]}_dROC_{self._earth.deltaR[1]}_dRMT_{self._earth.deltaR[2]}.png')
        plt.show()

        print(f'True earth mass is {self._earth.earth_mass/1000:e} kg, new conserved earth mass with kappa {self._earth.kappa} is {self._earth.new_earth_mass/1000:e} kg')
        print(f'True earth MOI is {self._earth.earth_moi/1E7:e} kg m^2, new conserved earth mass with kappa {self._earth.kappa} is {self._earth.new_earth_moi/1E7:e} kg m^2')