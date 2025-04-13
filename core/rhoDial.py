import numpy as np
def _calc_layer_inertia(r:np.ndarray, rho:np.ndarray):
    # r is in km, rho is in g/cm^3
    return np.trapz(8/3*np.pi*rho*(r*1E5)**4, r*1E5).item()
def _calc_layer_mass(r:np.ndarray, rho:np.ndarray):
    # r is in km, rho is in g/cm^3
    return np.trapz(4.*np.pi*rho*(r*1E5)**2, r*1E5).item()

def _poly_fit_layer(r:np.ndarray, rho:np.ndarray, power:int=3):
    # r is in km, rho is in g/cm^3
    coeff = np.polyfit(r[1:-1], rho[1:-1], power) # exclude the edges
    poly_func = np.poly1d(coeff)
    return poly_func

class newEarth:
    def __init__(self, cfg:dict, kappa:float, deltaR:list):
        '''
        Constructor for the newEarth class
        :param cfg:
            config file that contains the following:
                data: path to the original PREM model density file
                var_rhos: number of variation steps in density to consider the equilibriums
                layers: name of layers to be considered, e.g. center, IC, OC, MT, CR
                save_new: flag to save the new PREM model, default to false
                IC_upper_limit: upper limit of the IC layer density factor
                verbose: print out the process
                check_thre: threshold for agreement check of mass and MOI
        :param kappa:
            dial factor for the OC layer density (in fraction)
        :param deltaR:
            list of dial factors for the layer boundaries (in km)
        '''

        self._file = cfg['data']['filepath']
        self.var_rhos = int(cfg['dial'].get('var_rhos', 10))
        self.save_new = cfg['data'].get('save_new', False)
        self.verbose = cfg['dial'].get('verbose', False)
        self.threshold = cfg['dial'].get('check_thre', 0.005)
        self.kappa = kappa
        self.deltaR = deltaR

        self.IC_upper_kappa = cfg['dial'].get('IC_upper_limit', 1.)
        if not isinstance(self.deltaR, list):
            return RuntimeError('deltaR must be a list')

        self.earth_keys = cfg['earth']['layers']
        self.nlayers = len(self.earth_keys) - 1 # exclude the center
        assert len(self.deltaR) == self.nlayers, 'deltaR must have the same length as the number of layers.'
        self.earth_layers_R = {}
        self.earth_layers_Rshift = {}
        self.earth_layers_index = {}
        self.earth_layers_index_shift = {}
        self.earth_layers_poly = {}

        self.earth_raw_layer_mass = {}
        self.earth_raw_layer_MOI = {}

        self.r = None
        self.rho = None

        for i, key in enumerate(self.earth_keys):
            if key == 'center':
                self.earth_layers_R[key] = 0
                self.earth_layers_Rshift[key] = 0
            else:
                self.earth_layers_R[key] = cfg['earth'][key+'_R']
                self.earth_layers_Rshift[key] = cfg['earth'][key+'_R'] + self.deltaR[i-1]
                self.earth_layers_index[key] = None
                self.earth_layers_index_shift[key] = None
                self.earth_layers_poly[key] = None

        self.raw_PREM = None
        self.mod_PREM = None
        self.earth_mass = None
        self.new_earth_mass = None
        self.earth_MOI = None

        self._load_data()
        self._fit_layers()
        self._set_rho_OC()
        self._set_compensation_rho()

    def _load_data(self):

        self.r = np.loadtxt(self._file)[:, 0]
        self.rho = np.loadtxt(self._file)[:, 1]
        assert self.r is not None, 'Radius data not loaded'
        assert self.rho is not None, 'Density data not loaded'

        print(f'File {self._file} is loaded.')

        self.raw_PREM = np.column_stack((self.r, self.rho))
        self.mod_PREM = self.raw_PREM.copy()

        for i, key in enumerate(self.earth_keys):
            if i == 0:
                continue
            self.earth_layers_index[key] = np.where((self.raw_PREM[:, 0] >= self.earth_layers_R[self.earth_keys[i-1]]) & (self.raw_PREM[:, 0] < self.earth_layers_R[key]))
            self.earth_layers_index_shift[key] = self.earth_layers_index[key]

    def _fit_layers(self):
        #print(self.raw_PREM[self.IC_index, 0])

        for i, key in enumerate(self.earth_keys):
            if i == 0:
                continue
            self.earth_layers_poly[key] = _poly_fit_layer(self.raw_PREM[self.earth_layers_index[key], 0].squeeze(), self.raw_PREM[self.earth_layers_index[key], 1].squeeze())

        if max(self.deltaR) > 0:
            # so consider the shift of layer boundaries
            for i, key in enumerate(self.earth_keys):
                if i == 0:
                    continue
                self.earth_layers_index_shift[key] = np.where((self.raw_PREM[:, 0] >= self.earth_layers_Rshift[self.earth_keys[i-1]]) & (self.raw_PREM[:, 0] < self.earth_layers_Rshift[key]))

            self._extrapolate_layers()

    def _extrapolate_layers(self):
        for i, key in enumerate(self.earth_keys):
            if i == 0:
                continue
            if i == self.nlayers:
                break
            self.mod_PREM[self.earth_layers_index_shift[key], 1] = np.expand_dims(self.earth_layers_poly[key](self.mod_PREM[self.earth_layers_index_shift[key], 0].squeeze()), axis = 0)

    def _set_rho_OC(self):
        if self.raw_PREM is None:
            self._load_data()
        #self.mod_PREM[self.OC_index_shift, 1] *= (1. + self.kappa)
        assert self.earth_layers_index_shift['OC'] is not None, 'OC layer index not set'
        self.mod_PREM[self.earth_layers_index_shift['OC'], 1] *= (1 + self.kappa)

    def _set_compensation_rho(self):

        for i, key in enumerate(self.earth_keys):
            if i == 0:
                continue
            self.earth_raw_layer_mass[key] = _calc_layer_mass(self.raw_PREM[self.earth_layers_index[key], 0], self.raw_PREM[self.earth_layers_index[key], 1])
            self.earth_raw_layer_MOI[key] = _calc_layer_inertia(self.raw_PREM[self.earth_layers_index[key], 0], self.raw_PREM[self.earth_layers_index[key], 1])

        earth_mod_layer_mass = {}
        earth_mod_layer_MOI = {}
        earth_mod_rho = {}

        for i, key in enumerate(self.earth_keys):
            if i == 0:
                continue
            earth_mod_rho[key] = self.mod_PREM[self.earth_layers_index_shift[key], 1]

        assert 'OC' in self.earth_keys, 'OC is not defined as a layer in the earth model.'
        earth_mod_layer_mass['OC'] = _calc_layer_mass(self.mod_PREM[self.earth_layers_index_shift['OC'], 0], self.mod_PREM[self.earth_layers_index_shift['OC'], 1])
        earth_mod_layer_MOI['OC'] = _calc_layer_inertia(self.mod_PREM[self.earth_layers_index_shift['OC'], 0], self.mod_PREM[self.earth_layers_index_shift['OC'], 1])

        # #ensure hydrostatic equivalence
        for i, key in enumerate(self.earth_keys):
            if i == 0 or key == 'MT':
                continue
            if i == self.nlayers:
                break
            if np.min(earth_mod_rho[self.earth_keys[i]]) < np.max(earth_mod_rho[self.earth_keys[i+1]]):
                if self.earth_keys[i] != 'OC':
                    self.mod_PREM[self.earth_layers_index_shift[self.earth_keys[i]], 1] += (np.max(earth_mod_rho[self.earth_keys[i+1]])-np.min(earth_mod_rho[self.earth_keys[i]]))
                else:
                    self.mod_PREM[self.earth_layers_index_shift[self.earth_keys[i+1]], 1] += (np.min(earth_mod_rho[self.earth_keys[i]])-np.max(earth_mod_rho[self.earth_keys[i+1]]))
                    self.mod_PREM[self.earth_layers_index_shift[self.earth_keys[i+1]], 1] = np.clip(earth_mod_rho[self.earth_keys[i+1]], 0, np.min(earth_mod_rho[self.earth_keys[i]]))

        self.earth_mass = sum(self.earth_raw_layer_mass.values())
        self.earth_moi = sum(self.earth_raw_layer_MOI.values())

        self.new_earth_mass = self.earth_mass + earth_mod_layer_mass['OC'] - self.earth_raw_layer_mass['OC']
        self.new_earth_moi = self.earth_moi + earth_mod_layer_MOI['OC'] - self.earth_raw_layer_MOI['OC']

        print(f'Default earth mass is {self.earth_mass/1000:e} kg, new non-conserved earth mass with kappa {self.kappa} is {self.new_earth_mass/1000:e} kg')
        print(f'Default earth MOI is {self.earth_moi/1E7:e} kg m^2, new non-conserved earth MOI with kappa {self.kappa} is {self.new_earth_moi/1E7:e} kg m^2')

        IC_lower_kappa = np.max(earth_mod_rho['OC'])/np.min(earth_mod_rho['IC']) - 1.
        MT_upper_kappa = np.min(earth_mod_rho['OC'])/np.max(earth_mod_rho['MT']) - 1.

        break_flag = False
        mass_diff_min = 1.E+12
        moi_diff_min = 1.E+12

        for comp_k in np.linspace(IC_lower_kappa, self.IC_upper_kappa, self.var_rhos):
            if break_flag:
                break
            earth_mod_rho['IC'] = self.mod_PREM[self.earth_layers_index_shift['IC'], 1] * (1 + comp_k)
            if np.min(earth_mod_rho['IC']) < np.max(earth_mod_rho['OC']):
                continue


            for comp_j in np.linspace(-0.999, MT_upper_kappa, self.var_rhos):
                earth_mod_rho['MT'] = self.mod_PREM[self.earth_layers_index_shift['MT'], 1] * (1 + comp_j)
                if np.min(earth_mod_rho['OC']) < np.max(earth_mod_rho['MT']) and self.verbose:
                    print(f"Mantle has larger density than OC at dial factor {comp_j}.")
                    break

                mass_diff = 0
                moi_diff = 0

                for key in self.earth_keys:
                    if key == 'center':
                        continue

                    earth_mod_layer_mass[key] = _calc_layer_mass(self.mod_PREM[self.earth_layers_index_shift[key], 0], earth_mod_rho[key])
                    earth_mod_layer_MOI[key] = _calc_layer_inertia(self.mod_PREM[self.earth_layers_index_shift[key], 0], earth_mod_rho[key])
                    mass_diff += (earth_mod_layer_mass[key] - self.earth_raw_layer_mass[key])
                    moi_diff += (earth_mod_layer_MOI[key] - self.earth_raw_layer_MOI[key])

                if any(earth_mod_layer_mass.values()) < 0 or any(earth_mod_layer_MOI.values()) < 0:
                    print('<0')
                    continue

                if np.abs(mass_diff) < np.abs(mass_diff_min):
                    mass_diff_min = mass_diff

                if np.abs(moi_diff) < np.abs(moi_diff_min):
                    moi_diff_min = moi_diff

                if (np.abs(mass_diff) / self.earth_mass <= self.threshold) and (np.abs(moi_diff) / self.earth_moi <= self.threshold):

                    print(f"Reached new equilibrium with IC rho factor {comp_k} and MT rho factor {comp_j}.")
                    self.mod_PREM[self.earth_layers_index_shift['IC'], 1] = earth_mod_rho['IC']
                    self.mod_PREM[self.earth_layers_index_shift['MT'], 1] = earth_mod_rho['MT']
                    break_flag = True
                    break


        print(f"Minimal fractional difference in mass is {mass_diff_min/np.abs(self.new_earth_mass - self.earth_mass):.2f}, minimal fractional difference in MOI is {np.abs(moi_diff_min)/np.abs(self.new_earth_moi - self.earth_moi):.2f}.")
        assert np.all(self.mod_PREM[:, 1] > 0) and break_flag, f'Cannot conserve earth mass with kappa = {self.kappa}.'

        self.new_earth_mass = sum([_calc_layer_mass(self.mod_PREM[self.earth_layers_index_shift[key], 0], self.mod_PREM[self.earth_layers_index_shift[key], 1]) for key in self.earth_keys if key != 'center'])