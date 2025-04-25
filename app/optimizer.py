from app.path import *
from app.dataManager import *
import math
import matplotlib.pyplot as plt

class optManager:
    """
    A class to evaluate the current configuration likelihood and control parameter optimization
    """
    def __init__(self, scale_to_HK=True, tag=None):

        self.outpath = base_dir_path+'../profile_LL_files/'
        self.det = 'SK'
        if scale_to_HK:
            self.det= 'HK'

        data_folder  = base_dir_path+'../data/unoscillated'
        binning_file = base_dir_path+'../data/sk_2023_BinInfo.txt'
        self.tag = tag

        if self.tag is not None:
            data_folder = base_dir_path+'../fake_data/FDS_'+tag+'/unoscillated/'
            binning_file= base_dir_path+'../fake_data/FDS_'+tag+'/binning.txt'
            self.tag = 'FDS_'+self.tag
        else:
            self.tag = 'Nominal'

        self.ana_master = AnaMaster(data_folder=data_folder, binning_file=binning_file, scale_to_HK=scale_to_HK)
        self.set_nominal_params()
        self.ana_master.osc_weight_all()
        self.ana_master.fill_histograms()
        self.nominal_counts = self.get_counts()
        self.set_nominal_params()
        self.set_params_grid()

    def set_nominal_params(self):
        self.nominal_t23 = 0.58
        self.nominal_t13 = 0.022
        self.nominal_mAtm = 2.4e-3
        self.nominal_dcp = 4.53

    def set_params_grid(self):
        # with these definitions all grids should evaluate the nominal point (where we expect LL=0).

        step_t23 = 0.02
        self.list_of_t23 = np.arange(0.35, 0.65+step_t23, step_t23)

        offset = -0.0031853100000001078
        self.list_of_dcp = np.array([-math.pi*.999] + list(np.arange(-3.1+0.45+offset,3.1,0.45)) + [math.pi*.999])

        step_mAtm = 0.1
        self.list_of_mAtm = np.arange(1.70, 3.1+step_mAtm, step_mAtm)*1e-3

        step_t13 = 0.0055
        self.list_of_t13 = np.arange(0., 0.075+step_t13, step_t13)

    def reset_to_nominal_params(self):
        self.ana_master.OscProb.t23 = self.nominal_t23
        self.ana_master.OscProb.t13 = self.nominal_t13
        self.ana_master.OscProb.mAtm = self.nominal_mAtm
        self.ana_master.OscProb.delta = self.nominal_dcp

    def get_counts(self):
        tot_counts = []
        for i in range(len(self.ana_master.samples)):
            s = self.ana_master.samples[i]
            if valid_sample_condition(s.name):
                tot_counts += list(get_counts_from_hist(self.ana_master.samples[i].data_hist))
        return tot_counts

    def calc_1D_profile_DeltaChi2(self, case):
        self.reset_to_nominal_params()

        def calc_Chi2():
            self.ana_master.osc_weight_all()
            self.ana_master.fill_histograms()
            return calc_chi2(self.nominal_counts, self.get_counts())

        list_of_chi2 = []
        outname = None

        if case == 't13':
            print('doing t13')
            for t13 in self.list_of_t13:
                print(t13)
                self.ana_master.OscProb.t13 = t13
                list_of_chi2.append(calc_Chi2())

            df = pd.DataFrame({'x': np.array(self.list_of_t13), 'chi2': np.array(list_of_chi2)})
            outname = self.outpath+self.det+'_'+self.tag+'_'+case+'_1D.csv'
            df.to_csv(outname, index=False)

        elif case == 't23':
            for t23 in self.list_of_t23:
                self.ana_master.OscProb.t23 = t23
                list_of_chi2.append(calc_Chi2())

            df = pd.DataFrame({'x': self.list_of_t23, 'chi2': list_of_chi2})
            outname = self.outpath+self.det+'_'+self.tag+'_'+case+'_1D.csv'
            df.to_csv(outname, index=False)

        elif case == 'mAtm':
            for mAtm in self.list_of_mAtm:
                self.ana_master.OscProb.mAtm = mAtm
                list_of_chi2.append(calc_Chi2())

            df = pd.DataFrame({'x': self.list_of_mAtm, 'chi2': list_of_chi2})
            outname = self.outpath+self.det+'_'+self.tag+'_'+case+'_1D.csv'
            df.to_csv(outname, index=False)

        elif case == 'dcp':
            for dcp in self.list_of_dcp:
                self.ana_master.OscProb.delta = dcp
                list_of_chi2.append(calc_Chi2())

            df = pd.DataFrame({'x': self.list_of_dcp, 'chi2': list_of_chi2})
            outname = self.outpath+self.det+'_'+self.tag+'_'+case+'_1D.csv'
            df.to_csv(outname, index=False)

        self.reset_to_nominal_params()
        return outname


    def plot_1D(self, case):

        df = pd.read_csv(base_dir_path+'../bin/'+self.det+'_'+case+'_1D.csv')
        df.plot(df.x, df.chi2, color='black', zorder=-1)
        df.scatter(df.x, df.chi2, marker='o', color='red')

        if case == 't13':
            plt.gca().set_xlabel('$\\sin^2 \\theta_{13}$ ')
        elif case == 't23':
            plt.gca().set_xlabel('$\\sin^2 \\theta_{23}$ ')
        elif case == 'mAtm':
            plt.gca().set_xlabel('$\\Delta m^2_{23}$ ')
        elif case == 'dcp':
            plt.gca().set_xlabel('$\\delta_{CP}$ ')

        plt.gca().set_ylabel('$\\Delta \\chi^2$')
        plt.tight_layout()

        plt.gca().set_ylim(0,16)
        plt.gca().set_xlim(0.,0.07)


    def calc_all_1D_profiles(self):
        for x in ['t13', 't23', 'mAtm', 'dcp']:
            print('Processing ', x)
            opt.calc_1D_profile_DeltaChi2(x)
