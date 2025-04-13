import os
import pandas as pd
import numpy as np
import ROOT
import copy
from app.path import *
from app import oscProb, oscMaster
import yaml
from app import fluxManager
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt


def valid_sample_condition(name):
    return 'Showering' not in name and 'Sub-GeV' not in name and 'pi^{0}' not in name and '*' not in name
    #return True

def plot_sampling(q_values):

    # Define the quantiles and the corresponding values
    quantiles = np.array([2.3, 15.9, 50, 84.1, 97.7])
    values = q_values
    
    # Step 1: Generate uniformly distributed data
    uniform_data = np.random.uniform(0, 1, 1000000)
    
    # Step 2: Map uniform data to empirical distribution
    # Convert uniform data to percentiles (0 to 100)
    uniform_percentiles = uniform_data * 100
    
    # Use linear interpolation to find corresponding empirical values
    pchip_interpolator = PchipInterpolator(quantiles, values)
    empirical_data = pchip_interpolator(uniform_percentiles)
    
    # Plotting the results
    plt.figure(figsize=(18, 6))
    
    plt.subplot(1, 3, 1)
    plt.hist(uniform_data, bins=50, density=True, color='cornflowerblue')
    plt.title('Uniform Distribution')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    
    # Generate points for plotting the spline
    x_new = np.linspace(min(quantiles), max(quantiles), 100)
    y_new = pchip_interpolator(x_new)
    
    plt.subplot(1, 3, 2)
    plt.plot(x_new, y_new, 'r-', label='Spline')
    plt.plot(quantiles, values, 'ko', label='Quantile Points')
    plt.title('Quantile Points of Empirical Distribution')
    plt.xlabel('Quantiles')
    plt.ylabel('Values')
    plt.legend()
    
    plt.subplot(1, 3, 3)
    plt.hist(empirical_data, bins=50, density=True, color='cornflowerblue')
    plt.title('Transformed Empirical Distribution')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    
    plt.show()

def load_data(file_path):
    data = pd.read_csv(file_path, delim_whitespace=True)
    return data

def load_data_and_create_bins(file_path):
    data_frame = load_data(file_path)

    ana_bins = []
    for index, row in data_frame.iterrows():
        ana_bin = AnaBin(row['sample_name'], row['log_p_lo'], row['log_p_hi'], row['cos_z_lo'], row['cos_z_hi'])
        ana_bins.append(ana_bin)

    return ana_bins

def calc_chi2(expected,observed):
    if len(observed) != len(expected):
        print("Length of observed does not match length of expected")
        return

    N = len(observed)
    chi2 = 0
    for i in range(N):
        E = expected[i]
        O = observed[i]
        chi2 += (E-O)+O*np.log(O/E)

    return chi2


def get_counts_from_hist(hist):
    counts = []
    for i in range(hist.GetXaxis().GetNbins()):
        counts.append(hist.GetBinContent(int(i+1)))
    return counts


class BinManager:
    """
    A class to manage bin meta data
    """
    def __init__(self, bin_info_path=base_dir_path+"../data/sk_2023_BinInfo.txt"):
        self.bin_info_df = load_data(bin_info_path)
        self.bin_names = self.bin_info_df['SampleName']

    def get_bin_names(self, bin_indices):
        names = []
        bin_subset =  self.bin_info_df.loc[bin_indices]
        for i in range(len(bin_subset)):
            names.append(str(bin_subset.iloc[i][['LogP_lo', 'LogP_hi', 'CosZ_lo', 'CosZ_hi']].to_numpy()))

        return names

    def get_indices_with_value_one(self, bin_names):
        indices = [i for i, s in enumerate(bin_names) if (s.strip().endswith('1.0]') \
                   or s.strip().endswith('-0.1 0.0]')) and not s.strip().endswith('-1.0 1.0]')]
        return indices

class AnaBin:
    """
    An analysis bin containing all its relevant information to allow reweight
    """
    def __init__(self, id, reaction, bin_data):
        self.id = id
        self.reaction = reaction
        self.counts = bin_data['Counts']
        self.energy_avg = bin_data['EnergyAvg']
        self.energy_rms = bin_data['EnergyRMS']
        self.cos_z_avg = bin_data['CosZAvg']
        self.cos_z_rms = bin_data['CosZRMS']
        self.weight = 1

        self.input_file_E_quantiles = []
        self.input_file_Z_quantiles = []

        for Q in ['2.3', '15.9', '50.0', '84.1', '97.7']:

            E = bin_data['EnergyQuantile'+Q+'Percent']
            Z = bin_data['CosZQuantile'+Q+'Percent']

            if Z<-1:
                Z=-1
            elif Z>1:
                Z=1

            self.input_file_E_quantiles.append(E)
            self.input_file_Z_quantiles.append(Z)

        # Define the quantiles and the corresponding values
        quantiles = np.array([2.3, 15.9, 50, 84.1, 97.7])
        #print(self.input_file_E_quantiles)
        self.pchip_E_interpolator = PchipInterpolator(quantiles, self.input_file_E_quantiles)
        self.pchip_Z_interpolator = PchipInterpolator(quantiles, self.input_file_Z_quantiles)

        # now we need to break the distribution in small pieces that contain the fraction of events
        # given by the difference in adjacent quantile values (that I am calling quantile edges)
        # and we will evaluate them in the mid point of the small piece.

        Npieces = 15
        self.quantile_edges = np.linspace(0,100, Npieces)
        # the width between consecutive quantile edges
        self.quantile_W = (self.quantile_edges[1:]-self.quantile_edges[:-1])

        # the mid points in those edges
        self.quantile_X = self.quantile_edges[:-1]+self.quantile_W/2

        self.generalized_E_quantiles = self.pchip_E_interpolator(self.quantile_X)
        self.generalized_Z_quantiles = self.pchip_Z_interpolator(self.quantile_X)

        self.logp_lo = None
        self.logp_hi = None
        self.cosz_lo = None
        self.cosz_hi = None

        self.E_bin = None
        self.cosZ_bin = None

        self.ignore_bin = False
        if self.counts == 0 or self.energy_avg < 0.1:
            self.ignore_bin = True


    def set_bin_meta_info(self, bin_data):
        self.logp_lo = bin_data['LogP_lo']
        self.logp_hi = bin_data['LogP_hi']
        self.cosz_lo = bin_data['CosZ_lo']
        self.cosz_hi = bin_data['CosZ_hi']

    def get_flux_bin(self, E, cosZ, E_bins, cosZ_bins):
        found_E_bin = -1
        found_cosZ_bin = -1
        for i, E_bin in enumerate(E_bins):
            if E >= E_bin[0] and E <= E_bin[1]:
                found_E_bin = i
                break
        for i, c_bin in enumerate(cosZ_bins):
            if cosZ >= c_bin[0] and cosZ <= c_bin[1]:
                found_cosZ_bin = i
                break
        return found_E_bin, found_cosZ_bin

    def flux_from_flux_index(self, flux_data, energy, E_bin, cosZ_bin):
        lo_Enu_flux_bin = flux_data[cosZ_bin][E_bin]
        hi_Enu_flux_bin = flux_data[cosZ_bin][E_bin + 1]
        flavors = ['NuMu', 'NuMubar', 'NuE', 'NuEbar']
        flux = {}
        for flav in flavors:
            x = [lo_Enu_flux_bin['Enu'], hi_Enu_flux_bin['Enu']]
            y = [lo_Enu_flux_bin[flav], hi_Enu_flux_bin[flav]]
            flux[flav] = np.interp(energy, x, y)
        return flux

    def set_flux_info(self, flux_manager, E_bins, cosZ_bins):
        if self.ignore_bin:
            return

        # # --- first the flux assuming average energy and cosZ in the bin ---
        self.E_bin, self.cosZ_bin = self.get_flux_bin(self.energy_avg, self.cos_z_avg, E_bins, cosZ_bins)
        self.flux = self.flux_from_flux_index(flux_manager.flux_data, self.energy_avg, self.E_bin, self.cosZ_bin)

    def get_quantile_based_flux(self, OscProb, nu, alpha, beta):

        nu_flav = {'e': 1, 'mu': 2, 'tau': 3}
        nu_type = {'nu': 1, 'nubar': -1}
        weighted_A_to_A_flux = 0
        weighted_B_to_A_flux = 0
        Q_weights = np.array(self.quantile_W)/100

        nominal_A_flux_str = ''
        nominal_B_flux_str = ''
        if alpha == 'e':
            nominal_A_flux_str += 'NuE'
            nominal_B_flux_str += 'NuMu'
        elif alpha == 'mu':
            nominal_A_flux_str += 'NuMu'
            nominal_B_flux_str += 'NuE'

        if nu == 'nubar':
            nominal_A_flux_str += 'bar'
            nominal_B_flux_str += 'bar'

        tot_A_flux = 0
        cnt = 0
        for i, E in enumerate(self.generalized_E_quantiles):
            for j, cosZ in enumerate(self.generalized_Z_quantiles):
                A_flux = self.flux[nominal_A_flux_str]
                B_flux = self.flux[nominal_B_flux_str]
                tot_A_flux += A_flux

                if beta != 'tau':
                    weighted_A_to_A_flux  += A_flux*Q_weights[i]*Q_weights[j]*OscProb.get_osc_prob(E, cosZ, nu_type[nu], nu_flav[alpha], nu_flav[alpha])
                    weighted_B_to_A_flux  += B_flux*Q_weights[i]*Q_weights[j]*OscProb.get_osc_prob(E, cosZ, nu_type[nu], nu_flav[beta], nu_flav[alpha])
                else:
                    weighted_B_to_A_flux  += A_flux*Q_weights[i]*Q_weights[j]*OscProb.get_osc_prob(E, cosZ, nu_type[nu], nu_flav[alpha], nu_flav[beta])

                cnt += 1
        return weighted_A_to_A_flux, weighted_B_to_A_flux, tot_A_flux/cnt


    def quantile_based_osc_reweight(self, OscProb):

        nu_flav = {'e': 1, 'mu': 2, 'tau': 3}
        nu_type = {'nu': 1, 'nubar': -1}

        self.weight = 1
        self.osc_weight = 1

        if self.ignore_bin == False:
            if self.reaction == 'MCNue':
                dis_flux, app_flux, nom_flux = self.get_quantile_based_flux(OscProb,'nu', 'e', 'mu')
                self.osc_weight = (dis_flux+app_flux)/nom_flux

            elif self.reaction == 'MCNueBar':
                dis_flux, app_flux, nom_flux = self.get_quantile_based_flux(OscProb,'nubar', 'e', 'mu')
                self.osc_weight = (dis_flux+app_flux)/nom_flux

            elif self.reaction == 'MCNumu':
                dis_flux, app_flux, nom_flux = self.get_quantile_based_flux(OscProb,'nu', 'mu', 'e')
                self.osc_weight = (dis_flux+app_flux)/nom_flux

            elif self.reaction == 'MCNumuBar':
                dis_flux, app_flux, nom_flux = self.get_quantile_based_flux(OscProb,'nubar', 'mu', 'e')
                self.osc_weight = (dis_flux+app_flux)/nom_flux

            elif self.reaction == 'MCNutau':
                _, numu_app_flux, numu_nom_flux       = self.get_quantile_based_flux(OscProb,'nu', 'mu', 'tau')
                _, nue_app_flux, _                   = self.get_quantile_based_flux(OscProb,'nu', 'e', 'tau')
                _, numubar_app_flux, numubar_nom_flux = self.get_quantile_based_flux(OscProb,'nubar', 'mu', 'tau')
                _, nuebar_app_flux, _                = self.get_quantile_based_flux(OscProb,'nubar', 'e', 'tau')
                self.osc_weight = (numu_app_flux+nue_app_flux+numubar_app_flux+nuebar_app_flux)/(numu_nom_flux+numubar_nom_flux)

        self.weight *= self.osc_weight

    def get_metrics(self):
        cosz_bin_width = self.cosz_hi-self.cosz_lo
        cosz_bin_center = self.cosz_lo + cosz_bin_width/2

        logp_bin_width  = 10**self.logp_hi - 10**self.logp_lo
        logp_bin_center = 10**self.logp_lo + logp_bin_width/2

        cosz_bias = (self.cos_z_avg-cosz_bin_center)/cosz_bin_width
        cosz_width_ratio = self.cos_z_rms/cosz_bin_width

        e_resolution = self.energy_rms/self.energy_avg

        return cosz_bias, cosz_width_ratio, e_resolution



    def resample(self):

        nToys = 10000
        uniform_Z_data = np.random.uniform(0, 100, nToys)
        uniform_E_data = np.random.uniform(0, 100, nToys)

        self.sampled_Zs = self.pchip_Z_interpolator(uniform_Z_data)
        self.sampled_Es = self.pchip_E_interpolator(uniform_E_data)

class AnaReaction:
    """
    A collection of bins sharing the same physics selection
    that can be reweighted
    """
    def __init__(self, reaction_file_path, bin_manager, OscProb):
        self.bin_manager = bin_manager
        self.ana_bins = []
        self.name = os.path.basename(reaction_file_path).replace('sk_2023_', '').replace('.txt', '')
        self.name = self.name.strip('NoOsc').strip('NO').strip('IO')
        print('Reaction Name: ' + self.name)

        reac_info_df = load_data(reaction_file_path)
        for index, row in reac_info_df.iterrows():
            bin = AnaBin(index, self.name, row)
            bin.set_bin_meta_info(bin_manager.bin_info_df.iloc[index])
            self.ana_bins.append(bin)
        self.ana_bins = np.array(self.ana_bins)
        self.nominal_ana_bins = copy.deepcopy(self.ana_bins)
        self.OscProb = OscProb

    def apply_osc_weights(self):
        for bin in self.ana_bins:
            #bin.central_value_based_osc_reweight(self.OscProb)
            bin.quantile_based_osc_reweight(self.OscProb)

    def get_bin_counts(self, bin_indices):
        bin_subset =  self.ana_bins[bin_indices]
        counts = np.array([bin.counts for bin in bin_subset])
        weights = np.array([bin.weight for bin in bin_subset])
        return counts, weights

    def study_bins_reliability(self):
        cosz_biases = []
        cosz_width_ratios = []
        e_resolutions = []
        for bin in self.ana_bins:
            cosz_bias, cosz_width_ratio, e_resolution = bin.get_metrics()
            cosz_biases.append(cosz_bias)
            cosz_width_ratios.append(cosz_width_ratio)
            e_resolutions.append(e_resolution)

        return cosz_biases, cosz_width_ratios, e_resolutions

class AnaSample:
    """
    A collection of bins sharing the same physics selection
    that can be visualized.
    """
    def __init__(self, name, bin_manager, nbins, sample_bins, reaction_names):
        self.name = name
        self.bin_manager = bin_manager
        self.nbins = nbins
        self.bin_indices = sample_bins
        self.bin_lines = []
        self.data_hist = None
        self.reaction_names = ['MCNC','MCNutau','MCNumuBar','MCNumu','MCNueBar','MCNue','MCTotal']

        # some naming exceptions
        self.special_samples = ["Sub-GeV  e-like 1 d.e.*",
                                "Sub-GeV #mu-like 2 d.e.*",
                                "1-ring #pi^{0}-like",
                                "2-ring #pi^{0}-like"]

        self.reaction_names_to_latex = {}
        for reac_name in self.reaction_names:
            if 'Nutau' in reac_name:
                self.reaction_names_to_latex[reac_name] = "#nu_{#tau} & #bar{#nu}_{#tau} CC"
            elif 'NC' in reac_name:
                self.reaction_names_to_latex[reac_name] = "NC"
            elif 'NueBar' in reac_name:
                self.reaction_names_to_latex[reac_name] = "#bar{#nu}_{e} CC"
            elif 'NumuBar' in reac_name:
                self.reaction_names_to_latex[reac_name] = "#bar{#nu}_{#mu} CC"
            elif 'Nue' in reac_name:
                self.reaction_names_to_latex[reac_name] = "#nu_{e} CC"
            elif 'Numu' in reac_name:
                self.reaction_names_to_latex[reac_name] = "#nu_{#mu} CC"
            else:
                self.reaction_names_to_latex[reac_name] = "Data"

        self.bin_names = self.bin_manager.get_bin_names(self.bin_indices)
        self.histograms = {}
        self.custom_palette = [ROOT.kGray, ROOT.kMagenta+3, ROOT.kCyan-8, \
                               ROOT.kAzure+5, ROOT.kRed-8, ROOT.kRed+1, ROOT.kBlack]
        self.create_histograms()
        self.create_stack()

    def get_stat_chi2(self):

        self.set_total_MC_as_data_hist()
        mc = []
        for i in range(len(self.bin_names)):
            mc.append(self.data_hist.GetBinContent(int(i+1)))

        self.set_SK_true_data_as_data_hist()
        dt = []
        for i in range(len(self.bin_names)):
            dt.append(self.data_hist.GetBinContent(int(i+1)))

        chi2 = 0
        for i in range(len(self.bin_names)):
            E = mc[i]
            O = dt[i]
            chi2 += (dt[i]-mc[i])**2/mc[i]

        return chi2

    def create_data_histogram(self, title="Data"):
        data_hist = ROOT.TH1F(self.name+title, title, len(self.bin_names), 0, len(self.bin_names))
        data_hist.GetYaxis().SetTitle("Counts")
        data_hist.SetMarkerStyle(20)
        data_hist.SetMarkerColor(ROOT.kBlack)
        data_hist.SetLineColor(ROOT.kBlack)
        data_hist.SetLineWidth(2)
        return data_hist

    def format_histogram_axis(self, hist):
        if self.name not in self.special_samples:
            positions = self.bin_manager.get_indices_with_value_one(self.bin_names)
            for bin_idx, bin_name in enumerate(self.bin_names):
                if bin_idx == 0 or bin_idx-1 in positions:
                    hist.GetXaxis().SetBinLabel(bin_idx+1, "-1")
                elif bin_idx in positions:
                    if "Up" in self.name:
                        hist.GetXaxis().SetBinLabel(bin_idx+1, "0")
                    else:
                        hist.GetXaxis().SetBinLabel(bin_idx+1, "1")
                else:
                    hist.GetXaxis().SetBinLabel(bin_idx+1, "")
        else:
            for bin_idx, bin_name in enumerate(self.bin_names):
                momentum = [float(x) for x in self.bin_names[bin_idx].strip('[]').split()]
                hist.GetXaxis().SetBinLabel(bin_idx+1, '{:.1f}'.format(momentum[0]))
            hist.GetXaxis().SetTitle("#log_{10} [p (MeV/c)]")

        hist.GetXaxis().SetLabelSize(0.03)

    def create_histograms(self):
        # one histogram per reaction
        for reac_idx, reac_name in enumerate(self.reaction_names):
            hist = ROOT.TH1F(self.name+"  "+reac_name, self.name+"  "+reac_name, len(self.bin_names), 0, len(self.bin_names))
            hist.GetYaxis().SetTitle("Counts")
            if 'Total' not in reac_name:
                hist.SetFillColor(self.custom_palette[reac_idx])
                hist.SetLineWidth(1)
                hist.SetLineColor(ROOT.kBlack)
                self.format_histogram_axis(hist)
                self.histograms[reac_name] = hist

        # one histogram for the data
        self.data_hist = self.create_data_histogram()
        self.format_histogram_axis(hist)

    def create_zenith_angle_edge(self, hist):
        cos_zenith_is_one_bins = self.bin_manager.get_indices_with_value_one(self.bin_names)
        # could be any of the histograms, we just need the binning
        for bin_number in cos_zenith_is_one_bins:
            x_bin_end = hist.GetBinLowEdge(bin_number) + hist.GetBinWidth(bin_number)+1
            line = ROOT.TLine(x_bin_end, ROOT.gPad.GetUymin(), x_bin_end, ROOT.gPad.GetUymax())
            line.SetLineColor(ROOT.kBlack)
            line.SetLineWidth(2)
            line.SetLineStyle(2)
            self.bin_lines.append(line)

    def add_momentum_lables(self):
        positions = self.bin_manager.get_indices_with_value_one(self.bin_names)
        for i in range(len(positions)):
            x_position = i*(positions[0]+1)+(1+positions[0])/2
            momentum = [float(x) for x in self.bin_names[i*(positions[0]+1)].strip('[]').split()]
            result_string = ('[10^{'+'{:.1f}'.format(momentum[0])+'}, 10^{'+'{:.1f}'.format(momentum[1])+'}]')

            text1 = ROOT.TLatex()
            #text1.SetTextSize(0.035)
            text1.SetTextAlign(22)
            text1.DrawLatex(x_position, ROOT.gPad.GetUymax()*0.95, result_string)

            text2 = ROOT.TLatex()
            #text2.SetTextSize(0.035)
            text2.SetTextAlign(22)
            text2.DrawLatex(x_position, ROOT.gPad.GetUymax()*0.89, 'MeV/c')

            x_axis = self.data_hist.GetXaxis() #list(self.histograms.values())[-1].GetXaxis()
            x_min = x_axis.GetXmin()
            x_max = x_axis.GetXmax()
            left_margin = ROOT.gPad.GetLeftMargin()
            right_margin = ROOT.gPad.GetRightMargin()
            ndc_x = left_margin + (x_position - x_min) / (x_max - x_min) * (1 - left_margin - right_margin)

            text3 = ROOT.TLatex()
            text3.SetNDC()
            text3.SetTextAlign(22)
            text3.DrawLatex(ndc_x, 0.03,  "cos#theta_{z}")

    def reset_histograms(self):
        for hist in self.histograms.values():
            hist.Reset()
        self.data_hist.Reset()
    def add_reaction_to_hist(self, reaction):
        counts, weights = reaction.get_bin_counts(self.bin_indices)
        for i in range(len(self.bin_names)):
            self.histograms[reaction.name].SetBinContent(int(i+1), float(counts[i]*weights[i]))

    def set_total_MC_as_data_hist(self):
        self.data_hist.Reset()
        for h in self.histograms.values():
            for i in range(len(self.bin_names)):
                self.data_hist.SetBinContent(int(i+1), h.GetBinContent(i+1) + self.data_hist.GetBinContent(i+1))

    def set_SK_true_data_as_data_hist(self):
        self.data_hist.Reset()
        df= pd.read_csv(base_dir_path+'../data/sk_data/sk_2023_Data.txt', delim_whitespace=True)
        for local_bin_idx, global_bin_idx in enumerate(self.bin_indices):
            self.data_hist.SetBinContent(int(local_bin_idx+1), int(df.iloc[global_bin_idx]))

    def create_stack(self):
        self.stack  = ROOT.THStack(self.name, self.name)
        self.legend = ROOT.TLegend(0.2, 0.87, 0.9, 0.91)

        # let's have the legend in order reversed to the stacked histogram
        for name, hist in self.histograms.items():
            self.stack.Add(hist)

        #self.legend.AddEntry(self.data_hist, "Data", "E1P")
        for name, hist in reversed(list(self.histograms.items())):
            self.legend.AddEntry(hist, self.reaction_names_to_latex[name], "f")

    def plot(self, out_hist=None):
        # Set font to 'Serif' for all labels and text
        ROOT.gStyle.SetLabelFont(42, "xyz")  # Set font for axis labels
        ROOT.gStyle.SetLabelFont(42, "title")  # Set font for axis titles
        ROOT.gStyle.SetLabelFont(42, "legend")  # Set font for legend
        ROOT.gStyle.SetTextFont(42)  # Set font for text, such as TLatex

        # # Apply 'Serif' font to all labels and text
        # ROOT.gStyle.SetLabelSize(0.04, "xyz")  # Set label size for axis labels
        # ROOT.gStyle.SetLabelSize(0.04, "title")  # Set label size for axis titles
        # ROOT.gStyle.SetLabelSize(0.04, "legend")  # Set label size for legend
        ROOT.gStyle.SetTextSize(0.04)  # Set text size for text, such as TLatex

        c = ROOT.TCanvas()
        c.cd()
        hist = self.data_hist
        if out_hist is not None:
            hist = out_hist

        c.SetBottomMargin(0.1)
        c.SetTopMargin(0.15)

        self.stack.Draw("hist")
        #hist.Draw("E1P same")
        # hist.Draw("E1P same")
        # Update the Y-axis range to fit the highest value
        max_y = max(
            hist.GetBinContent(hist.GetMaximumBin()) + hist.GetBinError(hist.GetMaximumBin()),
            self.stack.GetMaximum()
        )
        self.stack.SetMaximum(max_y * 1.2)  # Add some extra space at the top
        self.stack.SetMinimum(0)
        ROOT.gStyle.SetTitleY(.99)

        c.Update()
        c.SetLeftMargin(0.14)
        c.SetRightMargin(0.02)

        self.create_zenith_angle_edge(hist)
        if len(self.bin_lines)>0:
            for line in self.bin_lines[:-1]:
                line.Draw()
        self.legend.SetNColumns(7)
        self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)
        self.legend.Draw()
        self.stack.GetYaxis().SetTitle("Counts")
        self.stack.GetYaxis().CenterTitle(True)

        if self.name in self.special_samples:
            self.stack.GetXaxis().SetTitle("log_{10} [p (MeV/c)]")
            self.stack.GetXaxis().CenterTitle(True)
            self.stack.GetXaxis().SetTitleOffset(1.3)

        self.stack.GetYaxis().SetTitleSize(0.04)
        self.stack.GetYaxis().SetLabelSize(0.04)
        self.stack.GetXaxis().SetLabelSize(0.06)
        self.stack.GetXaxis().LabelsOption("h")


        self.add_momentum_lables()

        return c

    def calc_aggregated_hist(self, hist):
        if len(self.bin_names)<10:
            return hist

        new_hist = ROOT.TH1F('1D '+hist.GetName(),'1D '+hist.GetName(), 10, 0, 10)

        new_hist.SetFillColor(hist.GetFillColor())
        new_hist.SetLineWidth(hist.GetLineWidth())
        new_hist.SetLineColor(hist.GetLineColor())
        #self.format_histogram_axis(new_hist)

        folds = int(len(self.bin_names)/10)
        for i in range(10):
            for j in range(folds):
                new_hist.SetBinContent(i+1, hist.GetBinContent(folds*j+i+1)+new_hist.GetBinContent(i+1))

        return new_hist

    def plot_1D(self, out_hist=None):
        # Set font to 'Serif' for all labels and text
        ROOT.gStyle.SetLabelFont(42, "xyz")  # Set font for axis labels
        ROOT.gStyle.SetLabelFont(42, "title")  # Set font for axis titles
        ROOT.gStyle.SetLabelFont(42, "legend")  # Set font for legend
        ROOT.gStyle.SetTextFont(42)  # Set font for text, such as TLatex

        ROOT.gStyle.SetTextSize(0.04)  # Set text size for text, such as TLatex

        c = ROOT.TCanvas()
        c.cd()

        new_histograms = [self.calc_aggregated_hist(h) for h in list(self.histograms.values())]
        new_data_hist  = self.calc_aggregated_hist(self.data_hist)
        # for h in new_histograms:
        #     h.Draw("HIST")

        self.stack  = ROOT.THStack(self.name, self.name)
        self.legend = ROOT.TLegend(0.1, 0.87, 0.9, 0.91)

        # let's have the legend in order reversed to the stacked histogram
        for h in new_histograms:
            self.stack.Add(h)

        #self.legend.AddEntry(self.data_hist, "Data", "E1P")
        names = list(self.histograms.keys())
        cnt = 0
        for hist in reversed(new_histograms):
            self.legend.AddEntry(hist, self.reaction_names_to_latex[names[-1-cnt]], "f")
            cnt+=1

        c.Update()
        # hist = self.data_hist
        # if out_hist is not None:
        #     hist = out_hist
        #
        c.SetBottomMargin(0.1)
        c.SetTopMargin(0.15)

        self.stack.Draw("hist")
        new_data_hist.DrawCopy("E1P same")
        # # Update the Y-axis range to fit the highest value
        max_y = max(
            new_data_hist.GetBinContent(new_data_hist.GetMaximumBin()) + new_data_hist.GetBinError(new_data_hist.GetMaximumBin()),
            self.stack.GetMaximum()
        )
        self.stack.SetMaximum(max_y * 1.2)  # Add some extra space at the top
        self.stack.SetMinimum(0)
        ROOT.gStyle.SetTitleY(.99)

        c.Update()

        self.legend.SetNColumns(7)
        self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)
        self.legend.Draw()
        self.stack.GetYaxis().SetTitle("Counts")
        self.stack.GetYaxis().CenterTitle(True)

        if self.name in self.special_samples:
            self.stack.GetXaxis().SetTitle("log_{10} [p (MeV/c)]")
            self.stack.GetXaxis().CenterTitle(True)
            self.stack.GetXaxis().SetTitleOffset(1.3)

        self.stack.GetYaxis().SetTitleSize(0.04)
        self.stack.GetYaxis().SetLabelSize(0.04)
        self.stack.GetXaxis().SetLabelSize(0.06)
        self.stack.GetXaxis().LabelsOption("h")

        return c, new_histograms, new_data_hist

class AnaMaster:
    """
    A class to control how to load, and access the data.
    """
    def __init__(self, data_folder = base_dir_path+'../data/unoscillated', scale_to_HK=False, dont_replace_names=False, binning_file=base_dir_path+"../data/sk_2023_BinInfo.txt", mask_bins=True):
        self.scale_to_HK = scale_to_HK
        self.config_file=base_dir_path+'../config/chi2_config.yaml'
        #self.config_file='/afs/cern.ch/work/c/cjesus/private/EarthTomo/config/chi2_config.yaml'
        
        with open(self.config_file,'r') as f:
            cfg=yaml.safe_load(f)
            #self.OscProb = oscProb(cfg)
            self.OscProb = oscMaster(cfg)
            #self.OscProb._load_barger_prop()

        self.bin_manager = BinManager(binning_file)
        self.reactions = self.process_unoscillated_data(data_folder)
        self.reaction_names = [reaction.name for reaction in self.reactions]

        def unique_ordered_with_counts(arr):
            arr = np.asarray(arr)
            order = np.argsort(np.unique(arr, return_index=True)[1])
            return np.unique(arr)[order], np.unique(arr, return_counts=True)[1][order]

        #samples, size = np.unique(self.bin_manager.bin_names, return_counts=True)
        samples, size = unique_ordered_with_counts(self.bin_manager.bin_names)
        self.sample_bins = {}
        self.samples = []
        self.exposure = {}

        self.dont_replace_names = dont_replace_names
        self.sample_name_to_formated_name = self.replace_names(samples)


        for idx, sample in enumerate(samples):
            sample_bins = np.where(self.bin_manager.bin_names==sample)[0]
            self.sample_bins[sample] = sample_bins
            ana_sample = AnaSample(self.sample_name_to_formated_name[sample], \
                                   self.bin_manager, size[idx], sample_bins, self.reaction_names)
            self.samples.append(ana_sample)

        print('MASK BINS: ', mask_bins)
        if mask_bins:
            self.mask_bins()
        self.scale_statistics()
        self.fill_histograms()

        self.flux_manager = fluxManager.FluxManager()
        self.flux_df = self.flux_manager.df
        self.cosZ_to_fluxbin = self.flux_manager.cosZ_bin_dic
        self.fluxbin_to_cosZ = self.flux_manager.bin_cosZ_dic
        flux_E_values = sorted(np.unique(self.flux_df.Enu.to_numpy()))

        self.E_bins = []
        for i in range(len(flux_E_values)-1):
            self.E_bins.append([flux_E_values[i], flux_E_values[i+1]])

        cosZ_bins = list(self.cosZ_to_fluxbin.values())

        for reac in self.reactions:
            for ana_bin in reac.ana_bins:
                ana_bin.set_flux_info(self.flux_manager, self.E_bins, cosZ_bins)
    def replace_names(self, names):
        replaced_names = {}
        replacements = {
            "sk1-3_": "",
            "sk1-5_": "",
            "sk4-5_": "",
            "multigev": "Multi-GeV ",
            "subgev": "Sub-GeV ",
            "_elike": " e-like",
            "_mulike": "#mu-like",
            "_nuebarlike": "#bar{#nu}_{e}-like",
            "_nuelike": "#nu_{e}-like",
            "_numubarlike": "#bar{#nu}_{#mu}-like",
            "_numulike": "#nu_{#mu}-like",
            "_0decaye": " 0 d.e.",
            "_1decaye": " 1 d.e.",
            "_2decaye": " 2 d.e.",
            "upmu_thru": "upmu",
            "upmu": "Up-#mu",
            "_stop": " Stopping",
            "_thru": " Through-going",
            "pc": "PC ",
            "_showering": " Showering",
            "_nonshowering": " Non-Showering",
            "_multiring": "Multi-Ring ",
            "1ring_ncpi0": "1-ring #pi^{0}-like",
            "2ring_ncpi0": "2-ring #pi^{0}-like",
            "ncpi0": "#pi^{0}",
            "_1ring": "",
            "_other": " Other",
            "fc_": "",
            "_2ring": "2 Ring",
            "_1neutron": " 1 n",
            "_0neutron": " 0 n",
        }

        if self.dont_replace_names:
            replacements={}            

        for name in names:
            original_name = name
            for old, new in replacements.items():
                name = name.replace(old, new)
            if "sk1-3_" in original_name:
                name += '*'
                self.exposure[name] = 6510.9-3705.4
            elif "sk4-5_" in original_name:
                self.exposure[name] = 3705.4
            else:
                self.exposure[name] = 6510.9

            replaced_names[original_name] = name

        print(self.exposure)
        print(self.exposure.values())
        return replaced_names

    def fill_histograms(self):
        for sample in self.samples:
            sample.reset_histograms()
            for reaction in self.reactions:
                if 'Total' not in reaction.name:
                    sample.add_reaction_to_hist(reaction)

            sample.set_total_MC_as_data_hist()
            #sample.set_SK_true_data_as_data_hist()

    def osc_weight_all(self):
        for r in self.reactions:
            if 'Total' not in r.name:
                r.apply_osc_weights()

    def resample_from_quantiles(self):
        # the concept here has to be something like iterating in every reaction in every sample and then resample from the quantiles in every bin.
        # then we get back the collection of sampled E and theta and we combined them scaling them up by the counts on the bin.
        # probably the best is to store the full collection of toys in every SK bin and then post-process that into new SK-like files.
        # let's do something like a map from bin ID to a collection of [E,costheta].
        # and then we pass that to a function that takes care of creating the new bins (maybe?).

        for r in self.reactions:
            for s in self.samples:
                for b_idx in s.bin_indices:
                    r.ana_bins[b_idx].resample()  

    def mask_bins(self):
        for r in self.reactions:
            for s in self.samples:
                if not valid_sample_condition(s.name):
                    for b_idx in s.bin_indices:
                        r.ana_bins[b_idx].counts=0
                        r.ana_bins[b_idx].ignore_bin=True

    def scale_statistics(self):
        for r in self.reactions:
            for s in self.samples:
                exposure = 1
                HK_factor = 1
                if self.scale_to_HK:
                    exposure = 6510.9/self.exposure[s.name] # 20 years
                    HK_factor = 8                    

                for b_idx in s.bin_indices:
                    r.ana_bins[b_idx].counts=r.ana_bins[b_idx].counts*exposure*HK_factor
                    if r.ana_bins[b_idx].counts:
                        print(s.name, "----", self.exposure[s.name], "----", r.ana_bins[b_idx].counts)


    def process_unoscillated_data(self, data_folder):
        reactions = []
        for file in os.listdir(data_folder):
            if file.endswith('.txt'):
                reaction = AnaReaction(os.path.join(data_folder, file), self.bin_manager, self.OscProb)
                reactions.append(reaction)

        return reactions
