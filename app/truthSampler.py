import numpy as np

def get_bin_as_float(input_bin):
    return [float(x) for x in input_bin.strip('[').strip(']').split(' ')]

def extend_binning(bins):
    new_bins = []
    for i in range(len(bins)):
        new_bins.append(bins[i])
        if i<len(bins)-1:
            new_bins.append((bins[i]+bins[i+1])/2) # add one bin in the middle...
    return new_bins

def create_bins(bins_x, bins_y):
    new_bins = []
    for i in range(len(bins_x)-1):
        for j in range(len(bins_y)-1):
            new_bins.append([bins_x[i], bins_x[i+1], bins_y[j], bins_y[j+1]])

    return new_bins

def dict_to_structured_file(data_dict, output_file):
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"{'SampleName':<44} {'LogP_lo':>9} {'LogP_hi':>9} {'CosZ_lo':>9} {'CosZ_hi':>9}\n")
        
        for key, value_list in data_dict.items():
            for value in value_list:
                # Remove asterisk from key if present
                sample_name = key.rstrip('*')
                
                # Write each line
                f.write(f"{sample_name:<44} {value[0]:9.2f} {value[1]:9.2f} {value[2]:9.2f} {value[3]:9.2f}\n")


class NewBin:
    def __init__(self, boundaries):
        self.boundaries = boundaries
        self.tZs = []  # collection of true angles
        self.tEs = []  # collection of true energies
        self.rZs = []  # collection of true angles
        self.rEs = []  # collection of true energies
        self.Cs = []   # collection of bin counts
        self.b_indices = [] # just to know how many original 'bins' true distributions contribute

        self.cnts = None
        self.avg_E = None
        self.rms_E = None
        self.E_quantiles = None
        self.avg_Z = None
        self.rms_Z = None
        self.Z_quantiles = None

    def evt_is_in_the_bin(self, Z, E, debug=False):
        xlo = self.boundaries[0]
        xhi = self.boundaries[1]
        ylo = self.boundaries[2]
        yhi = self.boundaries[3]

        # print(10**xlo, 10**xhi, ylo, yhi)

        if debug:
            print(Z,E)
            print(10**xlo, 10**xhi, ylo, yhi, 10**xlo <= E and E <= 10**xhi and ylo <= Z and Z <= yhi)
            print('\n')

        # note that here I am converting the bin edges to log10 for E
        if 10**xlo <= E and E <= 10**xhi and ylo <= Z and Z <= yhi:
            return True

        return False

    def add(self, tZ, tE, rZ, rE, C, b_idx, nsampled_evts):
        self.tZs.append(tZ)
        self.tEs.append(tE)
        self.rZs.append(rZ)
        self.rEs.append(rE)
        self.Cs.append(C)
        self.b_indices.append(b_idx)
        self.nsampled_evts = nsampled_evts

    def show_1D_tZ(self):
        plt.hist(self.tZs, weights = self.Cs, bins=np.linspace(self.boundaries[2],  self.boundaries[3], 20))

    def show_1D_rZ(self):
        plt.hist(self.rZs, weights = self.Cs, bins=np.linspace(min(self.rZs),  max(self.rZs), 20))

    def show_1D_tE(self):
        plt.hist(self.tEs, weights = self.Cs, bins=np.linspace(min(self.tEs),  max(self.tEs), 20))
        plt.axvline(10**self.boundaries[0], c='r')
        plt.axvline(10**self.boundaries[1], c='r')

    def show_1D_rE(self):
        plt.hist(self.rEs, weights = self.Cs, bins=np.linspace(min(self.rEs),  max(self.rEs), 20))
        plt.axvline(10**self.boundaries[0], c='r')
        plt.axvline(10**self.boundaries[1], c='r')


    def get_Z_quantile(self):
        # Ensure the data is in the correct format
        tZs = np.array(self.tZs)
        Cs = np.array(self.Cs)

        # Sort the data and weights by the data values
        sorted_indices = np.argsort(tZs)
        sorted_data = tZs[sorted_indices]
        sorted_weights = Cs[sorted_indices]

        # Compute the cumulative sum of the weights and normalize it
        cumulative_weights = np.cumsum(sorted_weights)
        cumulative_weights /= cumulative_weights[-1]

        # Quantiles to find
        quantiles = [2.3, 15.9, 50, 84.1, 97.7]
        quantiles = np.array(quantiles) / 100  # Convert percentages to fractions

        # Interpolation to find the quantile values
        quantile_values = np.interp(quantiles, cumulative_weights, sorted_data)

        return quantile_values

    def get_E_quantile(self):
        # Ensure the data is in the correct format
        tEs = np.array(self.tEs)
        Cs = np.array(self.Cs)

        # Sort the data and weights by the data values
        sorted_indices = np.argsort(tEs)
        sorted_data = tEs[sorted_indices]
        sorted_weights = Cs[sorted_indices]

        # Compute the cumulative sum of the weights and normalize it
        cumulative_weights = np.cumsum(sorted_weights)
        cumulative_weights /= cumulative_weights[-1]

        # Quantiles to find
        quantiles = [2.3, 15.9, 50, 84.1, 97.7]
        quantiles = np.array(quantiles) / 100  # Convert percentages to fractions

        # Interpolation to find the quantile values
        quantile_values = np.interp(quantiles, cumulative_weights, sorted_data)

        return quantile_values

    def fill_zeros(self):
        self.cnts = 0.
        self.avg_E = 0.
        self.rms_E = 0.
        self.E_quantiles = [0., 0., 0., 0., 0.]
        self.avg_Z = 0.
        self.rms_Z  = 0.
        self.Z_quantiles = [0., 0., 0., 0., 0.]
        #self.print_filled_values()
    
    def fill_values(self):
        def weighted_stats(values, weights):
            """Calculate weighted mean and standard deviation."""
            average = np.average(values, weights=weights)
            
            # Weighted variance
            variance = np.average((values - average)**2, weights=weights)
            
            # Weighted standard deviation
            std_dev = np.sqrt(variance)
            
            return average, std_dev

        self.cnts = np.sum(self.Cs)/self.nsampled_evts
        self.avg_E, self.rms_E = weighted_stats(self.rEs, self.Cs)
        self.E_quantiles = self.get_E_quantile()

        self.avg_Z, self.rms_Z = weighted_stats(self.rZs, self.Cs)
        self.Z_quantiles = self.get_Z_quantile()

    def print_filled_values(self):
        print("Filled Values:")
        print(f"{'Counts:':<30} {self.cnts:.6f}")
        print(f"{'Average Energy (log10):':<30} {self.avg_E:.6f}")
        print(f"{'RMS Energy (log10):':<30} {self.rms_E:.6f}")
        print("Energy Quantiles (log10):")
        for i, q in enumerate([2.3, 15.9, 50.0, 84.1, 97.7]):
            print(f"  {q}%:{self.E_quantiles[i]:>20.6f}")
        print(f"{'Average CosZ (log10):':<30} {self.avg_Z:.6f}")
        print(f"{'RMS CosZ (log10):':<30} {self.rms_Z:.6f}")
        print("CosZ Quantiles:")
        for i, q in enumerate([2.3, 15.9, 50.0, 84.1, 97.7]):
            print(f"  {q}%:{self.Z_quantiles[i]:>20.6f}")

    def get_bin_values_for_table(self):
        return self.cnts, self.avg_E, self.rms_E, *self.E_quantiles, self.avg_Z, self.rms_Z, *self.Z_quantiles 

def create_bin_content_table(filename, list_of_new_bin_objects):
    header = [
        "Counts", "EnergyAvg", "EnergyRMS",
        "EnergyQuantile2.3Percent", "EnergyQuantile15.9Percent",
        "EnergyQuantile50.0Percent", "EnergyQuantile84.1Percent",
        "EnergyQuantile97.7Percent", "CosZAvg", "CosZRMS",
        "CosZQuantile2.3Percent", "CosZQuantile15.9Percent",
        "CosZQuantile50.0Percent", "CosZQuantile84.1Percent",
        "CosZQuantile97.7Percent"
    ]
    
    # Calculate the width for each column
    column_width = 26
    
    # Specify the number of spaces before the first column
    first_column_padding_header = 20
    first_column_padding = 14
    
    with open(filename, 'w') as f:
        # Write the header
        f.write(" " * first_column_padding)  # Add padding before the first column
        for i, item in enumerate(header):
            if i == 0:
                # For the first item, use left alignment and reduced width
                f.write(f"{item:<{column_width - first_column_padding_header}}")
            else:
                # For all other items, use right alignment and full width
                f.write(f"{item:>{column_width}}")
        f.write("\n")
        
        # Write the data for each bin
        for bin_object in list_of_new_bin_objects:
            bin_values = bin_object.get_bin_values_for_table()
            f.write(" " * first_column_padding)  # Add padding before the first column
            for i, value in enumerate(bin_values):
                if i == 0:
                    # For the first item, use left alignment and reduced width
                    f.write(f"{value:<{column_width - first_column_padding}.6e}")
                else:
                    # For all other items, use right alignment and full width
                    f.write(f"{value:>{column_width}.6e}")
            f.write("\n")
    
    print(f"File '{filename}' has been created.")


def make_all_binned_tables(anaMaster, binning_scheme, outpath, tag):
    for r in anaMaster.reactions:
        empty_cnt = 0
        print('----------------------')
        print(r.name)
        print('----------------------')
        sample_to_new_bin_objects = []
        for s_idx, s in enumerate(anaMaster.samples):
            #print(s.name)
            list_of_new_bin_objects = []
            #print(s.bin_indices)  # check that the bins are in consecutive order!!!
            bin_boundaries = binning_scheme[s.name]
            for b in bin_boundaries:
                nbin = NewBin(b)
                list_of_new_bin_objects.append(nbin)
    
            for b_idx in s.bin_indices:
                true_Zs = r.ana_bins[b_idx].sampled_Zs
                true_Es = r.ana_bins[b_idx].sampled_Es*1000
                Cs = np.ones(len(r.ana_bins[b_idx].sampled_Zs))*r.ana_bins[b_idx].counts
    
                reco_Zs = [apply_Z_smearing(x, tag) for x in true_Zs]
                reco_Es = [apply_E_smearing(x, tag) for x in true_Es]

                cnt_bad = 0
    
                for rZ, rE, tZ, tE, C in zip(reco_Zs, reco_Es, true_Zs, true_Es, Cs):
                    flag = False
                    for b in list_of_new_bin_objects:
                        if b.evt_is_in_the_bin(rZ,rE):
                            flag = True
                            nsampled_evts = len(true_Zs)
                            b.add(tZ, tE/1000, rZ, rE/1000, C, b_idx, nsampled_evts)
    
                    if not flag:
                        cnt_bad +=1
    
            for idx, b in enumerate(list_of_new_bin_objects):
                if len(b.Cs):
                    b.fill_values()
                else:
                    b.fill_zeros()
    
            # we need first to do things per sample to avoid cross-sample migration, then we flatten
            for item in list_of_new_bin_objects:
                sample_to_new_bin_objects.append(item)
    
        print('Empty bins: ', 100*empty_cnt/len(sample_to_new_bin_objects), ' %')
        fname = outpath+'sk_2023_'+r.name+'NoOsc.txt'
        create_bin_content_table(fname, sample_to_new_bin_objects)


def smear_costheta(costheta, smearing_degrees=10):
    # Convert costheta to a 3D unit vector (assuming y=0)
    theta = np.arccos(costheta)
    x = np.sin(theta)
    z = costheta
    vector = np.array([x, 0, z])

    # Generate a random rotation axis
    rotation_axis = np.random.randn(3)
    rotation_axis /= np.linalg.norm(rotation_axis)

    # Generate a random rotation angle from a normal distribution
    # The standard deviation is set to achieve approximately 10% smearing
    angle_std = smearing_degrees * np.pi / 180  # Convert to radians
    rotation_angle = np.random.normal(0, angle_std)

    # Create the rotation matrix using Rodrigues' rotation formula
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    rotation_matrix = (np.eye(3) + np.sin(rotation_angle) * K + 
                       (1 - np.cos(rotation_angle)) * np.matmul(K, K))

    # Apply the rotation to the vector
    rotated_vector = np.dot(rotation_matrix, vector)

    # Calculate the new costheta
    smeared_costheta = rotated_vector[2]  # z-component is costheta
    angle_between = np.arccos(np.dot(vector, rotated_vector))

    return smeared_costheta

def smear_energy(energy, energy_resolution=1):
    return energy*np.random.normal(1,energy_resolution)

def apply_Z_smearing(x, tag):
    if tag == 'B':
        return smear_costheta(x, 10)
    if tag == 'C':
        return smear_costheta(x, 20)
    if tag == 'F':
        return smear_costheta(x, 10)
    if tag == 'G':
        return smear_costheta(x, 20)
    if tag == 'H':
        return smear_costheta(x, 10)
    if tag == 'I':
        return smear_costheta(x, 10)
    if tag == 'J':
        return smear_costheta(x, 20)

    return x

def apply_E_smearing(x, tag):
    if tag == 'D':
        return smear_energy(x, 0.3)
    if tag == 'E':
        return smear_energy(x, 0.5)
    if tag == 'F':
        return smear_energy(x, 0.3)
    if tag == 'G':
        return smear_energy(x, 0.5)
    if tag == 'H':
        return smear_energy(x, 0.2)
    if tag == 'I':
        return smear_energy(x, 0.1)
    if tag == 'J':
        return smear_energy(x, 0.3)


    # if tag == 'F':
    #     return smear_energy(x, 0.3)
    # if tag == 'G':
    #     return smear_energy(x, 1)

    return x


def process_FDS(anaMaster, tag, binning_scheme):
    path = "/Users/cjesus/Documents/HKnuTomo/fake_data/FDS_" + tag + "/"
    fname_bins = path + "/binning.txt"
    dict_to_structured_file(binning_scheme, fname_bins)
    print(f"File '{fname_bins}' has been created.")
    fname_content = path + "unoscillated/"
    make_all_binned_tables(anaMaster, binning_scheme, fname_content, tag)


def create_bins_using_standard_binning(anaMaster):

    sample_to_new_bin_boundaries = {}
    for s in anaMaster.samples:
            positions = anaMaster.bin_manager.get_indices_with_value_one(s.bin_names)
        
            if positions:
                xmin = get_bin_as_float(s.bin_names[0])[0]
                ymin = get_bin_as_float(s.bin_names[0])[2]
                xmax = get_bin_as_float(s.bin_names[-1])[1]
                ymax = get_bin_as_float(s.bin_names[-1])[3]

                Zs = [ymin]
                for idx in range(positions[0]):
                    Zs.append(get_bin_as_float(s.bin_names[idx])[3])
                Zs.append(ymax)

                Es = [xmin]
                for idx in positions:
                    Es.append(get_bin_as_float(s.bin_names[idx])[1])

                sample_to_new_bin_boundaries[s.name] = create_bins(Es, Zs)

            else:
                Es = []
                for i in range(len(s.bin_names)):
                    Es.append(get_bin_as_float(s.bin_names[i])[0])
                Es.append(get_bin_as_float(s.bin_names[-1])[1])

                Zs = [-1,1]
                sample_to_new_bin_boundaries[s.name] = create_bins(Es, Zs)
    
    
    return sample_to_new_bin_boundaries

def create_bins_using_narrow_binning(anaMaster):

    def get_bin_values(x, name):
        if 'pi0' in name:
            return (2.,3., 3.5, 4.)

        if 'upmu_stop' in name:
            return (3.0, 3.6, 3.8, 4., 4.5, 5.5)

        if 'subgev' in name and 'nuelike' in name:
            return (2.6, 2.9, 3.0, 3.2, 3.5)

        if 'pc_thru' in name:
            return (3.0, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 5.5)

        if 'multigev_1ring_nuebarlike' in name:
            return (3.0, 3.3, 3.6, 4.5)

        if 'multigev_1ring_nuelike' in name:
            return (3.0, 3.6, 4.5)

        if 'subgev_1ring_mulike_2decaye' in name or 'subgev_1ring_elike_1decaye' in name:
            return (2.5, 2.9, 3.2, 4)

        if 'multigev_multiring_nuelike' in name or 'multigev_multiring_nuebarlike' in name:
            return (3.0, 3.5, 3.7, 5.0)

        if 'multigev_multiring_mulike' in name:
            return (3.0, 3.2, 3.3, 3.4, 3.5, 3.6, 3.8, 4.0, 5.0)

        if 'multigev_multiring_other' in name:
            return (3.3, 3.8, 4.0, 5.0)

        if 'showering' in name:
            if 'non' in name:
                return (3.0, 4.5, 8.0)
            else:
                return (3.0, 8.0)

        if 'multigev_1ring' in name and 'mu' in name:
            return (3.0, 3.4, 5.0)

        if 'pc_stop' in name:
            return (2.0, 3.4, 5.0)

        if 'subgev' in name and 'nuebarlike' in name:
            return (2.0, 2.3, 2.5, 2.7, 2.8, 3.0, 3.2, 4)

        return (2.0, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.5)

    sample_to_new_bin_boundaries = {}
    for s in anaMaster.samples:
            print(s.name)
            positions = anaMaster.bin_manager.get_indices_with_value_one(s.bin_names)
        
            if positions:
                xmin = get_bin_as_float(s.bin_names[0])[0]
                ymin = get_bin_as_float(s.bin_names[0])[2]
                xmax = get_bin_as_float(s.bin_names[-1])[1]
                ymax = get_bin_as_float(s.bin_names[-1])[3]

                # Zs = [ymin]
                # for idx in range(positions[0]):
                #     Zs.append(get_bin_as_float(s.bin_names[idx])[3])
                # Zs.append(ymax)

                Zs = np.linspace(-1,0,41)

                print('Zs: ', Zs)

                Es = [xmin]
                for idx in positions:
                    Es.append(get_bin_as_float(s.bin_names[idx])[1])

                Es = get_bin_values(tuple(Es), s.name)
                sample_to_new_bin_boundaries[s.name] = create_bins(Es, Zs)

            else:
                Es = []
                for i in range(len(s.bin_names)):
                    Es.append(get_bin_as_float(s.bin_names[i])[0])
                Es.append(get_bin_as_float(s.bin_names[-1])[1])
                Es = get_bin_values(tuple(Es), s.name)

                Zs = [-1,1]
                sample_to_new_bin_boundaries[s.name] = create_bins(Es, Zs)
                print('special z-binning: ', s.name)

    return sample_to_new_bin_boundaries
