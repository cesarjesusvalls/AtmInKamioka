from app.path import *
import pandas as pd
import re

class FluxManager:
    def __init__(self):
        self.cosZ_bin_dic = {}
        self.bin_cosZ_dic = {}
        self.bin_id = -1
        self.data = {'Enu': [], 'NuMu': [], 'NuMubar': [], 'NuE': [], 'NuEbar': [], 'cosZ': [], 'cosZ_bin': [], 'phi_Az': []}
        self.load_flux_data()
        self.preprocess_flux_data()
    @staticmethod
    def extract_cosz_phi_values(input_string):
        # Regular expression pattern to extract numerical values
        pattern = r"[-+]?\d*\.\d+|\d+"
        matches = re.findall(pattern, input_string)

        # Extracted values
        cosZ_range = [float(matches[0]), float(matches[1])]
        phi_Az_range = [float(matches[2]), float(matches[3])]

        return cosZ_range, phi_Az_range

    def load_flux_data(self, filename=base_dir_path+'../data/kam-nu-20-01-mtn.d'):
        with open(filename, 'r') as file:
            lines = file.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("average flux in"):
                    # Extract cosZ and phi_Az from the caption
                    cosZ, phi_Az = self.extract_cosz_phi_values(line)
                    i += 2  # Skip the next line as it contains column headers
                    self.bin_id += 1
                    self.cosZ_bin_dic[self.bin_id] = cosZ
                    self.bin_cosZ_dic[tuple(cosZ)] = self.bin_id
                else:
                    # Extract binned information
                    for j in range(101):
                        row = lines[i+j].strip().split()
                        self.data['Enu'].append(float(row[0]))
                        self.data['NuMu'].append(float(row[1]))
                        self.data['NuMubar'].append(float(row[2]))
                        self.data['NuE'].append(float(row[3]))
                        self.data['NuEbar'].append(float(row[4]))
                        self.data['cosZ'].append(cosZ)
                        self.data['cosZ_bin'].append(self.bin_id)
                        self.data['phi_Az'].append(phi_Az)
                    i += 101

        self.df = pd.DataFrame(self.data)


    def preprocess_flux_data(self):
        self.flux_data = {}

        # Group data by cosZ_bin
        grouped_by_cosZ = self.df.groupby('cosZ_bin')

        for cosZ_bin, group in grouped_by_cosZ:
            self.flux_data[cosZ_bin] = {}

            # Sort by Enu
            group_sorted = group.sort_values(by='Enu')

            # Iterate over Enu bins
            for i in range(len(group_sorted) - 1):
                lo_bin = group_sorted.iloc[i]
                hi_bin = group_sorted.iloc[i + 1]
                self.flux_data[cosZ_bin][i] = {'Enu': lo_bin['Enu'], **lo_bin.drop('Enu').to_dict()}
                self.flux_data[cosZ_bin][i + 1] = {'Enu': hi_bin['Enu'], **hi_bin.drop('Enu').to_dict()}
