import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, IntegrationWarning
from scipy.optimize import minimize
import warnings
from itertools import product
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

class EarthModel:
    def __init__(self, cfg:dict, change_type:str, change_value:float, verbose=False, core_mode='single', use_average_density=True):
        self.change_type = change_type
        self.change_value = change_value
        self.verbose = verbose
        self.core_mode = core_mode  # 'single' or 'independent'
        self.use_average_density = use_average_density
        
        self.R = 6371  # Earth radius in km
        self.original_boundaries = {
            'IC_OC': 1221.5,
            'OC_IM': 3480.0,
            'IM_OM': 5701.0
        }

        self.modified_boundaries = self.original_boundaries.copy()
        self.original_densities = {
            'IC': self.density_IC,
            'OC': self.density_OC,
            'IM': self.density_IM,
            'OM': self.density_OM
        }
        self.modified_densities = self.original_densities.copy()
        self.layer_factors = {}

        self.PREM_file = cfg['data']['filepath']

        self.load_data()

    def calculate_average_density(self, layer):
        boundaries = self.modified_boundaries
        if layer == 'IC':
            r_start, r_end = 0, boundaries['IC_OC']
        elif layer == 'OC':
            r_start, r_end = boundaries['IC_OC'], boundaries['OC_IM']
        elif layer == 'IM':
            r_start, r_end = boundaries['OC_IM'], boundaries['IM_OM']
        elif layer == 'OM':
            r_start, r_end = boundaries['IM_OM'], self.R
        else:
            raise ValueError(f"Invalid layer: {layer}")
        
        density_func = self.modified_densities[layer]
        
        # Calculate average density using numerical integration
        num_points = 20
        r_values = np.linspace(r_start, r_end, num_points)
        densities = [density_func(r) for r in r_values]
        average_density = np.mean(densities)
        
        return average_density

    def impose_hydrostatic_equilibrium(self, change_type):
        boundaries = list(self.modified_boundaries.values())
        layers = ['IC', 'OC', 'IM', 'OM']
        layer_idx = layers.index(change_type)

        for i in range(len(layers)):
            if i < layer_idx:
                if self.use_average_density:
                    upper_density = self.calculate_average_density(layers[i+1])
                    lower_density = self.calculate_average_density(layers[i])
                else:
                    upper_density = self.density(boundaries[i] + 1e-6)
                    lower_density = self.density(boundaries[i] - 1e-6)
                
                if upper_density > lower_density:
                    factor = upper_density / lower_density
                    self.change_layer_density(layers[i], factor)
            
            elif i > layer_idx:
                if self.use_average_density:
                    upper_density = self.calculate_average_density(layers[i])
                    lower_density = self.calculate_average_density(layers[i-1])
                else:
                    upper_density = self.density(boundaries[i-1] + 1e-6)
                    lower_density = self.density(boundaries[i-1] - 1e-6)
                
                if upper_density > lower_density:
                    factor = upper_density / lower_density
                    if self.verbose:
                        print(f"Layer {layers[i]}: Upper density = {upper_density:.4f}, Lower density = {lower_density:.4f}, Factor = {factor:.4f}")
                    self.change_layer_density(layers[i], 1/factor)

        if self.verbose:
            print("Hydrostatic equilibrium imposed.")

    def check_hydrostatic_constraint(self, change_type):    
        boundaries = list(self.modified_boundaries.values())
        layers = ['IC', 'OC', 'IM', 'OM']

        if change_type not in layers:
            return True

        layer_idx = layers.index(change_type)
        for i in range(len(layers)):
            if i < layer_idx:
                if self.use_average_density:
                    upper_density = self.calculate_average_density(layers[i+1])
                    lower_density = self.calculate_average_density(layers[i])
                else:
                    upper_density = self.density(boundaries[i] + 1e-6, use_modified=True)
                    lower_density = self.density(boundaries[i] - 1e-6, use_modified=True)
                if upper_density > lower_density:
                    return False
            elif i > layer_idx:
                if self.use_average_density:
                    upper_density = self.calculate_average_density(layers[i])
                    lower_density = self.calculate_average_density(layers[i-1])
                else:
                    upper_density = self.density(boundaries[i-1] + 1e-6, use_modified=True)
                    lower_density = self.density(boundaries[i-1] - 1e-6, use_modified=True)
                if lower_density < upper_density:
                    return False
        return True

    def density_IC(self, r):
        return 13.0885 - 8.8381 * (r/self.R)**2

    def density_OC(self, r):
        return 12.5815 - 1.2638*(r/self.R) - 3.6426*(r/self.R)**2 - 5.5281*(r/self.R)**3

    def density_IM(self, r):
        return 7.9565 - 6.4761*(r/self.R) + 5.5283*(r/self.R)**2 - 3.0807*(r/self.R)**3

    def density_OM(self, r):
        if r < 5771:
            return 5.3197 - 1.4836*(r/self.R)
        elif 5771 <= r < 5971:
            return 11.2494 - 8.0298*(r/self.R)
        elif 5971 <= r < 6151:
            return 7.1089 - 3.8045*(r/self.R)
        elif 6151 <= r < 6346.6:
            return 2.6910 + 0.6924*(r/self.R)
        elif 6346.6 <= r < 6356.0:
            return 2.900
        elif 6356.0 <= r < 6368.0:
            return 2.600
        elif 6368.0 <= r <= 6371.0:
            return 1.020
        else:
            return 0  # or raise an error

    def density(self, r, use_modified=True):
        densities = self.modified_densities if use_modified else self.original_densities
        boundaries = self.modified_boundaries if use_modified else self.original_boundaries
        
        if 0 <= r < boundaries['IC_OC']:
            return densities['IC'](r)
        elif boundaries['IC_OC'] <= r < boundaries['OC_IM']:
            return densities['OC'](r)
        elif boundaries['OC_IM'] <= r < boundaries['IM_OM']:
            return densities['IM'](r)
        elif boundaries['IM_OM'] <= r <= self.R:
            return densities['OM'](r)
        else:
            return 0  # or raise an error

    def change_layer_density(self, layer, factor):
        if layer not in self.modified_densities:
            raise ValueError(f"Invalid layer name: {layer}")
        
        original_func = self.original_densities[layer]
        self.layer_factors[layer] = factor
        self.modified_densities[layer] = lambda r: original_func(r) * factor

        # If in single core mode and changing IC or OC, apply the same factor to both
        if self.core_mode == 'single' and layer in ['IC', 'OC']:
            other_layer = 'OC' if layer == 'IC' else 'IC'
            self.layer_factors[other_layer] = factor
            self.modified_densities[other_layer] = lambda r: self.original_densities[other_layer](r) * factor

    def change_boundary(self, boundary, new_value):
        if boundary not in ['IC_OC', 'OC_IM', 'IM_OM']:
            raise ValueError("Can only modify IC_OC, OC_IM, or IM_OM boundaries")
        
        if boundary == 'IC_OC':
            if not 0 < new_value < self.modified_boundaries['OC_IM']:
                raise ValueError("IC_OC boundary must be between 0 and OC_IM boundary")
        elif boundary == 'OC_IM':
            if not self.modified_boundaries['IC_OC'] < new_value < self.modified_boundaries['IM_OM']:
                raise ValueError("OC_IM boundary must be between IC_OC and IM_OM boundaries")
        elif boundary == 'IM_OM':
            if not self.modified_boundaries['OC_IM'] < new_value < self.R:
                raise ValueError("IM_OM boundary must be between OC_IM boundary and Earth's radius")
        
        self.modified_boundaries[boundary] = new_value

    def impose_hydrostatic_equilibrium(self, change_type):
        boundaries = list(self.modified_boundaries.values())
        layers = ['IC', 'OC', 'IM', 'OM']
        layer_idx = layers.index(change_type)

        for i in range(len(layers)):
            if i < layer_idx:
                upper_density = self.density(boundaries[i] + 1e-6)
                lower_density = self.density(boundaries[i] - 1e-6)
                
                if upper_density > lower_density:
                    factor = upper_density / lower_density
                    self.change_layer_density(layers[i], factor)
            
            elif i > layer_idx:
                upper_density = self.density(boundaries[i-1] + 1e-6)
                lower_density = self.density(boundaries[i-1] - 1e-6)
                
                if upper_density > lower_density:
                    factor = upper_density / lower_density
                    if self.verbose:
                        print(f"Layer {layers[i]}: Upper density = {upper_density:.4f}, Lower density = {lower_density:.4f}, Factor = {factor:.4f}")
                    self.change_layer_density(layers[i], 1/factor)

        if self.verbose:
            print("Hydrostatic equilibrium imposed.")
    
    def grid_search_optimize(self, change_type, change_value, grid_size=0.00625, tolerance=0.01):
        original_mass = self.total_mass(use_modified=False)
        original_moi = self.moment_of_inertia(use_modified=False)

        if change_type in ['IC_OC', 'OC_IM', 'IM_OM']:
            self.change_boundary(change_type, change_value)
            # otherwise the grid is too difficult...
            self.core_mode = 'single'
            grid_size=0.05
            variable_layers = ['OC', 'IM', 'OM']
        elif change_type in ['IC', 'OC', 'IM', 'OM']:
            self.change_layer_density(change_type, change_value)
            options = ['IC', 'OC', 'IM', 'OM'] if self.core_mode != 'single' else ['OC', 'IM', 'OM']
            variable_layers = [layer for layer in options if layer != change_type]
            if self.core_mode == 'single' and change_type in ['IC', 'OC']:
                variable_layers = [layer for layer in variable_layers if layer not in ['IC', 'OC']]
            self.impose_hydrostatic_equilibrium(change_type)
        else:
            raise ValueError("Invalid change_type. Must be 'IC_OC', 'OC_IM', 'IM_OM', 'IC', 'OC', 'IM', or 'OM'.")

        print(variable_layers)

        layer_ranges = {
            'IC': (0.6, 1.2),
            'OC': (0.6, 1.2),
            'IM': (0.8, 1.5),
            'OM': (0.3, 1.2)
        }

        grid_values = [np.arange(layer_ranges[layer][0], layer_ranges[layer][1], grid_size) 
               for layer in variable_layers]

        grid = list(product(*grid_values))

        best_result = None
        best_objective = float('inf')

        for point in grid:
            for layer, scale in zip(variable_layers, point):
                self.change_layer_density(layer, scale)
            
            # Check hydrostatic equilibrium constraint
            if not self.check_hydrostatic_constraint(change_type):
                continue

            mass_diff = abs((self.total_mass(use_modified=True) - original_mass) / original_mass)
            moi_diff = abs((self.moment_of_inertia(use_modified=True) - original_moi) / original_moi)
            
            objective_value = max(mass_diff, moi_diff)

            if objective_value < best_objective:
                best_objective = objective_value
                best_result = point
                self.best_mass_diff = mass_diff
                self.best_moi_diff = moi_diff

        if best_result is not None and best_objective <= tolerance:
            if self.verbose:
                print("Grid search optimization successful!")
            for layer, scale in zip(variable_layers, best_result):
                if self.verbose:
                    print(f"{layer} density scaled by: {scale:.4f}")
                self.change_layer_density(layer, scale)
        else:
            if self.verbose:
                print("Grid search optimization failed to find a solution within the specified tolerance.")

        print(best_objective, best_result)
        print(self.best_mass_diff, self.best_moi_diff)
        final_mass_diff = abs((self.total_mass(use_modified=True) - original_mass) / original_mass)
        final_moi_diff = abs((self.moment_of_inertia(use_modified=True) - original_moi) / original_moi)
        if self.verbose:
            print(f"Final relative mass difference: {final_mass_diff:.6f}")
            print(f"Final relative MoI difference: {final_moi_diff:.6f}")

        return best_result is not None and best_objective <= tolerance

    def plot_results(self, show_viability=True, conditions=None, figname=None, viability_c='pink'):
        fig, ax = plt.subplots(figsize=(4, 3.5))
        
        r = np.linspace(0, self.R, 1000)
        original_density = [self.density(ri, use_modified=False) for ri in r]
        
        ax.plot(r, original_density, label='PREM', color='orangered', linewidth=1)
        
        viability_regions = {
            'IC': (0.7, 1.1),
            'OC': (0.7, 1.1),
            'IM': (0.85, 1.45),
            'OM': (0.35, 1.1)
        }

        # viability_regions = {
        #     'IC': (0.9, 7),
        #     'OC': (0.59, 1.04),
        #     'IM': (1-0.02*factor, 1+0.05*factor),
        #     'OM': (1-0.02*factor, 1+0.05*factor)
        # }
        
        layer_boundaries = {
            'IC': (0, self.original_boundaries['IC_OC']),
            'OC': (self.original_boundaries['IC_OC'], self.original_boundaries['OC_IM']),
            'IM': (self.original_boundaries['OC_IM'], self.original_boundaries['IM_OM']),
            'OM': (self.original_boundaries['IM_OM'], self.R)
        }
        
        # Plot viability region if enabled
        if show_viability:
            for layer in ['IC', 'OC', 'IM', 'OM']:
                layer_start, layer_end = layer_boundaries[layer]
                layer_mask = (r >= layer_start) & (r <= layer_end)
                r_layer = r[layer_mask]
                
                original_density_layer = [self.density(ri, use_modified=False) for ri in r_layer]
                
                lower_viability = [d * viability_regions[layer][0] for d in original_density_layer]
                upper_viability = [d * viability_regions[layer][1] for d in original_density_layer]
                
                lbl = 'Viability' if layer == 'IC' else None
                ax.fill_between(r_layer, lower_viability, upper_viability, alpha=1.0, color=viability_c, label=lbl)
        
        # Plot custom shades for each condition
        if conditions:
            for condition in conditions:
                label = condition['label']
                color = condition['color']
                
                for layer in ['IC', 'OC', 'IM', 'OM']:
                    layer_start, layer_end = layer_boundaries[layer]
                    layer_mask = (r >= layer_start) & (r <= layer_end)
                    r_layer = r[layer_mask]
                    
                    original_density_layer = [self.density(ri, use_modified=False) for ri in r_layer]
                    
                    lower_factor, upper_factor = condition[layer]
                    lower_factor = max(lower_factor, viability_regions[layer][0])
                    upper_factor = min(upper_factor, viability_regions[layer][1])
                    
                    lower_shade = [d * lower_factor for d in original_density_layer]
                    upper_shade = [d * upper_factor for d in original_density_layer]
                    
                    lbl = label if layer == 'IC' else None
                    ax.fill_between(r_layer, lower_shade, upper_shade, alpha=1., color=color, label=lbl)
        
        ax.set_ylim(1.,16)
        ax.set_xlim(0.,self.R*1.01)

        # # Plot original boundaries and add text annotations
        # for boundary, value in self.original_boundaries.items():
        #     ax.axvline(x=value, color='gray', linestyle=':', linewidth=2)
        #     ax.text(value+350, ax.get_ylim()[1]*0.75, boundary.replace('_', '-'), rotation=90, 
        #              va='bottom', ha='right', fontsize=12, color='gray')
        
        ax.set_xlabel('Radius (km)')
        ax.set_ylabel('Density (g/cm³)')
        
        # # Adjust subplot parameters to make room for the legend
        # plt.subplots_adjust(top=0.9)
        
        # Place legend above the plot
        ax.legend(loc='upper center', bbox_to_anchor=(0.47, 1.32), ncol=2, fontsize=12, frameon=False, columnspacing=0.8)
        
        plt.tight_layout()

        if figname is not None:
            plt.savefig(figname, bbox_inches='tight')
        
        plt.show()

    def plot_improved_density_profiles(self, figname, h_factor=0.7, use_legend=True):
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        fig, ax = plt.subplots(figsize=(4, 3))
        
        r = np.linspace(0, self.R, 1000)
        original_density = [self.density(ri, use_modified=False) for ri in r]
        modified_density = [self.density(ri, use_modified=True) for ri in r]
        
        ax.plot(r, original_density, label='PREM', color='navy', linewidth=2)
        ax.plot(r, modified_density, label='Modified', color='orangered', linestyle='--', linewidth=2)

        # Plot boundaries and add text annotations
        for boundary, value in self.original_boundaries.items():
            if self.core_mode == 'single':
                if boundary == 'IC_OC':
                    continue

            print(boundary)
            ax.axvline(x=value, color='gray', linestyle=':', linewidth=2)
            ax.text(value+100, ax.get_ylim()[1]*h_factor, boundary.replace('_', '-').replace('OC', 'core').replace('IM', 'lman').replace('OM','uman'), rotation=90, 
                    va='top', ha='left', fontsize=12, color='gray')

        if self.change_type == 'OC_IM':
            modified_value = self.modified_boundaries[self.change_type]
            ax.axvline(x=modified_value, color='cornflowerblue', linestyle=':', linewidth=2)
            ax.text(modified_value-100, ax.get_ylim()[1]*0.5, 'New ' + self.change_type.replace('_', '-').replace('OC', 'C'), rotation=90, 
                    va='top', ha='right', fontsize=12, color='cornflowerblue')
        
        ax.set_xlabel('Radius (km)')
        ax.set_ylabel('Density (g/cm³)')
        
        # Adjust y-axis limits
        ax.set_ylim(1.5)
        ax.set_xlim(0, 6800)
        
        # Add legend
        if use_legend:
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, fontsize=11, frameon=False)
            
        # Adjust subplot parameters to make room for the legend
        plt.subplots_adjust(top=0.85)
        
        plt.tight_layout()
        if figname is not None:
            plt.savefig(figname, bbox_inches='tight')
        
        plt.show()

    def mass_integrand(self, r, use_modified=True):
        r_m = r * 1000
        density_kg_m3 = self.density(r, use_modified) * 1000
        return 4 * np.pi * r_m**2 * density_kg_m3

    def moi_integrand(self, r, use_modified=True):
        r_m = r * 1000
        density_kg_m3 = self.density(r, use_modified) * 1000
        return (8/3) * np.pi * r_m**4 * density_kg_m3

    def total_mass(self, use_modified=True):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=IntegrationWarning)
            return quad(lambda r: self.mass_integrand(r, use_modified), 0, self.R)[0]

    def moment_of_inertia(self, use_modified=True):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=IntegrationWarning)
            return quad(lambda r: self.moi_integrand(r, use_modified), 0, self.R)[0]

    def optimize(self):
        print('change_type: ', self.change_type)
        if not self.change_type:
            self.evaluate_modified_density(False)
            return True

        success = self.grid_search_optimize(self.change_type, self.change_value, grid_size=0.00625)
        print('Success: ', success)
        self.evaluate_modified_density()
        return success

    def load_data(self):
        self.r = np.loadtxt(self.PREM_file)[:, 0]
        self.rho = np.loadtxt(self.PREM_file)[:, 1]
        assert self.r is not None, 'Radius data not loaded'
        assert self.rho is not None, 'Density data not loaded'
        self.raw_PREM = np.column_stack((self.r, self.rho))

    def evaluate_modified_density(self, use_modified=True):
        self.mod_PREM = np.column_stack((self.r, [self.density(ri, use_modified=use_modified) for ri in self.r]))

    def dump_prem(self, filename='mod_PREM.dat'):
        """Save the modified PREM data to a file."""
        self.evaluate_modified_density()
        output_path = os.path.join(self.outdir, filename)
        np.savetxt(output_path, self.mod_PREM, fmt='%.3f', delimiter='\t')
        if self.verbose:
            print(f"Modified PREM data saved to {output_path}")

    def plot_earth_circle(self, figname=None, max_size=4000):
        fig, ax = plt.subplots(figsize=(4, 4))
        
        colors = ['cornflowerblue', 'lavender']
        n_bins = 100
        cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)
        factor = 0.4
        
        max_density = max(self.density(r) for r in np.linspace(0, self.R, 1000))
        for r in np.linspace(0, self.R, 1000):
            density = self.density(r)
            color = cmap(density / max_density)
            circle = plt.Circle((0.5, 0.5), r/self.R*factor, color=color, fill=False, linewidth=2)
            ax.add_artist(circle)

        layers = ['OC_IM', 'IM_OM']
        for layer in layers:
            r = self.original_boundaries[layer]
            circle = plt.Circle((0.5, 0.5), r/self.R*factor, color='black', fill=False, linestyle='-', linewidth=1)
            ax.add_artist(circle)
        circle = plt.Circle((0.5, 0.5), factor*1.01, color='black', fill=False, linestyle='-', linewidth=1)
        ax.add_artist(circle)

        for txt,pos in zip(['Core', 'IM', 'OM'], [0.5, 0.785, 0.885]):
            ax.text(0.5, 1-pos, txt, color='black', fontweight='normal', ha='center', va='center', fontsize=12, font='Serif')

        # layer_names = ['IC', 'OC', 'IM', 'OM']
        # for i, (layer, name) in enumerate(zip(layers, layer_names)):
        #     r = self.original_boundaries[layer]
        #     print(name)
        #     #print(r/self.R*factor)
        #     circle = plt.Circle((0.5, 0.5), r/self.R*factor, color='black', fill=False, linestyle='-', linewidth=1)
        #     ax.add_artist(circle)
        #     ax.text(0.435+r/self.R*0.38, 0.5, name, color='black', fontweight='normal', ha='center', va='center', fontsize=14, font='Serif')

        # circle = plt.Circle((0.5, 0.5), factor, color='black', fill=False, linestyle='-', linewidth=1)
        # ax.add_artist(circle)
        
        # Set equal aspect ratio and remove axes
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Explicitly set the limits to show the whole circle and detector
        ax.set_xlim(0.5-factor-0.05, 0.5+factor+0.05)
        ax.set_ylim(0.5-factor-0.05, 0.9+0.05)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_density))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, label='Density (g/cm³)', shrink=0.7,  pad=0.0)
        cbar.ax.set_ylabel('Density (g/cm³)', rotation=270, labelpad=15, fontsize=14, font='Serif')
        cbar.ax.tick_params(labelsize=10)

        det_color = 'k'
        
        detector_width = 0.015
        detector_height = 0.015
        detector = Rectangle((0.5-detector_width/2, 0.895), detector_width, detector_height, 
                             fill=True, edgecolor=det_color,  facecolor=det_color, linewidth=2, zorder=99)
        
        ax.add_artist(detector)
        ax.text(0.5, 0.9 + 0.02, 'Detector', 
                ha='center', va='bottom', fontsize=14, color=det_color, font='Serif')
        
        arrow_col='black'
        
        earth_D = 0.785
        angle = 22
        x_shift = 0.5-earth_D*np.sin(angle*np.pi/360)
        y_shift = 0.1+(0.8-earth_D*np.cos(angle*np.pi/360))
        dX = 0.5-x_shift
        dY = 0.895-y_shift

        print('COSTHETA: ', 1-np.cos(np.arctan(dY/dX)))

        ax.plot([x_shift, 0.5], [y_shift, 0.895], color=arrow_col, linestyle='-', linewidth=1.5)
        ax.annotate('', xy=(x_shift+dX*0.1, y_shift+dY*0.1), xytext=(x_shift, y_shift),
                arrowprops=dict(arrowstyle='->,head_width=0.5,head_length=0.7', color=arrow_col, lw=1),
                xycoords='data', textcoords='data')
        ax.text(0.34, 0.07, '$\\nu$', ha='center', va='bottom', fontsize=24, rotation=68, color=arrow_col)
        ax.text(0.335+dX*0.45, 0.07+dY*0.45, 'cos $\\theta$=-0.8', ha='center', va='bottom', fontsize=12, rotation=78, color=arrow_col)
        
        plt.tight_layout()
        
        if figname:
            plt.savefig(figname, bbox_inches='tight')
        else:
            plt.savefig('EarthDiagram.pdf', bbox_inches='tight')