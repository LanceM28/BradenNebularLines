# importing packages
import numpy as np
import shutil
import os
import matplotlib.pyplot as plt
import emission
import astropy
import yt
from yt.units import dimensions
import copy
from scipy.special import voigt_profile
from astropy.cosmology import FlatLambdaCDM
from matplotlib.colors import LogNorm
import sys
from scipy.ndimage import gaussian_filter
from matplotlib.gridspec import GridSpec

'''
galaxy_visualization.py

Author: Braden Nowicki

Visualization and analysis routines for RAMSES-RT Simulations.

'''

'''
Projection, Slice Plot Routines
'''

# Cloudy Grid Run Bounds (log values)
# Umin, Umax, Ustep: -6.0 1.0 0.5
# Nmin, Nmax, Nstep: -1.0 6.0 0.5 
# Tmin, Tmax, Tstop: 3.0 6.0 0.1

'''
lines=["H1_6562.80A","O1_1304.86A","O1_6300.30A","O2_3728.80A","O2_3726.10A",
       "O3_1660.81A","O3_1666.15A","O3_4363.21A","O3_4958.91A","O3_5006.84A", 
       "He2_1640.41A","C2_1335.66A","C3_1906.68A","C3_1908.73A","C4_1549.00A",
       "Mg2_2795.53A","Mg2_2802.71A","Ne3_3868.76A","Ne3_3967.47A","N5_1238.82A",
       "N5_1242.80A","N4_1486.50A","N3_1749.67A","S2_6716.44A","S2_6730.82A"]

wavelengths=[6562.80, 1304.86, 6300.30, 3728.80, 3726.10, 
             1660.81, 1666.15, 4363.21, 4958.91, 5006.84, 
             1640.41, 1335.66, 1906.68, 1908.73, 1549.00, 
             2795.53, 2802.71, 3868.76, 3967.47, 1238.82, 
             1242.80, 1486.50, 1749.67, 6716.44, 6730.82]
'''

class VisualizationManager:

    def __init__(self, filename, lines, wavelengths):
        '''
        
        Parameters:
        filename (str): filepath to the RAMSES-RT output_*/info_*.txt file
        lines (List, strings): List of nebular emission lines
        wavelengths (List, floats): List of corresponding wavelengths

        file_dir (str): filepath to output directory
        output_file (str): output folder, e.g. output_00273
        sim_run (str): Time slice number for simulation eg. 00273
        info_file (str): Filename with info file appended '/info_00273.txt'
        directory (str): analysis output directory
        '''

        
        self.filename = filename
        self.file_dir = os.path.dirname(self.filename)
        self.lines = lines
        self.wavelengths = wavelengths
        self.output_file = self.file_dir.split('/')[-1]
        self.sim_run = self.output_file.split('_')[1]
        #self.info_file = f'{self.file_dir}/info_{self.sim_run}.txt'

        # Analysis directory for saving
        self.directory = f'analysis/{self.output_file}_analysis'

        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        print(f'Filename = {self.filename}')
        print(f'File Directory = {self.file_dir}')
        print(f'Output File = {self.output_file}')
        print(f'Simulation Run = {self.sim_run}')
        print(f'Analysis Directory = {self.directory}')

    
    # Find center of mass of star particles
    # TODO ds and ad as attributes of class?
    def star_center(self, ad):
        '''
        Locate the center of mass of star particles in code units.

        Parameters:
        ad: data object from RAMSES-RT output loaded into yt-project

        Returns:
        ctr_at_code (List, float): Coordinates (code units) of center of mass

        '''

        x_pos = np.array(ad["star", "particle_position_x"])
        y_pos = np.array(ad["star", "particle_position_y"])
        z_pos = np.array(ad["star", "particle_position_z"])
        x_center = np.mean(x_pos)
        y_center = np.mean(y_pos)
        z_center = np.mean(z_pos)
        x_pos = x_pos - x_center
        y_pos = y_pos - y_center
        z_pos = z_pos - z_center
        ctr_at_code = np.array([x_center, y_center, z_center])
        return ctr_at_code


    def proj_plot(self, ds, sp, width, center, field, weight_field):
        '''
        Projection Plot Driver.

        Parameters:
        ds: loaded RAMSES-RT data set
        sp: sphere data object to project within
        width (tuple, int and str): width in code units or formatted with 
            units, e.g. (1500, 'pc')
        center (List, float): center (array of 3 values) in code units
        field (tuple, str): field to project, e.g. ('gas', 'temperature')
        weight_field (tuple, str): field to weight project (or None if 
            unweighted)

        Returns:
        Projection Plot Object
        '''

        if weight_field == None:
            p = yt.ProjectionPlot(ds, "z", field,
                          width=width,
                          data_source=sp,
                          buff_size=(1000, 1000),
                          center=center)
        else:
            p = yt.ProjectionPlot(ds, "z", field,
                          width=width,
                          weight_field=weight_field,
                          data_source=sp,
                          buff_size=(1000, 1000),
                          center=center)
        return p


    def slc_plot(self, ds, width, center, field):
        '''
        Slice Plot Driver.

        Parameters:
        ds: load RAMSES-RT data set
        width (tuple, int and str): width in code units or formatted with 
            units, e.g. (1500, 'pc')
        center (List, float): center (array of 3 values) in code units
        field (tuple, str): field to project, e.g. ('gas', 'temperature')
        
        Returns:
        Slice Plot Object
        '''

        slc = yt.SlicePlot(
                        ds, "z", field,
                        center=center,
                        width=width,
                        buff_size=(1000, 1000))

        return slc
    

    # TODO change redshift to self
    def convert_to_plt(self, yt_plot, plot_type, field, width, 
                       redshift, title, lims=None):
        '''
        Convert a yt projection or slice plot to matplotlib.

        Parameters:
        yt_plot: Projection or Slice Plot Object
        plot_type (str): Type of plot (for filename) - 'proj' or 'slc'
        field (tuple, str): field to plot, e.g. ('gas', 'temperature')
        width (tuple, int and str): width in code units or formatted with 
            units, e.g. (1500, 'pc')
        redshift (float): redshift of current time slice
        title (str): Plot title
        lims (None or List): [vmin, vmax] fixed limits on colorbar values
            for image if desired; otherwise None

        Returns:
        p_img (ndarray, float): 2D numpy array containing the image data

        Saves desired figures with usable file naming scheme.
        '''

        lbox = width[0]
        length_unit = width[1]
        field_comma = field[1].replace('.', ',')

        plot_title = f'{self.output_file}_{lbox}{length_unit}_' + \
            f'{field_comma}_{plot_type}'

        fname = os.path.join(self.directory, plot_title)
        if lims != None:
            fname = fname + '_lims'

        plot_frb = yt_plot.frb
        # TODO check below
        #p_img = np.array(plot_frb['gas', field])
        p_img = np.array(plot_frb[field[0], field[1]])

        # Clip non-positive values to avoid log of zero or negative numbers
        if np.min(p_img) <= 0:
            print('Warning: Data contains non-positive values. Adjusting ' +
                  'for LogNorm.')
            
            # Clip values below 1e-10
            p_img = np.clip(p_img, a_min=1e-10, a_max=None)

        # Replace NaN with 0 and Inf with finite numbers
        if np.any(np.isnan(p_img)) or np.any(np.isinf(p_img)):
            print('Warning: Data contains NaN or Inf values. ' +
                  'Replacing with 0.')
            p_img = np.nan_to_num(p_img)

        # TODO
        #p_img = gaussian_filter(p_img, sigma=1)

        # Set the extent of the plot
        extent_dens = [-lbox / 2, lbox / 2, -lbox / 2, lbox / 2]

        # Define the color normalization based on the range of the data
        if lims == None:
            dens_norm = LogNorm(vmin=np.min(p_img), vmax=np.max(p_img))
        else:
            # Set fixed color normalization limits
            dens_norm = LogNorm(vmin=lims[0], vmax=lims[1])


        # Viridis, Inferno, Magma maps work - perceptually uniform
        # TODO figsize, dpi
        fig = plt.figure(figsize=(8, 6))
        #im = plt.imshow(p_img, norm=dens_norm, extent=extent_dens, 
        #                origin='lower', aspect='auto', 
        #                interpolation='bilinear', cmap='viridis')

        im = plt.imshow(p_img, norm=dens_norm, extent=extent_dens, 
                        origin='lower', aspect='equal', 
                        interpolation='nearest', cmap='viridis')

        plt.xlabel(f'X [{length_unit}]', fontsize=12)
        plt.ylabel(f'Y [{length_unit}]', fontsize=12)
        plt.title(title, fontsize=14)

        plt.xlim(-lbox / 2, lbox / 2)
        plt.ylim(-lbox / 2, lbox / 2)

        cbar = plt.colorbar(im)

        # Add redshift
        plt.text(0.05, 0.05, f'z = {redshift:.5f}', color='white', fontsize=9,
                 ha='left', va='bottom', transform=plt.gca().transAxes)
        # TODO font

        # TODO print

        # Save the figure
        plt.savefig(fname, dpi=300)
        plt.close()

        return p_img


    def plot_wrapper(self, ds, sp, width, center, field_list,
                     weight_field_list, title_list, proj=True, slc=True,
                     lims_dict=None, lims_titles=None):
        '''
        Wrapper for plotting a variety of fields simultaneously.

        Parameters:
        -----------
        ds: loaded RAMSES-RT data set
        sp: sphere data object to project within
        center (List, float): center (array of 3 values) in code units
        width (tuple, int and str): width in code units or formatted with 
            units, e.g. (1500, 'pc')
        field_list (List of tuple, str): list of fields to plot, e.g. 
            ('gas', 'temperature')
        weight_field_list (List of tuple, str): list of fields to weight 
            projections (or None if unweighted)
        title_list (List of str): list of titles associated with plots
        lims_dict (None or Dict): dictionary of [vmin, vmax] fixed limits on
            colorbar values for image if desired; otherwise None
        lims_titles (List, str): titles associated with lims_dict for
            corresponding fields

        Returns:
        --------
        p_img_arr (list, ndarray, float): list of 2D image arrays
        '''

        redshift = ds.current_redshift

        p_img_arr = []

        for i, field in enumerate(field_list):
            if proj:
                p = self.proj_plot(ds, sp, width, center, field, 
                                   weight_field_list[i])
                
                if lims_dict == None:
                    p_img = self.convert_to_plt(p, 'proj', field, width,
                                                redshift,
                                                'Projected ' + title_list[i])
                else:
                    p_img = self.convert_to_plt(p, 'proj', field, width,
                                                redshift,
                                                'Projected ' + title_list[i],
                                                lims_dict[lims_titles[i]])

            if slc:
                p = self.slc_plot(ds, width, center, field)
                
                if lims_dict == None:
                    p_img = self.convert_to_plt(p, 'slc', field, width,
                                                redshift,
                                                title_list[i])
                else:
                    p_img = self.convert_to_plt(p, 'slc', field, width,
                                                redshift, title_list[i],
                                                lims_dict[lims_titles[i]])
                    
            p_img_arr.append(p_img)
        
        return p_img_arr
                    

    def phase_plot(self, ds, sp, x_field, y_field, z_field, extrema,
                   x_label, y_label, z_label):
        '''
        Generate a phase plot.
        
        Parameters:
        -----------
        ds: loaded RAMSES-RT data set
        sp: sphere data object to project within
        x_field (tuple, str): field to plot on the x-axis, i.e.
            ('gas', 'my_H_nuclei_density')
        y_field (tuple, str): field to plot on the y-axis, i.e.
            ('gas', 'my_temperature')
        z_field (tuple, str): field to plot with colormap, i.e.
            ('gas', 'flux_H1_6562.80A')
        extrema (dict): Dictionary specifying the extrema of the plot, i.e.
            extrema = {('gas', 'my_H_nuclei_density'): (1e-4, 1e4), 
                            ('gas', 'my_temperature'): (1e3, 1e8)}
        x_label (str): label for x-axis
        y_label (str): label for y-axis
        z_label (str): label for colorbar

        Returns:
        --------
        phase_profile: profile associated with PhasePlot object
            Can extract attributes
        x_vals (ndarray, float): x values associated with phase plot
        y_vals (ndarray, float): y values associated with phase plot
        z_vals (ndarray, float): 2D z values associated with phase plot
        '''

        # TODO add z extrema

        plot_title = f'{self.output_file}_' + \
            f'{x_field[1]}_{y_field[1]}_{z_field[1]}_phase.png'

        fname = os.path.join(self.directory, plot_title)

        profile = yt.create_profile(
            sp,
            [x_field, y_field],
            #n_bins=[128, 128],
            fields=[z_field],
            weight_field=None,
            #units=units,
            extrema=extrema,
        )

        plot = yt.PhasePlot.from_profile(profile)

        phase_profile = plot.profile

        plot.set_colorbar_label(z_field, z_label)
        plot.render()

        x_vals = phase_profile.x
        y_vals = phase_profile.y
        z_vals = phase_profile[z_field]  # alternatively field_data attr
        #print(x_vals.shape)
        #print(y_vals.shape)
        #print(z_vals.shape)

        #p_img = np.reshape(z_vals, (len(x_vals)-1, len(y_vals)-1))
        # TODO may need to average in each bin

        # Get a reference to the matplotlib axes object for the plot
        ax = plot.plots[z_field[0], z_field[1]].axes
        fig = plot.plots[z_field[0], z_field[1]].figure
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        plot.save(fname)

        return (phase_profile, x_vals, y_vals, z_vals)
    
    '''
    def phase_with_profiles(self, ds, sp, phase_profile,
                            x_field, y_field, z_field,
                            x_vals, y_vals, z_vals, x_label, y_label,
                            z_label):
        
        Generate a phase plot with additional profile plots.

        Parameters:
        -----------
        ds: loaded RAMSES-RT data set
        sp: sphere data object to project within
        phase_profile: profile associated with PhasePlot object
            Can extract attributes
        x_field (tuple, str): field to plot on the x-axis, i.e.
            ('gas', 'my_H_nuclei_density')
        y_field (tuple, str): field to plot on the y-axis, i.e.
            ('gas', 'my_temperature')
        z_field (tuple, str): field to plot with colormap, i.e.
            ('gas', 'flux_H1_6562.80A')
        x_vals (ndarray, float): x values associated with phase plot
        y_vals (ndarray, float): y values associated with phase plot
        z_vals (ndarray, float): z values associated with phase plot
        p_img (ndarray, float): 2D array of z values associated with
            phase plot
        x_label (str): label for x-axis
        y_label (str): label for y-axis
        z_label (str): label for colorbar

        Returns:
        --------
        TODO
        

        plot_title = f'{self.output_file}_' + \
            f'{x_field[1]}_{y_field[1]}_{z_field[1]}_phase_profile.png'

        fname = os.path.join(self.directory, plot_title)

        #print(x_vals.shape)
        #print(y_vals.shape)
        #print(z_vals.shape)
        x_vals = np.log10(x_vals)
        y_vals = np.log10(y_vals)
        z_vals = np.log10(z_vals).transpose()

        # Find the location of the peak z value (max of z_vals)
        peak_z_idx = np.unravel_index(np.argmax(z_vals), z_vals.shape)[::-1]
        #print(peak_z_idx)
        peak_x = x_vals[peak_z_idx[0]]
        peak_y = y_vals[peak_z_idx[1]]
        peak_z = z_vals[peak_z_idx]

        # Create the figure and gridspec layout
        fig = plt.figure(figsize=(12, 8))
        gs = GridSpec(4, 4, figure=fig)

        # Central phase plot (imshow)
        ax0 = fig.add_subplot(gs[1:4, 0:3])
        cax = ax0.imshow(z_vals, origin="lower", aspect="auto",
                         extent=(min(x_vals), max(x_vals), min(y_vals),
                                 max(y_vals)))
        ax0.set_xlabel(x_label)
        ax0.set_ylabel(y_label)
        #ax0.set_title(title)
        ax0.scatter(peak_x, peak_y, color="red",
                    label=f"Peak ({peak_x:.2f}, {peak_y:.2f}, {peak_z:.2f})")
        ax0.legend(loc="upper right")

        # Profile plot at the top (z vs x)
        # averaging over y bins (along each x slice)
        ax1 = fig.add_subplot(gs[0, 0:3])
        avg_z_vals_x = np.mean(10 ** z_vals, axis=1)
        ax1.plot(x_vals, np.log10(avg_z_vals_x), color="blue")
        ax1.set_ylabel(z_label)

        # Profile plot on the right (z vs y)
        # Averaging over x bins (along each y slice)
        ax2 = fig.add_subplot(gs[1:4, 3])
        avg_z_vals_y = np.mean(10 ** z_vals, axis=0)
        ax2.plot(np.log10(avg_z_vals_y), y_vals, color="blue")
        ax2.set_xlabel(z_label)

        # Adjust layout and position the colorbar
        # Leave space on the right for colorbar
        fig.tight_layout(rect=[0, 0, 0.9, 1])

        # Add colorbar on the right side
        fig.colorbar(cax, ax=ax0, orientation='vertical', label=z_label)

        # Annotate the peak z value location
        
        ax0.annotate(f"Peak z: ({peak_x:.2f}, {peak_y:.2f}, {peak_z:.2f})", 
                     xy=(peak_x, peak_y), xycoords='data', 
                     xytext=(peak_x + 0.1, peak_y + 0.1), textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='red'),
                     fontsize=10, color="black")
        

        plt.savefig(fname, dpi=300)
        plt.close()
    '''

    def phase_with_profiles(self, ds, sp, phase_profile,
                        x_field, y_field, z_field,
                        x_vals, y_vals, z_vals, x_label, y_label,
                        z_label):
        plot_title = f'{self.output_file}_' + \
                    f'{x_field[1]}_{y_field[1]}_{z_field[1]}_phase_profile.png'

        fname = os.path.join(self.directory, plot_title)

        # Logarithmic scaling of the data
        x_vals = np.log10(x_vals)
        y_vals = np.log10(y_vals)
        z_vals = np.log10(z_vals).transpose()

        # Find the location of the peak z value
        peak_z_idx = np.unravel_index(np.argmax(z_vals), z_vals.shape)[::-1]
        peak_x = x_vals[peak_z_idx[0]]
        peak_y = y_vals[peak_z_idx[1]]
        peak_z = z_vals[peak_z_idx[1]][peak_z_idx[0]]

        # Create the figure and gridspec layout
        fig = plt.figure(figsize=(10, 8))
        gs = GridSpec(4, 4, figure=fig)

        # Central phase plot (imshow)
        ax0 = fig.add_subplot(gs[1:4, 0:3])  # Phase plot takes a larger area
        cax = ax0.imshow(z_vals, origin="lower", aspect="auto",
                         extent=(min(x_vals), max(x_vals), min(y_vals),
                                 max(y_vals)))
        ax0.set_xlabel(x_label)
        ax0.set_ylabel(y_label)
        ax0.scatter(peak_x, peak_y, color="red",
                    label=f"Peak ({peak_x:.2f}, {peak_y:.2f}, {peak_z:.2f})")
        ax0.legend(loc="upper right")

        # Profile plot at the top (z vs x)
        ax1 = fig.add_subplot(gs[0, 0:3])  # Profile plot touching the top border of phase plot
        avg_z_vals_x = np.mean(10 ** z_vals, axis=1)[::-1]
        avg_z_vals_x[avg_z_vals_x < 1e-30] = 1e-30
        ax1.plot(x_vals, np.log10(avg_z_vals_x), color="blue")
        # Add a red dot at the peak location on the top profile
        ax1.scatter(peak_x, np.log10(avg_z_vals_x[np.argmax(x_vals == peak_x)]), color='red', s=50)  # Single red dot at peak
        ax1.set_xlim(ax0.get_xlim())  # Ensure x-axis of top profile matches phase plot
        ax1.tick_params(axis='x', which='both', bottom=False, top=False)  # Remove x-axis ticks
        ax1.set_xticklabels([])

        # Profile plot on the right (z vs y)
        ax2 = fig.add_subplot(gs[1:4, 3])  # Profile plot touching the right border of phase plot
        avg_z_vals_y = np.mean(10 ** z_vals, axis=0)[::-1]
        avg_z_vals_y[avg_z_vals_y < 1e-30] = 1e-30
        ax2.plot(np.log10(avg_z_vals_y), y_vals, color="blue")
        # Add a red dot at the peak location on the right profile
        ax2.scatter(np.log10(avg_z_vals_y[np.argmax(y_vals == peak_y)]), peak_y, color='red', s=50)  # Single red dot at peak
        ax2.set_ylim(ax0.get_ylim())  # Ensure y-axis of right profile matches phase plot
        ax2.tick_params(axis='y', which='both', left=False, right=False)  # Remove y-axis ticks
        ax2.set_yticklabels([])

        # Adjust layout and position the colorbar
        fig.tight_layout(rect=[0, 0, 0.85, 1])  # Leave more space for the colorbar

        # Add colorbar on the right side of the profile plot
        #fig.colorbar(cax, ax=ax0, orientation='vertical', label=z_label, pad=0.05)
        cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.7])  # Position the colorbar outside the profile plot area
        fig.colorbar(cax, cax=cbar_ax, orientation='vertical', label=z_label)

        # Save the figure
        plt.savefig(fname, dpi=300)
        plt.close()

    
    def save_array_with_headers(self, filename, array, headers, delimiter=','):
        '''
        Saves a NumPy array to a text file, including column headers.

        Parameters:
        filename (str): The name of the file to save to.
        array (np.ndarray): The NumPy array to save.
        headers (List): A list of strings representing the column headers.
        delimiter (str, optional): The delimiter to use between values. 
            Defaults to ','.
        '''
        with open(filename, 'w') as file:
            file.write(delimiter.join(headers) + '\n')
            np.savetxt(file, array, delimiter=delimiter, fmt='%s')


    def calc_luminosities(self, sp):
        '''
        Agreggate luminosities for each emission line in sphere sp.

        Parameters:
        sp: Data sphere.

        Returns:
        luminosities (List, float): array of luminosities for corresponding
            emission lines.
        '''

        lum_file_path = os.path.join(self.directory, 
                                     f'{self.output_file}_line_luminosity.txt')

        luminosities = []

        for line in self.lines:
            luminosity=sp.quantities.total_quantity(
                ('gas', 'luminosity_' + line)
            )
            luminosities.append(luminosity.value)
            print(f'{line} Luminosity = {luminosity} erg/s')

        # TODO
        #emission_line_str = ", ".join(self.lines)
        #np.savetxt(lum_file_path, luminosities, delimeter=',', 
        #           header=emission_line_str)

        #self.save_array_with_headers(lum_file_path, luminosities, self.lines)
        
        self.luminosities = luminosities

        with open(lum_file_path, 'w') as file:
            for i, line in enumerate(self.lines):
                file.write(f'{line} Luminosity: {self.luminosities[i]}\n')

        return luminosities
    

    def save_sim_info(self, ds):
        '''
        Save simulation parameters/information.

        Parameters:
        ds: RAMSES data loaded into yt.
        '''

        self.current_time = ds.current_time
        self.domain_dimensions = ds.domain_dimensions
        self.domain_left_edge = ds.domain_left_edge
        self.domain_right_edge = ds.domain_right_edge
        self.cosmological_simulation = ds.cosmological_simulation
        self.current_redshift = ds.current_redshift
        self.omega_lambda = ds.omega_lambda
        self.omega_matter = ds.omega_matter
        self.omega_radiation = ds.omega_radiation
        self.hubble_constant = ds.hubble_constant


        file_path = os.path.join(self.directory, 
                                f'{self.output_file}_sim_info.txt')
        
        with open(file_path, 'w') as file:
            file.write(f'current_time: {self.current_time}\n')
            file.write(f'domain_dimensions: {self.domain_dimensions}\n')
            file.write(f'domain_left_edge: {self.domain_left_edge}\n')
            file.write(f'domain_right_edge: {self.domain_right_edge}\n')
            file.write(f'cosmological_simulation: ' +
                       f'{self.cosmological_simulation}\n')
            file.write(f'current_redshift: {self.current_redshift}\n')
            file.write(f'omega_lambda: {self.omega_lambda}\n')
            file.write(f'omega_matter: {self.omega_matter}\n')
            file.write(f'omega_radiation: {self.omega_radiation}\n')
            file.write(f'hubble_constant: {self.hubble_constant}\n')


        # TODO
        '''
        column_headers = ['current_time', 'domain_dimensions',
                          'domain_left_edge', 'domain_right_edge',
                          'cosmological_simulation', 'current_redshift',
                          'omega_lambda', 'omega_matter',
                          'omega_radiation', 'hubble_constant']
        
        info_arr = [self.current_time, self.domain_dimensions,
                    self.domain_left_edge, self.domain_right_edge,
                    self.cosmological_simulation, self.current_redshift,
                    self.omega_lambda, self.omega_matter,
                    self.omega_radiation, self.hubble_constant]

        file_path = os.path.join(self.directory, 
                                f'{self.output_file}_sim_info.txt')

        self.save_array_with_headers(file_path, info_arr, self.lines)
        '''

        # Copy information files from data folder to analysis
        # TODO logSFC
        sim_info_files = [
            os.path.join(self.file_dir, f'header_{self.sim_run}.txt'),
            os.path.join(self.file_dir, 'hydro_file_descriptor.txt'),
            os.path.join(self.file_dir, f'info_{self.sim_run}.txt'),
            os.path.join(self.file_dir, f'info_rt_{self.sim_run}.txt'),
            os.path.join(self.file_dir, 'namelist.txt')
        ]

        for sim_info_file in sim_info_files:
            shutil.copy2(sim_info_file, self.directory)


    def save_sim_field_info(self, ds, sp):
        '''
        Save min, max, mean, and aggregate of each field in fields array.

        Parameters:
        ds: RAMSES data loaded into yt.
        TODO ad
        '''

        fields = [
            ('gas', 'temperature'),
            ('gas', 'density'),
            #('gas', 'number_density'),
            ('gas', 'my_H_nuclei_density'),
            ('gas', 'my_temperature'),
            ('gas', 'ion_param'),
            ('gas', 'metallicity')
        ]

        for line in self.lines:
            fields.append(('gas', 'flux_'  + line))
            fields.append(('gas', 'luminosity_'  + line))

        # Calculate desired quantities for each field
        field_info = []

        for field in fields:
            min = sp.min(field).value
            print(f'{field}_min: {min}')
            max = sp.max(field).value
            print(f'{field}_min: {max}')
            mean = sp.mean(field).value
            print(f'{field}_min: {mean}')
            agg = sp.quantities.total_quantity(field).value
            print(f'{field}_min: {agg}')

            field_info.append((min, max, mean, agg))

        # Save data to a file
        file_path = os.path.join(self.directory, 
                                f'{self.output_file}_field_info.txt')
        
        with open(file_path, 'w') as file:
            for i, field in enumerate(fields):
                file.write(f'{field}_min: {field_info[i][0]}\n')
                file.write(f'{field}_max: {field_info[i][1]}\n')
                file.write(f'{field}_mean: {field_info[i][2]}\n')
                file.write(f'{field}_agg: {field_info[i][3]}\n')

        '''
        Reading the data file example:

        Regex Pattern for float: r'[-+]?\d*\.\d+([eE][-+]?\d+)?'
        Scientific Notation Possible

        import re

        with open('data.txt', 'r') as file:
            file_content = file.read()

        temp_min_pattern = fr'{field}_min: [-+]?\d*\.\d+([eE][-+]?\d+)?'

        temp_min = float(re.search(temp_min_pattern, file_content).group(1)) 
        '''

    def plot_cumulative_field(self, ds, sp, field, title, fname,
                              idx_lims=None):
        '''
        Flatten and order the values in an image of a field 

        Parameters:
        -----------
        ds: loaded RAMSES-RT data set
        sp: sphere data object to project within
        field (tuple, str): list of fields to plot, e.g. 
            ('gas', 'temperature')
        title (str): title
        fname (str): figure name
        idx_lims (tuple, int): range of indices/cells to plot
        '''

        pix = sp[field].value
        pix_sort = np.sort(pix, axis=None)[::-1]
        idxs = np.arange(0, len(pix_sort), 1)
        cum_val = np.cumsum(pix_sort) / np.sum(pix_sort)

        if idx_lims is not None:
            idxs = idxs[idx_lims[0]: idx_lims[1]]
            cum_val = cum_val[idx_lims[0]: idx_lims[1]]
        
        fig = plt.figure(figsize=(8, 6))
        plt.xlabel('Index')
        plt.ylabel('Cumulative Value')
        plt.title(f'{title} Cumulative Sum')
        plt.plot(idxs, cum_val)
        plt.grid(True)
        plt.savefig(
            os.path.join(self.directory,
                         f'output_{self.output_file}_{fname}'), dpi=300
        )
        plt.close()


    def spectra_driver(self, ds, resolving_power, noise_lvl,
                       lum_lims=None, flux_lims=None, linear=False):
        '''
        Generate spectra.

        Parameters:
        -----------
        ds: loaded RAMSES-RT data set
        resolving_power (float): resolving power R = lambda/delta_lambda
            for the observational system. i.e. R = 1000
        noise_lvl (float): noise level/lower floor on signal, i.e. 10e-25
        lum_lims (List, float): manual limits on the luminosity values,
            i.e. lum_lims=[32, 44]
        flux_lims (List, float): manual limits of the flux values
            i.e. flux_lims=[-24, -19]
        '''

        cosmo = FlatLambdaCDM(H0=70, Om0=self.omega_matter)  # around 0.3
        
        # Mpc to cm
        d_1 = cosmo.luminosity_distance(self.current_redshift)*3.086e24
        self.flux_arr = (self.luminosities / (4 * np.pi * d_1 ** 2)).value

        fname = os.path.join(self.directory, self.output_file)

        # Raw spectra values
        self.plot_spectra(noise_lvl, resolving_power, 1000,
                          fname + '_raw_spectra', sim_spectra=False,
                          redshift_wavelengths=False)

        # Sim spectra, not redshifted
        self.plot_spectra(noise_lvl, resolving_power, 1000,
                          fname + '_sim_spectra', sim_spectra=True,
                          redshift_wavelengths=False)

        # Sim spectra, redshifted
        self.plot_spectra(noise_lvl, resolving_power, 1000,
                          fname + '_sim_spectra_redshifted', sim_spectra=True,
                          redshift_wavelengths=True)


        # With limits for animation
        # Sim spectra, not redshifted
        self.plot_spectra(noise_lvl, resolving_power, 1000,
                          fname + '_sim_spectra', sim_spectra=True,
                          redshift_wavelengths=False,
                          lum_lims=lum_lims, flux_lims=flux_lims)

        # Sim spectra, redshifted
        self.plot_spectra(noise_lvl, resolving_power, 1000,
                          fname + '_sim_spectra_redshifted', sim_spectra=True,
                          redshift_wavelengths=True,
                          lum_lims=lum_lims, flux_lims=flux_lims)


    def plot_spectra(self, noise_lvl, resolving_power, pad, figname,
                     sim_spectra=False, redshift_wavelengths=False,
                     lum_lims=None, flux_lims=None, linear=False):
        '''
        Plot a spectrum with certain options.

        Parameters:
        -----------
        noise_lvl (float): noise level/lower floor on signal
        resolving_power (float): resolving power
        pad (float): pad on wavelengths around each voigt profile, i.e. 1000A
        figname (str): filename of figure
        sim_spectra (bool): option to simulate spectra with voigt profiles.
            If False values are plotted in a scatter plot.
        redshift_wavelengths (bool): option to account for redshift in
            wavelengths on the x-axis.
        lum_lims (List, float): manual limits on the luminosity values
        flux_lims (List, float): manual limits of the flux values
        linear (bool): option to plot on linear rather than log y-axis
        '''

        wavelengths = self.wavelengths

        # Display spectra at redshifted wavelengths
        # lambda_obs = (1+z)*lambda_rest
        if redshift_wavelengths:
            wavelengths = (1 + self.redshift) * np.array(wavelengths)
            pad *= 5

        line_widths = np.array(wavelengths) / resolving_power  # Angstroms

        if sim_spectra:
            x_range, y_vals_f = self.plot_voigts(wavelengths, self.flux_arr,
                                                 line_widths,
                                                 [0.0]*len(wavelengths),
                                                 noise_lvl, pad)
            
            fig, ax1 = plt.subplots(1)
            
            if not linear:
                ax1.plot(x_range, np.log10(y_vals_f), color='black')
            else:
                ax1.plot(x_range, y_vals_f, color='black')

            if flux_lims != None:
                if not linear:
                    ax1.set_ylim(flux_lims)
                else:
                    ax1.set_ylim([10**flux_lims[0], 10**flux_lims[1]])

            ax1.set_xlabel(r'Wavelength [$\AA$]')
            if not linear:
                ax1.set_ylabel(
                    r'Log(Flux) [$erg\: s^{-1}\: cm^{-2} \: \AA^{-1}]$'
                )
            else:
                ax1.set_ylabel(r'Flux [$erg\: s^{-1}\: cm^{-2} \: \AA^{-1}]$')

            flux_fname = figname + '_flux'
            plt.savefig(flux_fname)
            plt.close()

            fig, ax1 = plt.subplots(1)
            x_range, y_vals_l = self.plot_voigts(wavelengths, self.luminosities,
                                                 line_widths,
                                                 [0.0]*len(wavelengths),
                                                 noise_lvl, pad)
            
            if not linear:
                ax1.plot(x_range, np.log10(y_vals_l), color='black')
            else:
                ax1.plot(x_range, y_vals_l, color='black')

            if lum_lims != None:
                if linear == False:
                    ax1.set_ylim(lum_lims)
                else:
                    ax1.set_ylim([10**lum_lims[0], 10**lum_lims[1]])

            ax1.set_xlabel(r'Wavelength [$\AA$]')
            if not linear:
                ax1.set_ylabel(
                    r'Log(Luminosity) [$erg\: s^{-1} \: \AA^{-1}]$'
                )
            else:
                ax1.set_ylabel(r'Luminosity [$erg\: s^{-1} \: \AA^{-1}]$')

            lum_fname = figname + '_lum'
            plt.savefig(lum_fname)
            plt.close()

        else:
            fig, (ax1, ax2) = plt.subplots(2, sharex=True)
            ax1.plot(wavelengths, np.log10(self.flux_arr), 'o')
            ax2.plot(wavelengths, np.log10(self.luminosities), 'o')
            ax2.set_xlabel(r'Wavelength [$\AA$]')
            ax1.set_ylabel(r'Log(Flux) [$erg\: s^{-1}\: cm^{-2}$]')
            ax2.set_ylabel(r'Log(Luminosity) [$erg\: s^{-1}$]')
            plt.savefig(figname)
            plt.close()


    def plot_voigts(self, centers, amplitudes, sigmas, gammas,
                    noise_lvl, pad):
        '''
        Plot voigt profiles for spectral lines over a specified noise level.

        Parameters:
        -----------
        All lists must be of same length.

        centers (list, float): centers of voigt profiles
        amplitudes (list, float): corresponding amplitudes (i.e.,
            luminosities) for each profile
        sigmas (list, float): list of associated standard deviations of
            a normal distribution
        gammas (list, float): list of associated FWHM of Cauchy distribution
        noise_lvl (float): noise level/lower floor on signal
        pad (float): pad on wavelengths around each voigt profile, i.e. 1000A

        Returns:
        ----------
        x_range (array, float): array of x values (wavelengths)
        y_vals (array, float): array of accumulated y values (i.e., the
            sum of luminosities or fluxes from each voigt profile)
        '''

        # TODO noiseless profile, Poisson noise

        x_range = np.linspace(min(centers) - pad, max(centers) + pad, 1000)
        y_vals = np.zeros_like(x_range) + noise_lvl

        for amp, center, sigma, gamma in \
            zip(amplitudes, centers, sigmas, gammas):
            y_vals += (amp) * voigt_profile(x_range - center, sigma, gamma)

            #if amp > noise_lvl:
                #y_vals += (amp-noise_lvl)*voigt_profile(x_range - center,
                #   sigma, gamma) # - noise after no sub

        #y_vals += noise_lvl

        return x_range, y_vals




'''
Projection Plots of Ionization Parameter, Number Density, 
Mass Density, Temperature, Metallicity
'''

'''
# TODO update convert to plt to do both versions of plots with and without limits

# Wrapper for making diagnostic plots of given fields
def plot_diagnostics(ds, sp, data_file, center, width, lims_dict=None):
    redshift = ds.current_redshift

    # Ionization Parameter
    proj_ion_param = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'ion-param'), ('gas', 'number_density'))
    if lims_dict == None:
        convert_to_plt_2(data_file, proj_ion_param, 'proj', 'ion-param', width, redshift, 'Ionization Parameter')
    else:
        convert_to_plt_2(data_file, proj_ion_param, 'proj', 'ion-param', width, redshift, 'Ionization Parameter', lims_dict['Ionization Parameter'])

    # Number Density
    proj_num_density = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'number_density'), ('gas', 'number_density'))
    if lims_dict == None:
        convert_to_plt_2(data_file, proj_num_density, 'proj', 'number_density', width, redshift, r'Number Density [$cm^{-3}$]')
    else:
        convert_to_plt_2(data_file, proj_num_density, 'proj', 'number_density', width, redshift, r'Number Density [$cm^{-3}$]', lims_dict['Number Density'])

    # Mass Density
    proj_density = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'density'), ('gas', 'number_density'))
    if lims_dict == None: 
        convert_to_plt_2(data_file, proj_density, 'proj', 'density', width, redshift, r'Density [$g\: cm^{-3}$]')
    else:
        convert_to_plt_2(data_file, proj_density, 'proj', 'density', width, redshift, r'Density [$g\: cm^{-3}$]', lims_dict['Mass Density'])

    # Temperature
    proj_temp = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'temperature'), ('gas', 'number_density'))
    if lims_dict == None:
        convert_to_plt_2(data_file, proj_temp, 'proj', 'temperature', width, redshift, 'Temperature [K]')
    else:
        convert_to_plt_2(data_file, proj_temp, 'proj', 'temperature', width, redshift, 'Temperature [K]', lims_dict['Temperature'])

    # Metallicity
    proj_metallicity = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'metallicity'), ('gas', 'number_density'))
    if lims_dict == None:
        convert_to_plt_2(data_file, proj_metallicity, 'proj', 'metallicity', width, redshift, 'Metallicity')
    else:
        convert_to_plt_2(data_file, proj_metallicity, 'proj', 'metallicity', width, redshift, 'Metallicity', lims_dict['Metallicity'])


Visualizing Line Intensities


def plot_intensities(ds, sp, data_file, center, width, lims_dict=None):
    redshift = ds.current_redshift

    for line in lines:
        proj = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'intensity_' + line), None)

        if line == 'H1_6562.80A':
            line_title = r'H$\alpha$_6562.80A'
        else:
            line_title = line

        if lims_dict == None:
            convert_to_plt_2(data_file, proj, 'proj', 'intensity_' + line, width, redshift, \
                        'Projected ' + line_title.replace('_', ' ') + r' Flux [$erg\: s^{-1}\: cm^{-2}$]')
        else:
            convert_to_plt_2(data_file, proj, 'proj', 'intensity_' + line, width, redshift, \
                        'Projected ' + line_title.replace('_', ' ') + r' Flux [$erg\: s^{-1}\: cm^{-2}$]', lims_dict[line])

    proj_halpha_l = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'luminosity_H1_6562.80A'), None)
    convert_to_plt_2(data_file, proj_halpha_l, 'proj', 'luminosity_H1_6562.80A', width, redshift, r'Projected H$\alpha$ 6562.80A Luminosity [$erg\: s^{-1}$]')

    slc_halpha = slc_plot(ds, (width, 'pc'), center, ('gas', 'intensity_H1_6562.80A'))
    convert_to_plt_2(data_file, slc_halpha, 'slc', 'intensity_H1_6562.80A', width, redshift, r'H$\alpha$ 6562.80A Flux [$erg\: s^{-1}\: cm^{-2}$]')
'''





# TODO
'''
Plot Spectra at simulated wavelengths
'''

'''
def spectra_driver(ds, luminosities, data_file, linear=False):
    z = ds.current_redshift
    omega_matter = ds.omega_matter
    omega_lambda = ds.omega_lambda
    cosmo = FlatLambdaCDM(H0=70, Om0=omega_matter)#, Om0=0.3)
    d_l = cosmo.luminosity_distance(z)*3.086e24 # convert Mpc to cm
    flux_arr = luminosities/(4*np.pi*d_l**2)
    flux_arr = flux_arr.value

    # resolving power
    # R = lambda/delta_lambda
    R = 1000

    directory = 'analysis/' + data_file + '_analysis'

    if not os.path.exists(directory):
        os.makedirs(directory)

    # Raw spectra values
    plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_raw_spectra"), \
             sim_spectra=False, redshift_wavelengths=False)

    # Sim spectra, not redshifted
    plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_sim_spectra"), \
             sim_spectra=True, redshift_wavelengths=False)

    # Sim spectra, redshifted
    plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_sim_spectra_redshifted"), \
             sim_spectra=True, redshift_wavelengths=True)
    
    # With Limits for animation
    # Sim spectra, not redshifted
    plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_sim_spectra_lims"), \
             sim_spectra=True, redshift_wavelengths=False, lum_lims=[32, 44], flux_lims=[-24, -19])

    # Sim spectra, redshifted
    plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_sim_spectra_redshifted_lims"), \
             sim_spectra=True, redshift_wavelengths=True, lum_lims=[32, 44], flux_lims=[-24, -19])
    
    # Linear Scale
    # Sim spectra, not redshifted
    #plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_sim_spectra_lims_linear"), \
    #         sim_spectra=True, redshift_wavelengths=False, lum_lims=[32, 38], flux_lims=[-24, -20], linear=True)

    # Sim spectra, redshifted
    #plot_spectra(wavelengths, luminosities, flux_arr, z, 10e-25, R, fname=os.path.join(directory, data_file + "_sim_spectra_redshifted_lims_linear"), \
    #         sim_spectra=True, redshift_wavelengths=True, lum_lims=[32, 38], flux_lims=[-24, -20], linear=True)

    

def plot_spectra(wavelengths, luminosities, flux_arr, z, noise_lvl, R, \
                 fname, sim_spectra=False, redshift_wavelengths=False, lum_lims=None, flux_lims=None, linear=False):

    pad = 1000

    # Display spectra at redshifted wavelengths
    # lambda_obs = (1+z)*lambda_rest
    if redshift_wavelengths:
        wavelengths = (1+z)*np.array(wavelengths)
        pad = 5000

    line_widths = np.array(wavelengths)/R # Angstrom

    if sim_spectra:
        fig, ax1 = plt.subplots(1)
        x_range, y_vals_f = plot_voigts(wavelengths, flux_arr, line_widths, [0.0]*len(wavelengths), noise_lvl, pad)
        if linear == False:
            ax1.plot(x_range, np.log10(y_vals_f), color='black')
        else:
            ax1.plot(x_range, y_vals_f, color='black')
        if flux_lims != None:
            if linear == False:
                ax1.set_ylim(flux_lims)
            else:
                ax1.set_ylim([10**flux_lims[0], 10**flux_lims[1]])
        ax1.set_xlabel(r'Wavelength [$\AA$]')
        if linear == False:
            ax1.set_ylabel(r'Log(Flux) [$erg\: s^{-1}\: cm^{-2} \: \AA^{-1}]$')
        else:
            ax1.set_ylabel(r'Flux [$erg\: s^{-1}\: cm^{-2} \: \AA^{-1}]$')
        plt.savefig(fname)
        plt.close()

        fig, ax1 = plt.subplots(1)
        x_range, y_vals_l = plot_voigts(wavelengths, luminosities, line_widths, [0.0]*len(wavelengths), noise_lvl, pad)
        if linear == False:
            ax1.plot(x_range, np.log10(y_vals_l), color='black')
        else:
            ax1.plot(x_range, y_vals_l, color='black')
        if lum_lims != None:
            if linear == False:
                ax1.set_ylim(lum_lims)
            else:
                ax1.set_ylim([10**lum_lims[0], 10**lum_lims[1]])
        ax1.set_xlabel(r'Wavelength [$\AA$]')
        if linear == False:
            ax1.set_ylabel(r'Log(Luminosity) [$erg\: s^{-1} \: \AA^{-1}]$')
        else:
            ax1.set_ylabel(r'Luminosity [$erg\: s^{-1} \: \AA^{-1}]$')
        lum_fname = fname + '_lum'
        plt.savefig(lum_fname)
        plt.close()
    else:
        fig, (ax1, ax2) = plt.subplots(2, sharex=True)
        ax1.plot(wavelengths, np.log10(flux_arr), 'o')
        ax2.plot(wavelengths, np.log10(luminosities), 'o')
        ax2.set_xlabel(r'Wavelength [$\AA$]')
        ax1.set_ylabel(r'Log(Flux) [$erg\: s^{-1}\: cm^{-2}$]')
        ax2.set_ylabel(r'Log(Luminosity) [$erg\: s^{-1}$]')
        plt.savefig(fname)
        plt.close()

# Plot voigt profiles for spectral lines over a noise level
# sigma - stdev of normal dist
# gamma - FWHM of cauchy dist
# TODO noiseless profile, Poisson Noise
def plot_voigts(centers, amplitudes, sigmas, gammas, noise_lvl, pad):

    x_range = np.linspace(min(centers) - pad, max(centers) + pad, 1000)
    y_vals = np.zeros_like(x_range)+noise_lvl

    for amp, center, sigma, gamma in zip(amplitudes, centers, sigmas, gammas):
        y_vals += (amp)*voigt_profile(x_range - center, sigma, gamma)
        #if amp > noise_lvl:
            #y_vals += (amp-noise_lvl)*voigt_profile(x_range - center, sigma, gamma) # - noise after no sub

    #y_vals += noise_lvl

    return x_range, y_vals



Star + Gas Plot
Adapted from work by Sarunyapat Phoompuang


def star_gas_overlay(ds, ad, sp, data_file, center, width, field, lims_dict=None):
    redshift = ds.current_redshift

    lims = lims_dict["H1_6562.80A"] # TODO

    directory = 'analysis/' + data_file + '_analysis'

    if not os.path.exists(directory):
        os.makedirs(directory)

    fname = os.path.join(directory, data_file + '_' + str(width) + 'pc_stellar_dist')

    # Finding center of the data
    x_pos = np.array(ad["star", "particle_position_x"])
    y_pos = np.array(ad["star", "particle_position_y"])
    z_pos = np.array(ad["star", "particle_position_z"])
    x_center = np.mean(x_pos)
    y_center = np.mean(y_pos)
    z_center = np.mean(z_pos)
    x_pos = x_pos - x_center
    y_pos = y_pos - y_center
    z_pos = z_pos - z_center

    # Create a ProjectionPlot
    #p = proj_plot(ds, sp, (width, 'pc'), center, ('gas', 'intensity_H1_6562.80A'), None)
    p = yt.ProjectionPlot(ds, "z", field,
                      width=(width, 'pc'),
                      data_source=sp,
                      buff_size=(2000, 2000),
                      center=center)

    p_frb = p.frb  # Fixed-Resolution Buffer from the projection
    p_img = np.array(p_frb["gas", 'intensity_H1_6562.80A'])
    star_bins = 2000
    star_mass = np.ones_like(x_pos) * 10
    pop2_xyz = np.array(ds.arr(np.vstack([x_pos, y_pos, z_pos]), "code_length").to("pc")).T
    extent_dens = [-width/2, width/2, -width/2, width/2]
    #gas_range = (20, 2e4)
    norm1 = LogNorm(vmin=lims[0], vmax=lims[1])
    #norm1 = LogNorm(vmin = gas_range[0], vmax = gas_range[1])
    #plt.figure(figsize = (12,8))
    #plt.imshow(p_img, norm = norm1, extent = extent_dens, origin = 'lower', aspect = 'auto', cmap = 'inferno')
    #plt.colorbar(label = r"Number Density [1/cm$^3$]")
    # plt.scatter(pop2_xyz[:, 0], pop2_xyz[:, 1], s = 5, marker = '.', color = 'red')
    #plt.xlabel("X (pc)")
    #plt.ylabel("Y (pc)")
    #plt.title("Gas Density and Particle Positions")
    #plt.xlim(-width/2, width/2)
    #plt.ylim(-width/2, width/2)
    #plt.show()
	
    stellar_mass_dens, _, _ = np.histogram2d(pop2_xyz[:, 0],
                                             pop2_xyz[:, 1],
                                             bins = star_bins,
                                             weights = star_mass,
                                             range = [[-width / 2, width / 2],
                                                      [-width / 2, width / 2],
                                                     ],
    )
    stellar_mass_dens = stellar_mass_dens.T
    stellar_mass_dens = np.where(stellar_mass_dens <= 1, 0, stellar_mass_dens)
    stellar_range = [1, 1200]
    norm2 = LogNorm(vmin = stellar_range[0], vmax = stellar_range[1])
    plt.figure(figsize = (12, 8))
    lumcmap = "cmr.amethyst"
    plt.imshow(stellar_mass_dens, norm = norm2, extent = extent_dens, origin = 'lower', aspect = 'auto', cmap = 'winter_r')
    plt.colorbar(label = "Stellar Mass Density")
    plt.xlabel("X (pc)")
    plt.ylabel("Y (pc)")
    plt.title("Stellar Mass Density Distribution")
    plt.savefig(fname=fname)
    plt.close()
    	
    #print(np.min(p_img), np.max(p_img))  # Check for min/max values of p_img
    #print(np.min(stellar_mass_dens), np.max(stellar_mass_dens))  # Check for min/max values of stellar_mass_dens

    overlay_fname = fname + '_H1' #+ field
    fig, ax = plt.subplots(figsize = (12, 8))
    alpha_star = stellar_mass_dens
    alpha_star = np.where(stellar_mass_dens <= 1, 0.0, 1)

    print(alpha_star.shape)
    print(p_img.shape)

    img1 = ax.imshow(p_img, norm = norm1, extent = extent_dens, origin = 'lower', aspect = 'auto', cmap = 'inferno', alpha = 1, interpolation='bilinear')
    cbar1 = fig.colorbar(img1, ax = ax, orientation = 'vertical', pad = 0.02)
    cbar1.set_label('Projected ' + r'H$\alpha$ Flux [$erg\: s^{-1}\: cm^{-2}$]') # TODO
    img2 = ax.imshow(stellar_mass_dens, norm = norm2, extent = extent_dens, origin = 'lower', aspect = 'auto', cmap = 'winter_r', interpolation='bilinear')

    # Make sure alpha_star matches the image shape (it must be the same size as the image)
    if img2.get_array().shape != alpha_star.shape:
        print(f"Shape mismatch: Image shape {img2.get_array().shape} vs. alpha_star shape {alpha_star.shape}")
    
    img2.set_alpha(alpha_star)  # Apply alpha mask after plotting

    cbar2 = fig.colorbar(img2, ax = ax, orientation = 'vertical', pad = 0.02)
    cbar2.set_label("Stellar Mass Density")
    # ax.scatter(pop2_xyz[:, 0], pop2_xyz[:, 1], s=5, marker='.', color='black')
    ax.set_xlabel("X (pc)")
    ax.set_ylabel("Y (pc)")
    ax.set_title(r'H$\alpha$ Flux' + ' and Stellar Mass Density Distribution')
    ax.set_xlim(-width / 2, width / 2)
    ax.set_ylim(-width / 2, width / 2)

    plt.text(0.05, 0.05, f'z = {redshift:.5f}', color='white', fontsize=9, ha='left', va='bottom', \
             transform=plt.gca().transAxes)

    plt.savefig(fname=overlay_fname)
    plt.close()

'''

# sp.quantities.center_of_mass(use_gas=False, use_particles=True, particle_type="star")


# TODO 1D Profile Plots, Phase Plots
# Save header information to a text file
# Save lines and wavelengths
# Save mins, maxs, means
# TODO cloudy run
# TODO change z label on plots annotation, sig figs
# TODO docstrings
# FITS images
# TODO star density plot, star + gas
