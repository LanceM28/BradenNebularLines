import os
import sys
import copy

import numpy as np
import yt
from yt.frontends.ramses.field_handlers import RTFieldFileHandler

from emission import EmissionLineInterpolator
import galaxy_visualization

'''
galaxy_emission.py

Author: Braden Nowicki

Braden Nowicki with Dr. Massimo Ricotti
University of Maryland, College Park Astronomy Department

Script to visualize RAMSES-RT Simulations of high-redshift galaxies in a 
variety of metal lines.
Ionization Parameter, Number Density, and Temperature for each pixel are input
into an interpolator for each line; the interpolator is created via the
module 'emission.py'. An EmissionLineInterpolator object is instantiated
given a filepath for a Cloudy-generated line flux list/data table.

'''

# TODO docstrings
# TODO change to use mytemperature
# TODO new hydrofile versions

'''
-------------------------------------------------------------------------------
Setup fields in yt
-------------------------------------------------------------------------------
'''

# Specify output_* folder
#filename = sys.argv[1]
#print(f'Line List Filepath = {filename}')
filename = "/Users/bnowicki/Documents/Research/Ricotti/output_00273/info_00273.txt"

# f1 = "/Users/bnowicki/Documents/Research/Ricotti/output_00273/info_00273.txt"
# Cloudy Grid Run Bounds for this line list (log values)
# Umin, Umax, Ustep: -6.0 1.0 0.5
# Nmin, Nmax, Nstep: -1.0 6.0 0.5
# Tmin, Tmax, Tstop: 3.0 6.0 0.1

<<<<<<< HEAD
lines=["H1_6562.80A","O1_1304.86A","O1_6300.30A","O2_3728.80A","O2_3726.10A",
       "O3_1660.81A","O3_1666.15A","O3_4363.21A","O3_4958.91A","O3_5006.84A", 
       "He2_1640.41A","C2_1335.66A","C3_1906.68A","C3_1908.73A","C4_1549.00A",
       "Mg2_2795.53A","Mg2_2802.71A","Ne3_3868.76A","Ne3_3967.47A",
       "N5_1238.82A",
       "N5_1242.80A","N4_1486.50A","N3_1749.67A","S2_6716.44A","S2_6730.82A"]

wavelengths=[6562.80, 1304.86, 6300.30, 3728.80, 3726.10, 1660.81, 1666.15,
             4363.21, 4958.91, 5006.84, 1640.41, 1335.66,
             1906.68, 1908.73, 1549.00, 2795.53, 2802.71, 3868.76,
=======
# TODO H2
lines=["H1_6562.80A","O1_1304.86A","O1_6300.30A","O2_3728.80A","O2_3726.10A",
       "O3_1660.81A","O3_1666.15A","O3_4363.21A","O3_4958.91A","O3_5006.84A", 
       "He2_1640.41A","C2_1335.66A","C3_1906.68A","C3_1908.73A","C4_1549.00A",
       "Mg2_2795.53A","Mg2_2802.71A","Ne3_3868.76A","Ne3_3967.47A","N5_1238.82A",
       "N5_1242.80A","N4_1486.50A","N3_1749.67A","S2_6716.44A","S2_6730.82A"]

wavelengths=[6562.80, 1304.86, 6300.30, 3728.80, 3726.10, 1660.81, 1666.15, \
             4363.21, 4958.91, 5006.84, 1640.41, 1335.66, \
             1906.68, 1908.73, 1549.00, 2795.53, 2802.71, 3868.76, \
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
             3967.47, 1238.82, 1242.80, 1486.50, 1749.67, 6716.44, 6730.82]


cell_fields = [
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
<<<<<<< HEAD
=======
    # "dark_matter_density",
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
    "xHI",
    "xHII",
    "xHeII",
    "xHeIII",
]

epf = [
    ("particle_family", "b"),
    ("particle_tag", "b"),
    ("particle_birth_epoch", "d"),
    ("particle_metallicity", "d"),
]

# Ionization Parameter Field
# Based on photon densities in bins 2-4
# Don't include bin 1 -> Lyman Werner non-ionizing
def _ion_param(field, data):
    p = RTFieldFileHandler.get_rt_parameters(ds).copy()
    p.update(ds.parameters)

    cgs_c = 2.99792458e10     #light velocity

    # Convert to physical photon number density in cm^-3
    pd_2 = data['ramses-rt','Photon_density_2']*p["unit_pf"]/cgs_c
    pd_3 = data['ramses-rt','Photon_density_3']*p["unit_pf"]/cgs_c
    pd_4 = data['ramses-rt','Photon_density_4']*p["unit_pf"]/cgs_c

    photon = pd_2 + pd_3 + pd_4

    return photon/data['gas', 'number_density']

<<<<<<< HEAD

=======
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
def _my_temperature(field, data):
    #y(i): abundance per hydrogen atom
    XH_RAMSES=0.76 #defined by RAMSES in cooling_module.f90
    YHE_RAMSES=0.24 #defined by RAMSES in cooling_module.f90
    mH_RAMSES=yt.YTArray(1.6600000e-24,"g") #defined by RAMSES in cooling_module.f90
    kB_RAMSES=yt.YTArray(1.3806200e-16,"erg/K") #defined by RAMSES in cooling_module.f90

    dn=data["ramses","Density"].in_cgs()
    pr=data["ramses","Pressure"].in_cgs()
    yHI=data["ramses","xHI"]
    yHII=data["ramses","xHII"]
    yHe = YHE_RAMSES*0.25/XH_RAMSES
    yHeII=data["ramses","xHeII"]*yHe
    yHeIII=data["ramses","xHeIII"]*yHe
    yH2=1.-yHI-yHII
    yel=yHII+yHeII+2*yHeIII
    mu=(yHI+yHII+2.*yH2 + 4.*yHe) / (yHI+yHII+yH2 + yHe + yel)
    return pr/dn * mu * mH_RAMSES / kB_RAMSES


# TODO see if it works in emission.py
# Luminosity field
# Cloudy Intensity obtained assuming height = 1cm
# Return intensity values erg/s/cm**2
# Multiply intensity at each pixel by volume of pixel -> luminosity
def get_luminosity(line):
   def _luminosity(field, data):
      return data['gas', 'flux_' + line]*data['gas', 'volume']
   return copy.deepcopy(_luminosity)


<<<<<<< HEAD
#number density of hydrogen atoms
def _my_H_nuclei_density(field, data):
    dn=data["ramses","Density"].in_cgs()
    XH_RAMSES=0.76 #defined by RAMSES in cooling_module.f90
    YHE_RAMSES=0.24 #defined by RAMSES in cooling_module.f90
    mH_RAMSES=yt.YTArray(1.6600000e-24,"g") #defined by RAMSES in cooling_module.f90

    return dn*XH_RAMSES/mH_RAMSES

'''
def _pressure(field, data):
    # TODO change to ds.add_field?
    #if 'hydro_thermal_pressure' in dir(ds.fields.ramses):
    return data['ramses', 'hydro_thermal_pressure']

def _xHI(field, data):
    return data['ramses', 'hydro_xHI']

def _xHII(field, data):
    return data['ramses', 'hydro_xHII']

def _xHeII(field, data):
    return data['ramses', 'hydro_xHeII']

def _xHeIII(field, data):
    return data['ramses', 'hydro_xHeIII']
'''
    
=======
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
'''
Add derived fields.
'''

# Ionization parameter
yt.add_field(
<<<<<<< HEAD
    ('gas', 'ion_param'),
=======
    ('gas', 'ion-param'),
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
    function=_ion_param,
    sampling_type="cell",
    units="cm**3",
    force_override=True
)

yt.add_field(
    ("gas","my_temperature"),
    function=_my_temperature,
    sampling_type="cell",
<<<<<<< HEAD
    # TODO units
    #units="K",
    #units="K*cm**3/erg",
    units='K*cm*dyn/erg',
    force_override=True
)

yt.add_field(
    ("gas","my_H_nuclei_density"),
    function=_my_H_nuclei_density,
    sampling_type="cell",
    units="1/cm**3",
    force_override=True
)

yt.add_field(
    ("gas","number_density"),
    function=_my_H_nuclei_density,
    sampling_type="cell",
    units="1/cm**3",
    force_override=True
)

'''
yt.add_field(
    ("ramses","Pressure"),
    function=_pressure,
    sampling_type="cell",
    units="1",
    force_override=True
)

yt.add_field(
    ("ramses","xHI"),
    function=_xHI,
    sampling_type="cell",
    units="1",
    force_override=True
)

yt.add_field(
    ("ramses","xHII"),
    function=_xHII,
    sampling_type="cell",
    units="1",
    force_override=True
)

yt.add_field(
    ("ramses","xHeII"),
    function=_xHeII,
    sampling_type="cell",
    units="1",
    force_override=True
)

yt.add_field(
    ("ramses","xHeIII"),
    function=_xHeIII,
    sampling_type="cell",
    units="1",
    force_override=True
)
'''


=======
    units="K",
    force_override=True
)

>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
# Normalize by Density Squared Flag
dens_normalized = True
if dens_normalized: 
    units = '1/cm**6'
else:
    units = '1'

# Instance of EmissionLineInterpolator for line list at filename
line_list = os.path.join(os.getcwd(), 'nebular_lines_v2/linelist.dat')
emission_interpolator = EmissionLineInterpolator(line_list, lines)

# Add flux and luminosity fields for all lines in the list
for i, line in enumerate(lines):
    yt.add_field(
        ('gas', 'flux_' + line),
        function=emission_interpolator.get_line_emission(
            i, dens_normalized=dens_normalized
        ),
        sampling_type='cell',
        units=units,
        force_override=True
    )
    # TODO change get_line_emission to accept line not idx

    yt.add_field(
        ('gas', 'luminosity_' + line),
        function=emission_interpolator.get_luminosity(lines[i]),
        #function=get_luminosity(lines[i]),
        sampling_type='cell',
        units='1/cm**3',
        force_override=True
    )


'''
-------------------------------------------------------------------------------
Load Simulation Data
Run routines on data
-------------------------------------------------------------------------------
'''

ds = yt.load(filename, extra_particle_fields=epf)
ad = ds.all_data()

print(ds.field_list)

viz = galaxy_visualization.VisualizationManager(filename, lines, wavelengths)
star_ctr = viz.star_center(ad)
sp = ds.sphere(star_ctr, (3000, "pc"))
sp_lum = ds.sphere(star_ctr, (10, 'kpc'))
width = (1500, 'pc')

sim_run = filename.split('/')[-1]

field_list = [
<<<<<<< HEAD
    #('gas', 'temperature'),
    ('gas', 'density'),
    ('gas', 'my_H_nuclei_density'),
    ('gas', 'my_temperature'),
    ('gas', 'ion_param'),
=======
    ('gas', 'temperature'),
    ('gas', 'density'),
    ('gas', 'number_density'),
    ('gas', 'my_temperature'),
    ('gas', 'ion-param'),
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
    ('gas', 'metallicity')
]

weight_field_list = [
<<<<<<< HEAD
    #('gas', 'my_H_nuclei_density'),
    ('gas', 'my_H_nuclei_density'),
    ('gas', 'my_H_nuclei_density'),
    ('gas', 'my_H_nuclei_density'),
    ('gas', 'my_H_nuclei_density'),
    ('gas', 'my_H_nuclei_density')
]

title_list = [
    #'Temperature [K]',
    r'Density [g cm$^{-3}$]',
    r'H Nuclei Number Density [cm$^{-3}$]',
=======
    ('gas', 'number_density'),
    ('gas', 'number_density'),
    ('gas', 'number_density'),
    ('gas', 'number_density'),
    ('gas', 'number_density'),
    ('gas', 'number_density')
]

title_list = [
    'Temperature [K]',
    r'Density [$g\: cm^{-3}$]',
    r'Number Density [$cm^{-3}$]',
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
    'My Temperature [K]',
    'Ionization Parameter',
    'Metallicity'
]

for line in lines:
    if line == 'H1_6562.80A':
        line_title = r'H$\alpha$_6562.80A'
    else:
        line_title = line

    field_list.append(('gas', 'flux_'  + line))
    title_list.append(line_title.replace('_', ' ') + 
                      r' Flux [$erg\: s^{-1}\: cm^{-2}$]')
    weight_field_list.append(None)

    field_list.append(('gas', 'luminosity_'  + line))
    title_list.append(line_title.replace('_', ' ') + 
                      r' Luminosity [$erg\: s^{-1}$]')
    weight_field_list.append(None)


viz.save_sim_info(ds)
viz.plot_wrapper(ds, sp, width, star_ctr, field_list,
                     weight_field_list, title_list, proj=True, slc=False)

#viz.calc_luminosities(sp)



'''
Line Luminosities
'''

'''
# Save data to new directory
directory = 'analysis/' + sim_run + '_analysis'

if not os.path.exists(directory):
    os.makedirs(directory)

lum_file_path = os.path.join(directory, sim_run + "_line_luminosity.txt")


def calc_luminosities(sp, file_path):

    luminosities = []

    for line in lines:
        luminosity=sp.quantities.total_quantity(
            ('gas', 'luminosity_' + line)
        )
        luminosities.append(luminosity.value)
        print(f'{lines} Luminosity = {luminosity} erg/s')
    
    np.savetxt(file_path, luminosities, delimeter=',')

    return luminosities
    
#luminosities = calc_luminosities(sp_lum)

directory = 'analysis/' + sim_run + '_analysis'

if not os.path.exists(directory):
    os.makedirs(directory)


# TODO save mins and max's of fields


Create figures


# store in files JSON TODO
# Dictionary of manual limits for each line for animating
# galaxies over time slices (fix colorbar limits for visual consistency)
lims_00273 = {
    'Ionization Parameter': [10e-7, 10e-1],
    'Number Density': [10e-2, 10e5],
    'Mass Density': [10e-26, 10e-19],
    'Temperature': [10e1, 10e6],
    'Metallicity': [10e-2, 10e1],
    "H1_6562.80A": [10e-7, 10e2],
    "O1_1304.86A": [10e-9, 10e1],
    "O1_6300.30A": [10e-8, 10e-3],
    "O2_3728.80A": [10e-7, 10e-2],
    "O2_3726.10A": [10e-7, 10e-2],
    "O3_1660.81A": [10e-8, 10e-4],
    "O3_1666.15A": [10e-8, 10e-3],
    "O3_4363.21A": [10e-9, 10e-4],
    "O3_4958.91A": [10e-8, 10e-3],
    "O3_5006.84A": [10e-8, 10e-3], 
    "He2_1640.41A": [10e-10, 10e-3],
    "C2_1335.66A": [10e-8, 10e2],
    "C3_1906.68A": [10e-7, 10e-2],
    "C3_1908.73A": [10e-7, 10e-2],
    "C4_1549.00A": [10e-16, 10e-9],
    "Mg2_2795.53A": [10e-8, 10e2],
    "Mg2_2802.71A": [10e-8, 10e2],
    "Ne3_3868.76A": [10e-9, 10e-4],
    "Ne3_3967.47A": [10e-9, 10e-5],
    "N5_1238.82A": [10e-14, 10e-4],
    "N5_1242.80A": [10e-14, 10e-4],
    "N4_1486.50A": [10e-10, 10e-4],
    "N3_1749.67A": [10e-8, 10e-4],
    "S2_6716.44A": [10e-8, 10e-2],
    "S2_6730.82A": [10e-8, 10e-3]
}

lims_fiducial_00319 = {
    'Ionization Parameter': [10e-7, 10e0],
    'Number Density': [10e-3, 10e2],
    'Mass Density': [10e-27, 10e-20],
    'Temperature': [10e2, 10e10],
    'Metallicity': [10e-3, 10e0],
    "H1_6562.80A": [10e-9, 10e2],
    "O1_1304.86A": [10e-13, 10e-4],
    "O1_6300.30A": [10e-12, 10e-3],
    "O2_3728.80A": [10e-11, 10e-2],
    "O2_3726.10A": [10e-11, 10e-2],
    "O3_1660.81A": [10e-12, 10e-3],
    "O3_1666.15A": [10e-11, 10e-3],
    "O3_4363.21A": [10e-12, 10e-4],
    "O3_4958.91A": [10e-11, 10e-3],
    "O3_5006.84A": [10e-11, 10e-3], 
    "He2_1640.41A": [10e-9, 10e-4],
    "C2_1335.66A": [10e-12, 10e-3],
    "C3_1906.68A": [10e-11, 10e-2],
    "C3_1908.73A": [10e-10, 10e-3],
    "C4_1549.00A": [10e-13, 10e-7],
    "Mg2_2795.53A": [10e-11, 10e-2],
    "Mg2_2802.71A": [10e-11, 10e-2],
    "Ne3_3868.76A": [10e-12, 10e-4],
    "Ne3_3967.47A": [10e-12, 10e-4],
    "N5_1238.82A": [10e-8, 10e-3],
    "N5_1242.80A": [10e-9, 10e-4],
    "N4_1486.50A": [10e-12, 10e-3],
    "N3_1749.67A": [10e-11, 10e-4],
    "S2_6716.44A": [10e-11, 10e-3],
    "S2_6730.82A": [10e-11, 10e-3]
}

def sim_diagnostics(ds, sp, data_file, width):
    
    #galaxy_visualization.plot_diagnostics(ds, sp, data_file, star_ctr, width)
    #galaxy_visualization.plot_intensities(ds, sp, data_file, star_ctr, width)

    #galaxy_visualization.plot_diagnostics(ds, sp, data_file, star_ctr, width, lims_dict=lims_00273)
    #galaxy_visualization.plot_intensities(ds, sp, data_file, star_ctr, width, lims_dict=lims_00273)
    #galaxy_visualization.spectra_driver(ds, luminosities, data_file)
    
    galaxy_visualization.star_gas_overlay(ds, ad, sp, data_file, star_ctr, width, 'intensity_H1_6562.80A', lims_dict=lims_00273)

sim_diagnostics(ds, sp, sim_run, width)


#maxT = ad.max(('gas', 'temperature'))'
'''

<<<<<<< HEAD
# TODO phase plots
# TODO check other data saving
# TODO spectra
=======
# TODO checkerboard
>>>>>>> 0660d65 (Copy nebular_lines_v2 3-24-2025. Plan to make edits to temperature calculation.)
