import os
import sys
import copy

import numpy as np
import yt
from yt.frontends.ramses.field_handlers import RTFieldFileHandler
import matplotlib.pyplot as plt

from emission import EmissionLineInterpolator
import galaxy_visualization

filename = "/Users/lamoreau/python/ASpec/SimulationFiles/output_00273/info_00273.txt"

lines=["H1_6562.80A","O1_1304.86A","O1_6300.30A","O2_3728.80A","O2_3726.10A",
       "O3_1660.81A","O3_1666.15A","O3_4363.21A","O3_4958.91A","O3_5006.84A", 
       "He2_1640.41A","C2_1335.66A","C3_1906.68A","C3_1908.73A","C4_1549.00A",
       "Mg2_2795.53A","Mg2_2802.71A","Ne3_3868.76A","Ne3_3967.47A",
       "N5_1238.82A",
       "N5_1242.80A","N4_1486.50A","N3_1749.67A","S2_6716.44A","S2_6730.82A"]

wavelengths=[6562.80, 1304.86, 6300.30, 3728.80, 3726.10, 1660.81, 1666.15,
             4363.21, 4958.91, 5006.84, 1640.41, 1335.66,
             1906.68, 1908.73, 1549.00, 2795.53, 2802.71, 3868.76,
             3967.47, 1238.82, 1242.80, 1486.50, 1749.67, 6716.44, 6730.82]

cell_fields = [
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
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


#number density of hydrogen atoms
def _my_H_nuclei_density(field, data):
    dn=data["ramses","Density"].in_cgs()
    XH_RAMSES=0.76 #defined by RAMSES in cooling_module.f90
    YHE_RAMSES=0.24 #defined by RAMSES in cooling_module.f90
    mH_RAMSES=yt.YTArray(1.6600000e-24,"g") #defined by RAMSES in cooling_module.f90

    return dn*XH_RAMSES/mH_RAMSES

def OII_ratio(field, data):
    # TODO lum or flux?
    return data['gas', 'luminosity_O2_3728.80A']/data['gas', 'luminosity_O2_3726.10A']


def _pressure(field, data):
    if 'hydro_thermal_pressure' in dir(ds.fields.ramses): # and 
        #'Pressure' not in dir(ds.fields.ramses):
        return data['ramses', 'hydro_thermal_pressure']


def _xHI(field, data):
    if 'hydro_xHI' in dir(ds.fields.ramses): # and \
        #'xHI' not in dir(ds.fields.ramses):
        return data['ramses', 'hydro_xHI']


def _xHII(field, data):
    if 'hydro_xHII' in dir(ds.fields.ramses): # and \
        #'xHII' not in dir(ds.fields.ramses):
        return data['ramses', 'hydro_xHII']


def _xHeII(field, data):
    if 'hydro_xHeII' in dir(ds.fields.ramses): # and \
        #'xHeII' not in dir(ds.fields.ramses):
        return data['ramses', 'hydro_xHeII']


def _xHeIII(field, data):
    if 'hydro_xHeIII' in dir(ds.fields.ramses): # and \
        #'xHeIII' not in dir(ds.fields.ramses):
        return data['ramses', 'hydro_xHeIII']

'''
-------------------------------------------------------------------------------
Load Simulation Data
Add Derived Fields
-------------------------------------------------------------------------------
'''

ds = yt.load(filename, extra_particle_fields=epf)

ds.add_field(
    ("gas","number_density"),
    function=_my_H_nuclei_density,
    sampling_type="cell",
    units="1/cm**3",
    force_override=True
)


ds.add_field(
    ("ramses","Pressure"),
    function=_pressure,
    sampling_type="cell",
    units="1",
    #force_override=True
)

ds.add_field(
    ("ramses","xHI"),
    function=_xHI,
    sampling_type="cell",
    units="1",
    #force_override=True
)

ds.add_field(
    ("ramses","xHII"),
    function=_xHII,
    sampling_type="cell",
    units="1",
    #force_override=True
)

ds.add_field(
    ("ramses","xHeII"),
    function=_xHeII,
    sampling_type="cell",
    units="1",
    #force_override=True
)

ds.add_field(
    ("ramses","xHeIII"),
    function=_xHeIII,
    sampling_type="cell",
    units="1",
    #force_override=True
)

ds.add_field(
    ("gas","my_temperature"),
    function=_my_temperature,
    sampling_type="cell",
    # TODO units
    #units="K",
    #units="K*cm**3/erg",
    units='K*cm*dyn/erg',
    force_override=True
)

# Ionization parameter
ds.add_field(
    ('gas', 'ion_param'),
    function=_ion_param,
    sampling_type="cell",
    units="cm**3",
    force_override=True
)

ds.add_field(
    ("gas","my_H_nuclei_density"),
    function=_my_H_nuclei_density,
    sampling_type="cell",
    units="1/cm**3",
    force_override=True
)

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
    ds.add_field(
        ('gas', 'flux_' + line),
        function=emission_interpolator.get_line_emission(
            i, dens_normalized=dens_normalized
        ),
        sampling_type='cell',
        units=units,
        force_override=True
    )
    # TODO change get_line_emission to accept line not idx

    ds.add_field(
        ('gas', 'luminosity_' + line),
        function=emission_interpolator.get_luminosity(lines[i]),
        #function=get_luminosity(lines[i]),
        sampling_type='cell',
        units='1/cm**3',
        force_override=True
    )

ad = ds.all_data()
print(ds.field_list)


'''
-------------------------------------------------------------------------------
Run routines on data
-------------------------------------------------------------------------------
'''

viz = galaxy_visualization.VisualizationManager(filename, lines, wavelengths)
star_ctr = viz.star_center(ad)
sp = ds.sphere(star_ctr, (3000, "pc"))
sp_lum = ds.sphere(star_ctr, (10, 'kpc'))
width = (1500, 'pc')

sim_run = filename.split('/')[-1]

field_list = [
    #('gas', 'temperature'),
    ('gas', 'density'),
    ('gas', 'my_H_nuclei_density'),
    ('gas', 'my_temperature'),
    ('gas', 'ion_param'),
    ('gas', 'metallicity'),
]

weight_field_list = [
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
#viz.plot_wrapper(ds, sp, width, star_ctr, field_list,
#                     weight_field_list, title_list, proj=True, slc=False)


#viz.calc_luminosities(sp)

'''
viz.plot_cumulative_field(ds, sp, width, star_ctr, ('gas', 'luminosity_H1_6562.80A'),
                              None, 'luminosity_H1_6562.80A', 'luminosity_H1_6562,80A_cumulative')

p = viz.proj_plot(ds, sp, width, star_ctr, ('gas', 'flux_H1_6562.80A'), None)
        
plot_frb = p.frb
p_img = np.array(plot_frb['gas', 'flux_H1_6562.80A'])

dim1, dim2 = p_img.shape

A_pt = (width[0] * 3.086e18 / dim1) * (width[0] * 3.086e18 / dim2)

lum = np.sum(p_img * A_pt)
print(lum)
'''

extrema = {('gas', 'my_temperature'): (1e3, 1e8),
           ('gas', 'my_H_nuclei_density'): (1e-4, 1e6)}

line_title = r'H$\alpha$_6562.80A'

viz.phase_plot(ds, sp, x_field=('gas', 'my_temperature'),
               y_field=('gas', 'my_H_nuclei_density'), z_field=('gas', 'luminosity_O2_3726.10A'),
               extrema=extrema, x_label='Temperature [K]', 
               y_label=r'H Nuclei Number Density [cm$^{-3}$]', 
               z_label=line_title.replace('_', ' ') + 
                      r' Luminosity [erg s$^{-1}$]')


'''
pix = sp[('gas', 'luminosity_H1_6562.80A')].value
pix_sort = np.sort(pix, axis=None)[::-1]
idxs = np.arange(0, len(pix_sort), 1)
cum_val = np.cumsum(pix_sort) / np.sum(pix_sort)

fig = plt.figure(figsize=(8, 6))
plt.xlabel('Index')
plt.ylabel('Cumulative Value')
plt.title(f'Cumulative Sum')
plt.plot(idxs, cum_val)
plt.grid(True)
plt.show()
plt.close()
'''

# width 1500 pc




# TODO line ratio - diagnostic O3 electron density
# colorbar on phase
# phase diagram - two images, ratio; or ratio of total luminosities
# compare e- density to phase diagram

# 1D histogram
# normalize value in colorbar
# min and max - gradient in yellow part 


extrema = {('gas', 'my_temperature'): (1e3, 1e8),
           ('gas', 'my_H_nuclei_density'): (1e-4, 1e6)}

line_title = r'H$\alpha$_6562.80A'

phase_profile, x_vals, y_vals, z_vals = viz.phase_plot(ds, sp, x_field=('gas', 'my_temperature'),
               y_field=('gas', 'my_H_nuclei_density'), z_field=('gas', 'flux_H1_6562.80A'),
               extrema=extrema, x_label='Temperature [K]', 
               y_label=r'H Nuclei Number Density [cm$^{-3}$]', 
               z_label=line_title.replace('_', ' ') + 
                      r' Flux [erg s$^{-1}$ cm$^{-2}$]')

viz.phase_with_profiles(ds, sp, phase_profile, x_field=('gas', 'my_temperature'),
                        y_field=('gas', 'my_H_nuclei_density'),
                        z_field=('gas', 'flux_H1_6562.80A'),
                        x_vals=x_vals, y_vals=y_vals, z_vals=z_vals,
                        x_label='Temperature [K]',
                        y_label=r'H Nuclei Number Density [cm$^{-3}$]',
                        z_label=line_title.replace('_', ' ') + 
                            r' Flux [erg s$^{-1}$ cm$^{-2}$]')
