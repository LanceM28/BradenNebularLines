# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:36:14 2023

@author: garli
"""
import matplotlib as mlt
mlt.use('TkAgg')
import numpy as np
from scipy import constants, integrate, special
import matplotlib.pyplot as plt
import yt
from yt import derived_field
import scipy.constants as con
from unyt import unyt_array 
import dill
from yt.frontends.ramses.field_handlers import RTFieldFileHandler

cgs_c=2.99792458e10

width=unyt_array([0.5, 0.5], 'kpc')

axis='y'
digits=2
lambda0=unyt_array([656.28], 'nm')
bandpass_vel=10
lambda1=(bandpass_vel*(1e3)/con.c)*lambda0
bandpass=unyt_array([lambda0.value-lambda1.value, lambda0.value+lambda1.value], 'nm')

yt.toggle_interactivity()

cell_fields = [

    "Density",

    "x-velocity",

    "y-velocity",

    "z-velocity",

    "Pressure",

    "metallicity",

    "H2-fraction",

    "HII-fraction",

    "HeII-fraction",

    "HeIII-fraction", 
    
    "refinement_param"

    ]

epf = [("p1", "double"), ("p2", "double"),("p3", "double"),("p4", "double")]

def _x0(field, data): 
    return data['index', 'x']-(data['index', 'dx']/2)

def _x1(field, data): 
    return data['index', 'x'] + (data['index', 'dx']/2)

def _y0(field, data): 
    return data['index', 'y'] -(data['index', 'dy']/2)

def _y1(field, data): 
    return data['index', 'y'] + (data['index', 'dy']/2)

def _z0(field, data): 
    return data['index', 'z']-(data['index', 'dz']/2)

def _z1(field, data): 
    return data['index', 'z'] + (data['index', 'dz']/2)

def _radial_vel(field, data): 
    if axis=='x': 
        y=data['ramses', 'x-velocity']
    elif axis=='y': 
        y=data['ramses', 'y-velocity']
    else: 
        y=data['ramses', 'z-velocity']
    return y

def _lambda_bar(field, data):
    lambda_bar=(data['gas', 'radial-vel']/(con.c*u.m/u.s) + 1)*lambda0 #NOTE: check sign!!!
    return lambda_bar

def _e_frac(field, data): 
    return (data['ramses', 'HII-fraction'] + (data['ramses', 'HeII-fraction'] + 2*data['ramses', 'HeIII-fraction']))

def _p_frac(field, data): 
    return data['ramses', 'HII-fraction'] + data['ramses', 'HeIII-fraction']

def _e_density(field, data): 
    return data['gas', 'e-fraction']*data['gas', 'number_density']

def _p_density(field, data):
    return data['gas', 'p-fraction']*data['gas', 'number_density']

def _m_bar(field, data): 
    return data['gas', 'density']/data['gas', 'number_density']

def _h_alpha(field, data):
    #remember that this field is power integrated over all wavelengths
    alpha=(2.6e-13)*(u.cm)**3*(1/u.s)*(data['gas', 'temperature']*(1e-4)/u.K)**(-0.7) #NOTE: change to include temperature!! 
    h_nu=con.c/(656e-9)*con.h*1e7*u.erg
    constant=0.45*h_nu*alpha
    return constant*data['gas', 'p-density']*data['gas', 'e-density']*data['index', 'volume']

def _h_alpha2(field, data):
    #remember that this field is power integrated over all wavelengths
    alpha=(2.6e-13)*(u.cm)**3*(1/u.s)*(data['gas', 'temperature']*(1e-4)/u.K)**(-0.7) #NOTE: change to include temperature!! 
    h_nu=con.c/(656e-9)*con.h*1e7*u.erg
    constant=0.45*h_nu*alpha
    return constant*data['gas', 'p-density']*data['gas', 'e-density']


def _bandpass_l(field, data): 
    sig=np.sqrt(con.k*u.J/u.K*data['gas', 'temperature']/data['gas', 'm-bar'])*lambda0/(con.c*u.m/u.s)
    A=data['gas', 'h-alpha2']
    
    x2=(bandpass[1]-data['gas', 'lambda-bar'])/(sig*np.sqrt(2))
    x1=(bandpass[0]-data['gas', 'lambda-bar'])/(sig*np.sqrt(2))

    flux=(A/2)*(special.erf(x2.value)-special.erf(x1.value))

    return flux

def _h_alpha3(field, data):
    #luminosity density per wavelength per volume
    sigma=np.sqrt(con.k*u.J/u.K*data['gas', 'temperature']/data['gas', 'm-bar'])*lambda0/(con.c*u.m/u.s)
    factor = 1/(np.sqrt(2*np.pi)*sigma)

    return factor*data['gas', 'h-alpha2']

def _ion_param(field, data): 
    photon=data['ramses-rt', 'Photon_density_1']+data['ramses-rt', 'Photon_density_2'] + data['ramses-rt', 'Photon_density_3'] + data['ramses-rt', 'Photon_density_4']

    return photon/data['gas', 'number_density']
"""
def _p_dens_1(field, data): 
    return data['ramses-rt', 'Photon_density_1']*p['unit_pf']/cgs_c
def _p_dens_2(): 
    return data['ramses-rt', 'Photon_density_2']*p['unit_pf']/cgs_c
def _p_dens_3(): 
    return data['ramses-rt', 'Photon_density_3']*p['unit_pf']/cgs_c

yt.add_field( 
    ("gas", "p-dens-1"), 
    function=_p_dens_1, 
    units="cm**(-3)", 
    sampling_type="cell"
)
"""
yt.add_field(
    ('gas', 'ion-param'), 
    function=_ion_param, 
    sampling_type="cell", 
    units="cm**3", 
    force_override=True
)



yt.add_field(
    ('gas', 'h-alpha3'), 
    function=_h_alpha3, 
    sampling_type="cell", 
    units="ergs/(s*cm**(4))", 
    force_override=True
)
yt.add_field(
    ('gas', 'radial-vel'),
    function=_radial_vel, 
    sampling_type="cell",
    units="m/s",
    force_override=True
)


yt.add_field(
    ('gas', 'lambda-bar'), 
    function=_lambda_bar, 
    sampling_type="cell", 
    units='cm', 
    force_override=True
)


yt.add_field(
    ('gas', 'l-bandpass'),
    function=_bandpass_l, 
    sampling_type="cell",
    units="erg/s/(cm**3)", 
    force_override=True
)




yt.add_field(
    ('gas', 'radial-vel'),
    function=_radial_vel, 
    sampling_type="cell",
    units="m/s",
    force_override=True
)

yt.add_field(
    ('gas', 'lambda-bar'), 
    function=_lambda_bar, 
    sampling_type="cell", 
    units='cm', 
    force_override=True
)

yt.add_field(
    ('gas', 'e-fraction'),
    function=_e_frac, 
    sampling_type="cell",
    units="dimensionless",
    force_override=True
)

yt.add_field(
    ('gas', 'e-density'),
    function=_e_density, 
    sampling_type="cell",
    units="cm**(-3)",
    force_override=True
)

yt.add_field(
    ('gas', 'p-fraction'),
    function=_p_frac, 
    sampling_type="cell",
    units="dimensionless",
    force_override=True
)

yt.add_field(
    ('gas', 'p-density'),
    function=_p_density, 
    sampling_type="cell",
    units="cm**(-3)",
    force_override=True
)

yt.add_field(
    ('gas', 'm-bar'), 
    function=_m_bar, 
    sampling_type="cell", 
    units="auto", 
    force_override=True
    )

yt.add_field(
    ('gas', 'h-alpha'),
    function=_h_alpha, 
    sampling_type="cell",
    units="erg/s", 
    force_override=True
    )
yt.add_field(
    ('gas', 'h-alpha2'),
    function=_h_alpha2, 
    sampling_type="cell",
    units="erg/s/(cm**3)", 
    force_override=True
    )

yt.add_field(
    ('gas', 'bandpass-l'),
    function=_bandpass_l, 
    sampling_type="cell",
    units="erg/s/(cm**3)", 
    force_override=True
    )


yt.add_field(
    ('index', 'x0'), 
    function=_x0, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )

yt.add_field(
    ('index', 'x0'), 
    function=_x0, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )


yt.add_field(
    ('index', 'x1'), 
    function=_x1, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )
yt.add_field(
    ('index', 'y0'), 
    function=_y0, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )
yt.add_field(
    ('index', 'y1'), 
    function=_y1, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )
yt.add_field(
    ('index', 'z0'), 
    function=_z0, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )

yt.add_field(
    ('index', 'z1'), 
    function=_z1, 
    sampling_type="cell",
    units="auto",
    force_override=True
    )


def gauss_func(x, avg, std): 
    const=1/(np.sqrt(2*np.pi)*std)
    exp=np.exp(-(((x-avg)/(np.sqrt(2)*std))**2))
    return const*exp

def axis_to_no(axis):
    if axis=='x':
        no=0
    elif axis=='y': 
        no=1
    elif axis=='z':
        no=2
    else:
         raise ValueError('Axis can only equal one of the following: x, y, z')
    return no

def find_spectra_lim(axis_no, center, width, lim):
     spectra_lim=unyt_array(np.empty((3, 2)), 'cm')
     
     spectra_lim[axis_no, 0]=unyt_array.min(lim[axis_no, 0])
     spectra_lim[axis_no, 1]=unyt_array.max(lim[axis_no, 1])
     
     n=0
     
     for s in range(spectra_lim.shape[0]): 
         if s!=axis_no:
             spectra_lim[s, 0]=center[s]-width[n]
             spectra_lim[s, 1]=center[s]+width[n]
             n=n+1
    
    
     return spectra_lim
 
def find_frac(spectra_lim, lim):
    
    '''
    This function gives the fraction of the given "box" of
    volume from yt, and then it gives it the fraction of that
    box that falls within the given spectral limits
    '''
    frac=unyt_array(np.empty((3, lim.shape[-1])), 'dimensionless')
    
    for i in range(spectra_lim.shape[0]):
        
        c0=(lim[i, 1]<=spectra_lim[i, 0]) | (lim[i, 0]>=spectra_lim[i, 1])
        #c1=(lim[i, 0]>spectra_lim[i, 0]) & (spectra_lim[i, 1]>lim[i, 1])
        #c1=((spectra_lim[i, 0]<lim[i, 0]) & (lim[i, 1]<spectra_lim[i, 1])) | ((spectra_lim[i, 0]==lim[i, 0]) & (spectra_lim[i, 1]==lim[i, 1])) 
        c1=( (spectra_lim[i, 0]<lim[i, 0]) | (spectra_lim[i, 0] == lim[i, 0])) & ( (lim[i, 1]<spectra_lim[i, 1]) | (spectra_lim[i, 1]==lim[i, 1]))
        
        c4=( (lim[i, 0] < spectra_lim[i, 0]) | (lim[i, 0] == spectra_lim[i, 0]) )  &  ((spectra_lim[i, 1]<lim[i,1]) | (lim[i, 1] == spectra_lim[i, 1]) )
        
        c2= (~c0) & (~c4) & (lim[i, 0]<spectra_lim[i, 0])
        c3= (~c0) & (~c4) & (lim[i, 1]> spectra_lim[i, 1])
        
        
        frac[i, c0]=0
        frac[i, c1]=1
        frac[i, c2]=(lim[i, 1, c2]-spectra_lim[i, 0])/(lim[i, 1, c2]-lim[i, 0, c2])
        frac[i, c3]=(spectra_lim[i, 1]-lim[i, 0, c3])/(lim[i, 1, c3]-lim[i, 0, c3])
        frac[i, c4]=(spectra_lim[i, 1]-spectra_lim[i, 0])/(lim[i, 1, c4]-lim[i, 0, c4])
        

        test1= c0 | c1 | c2 | c3 | c4
        index=np.where(test1==False)
        print(index)
        #print(np.min( c0 | c1 | c2 | c3 | c4 ))
        #print( np.max( c0 & c1 ))
        #print(np.max( c1 & c2 ))
        
    
    vol_frac=frac[0]*frac[1]*frac[2]
    return vol_frac

def trim_array(arr, vol_frac):
    condition=np.array(vol_frac==0)
    return arr[~condition]

def trim_array2(arr, int_flux, vol_frac, vol): 
    #condition=np.array(h_alpha>1e-8*np.sum(h_alpha))
    
    l=int_flux*vol*vol_frac

    cumsum=np.cumsum(np.sort(l))

    index=np.argmin(abs(cumsum-(1e-4)*np.sum(l)))
    
    sorted_i=np.argsort(l)[index:]
    
    #NOTE: returns stuff sorted by luminosity
    return arr[sorted_i]

def rad_vel(axis_no, vel): 
    radial=vel[axis_no]
    return radial

def get_lambda_bar(radial, lambda0): 
    lambda_bar=(radial/(con.c*u.m/u.s) + 1)*lambda0 #NOTE: check sign!!!
    return lambda_bar

def get_spectra(lambda_bar, lambda0, temp, m_bar, lambda_arr, h_alpha, vol_frac): 
    sigma=np.sqrt(con.k*u.J/u.K*temp/m_bar)*lambda0/(con.c*u.m/u.s)
    spectra=gauss_func(lambda_arr, lambda_bar, sigma)*h_alpha*vol_frac
    return spectra, sigma

ds = yt.load("/home/jcottin1/research_summer2023/Job2.2.2/output_01410/info_01410.txt", extra_particle_fields=epf)
u=ds.units
#fields=RTFieldFileHandler.detect_fields(ds).copy()
#p=RTFieldFileHandler.get_rt_parameters(ds).copy()
#p.update(ds.parameters)
s=ds.all_data()
v, c = ds.find_max(('gas', 'density'))

temp=ds.all_data()['gas', 'temperature']
h_alpha=ds.all_data()['gas', 'h-alpha']
m_bar=ds.all_data()['gas', 'm-bar'] #NOTE: fixed ratio of H to He, so this doesn't have to be an array
vol=ds.all_data()['index', 'volume']
int_f=ds.all_data()['gas', 'bandpass-l']

lim=unyt_array(np.empty((3, 2, *temp.shape)), 'cm')
vel=unyt_array(np.empty((3, *temp.shape)), 'm/s')

lim[0, 0, :]=ds.all_data()['index', 'x0']
lim[0, 1, :]=ds.all_data()['index', 'x1']
lim[1, 0, :]=ds.all_data()['index', 'y0']
lim[1, 1, :]=ds.all_data()['index', 'y1']
lim[2, 0, :]=ds.all_data()['index', 'z0']
lim[2, 1, :]=ds.all_data()['index', 'z1']

vel[0, :]=ds.all_data()['ramses', 'x-velocity']
vel[1, :]=ds.all_data()['ramses', 'y-velocity']
vel[2, :]=ds.all_data()['ramses', 'z-velocity']

axis_no=axis_to_no(axis)


s_lim=find_spectra_lim(axis_no, c, width, lim)

print('Finding the fraction of each cell within line of sight')
vol_frac=find_frac(s_lim, lim)

print('Finding radial velocity...')
radial=rad_vel(axis_no, vel)

print('Trimming arrays...')
"""
vol=trim_array(vol, vol_frac)
radial=trim_array(radial, vol_frac)
h_alpha=trim_array(h_alpha, vol_frac)
m_bar=trim_array(m_bar, vol_frac)
temp=trim_array(temp, vol_frac)
vol_frac=trim_array(vol_frac, vol_frac)
"""
#h_alpha2=h_alpha
vol2=vol
vol_frac2=vol_frac

vol=trim_array2(vol, int_f, vol2, vol_frac2)
radial=trim_array2(radial, int_f, vol2, vol_frac2)
vol_frac=trim_array2(vol_frac, int_f, vol2, vol_frac2)
m_bar=trim_array2(m_bar, int_f, vol2, vol_frac2)
temp=trim_array2(temp, int_f, vol2, vol_frac2)
h_alpha=trim_array2(h_alpha, int_f, vol2, vol_frac2)

#del(vol_frac2)
#del(h_alpha2)

print('Finding average wavelength of spectra')
lambda_bar=get_lambda_bar(radial, lambda0)

lambda_arr=unyt_array(np.arange(655, 657, 0.01), 'nm')[:, np.newaxis]
print('Calculating spectra')

#spectra, sig =get_spectra(lambda_bar, lambda0, temp, m_bar, lambda_arr, h_alpha, vol_frac)
#spectra=spectra.to('erg/(nm*s)')
#summed_spectra=np.sum(spectra, axis=1)

#remember to make map of peak height! 
#needed_var = np.array(['lambda_arr', 'lambda0', 'width',
"""
'digits', 'lambda_bar', 'axis', 'c', 'vol_frac', 'np', 'dill', 'h_alpha',
'spectra', 'radial', 'cell_fields', 'epf', '_m_bar', '_e_frac',
'_p_frac', '_e_density', '_p_density', '_h_alpha', 's_lim', '_bandpass_l', 'bandpass', '_lambda_bar', '_radial_vel'])
"""
#all_stuff=list(locals())
"""
for i in needed_var: 
    if np.sum(all_stuff==i): 
        dill.dumps(i, 'processed_data/dir1/'+ i + '.pkl')
"""
"""
for i in all_stuff: 
    if np.sum(i==needed_var)==0:
            if not i.startswith('__'):
                del locals()[i]
        



del(lim)
del(temp)
del(vol)
del(vel)
del(spectra)
del(sig)
del(vol_frac)
"""

#dill.dump_session('processed_data/pickle1.pkl')
"""
for n in range(num): 
    plt.plot(vel_arr, spectra[:, j[n]], label='Block number ' + str(n))
    box_text=box_text + '\nBlock number ' + str(n) + '\nLuminosity: ' + sci(weights[j[n]].to('erg/s')) + '\nAverage: ' + sci(lambda_bar[j[n]]) + '\nStdev: ' + sci(sig[j[n]].to('nm'))  

box = mlt.offsetbox.AnchoredText(box_text, loc=2, prop=dict(fontsize=6))
ax.add_artist(box)
print(box_text)
plt.legend()
"""

#Histogram of radial velocity (weighted by strength of h_alpha line) 
unit='km/s'
plt.figure()
a=np.array(radial.to(unit))
total_vol=np.prod(s_lim[:, 1]-s_lim[:, 0])
weight_arr=h_alpha*vol_frac
weight_arr=weight_arr/np.sum(weight_arr)
plt.hist(a, bins=100, weights=weight_arr)
#plt.xlim([-100, 100])
plt.xlabel('Radial velocity (viewed down the ' + axis + ' axis) in ' + unit)
lambda0=unyt_array([656.28], 'nm')
plt.ylabel('Count')
plt.title('Velocity down the ' + axis + ' axis weighted by h alpha')
plt.yscale('log')


unit='nm' 
b=np.array(lambda_bar.to(unit))
plt.figure()
plt.hist(b, bins=100, weights=weight_arr)
plt.xlabel('Average wavelength (viewed down the ' + axis + ' axis) in ' + unit)
plt.ylabel('Count')
plt.title('Wavelength down the ' + axis + ' axis weighted by h alpha')
plt.yscale('log')

fig, ax = plt.subplots()
var3=(b/lambda0.value-1)*con.c*10**(-3)
unit='km/s'
plt.hist(var3, bins=100, weights=weight_arr)
plt.xlabel('Average wavelength (viewed down the ' + axis + ' axis) in ' + unit)
plt.ylabel('Count')
plt.title('Wavelength down the ' + axis + ' axis weighted by h alpha')
plt.yscale('log')
#.ticklabel_format(axis='x', style='sci')

"""
unit='nm' 
b=np.array(sig.to(unit))
unit='km/s'
c2=b/lambda0.value*con.c*10**(-3)
plt.figure()
plt.xlim([1, 1000])
plt.hist(c2, bins=100, weights=weight_arr, log=True)
plt.xlabel('Doppler broadening in ' + unit)
plt.ylabel('Count')
#plt.xscale('log')
plt.yscale('log')
plt.title('Distribution of doppler broadening weighted by luminosity of h-alpha line')
#plt.ticklabel_format(style='sci')
"""
fig, (ax1, cax1) =plt.subplots(1, 2)
fig2, (ax2, cax2) = plt.subplots(1, 2)
fig3, (ax3, cax3) = plt.subplots(1, 2)

fig4, (ax4, cax4) = plt.subplots(1, 2)

plot5=yt.ProjectionPlot(ds, axis, ('gas', 'h-alpha2'), buff_size=(800, 800), width=(1, 'kpc'), center=c)
plot6=yt.ProjectionPlot(ds, axis, ('ramses', 'HII-fraction'),
weight_field=('gas', 'density'), buff_size=(800, 800), width=(1, 'kpc'), center=c)


#ax1.colorbar(plot5)
#ax2.colorbar(plot6)
plot7=yt.ProjectionPlot(ds, axis, ('gas', 'l-bandpass'), buff_size=(800, 800), width=(1, 'kpc'), center=c)
plot8=yt.PhasePlot(ds.all_data(), ('gas', 'density'), ('gas', 'temperature'), [('gas', 'h-alpha2')], weight_field=('index', 'volume'))
plot8.save('h-alpha.png')
plot9=yt.PhasePlot(ds.all_data(), ('gas', 'density'), ('gas', 'temperature'), [('gas', 'ion-param')])
plot9.save('ion-param1.png')
plot10=yt.PhasePlot(ds.all_data(), ('gas', 'h-alpha2'), ('gas', 'temperature'), [('gas', 'ion-param')])
plot10.save('ion-param2.png')
#plot5.set_cmap(('gas', 'h-alpha2'), 'viridis')

#plot7.set_cmap(('gas', 'l'), 'viridis')

plot=plot5.plots["h-alpha2"]
plot2=plot6.plots["HII-fraction"]
plot3=plot7.plots["l-bandpass"]
plot4=plot8.plots["h-alpha2"]

plot.figure = fig
plot.axes = ax1
plot.cax=cax1

plot2.figure=fig2
plot2.axes=ax2
plot2.cax=cax2

plot3.figure=fig3
plot3.axes=ax3
plot3.cax=cax3

plot4.figure=fig4
plot4.axes=ax4
plot4.cax=cax4

plot5._setup_plots()
plot6._setup_plots()
plot7._setup_plots()
plot8._setup_plots()

def plot(plot, field): 
    
    plot2=plot.plots[field]
    
    fig, (ax1, cax1) =plt.subplots(1, 2)
    plot2.figure=fig
    plot2.axes=ax1
    plot2.cax=cax1

    plot._setup_plots()

plot7=yt.ProjectionPlot(ds, axis, ('gas', 'h-alpha3'), buff_size=(800, 800), width=(1, 'kpc'), center=c)
plot(plot7, 'h-alpha3')


plot12=yt.ProjectionPlot(ds, axis, ('gas', 'ion-param'), buff_size=(800, 800), width=(1, 'kpc'), center=c)
plot(plot12, 'ion-param')

n=1
plot5_frb = plot5.data_source.to_frb((n, "kpc"), [800, 800], center=c)['gas', 'h-alpha2']


plot5_img=np.flipud(np.array(plot5_frb))

plt.figure()
plt.imshow(plot5_img, norm='log', extent=[-n/2, n/2, -n/2, n/2])
plt.colorbar(label='Luminosity of h-alpha line (erg/(s cm^2))')
#plt.xlim([-n/2, n/2])
#plt.ylim([-n/2, n/2])
plt.ylabel('Distance on x axis (kpc)') #x axis
plt.xlabel('Distance on z axis (kpc)')



print('Waiting for plots to appear...')
plt.show(block=False)


