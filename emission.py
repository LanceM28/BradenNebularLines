import numpy as np
import matplotlib as plt
from scipy.interpolate import RegularGridInterpolator
import copy
import yt

# Read line emission data (line list, run params)
filename='CloudyFiles/linelist.dat'
minU,maxU,stepU,minN,maxN,stepN,minT,maxT,stepT=np.loadtxt(filename,unpack=True,dtype=float, max_rows=1, skiprows=5)
print(minU,maxU,stepU,minN,maxN,stepN,minT,maxT,stepT)

ll=np.loadtxt(filename,unpack=True,dtype=float,skiprows=7)
print(ll.shape)

titls=["H1 6562.80A","O1 1304.86A","O1 6300.30A","O2 3728.80A","O2 3726.10A","O3 1660.81A",
       "O3 1666.15A","O3 4363.21A","O3 4958.91A","O3 5006.84A", "He2 1640.41A","C2 1335.66A",
       "C3 1906.68A","C3 1908.73A","C4 1549.00A","Mg2 2795.53A","Mg2 2802.71A","Ne3 3868.76A",
       "Ne3 3967.47A","N5 1238.82A","N5 1242.80A","N4 1486.50A","N3 1749.67A","S2 6716.44A","S2 6730.82A"]

# Number of emission lines
ncols=len(titls)

# Set the line to visualize
lineidx=0

# Reconfigure linelist into a data cube
dimU=int((maxU-minU)/stepU)+1
dimT=int((maxT-minT)/stepT)+1
dimN=int((maxN-minN)/stepN)+1
print(dimU,dimN,dimT)

# the log values of U, n, T in the run/grid
logU=minU+np.arange(dimU)*stepU
logN=minN+np.arange(dimN)*stepN
logT=minT+np.arange(dimT)*stepT

# (Ionization Parameter, Density, Temperature)
# (U, density, T)
# d defines the cube dimensions
# 4D cube with ncols line strengths at each U, N, T coordinate
# cub[i] is the cube for a single emission line
# reshape the 1D array ll[i, :] of a certain line's strengths
# to U, N, T grid
d=(dimU,dimN,dimT)
cub=np.zeros((ncols,dimU,dimN,dimT))
for i in range(ncols):
    cub[i]=np.reshape(ll[i,:], d)
    #print(cub[i].shape)

# Interpolation
# list of interpolators (for each emission line)
interpolator = [None]*ncols

for i in np.arange(ncols):
  interpolator[i] = RegularGridInterpolator((logU, logN, logT), cub[i])

# normalize by the density squared
dens_normalized_cub = cub.copy()

for i in np.arange(dimN):
    for j in np.arange(ncols):
      dens_normalized_cub[j,:,i,:]=dens_normalized_cub[j,:,i,:]/10**(2*logN[i])
      #dens_normalized_cub[j,:,i,:]=dens_normalized_cub[j,:,i,:]

# Density Squared Normalized Interpolators
dens_normalized_interpolator = [None]*ncols

for i in np.arange(ncols):
  dens_normalized_interpolator[i] = RegularGridInterpolator((logU, logN, logT), cub[i])

# Get an interpolator which either returns the intensity erg cm^-2 s^-1
# or normalized by the density squared
def get_interpolator(lineidx, dens_normalized):
    if dens_normalized:
       return dens_normalized_interpolator[lineidx]
    return interpolator[lineidx]

# Returns a function for line emission of index idx;
# Allows for the batch creation of intensity fields
# for a variety of lines
def get_line_emission(idx, dens_normalized):
    def _line_emission(field, data):
        interpolator=get_interpolator(idx, dens_normalized)

        # Change to log values
        U_val = data['gas', 'ion-param'].value
        N_val = data['gas', 'number_density'].value
        T_val = data['gas', 'temperature'].value

        # Cut off negative temperatures
        T_val = np.where(T_val < 0.0, 10e-4, T_val)

        U = np.log10(U_val)
        N = np.log10(N_val)
        T = np.log10(T_val)

        # Adjust log values to within bounds supported by
        # interpolation table
        Uadj = np.where(U < minU, minU, U)
        Uadj = np.where(Uadj > maxU, maxU, Uadj)

        Nadj = np.where(N < minN, minN, N)
        Nadj = np.where(Nadj > maxN, maxN, Nadj)

        Tadj = np.where(T < minT, minT, T)
        Tadj = np.where(Tadj > maxT, maxT, Tadj)
    
        tup = np.stack((Uadj, Nadj, Tadj), axis=-1)

        # Return interpolated values weighted by metallicity
        # for non-Hydrogen and Helium lines
        interp_val = interpolator(tup)

        if idx not in [0, 10]:
           interp_val = interp_val*data['gas', 'metallicity']

        if dens_normalized:
           # make sure to weight field with number_density as well
           interp_val = interp_val*data['gas', 'number_density'] 
           #interp_val = interp_val*data['gas', 'number_density']**2

        return interp_val
    return copy.deepcopy(_line_emission)

