import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from matplotlib import colormaps
list(colormaps)

Zero=1e-40
SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Usage - select data to use
data=0

if data == 0:
  filename='nebula.linesPredictionO'
  titl=["H1 6562.80A","O1 1304.86A","O1 6300.30A","O2 3728.80A","O2 3726.10A","O3 1660.81A","O3 1666.15A","O3 4363.21A","O3 4958.91A","O3 5006.84A"]

if data == 1:
  filename='nebula.linesPredictionCN'
  titl=["He2 1640.41A","C2 1335.66A","C3 1906.68A","C3 1908.73A","C4 1549.00A","Mg2 2795.53A","Mg2 2802.71A","Ne3 3868.76A","Ne3 3967.47A","N5 1238.82A","N5 1242.80A","N4 1486.50A","N3 1749.67A","S2 6716.44A","S2 6730.82A"]

ncols=len(titl)
cols=np.arange(ncols)+2

# Read information about the run
# Read in the line list
minU,maxU,stepU,minN,maxN,stepN,minT,maxT,stepT=np.loadtxt(filename,unpack=True,dtype=float, max_rows=1)
ll=np.loadtxt(filename,unpack=True,dtype=float,usecols=cols,skiprows=1)
print(ll.shape)
print(minU,maxU,stepU)
print(minN,maxN,stepN)
print(minT,maxT,stepT)

#minU=-6.0
#maxU= 1.0
#stepU=0.5
#minN=-1.0
#maxN= 6.0 
#stepN=0.5
#minT= 3.0
#maxT= 6.0
#stepT=0.1

# U, T, N dimensions: number of indices
dimU=int((maxU-minU)/stepU)+1
dimT=int((maxT-minT)/stepT)+1
dimN=int((maxN-minN)/stepN)+1
print(dimU,dimN,dimT)

# the log values of U, n, T in the run/grid
logU=minU+np.arange(dimU)*stepU
logn=minN+np.arange(dimN)*stepN
logT=minT+np.arange(dimT)*stepT

# (Ionization Parameter, Density, Temperature)
# (U, density, T)
# d defines the cube dimensions
# 4D cube with ncols line strengths at each U, N, T coordinate
# cub[i] is the cube for a single emission line
# reshape the 1D array ll[i, :] of a certain line's strengths
# to U, n, T grid
d=(dimU,dimN,dimT)
cub=np.zeros((ncols,dimU,dimN,dimT))
for i in range(ncols):
  cub[i]=np.reshape(ll[i,:], d)
  #print(cub[i].shape)

# normalize by the density squared
for i in np.arange(dimN):
  for j in np.arange(ncols):
   cub[j,:,i,:]=cub[j,:,i,:]/10**(2*logn[i])
   #cub[j,:,i,:]=cub[j,:,i,:]

# Interpolation
# list of interpolators (for each emission line)
interpolator = [None]*ncols

for i in np.arange(ncols):
  interpolator[i] = RegularGridInterpolator((logU, logn, logT), cub[i])

# Fill a grid with interpolated values (double precision)
interpU = np.linspace(minU, maxU, dimU*2)
interpN = np.linspace(minN, maxN, dimN*2)
interpT = np.linspace(minT, maxT, dimT*2)
X, Y, Z = np.meshgrid(interpU, interpN, interpT, indexing='ij')

interpcub = np.zeros((ncols, dimU*2, dimN*2, dimT*2))

for i in np.arange(ncols):
  interpcub[i] = interpolator[i]((X, Y, Z))

print(interpcub.shape)
  
# Plotting routine
def implot(cub, logU, logN, stepU, stepN):
  fig,ax = plt.subplots(2,ncols, sharex=True,sharey='row', dpi=100, figsize=(12,4))
  #fig,ax = plt.subplots(2,ncols, dpi=100, figsize=(12,4))
  plt.subplots_adjust(left=0.05,
                      bottom=0.1,
                      right=0.99,
                      top=0.95,
                      wspace=0,
                      hspace=0)

  for i in range(ncols):
    ax[1][i].set_xlabel('log T')
    #  for j in range(2):
    #   num=(i+1)*(j+1)-1

    #slicel=np.log10(cub[i,:,0,:])
    #sliceh=np.log10(cub[i,:,14,:])
    #ratio=slicel-sliceh

    nidx = (logN-minN)/stepN
    uidx = (logU-minU)/stepU
    
    slicel=np.log10(cub[i,:,0,:]+Zero)
    sliceh=np.log10(cub[i,0,:,:]+Zero)
    ratio=slicel/sliceh

    vmax=np.max(slicel)
    #vmin=np.min(slicel)
    vmin=vmax-3.0

    vmax1=np.max(sliceh)
    #vmin1=np.min(sliceh)
    vmin1=vmax1-3.0
    
    ind=np.argmax(slicel,keepdims=True)
    ind = np.unravel_index(np.argmax(slicel, axis=None), slicel.shape)
    #print(f"{i:2d}, {titl[i]:14s}, Line strength= {vmax:.2e} logT= {logT[ind[1]]:.2f} log U= {logU[ind[0]]:.2f} {ind}")
    ex0=(minT-stepT/2.0, maxT+stepT/2.0,maxU+stepU/2.0,minU-stepU/2.0)
    ex1=(minT-stepT/2.0, maxT+stepT/2.0,maxN+stepN/2.0,minN-stepN/2.0)

    #cmap='gist_ncar'
    ax[0][i].imshow(slicel,extent=ex0,vmin=vmin,vmax=vmax)
    #ax[0][i].imshow(slicel,extent=ex0)

    ax[1][i].imshow(sliceh,extent=ex1,vmin=vmin1,vmax=vmax1)

    #ax[1][i].imshow(ratio,extent=ex)
    ax[0][0].set_ylabel('log U')
    ax[1][0].set_ylabel('log N')
    ax[0][i].set_title(titl[i])

#plt.plot(logT,cub1[2,0,:])
#plt.plot(logT,cub1[5,0,:])
#plt.plot(logT,cub1[10,0,:])
implot(cub, 0.0, -1.0, stepU, stepN)
implot(interpcub, 0.0, -1.0, stepU/2.0, stepN/2.0)

plt.show()

# top and bottom plots - what value to keep constant
# github
# ytproject
