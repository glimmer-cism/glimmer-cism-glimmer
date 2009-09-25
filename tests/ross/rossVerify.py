#!/usr/bin/env python
import pycdf
import sys
import matplotlib.pyplot as plt
import numpy
from scipy import interpolate
import matplotlib.colors as colors
import matplotlib.cm as cm
print interpolate
#Code to return a descritized version of a colormap
#From http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations 
def cmap_discretize(cmap, N):
     cdict = cmap._segmentdata.copy()
     # N colors
     colors_i = numpy.linspace(0,1.,N)
     # N+1 indices
     indices = numpy.linspace(0,1.,N+1)
     for key in ('red','green','blue'):
         # Find the N colors
         D = numpy.array(cdict[key])
         I = interpolate.interp1d(D[:,0], D[:,1])
         clrs = I(colors_i)
         # Place these colors at the correct indices.
         A = numpy.zeros((N+1,3), float)
         A[:,0] = indices
         A[1:,1] = clrs
         A[:-1,2] = clrs
         # Create a tuple for the dictionary.
         L = []
         for l in A:
             L.append(tuple(l))
         cdict[key] = tuple(L)
     # Return colormap object.
     cmap = colors.LinearSegmentedColormap('colormap',cdict,1024)
     return cmap



#Searches a sorted list for a value, and returns in a tuple:
#1. The index of the last value that is *less* than the given value
#2. An alpha value for interpolation
def searchList(l, value):
    for i in range(len(l) - 1):
        if l[i+1] > value:
            return i, (value - l[i]) / (l[i+1] - l[i]) 
    print "NOT FOUND!!!"
    print l
    print value


def interp(data, i, j, alpha1, alpha2):
    return (data[i, j]  *(1-alpha1) + data[i+1, j]  *alpha1) * (1-alpha2) + \
           (data[i, j+1]*(1-alpha1) + data[i+1, j+1]*alpha1) * alpha2

    

rows = 115
cols = 151
shape = (rows, cols)

delta_u = 30 #Assume constant error for the observed velocity, per MacAyeal '96

output_nc = pycdf.CDF("ross.out.nc")
raw_nc    = pycdf.CDF("ross-raw.nc")

griddatafile = open("111by147Grid.dat")

#Read the output normalized velocity
velnorm = output_nc.var("velnormhom")[0,0,2:-2,2:-2]
velraw = raw_nc.var("velo_mag")[0,2:-2,2:-2]
fake   = raw_nc.var("fake_shelf_mask")[0,2:-2,2:-2]
velraw = numpy.where(fake == 1, 0, velraw)


#Read the lat and lon
for i in range(4):
    griddatafile.readline()


#Read the row and column positions in the local coordinate system
row_positions = []
for i in range(111):
    row_positions.append(float(griddatafile.readline()))

for i in range(3):
    griddatafile.readline()

column_positions = []
for i in range(147):
    column_positions.append(float(griddatafile.readline()))

#Read the RIGGS data (thanks to Ed Bueler and PISM for providing this!)
riggsdata = open("riggs_clean.dat")
riggsLats = []
riggsLons = []
riggsVels = []
computedVels = []

t = 0
n = 0
for line in riggsdata:
    tokens = [float(s) for s in line.strip().split()]

    #Convert degrees, minutes, and seconds into lat and lon
    #in the local coordinate system
    lat = tokens[3] + tokens[4]/60 + tokens[5]/60**2
    lon = tokens[6] + tokens[7]/60 + tokens[8]/60**2
    
    lat = -lat
    lon *= -tokens[9]

    vel = tokens[10]

    if lat < row_positions[0] or lat > row_positions[-1] or lon < column_positions[0] or lon > column_positions[-1]:
        continue


    #Find a set of four points to interpolate to get the computed velocity
    #at this RIGGS station
    row, alpha1 = searchList(row_positions, lat)
    col, alpha2 = searchList(column_positions, lon)
    computedVel = interp(velnorm, row, col, alpha1, alpha2)
    rawVel = interp(velraw, row, col, alpha1, alpha2)


    #If our computed velocity is 0, assume that this is not in the domain
    if computedVel == 0.0:
        continue

    riggsLats.append(lat)
    riggsLons.append(lon)
    riggsVels.append(vel)
    computedVels.append(computedVel)

    n += 1
    t += (vel-computedVel)**2 / 30**2



print "X^2 =", t
print "Max Value =", numpy.max(velnorm)

print len(column_positions), len(row_positions)
print velnorm.shape

nbreaks = 12
breaks = numpy.linspace(0,1500, nbreaks)
norm = colors.Normalize(0,1500)
cmap = cm.get_cmap("gist_ncar")
cmap = cmap_discretize(cmap, nbreaks)
plt.subplot(1, 2, 1)

plt.contourf(column_positions, row_positions, velnorm, breaks, norm=norm, cmap=cmap)
plt.scatter(riggsLons, riggsLats, c=riggsVels, norm=norm, cmap=cmap)
plt.colorbar()

plt.subplot(1, 2, 2)
plt.scatter(riggsVels, computedVels)
plt.plot([0, 1200], [0, 1200])
plt.show()

