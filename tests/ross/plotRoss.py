#!/usr/bin/env python
import pycdf as pycdf
import numpy as np
from pycdf import NC as NC
from pylab import clf,colorbar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
from getopt import gnu_getopt


#Code to return a descritized version of a colormap
#From http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations 
def cmap_discretize(cmap, N):
     """Return a discrete colormap from the continuous colormap cmap.     
         cmap: colormap instance, eg. cm.jet. 
         N: Number of colors.
     
      Example
         x = resize(arange(100), (5,100))
         djet = cmap_discretize(cm.jet, 5)
         imshow(x, cmap=djet)
     """
     cdict = cmap._segmentdata.copy()
     # N colors
     colors_i = linspace(0,1.,N)
     # N+1 indices
     indices = linspace(0,1.,N+1)
     for key in ('red','green','blue'):
         # Find the N colors
         D = array(cdict[key])
         I = interpolate.interp1d(D[:,0], D[:,1])
         colors = I(colors_i)
         # Place these colors at the correct indices.
         A = zeros((N+1,3), float)
         A[:,0] = indices
         A[1:,1] = colors
         A[:-1,2] = colors
         # Create a tuple for the dictionary.
         L = []
         for l in A:
             L.append(tuple(l))
         cdict[key] = tuple(L)
     # Return colormap object.
     return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


optlist, args = gnu_getopt(sys.argv, "", ["nrows=", "log=", "invert"])
optdict = dict(optlist)

if "--nrows" in optdict:
    nrows = int(optdict["--nrows"])
else:
    nrows = 1

#Get a list of variables that need a logarithmic transform
if "--log" in optdict:
    logVars = optdict["--log"].split(',')
else:
    logVars = []
print logVars

#Whether we need to invert the vertical
invert = ("--invert", '') in optlist

infilename = args[1]
vars = args[2:]

ncols = len(vars) / nrows

#Set up the Output netCDF file
nc = pycdf.CDF(infilename, NC.NOWRITE)
nc.automode()

subplotNum = 0
fig = plt.figure()

for varname in vars:
    isLogarithmic = (varname in logVars)
    
    var  = nc.var(varname)

    #Determine whether this variable is 3D.  If so, take just the surface values
    #(TODO: Add cmd option to specify the slice?)
    is3d = var.inq_ndims() == 4

    if is3d:
        data = var[0,0,:]
    else:
        data = var[0,:]


    #Determine whether this variable is on a staggered grid.
    #If so, pad with zeroes (TODO: Average properly!)
    isStagger = 'x0' in [nc.dim(d).inq_name() for d in var.inq_dimid()]
    if isStagger:
        data2 = data
        data = np.zeros([i+1 for i in data2.shape])
        data[:-1, :-1] = data2


    if invert:
        data = np.flipud(data)
    #Compute the minimum value, maximum value, and contour intervals
    #for this data
    minval = np.min(data)
    maxval = np.max(data)
    print varname, np.min(data), np.max(data)

    subplotNum += 1
    plt.subplot(nrows, ncols, subplotNum)

    if isLogarithmic:
        minval = np.log(minval+1)/np.log(2)
        maxval = np.log(2*maxval)/np.log(2)
    
    contourInterval = (maxval-minval)/14
    contours = np.arange(minval, maxval, contourInterval) 
    
    if isLogarithmic:
        contours = 2**contours
        norm = colors.LogNorm()
    else:
        norm = None
    print contours

    #contours = [1,4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

    y = np.arange(data.shape[0])
    x = np.arange(data.shape[1])

    cs = plt.contourf(x,y,data,contours,
                    cmap=plt.cm.gist_ncar,
                    norm=norm,
                    alpha=0.75,
                    antialiasing=True)

    cb = colorbar(ticks=contours)
    cb.ax.set_yticklabels([str(int(x)) for x in contours])

    title = varname
    if isLogarithmic:
        title = "Log " + title

    plt.title(title) # add a title

plt.show()

nc.close()
