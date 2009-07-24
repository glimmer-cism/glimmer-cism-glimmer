#!/usr/bin/python

import glimcdf
import numpy
import sys
import math
from getopt import gnu_getopt

#Get command-line arguments that specify specifics of the experiment
optlist, args = gnu_getopt(sys.argv, '', ["smooth-beta", "dirichlet-center", "sloped"])

optdict = dict(optlist)

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(args[1])

#Create variables for topography and ice thickness
topg = glimcdf.setup_variable(nc, "topg")
thk = glimcdf.setup_variable(nc, "thk")
beta = glimcdf.setup_variable(nc, "beta", staggered=True)

if "--dirichlet-center" in optdict:
    uvelbc = glimcdf.setup_variable(nc, "uvelbc", staggered=True, useZ=True)
    vvelbc = glimcdf.setup_variable(nc, "vvelbc", staggered=True, useZ=True)

#Determine the total length of the domain
L = shape.nx*shape.dx

for i in range(shape.nx):
    x = shape.dx * (i-shape.nx/2)
    xhat = float(i)/(shape.nx-1)
    for j in range(shape.ny):
        y = shape.dy * (j-shape.ny/2)
        yhat = float(j)/(shape.ny-1)
        topg.put_1((0,j,i), -2000) #Everything at sea level
        
        #compute dist. from center
        r = math.sqrt((xhat - .5)**2 + (yhat - .5)**2) 
        
        #Kilometer-thick ice in a circle
        if r < .48:
            #Get down with the thickness
            thickness = 1000
            #If "--sloped" was passed as a command line arg, create a sloped shelf.
            #Otherwise, create a 1000-m thick "pancake"
            if "--sloped" in optdict:
                thickness *= 1-r

            thk.put_1((0,j,i), thickness)
        else:
            thk.put_1((0,i,j), 0)


        if i < shape.ny - 1 and j < shape.nx - 1:
            #Put a zero-slip boundary at the very center of the domain
            #This is required to produce a well-posed problem
            #Everywhere else should have 0 traction.
            if "--smooth-beta" in optdict:
                beta.put_1((0,j,i),1+1e10*math.exp(-(x**2 + y**2)/5e5))
            else:
                beta.put_1((0,j,i), 1)
            
    
            if "--dirichlet-center" in optdict:
                for k in range(shape.nz):
                    uvelbc.put_1((0,k,i,j),float("NaN"))
                    vvelbc.put_1((0,k,i,j),float("NaN"))

for i in range(-1,2):
    for j in range(-1,2):
        x = shape.nx/2-i
        y = shape.ny/2-j
        
        if "--smooth-beta" not in optdict:
            beta.put_1((0,x,y),1e10)
        
        if "--dirichlet-center" in optdict:
            for k in range(shape.nz):
                uvelbc.put_1((0,k,x,y),0)
                vvelbc.put_1((0,k,x,y),0)

nc.close()
