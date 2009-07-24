#!/usr/bin/python

import glimcdf
import numpy
import sys
import math

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(sys.argv[1])

#Create variables for topography and ice thickness
topg = glimcdf.setup_variable(nc, "topg")
thk = glimcdf.setup_variable(nc, "thk")
beta = glimcdf.setup_variable(nc, "beta", staggered=True)
uvelbc = glimcdf.setup_variable(nc, "uvelbc", staggered=True, useZ=True)
vvelbc = glimcdf.setup_variable(nc, "vvelbc", staggered=True, useZ=True)
#Determine the total length of the domain
L = shape.nx*shape.dx

for i in range(shape.nx):
    xhat = float(i)/(shape.nx-1)
    for j in range(shape.ny):
        yhat = float(j)/(shape.ny-1)
        topg.put_1((0,j,i), -2000) #Everything at sea level
        
        #Kilometer-thick ice in a circle
        if (i > 2 and i < shape.nx-2 and j > 2 and j < shape.ny-2):
            thk.put_1((0,j,i), 1000)
        else:
            thk.put_1((0,i,j), 0)


        if i < shape.ny - 1 and j < shape.nx - 1:
            #Put a zero-slip boundary at the very center of the domain
            #This is required to produce a well-posed problem
            #Everywhere else should have 0 traction.
            beta.put_1((0,j,i),1)
            for k in range(shape.nz):
                uvelbc.put_1((0,k,i,j),float("NaN"))
                vvelbc.put_1((0,k,i,j),float("NaN"))

for i in range(-2,3):
    for j in range(-2,3):
        x = shape.nx/2-i
        y = shape.ny/2-j
        #beta.put_1((0,shape.nx/2-i,shape.ny/2-j),1e10)
        beta.put_1((0,x,y),float("NaN"))
        for k in range(shape.nz):
            uvelbc.put_1((0,k,x,y),0)
            vvelbc.put_1((0,k,x,y),0)
nc.close()
