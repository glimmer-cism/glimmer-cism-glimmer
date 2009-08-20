#!/usr/bin/env python

import glimcdf
from math import sin,cos,tan,pi
import numpy
import sys

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(sys.argv[1])

#Create variables for topography and ice thickness
topg = glimcdf.setup_variable(nc, "topg")
thk = glimcdf.setup_variable(nc, "thk")
beta = glimcdf.setup_variable(nc, "beta", staggered=True)
#Determine the total length of the domain
L = shape.nx*shape.dx

for i in range(shape.nx):
    #Our surface is a uniform slope.
    z = 2000 - tan(.1*pi/180)*i*shape.dx
    for j in range(shape.ny):
        b = 1000 + 1000*sin(2*pi*float(i)/(shape.nx - 2)) * sin(2*pi*float(j)/(shape.ny - 2))
        
        #Write this data point
        topg.put_1((0,j,i), z - 1000)
        thk.put_1((0,j,i), 1000)
        if i < shape.ny - 1 and j < shape.nx - 1:
           beta.put_1((0,j,i),b)

nc.close()
