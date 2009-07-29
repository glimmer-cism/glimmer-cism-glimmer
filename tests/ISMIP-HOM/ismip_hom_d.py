#!/usr/bin/python

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
for j in range(shape.nx):
    #Our surface is a uniform slope.
    x = L*float(j)/(shape.nx - 1)
    z = 2000 - x*tan(.1*pi/180)
    b = 1000 + 1000*sin(2*pi*float(j)/(shape.nx - 1))
    for i in range(shape.ny):
        #Write this data point
        topg.put_1((0,i,j), z - 1000)
        thk.put_1((0,i,j), 1000)
        if j < shape.nx - 1 and i < shape.ny - 1:
           beta.put_1((0,i,j),b)

nc.close()
