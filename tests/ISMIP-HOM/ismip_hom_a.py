#!/usr/bin/env python

#Run with a command line argument that is the config file to read filenames and grid data from

import glimcdf
from math import sin,cos,tan,pi
import numpy
import sys

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(sys.argv[1])

#Create variables for topography and ice thickness
topg = glimcdf.setup_variable(nc, "topg")
thk = glimcdf.setup_variable(nc, "thk")

#Determine the total length of the domain

for i in range(shape.nx):
    #Our surface is a uniform slope.
    z = 4000 - tan(.5*pi/180)*i*shape.dx
    for j in range(shape.ny):
        #Bumpy bed
        thick = 1000 - 500*sin(2*pi*float(i)/(shape.nx - 2)) * sin(2*pi*float(j)/(shape.ny - 2))
        #Write this data point
        topg.put_1((0,j,i), z - thick)
        thk.put_1((0,j,i), thick)

nc.close()
