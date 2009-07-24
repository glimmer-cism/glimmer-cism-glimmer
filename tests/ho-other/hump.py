#!/usr/bin/python

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

#Set this on a flat surface
topg[:] = 0

#Create a hump of ice in the middle third of the domain, with a maximum heigt of 1000 m.
cenx = shape.nx / 2
ceny = shape.ny / 2

i = numpy.indices((shape.nx, shape.ny)).astype("float32")

#Create a field with the normalized distance from the center of the domain
r_squared = 1.0/8.0 - ((i[1,:] - cenx)/shape.nx)**2 - ((i[0,:] - ceny)/shape.ny)**2 

H = (2000.0 * numpy.sqrt( numpy.where(r_squared > 0, r_squared, 0)))

thk[0,:,:] = H

nc.close()
