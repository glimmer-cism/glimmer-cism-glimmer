#!/usr/bin/python

import glimcdf
import numpy
import sys
import math
import ConfigParser
from pycdf import NC

#Check whether the correct command line option for a config file was specified
if len(sys.argv) == 1:  
    print "ERROR: You must specify the configuration file to read from."
    sys.exit(1)

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(sys.argv[1])

#Grab from the config file whether or not we are running with periodic boundary conditions
#in the x direction.  If we are, then we will ignore the kinematic boundary condition
#along the vertical edges of the ice shelf (this allows us to simulate a simple 1D
#solution.
parser = ConfigParser.ConfigParser()
parser.read(sys.argv[1])
periodic_ew = int(parser.get("options","periodic_ew"))

#Set the flow law for this experiment in the configuration file and write the file back out
parser.set("parameters", "default_flwa", "4.6e-18")  
f=open(sys.argv[1], "w")
parser.write(f)
f.close()

#Create variables for topography and ice thickness
topg = glimcdf.setup_variable(nc, "topg")
thk = glimcdf.setup_variable(nc, "thk")
uvelbc = glimcdf.setup_variable(nc, "uvelhom", staggered=True, useZ=True)
vvelbc = glimcdf.setup_variable(nc, "vvelhom", staggered=True, useZ=True)
beta = glimcdf.setup_variable(nc, 'beta', staggered=True)
kinbcmask = glimcdf.setup_variable(nc, 'kinbcmask', staggered=True, type=NC.INT)

#Determine the total length of the domain
L = shape.nx*shape.dx

rho_i = 910.0
rho_w = 1028.0

for i in range(shape.nx):
    xhat = float(i)/(shape.nx-1)
    for j in range(shape.ny):
        loc = (0,j,i)
        yhat = float(j)/(shape.ny-1)
        topg.put_1(loc, -2000) #Everything at sea level
        
        if j < 2:
            thk.put_1(loc, 0)
        else:
            thk.put_1(loc,1000)

        if i < shape.nx - 1 and j < shape.ny - 1:
            beta.put_1(loc, 0)

uvelbc[:] = 0
vvelbc[:] = 0
kinbcmask[:] = 0
kinbcmask[:,shape.ny-2,:] = 1

nc.close()
