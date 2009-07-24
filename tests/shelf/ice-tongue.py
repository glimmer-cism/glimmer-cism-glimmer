#!/usr/bin/python

import glimcdf
import numpy
import sys
import math
import ConfigParser

#Python file to set up NetCDF input for computation of steady-state velocities
#in Van Der Veen's ice tongue analytical sln.

#Non-dimensional variable scaling parameters
#for use in the steady-state thickness equations
Z = 205.7426
U = 400.0
L = 274.32
n    = 3          #Power dependence of flow law
rhoi = 910.0      #Density of ice, kg/m^3
rhow = 1028.0     #Density of seawater, kg/m^3
g    = 9.8        #Grav accel. m/s^2
A    = 1.15e-17   #Glen's flow law exponent

#Steady-state thickness profile
def steady_state_thick(x, h0, u0, n):
    #x  - x coordinate of the point to compute for, with x=0 at grounding line.
    #h0 - Thickness at grounding line
    #u0 - Velocity at grounding line
    #n  - Flow law exponent (usually 3)

    #All args should be real numbers
    x  = float(x)
    h0 = float(h0)
    u0 = float(u0)
    n  = float(n)

    #Convert all arguments to non-dimensional forms
    h0 /= Z
    x  /= L
    u0 /= U

    #Get the ice flux at grounding line
    q0 = u0*h0

    #Evaluate the exact solution
    num   = q0**(n+1.0) * (h0**(-(n+1.0)) - 1.0)
    denom = (q0 + x)**(n+1.0)
    h     = (1.0 + num / denom)**(-1.0 / (n + 1.0))

    #Convert to non-dimensional values
    #Numbers copied from MacAyeal Lessons
    h *= Z

    return h

#Computes the analytical velocity for point n given the thicknesses H_n and H_{n-1}, and the velocity v_{n-1}.
#Uses trapezoidal integration.
#This is written to the NetCDF file to potentially use as an initial guess.
def analyticalVel(v_nm1, h_nm1, h_n, A, delta):
    def strainRate(h):
        return A * (.25 * rhoi * g * h * (1-rhoi/rhow) )**n

    return v_nm1 + .5 * (strainRate(h_nm1) + strainRate(h_n))*delta

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(sys.argv[1])


#Set arrhenius parameter in the config file
parser = ConfigParser.ConfigParser()
parser.read(sys.argv[1])
parser.set("parameters", "default_flwa", str(A))  
f=open(sys.argv[1], "w")
parser.write(f)
f.close()

#Create variables for topography and ice thickness
topg = glimcdf.setup_variable(nc, "topg")
thk  = glimcdf.setup_variable(nc, "thk")
beta = glimcdf.setup_variable(nc, "beta", staggered=True)
kinbcu = glimcdf.setup_variable(nc, "uvelbc",  staggered=True, useZ=True)
kinbcv = glimcdf.setup_variable(nc, "vvelbc",  staggered=True, useZ=True)
guessu = glimcdf.setup_variable(nc, "uvelhom", staggered=True, useZ=True)
guessv = glimcdf.setup_variable(nc, "vvelhom", staggered=True, useZ=True)

#Determine the total length of the domain
L = shape.nx*shape.dx

rho_i = 910.0
rho_w = 1028.0

q_0 = 4e5 #Ice flux at the grounding line, cannot be zero!!
h_0 = 1000
u_0 = q_0 / h_0
h_last = float("NaN")
v_last = float("NaN")

for i in range(shape.ny):
    y = i * shape.dy

    if i == 0:
        thickness = h_0
    elif i < shape.ny - 2:
        thickness = steady_state_thick(y, h_0, u_0, 3)
    else: #Put a spot of no ice at this edge of the domain
        thickness = 0
    for j in range(shape.nx):
        loc = (i,j)
        loc2d = (0,)+loc
        
        #Set the thickness
        thk.put_1(loc2d, thickness)

        #Make sure that the ice is never touching
        #CISM will enforce hydrostatic equilibrium
        topg.put_1(loc2d, -h_0)
        mybeta = 0.0
        if i == 0:
            kinematic_u = 0
            kinematic_v = u_0
            vel_u = 0
            vel_v = u_0
        else:
            kinematic_u = float("NaN")
            kinematic_v = float("NaN") 
            vel_u = 0
            vel_v = analyticalVel(v_last, h_last, thickness, A, shape.dy)

        #Beta and kinematic b.c. are on the staggered grid and should
        #only be placed if we're not on the last dimensios
        
        if (loc[1] < shape.nx - 1 and loc[0] < shape.ny - 1):
            beta.put_1(loc2d, mybeta)
            #Place the kinematic b.c. (This needs to be a 3D field)
            for k in range (shape.nz):
                loc3d = (0,k) + loc
                kinbcu.put_1(loc3d, kinematic_u)
                kinbcv.put_1(loc3d, kinematic_v)
                guessu.put_1(loc3d, vel_u)
                guessv.put_1(loc3d, vel_v)
    h_last = thickness
    v_last = vel_v

nc.close()
