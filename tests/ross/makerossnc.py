#!/usr/bin/python
import numpy
import pycdf
from pycdf import NC
import glimcdf
import sys
from getopt import gnu_getopt
import math
rhoi=910.0
rhow=1028.0

rows = 115
cols = 151
shape = (rows, cols)


#Parse command line arguments.
optlist, args = gnu_getopt(sys.argv, "", ["nlevels=", "fake-shelf"])
optdict=dict(optlist)

nlevels = int(optdict["--nlevels"])
fakeShelf = ("--fake-shelf",'') in optlist

#Reads the data in the input file into a dictionary mapping to lists of lists.
#The dictionary keys are the comments that begin each section.
def readData(inputfile):
    currentKey = None
    currentList = []
    d={}
    for line in inputfile:
        line = line.strip()
        if line.startswith("#"):
            if currentKey:
                d[currentKey] = currentList
            currentList = []
            currentKey = line[1:].strip().lower()
        elif line: #Non-empty non-header line, assume it has numeric data
            currentList.append([float(s) for s in line.split()])
    d[currentKey] = currentList
    return d

#This function puts a border around the dataset, since the Pattyn model at least needs
#this to work.
def newfield(inputData, dtype, default=0):
    global shape
    field = numpy.zeros(shape, dtype)
    field[:,:] = default
    field[2:-2,2:-2] = numpy.array(inputData, dtype=dtype)
    return field

inputfile = open("111by147Grid.dat")

#Read the number of rows and columns from the input file
line = inputfile.readline().strip().split()
d = readData(inputfile)
inputfile.close()

print d.keys()

#Create data fields
existency = newfield(d["existency table:"], dtype='i')
ice_vel_azimuth = newfield(d["ice velocity azimuth grid"], dtype='d')
ice_vel_magnitude = newfield(d["ice velocity magnitude"], dtype='d')
thickness = newfield(d["thickness"], dtype='d')
vel_is_reliable = newfield(d["reliable velocity obs"], dtype='i')
seabed_depth = newfield(d["seabed depth"], dtype='d', default=5000.0)
fake_shelf_mask = newfield(d["fake ice shelf region"], dtype='i')
accumulation = newfield(d["surface accumulation"], dtype='d')
bbar = newfield(d["flowlaw"], dtype='d')
surface_temp = newfield(d["surface temperature"], dtype='d')

#Read the kinematic boundary mask file.  This is a set of (i,j) coordinates that
#specifies where the velocity read from the data file needs to be held as a bc
#in the model.
kinematic_bc_mask = numpy.zeros(shape, dtype='i')

#HACK: Get rid or Roosovelt
#existency[22:57, 22:57] = 1

#Create a "land ice mask" that specifies where ice is grounded. 
#These points are masked out by "existency", as this intercomparison considers
#floating ice only.
#We need this to determine whether a boundary condition that is unspecified
#should be treated as a 0-kinematic boundary or as a shelf front.
land_ice_mask = numpy.where(rhoi/rhow*thickness > seabed_depth, 1, 0)

shelf_front_mask = newfield(0, dtype='i')

#Remove parts of the "fake shelf" that don't correspond to a point on the existency mask
fake_shelf_mask = numpy.where(existency == 1, fake_shelf_mask, 0)

#Cheat some fake shelf points
fake_shelf_mask[38, 23] = 0

#Derive the "real shelf" as all points where ice exists that aren't part of the fake shelf
real_shelf_mask = existency - fake_shelf_mask

#If we are going to use the fake shelf in our model, we need to identify as shelf front
#all points that are on the front of the fake shelf mask
#Otherwise, we identify as shelf front all points that are adjacent to a fake shelf point
#Note that these points are one *past* the shelf front!!
if fakeShelf:
    for i in range(cols):
        if existency[2,i] == 1:
            shelf_front_mask[1,i] = 1
else:
    for i in range(1,rows-1):
        for j in range(1,cols-1):
            #If this point is on the fake shelf (meaning real_shelf_mask[i,j] will be 0
            #necessarily) and at least one surrounding point is on the real shelf,
            #then this point is the shelf front.
            if fake_shelf_mask[i,j] and numpy.sum(real_shelf_mask[i-1:i+2,j-1:j+2]) > 0:
                shelf_front_mask[i,j] = 1
    #HACK: the right side of the shelf front doesn't get identified using the above,
    #so I'm just going to hack it in nice and ugly here.
    for i in range(cols):
        if real_shelf_mask[2,i]:
            shelf_front_mask[1,i] = 1


shelf_front_mask[1,101] = 1

inputfile = open("kbc.dat")
for line in inputfile:
    #Read the coordinates of the kinematic bcs from the 
    i,j = line.strip().split()
    i=int(i)
    j=int(j)
    kinematic_bc_mask[i,j] = 1
inputfile.close()

#Read in the inlets.dat file, which specifies additional dirichlet conditions
inputfile = open("inlets.dat")
for line in inputfile:
    i, j, mag, azimuth = line.strip().split()
    kinematic_bc_mask[i,j] = 1
    ice_vel_azimuth[i,j] = azimuth
    ice_vel_magnitude[i,j] = mag



#Create a NetCDF with the raw data fields.  This will be a useful
#debugging tool
nc = pycdf.CDF("ross-raw.nc", NC.WRITE | NC.CREATE | NC.TRUNC)
nc.automode()
glimcdf.setup_dimensions(nc, cols, rows, 1, 6822, 6822)

existency_var = glimcdf.setup_variable(nc, "existency", type=NC.INT)
azimuth_var = glimcdf.setup_variable(nc, "velo_azimuth", type=NC.DOUBLE)
magnitude_var = glimcdf.setup_variable(nc, "velo_mag", type=NC.DOUBLE)
thck_var = glimcdf.setup_variable(nc, "thk", type=NC.DOUBLE)
topg_var = glimcdf.setup_variable(nc, "topg", type=NC.DOUBLE)
reliable_var = glimcdf.setup_variable(nc, "velo_reliable", type=NC.INT)
fake_shelf_var = glimcdf.setup_variable(nc, "fake_shelf_mask", type=NC.INT)
real_shelf_mask_var = glimcdf.setup_variable(nc, "real_shelf_mask", type=NC.INT)
acab_var = glimcdf.setup_variable(nc, "acab", type=NC.DOUBLE)
bbar_var = glimcdf.setup_variable(nc, "bbar", type=NC.DOUBLE)
surface_temp_var = glimcdf.setup_variable(nc, "surface_temp", type=NC.DOUBLE)
kinbc_var = glimcdf.setup_variable(nc, "kinematic_bc_mask", type=NC.INT)
land_ice_mask_var = glimcdf.setup_variable(nc, "land_ice_mask", type=NC.INT)
shelf_front_mask_var = glimcdf.setup_variable(nc, "shelf_front_mask", type=NC.INT)

existency_var[0,:,:] = existency
azimuth_var[0,:,:] = ice_vel_azimuth
magnitude_var[0,:,:] = ice_vel_magnitude
thck_var[0,:,:] = thickness
topg_var[0,:,:] = seabed_depth
fake_shelf_var[0,:,:] = fake_shelf_mask
acab_var[0,:,:] = accumulation
bbar_var[0,:,:] = bbar
surface_temp_var[0,:,:] = surface_temp
kinbc_var[0,:,:] = kinematic_bc_mask
land_ice_mask_var[0,:,:] = land_ice_mask
shelf_front_mask_var[0,:,:] = shelf_front_mask
real_shelf_mask_var[0,:,:] = real_shelf_mask
reliable_var[0,:,:] = vel_is_reliable
nc.close()

#Now, we can actually start processing the data so it's suitable for use
#in CISM.

#Convert velocity from polar specification (azimuth/mag) to rectangular
#specification (u/v)
uvel = numpy.sin(ice_vel_azimuth*math.pi/180) * ice_vel_magnitude
vvel = numpy.cos(ice_vel_azimuth*math.pi/180) * ice_vel_magnitude

#We need to get rid of the "fake shelf" part of the existency mask; CISM's ice shelf
#models don't need it.
if not fakeShelf:
    existency = real_shelf_mask

#Create fields that hold only the velocity components where they need to be
#specified as kinematic boundary conditions.
uvelbc = numpy.where(kinematic_bc_mask == 1, uvel, float("NaN"))
vvelbc = numpy.where(kinematic_bc_mask == 1, vvel, float("NaN"))

existency2 = existency



#Identify points as being next to the ice shelf.  If a shelf front wasn't
#identified at that point, extend the ice shelf onto that point artificialy.
for i in range(1, rows-1):
    for j in range(1, cols-1):
        #This is a point right off of the ice shelf
        if not existency[i,j] and numpy.sum(existency[i-1:i+2,j-1:j+2]) > 0: 
            #This point is on not the shelf front
            if not shelf_front_mask[i,j]:
                #if uvel and vvel are both NaN, set them to 0
                if uvelbc[i,j] != uvelbc[i,j] and vvelbc[i,j] != vvelbc[i,j]:
                    uvelbc[i,j] = 0
                    vvelbc[i,j] = 0
                    uvel[i,j] = 0
                    vvel[i,j] = 0
                #Artificially extend the shelf out to here.
                localexist = existency[i-1:i+2, i-1:i+2]
                localthickness = numpy.where(localexist == 1, thickness[i-1:i+2, i-1:i+2], 0)
                avgthck = numpy.sum(localthickness)/numpy.sum(localexist)
                if avgthck != avgthck:
                    avgthck = 0
                #thickness[i,j] = avgthck
                existency2[i,j] = 1

existency = existency2

#Use the new existency to set the velocity and the thickness to zero.
thickness = numpy.where(existency == 1, thickness, 0) 


#Move the boundary conditions onto the staggered grid, where CISM expects them
#For now I'll just discard the last row and col of data
uvelbc = uvelbc[0:-1, 0:-1]
vvelbc = vvelbc[0:-1, 0:-1]

#Stack the kinematic boundary condition 'nlevels' times so that it is in the 3D format that CISM
#wants
uvelbc = numpy.array([uvelbc]*nlevels)
vvelbc = numpy.array([vvelbc]*nlevels)

#write the netcdf file in a format that CISM can use
nc = pycdf.CDF("ross.nc", NC.WRITE|NC.CREATE|NC.TRUNC)
nc.automode()
glimcdf.setup_dimensions(nc, cols, rows, nlevels, 6822, 6822)

#Thickness
thk_var = glimcdf.setup_variable(nc, "thk", type=NC.DOUBLE)
thk_var[0,:,:] = thickness

#Topography
#This is specified as "depth" in the input rather that "elevation", so
#we need to negate it.
topg_var = glimcdf.setup_variable(nc, "topg", type=NC.DOUBLE)
topg_var[0,:,:] = -seabed_depth
topg_var[0,:,:] = -10000

#Kinematic boundary conditions
uvelbc_var = glimcdf.setup_variable(nc, "uvelbc", type=NC.DOUBLE, staggered=True, useZ=True)
vvelbc_var = glimcdf.setup_variable(nc, "vvelbc", type=NC.DOUBLE, staggered=True, useZ=True)

uvelbc_var[0,:,:,:] = uvelbc
vvelbc_var[0,:,:,:] = vvelbc

#Run everything on 0 basal traction
beta_var = glimcdf.setup_variable(nc, "beta", staggered=True)
beta_var[:,:,:] = 0

#Write the velocities as an initial guess
uvel_var = glimcdf.setup_variable(nc, "uvelhom", type=NC.DOUBLE, staggered=True, useZ=True)
vvel_var = glimcdf.setup_variable(nc, "vvelhom", type=NC.DOUBLE, staggered=True, useZ=True)

uvel[:rows-1, :cols-1] = numpy.where(uvelbc[0,:,:] == uvelbc[0,:,:], uvel[:rows-1, :cols-1], 0)
vvel[:rows-1, :cols-1] = numpy.where(vvelbc[0,:,:] == vvelbc[0,:,:], vvel[:rows-1, :cols-1], 0)

uvel_var[0,:,:,:] = numpy.array([uvel[:rows-1, :cols-1]]*nlevels)
vvel_var[0,:,:,:] = numpy.array([vvel[:rows-1, :cols-1]]*nlevels)

kinbcmask_var = glimcdf.setup_variable(nc, "kinbcmask", type=NC.INT)
kinbcmask_var[:] = 1
kinbcmask_var[0, :rows-1, :cols-1] = numpy.where(uvelbc[0,:,:] == uvelbc[0,:,:], 1, 0)
