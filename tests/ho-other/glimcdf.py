import pycdf
import numpy
from pycdf import NC
import glimcdf
#import scipy.interpolate as interpolate
import scipy.ndimage as sp
import scipy

from ConfigParser import ConfigParser

class Shape:
    def __init__(self, nx, ny, nz, dx, dy):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.xlen = nx*dx
        self.ylen = ny*dy


def dataStream(openfile):
    for line in openfile:
        if line[0] == '#':
            continue
        for token in line.strip().split():
            yield float(token)

def readData(file,header = True):
    print "reading data from ", file.filename
    openfile = open(file.filename)
    if(header):
        for line in openfile:
            tokens = []
            for token in line.strip().split():
                tokens.append(token)   
            if tokens[0] == '#ncols':
                file.ncols = int(tokens[2])
            if tokens[0] == '#nrows':
                file.nrows = int(tokens[2])
            if tokens[0] == '#xllcorner':
                file.xllcorner = float(tokens[2])
            if tokens[0] == '#yllcorner':
                file.yllcorner = float(tokens[2])
            if tokens[0] == '#gridspacing':
                file.gridSpacing = float(tokens[2])
            if tokens[0] == '#nodata':
                file.nodata = float(tokens[2])
            if tokens[0] == '#source':
                file.source = tokens[2 :]
            if tokens[0][1] == '=':
                break
        
    w = file.ncols
    h = file.nrows        
    stream = dataStream(openfile)
    data = numpy.zeros((w, h))
    for j in range(h):
        for i in range(w):
            data[i,j] = stream.next()
    openfile.close()
    #Flip data
    #data = data[:,::-1]
    # Make default type float32 and transpose it
    data=numpy.asarray(data.T,dtype='float32')
    #data=numpy.asarray(data,dtype='float32')
    return data


    
class file_reader():
    ncols = 0
    nrows = 0
    xllcorner = 0
    yllcorner = 0
    gridSpacing = 0
    nodata = -9999
    def __init__(self, file, w = 0, h = 0, header = True):
        self.filename = file
        if w != 0:
            self.ncols = w
        if h != 0:
            self.nrows = h
        self.header = header
        self.data = glimcdf.readData(self, self.header)
        
#Returns an empty NetCDF file from the cf_output
def nc_from_config(configFilename):
    parser = ConfigParser()
    parser.read(configFilename)

    nx = int(parser.get("grid","ewn"))
    ny = int(parser.get("grid","nsn"))
    nz = int(parser.get("grid","upn"))
    
    dx = float(parser.get("grid","dew"))
    dy = float(parser.get("grid","dns"))
    fname = parser.get("CF input", "name")
    
    print "Writing to", fname

    nc = pycdf.CDF(fname, NC.WRITE | NC.CREATE | NC.TRUNC)
    nc.automode()
    setup_dimensions(nc, nx, ny, nz, dx, dy)
    shape = Shape(nx, ny, nz, dx, dy)
    return nc, shape


def setup_dimensions(nc, nx, ny, nz, deltax, deltay):

    xgrid = numpy.vstack( [numpy.arange(nx)*deltax] * ny ).transpose()
    ygrid = numpy.vstack( [numpy.arange(ny)*deltay] * nx )

    setup_dimensions_gridded(nc, xgrid, ygrid, nlevels=nz)

def setup_dimensions_gridded(nc,xgrid,ygrid, nlevels=1, temps = None, sealevel = None, oisotopes = None):
    s = xgrid.shape
    nx = s[0]
    ny = s[1]
    timedim = nc.def_dim("time", 1)
    leveldim = nc.def_dim("level", nlevels)
    #Setup staggered coordinates
    x0dim = nc.def_dim("x0",nx-1)
    y0dim = nc.def_dim("y0",ny-1)
    #Setup unstaggered coordinates
    x1dim = nc.def_dim("x1",nx)
    y1dim = nc.def_dim("y1",ny)
    
    #Dimensions for the temperature, sealevel and oisotopes series, set up the variables and the values
    if temps != None:
        temptimes = nc.def_dim("temptimes",len(temps))
        tempvar = nc.def_var("temptimes",NC.FLOAT, temptimes)
        glimcdf.set_variable_attribute(nc,tempvar, "long_name","Temperature_times")
        glimcdf.set_variable_attribute(nc,tempvar, "units","years")
        tempvar[:] = numpy.asarray(temps, dtype = "float32")
    
    if sealevel != None:
        sealeveltimes = nc.def_dim("sealeveltimes",len(sealevel))
        sealevelvar = nc.def_var("sealeveltimes",NC.FLOAT, sealeveltimes)
        glimcdf.set_variable_attribute(nc,sealevelvar, "long_name","Sealevel_times")
        glimcdf.set_variable_attribute(nc,sealevelvar, "units","years")
        sealevelvar[:] = numpy.asarray(sealevel, dtype = "float32")
    
    if oisotopes != None:
        oisotopestimes = nc.def_dim("oisotopestimes",len(oisotopes))
        oisotopesvar = nc.def_var("oisotopestimes",NC.FLOAT, oisotopestimes)
        glimcdf.set_variable_attribute(nc,oisotopesvar, "long_name", "Oisotope_times")
        glimcdf.set_variable_attribute(nc,oisotopesvar, "units","years")
        oisotopesvar[:] = numpy.asarray(oisotopes, dtype = "float32")

    timevar = nc.def_var("time",NC.FLOAT,timedim)
    glimcdf.set_variable_attribute(nc,timevar, "long_name","Model time")
    glimcdf.set_variable_attribute(nc,timevar, "standard_name","time")
    glimcdf.set_variable_attribute(nc,timevar, "units","year since 1-1-1 0:0:0")
    glimcdf.set_variable_attribute(nc,timevar, "calender","none")
    x0var = nc.def_var("x0",NC.FLOAT,x0dim)
    glimcdf.set_variable_attribute(nc,x0var, "long_name","Cartisian x-coordinate, velocity grid")
    glimcdf.set_variable_attribute(nc,x0var, "standard_name","projection_x_coordinate")
    glimcdf.set_variable_attribute(nc,x0var, "units","meter")
    y0var = nc.def_var("y0",NC.FLOAT,y0dim)
    glimcdf.set_variable_attribute(nc,y0var, "long_name","Cartisian y-coordinate, velocity grid")
    glimcdf.set_variable_attribute(nc,y0var, "standard_name","projection_y_coordinate")
    glimcdf.set_variable_attribute(nc,y0var, "units","meter")
    
    x1var = nc.def_var("x1",NC.FLOAT,x1dim)
    glimcdf.set_variable_attribute(nc,x1var, "long_name","Cartisian y-coordinate")
    glimcdf.set_variable_attribute(nc,x1var, "standard_name","projection_x_coordinate")
    glimcdf.set_variable_attribute(nc,x1var, "units","meter")
    y1var = nc.def_var("y1",NC.FLOAT,y1dim)
    glimcdf.set_variable_attribute(nc,y1var, "long_name","Cartisian y-coordinate")
    glimcdf.set_variable_attribute(nc,y1var, "standard_name","projection_y_coordinate")
    glimcdf.set_variable_attribute(nc,y1var, "units","meter")
    
    # Make the assignments directly
    timevar[0] = 0.

    # Note type coersion
    x0var[:] = numpy.asarray((xgrid[0:-1,0] + xgrid[1:,0]) / 2.,dtype='float32')
    y0var[:] = numpy.asarray((ygrid[0,0:-1] + ygrid[0,1:]) / 2.,dtype='float32')

    x1var[:] = numpy.asarray(xgrid[:,0],dtype='float32')
    y1var[:] = numpy.asarray(ygrid[0,:],dtype='float32')

def set_global_attribute(nc, name, text):
    attr = nc.attr(name)
    attr.put(NC.CHAR,text)

def set_variable_attribute(nc,variableInstance, name, text, type = NC.CHAR):
    attr = variableInstance.attr(name)
    attr.put(type,text)

def setup_variable(nc, name, type = NC.FLOAT, f=None, staggered=False, other=None, useZ=None):
    if other != None:
        if other == "temp":
            return nc.def_var(name, type, nc.dim("temptimes"))
        if other == "sealevel":
            return nc.def_var(name, type, nc.dim("sealeveltimes"))
        if other == "oisotopes":
            return nc.def_var(name, type, nc.dim("oisotopestimes"))
    
    if staggered:
        xvar = nc.var("x0")
        yvar = nc.var("y0")
        xdim = nc.dim("x0")
        ydim = nc.dim("y0")
    else:
        xvar = nc.var("x1")
        yvar = nc.var("y1")
        xdim = nc.dim("x1")
        ydim = nc.dim("y1")
    leveldim = nc.dim("level")
    timedim = nc.dim("time")

    dims = (timedim,)
    if useZ:
        dims = dims + (leveldim,)
    dims = dims + (ydim, xdim)
    
    return nc.def_var(name, type, dims) 

def dataStream(openfile):
    for line in openfile:
        for token in line.strip().split():
            yield float(token)

def readData(filename, w, h):
    print "reading data from ", filename
    openfile = open(filename)
    stream = dataStream(openfile)
    data = numpy.zeros((w, h))
    for j in range(h):
        for i in range(w):
            data[i,j] = stream.next()
    openfile.close()
    # Flip data
    data = data[:,::-1]
    # Make default type float32 and transpose it
    data=numpy.asarray(data.T,dtype='float32')
    return data

def regridDataSimple(data, oldw, oldh, regridFactor):
    neww = int(oldw / regridFactor)
    newh = int(oldh / regridFactor)
    newData = numpy.zeros((neww,newh))
    for i in range(neww):
        for j in range(newh):
            newData[i,j] = data[i*regridFactor, j*regridFactor]
    newData = numpy.asarray(newData, dtype='float32')
    return newData

def regridData(data,xold,yold,xnew,ynew,xlen,ylen,oldGridSpacing,regridFactor = 1.0,cval = 0.0):

    
    #This assumes regularly gridded data input and outputs
    firstCornerX = (xnew - xold)/oldGridSpacing
    firstCornerY = (ynew - yold)/oldGridSpacing
    
    print firstCornerX
    print firstCornerY
    
    [x,y] = numpy.mgrid[firstCornerY:(ylen + firstCornerY):regridFactor, firstCornerX:(xlen + firstCornerX):regridFactor]
   
    newData = sp.map_coordinates(data, [x,y], order=3, cval=cval, prefilter = False)
    
    return newData


#Copies dimensions from one NetCDF file to another NetCDF file.
#ExceptionsCallback is a function that takes as its input the dimension name and length
#and returns as its output a new name, length tuple, or none if the dimension should not
#be copied (if the dimension is unchanged, this should just be a passthrough)
def copyDimensions(inCDF, outCDF, exceptionsCallback=None):
    for dimname, dimlen in inCDF.dimensions().iteritems():
        if exceptionsCallback:
            tup = exceptionsCallback(dimname, dimlen)
            if not tup:
                continue
            else:
                dimname, dimlen = tup
        outCDF.def_dim(dimname, dimlen)


#Copies variables from one NetCDF file to another NetCDF file.
#ExceptionsCallback is a function that takes as its input the following info about the
#variables:
#variable_name, dimension name tuple, variable shape tuple, variable type, and variable data
#It should return a tuple of the above information in that order, or "None" if the variable
#is not to be copied.
#When a variable is called with the callback, the dimensions will have been rotated such that
#x, y are the initial dimensions
def copyVariables(inCDF, outCDF, exceptionsCallback=None):
    for varname, (dimnames, shape, type, index) in inCDF.variables().iteritems():
        print varname
        data = inCDF.var(varname)[:]

        nrotates = 0
        firstTwoSwapped = False
        if ('x0' in dimnames or 'x1' in dimnames) and ('y0' in dimnames or 'y1' in dimnames):
            #Rotate the variable so that either x or y (whichever comes first) is the first dimension
            while dimnames[0][0] != "x" and dimnames[0][0] != "y":
                #Rotate the first dimension to the end
                dimnames = dimnames[1:] + dimnames[:1]
                shape = shape[1:] + shape[:1]
                data = numpy.rollaxis(data, 0, len(data.shape))
                nrotates += 1
            #If y comes first, swap the first two dimensions so that x is guaranteed to come first
            if dimnames[0][0] == "y":
                dimnames = dimnames[1:2] + dimnames[:1] + dimnames[2:]
                shape = shape[1:2] + shape[:1] + shape[2:]
                data = numpy.rollaxis(data, 1, 0)
                firstTwoSwapped = True


        if exceptionsCallback:
            tup = exceptionsCallback(varname, dimnames, shape, type, data)
            if not tup:
                continue
            else:
                varname, dimnames, shape, type, data = tup
        
        #Undo swapping if needed
        if firstTwoSwapped:
            dimnames = dimnames[1:2] + dimnames[:1] + dimnames[2:]
            shape = shape[1:2] + shape[:1] + shape[2:]
            data = numpy.rollaxis(data, 1, 0)
             
        for i in range(nrotates):
            dimnames = dimnames[-1:] + dimnames[:-1]
            shape = shape[-1:] + shape[:-1]
            data = numpy.rollaxis(data, len(data.shape) - 1, 0)

        #Create the new version of the variable in the output file
        outvar = outCDF.def_var(varname, type, dimnames)

        #Assign the new data to this variable
        print outvar.shape()
        print data.shape
        
        outvar[:] = data
       
        
