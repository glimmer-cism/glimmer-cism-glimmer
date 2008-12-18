import pycdf
from pycdf import NC
from ConfigParser import ConfigParser

class Shape:
    def __init__(self, nx, ny, nz, dx, dy):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy

def init_dimension_var(var,min,step,n=None):
    if n == None:
        n = var.shape()[0]

    print var.inq_name(),var.shape()
    current = min
    for i in range(n):
        var[i] = current
        current += step

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
    timedim = nc.def_dim("time", 1)
    leveldim = nc.def_dim("level", nz)
    #Setup staggered coordinates
    x0dim = nc.def_dim("x0",nx-1)
    y0dim = nc.def_dim("y0",ny-1)
    #Setup unstaggered coordinates
    x1dim = nc.def_dim("x1",nx)
    y1dim = nc.def_dim("y1",ny)
    
    #Set up the variables that determine the values of each
    #dimension coordinate
    timevar = nc.def_var("time",NC.FLOAT,timedim)
    x0var = nc.def_var("x0",NC.FLOAT,x0dim)
    y0var = nc.def_var("y0",NC.FLOAT,y0dim)
    x1var = nc.def_var("x1",NC.FLOAT,x1dim)
    y1var = nc.def_var("y1",NC.FLOAT,y1dim)
    
    init_dimension_var(timevar,0,1,1)
    init_dimension_var(x0var, deltax/2, deltax)
    init_dimension_var(y0var, deltay/2, deltay)
    init_dimension_var(x1var, 0, deltax)
    init_dimension_var(y1var, 0, deltay)

def set_attribute(nc, name, text):
    attr = nc.attr(name)
    attr.put(NC.CHAR,text)

def setup_variable(nc, name, staggered=False):
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
    nx = xvar.shape()[0]
    ny = yvar.shape()[0]
    maxx = xvar[nx - 1]
    maxy = yvar[ny - 1]
    
    var = nc.def_var(name, NC.FLOAT, (nc.dim("time"), ydim, xdim))
    return var
