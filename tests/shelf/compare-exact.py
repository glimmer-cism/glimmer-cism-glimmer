#!/usr/bin/python
import glimcdf
import pycdf
from pycdf import NC
from getopt import gnu_getopt
import sys
import matplotlib.pyplot as plt
import math
import numpy
import ConfigParser

#This script compares an ice shelf solution to the Weertman analytical 1D solution

#Given a set of levels (such that levels[i+1] > levels[i]),
#returns an array containing the spacing between the leels.
def levelsToSpacings(levels):
    spacings = []
    for i in range(len(levels) - 1):
        spacings.append(levels[i+1] - levels[i])
    return spacings

#Integrate a function given a potentially uneven spacing
#Note that I am not using scipy.integrate because those
#routines require a function that is defined at all values
#of f, not a discretized function.
def integrate(values, spacings, initialValue = 0):
    integratedValues = [initialValue]
    tot = initialValue
    for i in range(len(values) - 1):
        tot += (values[i] + values[i+1])/2 * spacings[i]
        integratedValues.append(tot)
    return integratedValues

#Returns start, stop
def getArrayBounds(thcks):
    firstIndex = -1
    lastIndex = -1

    for i, H in enumerate(thcks):
        if H != 0:
            lastIndex = i
            if firstIndex == -1:
                firstIndex = i
    return firstIndex, lastIndex
        
#Vertically averages a value down the ice column
#Note that, because of the rescaled vertical coordinate, the "height" is 1
#(hence the lack of division in our averaging!)
def verticallyAverageColumn(f, levels):
    v =  integrate(f, levelsToSpacings(levels))
    return v[len(f) - 1]

#Model Parameters
n    = 3          #Power dependence of flow law
rhoi = 910.0      #Density of ice, kg/m^3
rhow = 1028.0     #Density of seawater, kg/m^3
g    = 9.8        #Grav accel. m/s^2

#Parse command line arguments
optlist, args = gnu_getopt(sys.argv, '', ["compare-mode=", "exp=", "tstep=", "ignore="])

#Load an options dictionary with defaults
optdict = {}
optdict["--compare-mode"] = "averaged"
optdict["--exp"] = "confined"
optdict["--tstep"] = 0
optdict["--ignore"] = 0
if len(args) == 1:
    print "ERROR: Must specify a Glimmer/CISM configuration file!"
    sys.exit(1)

#Read the configuration file.
parser = ConfigParser.ConfigParser()
parser.read(args[1])
#Get the value of the arrhenius relationship as specified in the config file.
A = float(parser.get("parameters", "default_flwa"))
#Get the name of the output file
inputFileName = parser.get("CF output", "name")


#Override defaults with anything specified on the command line
for arg, val in optlist:
    optdict[arg] = val
time = int(optdict["--tstep"])
#Open the input file.  This should be the output of a circular shelf run
try:
    inputFile = pycdf.CDF(inputFileName)
except:
    print "Could not open the input file", inputFileName
    print "Make sure that you ran Glimmer/CISM using the configuration file you specified before running this script."
    sys.exit(1)

#Get the variables from the input file that we will use.
uvel = inputFile.var("uvelhom")
vvel = inputFile.var("vvelhom")
thck = inputFile.var("thk")

#Get the number of grid points and the grid spacing from the NetCDF file
nx = inputFile.dim("x0").inq_len()
ny = inputFile.dim("y0").inq_len()
nz = inputFile.dim("level").inq_len()
dx = inputFile.var("x0")[1] - inputFile.var("x0")[0]
dy = inputFile.var("y0")[1] - inputFile.var("y0")[0]
levels = inputFile.var("level").get().squeeze()

#Based on the type of experiment being conducted, set the starting point and the direction for
#the flowline
if optdict["--exp"] == "circle90": #Take from a circular shelf solution parallel to the grid spacing
    startx = nx / 2
    starty = ny / 2
    deltax = 1
    deltay = 0
    spacing = dx
elif optdict["--exp"] == "circle45": #Take from a circular shelf solution perpendicular to the grid spacing
    startx = nx / 2
    starty = ny / 2
    deltax = 1
    deltay = 1
    spacing = math.sqrt(dx**2 + dy**2)
elif optdict["--exp"] == "confined": #Take from the confined shelf solution
    starty = ny - 1 - int(optdict["--ignore"])
    startx = nx / 2
    deltax = 0
    deltay = -1
    spacing = dy
else:
    print "Invalid experiment type specified"
    sys.exit(1)

#Extract a flowline.  Note that we will need to average thickness onto the staggered
#grid to do an intercomparison.
vels = []
thcks = []
x = startx
y = starty

while x >= 0 and x < nx and y >= 0 and y < ny:
    ucolumn = uvel.get(start=(time, 0, y, x), count=(1,nz,1,1)).squeeze() 
    vcolumn = vvel.get(start=(time, 0, y, x), count=(1,nz,1,1)).squeeze() 
     
    #Get the velocity components
    u = verticallyAverageColumn(ucolumn, levels)
    v = verticallyAverageColumn(vcolumn, levels)
 
    #Compute the velocity magnitude
    vel = math.sqrt(u**2 + v**2)
    vels.append(vel)

    thcks.append(thck.get_1((time, y, x)))

    x += deltax
    y += deltay

#Convert thicknesses to numpy array so that we can vectorize operations
thcks = numpy.array(thcks)
analyticalStrainRates = A * (.25 * rhoi * g * thcks * (1-rhoi/rhow))**n

print A
print analyticalStrainRates

#Integrate the analytical solution
analyticalVels = integrate(analyticalStrainRates,[spacing]*len(thcks))
#The above integration assumed an integration constant v_0 = 0.
#This may not be the case if a different velocity was specified as a dirichlet condition
#To handle this, we assume that v0 is the model's computed initial velocity.
v0 = vels[0]
analyticalVels = [v0 + v for v in analyticalVels]

#Chop of the parts of the array that have 0 thickness
i,j = getArrayBounds(thcks)
vels = vels[i:j]
analyticalVels = analyticalVels[i:j]


#Print out a numeric comparison
print "x thickness exact analytical % difference"
for i, (exact, approx) in enumerate(zip(analyticalVels, vels)):
    print i*spacing, thcks[i], exact, approx, 100*(approx-exact)/(exact+.00001)

#Plot the analytical solution and the numeric solution.
xs = numpy.arange(0, len(vels)) * spacing
plt.plot(xs, vels, 'b', label="Model Velocities")
plt.plot(xs, analyticalVels, 'g', label="Analytical Velocities")
plt.legend(loc="upper left")
plt.show()
