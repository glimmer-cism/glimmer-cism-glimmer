#!/usr/bin/python
import math
import numpy
import pylab
import glimcdf

#Create a vector of (x,y,bed,thk) coordinates the represent the coordinates of each point generated for the 
#bed of ismip-hom f.  This is a funky way to do it, but we do it this way because additional
#transformation will be needed before these are ready to use
def createBumpPoints(nx, ny, H0=1000):
    #Implement constants in test requirement 1
    sigma = 10.0*H0
    a0 = 0.1*H0

    #Implement domain size per test requirement 2, domain center per test requirement 1
    domainSize = 100*H0
    deltaX = domainSize / nx
    deltaY = domainSize / ny
    minX = -domainSize / 2
    minY = -domainSize / 2

    points = []

    for i in range(nx):
        for j in range(ny):
            x = float(minX + deltaX*i)
            y = float(minY + deltaY*j)
            z = -H0 + a0*( math.exp( -(x**2 + y**2) / (sigma**2) ) )
            thk = -z
            points.append((x,y,z))

    return points

#Rotates a list of (x,y,z) points by the given angle about the y axis (so that all the slope is in
#the x axis and none in the y)
def rotateBumpPoints(points, theta=.3):

    outputPoints = []
    for point in points:
        x,y,z = point
        #Perform rotation about y axis, this was derived from the product of (x,y,z) with the rotation matrix
        xprime = x*math.cos(theta) - z*math.sin(theta)
        yprime = y
        zprime = x*math.sin(theta) + z*math.cos(theta)
        outputPoints.append((xprime, yprime, zprime))

    return outputPoints

#Collects points into lists where each list has the same y coordinate.  This is doable because the rotation
#only changed the x and z coordinates. We do this because we need to interpolate between x coordinates
#but all on the same line of y coordinates
#Lines are collected into a dictionary 
def collectYCoordinateLines(points):
    collectedLines = {}
    for point in points:
        y = point[1]
        if y not in collectedLines:
            collectedLines[y] = []
        collectedLines[y].append(point)
    return collectedLines

#Sorts an array of (x,y,z) corrdinates by the x value and returns a new sorted array
def sortCoordinateLine(points):
    return sorted( points, lambda l1, l2:cmp(l1[0], l2[0]) )

def interpolateCoordinateLine(points, oldPoints):
    extrapolationSlope = (points[1][2] - points[0][2]) / (points[1][0] - points[0][0])
    initialX = points[0][0]
    initialZ = points[0][2]

    finalX = points[len(points) - 1][0]
    finalZ = points[len(points) - 1][2]

    i = 0
    newPoints = []

    currentInterpIndex = 0

    #Extrapolate with a sloped line up until the first place where we have data.
    for x,y,z in oldPoints:
        #Extrapolate with a sloped line up until the first place where we have data.
        if x < initialX:
            zprime = initialZ - extrapolationSlope*(initialX - x)
            #zprime = 0
        elif x > finalX:
            zprime = finalZ   + extrapolationSlope*(x - finalX)
            #zprime = 0
        else: #Interpolate
            #Because these points are ordered by x coordinate in both lists, we can be sure that
            #the starting point for the interpolation will never be *before* the starting point
            #of the last interpolation.  This cuts our searching down to an average of constant time
            #for the whole run.
            while points[currentInterpIndex + 1][0] < x:
                currentInterpIndex += 1
            #We stopped at the first new point that is after the x value we want.  Therefore, the x value
            #we want falls between the two indices we've identified
            alpha = (x - points[currentInterpIndex][0]) / (points[currentInterpIndex+1][0] - points[currentInterpIndex][0])
            
            zprime = alpha*points[currentInterpIndex+1][2] + (1-alpha)*points[currentInterpIndex][2]
        newPoints.append((x,y,zprime))
    return newPoints

def xzPlot(points):
    xs = []
    zs = []
    for x,y,z in points:
        xs.append(x)
        zs.append(z)
    pylab.show(pylab.plot(xs,zs))

if __name__ == "__main__":

    nc, shape = glimcdf.nc_from_config("ishom.f.config")
    topg = glimcdf.setup_variable(nc, "topg")
    thk  = glimcdf.setup_variable(nc, "thk" )

    nx = shape.nx
    ny = shape.ny

    points = createBumpPoints(nx, ny)

    lineDict = collectYCoordinateLines(points)
    lines = sorted(list(lineDict.items()), lambda l1, l2:cmp(l1[0], l2[0]))

    j = 0

    topg_field = numpy.zeros((nx, ny))
    thck_field = numpy.zeros((nx, ny))
    
    for y, points in lines:
        newpoints = rotateBumpPoints(points,3*math.pi/180)
        newpoints = sortCoordinateLine(newpoints)
        newpoints = interpolateCoordinateLine(newpoints, points)
        
        #Get a point 1000 meters above the first point to represent the (smooth) surface
        #elevation.  We'll keep adding a constant slope to that to give us the surface, which
        #in turn easily gives us the thickness.
        currentSurf = 1000 + newpoints[0][2]
        slope = 1000 * math.tan(3.0*math.pi/180)

        for i, point in enumerate(newpoints):
            topg_field[i,j] = point[2]
            thck_field[i,j] = currentSurf - point[2]
            currentSurf += slope

        j += 1

    for i in range(nx):
        for j in range(ny):
            topg.put_1((0,i,j), topg_field[i,j])
            thk.put_1( (0,i,j), thck_field[i,j])

