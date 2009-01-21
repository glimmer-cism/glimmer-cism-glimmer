#!/usr/bin/python
import math
import numpy
import pylab
import glimcdf

#Create a vector of (x,y,z) coordinates the represent the coordinates of each point generated for the 
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
        for j in range(ny/2, ny/2+1):
            x = float(minX + deltaX*i)
            y = float(minY + deltaY*j)
            z = -H0 + a0*( math.exp( -(x**2 + y**2) / (sigma**2) ) )

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
        print collectedLines
    return collectedLines

#Sorts an array of (x,y,z) corrdinates by the x value and returns a new sorted array
def sortCoordinateLine(points):
    return sorted( points, lambda l1, l2:cmp(l1[0], l2[0]) )

def interpolateCoordinateLine(points, startX, deltaX):
    if points[0][0] > startX:
        raise Exception("The starting x is before the first x provided (you are asking for extrapolation, not interpolation)")

    newPoints = []
    x = startX
    maxX = points[len(points) - 1][0]

    xzPlot(points)

    currentOldPoint = 0
    while x < maxX:
        #See if we need to consider the next interval.  This happens if the x value we are interpolating onto
        #is greater than the end point of the current interval
        while points[currentOldPoint + 1][0] < x:
            currentOldPoint += 1
        
        p1 = points[currentOldPoint]
        p2 = points[currentOldPoint + 1]
        
        print p1[0],x,p2[0]
        
        alpha = (x - p1[0])/(p2[0] - p1[0])
        z = p1[2]*alpha + p2[2]*(1-alpha)
        newPoints.append((x, p1[1], z))
        x += deltaX
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
    nx = shape.nx
    ny = shape.ny

    points = createBumpPoints(nx, ny)
    points = rotateBumpPoints(points, 3*math.pi/180) #Rotate points so average slope is 3 degrees

    lineDict = collectYCoordinateLines(points)
    lines = sorted(list(lineDict.items()), lambda l1, l2:cmp(l1[0], l2[0]))

    j = 0

    output = numpy.zeros((nx, ny))

    for y, points in lines:
        points = sortCoordinateLine(points)
        points = interpolateCoordinateLine(points, 0, 1000*100/nx)
        for i, point in enumerate(points):
            output[i,j] = point[2]

        j += 1

    for i in range(nx):
        for j in range(ny):
            topg.put_1((0,i,j), output[i,j])

    print output
