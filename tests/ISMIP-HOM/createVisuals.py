import os
import os.path
import math
from matplotlib import pyplot

#Lists of model classifications
fullStokes = ["aas1", "aas2", "cma1", "fpa2", "ghg1", "jvj1", "mmr1", "oga1", "rhi1", "rhi3", "spr1", "ssu1", "yko1"]

def getDataFiles(root, experiment, domainSizeKm, isFullStokes):    
    fnameSuffix = experiment + "%03d"%domainSizeKm +".txt"
    files = []
    for root, ds, fs in os.walk(root):
        for file in fs:
            if file.endswith(fnameSuffix):
                model = file[:4]
                if isFullStokes and model in fullStokes or (not isFullStokes) and model not in fullStokes:
                    files.append(os.path.join(root, file))
    return files


class Flowline(object):
    def __init__(self):
        self.__dataPoints = {}
        self.__size = None

    def addDataPoint(self, indVar, depVars):
        self.__dataPoints[indVar] = depVars
        #Make sure that the collections of variables are consistent in length
        if self.__size == None:
            self.__size = len(depVars)
        elif self.__size != len(depVars):
            raise Exception("A list of dependant variables was added whose length is different than the established data size in the flowline")

    def getDataPoint(self, indVar):
        return self.__dataPoints[indVar]

    def getPointLocations(self):
        return sorted(self.__dataPoints.keys())

    #Returns a list of values for one of the stored dependat variables identified by position
    #in the dep. variable list
    def getDependantVariable(self, varPositionNumber):
        l = []
        for x in self.getPointLocations():
            l.append(self.getDataPoint(x)[varPositionNumber])
        return l

    def getDataSize(self):
        if self.__size == None:
            raise Exception("I don't know the size of the data yet!")
        else:
            return self.__size

    #A hack in case the endpoints are not exactly 0 or 1.  In this case,
    #we'll just extend the range out according to a constant extrapolation.
    #TODO: figure out something better to do in this case
    def fixRange(self):
        locs = self.getPointLocations()
        if locs[0] > 0:
            self.addDataPoint(0.0, self.getDataPoint(locs[0]))
        if locs[len(locs) - 1] < 1:
            self.addDataPoint(1.0, self.getDataPoint(locs[len(locs) - 1]))

    #Interpolates this flowline onto a list of alternate independant variable values
    #This is done linearly for now....
    def interpolateOntoValues(self, newPoints):
        #If passed another flowline, interpolate onto that flowline's independant variable locs.
        if type(newPoints) == Flowline:
            newPoints = newPoints.getPointLocations()

        interpolatedFlowline = Flowline()
        #Get the data points that we have
        myPoints = self.getPointLocations()
        #The data point of this flowline that is the lower value in the interpolation.
        currentPoint = 0
        for x in newPoints:
            while myPoints[currentPoint+1] < x:
                currentPoint += 1
            #Invariant: myPoints[currentPoint] <= x, myPoints[currentPoint+1] > x
            alpha = (x - myPoints[currentPoint]) / (myPoints[currentPoint+1] - myPoints[currentPoint])
            data1 = self.getDataPoint(myPoints[currentPoint])
            data2 = self.getDataPoint(myPoints[currentPoint+1])

            interp = [i*(1-alpha) + j*alpha for i,j in zip(data1, data2)]
            interpolatedFlowline.addDataPoint(x,interp)
        return interpolatedFlowline

#Represents a plan view of a 3D experiment.  Because of the display requirements,
#the plan view is internally stored as a list of flowlines, with each flowline
#associated with the y coordinate that it runs along.
class Planview(object):
    def __init__(self):
        self.__flowlines = {}

    def addDataPoint(self, x, y, depVars):
        if y not in self.__flowlines:
            self.__flowlines[y] = Flowline()
        self.__flowlines[y].addDataPoint(x, depVars)

    def extractFlowline(self, y):
        #If the requested y value is represented directly, grab it
        if y in self.__flowlines:
            return self.__flowlines[y]
        
        #Otherwise, we need to grab the flowline on either side of the
        #desired y value and interpolate.

        #For starters, figure out which y values we want.
        belowFl = 0
        aboveFl = 1
        for fl in self.__flowlines.iterkeys():
            if y - fl < y - belowFl:
                belowFl = fl
            if fl - y < aboveFl - y:
                aboveFl = fl

        flowline1 = self.__flowlines[belowFl]
        flowline2 = self.__flowlines[aboveFl]

        #Figure out the interpolation value
        alpha = float(y - belowFl) / (aboveFl - belowFl)

        interpolatedFlowline = Flowline()
        #Loop through each flowline and interpolate
        for x in flowline1.getPointLocations():
            values1 = flowline1.getDataPoint(x)
            values2 = flowline2.getDataPoint(x)
            interp = [i*(1-alpha) + j*alpha for i,j in zip(values1, values2)]
            interpolatedFlowline.addDataPoint(x, interp)
        return interpolatedFlowline

#Given a list of flowlines and a set of x coordinates, interpolates the flowlines
#onto the given x coordinates and returns a tuple of two flowlines with the means
#and standard deviations respectively.
def aggregateExperimentFlowlines(flowlines, xPoints):
    interpolatedFlowlines = [fl.interpolateOntoValues(xPoints) for fl in flowlines]
    n = len(interpolatedFlowlines)
    avgFlowline = Flowline()
    stdevFlowline = Flowline()

    for x in xPoints:
        sum = [0] * interpolatedFlowlines[0].getDataSize()
        sumOfSquares = [0] * interpolatedFlowlines[0].getDataSize()
        
        for data in [fl.getDataPoint(x) for fl in interpolatedFlowlines]:
            for i, d in enumerate(data):
                sum[i] = sum[i] + d
                sumOfSquares[i] = sumOfSquares[i] + d**2
        avg = [s/n for s in sum]
        stdev = [math.sqrt(s/n - m**2) for s,m in zip(sumOfSquares, avg)]
        avgFlowline.addDataPoint(x, avg)
        stdevFlowline.addDataPoint(x, stdev)
    return avgFlowline, stdevFlowline

def readFlowlineFile(filename):
    f = open(filename)
    flowline = Flowline()
    for line in f:
        data = [float(i) for i in line.strip().split()]
        if checkForNan(data):
            raise Exception("NaN found in input data!")
        x = data[0]
        data = data[1:]
        flowline.addDataPoint(x,data)
    f.close()
    return flowline

def readPlanviewFile(filename):
    f = open(filename)
    plan = Planview()
    for line in f:
        data = [float(i) for i in line.strip().split()]
        if checkForNan(data):
            raise Exception("NaN found in input data!")
        x = data[0]
        y = data[1]
        data = data[2:]
        plan.addDataPoint(x,y,data)
    f.close()
    return plan

def checkForNan(numberList):
    for n in numberList:
        if not n == n:
            return True
    return False
    

#Returns a new flowline with one data member: the norm of the surface velocity.
#Assumes that surface velocity components are the first two data members
def normSurfaceVelocity(flowline):
    result = Flowline()
    for x in flowline.getPointLocations():
        data = flowline.getDataPoint(x)
        result.addDataPoint(x, [math.sqrt(data[0]**2 + data[1]**2)])
    return result

def plotAggregatedFlowlines(axis, meanFlowline, stdevFlowline, dataMember, color, labelPrefix):
    xs = meanFlowline.getPointLocations()
    means = meanFlowline.getDependantVariable(dataMember)
    stdevs = stdevFlowline.getDependantVariable(dataMember)

    #Make sure the color is a tuple, not a list
    color = tuple(color)

    #Plot the mean
    axis.plot(xs, means, ":", color=color, label=labelPrefix + " Mean")
    #Shade the standard deviation around the mean
    #Set up a list of points on the polygon to shade.
    #First, we get a list of the upper and lower lines
    meanPlusSd =  [m+s for m,s in zip(means,stdevs)]
    meanMinusSd = [m-s for m,s in zip(means,stdevs)]
    
    #axis.plot(xs,meanPlusSd)
    #axis.plot(xs,meanMinusSd)
    
    #Now, we set up the list of points that form the polygon
    #This funky list reversing is so that we go around the polygon counter-clockwise;
    #that is, we specify in one list the upper bound going left-to-right, then the
    #lower bound going right-to-left.
    polyX = xs + list(reversed(xs))
    polyY = meanPlusSd + list(reversed(meanMinusSd))
    #Plot the polygon as the same color but at only 1/4 opacity
    axis.fill(polyX, polyY, facecolor=color,edgecolor=color,alpha=.25, label=labelPrefix + " Std. Dev.")

def createPlot(experiment, domainSizeKm):
    fig = pyplot.figure()
    axis = fig.add_subplot(111, xlabel="Rescaled X coordinate", ylabel="Velocity (m/a)", 
           title="ISMIP-HOM Experiment " + experiment.upper() + ", " + str(domainSizeKm) + " km")
    print axis.__class__
    for isFullStokes in [True, False]:
        models = getDataFiles("ismip_all", experiment, domainSizeKm, isFullStokes)
        print models
        
        #Certain ISMIP-HOM experiments are 3-D experiments that produce a plan view of the surface
        #velocity as their output.  If this is the case, we need to read in a plan view format 
        #file, extract a flowline at a standard location (y=L/4), and find the norm of the surface
        #velocity.  Other experiments are flowline experiments; that is, they only occur in the x and
        #z axis.  If this is the case, we can just read in the flowline file and use it as-is
        isFlowlineExp = experiment in ["b","d","e"] 
        if isFlowlineExp: #Experiment has flowline output already
            flowlines = []
            for f in models:
                try:
                    flowlines.append(readFlowlineFile(f))
                except Exception, e:
                    print f, "failed to load:",e
        else:
            plans = []
            for f in models:
                try:
                    plans.append(readPlanviewFile(f))
                except Exception, e:
                    print f, "failed to load:",e
            #Extract a flowline from the plan view at y=L/4
            flowlines = [pv.extractFlowline(.25) for pv in plans]
            #For each flowline, create another flowline with the norm of the
            #surface velocity.
            flowlines = [normSurfaceVelocity(fl) for fl in flowlines]

        #Simple fix if a flowline has a range that is less than [0..1].  In this case,
        #we are at the moment simply extending the boundaries out.
        for fl in flowlines:
            fl.fixRange()

        #Compute the mean and standard deviations of the experiments
        mean, stdev = aggregateExperimentFlowlines(flowlines, flowlines[0].getPointLocations())
        
        #Plot the mean and std. dev.
        if isFullStokes:
            plotAggregatedFlowlines(axis, mean, stdev, 0, (1,0,0), "Full Stokes")
        else:
            plotAggregatedFlowlines(axis, mean, stdev, 0, (0,0,1), "First Order")
    pyplot.legend(loc='center right')

    filename = experiment + "_" + str(domainSizeKm) + "km.png" 

    pyplot.savefig(filename)

if __name__ == "__main__":
    for experiment in ["a","b","c","d"]:
        for domainSize in [5, 10, 20, 40, 80, 160]:
            createPlot(experiment, domainSize)
