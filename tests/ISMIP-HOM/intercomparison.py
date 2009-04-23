# -*- coding: utf-8 -*-
#Library that contains routines for loading and manipulating data
#from the intercomparison experiments
import os
import math

#Lists of model classifications
fullStokes = ["aas1", "aas2", "cma1", "fpa2", "ghg1", "jvj1", "mmr1", "oga1", "rhi1", "rhi3", "spr1", "ssu1", "yko1"]

lmla = ["ahu1", "ahu2", "bds1", "fpa1", "mbr1", "rhi2", "tpa1"]

def readFlowlineFile(filename):
    f = open(filename)
    flowline = Flowline()
    for line in f:
        data = [float(i) for i in line.strip().split()]
        fixNan(data)
        x = data[0]
        data = data[1:]
        flowline.addDataPoint(x,data)
    f.close()
    return flowline

def readPlanviewFile(filename, transpose=False):
    f = open(filename)
    plan = Planview()
    for line in f:
        data = [float(i) for i in line.strip().split()]
        fixNan(data)
        x = data[0]
        y = data[1]
        data = data[2:]
        if not transpose:
            plan.addDataPoint(x,y,data)
        else:
            plan.addDataPoint(y,x,data)
    f.close()
    return plan

def checkForNan(numberList):
    for n in numberList:
        if not n == n:
            return True
    return False
    
def fixNan(data):
    for i,n in enumerate(data):
        if not n == n:
            data[i] = 0

def getDataFiles(root, experiment, domainSizeKm, modelType):    
    fnameSuffix = experiment + "%03d"%domainSizeKm +".txt"
    files = []
    for root, ds, fs in os.walk(root):
        for file in fs:
            if file.endswith(fnameSuffix):
                model = file[:4]
                if modelType == "full-stokes"  and model in fullStokes or \
                   modelType == "partial-stokes" and model not in fullStokes or \
                   modelType == "lmla" and model in lmla:
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
        if self.__size != len(depVars):
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
    
    def length(self):
        return len(self.__dataPoints)

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

    def getPointLocations(self):
        locations = []
        for y,fl in self.__flowlines.iteritems():
            for x in fl.getPointLocations():
                locations.append((x,y))
        return locations

    #surface velocity.
    def getDataPoint(self, x, y):
        return self.__flowlines[y].getDataPoint(x)

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
            if fl <= y and y - fl < y - belowFl:
                belowFl = fl
            if fl >= y and fl - y < aboveFl - y:
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

    def plot(self):
        xes = self.__flowlines.values()[0].getPointLocations()
        yes = sorted(list(self.__flowlines.keys()))
        nx = len(xes)
        ny = len(yes)
        arr = numpy.zeros((nx, ny))
        for i in range(nx):
            for j in range(ny):
                arr[i,j] = self.getDataPoint(xes[i], yes[j])[0]

        pyplot.imshow(arr)
        pyplot.show()

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

#Returns a new flowline with one data member: the norm of the surface velocity.
#Assumes that surface velocity components are the first two data members
def normSurfaceVelocity(flowline):
    result = Flowline()
    for x in flowline.getPointLocations():
        data = flowline.getDataPoint(x)
        result.addDataPoint(x, [math.sqrt(data[0]**2 + data[1]**2)])
    return result

def normSurfaceVelocityPlan(planview):
    result = Planview()
    for x, y in planview.getPointLocations():
        data = planview.getDataPoint(x,y)
        result.addDataPoint(x,y,[math.sqrt(data[0]**2 + data[1]**2)])
    return result



#Certain ISMIP-HOM experiments are 3-D experiments that produce a plan view of the surface
#velocity as their output.  If this is the case, we need to read in a plan view format 
#file, extract a flowline at a standard location (y=L/4), and find the norm of the surface
#velocity.  Other experiments are flowline experiments; that is, they only occur in the x and
#z axis.  If this is the case, we can just read in the flowline file and use it as-is
def grabFlowline(filename):
    experiment = filename[-8]
    isFlowlineExp = experiment in ["b","d","e"]
    if isFlowlineExp:
        return readFlowlineFile(filename)
    else:
        plan = readPlanviewFile(filename, transpose=True)
        plan = normSurfaceVelocityPlan(plan)
        return plan.extractFlowline(.25)

#STANDARD OPTION PARSING

#Lists of options to use with getopt.  Append these to whatever other options you want
#or use them straight out.
ExperimentOptionsShort = "abcd"
ExperimentOptionsLong  = ["5km","10km","20km","40km","80km","160km"]

#Parses a list of arguments extracted by getopt.  The user may specify specific experiments and/or domain
#sizes to run on the command line.  For example, -ab --5km --10km runs only experiments a and b for 5km and
#10km domains.  If nothing is specified for one of the lists, everything is run.  For example, -ac runs
#experiments a and c for all domain sizes, and --40km runs all experiments for a 40 km domain.
#Returns an (experimentList, domainSizes) tuple
def getExperimentsToRun(optlist):
    domainSizesDict = {"5km":5000, "10km":10000, "20km":20000, "40km":40000, "80km":80000, "160km":160000}
    experimentList = ['a','b','c','d']
    
    experimentDefaults = True
    experiments = []
    domainSizeDefaults = True
    domainSizes = []

    for option, arg in optlist:
        option = option.replace("-","")
        if option in experimentList:
            experimentDefaults = False
            experiments.append(option)
        elif option in domainSizesDict:
            domainSizeDefaults = False
            domainSizes.append(domainSizesDict[option])

    if experimentDefaults:
        experiments = experimentList
    if domainSizeDefaults:
        domainSizes = sorted(domainSizesDict.values())
    
    return experiments, domainSizes



