#!/usr/bin/python
import intercomparison

#Given a flowline of values and the average and standard dev. flowlines returned from
#aggregateExperimentFlowlines, returns a tuple containing:
#1. The maximum observed error from the average.
#2. The maximum observed percent error from the average.
#3. Boolean that's true if the value was ever outside one standard deviation
def verifyFlowlineNumerically(flowline, avgFlowline, stdevFlowline, varPositionNumber):
    maxError = float("-inf")
    maxPercentError = float("-inf")
    outsideStdev = False
    for i, (val, avg, stdev) in enumerate(zip(flowline.getDependantVariable(varPositionNumber), 
                               avgFlowline.getDependantVariable(varPositionNumber), 
                               stdevFlowline.getDependantVariable(varPositionNumber))):
        error = abs(val-avg)
        percentError = 100 * error/(avg)
        maxError = max(maxError, error + 1e-10) #Regularization
        maxPercentError = max(maxPercentError, percentError)
        if val > avg + stdev or val < avg - stdev:
            outsideStdev = True
            print "Position", str(i), ":", val, "not in", "[" + str(avg - stdev) + ", " + str(avg+stdev) + "]"
    return maxError, maxPercentError, outsideStdev

def performRegressionTest(inputFile):
    experiment = inputFile[4]
    domainSizeKm = int(inputFile[5:8])

    #Read the data file from the Glimmer run
    glimFlowline = intercomparison.grabFlowline(inputFile)
   
    #Read the data files from the intercomparison experiment
    models = intercomparison.getDataFiles("ismip_all", experiment, domainSizeKm, "partial-stokes")
    flowlines = []
    for f in models:
        try:
            fl = intercomparison.grabFlowline(f)
            fl.fixRange()
            flowlines.append(fl)
        except Exception, e:
            print f, "failed to load:", e

        #Compute the mean and standard deviations of the experiments
        meanFlowline, stdevFlowline = intercomparison.aggregateExperimentFlowlines(flowlines, glimFlowline.getPointLocations())
 
    maxError, maxPercentError, outsideStdev = verifyFlowlineNumerically(glimFlowline, meanFlowline, stdevFlowline, 0)
    if outsideStdev:
        print "FAILED:",
    else:
        print "PASSED:",
    print "ISMIP-HOM", experiment, domainSizeKm,"km is within one standard deviation of other model outputs"
    print "Error = ", maxError, "% Error =", maxPercentError

if __name__ == "__main__":
    import sys

    performRegressionTest(sys.argv[1])
