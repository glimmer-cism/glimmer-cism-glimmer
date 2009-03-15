#!/usr/bin/python
from intercomparison import *

#Given a flowline of values and the average and standard dev. flowlines returned from
#aggregateExperimentFlowlines, returns a tuple containing:
#1. The maximum observed error from the average.
#2. The maximum observed percent error from the average.
#3. Boolean that's true if the value was ever outside one standard deviation
def verifyFlowlineNumerically(flowline, avgFlowline, stdevFlowline, varPositionNumber):
    maxError = float("-inf")
    maxPercentError = float("-inf")
    outsideStdev = False
    for val, avg, stdev in zip(flowline.getDependantVariable(varPositionNumber), 
                               avgFlowline.getDependantVariable(varPositionNumber), 
                               stdevFlowline.getDependantVariable(varPositionNumber)):
        error = abs(val-avg)
        percentError = 100 * error/(avg)
        maxError = max(maxError, error + 1e-10) #Regularization
        maxPercentError = max(maxPercentError, percentError)
        if val > avg + stdev or val < avg - stdev:
            outsideStdev = True
    
    return maxError, maxPercentError, outsideStdev
def performRegressionTest(inputFile, errorTolerance):
    experiment = inputFile[4]
    domainSizeKm = int(inputFile[5:8])

    #bla bla bla we get avg and stdev and values somehow
    #also, we need position = component we are verifying
    maxError, maxPercentError, outsideStdev = verifyFlowlineNumerically(flowline, avgFlowline, stdevFlowline, position)
    if outsideStdev:
        print "FAILED:",
    else:
        print "PASSED:",
    print "ISMIP-HOM", experiment, domainSizeKm,"km is within one standard deviation of other model outputs"
    print "Error = ", maxError, "% Error =", maxPercentError
