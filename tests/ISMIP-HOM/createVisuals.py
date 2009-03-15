#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import os.path
import math
from matplotlib import pyplot
import matplotlib.lines
import matplotlib.patches
import matplotlib.figure
import matplotlib.font_manager
import numpy
from getopt import gnu_getopt

from intercomparison import *


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

#Returns a tuple with two tuples suitable for creating a legend.  The first tuple contains
#the lines and patches, the second contains the names
def createPlot(experiment, domainSizeKm, fig, subplotNum, notFullStokesModelType):
    isFlowlineExp = experiment in ["b","d","e"] 
    
    axis = fig.add_subplot(2,3,subplotNum)
    axis.set_title(str(domainSizeKm) + " km", size="medium")
    #Read the data that came from the Glimmer run
    glimFileName = glimPrefix + experiment + "%03d"%domainSizeKm + ".txt"
    if isFlowlineExp:
        glimFlowline = readFlowlineFile(glimFileName)
    else:
        glimPlan = readPlanviewFile(glimFileName, transpose=True)
        glimPlan = normSurfaceVelocityPlan(glimPlan)
        #glimPlan.plot()
        glimFlowline = glimPlan.extractFlowline(.25)
    #DEBUG!!
    #m = -9999999999999999999
    #for x in glimFlowline.getPointLocations():
    #    m = max(m, glimFlowline.getDataPoint(x)[0])
    #print m
    #End debug

    axis.plot(glimFlowline.getPointLocations(), glimFlowline.getDependantVariable(0), color=(0,0,0))

    for isFullStokes in [True, False]:
        if isFullStokes:
            modelType = "full-stokes"
        else:
            modelType = notFullStokesModelType
        models = getDataFiles("ismip_all", experiment, domainSizeKm, modelType)
        
        #Certain ISMIP-HOM experiments are 3-D experiments that produce a plan view of the surface
        #velocity as their output.  If this is the case, we need to read in a plan view format 
        #file, extract a flowline at a standard location (y=L/4), and find the norm of the surface
        #velocity.  Other experiments are flowline experiments; that is, they only occur in the x and
        #z axis.  If this is the case, we can just read in the flowline file and use it as-is
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
            plans = [normSurfaceVelocityPlan(pv) for pv in plans]
            #Extract a flowline from the plan view at y=L/4
            flowlines = [pv.extractFlowline(.25) for pv in plans]
        
        #Simple fix if a flowline has a range that is less than [0..1].  In this case,
        #we are at the moment simply extending the boundaries out.
        for fl in flowlines:
            fl.fixRange()

        #Compute the mean and standard deviations of the experiments
        mean, stdev = aggregateExperimentFlowlines(flowlines, glimFlowline.getPointLocations())
        
        #Plot the mean and std. dev.
        if isFullStokes:
            plotAggregatedFlowlines(axis, mean, stdev, 0, (1,0,0), "Full Stokes")
        else:
            plotAggregatedFlowlines(axis, mean, stdev, 0, (0,0,1), "First Order")

        #Modify the axes tick lables to have smaller fonts
        for tick in axis.xaxis.get_major_ticks():
            tick.label1.set_fontsize("xx-small")
        for tick in axis.yaxis.get_major_ticks():
            tick.label1.set_fontsize("xx-small")
        

if __name__ == "__main__":
    optlist, args = gnu_getopt(sys.argv[1:], "abcd", ["lmla", "prefix=", "subtitle="])
    optdict = dict(optlist)
    
    if "--lmla" in optdict:
        notFullStokesModelType = "lmla"
    else:
        notFullStokesModelType = "partial-stokes"

    if "--prefix" in optdict:
        glimPrefix = optdict["--prefix"]
    else:
        glimPrefix = "glm1"

    if "--subtitle" in optdict:
        subtitle = optdict["--subtitle"]
    else:
        subtitle = None

    #Get the list of experiments
    allTheExperiments = ['a','b','c','d']
    experiments = []
    
    for opt in optdict.keys():
        opt=opt.replace("-","")
        if opt in allTheExperiments:
            experiments.append(opt)
    if not experiments:
        experiments = allTheExperiments

    for experiment in experiments:
        print "ISMIP-HOM",experiment.upper()

        fig = pyplot.figure(subplotpars=matplotlib.figure.SubplotParams(top=.85,bottom=.15))
        for i, domainSize in enumerate([5, 10, 20, 40, 80, 160]):
            createPlot(experiment, domainSize, fig, i+1, notFullStokesModelType)

        #Create the legend!  This is overly complicated because I want the legend to be
        #in multiple columns at the bottom of the visual, and that's impossible with my
        #version of matplotlib (it's apparently coming though?)
        #So, I am creating a different legend for each column
        #l.draw_frame(False) turns off the borders.
        l = fig.legend([matplotlib.lines.Line2D([1,2],[1,2],color=(0,0,0))],
                   ['Model Output'] ,
                   loc=(.1, .05),prop=matplotlib.font_manager.FontProperties(size="x-small"))
        l.draw_frame(False)
        l = fig.legend([matplotlib.lines.Line2D([1,2],[1,2],ls=':',color=(1,0,0)),
                    matplotlib.patches.Patch(edgecolor=None, facecolor=(1,0,0), alpha=.25)],
                    ['Full Stokes Mean', 'Full Stokes Std. Dev.'], 
                   loc=(.3, .02),prop=matplotlib.font_manager.FontProperties(size="x-small"))
        l.draw_frame(False)
        l = fig.legend([matplotlib.lines.Line2D([1,2],[1,2],ls=':',color=(0,0,1)),
                    matplotlib.patches.Patch(edgecolor=None, facecolor=(0,0,1), alpha=.25)],
                    ['First Order Mean', 'First Order Std. Dev.'],
                   loc=(.55, .02),prop=matplotlib.font_manager.FontProperties(size="x-small"))
        l.draw_frame(False)

        #Add an overall title to the figure
        fig.text(.5, .92, "ISMIP-HOM Experiment " + experiment.upper(), horizontalalignment='center', size='large')
        if subtitle:
            fig.text(.5, .89, subtitle, horizontalalignment='center', size='small')
        
        #Add one axis label per side (the labels are the same because of the small multiple format,
        #so we only want to repeat them once in the whole figure!)
        fig.text(.5, .1, "Normalized X coordinate", horizontalalignment='center',size='small')
        fig.text(.06, .5, "Ice Speed (m/a)", rotation="vertical", verticalalignment='center') 
        #Save the figure
        filename = "ISMIP-HOM-" + experiment.upper() + "-" + glimPrefix + ".png" 
        pyplot.savefig(filename)

