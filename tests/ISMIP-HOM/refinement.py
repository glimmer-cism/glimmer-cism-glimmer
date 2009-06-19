#!/usr/bin/python
#The purpose of this script is to run the ISMIP-HOM suite for various grid refinements and solver options
#to study their impact on the model output.
import os
import sys

def runExperiment(gridSizeArg="40", vertGridArg="40", vertGridSpacing="even", diagnosticScheme="2"):
    global k
    global guideFile

    gridSizeArg= str(gridSizeArg)
    vertGridArg= str(vertGridArg)
    diagnosticScheme= str(diagnosticScheme)

    strDiagnostic = {"1":"unstaggered", "2":"staggered"}[diagnosticScheme]

    print >>guideFile, k, gridSizeArg, vertGridArg, vertGridSpacing, strDiagnostic
    prefix = "gl"+str(k)
            
    args = []
    args.append("--horiz-grid-size=" +gridSizeArg)
    args.append("--vert-grid-size=" +vertGridArg)
    args.append("--vert-grid-spacing=" +vertGridSpacing)
    args.append("--prefix="+prefix)
    args.append("--diagnostic-type="+diagnosticScheme)

    verifyCommand = sys.executable + " verify.py --exp=a,c " + " ".join(args) + " >>refinement-log.txt"

    caption = gridSizeArg + "x" + gridSizeArg + "x" + vertGridArg + ", " + vertGridSpacing + " spacing, " + strDiagnostic

    visualsCommand = sys.executable + " createVisuals.py -abc --prefix=" + prefix + " '--subtitle=" + caption + "'"

    print "RUNNING:",verifyCommand
    os.system(verifyCommand)
    print "RUNNING:",visualsCommand
    os.system(visualsCommand)            
    k += 1
	

def fullCombinatorics():
    for gridSizeArg in ["20","40","60","80"]:
        for vertGridArg in ["20","40",,"60","80"]:
            for vertGridSpacing in ["glimmer","pattyn","even"]:
                for diagnosticScheme in ["1","2"]:
                    runExperiment(gridSizeArg, vertGridArg, vertGridSpacing, diagnosticScheme)

def isolatedVariables():
    #First do all defaults
    runExperiment()
    #Try varying horizontal grid spacing
    runExperiment(gridSizeArg=20)
    runExperiment(gridSizeArg=60)
    runExperiment(gridSizeArg=80)

    #Try varying vertical grid spacing
    runExperiment(vertGridArg=20)
    runExperiment(vertGridARg=60)
    runExperiment(vertGridArg=80)

    #Try different grid spacings
    runExperiment(vertGridSpacing = "pattyn")
    runExperiment(vertGridSpacing = "glimmer")

    #Try a colocated grid
    runExperiment(diagnosticScheme = 1)

k = 10
guideFile = open("refinement-experiments.txt","w")


#Uncomment to run the fully combinatoric version of the experiment
#fullCombinatorics()

isolatedVariables()

guideFile.close()

