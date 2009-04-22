#!/usr/bin/python
import os

k = 10
guideFile = open("refinement-experiments.txt","w")
for gridSizeArg in ["20","40","80"]:
    for vertGridArg in ["20","40","80"]:
        for vertGridSpacing in ["glimmer","pattyn","even"]:
            for diagnosticScheme, strDiagnostic in [("1","unstaggered"),("2","staggered")]:
                
                print >>guideFile, k, gridSizeArg, vertGridArg, vertGridSpacing, strDiagnostic
                prefix = "gl"+str(k)
            
                args = []
                args.append("--horiz-grid-size=" +gridSizeArg)
                args.append("--vert-grid-size=" +vertGridArg)
                args.append("--vert-grid-spacing=" +vertGridSpacing)
                args.append("--prefix="+prefix)
                args.append("--diagnostic-type="+diagnosticScheme)

                verifyCommand = sys.executable + " verify.py -abc " + " ".join(args) + " >>refinement-log.txt"

                caption = gridSizeArg + "x" + gridSizeArg + "x" + vertGridArg + ", " + vertGridSpacing + " spacing, " + strDiagnostic

                visualsCommand = sys.executable + " createVisuals.py -abc --prefix=" + prefix + " '--subtitle=" + caption + "'"

                print "RUNNING:",verifyCommand
                os.system(verifyCommand)
                print "RUNNING:",visualsCommand
                os.system(visualsCommand)            
                k += 1

guideFile.close()
            
