#!/usr/bin/python
#This python script replaces the Makefile that used to be in this directory.
#The makefile is replaced because it was growing difficult to avoid a lot of
#code duplication especially with the requirement that ISMIP-HOM A through D
#be run for different domain sizes.

import os
import sys
from getopt import gnu_getopt

import changeDomainSize

import intercomparison

#Parse command line arguments
optlist, args = gnu_getopt(sys.argv, intercomparison.ExperimentOptionsShort, 
         intercomparison.ExperimentOptionsLong + 
         ["format-only","horiz-grid-size=","vert-grid-size=","vert-grid-spacing=","prefix=","diagnostic-type=", "no-visuals", "glide="])

optdict = dict(optlist)

#Determine whether to actually run Glimmer or just prepare the input files
formatOnly = ("--format-only",'') in optlist
createVisuals = ("--no-visuals",'') not in optlist

#Determine which experiments and domain sizes to run
experiments, domainSizes = intercomparison.getExperimentsToRun(optlist)

#Determine whether the grid sizes specified in the input config files have been overridden by a command line argument
if "--horiz-grid-size" in optdict:
    horizontalGridOverride = int(optdict["--horiz-grid-size"])
else:
    horizontalGridOverride = None

if "--vert-grid-size" in optdict:
    verticalGridOverride = int(optdict["--vert-grid-size"])
else:
    verticalGridOverride = None

if "--vert-grid-spacing" in optdict:
    arg = optdict["--vert-grid-spacing"]
    vertGridOptions = {"glimmer":0, "even":1, "pattyn":2}
    verticalSpacingSchemeOverride = vertGridOptions[arg]
else:
    verticalSpacingSchemeOverride = None

if "--diagnostic-type" in optdict:
    diagnosticSchemeOverride = optdict["--diagnostic-type"]
else:
    diagnosticSchemeOverride = None

if "--glide" in optdict:
    glidePath = optdict["--glide"]
else:
    glidePath = "simple_glide"

#Determine the prefix to use for the ISMIP-HOM format output files
if "--prefix" in optdict:
    ismipHomOutputPrefix = optdict["--prefix"]
else:
    ismipHomOutputPrefix = "glm1"

for experiment in experiments:
    #Name of the configuration file for this experiment
    filename = "ishom." + experiment + ".config"
    #Name of the netcdf generation script for this experiment
    ncScript = sys.executable + " ismip_hom_" + experiment + ".py"
    for domain in domainSizes:
        #Create a new configuration file that reflects the new domain size and (possibly) grid spacing
        newConfigFile = changeDomainSize.changeDomainSize(filename, domain, 
                                          horizontalGridOverride, verticalGridOverride, 
                                          verticalSpacingSchemeOverride, diagnosticSchemeOverride)

        #Figure out the name of the netCDF file that the output will be written to
        #(this is a heuristic hack, it should probably really read the config file!!!!)
        ncOutputFilename = newConfigFile.replace("config","out.nc")
        
        #Name of the file that will be written to, using the ISMIP-HOM intercomparison standard
        intercompareOutputFilename = ismipHomOutputPrefix + experiment + "%03d"%int(domain/1000) + ".txt"
        
        #Create the domain size-specific config file
        #Generate the netcdf file for this config file
        os.system(ncScript + " " + newConfigFile)
        if not formatOnly:
            #Run Glimmer
            os.system("echo " + newConfigFile + "|" + glidePath)
            #Reformat the output as the ISMIP-HOM format
            os.system(sys.executable + " formatData.py " + experiment + " " + ncOutputFilename + " " + intercompareOutputFilename)

#If intercomparisons were run, create the output visuals
if (not formatOnly) and createVisuals:
    os.system(sys.executable + " createVisuals.py " + " ".join(sys.argv[1:]))
