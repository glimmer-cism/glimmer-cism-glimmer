#!/usr/bin/python
#This python script replaces the Makefile that used to be in this directory.
#The makefile is replaced because it was growing difficult to avoid a lot of
#code duplication especially with the requirement that ISMIP-HOM A through D
#be run for different domain sizes.

import os
import sys
from getopt import gnu_getopt

import changeDomainSize

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


#Parse command line arguments
optlist, args = gnu_getopt(sys.argv, 'abcd', 
    ["5km","10km","20km","40km","80km","160km",
     "format-only","horiz-grid-size=","vert-grid-size=","vert-grid-spacing=","ismip-hom-prefix="])
optdict = dict(optlist)

#Determine whether to actually run Glimmer or just prepare the input files
formatOnly = ("--format-only",'') in optlist

#Determine which experiments and domain sizes to run
experiments, domainSizes = getExperimentsToRun(optlist)

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

#Determine the prefix to use for the ISMIP-HOM format output files
if "--ismip-hom-prefix" in optdict:
    ismipHomOutputPrefix = optdict["--ismip-hom-prefix"]
else:
    ismipHomOutputPrefix = "glm1"

for experiment in experiments:
    #Name of the configuration file for this experiment
    filename = "ishom." + experiment + ".config"
    #Name of the netcdf generation script for this experiment
    ncScript = "python ismip_hom_" + experiment + ".py"
    for domain in domainSizes:
        #Create a new configuration file that reflects the new domain size and (possibly) grid spacing
        newConfigFile = changeDomainSize.changeDomainSize(filename, domain, 
                                          horizontalGridOverride, verticalGridOverride, verticalSpacingSchemeOverride)

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
            os.system("echo " + newConfigFile + "|simple_glide")
            #Reformat the output as the ISMIP-HOM format
            os.system("./formatData.py " + experiment + " " + ncOutputFilename + " " + intercompareOutputFilename)
