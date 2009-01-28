#!/usr/bin/python
#This python script replaces the Makefile that used to be in this directory.
#The makefile is replaced because it was growing difficult to avoid a lot of
#code duplication especially with the requirement that ISMIP-HOM A through D
#be run for different domain sizes.

import os
import sys
from getopt import gnu_getopt

import changeDomainSize

optlist, args = gnu_getopt(sys.argv, '', ["short","format-only"])

formatOnly = ("--format-only",'') in optlist

if ("--short", '') in optlist:
    domain_sizes = [40000]
else:
    domain_sizes = [5000, 10000, 20000, 40000, 80000, 160000]

#Run ISMIP-HOM A-D
for experiment in ['a', 'b', 'c', 'd']:
    #Name of the configuration file for this experiment
    filename = "ishom." + experiment + ".config"
    #Name of the netcdf generation script for this experiment
    ncScript = "python ismip_hom_" + experiment + ".py"
    for domain in domain_sizes:
        if not formatOnly:
            #Create the domain size-specific config file
            newConfigFile = changeDomainSize.changeDomainSize(filename, domain)
            #Generate the netcdf file for this config file
            os.system(ncScript + " " + newConfigFile)
            #Run Glimmer
            os.system("echo " + newConfigFile + "|simple_glide")
            #Reformat the output as the ISMIP-HOM format
            ncOutputFilename = newConfigFile.replace("config","out.nc")
            intercompareOutputFilename = "glm1" + experiment + "%03d"%int(domain/1000) + ".txt"
        os.system("./formatData.py " + experiment + " " + ncOutputFilename + " " + intercompareOutputFilename)
