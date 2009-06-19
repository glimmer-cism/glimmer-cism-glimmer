#!/usr/bin/python
#This python script replaces the Makefile that used to be in this directory.
#The makefile is replaced because it was growing difficult to avoid a lot of
#code duplication especially with the requirement that ISMIP-HOM A through D
#be run for different domain sizes.

import os
import sys
from optparse import OptionParser

import changeDomainSize
import intercomparison


parser=OptionParser()
intercomparison.standardOptions(parser)
parser.add_option("-f", "--format_only", dest="format_only", action="store_true", help = "If specified will stop after generating the NetCDF input")
parser.add_option("-o", "--horiz-grid-size", dest="horiz_grid_size", type="int", help="Alternate horizontal grid size to use (overrides config file)")
parser.add_option("-v", "--vert-grid-size", dest="vert_grid_size", type="int", help="Alternate vertical grid size to use (overrides config file)")
parser.add_option("-g", "--vert-grid-spacing", dest="vert_grid_spacing", help="Alternate vertical grid spacing to use (overrides config file)")
parser.add_option("-d", "--diagnostic-type", dest="diagnostic_type", help="Alternate diagnostic computation to use (overrides config file)")
parser.add_option("-n", "--no-visuals", dest="generate_visuals", action="store_false", help="Suppresses the creation of intercomparison visuals")
parser.add_option("-i", "--glide", dest="glide", default="simple_glide", help="Path to the GLIMMER runtime to use (defaults to simple_glide)")

options, args = parser.parse_args()
#If the user didn't specify a list of experiments or domain sizes, run the whole suite
if not options.experiments:
    options.experiments = ["a","b","c","d"]
if not options.sizes:
    options.sizes = ["5", "10", "20", "40", "80", "160"]

for experiment in options.experiments:
    #Name of the configuration file for this experiment
    filename = "ishom." + experiment + ".config"
    #Name of the netcdf generation script for this experiment
    ncScript = sys.executable + " ismip_hom_" + experiment + ".py"
    for domain in options.sizes:
        #Create a new configuration file that reflects the new domain size and (possibly) grid spacing
        newConfigFile = changeDomainSize.changeDomainSize(filename, int(domain)*1000, 
                                          options.horiz_grid_size, options.vert_grid_size, 
                                          options.vert_grid_spacing, options.diagnostic_type)

        #Figure out the name of the netCDF file that the output will be written to
        #(this is a heuristic hack, it should probably really read the config file!!!!)
        ncOutputFilename = newConfigFile.replace("config","out.nc")
        
        #Name of the file that will be written to, using the ISMIP-HOM intercomparison standard
        intercompareOutputFilename = options.prefix + experiment + "%03d"%int(domain) + ".txt"
        
        #Create the domain size-specific config file
        #Generate the netcdf file for this config file
        os.system(ncScript + " " + newConfigFile)
        if not options.format_only:
            #Run Glimmer
            os.system("echo " + newConfigFile + "|" + options.glide)
            #Reformat the output as the ISMIP-HOM format
            os.system(sys.executable + " formatData.py " + experiment + " " + ncOutputFilename + " " + intercompareOutputFilename)

#If intercomparisons were run, create the output visuals
if (not options.format_only) and options.generate_visuals:
    os.system(sys.executable + " createVisuals.py " + " ".join(sys.argv[1:]))
