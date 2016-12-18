import sys
import os
import random
import string
import shutil
import ConfigParser
import optparse
import json
import rlcompleter
import pdb

pdb.Pdb.complete = rlcompleter.Completer(locals()).complete

###############################################################################
##########################                             ########################
##########################   Func: parse_commandline   ########################
##########################                             ########################
###############################################################################
# Definite Command line arguments here

def parse_commandline():
    """Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()
    parser.add_option("--inifile", help="Name of ini file of params")
    parser.add_option("--eventTime", type=float,help="Trigger time of the glitch")
    parser.add_option("--uniqueID", action="store_true", default=False,help="Is this image being generated for the GravitySpy project, is so we will create a uniqueID strong to use for labeling images instead of GPS time")
    parser.add_option("--ID", default='',help="Already supplying an ID? If not then ignore this flag. Only to be used in conjunction with --uniqueID")
    parser.add_option("--outDir", help="Outdir of omega scan and omega scan webpage (i.e. your html directory)")
    parser.add_option("--NSDF", action="store_true", default=False,help="No framecache file available want to use NSDF server")
    parser.add_option("--condor", action="store_true", default=False,help="Want to run as condor job?")
    parser.add_option("--plot-whitened-timeseries", action="store_true", default=False,help="Plot whitened timeseries")
    parser.add_option("--plot-highpassfiltered-timeseries", action="store_true", default=False,help="Plot high pass filtered timeseries")
    parser.add_option("--plot-raw-timeseries", action="store_true", default=False,help="Plot raw timeseries")
    parser.add_option("--plot-eventgram", action="store_true", default=False,help="Plot eventgram")
    parser.add_option("--runML", action="store_true", default=False,help="Run the ML classifer on the omega scans")
    parser.add_option("--verbose", action="store_true", default=False,help="Run in Verbose Mode")
    opts, args = parser.parse_args()


    return opts

if __name__ == '__main__':
    main()
