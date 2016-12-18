import sys
import os
import random
import string
import shutil
import ConfigParser
import optparse
import argparse
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
       Some commandine argument are optional and some are not
    """
    parser = argparse.ArgumentParser(description='Process arguments needed to run XPypeline.')
    parser.add_argument("--condor", action="store_true", default=False,help="Want to run X thorugh condor?") 
    parser.add_argument("--eventTime", type=float,help="Trigger time of the glitch")
    parser.add_argument('--ifos', nargs='+', default=['H1', 'L1', 'V1'], help='IFOs for the analysis.')
    parser.add_argument("--inifile", help="Name of ini file containing the parameters of the X-Pypeline run")
    parser.add_argument("--name", default="TESTRUN",help="The name of the event being followed-up")
    parser.add_argument("--NSDF", action="store_true", default=False,help="No framecache file available want to use NSDF server to retrieve your data.")
    parser.add_argument("--outDir", help="output files containing the triggers found from the search will go here")
    parser.add_argument('--ra', default=None, help='RA of trigger.')
    parser.add_argument('--dec', default=None, help='Dec of trigger')
    parser.add_argument("--verbose", action="store_true", default=False,help="Run in Verbose Mode")

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    pdb.set_trace()
    args = parse_commandline()
