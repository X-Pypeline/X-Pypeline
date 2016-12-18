#!/usr/bin/env python

# ---- Import standard modules to the python path.

from __future__ import division

import sys
import os
import random
import string
import shutil
import json
import rlcompleter
import pdb

import numpy as np

from scipy import signal
from scipy.interpolate import InterpolatedUnivariateSpline


from gwpy.timeseries import TimeSeries

pdb.Pdb.complete = rlcompleter.Completer(locals()).complete

args = parse_commandline()

if __name__ == '__main__':
    ###########################################################################
    #                                   Parse Ini File                        #
    ###########################################################################

    # ---- Create configuration-file-parser object and read parameters file.
    cp = ConfigParser.ConfigParser()
    cp.read(opts.inifile)

    # ---- Read needed variables from [parameters] and [channels] sections.
    # ---- Window for FFTs.
    windowType = 'modifiedhann';
# ---- Path to injection file.
injectionFileName = '';
# ---- Scale factor applied to all injections (including MDC channels).
injectionScales = 0;
# ---- If nonempty replace frame data with simulated noise of the indicated type.
makeSimulatedNoise = '';
# ---- Path to file listing MDC channels for simulations.
mdcChannelFileName = '';
# ---- Path to file listing parameters used to miscalibrate injections.
misCalibFileName = '';
# ---- If nonzero take coincidence of clusters with injection list 
#      and report coincidences instead of all clusters.
postProcessInjections = 0;
# ---- Rescale injections by sum of Fp^2 over detectors (use with care).
rescaleByAntennaResponse = 0;
# ---- Compute sky map quantities at exact locations of injected signals.
testSourcePosition = 0;
# ---- 1 for extra verbosity.
verboseFlag = 1;
# ---- Length of whitening filters to be made by xcondition.  Should only 
#      change from 1 with changes to xcondition.
whiteningTime = 1;
# ---- Path to directory containing cataloged waveforms.
catalogDirectory = '';
# ---- Flag indicating whether to apply corrections to known calibration
#      errors
applyCalibCorrection = 0;
# ---- Use xfindchirpcondition instead of xcondition for data conditioning.
#      THIS OPTION IS NOW DEPRECATED - the value is ignored.
useXFindChirpCondition = 0;
# ---- Use xbayesiantimefrequencymap instead of xtimefrequencymap.
usexbayesiantimefrequencymap = 0;
# ---- Trigger rate used for per analysis time, sky position decimate in Hz
decimateRate = 1;
# ---- Trigger rate used for final decimation after super clustering in Hz
superDecimateRate = 0.25;
# ---- Black pixel percentile level used to determine which fraction of
#      pixels to discard from cluster construction
blackPixelPrctile = 99;
# ---- Threshold for generalized clustering, if above 0 then clusters
#      with less energy than this times the black pixel threshold are
#      rejected and not used in generalized clustering (connectivity 24).
genClusteringThresh = 0;
# ---- detection statistic to use when comparing triggers
detectionStat = [];
# ---- time step to use for circular time slides, zero if no circular time
# slides 
circTimeSlideStep = 0;
seedlessParams = '';
###########################################################################
#                 ensure output arguments are defined                     #
###########################################################################

# ---- Make sure all possible outputs are defined.

# ---- Output argument list for interactive use.
skyPositions = [];
likelihoodMap = [];
skyPositionIndex = [];
sourcePositions = [];
sourceLikelihoodMap = [];
sourcePositionIndex = [];
spectrogram = [];


        case 'analysismode',
            analysisMode = paramValStr{1};
        case 'analysistimes',
            analysisTimes = str2num(paramValStr{1});
        case 'applycalibcorrection',
            applyCalibCorrection = str2num(paramValStr{1});
        case 'blackpixelprctile',
            blackpixelprctile = str2num(paramValStr{1});
        case 'blocktime',
            blockTime = str2num(paramValStr{1});
        case 'catalogdirectory'
            catalogDirectory = paramValStr{1};
        case 'channelfilename',
            channelFileName = paramValStr{1};
        case 'circtimeslidestep',
            circTimeSlideStep = str2num(paramValStr{1});
        case 'decimaterate'
            decimateRate = str2num(paramValStr{1});
        case 'detectionstat',
            detectionStat = paramValStr{1};
        case 'donotinject'
            doNotInject = str2num(paramValStr{1});
        case 'eventfilename',
            eventFileName = paramValStr{1};
        case 'extractparamfilename'
            extractParamFileName = paramValStr{1};
        case 'filterfile'
            filterFile = paramValStr{1};
        case 'framecachefile',
            frameCacheFile = paramValStr{1};
        case 'frequencybands',
            frequencyBands = str2num(paramValStr{1});
        case 'genclusteringthresh',
            genClusteringThresh = str2num(paramValStr{1});
        case 'injectionfilename',
            injectionFileName = paramValStr{1};
        case {'injectionscale','injectionscales'},
            injectionScales = str2num(paramValStr{1});
        case 'likelihoodtype',
            likelihoodType = dataread('string',paramValStr{1},...
                            '#s','delimiter',',');
        case 'makesimulatednoise',
                makeSimulatedNoise = paramValStr{1};
        case 'maximumfrequency',
            maximumFrequency = str2num(paramValStr{1});
        case 'maximumtimingerror',
            maximumTimingError = str2num(paramValStr{1});
        case 'mdcchannelfilename',
            mdcChannelFileName = paramValStr{1};
        case 'minimumfrequency',
            minimumFrequency = str2num(paramValStr{1});
        case 'miscalibfilename',
            misCalibFileName = paramValStr{1};
        case 'multipleinjections'
            multipleInjections = str2num(paramValStr{1});
        case 'offsetfraction',
            offsetFraction = str2num(paramValStr{1});
        case 'outputtype',
            outputType = paramValStr{1};
        case 'postprocessinjections',
            postProcessInjections = str2num(paramValStr{1});
        case 'rescalebyantennaresponse',
            rescaleByAntennaResponse = str2num(paramValStr{1});
        case 'samplefrequency',
            sampleFrequency = str2num(paramValStr{1});
        case 'savevariables',
            saveVariables = dataread('string',paramValStr{1},...
                            '#s','delimiter',',');
        case 'seed',
            seed = paramValStr{1};
        case 'seedlessparams',
            seedlessParams = paramValStr{1};
        case 'skycoordinatesystem',
            skyCoordinateSystem = paramValStr{1};
        case 'skypositionlist',
            skyPositionList = paramValStr{1};
        case 'sphradparameterfile',
            sphradParameterFile = paramValStr{1};
        case 'superdecimaterate'
            superDecimateRate = str2num(paramValStr{1});
        case 'testsourceposition',
            testSourcePosition = str2num(paramValStr{1});
        case 'usexfindchirpcondition',
            useXFindChirpCondition = str2num(paramValStr{1});
        case 'usexbayesiantimefrequencymap',
            usexbayesiantimefrequencymap = str2num(paramValStr{1});
        case 'verboseflag',
            verboseFlag = str2num(paramValStr{1});
        case 'whiteningtime',
            whiteningTime = str2num(paramValStr{1});
        case 'windowtype',
            windowType = paramValStr{1};
        case {'onsourcebeginoffset','onsourceendoffset'},
        case{'doasymmetricbackground', 'backgroundasymmetryfactor'},

# ---- If no detection statistic specified, use the first specified
#      likelihood. This is backward compatible with xpipeline before r3470.
if isempty(detectionStat)
  detectionStat = likelihoodType{1};
elseif sum(strcmp(detectionStat,likelihoodType)) ~= 1
  likelihoodType
  error(['Detection statistic: ' detectionStat ' does not appear exactly once ' ...
         'in the list of computed likelihoods.']);
end

#----- Numeric value of injection number(s).
injectionNumbers = jobstr2num(injectionNumberString);
if length(injectionNumbers)>1
    warning('LINEAR test code: length(injectionNumbers)>1.');
end


