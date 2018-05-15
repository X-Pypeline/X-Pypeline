#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np

from gwpy.timeseries import TimeSeriesDict
from gwpy.table import EventTable

from .xdetector import Detector


class XEvent(object):
    def retrieve_data(self, channel_names=[], frame_types=[], verbose=False):
        """Obtain data for either on source, off source, or injections.

        This uses the gwpy `TimeSeriesDict.get` method

        Parameters
        ----------
        channel_names (`list`) :
            required data channels.

        frame_types : `str`, optional
            name of frametype in which this channel is stored, by default
            will search for all required frame types

        verbose : `bool`, optional
            print verbose output about NDS progress.

        Returns:

            `TimeSeriesDict` :
        """
        if hasattr(self, 'channel_names'):
            channel_names = self.channel_names
        if hasattr(self, 'frame_types'):
            frame_types = self.frame_types

        #----- Start and stop time for this event.
        startTime = self.event_time - self.blocktime / 2;
        stopTime = self.event_time + self.blocktime / 2;

        # zip frameTypes and detectors, and channel names and detectors
        data = TimeSeriesDict.get(
                        channel_names,
                        startTime, stopTime, verbose=verbose
                       )

        for (idet, iseries) in data.items():
            # resample data
            if iseries.sample_rate.decompose().value != self.samplefrequency:
                data[idet] = iseries.resample(self.samplefrequency)

        self.timeseries = data
        return data


    def whiten(self):
        """White this `TimeSeries` against its own ASD
            
            Parameters
            ----------
            fft_length : `float`
                number of seconds in single FFT
        """
        whitened_timeseries = TimeSeriesDict()
        for (idet, iseries) in self.timeseries.items():
            asd = iseries.asd(self.whiteningtime,
                              self.whiteningtime/2.,
                              method='lal_median_mean')
            whitened = iseries.whiten(self.whiteningtime,
                                      self.whiteningtime/2.,
                                      asd=asd)
            whitened_timeseries.append(
                                       {idet : whitened}
                                       )

        self.whitened_timeseries = whitened_timeseries
        return whitened_timeseries


    def compute_antenna_patterns(self):
        # dimensions are number of skypostions by detectors
        # currently 6 patterns are calculated
        Fp = np.zeros([len(self.phi), len(self.detectors.keys())]) 
        Fc = np.zeros([len(self.phi), len(self.detectors.keys())]) 
        Fb = np.zeros([len(self.phi), len(self.detectors.keys())]) 
        FL = np.zeros([len(self.phi), len(self.detectors.keys())]) 
        F1 = np.zeros([len(self.phi), len(self.detectors.keys())]) 
        F2 = np.zeros([len(self.phi), len(self.detectors.keys())]) 
        for idet in self.detectors.items():
            idet = 1
        


class XCreateEventFromFile(XEvent):
    def __init__(self, paramsfile, eventNumber):
        """Initialize an XEvent on-source off-source or injection with pfile

        Parameters
        ----------
        paramsfile (str):
            a xpipeline param file
        eventNumber (int):
            an integer refering to what events from the
            input/event_off/on/inj.txt to grab for processing

        Returns:

            `XEvent`
        """
        print("You are generating an xevent by supplying a "
              "a xpipeline params file, this will overwite the defaults")
        with open(paramsfile, 'r') as f:
            for line in f.readlines():
                parsed_text = line.split('\n')[0].split(':')
                # check if param is also command separated
                try:
                    parsed_text[1].split(',')[1]
                    setattr(self, parsed_text[0], parsed_text[1].split(','))
                except:
                    setattr(self, parsed_text[0], parsed_text[1])

        self.phi = list(EventTable.read(self.skyPositionList,
                        format='ascii')['col2'])
        self.theta = list(EventTable.read(self.skyPositionList,
                          format='ascii')['col1'])
        self.event_time = list(EventTable.read(self.eventFileName,
                               format='ascii')['col1'])[eventNumber]

        for key, item in vars(self).items():
            try:
                setattr(self, key, float(item))
            except:
                pass

        analysistimes = [float(i) for i in self.analysistimes]
        self.analysistimes = analysistimes
        channel_names = []
        frame_types = []
        detectors = {}
        with open(self.channelFileName, 'r') as f:
            for det in f.readlines():
                detector_name = det.split(' ')[0].split(':')[0]
                channel_names.append(det.split(' ')[0])
                frame_types.append(det.split(' ')[1])
                detectors[detector_name] = Detector(detector_name)

        self.channel_names = channel_names
        self.frame_types = frame_types
        self.detectors = detectors


class XCreateEvent(XEvent):
    def __init__(self, event_time, phi, theta,
                 likelihoodtype=['loghbayesiancirc', 'standard', 'circenergy',
                                 'circinc', 'circnullenergy', 'circnullinc',
                                 'nullenergy', 'nullinc', 'powerlaw',
                                 'energyitf1', 'energyitf2', 'energyitf3',
                                 'skypositiontheta', 'skypositionphi'],
                 analysistimes=[1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125,
                                0.015625,0.0078125],
                 blocktime=256, minimumfrequency=64, maximumfrequency=500,
                 offsetfraction=0.5, outputtype='clusters',
                 samplefrequency=1024, whiteningtime=1.0,
                 seed=1235, applycalibcorrection=False,
                 onsourcebeginoffset=-600, onsourceendoffset=60,
                 makesimulatednoise=False,
                 circtimeslidestep=3, verboseflag=True):
        """Initialize an XEvent on-source off-source or injection.

        Parameters
        ----------
        event_time (float):
        block_time (int, optional):
        sample_frequency (optional, float): 
        rightascension (float, optional):
        declination (float, optional):

        Returns:

            `XEvent`
        """
        self.event_time = event_time
        self.phi = phi
        self.theta = theta
        self.likelihoodtype = likelihoodtype
        self.analysistimes = analysistimes
        self.blocktime = blocktime
        self.minimumfrequency = minimumfrequency
        self.maximumfrequency = maximumfrequency
        self.offsetfraction = offsetfraction
        self.outputtype = outputtype
        self.samplefrequency = samplefrequency
        self.verboseflag = verboseflag
        self.whiteningtime = whiteningtime
        self.seed = seed
        self.applycalibcorrection = applycalibcorrection
        self.onsourcebeginoffset = onsourcebeginoffset
        self.onsourceendoffset = onsourceendoffset
        self.makesimulatednoise = makesimulatednoise
        self.circtimeslidestep = circtimeslidestep
