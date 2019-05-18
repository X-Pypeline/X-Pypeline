#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np

from gwpy.timeseries import TimeSeriesDict
from gwpy.table import EventTable

from .xdetector import Detector


class XEvent(object):
    def __init__(self):
        self.event = 0


class XCreateEventFromFile(XEvent):
    def __init__(self, paramsfile, event_number):
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
        parameters = dict(line.rstrip('\n').split(':', 1) for line in open(args.parameter_file))
import pdb
pdb.set_trace()
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
                               format='ascii')['col1'])[event_number]

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
