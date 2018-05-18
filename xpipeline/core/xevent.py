#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np

from gwpy.timeseries import TimeSeriesDict
from gwpy.table import EventTable

from .xdetector import Detector


class XEvent(object):
    def compute_antenna_patterns(self):
        """White this `TimeSeries` against its own ASD

            Parameters
            ----------
            fft_length : `float`
                number of seconds in single FFT
        """
        Fp = {}
        Fc = {}
        Fb = {}
        FL = {}
        F1 = {}
        F2 = {}
        for (idet_name, idet) in self.detectors.items():
            [Fptmp, Fctmp, Fbtmp, FLtmp, F1tmp, F2tmp] = \
                idet.compute_antenna_response(self.phi, self.theta)
            Fp[idet_name] = Fptmp
            Fc[idet_name] = Fctmp
            Fb[idet_name] = Fbtmp
            FL[idet_name] = FLtmp
            F1[idet_name] = F1tmp
            F2[idet_name] = F2tmp

        self.Fp = Fp
        self.Fc = Fc
        self.Fb = Fb
        self.FL = FL
        self.F1 = F1
        self.F2 = F2
        return Fp, Fc, Fb, FL, F1, F2


    def compute_time_delays(self):
        """White this `TimeSeries` against its own ASD

            Parameters
            ----------
            fft_length : `float`
                number of seconds in single FFT
        """
        time_delays = {}
        for (idet_name, idet) in self.detectors.items():
            self.reference_detector = idet_name
            time_delays[idet_name] = \
                idet.time_delay_from_earth_center_phi_theta(self.phi, self.theta)
        if len(self.detectors) == 1:
            self.time_delays = time_delays
            return time_delays
        else:
            # shift everything in relation to the first detector
            for idet, itime_delay in time_delays.items():
                time_delays[idet] = itime_delay - time_delays[self.reference_detector]
        self.time_delays = time_delays
        return time_delays


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
