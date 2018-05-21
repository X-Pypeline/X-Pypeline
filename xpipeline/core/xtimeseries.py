#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np

from gwpy.timeseries import TimeSeriesDict
from .xfrequencyseries import XFrequencySeriesDict
from .xtimefrequencymap import XTimeFrequencyMapDict, XTimeFrequencyMap
from collections import OrderedDict

class XTimeSeries(TimeSeriesDict):
    @classmethod
    def retrieve_data(cls, event_time, block_time,
                      channel_names, sample_frequency,
                      frame_types=[], verbose=False):
        """Obtain data for either on source, off source, or injections.

        This uses the gwpy `TimeSeriesDict.get` method

        Parameters
        ----------
        event_time : (`float)
            trigger time of event to be processing

        blcok_time : (`int`)
            length of data to be processed

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
        #----- Start and stop time for this event.
        startTime = event_time - block_time / 2;
        stopTime = event_time + block_time / 2;

        # Retrieve data and then resample and set epoch
        data = cls.get(
                        channel_names,
                        startTime, stopTime, verbose=verbose
                       )

        data.resample(sample_frequency)

        for (idet, iseries) in data.items():
            # set epoch of timeseries to the event_time
            iseries.epoch = event_time

        # set one of the detectors to be the reference detecotr
        # for any future coherent combinations

        return XTimeSeries(data)


    def asd(self, fftlength, **kwargs):
        """Obtain the asd of items in this dict.

        Parameters
        ----------
        fftlength : `dict`, `float`
            either a `dict` of (channel, `float`) pairs for key-wise
            asd calc, or a single float/int to computer as of all items.

        **kwargs
             other keyword arguments to pass to each item's asd
             method.
        """
        asds = XFrequencySeriesDict()
        if not isinstance(fftlength, dict):
            fftlength = dict((c, fftlength) for c in self)

        for key, fftlen in fftlength.items():
            asds[key] = self[key].asd(fftlen,
                                      fftlen/2.,
                                      method='lal_median_mean'
                                     )

        return asds


    def whiten(self, asds):
        """White this `XTimeSeries` against its own ASD

            Parameters
            ----------
            asds : `dict`
                a `dict` of (channel, `XFrequencySeries`) pairs for key-wise
                            whitened of the timeseries
        """
        if not isinstance(asds, dict):
            raise ValueError("asds must be supplied in the form of a dict")

        whitened_timeseries = XTimeSeries()
        for (idet, iseries) in self.items():
            whitened = iseries.whiten(asds[idet].dx,
                                      asds[idet].dx/2.,
                                      asd=asds[idet])
            whitened_timeseries.append(
                                       {idet : whitened}
                                       )

        return whitened_timeseries


    def spectrogram(self, fftlength):
        """Obtain the specotrograms of items in this dict.

        Parameters
        ----------
        fftlength : `dict`, `float`
            either a `dict` of (channel, `float`) pairs for key-wise
            asd calc, or a single float/int to computer as of all items.

        **kwargs
             other keyword arguments to pass to each item's asd
             method.
        """
        tfmaps = XTimeFrequencyMapDict()
        if not isinstance(fftlength, dict):
            fftlength = dict((c, fftlength) for c in self)

        for (idet, iseries) in self.items():
            tfmaps[idet] = XTimeFrequencyMap(iseries.spectrogram(
                                             stride=fftlength[idet],
                                             fftlength=fftlength[idet],
                                             overlap=0.5,
                                             window='hann'))
        return tfmaps
