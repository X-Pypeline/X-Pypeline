# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017-)
#
# This file is part of the XPypeline python package.
#
# hveto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hveto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hveto.  If not, see <http://www.gnu.org/licenses/>.

# ---- Import standard modules to the python path.
import numpy

from gwpy.timeseries import TimeSeriesDict
from gwpy.timeseries import TimeSeries
from gwpy.signal import filter_design
from .xfrequencyseries import XFrequencySeriesDict
from .xtimefrequencymap import XTimeFrequencyMapDict, XTimeFrequencyMap
from collections import OrderedDict

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['XTimeSeries']

class XTimeSeries(TimeSeriesDict):
    @classmethod
    def retrieve_data(cls, event_time, block_time,
                      channel_names, sample_frequency,
                      **kwargs):
        """Obtain data for either on source, off source, or injections.

        This uses the gwpy `TimeSeriesDict.get` method

        Parameters
        ----------
        event_time : (`float)
            trigger time of event to be processing

        block_time : (`int`)
            length of data to be processed

        channel_names (`list`) :
            required data channels.

        sample_frequency (`int`):
            sample rate of the data desired

        frame_types : `dict`, optional
            key channel and value frametype in whcih the channel is stored

        verbose : `bool`, optional
            print verbose output about NDS progress.

        Returns:

            `TimeSeriesDict` :
        """
        # Default is you do not want a timeslided dictionary
        # of timeseries
        time_slides = kwargs.pop('time_slides', {det : 0 for det in channel_names})
        if not isinstance(time_slides, dict):
            raise ValueError("time_slides must be supplied "
                              "in the form of a dict "
                              "with key channel and value integer "
                              "time to slide start and stop time")

        # If frametype is None then make dictionary
        # of key channel and value None
        frame_types = kwargs.pop('frame_types', {det : None for det in channel_names})
        if not isinstance(time_slides, dict):
            raise ValueError("frame_types must be supplied "
                              "in the form of a dict "
                              "with key channel and value frametype "
                              "assosciated with that channel")

        #----- Start and stop time for this event.
        start_time = event_time - block_time / 2;
        stop_time = event_time + block_time / 2;

        try:
            data = cls()
            for det in channel_names:
                data.append({det : TimeSeries.get(det,
                                           start_time + time_slides[det],
                                           stop_time + time_slides[det],
                                           frametype=frame_types[det],
                                           **kwargs
                                          )
                             })
        except:
            # Retrieve data useing the get classmethod 
            data = cls()
            for det in channel_names:
                data.append({det : TimeSeries.fetch_open_data(
                                                              det.split(':')[0],
                                           start_time + time_slides[det],
                                           stop_time + time_slides[det],
                                           sample_rate=sample_frequency,
                                           **kwargs
                                          )
                                 })

        for key, series in data.items():
            if series.sample_rate.decompose().value != sample_frequency:
                data[key] = series.resample(sample_frequency)

        # set one of the detectors to be the reference detecotr
        # for any future coherent combinations

        return XTimeSeries(data)


    @classmethod
    def generate_data(cls, event_time, block_time,
                      channel_names, sample_frequency,
                      verbose=False):
        """Obtain data for either on source, off source, or injections.

        This uses the gwpy `TimeSeriesDict.get` method

        Parameters
        ----------
        event_time : (`float)
            trigger time of event to be processing

        block_time : (`int`)
            length of data to be processed

        channel_names (`list`) :
            required data channels.

        sample_frequency (`int`):
            sample rate of the data desired

        verbose : `bool`, optional
            print verbose output about NDS progress.

        Returns:

            `TimeSeriesDict` :
        """
        #----- Start and stop time for this event.
        start_time = event_time - block_time / 2;
        stop_time = event_time + block_time / 2;

        # Retrieve data and then resample and set epoch
        data = {}
        for det in channel_names:
            data[det] = TimeSeries(numpy.random.normal(scale=.1,
                                        size=sample_frequency*block_time),
                                   sample_rate=sample_frequency)

        data = XTimeSeries(data)

        for (idet, iseries) in data.items():
            # set epoch of timeseries to the event_time
            iseries.t0 = start_time

        # set one of the detectors to be the reference detecotr
        # for any future coherent combinations

        return data

    @classmethod
    def generate_noise_from_file(cls, file_name, event_time, block_time,
                                 channel_names, sample_frequency,
                                 verbose=False):
        """Obtain data for either on source, off source, or injections.

        This uses the gwpy `TimeSeriesDict.get` method

        Parameters
        ----------
        event_time : (`float)
            trigger time of event to be processing

        block_time : (`int`)
            length of data to be processed

        channel_names (`list`) :
            required data channels.

        sample_frequency (`int`):
            sample rate of the data desired

        verbose : `bool`, optional
            print verbose output about NDS progress.

        Returns:

            `TimeSeriesDict` :
        """
        #----- Start and stop time for this event.
        start_time = event_time - block_time / 2;
        stop_time = event_time + block_time / 2;

        # Retrieve data and then resample and set epoch
        data = {}
        for det in channel_names:
            data[det] = TimeSeries.read(file_name,
                                        path='/noise/aLIGOZeroDetHighPower')

        data = XTimeSeries(data)

        for (idet, iseries) in data.items():
            # set epoch of timeseries to the event_time
            iseries.t0 = start_time

        # set one of the detectors to be the reference detecotr
        # for any future coherent combinations

        return data

    def asd(self, fftlength, **kwargs):
        """Obtain the asd of items in this dict.

        Parameters
        ----------
        fftlength : `dict`, `float`
            either a `dict` of (channel, `float`) pairs for key-wise
            asd calc, or a single float/int to compute as of all items.

        **kwargs
             other keyword arguments to pass to each item's asd
             method.
        """
        asds = XFrequencySeriesDict()
        if not isinstance(fftlength, dict):
            fftlength = dict((c, fftlength) for c in self)

        for key, fftlen in fftlength.items():
            asds[key] = self[key].asd(fftlength=fftlen,
                                      overlap=fftlen*0.5,
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
            dt = 1./asds[idet].dx
            whitened = iseries.whiten(fftlength=dt,
                                      overlap=dt*0.5,
                                      asd=asds[idet])
            whitened_timeseries.append(
                                       {idet : whitened}
                                       )

        return whitened_timeseries

    def highpass(self, frequency, corruption=2):
        """Design a high-pass filter.
        
        Parameters
        ----------
        fhigh : `float`
            high-pass corner frequency
        """
        for k, v in self.items():
            ts_hp = v.highpass(frequency)
            self[k] = ts_hp.crop(*ts_hp.span.contract(corruption))

        return self

    def spectrogram(self, fftlength):
        """Obtain the spectrograms of items in this dict.

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
            tfmaps[idet] = XTimeFrequencyMap(iseries.spectrogram2(
                                             fftlength=fftlength[idet],
                                             overlap=0.5*fftlength[idet],
                                             window='hann'))
        return tfmaps


    def fftgram(self, fftlength):
        """Obtain the spectrograms of items in this dict.

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
            tfmaps[idet] = XTimeFrequencyMap(iseries.fftgram(
                                             fftlength=fftlength[idet],
                                             overlap=0.5*fftlength[idet],
                                             window='hann'))
        return tfmaps


    def inject(self, injection_data, **kwargs):
        """Take an injection time series and inject into `XTimeSeries`

        This is essentially a wrapper arounf the very useful
        `gwpy.timeseries.TimeSeries.inject` method

        Parameters
        ----------
        injection_data : `XTimeSeries`,
            A `XtimeSeries` containing the same keys as
            the `XTimeSeries` you are injecting into
            i.e. if you are injecting into an
            `XTimeSeries` with 'H1', 'L1' and 'V1'
            data then your injeciton data
            should have the same keys

        **kwargs
             other keyword arguments to pass to each item's asd
             method.
        """
        injection_timeseries = XTimeSeries()
        for (idet, iseries) in self.items():
            injection_timeseries.append(
                {idet : iseries.inject(injection_data[idet])}
                )

        return injection_timeseries
