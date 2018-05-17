#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np

from gwpy.timeseries import TimeSeriesDict


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

        return XTimeSeries(data)


    def whiten(self, fftlength):
        """White this `TimeSeries` against its own ASD

            Parameters
            ----------
            fft_length : `float`
                number of seconds in single FFT
        """
        whitened_timeseries = XTimeSeries()
        asd_frequency_series = XTimeSeries()
        for (idet, iseries) in self.items():
            asd = iseries.asd(fftlength,
                              fftlength/2.,
                              method='lal_median_mean')
            whitened = iseries.whiten(fftlength,
                                      fftlength/2.,
                                      asd=asd)
            whitened_timeseries.append(
                                       {idet + '-whitened' : whitened}
                                       )
            asd_frequency_series.append(
                                       {idet + '-ASD' : asd}
                                       )

        return whitened_timeseries, asd_frequency_series


    def normalize_whitened_data(self, fftlength)
        return


    def shift_by_sky_location(self, phi, theta):
        return
