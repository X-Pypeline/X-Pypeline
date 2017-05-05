#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np

from glue import datafind

from gwpy.timeseries import TimeSeries

def main(centerTime,duration,frameTypes,channelNames,detectors,rightascension,declination,FFTlength,sampleFrequency):

    #----- Start and stop time for this event.
    startTime = centerTime - duration / 2;
    stopTime = centerTime + duration / 2;

    # zip frameTypes and detectors, and channel names and detectors
    frameType   = dict(zip(detectors,frameTypes))
    channelName = dict(zip(detectors,channelNames))
    data        = dict()
    white_data  = dict()

    # Read in the data
    for iDet in detectors:
        connection = datafind.GWDataFindHTTPConnection()
        cache      = connection.find_frame_urls(iDet.strip('1'), frameType[iDet], startTime, stopTime, urltype='file')
        data[iDet] = TimeSeries.read(cache,channelName[iDet], format='gwf',start=startTime,end=stopTime)

    for (iDet,iSeries) in data.iteritems():
        # resample data
        if iSeries.sample_rate.decompose().value != sampleFrequency:
            iSeries = iSeries.resample(sampleFrequency)
        asd = iSeries.asd(FFTlength, FFTlength/2., method='median-mean')
        # Apply ASD to the data to whiten it
        whitened = iSeries.whiten(FFTlength, FFTlength/2., asd=asd) 
        white_data[iDet] = whitened.fft()

if __name__ == '__main__':
   main()
