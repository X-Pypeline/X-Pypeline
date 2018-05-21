.. _examples:

########################
xpipeline an explanation
########################

============
Introduction
============
The story of one “trial” of core X-Pipeline,

There is some timestamp considered as an “event time”

Some block of data from N detectors if obtain around this event time.

A time-frequency is made for a given fftlength, and based on time delays between detectors which is based on a potential sky location of the “source”, the tfmaps of N-1 detectors is phase shifted appropriately.

From these individual maps loud pixels are identified and then a “coherent” set of pixels is calculated from the overlap of these individually loud pixels.

a coherent energy is calculated and a threshold on loud coherent pixels is used,

These final loud enough coherent pixels are then turned from individual pixels into clusters using whatever clustering algorithm

these clusters then are assigned a variety of *likelihoods* these *likelihoods* for each cluster are based on the concept of translating the data in the dominant polarization frame i.e. plus and cross polarized time frequency maps instead of the standard time frequency maps considered above.

So the final important stat calculated is energy_of_cluster * likelihood_of_cluster for a given likelihood.

These likelihood rely largely on the ability to project the N stream of gravitational wave data in the relative antenna response pattern weighted time frequency space (so fplus fcross, fscalar, etc)

In order to calculate these antenna weighted tfmaps we need to re-weight the ASDs of the data streams appropriately. This re-weighted ASD is then multiplied with the time frequency map obtained normally above.

The Time-Frequency Map
----------------------

.. ipython::

    In [1]: from xpipeline.core.xtimeseries import XTimeSeries

    In [2]: data = XTimeSeries.retrieve_data(event_time=1135223174.0, block_time=256, sample_frequency=1024, channel_names=['H1:DCS-CALIB_STRAIN_C01','L1:DCS-CALIB_STRAIN_C01'], verbose=True)

    In [3]: asds = data.asd(1.0)

    In [4]: data.whiten(asds)

    In [5]: whitened_timeseries = data.whiten(asds)

    In [6]: print(whitened_timeseries)

The Dominant Polarization Frame
-------------------------------

Coherent Data
-------------

xpipeline likelihoods
---------------------
