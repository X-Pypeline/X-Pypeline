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

    In [2]: data = XTimeSeries.read('examples/example.gwf', channels=['H1:DCS-CALIB_STRAIN_C01','L1:DCS-CALIB_STRAIN_C01']) 

    In [3]: asds = data.asd(1.0)

    In [4]: data.whiten(asds)

    In [5]: whitened_timeseries = data.whiten(asds)

    In [6]: print(whitened_timeseries)

    In [7]: tfmaps = whitened_timeseries.spectrogram(1.0)

    In [8]: print(tfmaps)

The Dominant Polarization Frame
-------------------------------

.. ipython::

    In [1]: from xpipeline.core.xtimeseries import XTimeSeries

    In [2]: data = XTimeSeries.read('examples/example.gwf', channels=['H1:DCS-CALIB_STRAIN_C01','L1:DCS-CALIB_STRAIN_C01'])

    In [3]: asds = data.asd(1.0)

    In [4]: data.whiten(asds)

    In [5]: whitened_timeseries = data.whiten(asds)

    In [6]: tfmaps = whitened_timeseries.spectrogram(1.0)

    In [7]: from xpipeline.core.xdetector import compute_antenna_patterns

    In [8]: phi = 0.7728; theta = 1.4323

    In [9]: antenna_patterns = compute_antenna_patterns(['H1', 'L1'], phi, theta, antenna_patterns=['f_plus', 'f_cross', 'f_scalar'])

    In [10]: print(antenna_patterns)

    In [11]: projected_asds = asds.project_onto_antenna_patterns(antenna_patterns)

    In [12]: print(projected_asds)

    In [13]: projected_tfmaps = tfmaps.to_dominant_polarization_frame(projected_asds)

    in [14]: print(projected_tfmaps)

Coherent Data
-------------

xpipeline likelihoods
---------------------
