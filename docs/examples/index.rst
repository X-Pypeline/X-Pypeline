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
The following is all open data obtained via `LOSC <https://losc.ligo.org/>`_

.. ipython::

    In [1]: from xpipeline.core.xtimeseries import XTimeSeries

    In [2]: data = XTimeSeries.read('examples/GW150914.gwf', channels=['H1:GDS-CALIB_STRAIN','L1:GDS-CALIB_STRAIN'])

    In [3]: asds = data.asd(1.0)

    In [4]: whitened_timeseries = data.whiten(asds)

    In [5]: tfmaps = whitened_timeseries.spectrogram(1. /64)

    In [6]: gps = 1126259462.427

    In [7]: plot = tfmaps.plot(figsize=[8, 4])

    In [8]: for ax in plot.axes:
       ...:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ...:     ax.set_epoch(gps)
       ...:     ax.set_xlabel('Time [milliseconds]')
       ...:     ax.set_ylim(20, 500)
       ...:

    @savefig plot-time-frequency-map.png
    In [9]: plot

The Coherent Time-Frequency Map
-------------------------------
In this example we use a sky location chosen from the sky map assocaited with GW150914
to illustrate what a coherent anlaysis might look like if say you had a Gamma-Ray-Burst
or Supernova counterpart you were following up.

The way we can accomplish this is by either physically shifting the data of N-1 detectors
relative to a baseline detector some delta T amount or we can phase shift the pixels
of the timefrequencymap, here we do each of these methods.

.. ipython::

    In [10]: from xpipeline.core.xdetector import Detector

    In [11]: hanford = Detector('H1')

    In [12]: livingston = Detector('L1')

    In [13]: phi = -0.3801; theta = 2.7477 # Earth fixed coordinates

    In [14]: time_shift = livingston.time_delay_from_earth_center_phi_theta([phi], [theta]) - hanford.time_delay_from_earth_center_phi_theta([phi], [theta])

    In [15]: whitened_timeseries['L1:GDS-CALIB_STRAIN'].shift(time_shift[0]) # In place shift

    In [16]: tfmaps_ts_shifted = whitened_timeseries.spectrogram(1. /64)

    In [17]: plot = tfmaps_ts_shifted.plot(figsize=[8, 4])

    In [18]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-ts-shifted.png
    In [19]: plot

    In [20]: plot = tfmaps_ts_shifted.to_coherent().plot(figsize=[8, 4])

    In [21]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-time-shifted-coherent.png
    In [22]: plot

    In [23]: tfmaps['L1:GDS-CALIB_STRAIN'] = tfmaps['L1:GDS-CALIB_STRAIN'].phaseshift(time_shift[0]).abs()

    In [24]: plot = tfmaps.plot(figsize=[8, 4])

    In [25]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-phase-shifted.png
    In [26]: plot

    In [27]: plot = tfmaps.to_coherent().plot(figsize=[8, 4])

    In [28]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-coherent-phase-shifted.png
    In [26]: plot

The Dominant Polarization Frame
-------------------------------
Now the we have a sky location assosciated with the event we can project every time-freqeuncy pixel
into the Dominant Polarization Frame (DPF). What this means is the is we assume the GW has a plus and cross
polarization there is some orthoganal projection of the pixels onto the plus-cross plane for 2 or more detectors

.. ipython::

    In [13]: from xpipeline.core.xdetector import compute_antenna_patterns

    In [14]: import numpy as np

    In [14]: phi = -0.3801; theta = 2.7477 # Earth fixed coordinates

    In [15]: antenna_patterns = compute_antenna_patterns(['H1', 'L1'], phi, theta, antenna_patterns=['f_plus', 'f_cross', 'f_scalar'])

    In [16]: frequencies = np.in1d(asds['L1:GDS-CALIB_STRAIN'].xindex.to_value(),tfmaps['L1:GDS-CALIB_STRAIN'].yindex.to_value())

    In [17]: sliced_asds = asds.slice_frequencies(frequencies) 

    In [18]: projected_asds = sliced_asds.project_onto_antenna_patterns(antenna_patterns, to_dominant_polarization_frame=True)

    In [19]: projected_tfmaps = tfmaps.to_dominant_polarization_frame(projected_asds)

    In [20]: plot = projected_tfmaps['f_plus'].plot(figsize=[8, 4])

    In [21]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-dpf-plus.png
    In [22]: plot

xpipeline likelihoods
---------------------
Now we have a basis to determine whether or not a particular cluster of pixels
can be considered likely was a gravitational wave

A gravitational wave not only should be coherent between the multiple data streams
but if it originated from a certain part of the sky the projection of the cluster onto
the plus and cross polarization plane (i.e. `projected_tfmaps` should also be large.

.. ipython::

    In [21]: from xpipeline.core.xlikelihood import XLikelihood

    In [22]: mpp = projected_asds['f_plus'].to_m_ab()

    In [23]: mcc = projected_asds['f_cross'].to_m_ab()

    In [24]: wfptimefrequencymap = projected_tfmaps['f_plus'].to_coherent()

    In [25]: wfctimefrequencymap = projected_tfmaps['f_cross'].to_coherent()

    In [26]: likelihood_map_standard = XLikelihood.standard(mpp, mcc, wfptimefrequencymap, wfctimefrequencymap)

    In [27]: likelihood_map_circenergy = XLikelihood.circenergy(mpp, mcc, wfptimefrequencymap, wfctimefrequencymap)

    In [28]: likelihood_map_circinc = XLikelihood.circinc(tfmaps, mpp, mcc, projected_asds)

    In [29]: likelihood_map_circnullinc = XLikelihood.circnullinc(tfmaps, mpp, mcc, projected_asds)

    In [30]: likelihood_map_circnullenergy = XLikelihood.circnullenergy(mpp, mcc, wfptimefrequencymap, wfctimefrequencymap)

    In [31]: plot = likelihood_map_standard.plot(figsize=(12,8), label='standard')

    In [32]: plot.add_spectrogram(likelihood_map_circinc, newax=True, label='circinc')

    In [33]: plot.add_spectrogram(likelihood_map_circnullenergy, newax=True, label='circnullenergy')

    In [34]: plot.add_spectrogram(likelihood_map_circnullinc, newax=True, label='circnullinc')

    In [35]: plot.add_spectrogram(likelihood_map_circenergy, newax=True, label='circenergy')

    In [31]: for ax in plot.axes:
       ....:     plot.add_colorbar(ax=ax)
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-likelihood-maps.png
    In [32]: plot
