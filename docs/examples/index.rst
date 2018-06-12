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

    In [7]: plot = tfmaps.plot(figsize=[ 12, 6])

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

    In [15]: whitened_timeseries['L1:GDS-CALIB_STRAIN'].shift(-time_shift[0]) # In place shift

    In [16]: tfmaps_ts_shifted = whitened_timeseries.spectrogram(1. /64)

    In [17]: plot = tfmaps_ts_shifted.plot(figsize=[ 12, 6])

    In [18]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-ts-shifted.png
    In [19]: plot

    In [20]: plot = tfmaps_ts_shifted.to_coherent().plot(figsize=[ 12, 6])

    In [21]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-time-shifted-coherent.png
    In [22]: plot

    In [23]: tfmaps['L1:GDS-CALIB_STRAIN'] = tfmaps['L1:GDS-CALIB_STRAIN'].phaseshift(time_shift[0]).abs()

    In [24]: plot = tfmaps.plot(figsize=[ 12, 6])

    In [25]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-phase-shifted.png
    In [26]: plot

    In [27]: plot = tfmaps.to_coherent().plot(figsize=[ 12, 6])

    In [28]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-coherent-phase-shifted.png
    In [26]: plot

The Dominant Polarization Frame
-------------------------------
Now the we have a sky location assosciated with the event we can proclustersum.clustersum_wrapperject every time-freqeuncy pixel
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

    In [20]: plot = projected_tfmaps['f_plus'].plot(figsize=[ 12, 6])

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


Clustering Pixels
-----------------
There are a few ways to speed up the processing of the map. Many of the pixels
are going to not be significant, so we can threhold on what pixels we want
(say the loudest 1 percent of pixels) and then employ a method to group the pixels
together in what are referred to as `clusters`. These `clusters` become our possible
gravitational wave `triggers` on which we evaluate the likelihoods described above

column 0: minimum time of cluster
column 1: weighted center time of cluster
column 2: maximum time of cluster
column 3: minimum frequency of cluster
column 4: weighted center frequency of cluster
column 5: maximum frequency of cluster
column 6: number of pixels in cluster
column 7-?: sum-over-cluster map values for each likelihood

.. ipython::

    In [34]: from xpipeline.cluster import nearestneighbor

    In [34]: import numpy

    In [35]: tfmaps_ts_shifted = tfmaps_ts_shifted.blackout_pixels(99)

    In [35]: coh_map = tfmaps_ts_shifted.to_coherent() 

    In [35]: pixels = numpy.argwhere(coh_map).T

    In [37]: coord_dim_array = coh_map.shape

    In [38]: npixels = pixels.shape[1]; connectivity = 8;

    In [39]: labelled_map = nearestneighbor.fastlabel_wrapper(pixels + 1, coord_dim_array, connectivity, npixels).astype(int)

    In [40]: print(labelled_map)

Now the we have labelled are remaining pixels (the non-zeroed out pixels), let's extract
some fo the cluster properites of these clusters. i.e. how many piels are in the clsuter
the bounding box of the cluster (i.e. [[min-time, max-time], [min-freq, max-freq]] and the
sum of energy over the cluster.

Specifically the function `clusterproperities` outputs the following information

column 0: minimum time of cluster
column 1: weighted center time of cluster
column 2: maximum time of cluster
column 3: minimum frequency of cluster
column 4: weighted center frequency of cluster
column 5: maximum frequency of cluster
column 6: number of pixels in cluster
column 7-?: sum-over-cluster map values for each likelihood

.. ipython::

    In [41]: from xpipeline.cluster import clusterproperties

    In [41]: from gwpy.table import EventTable

    In [41]: total_energy = coh_map[pixels[0,:], pixels[1,:]] 

    In [43]: dim_array = numpy.array([total_energy.shape[0], 1, 2.0])

    In [42]: cluster_array = clusterproperties.clusterproperities_wrapper(dim_array, labelled_map, total_energy, pixels[0,:] + 1, pixels[1,:] + 1).T

    In [43]: cluster_array[:,0:3] = cluster_array[:,0:3]  * coh_map.dx + coh_map.t0

    In [44]: cluster_array[:,3:6] = cluster_array[:,3:6] * coh_map.dy + coh_map.y0

    In [67]: clusters = EventTable(cluster_array,
       ....:                       names=['min_time_of_cluster',
       ....:                              'weighted_center_time', 'max_time_of_cluster',
       ....:                              'min_frequency_of_cluster',
       ....:                              'weighted_center_frequency',
       ....:                              'max_frequency_of_cluster',
       ....:                              'number_of_pixels', 'energy_of_cluster'])

    In [47]: print(clusters)

    In [47]: loudest_cluster_idx = clusters['energy_of_cluster'].argmax()

    In [48]: min_time = clusters['min_time_of_cluster'][loudest_cluster_idx]; max_time = clusters['max_time_of_cluster'][loudest_cluster_idx]; weighted_center_time = clusters['weighted_center_time'][loudest_cluster_idx]; min_freq = clusters['min_frequency_of_cluster'][loudest_cluster_idx]; max_freq = clusters['max_frequency_of_cluster'][loudest_cluster_idx];

    In [50]: plot = coh_map.plot()

    In [51]: for ax in plot.axes:
       ....:     ax.set_xlim(min_time, max_time)
       ....:     ax.set_epoch(weighted_center_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(min_freq, max_freq)

    @savefig loudest-cluster-gw150914.png
    In [32]: plot


Alright, we now have a labelling of all pixels into clusters and likelihood maps.
So, let us calculated the likelihood of the clusters

.. ipython::

    In [41]: from xpipeline.cluster import clustersum

    In [42]: clustersum.clustersum_wrapper

    In [43]: likelihood = clustersum.clustersum_wrapper(labelled_map, likelihood_map_standard[pixels[0,:], pixels[1,:]])


The Waveform
------------
In order to train these likelihoods so we can understand what values to expect from
gravitational wave clusters instead of random noise fluctations or `glitches` we must
inject a number of fake gravitational wave like signals.

This involves to steps, generating a gravitational-wave like waveform on the fly
and then injecting that signal into a stretch of data.

The parametrs that go into xmakewaveform are the `family` of waveform, a set of parameters specific for that
waveform. In this case, the hrss is the quadrature sum of the RSS amplitudes of the plus and cross
polarizations, tau is the duration, f0 is the central frequency, alpha is
the chirp parameter, and delta is the phase at the peak of the envelope.

.. ipython::

    In [40]: from xpipeline.waveform import xwaveform

    In [41]: from gwpy.plotter import TimeSeriesPlot

    In [42]: t, hp, hc, hb = xwaveform.xmakewaveform(family='chirplet', parameters=[1e-22, 0.0033, 300.0, 0, 0, 1], T=513, T0=256.6161, fs=1024)

    In [43]: plot = TimeSeriesPlot(hp, hc)

    In [44]: plot.set_epoch(256.6161)

    In [45]: plot.set_xlim([256.6161 - 0.05, 256.6161 + 0.05])

    @savefig chirplet.png
    In [46]: plot

Now let's say this is not an analytical waveform and instead an hplus and hcross
from say a supernova simulation. We can also handles that, tracked by `git-lfs`,
the waveforms folder of X-Pypeline repository houses a number of hdf5 files
full of pregenerated waveforms.

.. ipython::

    In [40]: from xpipeline.waveform import xwaveform

    In [41]: from gwpy.plotter import TimeSeriesPlot

    In [42]: t, hp, hc, hb = xwaveform.xmakewaveform(family='o1snews',
       ....:     parameters=[1e-21, 1e-21, 'R4E1FC_L_theta2.094_phi2.094'],
       ....:     T=1, T0=0, fs=16384, catalogdirectory='../waveforms/')

    In [43]: plot = TimeSeriesPlot(hp, hc)

    In [44]: plot.set_xlim([0, 0.1])

    @savefig supernova-R4E1FC_L_theta2.094_phi2.094.png
    In [45]: plot


The Injection
-------------

In a coherent search it is not enough to simply inject any old signal.
You must take in a set of sky coordinates and project an individual
signal with its antenna pattern (for example Fp*hp and Fc*hc)
just like we do for the data.

.. ipython::

    In [1]: from xpipeline.waveform import xinjectsignal

    In [2]: start_time = 1156609396.0; block_time = 256; channels = ['H1', 'L1', 'V1']; sample_rate = 1024; injection_file_name ='examples/injection_sgc300.txt'; injection_number=0; catalogdirectory='';

    In [3]: [injection_data, gps_s, gps_ns, phi, theta, psi] = xinjectsignal.xinjectsignal(start_time=start_time, block_time=block_time, channels=channels, injection_file_name=injection_file_name, injection_number=injection_number, sample_rate= sample_rate, catalogdirectory=catalogdirectory)

    In [4]: print(gps_s, gps_ns, phi, theta, psi)

    In [7]: peak_time = injection_data['H1'].peak

    In [5]: for det, series in injection_data.items():
       ...:     injection_data[det] = series * 4.87

    In [5]: plot = injection_data.plot()

    In [6]: plot.add_legend()

    In [8]: plot.set_epoch(peak_time)

    In [9]: plot.set_xlim([peak_time - 0.1, peak_time + 0.1])

    @savefig chirplet-h1-l1-v1.png
    In [10]: plot

Now let's inject this into some data, we could use real data but let's just generate
some data and scale it to an amplitude where we would expect this waveform to show up.

.. ipython::

    In [11]: event_time = 1156609524; block_time = 256; channel_names = ['H1', 'L1', 'V1']; sample_frequency = 1024

    In [12]: data = XTimeSeries.generate_data(event_time=event_time,
       ....:                                  block_time=block_time,
       ....:                                  channel_names=channel_names,
       ....:                                  sample_frequency=sample_frequency)
       ....:

    In [13]: for det, series in data.items():
       ....:     data[det] = series * 1e-21

    In [14]: injection_series = data.inject(injection_data=injection_data)

    In [15]: injection_series.plot()

    In [16]: plot = injection_series.plot()

    In [17]: plot.add_legend()

    @savefig chirplet-h1-l1-v1-in-data.png
    In [17]: plot

Now you can see where the injection went in terms of the entire length of data
we are analyzing (a 256 second block) but let us zoom in a bit.

.. ipython::

    In [18]: plot.set_epoch(peak_time)

    In [19]: plot.set_xlim([peak_time - 0.1, peak_time + 0.1])

    @savefig chirplet-h1-l1-v1-in-data-zoom.png
    In [20]: plot

You will notice that just like int he case where we read in the data surrounding GW150914
we now has a variable TimeSeries that is bascially the same as above, except it has
an injected signal in there. Well let us look at what the likelihoods look like for this waveform


.. ipython::

    In [3]: asds = injection_series.asd(1.0)

    In [4]: whitened_timeseries = injection_series.whiten(asds)

    In [5]: tfmaps = whitened_timeseries.spectrogram(1. /64)

    In [7]: plot = tfmaps.plot(figsize=[12, 6])

    In [8]: for ax in plot.axes:
       ...:     ax.set_xlim(peak_time - 0.05, peak_time + 0.05)
       ...:     ax.set_epoch(peak_time)
       ...:     ax.set_xlabel('Time [milliseconds]')
       ...:     ax.set_ylim(20, 500)
       ...:

    @savefig chirplet-time-frequency-map.png
    In [9]: plot

    In [10]: from xpipeline.core.xdetector import Detector

    In [11]: hanford = Detector('H1')

    In [12]: livingston = Detector('L1')

    In [13]: virgo = Detector('V1')

    In [14]: time_shift_livingston = livingston.time_delay_from_earth_center_phi_theta([phi], [theta]) - hanford.time_delay_from_earth_center_phi_theta([phi], [theta])

    In [14]: time_shift_virgo = virgo.time_delay_from_earth_center_phi_theta([phi], [theta]) - hanford.time_delay_from_earth_center_phi_theta([phi], [theta])

    In [15]: whitened_timeseries['L1'].shift(-time_shift_livingston[0]) # In place shift

    In [15]: whitened_timeseries['V1'].shift(-time_shift_virgo[0]) # In place shift

    In [16]: plot = whitened_timeseries.plot()

    In [17]: plot.add_legend()

    In [19]: plot.set_xlim([peak_time - 0.05, peak_time + 0.05])

    @savefig plot-chirplet-wts-shifted.png
    In [22]: plot

    In [16]: tfmaps_ts_shifted = whitened_timeseries.spectrogram(1. /64)

    In [17]: plot = tfmaps_ts_shifted.plot(figsize=[ 12, 6])

    In [18]: for ax in plot.axes:
       ....:     ax.set_xlim(peak_time - 0.05, peak_time + 0.05)
       ....:     ax.set_epoch(peak_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-chirplet-time-frequency-map-ts-shifted.png
    In [19]: plot

    In [13]: from xpipeline.core.xdetector import compute_antenna_patterns

    In [14]: import numpy as np

    In [15]: antenna_patterns = compute_antenna_patterns(['H1', 'L1', 'V1'],
       ....:     phi, theta, antenna_patterns=['f_plus', 'f_cross', 'f_scalar'])

    In [16]: frequencies = np.in1d(asds['L1'].xindex.to_value(), tfmaps_ts_shifted['L1'].yindex.to_value())

    In [17]: sliced_asds = asds.slice_frequencies(frequencies)

    In [18]: projected_asds = sliced_asds.project_onto_antenna_patterns(antenna_patterns, to_dominant_polarization_frame=True)

    In [19]: projected_tfmaps = tfmaps_ts_shifted.to_dominant_polarization_frame(projected_asds)

    In [20]: plot = projected_tfmaps['f_plus'].plot(figsize=[12, 6])

    In [21]: for ax in plot.axes:
       ....:     ax.set_xlim(peak_time - 0.05, peak_time + 0.05)
       ....:     ax.set_epoch(peak_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-chirplet-time-frequency-map-dpf-plus.png
    In [22]: plot

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
       ....:     ax.set_xlim(peak_time - 0.05, peak_time + 0.05)
       ....:     ax.set_epoch(peak_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-chirplet-time-frequency-map-likelihood-maps.png
    In [32]: plot
