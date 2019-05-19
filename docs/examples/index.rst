.. _examples:

########################
xpipeline an explanation
########################

============
Introduction
============
Xpipeline is an unmodelled Gravitational Wave Burst analysis algorithm aimed at
using coherent detection statistics as well as coherent conistency checks to eliminate
surperfluous noise transients and extact real gravitational wave events from the data.

In order to build the statistics necessary to claim a detection
and to calculate the coherent statistics needed to eliminate excess back ground noise,
xpipeline conists of performing the same statistical "trials" over stretches of
data that couldhave a gravitational wave event and data which it is
either not believed or not possible for it to contain a graviational wave event.
In addition to these two stretches of data, an analysis is also performed on a
variety of simulated signals which is injected
into the data that is believed to contain a graviational wave in order that we
may understand what to expect from an actual gravitational wave
if it were to appear in the data during the stretch fo time believed to contain and gravitational wave.

Below we highlight what one of these "statistical trials" look like (either for the time
that does or does not contain a graviational wave). First, we use time from a known graviational wave event
to hightlight what the trial would look like for it, and then second we use a known noise transient
to see what it would look like for it.

In this breakdown, we imagine we are running over a stretch of data due to an
alert from a electromagnetic countrpart, i.e. an event that is believed to create
photometric and graviational wave signatures.

Through the alert, there is some timestamp considered as an “event time”, and around said
event time is some stretch of probable time that could contain the graviational wave counterpart to the
observation.

The data corresponding to the appropriate stretch of plausible time when
the gravitational wave could have passed through earth is obtained for N detectors
that were operating during that time.

The data from these N detectors are then FFT'ed and taken from time-amplitude
into time-frequency space. Since we have a plausible sky location of the source,
not only can we shrink the amount of data we search for a graviational wave over,
we can combine the data streams from the N detectors based on the known travel time of
graviational wave (i.e. the speed of light).

The tfmaps of N-1 detectors are time/phase shifted appropriately with the sky
location (locations if there is some error box).

From these individual maps loud time-frequency pixels are identified and
a “coherent” set of pixels is calculated from the overlap of these individually loud pixels.

After these coherent (overlapping) pixels are calculted a chosen "clustering" algorithm is
applied to the pizels that group near by loud pixels into possible graviational wave candidate events.

So essentially, our possible grviational wave events are simply a cluster of loud time frequency pixels upon
which numerous fancy statistics will be applied (these statistics often are called coherent statistics as
they are most useful when N streams of data can be projected combined together).

Specificlly, these clusters then are assigned a variety of *likelihoods* these *likelihoods* for each cluster are based on the concept of translating the data in the dominant polarization frame i.e. plus and cross polarized time frequency maps instead of the standard time frequency maps considered above.

So the final important stat calculated is energy_of_cluster * likelihood_of_cluster for a given likelihood.

These likelihood rely largely on the ability to project the N stream of gravitational wave data in the relative antenna response pattern weighted time frequency space (so fplus fcross, fscalar, etc)

In order to calculate these antenna weighted tfmaps we need to re-weight the ASDs of the data streams appropriately. This re-weighted ASD is then multiplied with the time frequency map obtained normally above.

So let us see what the above looks like if we had had an counterpart signal to the very
first graviational wave detection.


The Time-Frequency Map
----------------------
The following is all open data obtained via `LOSC <https://losc.ligo.org/>`_

Data from Hanford and Livingston is obtained, and then whitened (one can think of this as a normalization process).
The realty is that our detectors are more or less senstive to different frequencies and therefore in order
to detect excess noise we must first "whiten that data."

As I know the event time I have zoomed in close on the over all time frequnecy map
so that the signal is quite clear even without any fancy statistics. Nonetheless,
let us go through the whole process.

.. ipython::
    :okwarning:

    In [1]: from xpipeline.core.xtimeseries import XTimeSeries

    In [2]: data = XTimeSeries.read('examples/GW150914.gwf', channels=['H1:GDS-CALIB_STRAIN','L1:GDS-CALIB_STRAIN'])

    In [3]: asds = data.asd(1.0)

    In [4]: whitened_timeseries = data.whiten(asds)

    In [5]: fft_maps = whitened_timeseries.fftgram(1. /64)

    In [6]: energy_maps = fft_maps.abs()

    In [6]: gps = 1126259462.427

    In [7]: plot = energy_maps.plot(figsize=[ 12, 6])

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
In the above example you may notice that the data streams from both detectors do not seem to be
in sync. Well this is because we have not utilized the most important concept in this
unmodeled gravitational wave search analysis, the coherent combining of data streams, based upon
an a sky location value known ahead of time.

In this example we use a sky location chosen from the sky map assocaited with GW150914
to illustrate what a coherent anlaysis might look like if say you had a Gamma-Ray-Burst
or Supernova counterpart you were following up.

The way we can accomplish this is by either physically shifting the data of N-1 detectors
relative to a baseline detector some delta T amount or we can phase shift the pixels
of the timefrequencymap, here we physically shift the livingston data.

.. ipython::

    In [10]: from xpipeline.core.xdetector import Detector; import numpy

    In [11]: hanford = Detector('H1')

    In [12]: livingston = Detector('L1')

    In [13]: phi = -0.3801; theta = 2.7477 # Earth fixed coordinates

    In [14]: time_shift = numpy.round((livingston.time_delay_from_earth_center_phi_theta([phi], [theta]) - hanford.time_delay_from_earth_center_phi_theta([phi], [theta]))*data['H1:GDS-CALIB_STRAIN'].sample_rate)

    In [15]: whitened_timeseries['L1:GDS-CALIB_STRAIN'] = numpy.roll(whitened_timeseries['L1:GDS-CALIB_STRAIN'], -time_shift.astype(int))

    In [16]: fft_grams = whitened_timeseries.fftgram(1. /64)

    In [17]: energy_maps = fft_grams.abs()

    In [17]: plot = energy_maps.plot(figsize=[ 12, 6])

    In [18]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-ts-shifted.png
    In [19]: plot

    In [20]: coh_energy_maps = energy_maps.to_coherent()

    In [20]: plot = coh_energy_maps.plot(figsize=[ 12, 6])

    In [21]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-time-frequency-map-time-shifted-coherent.png
    In [22]: plot


Clustering Pixels
-----------------
There are a few ways to speed up the processing of the map. Many of the pixels
are not going to be significant, so we can threhold on what pixels we want
(say the loudest 1 percent of pixels) and then employ a method to group the pixels
together in what are referred to as `clusters`. These `clusters` become our possible
gravitational wave `triggers` on which we will later perform statistics on to determine
whether they originate from gravitational wave source or not.

In order to retain the visual key of a time frequency map pixels being grouped
and added together, but also perform algorithmic opertaions quickly we
utilize a subclass of the `scipy.sparse.csc_matrix` class which is designed
to efficiently perform operation on 2D matrices that are mostly zeroes
(which is what happens when we set 99 percent of the pixels to zero).

The algorithm employed to label the remaining pixels in our map into groups
is a nearest neighbor cpp wrapped algorithm called fastlabel.

.. ipython::

    In [35]: energy_map_zeroed = energy_maps.blackout_pixels(99)

    In [20]: plot = energy_map_zeroed.plot(figsize=[ 12, 6])

    In [21]: for ax in plot.axes:
       ....:     ax.set_xlim(gps - 0.15, gps + 0.05)
       ....:     ax.set_epoch(gps)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(20, 500)

    @savefig plot-sparse-ind-tfmaps.png
    In [22]: plot

    In [35]: coh_map = energy_map_zeroed.to_coherent()

    In [36]: tf_indices = coh_map.nonzero() # find the union of pixels from the sparse in coherent maps

    In [41]: tindex = {k : tf_indices[0]
       ....:           for k in energy_maps}

    In [42]: findex = {k : tf_indices[1]
       ....:           for k in energy_maps}

    In [43]: energy_map_zeroed = energy_maps.to_sparse(tindex, findex)

    In [44]: clusters = energy_map_zeroed.cluster() # cluster over the coherent pixels

    In [40]: print(clusters)


Now the we have labelled our remaining pixels (the non-zeroed out pixels), let's extract
some of the cluster properites of these clusters. i.e. how many pixels are in the cluster
the bounding box of the cluster (i.e. [[min-time, max-time], [min-freq, max-freq]] and the
sum of energy over the cluster.

Specifically, the cpp wrapped function `clusterproperities` outputs the following information

    * column 0: minimum time of cluster
    * column 1: weighted center time of cluster
    * column 2: maximum time of cluster
    * column 3: minimum frequency of cluster
    * column 4: weighted center frequency of cluster
    * column 5: maximum frequency of cluster
    * column 6: number of pixels in cluster
    * column 7: sum-over-cluster map values for each likelihood

Before we do this though, we must re-make our sparse maps.
for we have zeroed out some pixels in either map that are now part of our
clusters. i.e. some pixels may have been in the top 1 percent of one but not
all maps.

.. ipython::

    In [48]: min_time = clusters['min_time_of_cluster'].iloc[0]; max_time = clusters['max_time_of_cluster'].iloc[0]; weighted_center_time = clusters['weighted_center_time'].iloc[0]; min_freq = clusters['min_frequency_of_cluster'].iloc[0]; max_freq = clusters['max_frequency_of_cluster'].iloc[0];

    In [50]: plot = energy_map_zeroed.to_xtimefrequencymapdict().to_coherent().plot()

    In [51]: for ax in plot.axes:
       ....:     ax.set_xlim(min_time, max_time)
       ....:     ax.set_epoch(weighted_center_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(min_freq, max_freq)

    @savefig loudest-cluster-gw150914.png
    In [32]: plot

The Dominant Polarization Frame
-------------------------------
Now the we have a sky location assosciated with the event we can project every time-freqeuncy pixel
into the Dominant Polarization Frame (DPF). What this means is that we assume the GW has a plus and cross
polarization there is some orthoganal projection of the pixels onto the plus-cross plane for 2 or more detectors

.. ipython::

    In [13]: from xpipeline.core.xdetector import compute_antenna_patterns

    In [14]: phi = -0.3801; theta = 2.7477 # Earth fixed coordinates

    In [15]: antenna_patterns = compute_antenna_patterns(['H1', 'L1'], phi, theta, antenna_patterns=['f_plus', 'f_cross', 'f_scalar'])

    In [18]: projected_asds = asds.project_onto_antenna_patterns(antenna_patterns, to_dominant_polarization_frame=True)

    In [19]: projected_fftmaps = fft_grams.to_dominant_polarization_frame(projected_asds)

Now that we have projected each pixels onto the plus and cross phase + amplitude space
Let's see what it looks like if we simply take these projections and plot them.

.. ipython::

    In [19]: sparse_projected_fftmaps = {k : v.to_sparse(tindex, findex) for k, v in projected_fftmaps.items()}

    In [19]: plot = sparse_projected_fftmaps['f_plus'].to_xtimefrequencymapdict().to_coherent().abs().plot()

    In [51]: for ax in plot.axes:
       ....:     ax.set_xlim(min_time, max_time)
       ....:     ax.set_epoch(weighted_center_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(min_freq, max_freq)

    @savefig fplus-gw150914.png
    In [32]: plot

    In [20]: plot = sparse_projected_fftmaps['f_cross'].to_xtimefrequencymapdict().to_coherent().abs().plot()

    In [51]: for ax in plot.axes:
       ....:     ax.set_xlim(min_time, max_time)
       ....:     ax.set_epoch(weighted_center_time)
       ....:     ax.set_xlabel('Time [milliseconds]')
       ....:     ax.set_ylim(min_freq, max_freq)

    @savefig fcross-gw150914.png
    In [32]: plot

Likelihoods
-----------
So now that we have possible graviational wave candidates in the form
of clusters of loud pixels and the projected values of those pixels,
how do we "rank" these clusters in order or more or less likely to have orginiated from
graviational wave origins.

The Waveform
------------
In order to train these likelihoods so we can understand what values to expect from
gravitational wave clusters instead of random noise fluctations or `glitches` we must
inject a number of fake gravitational wave like signals.

This involves to steps, generating a gravitational-wave like waveform on the fly
and then injecting that signal into a stretch of data.

The parameters that go into xmakewaveform are the `family` of waveform, a set of parameters specific for that
waveform. In this case, the hrss is the quadrature sum of the RSS amplitudes of the plus and cross
polarizations, tau is the duration, f0 is the central frequency, alpha is
the chirp parameter, and delta is the phase at the peak of the envelope.

.. ipython::

    In [40]: from xpipeline.waveform import xwaveform

    In [41]: from gwpy.plot import Plot

    In [42]: t, hp, hc, hb = xwaveform.xmakewaveform(family='chirplet', parameters=[1e-22, 0.0033, 300.0, 0, 0, 1], T=513, T0=256.6161, fs=1024)

    In [43]: fig = Plot(hp, hc)

    In [43]: ax = fig.gca()

    In [44]: ax.set_epoch(256.6161)

    In [45]: ax.set_xlim([256.6161 - 0.05, 256.6161 + 0.05])

    @savefig chirplet.png
    In [46]: plot

Now let's say this is not an analytical waveform and instead an hplus and hcross
from say a supernova simulation. We can also handles that, tracked by `git-lfs`,
the waveforms folder of X-Pypeline repository houses a number of hdf5 files
full of pregenerated waveforms.

.. ipython::

    In [40]: from xpipeline.waveform import xwaveform

    In [41]: from gwpy.plot import Plot

    In [42]: t, hp, hc, hb = xwaveform.xmakewaveform(family='Morozova2018',
       ....:     parameters=[8, 'M10_LS220'],
       ....:     T=1, T0=0, fs=16384, catalogdirectory='../waveforms/xpipeline-waveforms')

    In [43]: fig = Plot(hp, hc)

    In [43]: ax = fig.gca()

    In [44]: ax.set_xlim([0, 1.0])

    @savefig supernova-M10_LS220.png
    In [45]: plot

The Injection
-------------

In a coherent search it is not enough to simply inject any old signal.
You must take in a set of sky coordinates and project an individual
signal with its antenna pattern (for example Fp*hp and Fc*hc)
just like we do for the data. Let us say we have the SN waveform from above
now we will assume this SN signal occurred at the same earth fixed coordinates
of GW150914 from above, but since this is a simulation let us imagine VIRGO was
on at the time too.

.. ipython::

    In [1]: from xpipeline.waveform import xinjectsignal

    In [3]: start_time = 1156609396.0; block_time = 256; channels = ['H1', 'L1']; sample_rate = 16384; injection_file_name ='examples/injection_Morozova.txt'; injection_number=2; catalogdirectory='../waveforms/xpipeline-waveforms';

    In [4]: [injection_data, gps_s, gps_ns, phi, theta, psi] = xinjectsignal.xinjectsignal_fromfile(start_time=start_time, block_time=block_time, channels=channels, injection_file_name=injection_file_name, injection_number=injection_number, sample_rate= sample_rate, catalogdirectory=catalogdirectory)

    In [6]: from gwpy.plot import Plot

    In [7]: fig = Plot()

    In [7]: ax = fig.gca()

    In [8]: for k, v in injection_data.items():
       ...:     ax.plot(v, label=k,)

    In [9]: ax.set_xlim([gps_s, gps_s + 2])

    In [10]: fig.legend()

    @savefig Morozova2018-h1-l1.png
    In [10]: plot
