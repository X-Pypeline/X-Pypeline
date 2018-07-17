from xpipeline.core.xdetector import Detector
from .xwaveform import xmakewaveform
from gwpy.timeseries import TimeSeries
from xpipeline.core.xtimeseries import XTimeSeries

import numpy
import math


def xinjectsignal(family, start_time, block_time, channels, sample_rate,
                  phi, theta, psi=0, **kwargs):

    # ---- Check for optional arguments.
    catalogdirectory = kwargs.pop('catalogdirectory', '')
    parameters = kwargs.pop('parameters', [])

    # ----- Speed of light (m/s).
    c = 299792458

    # ----- Parse channel list and load info on corresponding detectors.
    injection_data = XTimeSeries()
    for det in channels:

        #----- Make timeseries data for given detector.
        injection_data[det] = TimeSeries(numpy.zeros(block_time * sample_rate),
                                         dx= 1./sample_rate,
                                         name=det)

        # create detector instance
        detector = Detector(det)

        #----- Time delay for incoming signal wrt center of Earth
        delay = detector.time_delay_from_earth_center_phi_theta(
                                       [phi], [theta])

        #----- Compute antenna response, including polarization angle.
        [Fp, Fc, Fb, F1, F2, FL] = detector_dict.compute_antenna_response(
                                       [phi], [theta], [psi])

        snippet_pad = block_time
        snippet_duration = 2*snippet_pad + 1

        #----- Peak Time of signal relative to start_time.
        peak_time = start_time + delay
        peak_time_for_waveform = snippet_pad + peak_time - math.floor(peak_time)

        #----- Waveform plus and cross polarizations, in wave frame.
        [t,hp,hc,hb] = xmakewaveform(family=family,
                                     parameters=parameters,
                                     T=snippet_duration,
                                     T0=peak_time_for_waveform[0],
                                     fs=sample_rate,
                                     catalogdirectory=catalogdirectory)

        # ---- Apply taper to start and end of snippet to avoid
        #      discontinuous turn-on/off of very long signals.
        # Only taper if there is an signal to taper
        if numpy.where(hp != 0)[0].size:
            hp.taper()
        if numpy.where(hc != 0)[0].size:
            hc.taper()
        if numpy.where(hb != 0)[0].size:
            hb.taper()

        #----- Reset t to GPS time.
        t = t + math.floor(peak_time) - snippet_pad

        #----- Sample indices.
        injection_samples = numpy.arange(0, t.size)

        #----- keep only samples which fall inside the desired time interval.
        k = numpy.where( (t>=start_time) & (t<start_time + block_time) )[0]
        injection_samples = injection_samples[k]
        t = t[k]
        hp = hp[k]
        hc = hc[k]
        hb = hb[k]

        #----- Combine with antenna response to give signal.
        injection_samples = (injection_samples +
                             numpy.round(
                             (math.floor(peak_time) - snippet_pad - start_time)
                             *sample_rate)
                            ).astype(int)

        injection_data[det][injection_samples] = Fp*hp + Fc*hc + Fb*hb
        injection_data[det].t0 = start_time
        injection_data[det].peak = peak_time
