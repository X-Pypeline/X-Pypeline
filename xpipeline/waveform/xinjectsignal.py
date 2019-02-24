from xpipeline.core.xdetector import Detector
from .xwaveform import xmakewaveform
from gwpy.timeseries import TimeSeries
from xpipeline.core.xtimeseries import XTimeSeries
from gwpy.table import Table
from astropy.table import hstack

import numpy
import math

def xinjectsignal(family, start_time, block_time, detectors, sample_rate,
                  phi, theta, psi=0, **kwargs):
    """Create simulated signals using a file defining the parameters

    Parameters:
        family (str):
            waveform family name to be generated

        start_time (int):
            Start time of the data being analysed.

        block_time (int):
            Duration (s) of the data being analysed.

        channel (str)
            channel[j] holds the name of the
            jth detector channel being analysed.  The detector name is
            inferred from the first character of the channel name, and
            must be one of the names recognized by LoadDetectorData.

        sample_rate (int):
            Sample rates corresponding to channels.

        phi (float):
            Earth fixed coordinate

        theta (float):
            Earth fixed coordinate

        **kwargs:
            catalog_directory (optional, str):
                Path to directory containing pregenerated waveforms

    Returns:
        injection_data:
            `XTimeSeries`
        
    Notes:

        The coordinate system is explained in ComputeAntennaResponse. It is
        assumed that the sky positions in the injection file are supplied in this
        coordinate system.  For glitch injections the reported times and angles
        are those for the injection into the first detector.

        Any portion of any signal which lies outside the interval
        [start_time, start_time+block_time] will be clipped. Otherwise, XINJECTSIGNAL
        will inject any signal of duration up to block_time without clipping.
        The edges of long-duration signals will be tapered with
        a hann window to avoid turn-on/off transients.
    """
    # ---- Check for optional arguments.
    catalogdirectory = kwargs.pop('catalogdirectory', '')
    parameters = kwargs.pop('parameters', [])

    # ----- Speed of light (m/s).
    c = 299792458

    # ----- Parse channel list and load info on corresponding detectors.
    injection_data = XTimeSeries()
    for det in detectors:

        #----- Make timeseries data for given detector.
        injection_data[det] = TimeSeries(numpy.zeros(int(block_time * sample_rate)),
                                         dx= 1./sample_rate,
                                         name=det)

        # create detector instance
        detector = Detector(det)

        #----- Time delay for incoming signal wrt center of Earth
        delay = detector.time_delay_from_earth_center_phi_theta(
                                       [phi], [theta])

        #----- Compute antenna response, including polarization angle.
        [Fp, Fc, Fb, F1, F2, FL] = detector.compute_antenna_response(
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
            hp = hp.taper('right')
        if numpy.where(hc != 0)[0].size:
            hc = hc.taper('right')
        if numpy.where(hb != 0)[0].size:
            hb = hb.taper('right')

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

    return injection_data

def xinjectsignal_fromfile(start_time, block_time, channels, sample_rate,
                  injection_file_name, rescale_by_antenna_response=False,
                  injection_number=0, **kwargs):
    """Create simulated signals using a file defining the parameters

    Parameters:
        start_time (int):
            Start time of the data being analysed.

        block_time (int):
            Duration (s) of the data being analysed.

        channel (str)
            channel[j] holds the name of the
            jth detector channel being analysed.  The detector name is
            inferred from the first character of the channel name, and
            must be one of the names recognized by LoadDetectorData.

        sample_rate (int):
            Sample rates corresponding to channels.

        injection_file_name (str):
            Name of file specifying simulated signals to be
            injected.

        injection_number
            Array (positive integer).  Optional.  If given then only
            the injections specified on rows "injection_number" of the
            injection file are injected.  This over-rides the default
            behaviour, which is to inject all injections in the file
            (including those outside the interval [start_time,start_time+block_time]!).

        **kwargs:
            catalog_directory (optional, str):
                Path to directory containing pregenerated waveforms

    Returns:
        injection_data:
            `XTimeSeries`
        
        gps_s:
            Vector of peak time (GPS seconds, integer part) of each
            injection at the center of the Earth.

        gps_ns:
            Vector of peak time (GPS seconds, nanosecond part) of each
            injection at the center of the Earth.

        phi:
            Vector of azimuthal sky coordinate of each injection
            (radians, Earth-fixed coordinates).

        theta:
            Vector of polar sky coordinate of each injection

        psi:
            Vector of polarization angle of each injection
            (radians Earth-fixed coordinates).

    Notes:

        The coordinate system is explained in ComputeAntennaResponse. It is
        assumed that the sky positions in the injection file are supplied in this
        coordinate system.  For glitch injections the reported times and angles
        are those for the injection into the first detector.

        Any portion of any signal which lies outside the interval
        [start_time, start_time+block_time] will be clipped. Otherwise, XINJECTSIGNAL
        will inject any signal of duration up to block_time without clipping.
        The edges of long-duration signals will be tapered with
        a hann window to avoid turn-on/off transients.
    """
    ifos = [ifo.split(':')[0] for ifo in channels]

    injection_number = int(injection_number)

    # ---- Check for optional arguments.
    do_not_inject = kwargs.pop('do_not_inject', 0)
    catalogdirectory = kwargs.pop('catalogdirectory', '')

    # ----- Speed of light (m/s).
    c = 299792458

    # ----- Parse channel list and load info on corresponding detectors.
    detector_dict = {}
    for det in ifos:
        detector_dict[det] = Detector(det)

    # We must check if there are detector specific injection
    # parameters
    injection_file_parameters = Table.read(injection_file_name, format='ascii',)
    _injection_file_columns = ['gps_s','gps_ns', 'phi', 'theta', 'psi',
                               'gwb_type', 'gwb_params']
    if len(injection_file_parameters.keys()) > 7:
        # ----- Create appropriate columns for reading the injection
        #       file based on detectors supplied above
         
        injection_file_columns = []
        for det in ifos:
            for columns in _injection_file_columns:
                injection_file_columns.append(columns + '_' + det)
        #----- Read injection file and extract parameters of injections
        injection_file_parameters = Table.read(injection_file_name, format='ascii',
                                               names=injection_file_columns)
    else:
        # There is a single set of parameters for ever detector
        injection_file_parameters = Table()
        for det in ifos:
            injection_file_columns = []
            for columns in _injection_file_columns:
                injection_file_columns.append(columns + '_' + det)
            injection_file_parameters = hstack((injection_file_parameters, Table.read(injection_file_name, format='ascii',names=injection_file_columns)))


    #----- If specific injections are specified then keep only those injections.
    injection_file_parameters = injection_file_parameters[injection_number]

    #----- Compute sum over detectors of Fp^2 for each injection, if desired.
    """
    if (rescale_by_antenna_response)
        totalInjectedPower = zeros(nInjection,1)
        #----- Loop over simulated signals and record peak time, sky angles.
        for jInjection=1:nInjection
            #----- Use angles for first detector.
            phi    = str2num(currentParameters{3}{jInjection})
            theta  = str2num(currentParameters{4}{jInjection})
            psi    = str2num(currentParameters{5}{jInjection})
            #----- Compute Fp for each detector.
            for jDetector=1:nDetector
                Fp = ComputeAntennaResponse(phi,theta,psi,detector{jDetector}.d)
                totalInjectedPower(jInjection) = ...
                    totalInjectedPower(jInjection) + Fp^2
            end
        end
    end
    """

    #----- Loop over detectors and construct simulated signals.
    injection_data = {}
    for channel, det in zip(channels, ifos):
        #----- Make timeseries data for given detector.
        injection_data[channel] = TimeSeries(numpy.zeros(int(block_time * sample_rate)),
                                         dx= 1./sample_rate,
                                         name=det)
        #---------- Extract parameters for current injection.
        gps_s = injection_file_parameters['gps_s_' + det]
        gps_ns = injection_file_parameters['gps_ns_' + det]
        phi = injection_file_parameters['phi_' + det]
        theta = injection_file_parameters['theta_' + det]
        psi = injection_file_parameters['psi_' + det]
        gwb_type = injection_file_parameters['gwb_type_' + det]
        gwb_params = injection_file_parameters['gwb_params_' + det]


        #----- Time delay for incoming signal wrt center of Earth
        delay = detector_dict[det].time_delay_from_earth_center_phi_theta(
                                       [phi], [theta])

        #----- Compute antenna response, including polarization angle.
        [Fp, Fc, Fb, F1, F2, FL] = detector_dict[det].compute_antenna_response(
                                       [phi], [theta], [psi])

        """
        #----- Rescale amplitude by sum_over_detectors of Fp^2, if desired:
        if (rescale_by_antenna_response):
            Fp = Fp / sqrt(totalInjectedPower(jInjection))
            Fc = Fc / sqrt(totalInjectedPower(jInjection))
        end
        """
        #----- Make injection data in a several-second snippet.  This will
        #      work for any signal with duration less than snippet_pad for
        #      longer signals some clipping may occur.  Clipping will
        #      definitely occur for injections longer than snippet_duration.
        snippet_pad = block_time
        snippet_duration = 2*snippet_pad + 1
        #----- Peak Time of signal relative to start_time.
        peak_time = gps_s + 1e-9 * gps_ns + delay

        peak_time_for_waveform = snippet_pad + peak_time - math.floor(peak_time)
        #----- Waveform plus and cross polarizations, in wave frame.
        [t,hp,hc,hb] = xmakewaveform(family=gwb_type,
                                     parameters=gwb_params,
                                     T=snippet_duration,
                                     T0=peak_time_for_waveform[0],
                                     fs=sample_rate,
                                     catalogdirectory=catalogdirectory)

        # ---- Apply taper to start and end of snippet to avoid
        #      discontinuous turn-on/off of very long signals.
        # Only taper if there is an signal to taper
        #if numpy.where(hp != 0)[0].size:
        #    hp.taper()
        #if numpy.where(hc != 0)[0].size:
        #    hc.taper()
        #if numpy.where(hb != 0)[0].size:
        #    hb.taper()

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

        injection_data[channel][injection_samples] = Fp*hp + Fc*hc + Fb*hb
        injection_data[channel].t0 = start_time

    #----- Record for output the peak time and sky angles for each of the
    #      injections.  For glitches record only the parameters for the first
    #      detector.
    #----- Loop over simulated signals and record peak time, sky angles.
    #---------- Extract parameters for current injection.
    gps_s = injection_file_parameters['gps_s_' + ifos[0]]
    gps_ns = injection_file_parameters['gps_ns_' + ifos[0]]
    phi = injection_file_parameters['phi_' + ifos[0]]
    theta = injection_file_parameters['theta_' + ifos[0]]
    psi = injection_file_parameters['psi_' + ifos[0]]

    #----- Done
    return XTimeSeries(injection_data), gps_s, gps_ns, phi, theta, psi
