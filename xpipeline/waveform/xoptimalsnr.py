import warnings
from scipy.interpolate import interp1d
import numpy


def xoptimalsnr(h, t0, fs, S=None, F0=None, dF=None, Fmin=None, Fmax=None):
    """Compute SNR and other properties of a signal in noise.

    XOPTIMALSNR - Compute the SNR and time-frequency measures of a waveform
    in a specified noise background.  Simpler noise-independent measures of
    the wave amplitude are also provided.

       [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = ...
           xoptimalsnr(h,t0,fs,S,F0,dF,Fmin,Fmax)

       h      Array.  Waveform timeseries data.  Each column holds the
              timeseries for one of the GW polarizations.  (There may be any
              number of polarisations, for the non-GR case.  The order of
              polarizations does not matter.)  The timeseries duration should
              be a power of 2.
       t0     Scalar.  Time at which the first waveform data point h(1,:) is
              sampled.
       fs     Scalar.  Sampling rate (Hz) of waveform data.  This should be a
              power of 2.
       S      Vector (optional).  Noise background one-sided POWER (not
              amplitude) spectrum.  If supplied, it must cover at least the
              range [Fmin, Fmax], and must be sampled in ascending
              order of frequency.  For noise-independent signal measures use
              S=f0=df=[].  In this case the constant spectrum S(f)=2 will be
              used (corresponding to a two-sided noise spectrum of unity).
       F0     Scalar or vector (optional).  Scalar: Frequency at which S(1) is
              sampled.  Vector: frequencies at which S is sampled.
       dF     Scalar (optional).  Frequency spacing of the noise spectrum S.
              If F0 is a vector then dF = [] must be used.
       Fmin   Scalar (optional).  Minimum frequency (Hz) to include in
              frequency-domain calculations.  Defaults to 0 if no noise
              spectrum specified, or minimum frequency of noise spectrum.
      Fmax   Scalar (optional).  Maximum frequency (Hz) to include in
              frequency-domain calculations. Defaults to Nyquist if no noise
              spectrum specified, or maximum frequency of noise spectrum.

     Computations are done in the time domain (TD) and frequency
     domain (FD) using the energy distributions
          p_TD = h(:,1).^2 + h(:,2).^2
          p_FD = 2(|FFT(h(:,1))|.^2 + |FFT(h(:,2))|.^2)  # for f>=0
     With these conventions, the output is
       SNR    Signal-to-noise ratio of h in the given noise background,
              defined as
                SNR = (2 \int_Fmin^Fmax df p_FD./S).^0.5
       h_rss  The root-sum-square amplitude (Hz^-0.5) of the waveform h:
                h_rss = \int_Fmin^Fmax df p_FD
       h_peak The maximum absolute amplitude of the waveform h:
                h_peak = max(p_TD).^0.5
       Fchar  Characteristic frequency (Hz) of h in the given noise
              background, defined as
                        \int_Fmin^Fmax df f p_FD./S
                Fchar = ----------------------------------------
                         \int_Fmin^Fmax df p_FD./S
              where Fmin = max(f) and \tilde(h) is the FFT of h.
       bw     Effective bandwidth (Hz) of h in the given noise
              background, defined as
                     \int_Fmin^Fmax df (f-Fchar).^2 p_FD./S
                bw = ---------------------------------------------------
                      \int_Fmin^Fmax df p_FD./S
              where Fmin = max(f) and \tilde(h) is the FFT of h.
       Tchar  Characteristic time at which the event occurs, defined as
                        \int dt t p_TD(t)
                Tchar = -----------------
                         \int dt p_TD(t)
       dur    Duration of the signal, defined as
                        \int dt (t-Tchar).^2 p_TD(t)
                Tchar = ----------------------------
                         \int dt p_TD(t)

     Notes:

     The power spectrum S is interpolated (if needed) to the frequencies of
     fft(h).

     Note that no windowing is used for computing FFTs windowing and
     averaging is not really sensible for a transient signal, since by
     definition the statistical properties of the signal are changing over
     time.  There may be problems for, e.g., band-passed noise bursts.

     The SNR and h_rss measures are restricted to the frequency interval
     [Fmin,Fmax], while h_peak is evaluated in the time domain (i.e., using
     data from the full frequency range.

     See FFTCONVENTIONS for information on the conventions used for Fourier
     transforms.
    """


    ###########################################################################
    #   Checks.
    ###########################################################################

    # ---- Optional arguments.
    t0 = int(t0)
    # ---- Is there an input noise spectrum?
    if (not S) or (not F0):
        # ---- Set flag to make dummy noise spectrum.
        noise = 0
    else:
        noise = 1
        # ---- Validate noise spectrum data.
        if type(S) is np.ndarray:
            raise ValueError('Spectrum S be a vector or empty array.')

        if not ((type(dF) == int and  dF > 0)  or not dF):
            raise ValueError('Spectrum frequency resolution dF '
                             'must be a positive scalar (if F0 is a scalar) '
                             'or empty (if F0 is a vector).')

        if not ( (type(F0) == int and F0 >= 0) or 
                 (type(F0) == np.ndarray and len(S) == len(F0))
                ):
            raise ValueError('Spectrum lowest frequency F0 must be a '
                             'non-negative scalar, '
                             'or a vector of equal length to S')

        if type(F0) == int and not dF:
            raise ValueError('Frequency sampling step not provided and F0 is a scalar')

        if not (type(F0) ==int) and dF:
            raise ValueError('Frequency sampling step provided, but F0 is a vector')

        # ---- Vector of sampled noise frequencies.
        if type(F0) == int:
            F = F0 + numpy.arange(0, len(S)) * dF
        else:
            if any(diff(F0)<=0):
                raise ValueError('F0 is not an ascending frequency series')
            F = F0

    # ---- Frequency range of analysis.
    if not Fmax:
        if noise:
            # ---- Default to highest frequency in supplied spectrum.
            Fmax = F[-1]
        else:
            # ---- Default to Nyquist.
            Fmax = numpy.floor(fs/2).astype(int)

    if not Fmin:
        if noise:
            # ---- Default to lowest frequency in supplied spectrum.
            Fmin = F[0]
        else:
            # ---- Default to DC.
            Fmin = 0

    # ---- Error checks.
    
    if (not (type(Fmax) in [int, numpy.int64, float, numpy.float64])
        or Fmax<=0):
        raise ValueError('Fmax must be positive and either a float or int')
    if fs<=0:
        raise ValueError('Sampling rate fs must be positive.')

    if not(type(t0) == int):
        raise ValueError('Timeseries start time t0 must be a scalar.')

    if (not (type(Fmin) in [int, numpy.int64, float, numpy.float64])
        or Fmin<0):
        raise ValueError('Frequency limit Fmin must be a non-negative scalar.')

    if Fmin >= Fmax:
        raise ValueError('Frequency limits must satisfy Fmin<Fmax.')

    # ---- Require positive sampling rate.
    if fs<=0:
        raise ValueError('Sampling rate fs must be positive.')

    if not(type(t0) == int):
        raise ValueError('Timeseries start time t0 must be a scalar.')


    ###########################################################################
    #   Preparations.
    ###########################################################################
    # ---- Number of data points.  Force to be even.
    N = h.shape[1]
    if not (N % 2 == 0):
        warnings.warn('Function not guaranteed for waveforms with odd number '
                      'of samples.  Dropping last sample')
        h = h[:, :-2]
        N = h.shape[1]    

    # ---- Duration of timeseries [sec].
    T = N/fs

    # ---- Verify that T is a power of 2.
    if T != 2**numpy.round(numpy.log2(T)):
        warnings.warn('Function is not guaranteed for timeseries durations '
                       'that are not powers of 2.')

    # ---- Sample times.
    t = t0 + numpy.arange(0, N) * 1/fs

    # ---- If no noise spectrum is supplied then make a dummy noise vector
    #      covering [0,Nyquist] Hz.  This will allow us to assume henceforth
    #      that S, F, F0, and dF are defined.
    if not noise:
        dF = 1./T
        F0 = 0
        S = 2*numpy.ones(int(N / 2))
        F = F0 + numpy.arange(0, S.size)*dF

    ###########################################################################
    #   Frequency-domain calculations.
    ###########################################################################

    # ---- Signal.

    # ---- Standard FFT works by column.
    hf = 1/fs*numpy.fft.fft(h)
    # ---- Vector holding frequencies in usual screwy FFT order:
    #        vector element:   [ 1  2  ...  N/2-1  N/2    N/2+1            N/2+2   ... N-1  N ]
    #        frequency (df):   [ 0  1  ...  N/2-2  N/2-1  (N/2 or -N/2)   -N/2+1  ... -2   -1 ]
    f_one_side = numpy.arange(0, N/2 + 1) / T
    # ---- Keep only frequency components in [Fmin,Fmax].  Note that most of
    #      the frequency-domain formulas include a factor of 2 to account for
    #      negative frequencies.
    index = numpy.where((f_one_side>=Fmin) & (f_one_side<=Fmax))[0]
    f = f_one_side[index]
    # ---- Distribution of signal energy in frequency.
    p_FD = numpy.sum(hf[:, index] * hf[:, index].conj(), axis=0)

    # ---- f=0 and f=Nyq bins should count for 1/2.
    if (f[0] == 0):
        p_FD[0] = 0.5 * p_FD[0]

    if (f[-1] == fs/2):
        p_FD[-1] = 0.5 * p_FD[-1]

    # ---- Noise.

    # ---- Does vector of sampled noise frequencies cover requested range?
    if ( (F[0] > f[0]) or (F[-1] < f[-1]) ):
        raise ValueError('Noise spectrum does not cover desired frequency range.')

    # ---- Force interpolation of S from sampled frequencies F to data
    #      frequencies f.
    interp_function = interp1d(F, S)
    S = interp_function(f)
    F = f
    dF = 1./T
    F0 = F[0]


    # ---- All set to do calculations.  Assured that f, F interchangable.

    # ---- SNR^2 versus frequency.
    SNR2f = 2*p_FD / S

    # ---- SNR on [Fmin,Fmax].
    SNR = (dF*sum(SNR2f))**0.5

    # ---- Characteristic frequency.
    Fchar = sum(f * SNR2f) / sum(SNR2f)

    # ---- Characteristic bandwidth.
    bw = (sum((f-Fchar)**2 * SNR2f) / sum(SNR2f))**0.5

    # ---- RSS amplitude.
    h_rss = (dF * sum(p_FD))**0.5


    ###########################################################################
    #   Time-domain calculations.
    ###########################################################################

    # ---- Distribution of signal energy in time.
    p_TD = (h**2).sum(0)

    # ---- Peak amplitude.
    h_peak = max(p_TD)**0.5

    # ---- Characteristic time.
    Tchar = sum(t * p_TD) / sum(p_TD)

    # ---- Duration.
    dur = (sum((t-Tchar)**2 * p_TD) / sum(p_TD))**0.5

    # ---- Done

    return SNR, h_rss, h_peak, Fchar, bw, Tchar, dur
