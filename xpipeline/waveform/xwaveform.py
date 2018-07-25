from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from gwpy.timeseries import TimeSeries
from .xoptimalsnr import xoptimalsnr

import numpy
import scipy.signal
import math
import os
import warnings

###########################################################################
#   Physical constants.
###########################################################################

#-- 1 parsec (m)
pc = 3.08567756707701e+016

#-- 1 Mega-parsec (m)
mpc = 1e6*pc

#-- speed of light (m/s)
c = 299792458

newtong = 6.67e-11
speedoflight_meterspersecond = 299792458
megaparsec_meters = 1e6*3.26*speedoflight_meterspersecond*365.25*86400
solarmass_kilograms = 1.99e30
solarmass_seconds = newtong*solarmass_kilograms/speedoflight_meterspersecond**3
solarmass_meters = newtong*solarmass_kilograms/speedoflight_meterspersecond**2
sqrt_of_neg_1 = numpy.sqrt(numpy.array([-1], dtype=complex))

def xmakewaveform(family, parameters, T, T0, fs, **kwargs):
    """Create timeseries data containing a simulated GW.

       Parameters:

           family (str):

               A char or cell array containing one of (not case sensitive):
                  adi-a,...,-e  accretion disk
                                instability waveforms as used in STAMP searches
                  AORDDMJ10  Abdikamalov et al. accretion-induced collapse waveforms
                             [PRD81 044012 2010],
                  BLDOM07  Burrows et al. acoustic mechanism supernova waveforms
                             [ApJ655 416 2007],
                  CG  cosine-Gaussians (special case of 'chirplet'),
                  chirplet  chirping sine/cosine-Gaussian,
                  cusp  waveform from a cosmic string cusp (see cuspgw.m),
                  DFM  NOT FUNCTIONAL. Dimmelmeier-Font-Mueller 2002 supernova waveforms,
                             relativistic simulations only [A&A393 523 2002],
                             WARNING: Normalization may be incorrect.
                  DFMc  NOT FUNCTIONAL. Conditioned Dimmelmeier-Font-Mueller supernova
                             waveforms.  Same as DFM except detrended.
                             WARNING: Normalization may be incorrect.
                  DFM02  Dimmelmeier-Font-Mueller 2002 supernova waveforms,
                             relativistic & Newtonian simulations, conditioned
                             and detrended [A&A393 523 2002].
                             WARNING: Normalization may be incorrect.
                  DOJMM07  Dimmelmeier et al. 2007 supernova waveforms [PRL98
                             251101 2007],
                  DOMJ08  Dimmelmeier et al. 2008 supernova waveforms [PRD78
                             064056 2008]
                  DS  damped sinusoids,
                  DS2P  2-polarization damped sinusoids with smooth turn-on,
                  ebbh-a,...,-e  eccentric binary black-hole waveforms as used in STAMP
                                 searches
                  G  Gaussians (special case of 'chirplet'),
                  inspiral    post-Newtonian inspiral,
                  inspiralgated  As inspiralsmooth, with a gate applied centered 1 sec
                             before coalescence time. Gate is 0.5 sec hann roll-off,
                             0.3 sec zeros, 0.5 sec hann roll-on.
                  inspiralsmooth  post-Newtonian inspiral, tapered at the end
                  iscochirp*  Broadband chirps from ISCO waves around rotating black holes,
                             from https://arxiv.org/abs/1602.03634 (Van Putten)
                  kotake09  Kotake et al. supernova waveform -- NEED INFO,
                  lalinspiral  post-Newtonian inspiral, generated using LALsuite
                  Lazarus  Lazarus project black-hole merger waveforms (analytic
                             approximation) [PRD65 124012 2002],
                  longbar  long-lived bar mode,
                  magnetar  GW emission by magnetar undergoing spin down caused by EM
                             and/or EM emission,
                  marek09  Marek et al. 2009 supernova waveforms -- NEED INFO,
                  methodspaper  One of four numerical supernova waveforms from
                             http://arxiv.org/abs/1110.5107v2.
                  mono       Monochromatic waveforms used by STAMP-AS (including mono-*,
                             line-*, quad-*) see https://wiki.ligo.org/Main/O1analysisInjections,
                  mob09      Murphy-Ott-Burrows 2009 supernova waveforms
                             [arXiv:0907.4762],
                  murphy09  Murphy+09 waveforms
                  mvp       van Putten based waveform model -- NEED INFO,
                  NCSA-CAM  eccentric binary neutron star mergers, generated with the
                             NCSA-CAM waveform model. The waveforms start at 15Hz and
                             continue through merger and ringdown. See
                             https://dcc.ligo.org/DocDB/0140/G1700265/001/ecc_presentation.pdf
                  OBDL06     Ott, Burrows, Dessart, & Livne supernova waveforms
                             [PRL 96 201102 2006]
                  OBLW04     Ott-Burrows-Livne-Walder 2004 supernova waveforms
                             [ApJ600 834 2004],
                  onecyclesine  single cycle (phase -pi->0->pi) of a sine wave,
                  osnsearch  test waveform used in S5/S6 optical supernova search,
                  o1snews  Waveforms for O1 X SNEWS Search.
                  O09  Ott SASI supernova waveforms [CQG26 063001 2009],
                  O10  Ott 2010 supernova waveforms -- NEED INFO,
                  PCA  Principle Component Analysis Wave Forms
                  piro  NOT FUNCTIONAL. Piro-based waveform model -- NEED INFO,
                  scalarchirplet  scalar-mode chirping cosine-Gaussian,
                  pmns_pca  A large number of precalculated Post-Merger Remnant Signals for a
                             variety of EOS mass configurations. some include
                             'nl3_1515' 'tma_1215' 'tma_1212' 'dd2_1314'
                             'ls375_1215' 'dd2_1215' 'dd2_165165' 'tma_1515'
                             'ls375_1212' 'ls220_1414'
                  scalarsn  Model A follows solution number 4 from Table 1 for a 1.5
                             solar mass star undergoing spontaneous scalarization. See
                             http://journals.aps.org/prd/abstract/10.1103/PhysRevD.58.064019
                             Model B follows the collpase of aspherical NS star to
                             blackhole following solution B in table 1.
                             http://journals.aps.org/prd/abstract/10.1103/PhysRevD.57.4789
                             Models C, D, E are Novak and Ibanez scalar supernovae waveforms.
                            http://stacks.iop.org/0004-637X/533/i=1/a=392
                  scheide    Scheidegger et al. 2010 equatorial waveforms -- NEED INFO,
                  scheidp     Scheidegger et al. 2010 for polar waveforms -- NEED INFO,
                         SG  sine-Gaussians (special case of 'chirplet'),
                  snmp  Waveforms used for SN Methods Paper that are not included in
                         osnsearch.
                  SNN94  Shibata, Nakao, and Nakamura scalar-mode GWs from
                             gravitational collapse [PRD 50 7304 1994],
                  stamp_*  STAMP-AS waveforms for O1 long-duration search see
                             https://wiki.ligo.org/Bursts/Allsky20160425
                  WNB  windowed, band-limited white-noise burst (2 independent
                             polarizations).
                  Yakunin2010  Yakunin et al. http://arxiv.org/abs/1005.0779
                  zero  null signal (hp=0=hc),
                  ZM  Zwerger-Mueller supernova waveforms [A&A320 209 1997].
                             WARNING: Normalization may be incorrect see eqn (21) of
                             the reference.
           params (Array, or tilde-delimited string):

               containing parameters for the specified waveform type.  The format and
               definition of the parameters depend on the waveform type.  For
               details on a specific waveform, e.g. WNB, do xmakewaveform('WNB').

           T (`int`):
               Duration of the waveform.  Must have T*fs = integer.
               Recommend T>=1 for astrophysical waveforms.

           T0 (`int`):
               Desired peak time of the waveform, as measured in hrss
               for (hp.^+hc.^2).^(0.5).  Must have 0<T0<T.  See "Notes" below
               for caveats.

           fs (`int`):
               Sampling rate of the time series.  T*fs must be an
           integer.

       Returns:
           t:
               Times at which the waveform is sampled, starting from zero.
           hp:
               Plus polarization waveform [strain].
           hc:
               Cross polarization waveform [strain].
           hb:
               Breathing (scalar) polarization waveform [strain].  Identically
               zero for all waveform types unless explicitly stated otherwise.

       Notes:

            The waveform is interpolated so that the RSS-weighted center time,
            [\int dt t (hp.^+hc.^2+hb.^2)].^(0.5), is at T0.  Note that for
            2-polarization waveforms that include an inclination angle (eg, DS2P)
            different inclination angles will cause hp, hc to be weighted differently
            and therefore will produce different shifts.

            The cusp waveform (generated by P. Shawhan's cuspgw function) and most
            astrophysical waveforms that are loaded from files are always 1 sec long
            with sampling rate 16384 Hz when created/loaded. These are resampled
            linearly and truncated or zero-padded to the desired length.  One should
            be cautious in using T<1 for these waveforms.
    """

    ###########################################################################
    #   Check inputs.
    ###########################################################################

    inclinationtype = kwargs.pop('inclinationtype', 'total')
    if inclinationtype == 'total':
        use_total_inclination = True
    elif inclinationtype == 'orbital':
        use_total_inclination = False

    hrss = kwargs.pop('hrss', -1)
    filedir = kwargs.pop('catalogdirectory',
                         '/home/scoughlin/PhD/waveforms/')

    ###########################################################################
    #   Other preparatory.
    ###########################################################################

    #----- Boolean: If true then waveform "type" specifies a waveform made
    #      outside of this function (loaded from file, or using another
    #      function) and which therefore may need to be re-sampled, shifted in
    #      time, zero-padded, or truncated.
    pregen = 1      #-- By default perform conditioning on everything.
    pregen_fs = fs  #-- sampling rate of pregenerated waveform.
    pregen_T = T    #-- duration of pregenerated waveform.

    #----- Time in column vector.
    t = numpy.arange(0, T, 1./fs)

    #----- If parameters are supplied as tilde-delimited string, then convert
    #      to cell or double array.
    if type(parameters) in [str, numpy.str_, numpy.str, numpy.string_]:
        parameters_tmp = parameters.split('~')
        parameters = []
        for param in parameters_tmp:
            try:
                parameters.append(float(param))
            except:
                parameters.append(param)
        

    if family.lower() in ['aorddmj10', 'dojmm07', 'domj08', 'mob09', 'obdl06', 'oblw04']:

        # ----  Entry from one of various supernova waveform catalogs
        #       (pregenerated).

        # ---- Find specified waveform type.
        name = parameters[1]

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))
        pregen = 1
        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        hp = 10./parameters[0] * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value)

    elif family.lower() in ['murphy09']:

        # ----  Entry from one of various supernova waveform catalogs
        #       (pregenerated).

        # ---- Find specified waveform type.
        name = parameters[1]
        hrss = parameters[0]

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))
        pregen = 1
        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Find specified waveform type.
        hp = hp / hrss
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value)

    elif family.lower() in ['adi-a', 'adi-b', 'adi-c',
                            'adi-d','adi-e','ebbh-a','ebbh-d','ebbh-e']:

        familyname = family.split('-')[0]
        hp = TimeSeries.read(os.path.join(filedir,
                                          familyname + '.hdf5'),
                             path='/{0}/{1}/hp'.format(familyname, family))
        hc = TimeSeries.read(os.path.join(filedir,
                                          familyname + '.hdf5'),
                             path='/{0}/{1}/hc'.format(familyname, family))

        # ---- Parse parameters.
        distance = parameters[0]  #-- mpc
        ciota    = parameters[1]  #-- inclination

        # ---- Extract +/x polarizations. Default distance is 1 mpc. Assume
        #      l=m=2 mode dominates when applying inclination.
        hp = 0.5 * (1 + ciota**2) * hp / distance
        hc = ciota * hc / distance

        # ----  Set "pregenerated" variables.
        pregen = 1
        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

    elif family.lower() in ['bldom07']:

        # ----  acoustic.
        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        distance = 0.1 * parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        hp = (1./distance) * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value)

    elif family.lower() in ['cg']:

        # ---- cosine-Gaussians.
        h_rss = parameters[0]
        Q = parameters[1]
        f0 = parameters[2]
        h_peak = h_rss * ((4 * math.pi**0.5 * f0)/(Q * (1 + math.exp(-Q**2))))**0.5
        tau = Q / (2**0.5 * math.pi * f0)
        hp = TimeSeries(h_peak * numpy.cos(2 * math.pi * (t-T0) * f0) * numpy.exp(-(t-T0)**2/tau**2),
                        name=family, dx = 1./fs)
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value)

        # ---- Turn off default interpolation (symmetric ad hoc waveform).
        pregen = 0

    elif family.lower() in ['chirplet']:

        # ---- Chirplet - Gaussian-modulated sinusoid with frequency
        #      changing linearly with time.  Put chirping cosine-Gaussian
        #      in plus polarization, chirping sine-Gaussian in cross.

        # ---- Required parameters.
        h_rss = parameters[0]
        tau = parameters[1]
        f0 = parameters[2]

        # ---- Optional parameters.
        alpha = 0
        delta = 0
        ciota = 1
        if len(parameters) >= 4:
            alpha = parameters[3]

        if len(parameters) >= 5:
            delta = parameters[4]

        if len(parameters) >= 6:
            ciota = parameters[5]

        # ---- Waveform.
        h = h_rss * (
                     numpy.exp(
                              (-1 + sqrt_of_neg_1 * alpha) * (t-T0)**2 /
                              (4 * tau**2)
                              + sqrt_of_neg_1 * 2 * numpy.pi * (t-T0) *f0
                              + sqrt_of_neg_1 * delta)
                     / 
                     (2 * math.pi * tau**2)**(1/4))
        h = TimeSeries(h, dx=1./fs, name='chirplet')
        hp = 1./2 * (1 + ciota**2) * h.real
        hc = ciota * h.imag

        # ---- Turn off default interpolation (ad hoc waveform is designed
        #      to produce desired T0).
        pregen = 0


    elif family.lower() in ['dfm02']:

        # ----  Dimmelmeier, Font, and Mueller supernova waveform (pregenerated).
        distance = 0.1 * parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        hp = (1./distance) * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)

    elif family.lower() in ['ds']:

        # ---- Damped sinusoids with sharp turn-on.
        h_peak = parameters[0]
        tau = parameters[1]
        f0 = parameters[2]

        hp = TimeSeries(numpy.zeros(t.size), dx=1/fs,
                        name=family)
        hc = TimeSeries(numpy.zeros(t.size), dx=1/fs,
                        name=family)
        k = numpy.where(t>=T0)[0]
        hp[k] = h_peak * (numpy.cos(2 * numpy.pi * 
                                   (t[k] - T0) * f0)
                         * numpy.exp(-(t[k] - T0) / tau))
        hc[k] = h_peak * (numpy.sin(2 * numpy.pi * 
                                    (t[k] - T0) * f0)
                         * numpy.exp(-(t[k] - T0) / tau))

        # ---- Force interpolation to get correct T0, since waveform is
        #      asymmetric.
        pregen = 1

    elif family.lower() in ['ds2p']:

        # ---- Damped sinusoids, 2-polarization, with smooth turn-on.
        h_peak = parameters[0]
        tau = parameters[1]
        f0 = parameters[2]
        delta = parameters[3]
        ciota = parameters[4]

        # ---- Smooth turn on: prepend time-reversed damped sinusoid with
        #      same parameters, except tau->tau/10.
        tau_vec = tau * numpy.ones(t.size)
        k = numpy.where(t<T0)[0]
        tau_vec[k] = tau_vec[k] / 10.
        hp = h_peak * 0.5 * (1+ciota**2) * (
                                            numpy.cos(2*pi * (t-T0)*f0 + delta)
                                            * numpy.exp(-abs(t-T0) / tau_vec)
                                           )
        hc = h_peak * ciota * (
                               numpy.sin(2 * pi * (t-T0)*f0 + delta)
                               *exp(-abs(t-T0) / tau_vec)
                              )
        hp = TimeSeries(hp, name=family, dx=1./fs)
        hc = TimeSeries(hc, name=family, dx=1./fs)

        # ---- Force interpolation to get correct T0, since waveform is
        #      asymmetric.
        pregen = 1

    elif family.lower() in ['g']:

        # ---- Gaussians.
        h_rss = parameters[0]
        tau = parameters[1]

        h_peak = h_rss * (2/(math.pi * tau**2))**(1/4)
        hp = TimeSeries(h_peak * numpy.exp(-(t-T0)**2 / tau*2), dx=1./fs,
                        name=family)
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)

        # ---- Turn off default interpolation (symmetric ad hoc waveform).
        pregen = 0

    elif family.lower() in ['inspiral']:
        """
        import lal
        import lalsimulation

        # ---- Post-Newtonian inspiral waveform, using INSPIRAL2PN.
        mass1 = parameters[0]
        mass2 = parameters[1]
        iota  = numpy.arccos(parameters[2])
        dist  = parameters[3]

        lalsimulation.SimInspiralChooseTDWaveform(m1=mass1, m2=mass2,
                                                  s1x=0, s1y=0, s1z=0, s2x=0,
                                                  s2y=0, s2z=0, distance=dist,
                                                  inclination=iota, phiRef=0,
                      longAscNodes=0, eccentricity=0, meanPerAno=0,
                      deltaT=1./fs, f_min =20, f_ref=20,
                      params=lal.CreateDict(),
                      approximant=lalsimulation.TaylorF2) 
        [hp, hc] = inspiral2pn(mass1,mass2,iota,dist,fs,T,T0)

        # ---- Turn off default interpolation (leave coalescence time
        #      intact).
        pregen = 0
        """
    elif family.lower() in ['inspiralsmooth', 'inspiralgated']:

        """
        # ---- Post-Newtonian inspiral waveform, using INSPIRAL2PN.
        mass1 = parameters(1)
        mass2 = parameters(2)
        iota  = acos(parameters(3))
        dist  = parameters(4)
        [hp, hc, tpn, tinspiral] = inspiral2pn(mass1,mass2,iota,dist,fs,T,T0)
        inspiralEnd = length(tinspiral)
        taperLength = round(fs/512)
        taper = hann(2*taperLength+1)
        if(length(hp)<10*taperLength)
          error('Inspiral waveform is too short for smoothing')
        end
        hp(inspiralEnd-taperLength:inspiralEnd) = ...
            hp(inspiralEnd-taperLength:inspiralEnd).*taper(end-taperLength:end)
        hc(inspiralEnd-taperLength:inspiralEnd) = ...
            hc(inspiralEnd-taperLength:inspiralEnd).*taper(end-taperLength:end)

        if strcmpi(lower(type),'inspiralgated')
            # ---- Apply gate to waveform centered 1 sec before coalescence
            #      time. Gate is 0.5 sec hann roll-off, 0.3 sec zeros,
            #      0.5 sec hann roll-on. So gate starts 1+0.15*0.5 sec = 1.65 sec
            #      before coalescence.
            taper = hann(fs)
            gate = [taper(end/2+1:end) zeros(round(0.3*fs),1) taper(1:end/2)]
            gatelength = length(gate)
            [val,idx] = min(abs(t-(T0-1.65)))
            hp((idx+1):(idx+gatelength)) =  gate .* hp((idx+1):(idx+gatelength))
            hc((idx+1):(idx+gatelength)) =  gate .* hc((idx+1):(idx+gatelength))
        end

        # ---- Turn off default interpolation (leave coalescence time
        #      intact).
        pregen = 0
        """
    elif family.lower() in ['kotake09']:

        # ----  Dimmelmeier, Font, and Mueller supernova waveform (pregenerated).
        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        hc = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value


        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        hp = 10. / distance * hp
        hc = 10. / distance * hc

    elif family.lower() in ['lalinspiral']:

        """
        # --- add sampling frequency to parameters
        if iscell(parameters)
          parameters = {fs, parameters{:}}
        else
          parameters = num2cell([fs parameters(:)])
        end

        # --- realign binary so that inclination referes to angle between
        # total angular momentum and line of sight
        if use_total_inclination
          parameters = alignTotalSpin(parameters)
        end

        # ---- Post-Newtonian inspiral waveform, using LALGenerateInspiral
        [tLAL, hpLAL, hcLAL, freq, a1LAL, a2LAL] = lalinspiral(parameters)

        # ---- Align coalescence time with T0 in the requested [0,T]
        #      window, chop or fill with zeros to complete the waveform
        hp=zeros(size(t))
        hc=zeros(size(t))
        nTimeBins = min(length(tLAL),sum(t<=T0))-1
        copyRange = (length(tLAL)-nTimeBins):length(tLAL)
        pasteRange = find(t<=T0)
        pasteRange = pasteRange((end-nTimeBins):end)
        hp(pasteRange) = hpLAL(copyRange)
        hc(pasteRange) = hcLAL(copyRange)
        a1(pasteRange) = a1LAL(copyRange)
        a2(pasteRange) = a2LAL(copyRange)

        # ---- Turn off default interpolation (leave coalescence time
        #      intact).
        pregen = 0
        """

    elif family.lower() in ['lazarus']:
        """

        # ---- Lazarus Project black-hole mergers.  Use analytic
        #      approximation from Baker, Campanelli, Lousto, and Takahashi,
        #      PRD 65 124012 2002.

        # ---- Input parameters.
        mass_solarMasses = parameters[0]
        distance_megaParsecs = parameters[1]
        ciota = parameters[2]

        # ---- Constants of analytic approximation.
        alpha0 = 0.0085
        omega0 = 0.2894
        tomega0 = 33.00
        tomega1 = 67.05
        sigma0 = 0.2192
        tsigma0 = 32.54
        tsigma1 = 62.92
        omegaQNM = 0.55
        sigmaQNM = -0.073
        phi0 = -4.40
        a0 = -6.31

        # ---- This range of t includes the entire range where the envelope
        #      A > 1e-4*max(A).
        dt = 0.1  #-- fine enough resolution (in units of system mass)
        tmin = 0
        tmax = 200 #-- range over which analytic approximation has support
        t = numpy.arange(tmin, tmax, dt)

        # ---- Construct phase phi(t) piecewise in time.
        phi = numpy.zeros(t.size)
        k = numpy.where(t<tomega0)[0]
        phi[k] = phi0 + omega0*(t[k] - tomega0) + alpha0/2*(t[k]-tomega0)**2
        k = numpy.where( (t>=tomega0) & (t<tomega1) )[0]
        phi[k] = phi0 + omega0*(t[k] - tomega0) ...
            + (omega0 - omegaQNM)/(tomega0 - tomega1) / 2*(t[k] - tomega0)**2
        k = numpy.where(t>=tomega1)[0]
        phi[k] = phi0 + 1./2 * ((omega0 + omegaQNM) * (tomega1-tomega0)
            + omegaQNM*(t[k] - tomega1))

        # ---- Construct amplitude A(t) piecewise in time.
        lnA = numpy.zeros(t.size)
        k = numpy.where(t<tsigma0)[0]
        lnA[k] = a0 + sigma0*(t[k]-tsigma0)
        k = numpy.where( (t>=tsigma0) & (t<tsigma1) )[0]
        lnA[k] = a0 + sigma0*(t[k]-tsigma0) ...
            + 1./2*(sigma0-sigmaQNM)/(tsigma0-tsigma1) * (t[k]-tsigma0)**2
        k = numpy.where(t>=tsigma1)[0]
        lnA[k] = a0 + 1./2*(sigma0+sigmaQNM)*(tsigma1-tsigma0) ...
            + sigmaQNM*(t[k]-tsigma1)
        A = numpy.exp(lnA)

        # ---- Construct r*Psi_4 variable.
        rPsi4 = A * numpy.exp(-srt_of_neg_1 * phi)

        # ---- Integrate twice in time to get hplus, hcross.
        ddhp = real(2*rPsi4)
        ddhc = imag(2*rPsi4)
        dhp = dt*real(cumsum(2*rPsi4))
        dhc = dt*imag(cumsum(2*rPsi4))
        hp = dt^2*real(cumsum(cumsum(2*rPsi4)))
        hc = dt^2*imag(cumsum(cumsum(2*rPsi4)))
        # ---- Remove linear trend in hplus.
        tint = t(end)-hp(end)/dhp(end)
        m = dhp(end)
        b = hp(end) - dhp(end)*t(end)
        hp_trend = zeros(size(t))
        k = find(t>tint)
        hp_trend(k) = m*t(k)+b
        hp = hp - hp_trend
        # ---- Remove linear trend in hcross.
        tint = t(end)-hc(end)/dhc(end)
        m = dhc(end)
        b = hc(end) - dhc(end)*t(end)
        hc_trend = zeros(size(t))
        k = find(t>tint)
        hc_trend(k) = m*t(k)+b
        hc = hc - hc_trend

        # ---- Rescale to seconds and megaParsecs.
        t = t * mass_solarMasses * solarmass_seconds
        hp = hp * ( mass_solarMasses * solarmass_meters ) ...
            / ( distance_megaParsecs * megaparsec_meters )
        hc = hc * ( mass_solarMasses * solarmass_meters ) ...
            / ( distance_megaParsecs * megaparsec_meters )

        # ---- Rescale waveforms according to inclination angle.
        hp = hp*(0.5*(1+ciota^2))
        hc = hc*ciota

        # ---- Set pregen flag so we will interpolate to desired sampling
        #      rate, etc.
        pregen = 1
        pregen_T = length(t) * dt * mass_solarMasses * solarmass_seconds
        pregen_fs = length(t)/pregen_T

        """
    elif family.lower() in ['longbar']:

        """
        # ---- Figure out parameters and call helper function.
        # [D,M,L,R,f,tau,phi,theta]
        D  = parameters(1)
        M = parameters(2)
        L = parameters(3)
        R  = parameters(4)
        f  = parameters(5)
        tau  = parameters(6)
        phi  = parameters(7)
        theta  = parameters(8)
        [dummy, hp, hc] = longbar(D, M, L, R, f, tau, 0, 0)
        hb = zeros(size(hp))

        # ---- Long-lived bar mode waveform (pregenerated).
        pregen = 1
        pregen_fs = fs
        pregen_T = length(hp)/pregen_fs

        # ---- Time vector for pregenerated waveform.
        pregen_t = [0:1/pregen_fs:pregen_T-1/pregen_fs].T
        # ---- Measure peak/characteristic time of pregenerated waveform.
        [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = xoptimalsnr( ...
            [hp,hc,hb],pregen_t(1),pregen_fs,[],[],[], ...
            1/pregen_T,0.5*pregen_fs-1/pregen_T ...
        )
        # ---- Shift pregenerated time vector by required amount.
        pregen_t = pregen_t + (T0 - Tchar)
        # ---- Interpolate to correct sampling rate and peak time.
        # ---- Time in column vector.
        t = [0:1/fs:T-1/fs].T
        # ---- Find desired times which overlap pregen_t.
        k = find(t>=pregen_t(1) & t<=pregen_t(end))
        # ---- Interpolate, with zero padding if pregenerated waveform is too
        #      short.  (Truncation of long waveforms handled automatically by
        #      specifying vector "t".)
        hp_interp = zeros(size(t))
        hc_interp = zeros(size(t))
        hb_interp = zeros(size(t))
        hp_interp(k) = interp1(pregen_t,hp,t(k))
        hc_interp(k) = interp1(pregen_t,hc,t(k))
        hb_interp(k) = interp1(pregen_t,hb,t(k))
        hp = hp_interp
        hc = hc_interp
        hb = hb_interp
        hp = 0.5 * (1 + (cos(theta))^2) * hp
        hc = cos(theta) * hc
        pregen = 0
        """

    elif family.lower() in ['magnetar']:

        # ----  Lasky et al. magnetar waveform see https://wiki.ligo.org/Main/LongDurationRemnantWaveforms .

        # ---- Parse parameters.
        hrss  = parameters[0]  #-- Hz^-0.5
        ciota = parameters[1]  #-- inclination

        # ---- Load waveform structure, if not already loaded or supplied.
        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, family))


        hc = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hc'.format(family, family))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ----  Set "pregenerated" variables.
        pregen = 1

        # ---- Extract +/x polarizations and rescale to desired hrss at
        #      OPTIMAL inclination.
        hrss_default = numpy.sum((hp**2+hc**2)/hp.sample_rate.value)**0.5
        hp = hp * hrss / hrss_default
        hc = hc * hrss / hrss_default

        # ---- Applying inclination factors, assuming l=|m|=2 mode
        #      dominates the emission.
        hp = 0.5 * (1+ciota**2) * hp
        hc = ciota * hc

    elif family.lower() in ['marek09']:

        # ----  Marek 2009 supernova matter waveform (pregenerated).
        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        hp = 10. / distance * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)


    elif family.lower() in ['methodspaper']:

        # ---- A set of four numerical waveforms to be included in the
        #      supernova methods paper. See also Table 3 of the review
        #      http://arxiv.org/abs/1110.5107v2. The names are
        #      gw_s15s7b2.dat, sasi3DL68randgwraynew_thphcPole,
        #      polar_R4E1FC_L and TK11_B12X1B01.
        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        hc = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value


        # ---- Waveforms have a default distance of 10 kpc.  Rescale if
        #      desired.  For historical reasons, distance = 0 is treated as
        #      the default of 10 kpc.
        if distance != 0:
            hp = 10. / distance * hp
            hc = 10. / distance * hc


    elif family.lower() in ['morozova2018']:

        # ----  Burrows 2018 2D SN waveform (pregenerated).
        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Read waveform, which is defined at a range of 1cm.
        #      therefore we convert it first to 1 kpc then divide
        #      through by provided distance
        #      Make sure h is same type of vector (column) as t.
        hp = 1. / distance * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)


    elif family.lower() in ['mono', 'mono-a', 'mono-b', 'mono-c',
                            'line-a', 'line-b', 'quad-a', 'quad-b']:

        """
        # ---- Long-duration sine waves used by STAMP-AS. See
        #      https://wiki.ligo.org/Main/O1analysisInjections (date: 2016.03.05).
        #        parameters(1) - hrss
        #        parameters(2) - tau (duration) - OPTIONAL
        #        parameters(3) - f0 (start frequency) - OPTIONAL
        #        parameters(4) - fp0 (frequency derivative) - OPTIONAL
        #        parameters(5) - fpp0 (frequency second derivative) - OPTIONAL
        h_rss = parameters[0]
        # ---- Default parameters: zero frequency derivatives.
        fp0  = 0
        fpp0 = 0
        if family.lower() == 'mono':
            tau = parameters[1]
            f0 = parameters[2]
            if len(parameters)>=4:
                fp0 = parameters[3]

            if len(parameters)>=5:
                fpp0 = parameters[4]

        # ---- Special cases with hardwired parameters.
        elif family.lower() == 'mono-a':
            tau = 150
            f0  = 90
        elif family.lower() == 'mono-b':
            tau = 250
            f0  = 550
        elif family.lower() == 'mono-c':
            tau = 250
            f0  = 405
        elif family.lower() == 'line-a':
            tau = 250
            f0  = 50
            fp0 = 0.6
        elif family.lower() == 'line-b':
            tau = 100
            f0  = 900
            fp0 = -2
        elif family.lower() == 'quad-a':
            tau =  30
            f0  =  50
            fp0 =  0
            fpp0 = 1./3
        elif family.lower() == 'quad-b':
            tau =  70
            f0  =  500
            fp0 =  0
            fpp0 = 0.041

        # ---- Sanity checks: these are meant to be long-duration injections.
        if tau < 2:
            warnings.warn('Injection shorter than 2 seconds. '
                          'Tapering may give unexpected behaviour.')

        # ---- Linearly polarised sine wave with frequency f(t) starting at time
        #      t=0. The frequency is defined as
        #        f(t)   = f0 + (df/dt)*t + (1/2)*(d^2f/dt^2)*t^2
        #      t=0. Therefore the phase will be
        #        Phi(t) = 2 * pi * \int_0^t dt f(t)
        #               = 2 * pi * \int_0^t dt [f0 + (df/dt)*t + (1/2)*(d^2f/dt^2)*t^2]
        #               = 2 * pi * [f0 t + (df/dt)*t^2/2 + (1/2)*(d^2f/dt^2)*t^3/3]
        t_temp = numpy.arange(0, tau, 1./fs)
        hp_temp = numpy.sin(2 * numpy.pi * (f0 * t_temp +
                                            1./2 * fp0 * t_temp**2 +
                                            1./6 * fpp0 * t_temp**3))

        # ---- Use Hann window to taper the ends of the waveform.
        window = scipy.signal.hanning(2*fs)
        hp_temp(0:fs) = hp_temp([1:1:fs]) .* window(1:fs)
        hp_temp([end-fs+1:end]) = hp_temp([end-fs+1:1:end]) .* window(fs+1:end)

        # ---- Rescale to the desired h_rss.
        hrss_norm = (numpy.sum(hp_temp**2)/fs)**0.5
        hp = hp_temp * h_rss/hrss_norm

        # ---- Define cross and scalar-mode components to be zero.
        hc = zeros(size(hp))
        hb = zeros(size(hp))

        # ---- Turn on interpolation to handle time-shifting of waveform to
        #      desired peak time.
        pregen = 1
        pregen_fs = fs
        pregen_T = tau
        """

    elif family.lower() in ['mvp']:

        """
        # ----  van Putten based waveform model (pregenerated).
        pregen = 1
        pregen_fs = 16384
        pregen_T = 30

        # ---- Load waveform structure, if not already loaded or supplied.
        if (isempty(catalog))
            # disp('Loading MvP catalog.')
            load([filedir 'MvPcat.mat'])
        end

        # ---- Find specified waveform type.
        NAME = parameters{1,2}
        for k=1:length(MvPcat)
            if (strcmp(MvPcat(k).name,NAME))
                wf_number = k
            end
        end

        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        distance = 0.1*parameters{1}
        # note: if these are read in from column ascii files, no need to invert vector dims
        hp = 1./distance*MvPcat(wf_number).hp
        hc = 1./distance*MvPcat(wf_number).hc
        """

    elif family.lower() in ['ncsa-cam']:

        # ---- These waveforms are generated with the NCSA-CAM waveform
        #      model. The main features of the model are summarized in the
        #      DCC document https://dcc.ligo.org/DocDB/0140/G1700265/001/ecc_presentation.pdf
        #      The files present waveforms describing BNS systems with mass
        #      combinations [1.4, 1.4] and [3, 3] solar masses,
        #      eccentricities of 0.2, 0.4, or 0.6, and located at a
        #      luminosity distance of 100 mpc. The waveforms are generated
        #      from an initial gravitational wave frequency of 15Hz.

        # ---- Each cell array entry holds an Nx2 array with columns
        #      [h_plus(t), h_cross(t)] for one of the six waveforms. Each
        #      waveform is sampled at 16384 samples per second. The
        #      waveform parameters are:
        #       1  m1=m2=1.4Msun, ecc=0.2
        #       2  m1=m2=1.4Msun, ecc=0.4
        #       3  m1=m2=1.4Msun, ecc=0.6
        #       4  m1=m2=3Msun, ecc=0.2
        #       5  m1=m2=3Msun, ecc=0.4
        #       6  m1=m2=3Msun, ecc=0.6

        # ---- Parameters.
        mass  = parameters[0]
        ecc   = parameters[1]
        ciota = parameters[2]
        dist  = parameters[3]

        # ---- Select desired waveform.
        if [mass, ecc] == [1.4, 0.2]:
            name = 'mass1.4_ecc0.2'
        elif [mass, ecc] == [1.4, 0.4]:
            name = 'mass1.4_ecc0.4'
        elif [mass, ecc] == [1.4, 0.6]:
            name = 'mass1.4_ecc0.6'
        elif [mass, ecc] == [3, 0.2]:
            name = 'mass3.0_ecc0.2'
        elif [mass, ecc] == [3, 0.4]:
            name = 'mass3.0_ecc0.4'
        elif [mass, ecc] == [3, 0.6]:
            name = 'mass3.0_ecc0.6'
        else:
            raise ValueError('mass, ecc values not an allowed pair')

        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        hc = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Extract waveforms and rescale as needed.
        hp = 100./dist * hp * (1+ciota**2) / 2
        hc = 100./dist * hc * ciota


    elif family.lower() in ['onecyclesine']:

        # ---- One cycle of a sine wave.
        h_peak = parameters[0]
        f0 = parameters[1]
        hp = TimeSeries(h_peak * numpy.sin(2 * pi *(t-T0) * f0),
                        dx = 1./fs,
                        name='onecyclesine'
                       )
        k = numpy.where(numpy.abs(t-T0) >= 1./numpy.abs(2*f0))
        hp[k] = 0
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)

        # ---- Turn off default interpolation (symmetric ad hoc waveform).
        pregen = 0

    elif family.lower() in ['osnsearch']:

        # ----  Entry from one of various waveform catalogs (pregenerated).
        #       The waveforms in this catalog are of different lengths, but
        #       have the same sampling rate.
        pregen = 1

        # ---- Find specified waveform type.
        hrssp = parameters[0]
        hrssc = parameters[1]
        name  = parameters[2]

        # ---- Read waveform.  Each has its own defined hrss for each
        #      polarisation.  First check for consistent hrss values.
        # Read hplus
        hp = hrssp * TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                                     path='/{0}/{1}/hp'.format(family, name))


        hc = hrssc * TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                                     path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value


    elif family.lower() in ['o1snews']:

        # ----  Entry from one of various waveform catalogs (pregenerated).
        #       The waveforms in this catalog are of different lengths, but
        #       have the same sampling rate.
        pregen = 1

        # ---- Find specified waveform type.
        hrssp = parameters[0]
        hrssc = parameters[1]
        name  = parameters[2]

        # ---- Read waveform.  Each has its own defined hrss for each
        #      polarisation.  First check for consistent hrss values.
        # Read hplus
        hp = hrssp * TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                                     path='/{0}/{1}/hp'.format(family, name))


        hc = hrssc * TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                                     path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

    elif family.lower() in ['ott09']:

        # ----  Marek 2009 supernova matter waveform (pregenerated).
        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        hp = 10. / distance * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)


    elif family.lower() in ['ott10']:

        # ----  Dimmelmeier, Font, and Mueller supernova waveform (pregenerated).
        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        hc = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value


        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        hp = 10. / distance * hp
        hc = 10. / distance * hc

    elif family.lower() in ['pca']:

        # ----  Marek 2009 supernova matter waveform (pregenerated).
        hrss = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        # ---- Read waveform, which is defined at a range of 10kpc.
        #      Make sure h is same type of vector (column) as t.
        hp = 10. / distance * hp
        hc = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)

    elif family.lower() in ['pmns_pca']:


        distance = parameters[0]
        name = parameters[1]
        pregen = 1

        # Read hplus
        hp = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hp'.format(family, name))


        hc = TimeSeries.read(os.path.join(filedir, family + '.hdf5'),
                             path='/{0}/{1}/hc'.format(family, name))

        pregen_fs = hp.sample_rate.value
        pregen_T = hp.duration.value

        hp = hp / (distance/20.)
        hc = hc / (distance/20.)

    elif family.lower() in ['scalarchirplet']:

        h_rss = parameters[0]
        tau = parameters[1]
        f0 = parameters[2]

        # ---- Optional parameters.
        alpha = 0
        delta = 0
        ciota = 1
        if len(parameters) >= 4:
            alpha = parameters[3]

        if len(parameters) >= 5:
            delta = parameters[4]

        if len(parameters) >= 6:
            ciota = parameters[5]

        # ---- Waveform.
        h = h_rss * (
                     numpy.exp(
                              (-1 + sqrt_of_neg_1 * alpha) * (t-T0)**2 /
                              (4 * tau**2)
                              + sqrt_of_neg_1 * 2 * numpy.pi * (t-T0) *f0
                              + sqrt_of_neg_1 * delta)
                     /
                     (2 * math.pi * tau**2)**(1/4))

        hb = TimeSeries(h, dx=1./fs, name=family.lower())

        hc = TimeSeries(numpy.zeros(hb.size), dx=1/hp.sample_rate.value,
                        name=hb.name)
        hp = TimeSeries(numpy.zeros(hb.size), dx=1/hp.sample_rate.value,
                        name=hb.name)

        # ---- Turn off default interpolation (ad hoc waveform is designed
        #      to produce desired T0).
        pregen = 0

        """
        elif family.lower() in ['scalarsn']:
                # ---- Ibanez and Novak scalar-mode gravitational
                #      collapse waveform (pregenerated).

                # ---- Model type and parameters.
                model    = parameters[0]
                distance = parameters[1]
                newalpha    = parameters[2]

                # ---- Pre-sampled waveform.  Columns: time [s], strain [unitless].
                #      For model details on each model look Tabel 1 of
                #      http://stacks.iop.org/0004-637X/533/i=1/a=392.
                if strcmpi(model,'A')
                    alphadefault = 0.01
                    data = [ ...
                        2.0038910505836576 -92.51101321585884...
                        2.953307392996109 -88.10572687224658...
                        3.1089494163424125 -88.10572687224658...
                        3.264591439688716 -79.29515418502183...
                        3.3813229571984436 -79.29515418502183...
                        3.498054474708171 -66.07929515418482...
                        3.583657587548638 -52.86343612334804...
                        3.6692607003891053 -35.24229074889854...
                        3.7237354085603114 -17.62114537444927...
                        3.7782101167315174 4.405286343612488...
                        3.8326848249027234 35.24229074889877...
                        3.8949416342412455 79.29515418502206...
                        3.9338521400778212 114.53744493392082...
                        3.9961089494163424 189.4273127753304...
                        4.073929961089494 308.3700440528635...
                        4.143968871595331 449.3392070484581...
                        4.198443579766537 607.9295154185022...
                        4.276264591439689 828.193832599119...
                        4.33852140077821 1026.431718061674...
                        4.392996108949417 1149.7797356828194...
                        4.439688715953308 1220.2643171806167...
                        4.470817120622568 1251.101321585903...
                        4.517509727626459 1277.533039647577...
                        4.571984435797665 1290.748898678414...
                        4.626459143968872 1290.748898678414...
                        4.696498054474708 1281.9383259911895...
                        4.735408560311284 1264.3171806167402...
                        4.797665369649805 1255.5066079295154...
                        4.875486381322958 1246.6960352422907...
                        4.976653696498055 1246.6960352422907...
                        5.085603112840467 1246.6960352422907...
                        5.186770428015564 1237.8854625550662...
                        5.373540856031129 1242.2907488986784...
                        5.575875486381324 1237.8854625550662...
                        5.747081712062257 1246.6960352422907...
                        5.8404669260700395 1242.2907488986784]

                elseif strcmpi(model,'B')
                    alphadefault = 0.004
                    data = [ ...
              TimeSeries          0.017094017094017144 509.0909090909091...
                        0.8974358974358977 509.0909090909091...
                        1.0683760683760686 509.0909090909091...
                        1.2136752136752138 509.0909090909091...
                        1.3162393162393164 506.06060606060606...
                        1.4615384615384617 503.030303030303...
                        1.5384615384615388 500.0...
                        1.6495726495726495 490.9090909090909...
                        1.717948717948718 481.8181818181818...
                        1.7863247863247866 469.69696969696975...
                        1.8461538461538463 451.5151515151515...
                        1.9145299145299148 421.21212121212125...
                        1.9572649572649574 396.969696969697...
                        2.0085470085470085 363.6363636363636...
                        2.0512820512820515 306.06060606060606...
                        2.0854700854700856 263.6363636363636...
                        2.1196581196581197 203.03030303030306...
                        2.1623931623931627 133.33333333333337...
                        2.1880341880341883 81.81818181818187...
                        2.213675213675214 30.30303030303037...
                        2.247863247863248 -15.151515151515127...
                        2.2820512820512824 -48.4848484848485...
                        2.3162393162393164 -63.636363636363626...
                        2.3504273504273505 -57.57575757575751...
                        2.3931623931623935 -39.39393939393938...
                        2.4444444444444446 -18.18181818181813...
                        2.4957264957264957 -6.0606060606060055...
                        2.5384615384615388 3.0303030303030027...
                        2.717948717948718 3.0303030303030027...
                        2.7521367521367526 3.0303030303030027]

                elseif strcmpi(model,'C')
                    alphadefault = 0.01
                    data = [...
                        69.2929292929293 21.333333333333336  ...
                        69.4949494949495 21.333333333333336  ...
                        69.72222222222223 21.333333333333336  ...
                        69.97474747474747 21.333333333333336  ...
                        70.25252525252526 21.333333333333336  ...
                        70.5050505050505 21.393939393939394  ...
                        70.8838383838384 21.393939393939394  ...
                        71.23737373737374 21.393939393939394  ...
                        71.74242424242425 21.393939393939394  ...
                        72.2979797979798 21.515151515151516  ...
                        72.77777777777779 21.515151515151516  ...
                        73.15656565656566 21.575757575757578  ...
                        73.66161616161617 21.636363636363637  ...
                        74.11616161616162 21.636363636363637  ...
                        74.64646464646465 21.6969696969697  ...
                        75.22727272727273 21.75757575757576  ...
                        75.85858585858587 21.81818181818182  ...
                        76.5909090909091 21.93939393939394  ...
                        77.52525252525253 22.060606060606062  ...
                        78.51010101010101 22.303030303030305  ...
                        79.57070707070707 22.666666666666668  ...
                        80.42929292929293 23.272727272727273  ...
                        80.83333333333334 23.6969696969697  ...
                        81.03535353535354 24.121212121212125  ...
                        81.23737373737374 24.72727272727273  ...
                        81.43939393939394 25.63636363636364  ...
                        81.64141414141415 27.090909090909093  ...
                        81.76767676767678 28.969696969696972  ...
                        81.86868686868688 31.212121212121215  ...
                        81.96969696969697 33.45454545454545  ...
                        81.9949494949495 36.0  ...
                        82.04545454545455 39.333333333333336  ...
                        82.0959595959596 42.24242424242424  ...
                        82.12121212121212 43.27272727272727  ...
                        82.24747474747475 46.909090909090914  ...
                        82.42424242424244 42.54545454545455  ...
                        82.57575757575758 40.0  ...
                        82.70202020202021 38.78787878787879  ...
                        82.85353535353536 39.57575757575758  ...
                        82.95454545454547 39.93939393939394  ...
                        83.13131313131314 39.45454545454545  ...
                        83.28282828282829 39.21212121212122  ...
                        83.56060606060606 38.484848484848484  ...
                        83.7878787878788 38.66666666666667  ...
                        83.88888888888889 38.60606060606061  ...
                        84.04040404040404 38.121212121212125  ...
                        84.1919191919192 38.0  ...
                        84.31818181818183 37.87878787878788  ...
                        84.4949494949495 37.81818181818182  ...
                        84.62121212121212 37.93939393939394  ...
                        84.84848484848486 37.75757575757576  ...
                        85.02525252525253 37.6969696969697  ...
                        85.12626262626263 37.75757575757576  ...
                        85.30303030303031 37.515151515151516  ...
                        85.8838383838384 37.3939393939394  ...
                        86.41414141414143 37.151515151515156  ...
                        86.61616161616162 37.27272727272727  ...
                        86.84343434343435 36.96969696969697  ...
                        87.14646464646465 37.151515151515156  ...
                        87.2979797979798 37.151515151515156  ...
                        87.57575757575758 37.151515151515156  ...
                        87.87878787878788 36.909090909090914]

                elseif strcmpi(model,'D')
                    alphadefault = 0.01
                    data = [ ...
                        76.26666666666667 21.627906976744185  ...
                        78.26666666666667 22.093023255813947  ...
                        79.48333333333333 22.558139534883722  ...
                        80.31666666666666 23.023255813953483  ...
                        80.75 23.720930232558132  ...
                        81.01666666666667 24.186046511627907  ...
                        81.26666666666667 24.65116279069767  ...
                        81.5 25.581395348837205  ...
                        81.7 27.441860465116278  ...
                        81.81666666666666 30.0  ...
                        81.9 33.25581395348837  ...
                        81.98333333333333 40.69767441860465  ...
                        82.06666666666666 52.790697674418595  ...
                        82.13333333333334 41.16279069767441  ...
                        82.15 31.860465116279073  ...
                        82.18333333333334 23.720930232558132  ...
                        82.25 30.69767441860465 ...
                        82.28333333333333 35.81395348837209  ...
                        82.31666666666666 44.18604651162791  ...
                        82.33333333333334 52.09302325581395  ...
                        82.38333333333334 63.72093023255813  ...
                        82.43333333333334 73.02325581395348  ...
                        82.5 80.69767441860465  ...
                        82.56666666666666 90.93023255813952  ...
                        82.66666666666667 101.62790697674419  ...
                        82.73333333333333 108.13953488372093  ...
                        82.81666666666666 112.55813953488372  ...
                        82.88333333333334 115.81395348837209  ...
                        82.95 118.37209302325581  ...
                        83.1 120.23255813953487  ...
                        83.23333333333333 120.93023255813954  ...
                        83.51666666666667 119.76744186046511  ...
                        83.71666666666667 118.13953488372093  ...
                        83.95 116.9767441860465  ...
                        84.23333333333333 115.81395348837209  ...
                        84.5 115.34883720930232  ...
                        84.9 113.72093023255813  ...
                        85.38333333333334 112.55813953488372  ...
                        85.88333333333334 111.62790697674419  ...
                        86.35 110.93023255813952  ...
                        86.66666666666667 110.69767441860465  ...
                        86.96666666666667 110.23255813953489  ...
                        87.13333333333334 110.23255813953489  ...
                        87.71666666666667 110.0  ...
                        88.15 109.30232558139534  ...
                        88.43333333333334 109.30232558139534]

                elseif strcmpi(model,'E')
                    alphadefault = 0.001
                    data = [ ...
                        69.33920704845815 0.0  ...
                        71.01321585903084 0.0  ...
                        72.99559471365639 0.0  ...
                        75.0 0.0  ...
                        76.98237885462555 0.0  ...
                        78.98678414096916 0.0  ...
                        80.99118942731278 0.0  ...
                        81.87224669603525 3.9603960396040065 ...
                        81.98237885462555 5.940594059405953  ...
                        82.0704845814978 13.861386138613852  ...
                        82.11453744493393 25.74257425742576  ...
                        82.18061674008811 39.60396039603961  ...
                        82.31277533039648 99.00990099009903  ...
                        82.5330396475771 215.8415841584159  ...
                        82.73127753303964 314.8514851485148  ...
                        82.9295154185022 481.1881188118812  ...
                        83.10572687224669 671.2871287128713  ...
                        83.17180616740089 752.4752475247525  ...
                        83.23788546255507 807.9207920792079  ...
                        83.32599118942731 829.7029702970298  ...
                        83.50220264317181 819.8019801980198  ...
                        83.6784140969163 841.5841584158416  ...
                        83.89867841409692 819.8019801980198  ...
                        84.11894273127754 833.6633663366337  ...
                        84.33920704845815 819.8019801980198  ...
                        84.49339207048459 827.7227722772277  ...
                        84.75770925110132 817.8217821782179 ...
                        84.86784140969164 825.7425742574258  ...
                        85.11013215859032 817.8217821782179  ...
                        85.48458149779736 817.8217821782179]

                else

                    error(['Model ' model ' unknown.'])

                end

            # ---- Catalog time is in milliseconds convert to seconds.
            time = data(:,1)/1000
            # Normalize to 1Kpc then let user specify distance (in Kpc)
            # Also multiply the strain by 1 over user specified omega coupling
            # paramters
            strain = data(:,2)*(1000/megaparsec_meters)*(1/distance)*...
            (newalpha/alphadefault)

            # ---- Attached smooth tapers to start and end.
            pregen = 1
            pregen_fs = 16384
            pregen_T = 1.0
            Nsamp = round(pregen_fs * pregen_T)
            t_all = [0:(Nsamp-1)].T * 1/pregen_fs
            time =  time + (0.5*pregen_T-time(1))
            # ---- Time samples before, during, and after pre-generated waveform.
            index_start = find(t_all<time(1))
            index_mid   = find(t_all>=time(1) & t_all<=time(end))
            index_end   = find(t_all>time(end))

            Fc = 1/pregen_T    # Hertz of windowing
            window = zeros(size(t_all))
            # ---- Interpolate pre-generated waveform to regularly sampled time
            #      samples.
            window(index_mid) = interp1(time,strain,t_all(index_mid))
            # ---- Attach smoothly varying sinusoid to start to avoid discontinuities.
            window(index_start) = (1 - cos(2*pi*Fc*t_all(index_start)))
            window(index_start) = window(index_start) * (window(index_mid(1))/window(index_start(end)))
            # ---- Attach smoothly varying sinusoid to end to avoid discontinuities.
            window(index_end) = (1 - cos(2*pi*Fc*t_all(index_end)))
            window(index_end) = window(index_end) * (window(index_mid(end))/window(index_end(1)))

            hb = window
            t_resamp = t_all
            hp = zeros(size(hb))
            hc = zeros(size(hb))

        elif family.lower() in ['scheide', 'scheidp']:

                # ---- This is Scheidegger 2010 waveform for equatorial (e) or
                #      polar (p) observer.  Read pregenerated supernova catalog.
                pregen = 1
                pregen_fs = 8192
                pregen_T = 0.5

                # ---- Load waveform structure, if not already loaded or supplied.
                if (isempty(catalog))
                    load([filedir 'Scheidcat.mat'])
                end

                # ---- Find specified waveform type.
                NAME = parameters{1,2}
                for k=1:length(Scheidcat)
                    if (strcmp(Scheidcat(k).name,NAME))
                        index = k
                    end
                end

                # ---- Read waveform, which is defined at a range of 10kpc.
                #      Make sure h is same type of vector (column) as t.
                switch lower(type)
                    case 'scheide'
                        hp = 10/parameters{1} * Scheidcat(index).hpe
                        hc = 10/parameters{1} * Scheidcat(index).hce
                    case 'scheidp'
                        hp = 10/parameters{1} * Scheidcat(index).hpp
                        hc = 10/parameters{1} * Scheidcat(index).hcp
                end

        elif family.lower() in ['sg']:

                # ---- sine-Gaussians
                h_rss = parameters(1)
                Q = parameters(2)
                f0 = parameters(3)
                h_peak = h_rss * ((4 * pi^0.5 * f0)/(Q * (1-exp(-Q^2))))^0.5
                tau = Q / (2^0.5 * pi * f0)
                hp = h_peak*sin(2*pi*(t-T0)*f0).*exp(-(t-T0).^2./tau.^2)
                hc = zeros(size(hp))

                # ---- Turn off default interpolation (symmetric ad hoc waveform).
                pregen = 0

        elif family.lower() in ['snmp']:

                # ----  Entry from one of various waveform catalogs (pregenerated).
                #       The waveforms in this catalog are of different lengths, but
                #       have the same sampling rate.
                pregen = 1
                pregen_fs = 16384

                # ---- Load waveform structure, if not already loaded or supplied.
                if (isempty(catalog))
                    catalog = load([filedir lower(type) '.mat'])
                end

                # ---- Find specified waveform type.
                hrssp = parameters{1,1}
                hrssc = parameters{1,2}
                NAME  = parameters{1,3}
                index = find(strcmpi(catalog.name,NAME))

                # ---- Read waveform.  Each has its own defined hrss for each
                #      polarisation.  First check for consistent hrss values.
                if (catalog.hrssp(index)==0 & hrssp~=0)
                    error('Requested waveform has zero amplitude for plus polarisation.')
                end
                if (catalog.hrssc(index)==0 & hrssc~=0)
                    error('Requested waveform has zero amplitude for cross polarisation.')
                end
                if (catalog.hrssp(index)>0)
                    hp = (hrssp/catalog.hrssp(index))*catalog.hp{index}
                else
                    hp = catalog.hp{index}  #-- zero waveform
                end
                if (catalog.hrssc(index)>0)
                    hc = (hrssc/catalog.hrssc(index))*catalog.hc{index}
                else
                    hc = catalog.hc{index}  #-- zero waveform
                end

                pregen_T = length(hp)/pregen_fs

        elif family.lower() in ['snn94']:

                # ---- Shibata, Nakao, and Nakamura scalar-mode gravitational
                #      collapse waveform (pregenerated).

                # ---- Model type and parameters.
                model = parameters{1}
                mass = parameters{2}
                distance = parameters{3}
                alphanew = parameters{4}

                # ---- Pre-sampled waveform.  Columns: time [s], strain [unitless].
                #      This data if for a 10 Msun BH at 10 mpc.
                if strcmpi(model,'a1')
                    alphadefault = 0.0316
                    data = np.array([ ...
                        0, 0,
                        0.00327733, 0  ...
                        0.00579835, 0  ...
                        0.00593841, 3.45745e-24  ...
                        0.00613449, 1.40957e-23  ...
                        0.00630256, 2.73936e-23  ...
                        0.00644261, 3.96277e-23  ...
                        0.00661068, 5.50532e-23  ...
                        0.00672273, 6.70213e-23  ...
                        0.00686278, 7.84574e-23  ...
                        0.00705886, 9.09574e-23  ...
                        0.00728295, 9.94681e-23  ...
                        0.007395, 1.01596e-22  ...
                        0.00778716 9.9734e-23  ...
                        0.00820733 9.78723e-23  ...
                        0.00871153 9.49468e-23  ...
                        0.00927176 9.06915e-23  ...
                        0.00974795 8.51064e-23  ...
                        0.0100561 7.95213e-23  ...
                        0.0102522 7.20745e-23  ...
                        0.0104482 6.2234e-23  ...
                        0.0106163 5.13298e-23  ...
                        0.0107564 4.12234e-23  ...
                        0.0108404 3.16489e-23  ...
                        0.0109524 2.12766e-23  ...
                        0.0110365 1.2234e-23  ...
                        0.0111485 4.78723e-24  ...
                        0.0112606 -3.19149e-24  ...
                        0.0113446 -8.51064e-24  ...
                        0.0115127 -1.30319e-23  ...
                        0.0115967 -1.43617e-23  ...
                        0.0118208 -1.38298e-23  ...
                        0.0120449 -1.03723e-23  ...
                        0.012325 -6.38298e-24  ...
                        0.0125771 -3.7234e-24  ...
                        0.0129132 -1.8617e-24  ...
                        0.0132214 -1.06383e-24  ...
                        0.0137536 -1.06383e-24  ...
                        0.0144259 -7.97872e-25  ...
                        0.0154343 -2.65957e-25  ...
                        0.0169749 0  ...
                        0.0194679 0 ]

                else

                    error(['Model ' model ' unknown.'])

                end

                # ---- Rescale waveform by user-specified distance and system mass.
                t_pregen = data(:,1) * mass / 10
                h_pregen = data(:,2) * (mass / distance) * (alphanew/alphadefault)

                # ---- Set first time stamp to zero, and resample to regular
                #      timestamps.
                t_pregen = t_pregen - t_pregen(1)
                pregen = 1
                pregen_fs = 16384
                t_resamp = [0:1/pregen_fs:t_pregen(end)].T
                h_resamp = interp1(t_pregen,h_pregen,t_resamp)
                pregen_T = length(h_resamp)/pregen_fs
                #
                hb = h_resamp
                hp = zeros(size(hb))
                hc = zeros(size(hb))

        elif family.lower() in ['stamp_pt_a_tapered','stamp_adi_b_tapered','stamp_adi_e_tapered','stamp_pt_b_tapered','stamp_adi_c_tapered','stamp_adi_a_tapered','stamp_adi_d_tapered','stamp_maxgnetarf_tapered','stamp_maxgnetarg_tapered']:

                # ---- Load waveform structure, if not already loaded or supplied.
                if (isempty(catalog))
                    disp(['Loading ' lower(type) '.mat catalog.'])
                    data = load([filedir lower(type) '.mat'])
                end

                # ---- Parse parameters.
                distance = parameters(1)  #-- mpc
                ciota    = parameters(2)  #-- inclination

                # ---- Extract +/x polarizations. Default distance is 1 mpc. Assume
                #      l=m=2 mode dominates when applying inclination.
                hp = 0.5*(1+ciota^2) * data.hp / distance
                hc =           ciota * data.hc / distance

                # ----  Set "pregenerated" variables.
                pregen = 1
                pregen_fs = data.fs
                pregen_T = length(hp)/data.fs

        elif family.lower() in ['stamp_monoa','stamp_quadb','stamp_linea','stamp_monoc','stamp_sga','stamp_wnba','stamp_wnbc','stamp_lineb','stamp_quada','stamp_sgc','stamp_wnbb']:

                # ---- Load waveform structure, if not already loaded or supplied.
                if (isempty(catalog))
                    disp(['Loading ' lower(type) '.mat catalog.'])
                    data = load([filedir lower(type) '.mat'])
                end

                # ---- Parse parameters. Note: don't call the requested hrss "hrss"
                #      because this name will trigger the hrss rescaling code at
                #      the end of this function.
                hrss_requested = parameters(1)  #-- mpc
                ciota          = parameters(2)  #-- inclination

                # ---- Extract +/x polarizations. Default hrss for STAMP-AS
                #      injection files is 1e-20. Assume l=m=2 mode dominates
                #      when applying inclination.
                hp = 0.5*(1+ciota^2) * data.hp * hrss_requested/1e-20
                hc =           ciota * data.hc * hrss_requested/1e-20

                # ----  Set "pregenerated" variables.
                pregen = 1
                pregen_fs = data.fs
                pregen_T = length(hp)/data.fs

        elif family.lower() in ['wnb']:

                # ---- Turn off default interpolation.
                #      (Might actually be useful here, but don't want to alter
                #      frequency content of noise by low-pass filtering).
                pregen = 0

                # ---- Gaussian-modulated noise burst, white over specified band.
                h_rss = parameters(1)
                fc = parameters(2)
                df = parameters(3)
                dt = parameters(4)
                # ---- If fifth parameter is specified, use it to set the seed for
                #      the random number generator.
                #      KLUDGE: injections will be coherent only if all the
                #      detectors have the same sampling frequency.
                if(length(parameters)>=5)
                    randn('state',parameters(5))
                end
                # ---- Gaussian envelope
                env = exp(-(t-T0).^2/2/dt.^2)
                # ---- Band-limited noise (independent for each polarization)
                x = BLWNB(max(fc-df,0),2*df,T,fs)
                x = x.T
                x = circshift(x,round((T0-1.5)*fs))
                hp = env.*x
                hp = hp*h_rss/(hp.T*hp/fs).^0.5
                x = BLWNB(max(fc-df,0),2*df,T,fs)
                x = x.T
                x = circshift(x,round((T0-1.5)*fs))
                hc = env.*x
                hc = hc*h_rss/(hc.T*hc/fs).^0.5

        elif family.lower() in ['yakunin2010']:

                # ----  long-lived bar mode supernova waveform (pregenerated).
                pregen = 1
                pregen_fs = 16384
                pregen_T = 1

                # ---- Load waveform structure, if not already loaded or supplied.
                if (isempty(catalog))
                    load([filedir 'Yakunin2010cat.mat'])
                end

                # ---- Find specified waveform type.
                NAME = parameters{1,2}
                for k=1:length(Yakunin2010cat)
                    if (strcmp(Yakunin2010cat(k).name,NAME))
                        wf_number = k
                    end
                end

                # ---- Read waveform, which is defined at a range of 10kpc.
                #      Make sure h is same type of vector (column) as t.
                distance = 0.1*parameters{1}
                # note: if these are read in from column ascii files, no need to invert vector dims
                hp = 1./distance*Yakunin2010cat(wf_number).hp
                hc = 1./distance*Yakunin2010cat(wf_number).hc

        elif family.lower() in ['zero']:

                # ---- Trivial (null) signal.  Useful to avoid scripts crashing.
                hp = zeros(size(t))
                hc = zeros(size(t))

                # ---- Turn off default interpolation.
                pregen = 0

        elif family.lower() in ['zm']:

                # ---- Zwerger-Mueller supernova waveform (pregenerated).
                pregen = 1
                pregen_fs = 16384
                pregen_T = 1

                # ---- Load waveform structure, if not already loaded or supplied.
                if (isempty(catalog))
                    # disp('Loading Zwerger-Muller catalog.')
                    load([filedir 'ZMcat.mat'])
                end

                # ---- Find specified waveform type.
                NAME = parameters{1,2}
                for k=1:length(ZMcat)
                    if (strcmp(ZMcat(k).name,NAME))
                        ZM_number = k
                    end
                end

                # ---- Read waveform, which is defined at a range of 1mpc.
                #      Make sure h is same type of vector (column) as t.
                distance = 1e-3*parameters{1}
                h = 1./distance*ZMcat(ZM_number).hoft.T
                hp = h(:)
                hc = zeros(size(hp))

        """
    else:
       raise ValueError("Waveform family not recognized")

    ###########################################################################
    #   Process pre-generated waveforms.
    ###########################################################################

    # ---- Ensure breathing mode output is defined, even if not requested.
    #      (The interpolation and normalisation code that follows is simpler
    #      if this variable is defined.)
    if 'hb' not in locals():
        hb = TimeSeries(numpy.zeros(hp.size), dx=1/hp.sample_rate.value,
                        name=hp.name)

    # ---- Resample, time-shift, truncate, and/or zero-pad pregenerated
    #      waveform as needed.
    # ---- KLUDGE: We should issue warnings if any significant amount of the
    #      waveform (e.g., of hrss^2) is lost due to the shifting.
    if pregen:

        # ---- Time vector for pregenerated waveform.
        pregen_t = numpy.arange(0, pregen_T, 1./pregen_fs)
        if ((pregen_t.size != hp.size) or (pregen_t.size != hc.size) or
           (pregen_t.size != hb.size)):
            raise ValueError('Length of pregen_t '
                             'does not equal one of hb, hp, hc')

        # ---- Measure peak/characteristic time of pregenerated waveform.
        [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = xoptimalsnr(
            h = numpy.array([hp.value, hc.value, hb.value]), t0=pregen_t[0],
            fs=pregen_fs, S=None, F0=None, dF=None,
            Fmin=1./pregen_T, Fmax=0.5*pregen_fs-1./pregen_T)

        # ---- Shift pregenerated time vector by required amount.
        pregen_t = pregen_t + (T0 - Tchar)
        # ---- Interpolate to correct sampling rate and peak time.
        # ---- Time in column vector.
        t = numpy.arange(0, T, 1./fs)
        # ---- Find desired times which overlap pregen_t.
        k = numpy.where((t>=pregen_t[0]) & (t<=pregen_t[-1]))[0]
        # ---- Interpolate, with zero padding if pregenerated waveform is too
        #      short.  (Truncation of long waveforms handled automatically by
        #      specifying vector "t".)  Use spline rather than linear
        #      interpolation the latter effectively acts as a low-pass-filter
        #      which causes problems with high-frequency waveforms.
        hp_interp = TimeSeries(numpy.zeros(t.size), dx=1/fs,
                               name=hp.name)
        hc_interp = TimeSeries(numpy.zeros(t.size), dx=1/fs,
                               name=hp.name)
        hb_interp = TimeSeries(numpy.zeros(t.size), dx=1/fs,
                               name=hp.name)

        hp_interp_function = CubicSpline(pregen_t, hp)
        hc_interp_function = CubicSpline(pregen_t, hc)
        hb_interp_function = CubicSpline(pregen_t, hb)

        hp_interp[k] = hp_interp_function(t[k])
        hc_interp[k] = hc_interp_function(t[k])
        hb_interp[k] = hb_interp_function(t[k])

        hp = hp_interp
        hc = hc_interp
        hb = hb_interp

    # ---- Reset hrss amplitude to specified value, if any.
    #      For 2-polarization waveforms use hp^2+hc^2.
    #      Note: If the time-shifting moves a pre-generated waveform out of the
    #      interval of interest then this will screw up the amplitude!
    #if ( (hrss>0) and (hrss)):
    #    normalization = ((hp*hp + hc*hc + hb*hb) / fs)**0.5
    #    hp = hp * hrss/normalization
    #    hc = hc * hrss/normalization
    #    hb = hb * hrss/normalization

    return t, hp, hc, hp 
