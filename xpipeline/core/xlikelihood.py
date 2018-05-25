#!/usr/bin/env python
  
# ---- Import standard modules to the python path.
from gwpy.spectrogram import Spectrogram
import numpy as np

class XLikelihood(object):
    @classmethod
    def standard(cls, mpp, mcc, wfptimefrequencymap, wfctimefrequencymap):
        """Standard likelihood detection statistic.

        Parameters:
            mpp (array):
                sum(asds_plus**2, 1)
            mcc (array):
                sum(asds_cross**2, 1)
            wfptimefrequencymap (array):
                plus weighted time-frequency-map
            wfctimefrequencymap (array):
                cross weighted time-frequency-map

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        likelihood = (
                     (mpp**(-1) * (wfptimefrequencymap.real**2 +
                                   wfptimefrequencymap.imag**2))
                     +
                     (mcc**(-1) * (wfctimefrequencymap.real**2 + 
                                   wfctimefrequencymap.imag**2)))
        return likelihood

    @classmethod
    def circenergy(cls, mpp, mcc, wfptimefrequencymap, wfctimefrequencymap):
        """Circular energy likelihood detection statistic.

        Compute larger of left/rightcircularly
        polarized GWB.  (-) sign for last term SHOULD
        correspond to hp = cos(), hc = sin(), or inspiral with
        iota=0.

        Parameters:
            mpp (array):
                sum(asds_plus**2, 1)
            mcc (array):
                sum(asds_cross**2, 1)
            wfptimefrequencymap (array):
                plus weighted time-frequency-map
            wfctimefrequencymap (array):
                cross weighted time-frequency-map

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        # ---- Compute larger of left/rightcircularly
        #      polarized GWB.  (-) sign for last term SHOULD
        #      correspond 
        #      to hp = cos(), hc = sin(), or inspiral with
        #      iota=0.
        likelihood1 = (
            wfptimefrequencymap.real**2 +
            wfptimefrequencymap.imag**2
            + wfctimefrequencymap.real**2 +
            wfctimefrequencymap.image**2
            - 2 * (wfptimefrequencymap.conj *
                        wfctimefrequencymap).imag)

        likelihood1 = likelihood1 / (Mpp+Mcc);
        likelihood2 = (
            wfptimefrequencymap.real**2 +
            wfptimefrequencymap.imag**2
            + wfctimefrequencymap.real**2 +
            wfctimefrequencymap.image**2
            + 2 * (wfptimefrequencymap.conj *
                        wfctimefrequencymap).imag)

        likelihood2 = likelihood2 / (Mpp+Mcc);
        likelihood=max(likelihood1,likelihood2);
        return likelihood

    @classmethod
    def circinc(cls, tfmaps, mpp, mcc, projected_asds):
        """Circular energy likelihood detection statistic.

        Compute larger of left/rightcircularly
        polarized GWB.  (-) sign for last term SHOULD
        correspond to hp = cos(), hc = sin(), or inspiral with
        iota=0.

        Parameters:
            mpp (array):
                sum(asds_plus**2, 1)
            mcc (array):
                sum(asds_cross**2, 1)
            wfptimefrequencymap (array):
                plus weighted time-frequency-map
            wfctimefrequencymap (array):
                cross weighted time-frequency-map

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        # ---- For a general network (aligned or not):
        for channelNumber in range(numberOfChannels):
            likelihood = (likelihood +
              (wFpDP[:,channelNumber]**2 + wFcDP[:,channelNumber]**2)
              / (Mpp + Mcc) *
               pixEnergyMap[channelNumber])

    @classmethod
    def circnullenergy(cls, mpp, mcc, projected_tfmaps):
        """
        % ---- Compute larger of left/right
        % circularly
        %      polarized GWB.  (-) sign for last term SHOULD
        %      correspond 
        %      to hp = cos(), hc = sin(), or inspiral with
        %      iota=0.
        """
        irat=sqrt(-1);
        nMcc = Mcc/sqrt(Mpp*Mcc)
        nMpp = Mpp/sqrt(Mpp*Mcc)
        likelihood1 = abs(nMcc*wfptimefrequencymap + irat*nMpp*wfctimefrequencymap)**2/\
            (Mpp+Mcc)
        likelihood2 = abs(nMcc*wfptimefrequencymap - irat*nMpp*wfctimefrequencymap)**2/\
            (Mpp+Mcc)
        likelihood = min(likelihood1,likelihood2)

    @classmethod
    def circnullinc(cls, tfmaps, mpp, mcc, projected_asds):
        """
         ---- For a general network (aligned or not):
        """
        for channelNumber in range(numberOfChannels):
          likelihood = likelihood + ( \
              ((Mcc*wFpDP[:,channelNumber])**2 + (Mpp*wFcDP[:,channelNumber])**2)\
              / ((Mpp + Mcc)*Mpp*Mcc)) * \
              pixEnergyMap[channelNumber]
