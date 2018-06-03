# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017-)
#
# This file is part of the XPypeline python package.
#
# hveto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hveto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hveto.  If not, see <http://www.gnu.org/licenses/>.
 
# ---- Import standard modules to the python path.
from gwpy.spectrogram import Spectrogram
import numpy as np

sqrt_of_neg_1 = np.sqrt(np.array([-1],dtype=complex))

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['XLikelihood']

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
            wfctimefrequencymap.imag**2
            - 2 * (wfptimefrequencymap.conj() *
                        wfctimefrequencymap).imag)

        likelihood1 = likelihood1 / (mpp + mcc)
        likelihood2 = (
            wfptimefrequencymap.real**2 +
            wfptimefrequencymap.imag**2
            + wfctimefrequencymap.real**2 +
            wfctimefrequencymap.imag**2
            + 2 * (wfptimefrequencymap.conj() *
                        wfctimefrequencymap).imag)

        likelihood2 = likelihood2 / (mpp+mcc)
        likelihood = np.maximum(likelihood1, likelihood2)

        return likelihood

    @classmethod
    def circinc(cls, tfmaps, mpp, mcc, projected_asds):
        """Circular incoherent likelihood detection statistic.

        Compute larger of left/rightcircularly
        polarized GWB.  (-) sign for last term SHOULD
        correspond to hp = cos(), hc = sin(), or inspiral with
        iota=0.

        Parameters:
            mpp (array):
                sum(asds_plus**2, 1)

            mcc (array):
                sum(asds_cross**2, 1)

            tfmaps (`XTimeFrequencyMapDict`):
                time-frequency map

            projected_asds (`XFrequencySeriesDict`):
                Dominant Polarization Frame projected asds

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        # ---- For a general network (aligned or not):
        likelihood = 0
        for idet in tfmaps:
            likelihood =+ tfmaps[idet] * (
                              (projected_asds['f_plus'][idet]**2 +
                               projected_asds['f_cross'][idet]**2/
                               (mpp + mcc)
                              ))

        return likelihood

    @classmethod
    def circnullenergy(cls, mpp, mcc, wfptimefrequencymap, wfctimefrequencymap):
        """Compute larger of left/right circularly polarized GWB.

        (-) sign for last term SHOULD
        correspond to hp = cos(), hc = sin(), or inspiral with iota=0.

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
        n_mcc = mcc / np.sqrt(mpp*mcc)
        n_mpp = mpp / np.sqrt(mpp*mcc)

        likelihood1 = (n_mcc * wfptimefrequencymap +
                       sqrt_of_neg_1 * n_mpp * wfctimefrequencymap).abs()**2/\
                      (mpp + mcc)

        likelihood2 = (n_mcc * wfptimefrequencymap -
                       sqrt_of_neg_1 * n_mpp * wfctimefrequencymap).abs()**2/\
                      (mpp + mcc)

        likelihood = np.minimum(likelihood1, likelihood2)

        return likelihood

    @classmethod
    def circnullinc(cls, tfmaps, mpp, mcc, projected_asds):
        """Incoherent circular null energy

        Parameters:
            mpp (array):
                sum(asds_plus**2, 1)

            mcc (array):
                sum(asds_cross**2, 1)

            tfmaps (`XTimeFrequencyMapDict`):
                time-frequency map

            projected_asds (`XFrequencySeriesDict`):
                Dominant Polarization Frame projected asds

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        likelihood = 0
        for idet in tfmaps:
            likelihood =+ tfmaps[idet] * (
                              ((mcc * projected_asds['f_plus'][idet])**2 +
                               (mpp * projected_asds['f_cross'][idet])**2/
                               ((mpp + mcc) * mpp * mcc)
                              ))

        return likelihood
