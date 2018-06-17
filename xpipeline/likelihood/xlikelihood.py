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
from ..core.xtimefrequencymap import XTimeFrequencyMap
import numpy
import math

sqrt_of_neg_1 = numpy.sqrt(numpy.array([-1], dtype=complex))

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['XLikelihood']

class XLikelihoodMap(object):
    def __init__(self, mpp, mcc, tfmaps, projected_asds):
        """Initiated a XLikelihoodMap Class

        This class allows for calculating the probability
        a given pixels in your tfmap is likely to be from a
        GW or not for a variety of likelihoods.

        Specifically, by initializing this class
        with `mpp` `mcc` `tfmaps` and `projected_asds`
        you are able to calculate the rest of the necessary
        variables needed to calculate most likelihoods.

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
            `XLikelihoodMap`
        """
        self.mpp = mpp
        self.mcc = mcc
        self.tfmaps = tfmaps
        self.energy_maps = tfmaps.abs()
        self.projected_asds = projected_asds

        # take the tfmaps and projected them using the projected_asds
        projected_tfmaps = tfmaps.to_dominant_polarization_frame(projected_asds)

        # here we pre calculate a number of the values that are needed
        # in the likelihoods below. This is because many of the likelihood
        # share similar calculations.
        self.wfptimefrequencymap = projected_tfmaps['f_plus'].to_coherent()

        self.wfctimefrequencymap = projected_tfmaps['f_cross'].to_coherent()

        self.wfptimefrequencymap_squared = self.wfptimefrequencymap.abs()**2

        self.wfctimefrequencymap_squared = self.wfctimefrequencymap.abs()**2

        # Here we create an empty time-frequency map for the
        # incoherent statistics 
        times = self.wfctimefrequencymap.xindex
        frequencies = self.wfctimefrequencymap.yindex
        map_shape = self.wfctimefrequencymap.shape

        self.empty_map = XTimeFrequencyMap(numpy.zeros(map_shape),
                                           yindex=frequencies, xindex=times)


    def standard(self):
        """Standard likelihood detection statistic.

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        return self.plusenergy() + self.crossenergy()

    def plusenergy(self):
        """Plusenergy likelihood detection statistic.

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        return self.wfptimefrequencymap_squared * self.mpp**(-1)

    def crossenergy(self):
        """Standard likelihood detection statistic.

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        return self.wfctimefrequencymap_squared * self.mcc**(-1)

    def plusinc(self):
        """Standard likelihood detection statistic.

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        inc_energy = self.empty_map

        for (detector, tfmap) in self.energy_maps.items():
            inc_energy += ((self.projected_asds['f_plus'][detector]**2 /
                           self.mcc) * tfmap)
        return inc_energy

    def crossinc(self):
        """Standard likelihood detection statistic.

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        inc_energy = self.empty_map

        for (detector, tfmap) in self.energy_maps.items():
            inc_energy += ((self.projected_asds['f_cross'][detector]**2 /
                           self.mcc) * tfmap)
        return inc_energy

    def loghbayesiancirc(self, strain_sigma):
        """circularly polarization detection statistic.

        For more details see
        `Equation 9 <https://journals.aps.org/prd/pdf/10.1103/PhysRevD.86.022003#temp%3Aintralink-d9>`_

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        circular_response = (self.mpp + self.mcc) /2

        return self.mcc**(-1) * self.wfctimefrequencymap_squared


    def powerlaw(self):
        """powerlaw detection statistic.

        Returns:
            `likelihood` : `array`
                The standard likelihood assosciated with your time-
                frequency map and/or sub clusters of tf-pixels
        """
        # -- If clusters already defined compute likelihood only
        #    for those pixels
        if len(tfmaps) == 2:
            likelihood = (
                pixEnergyMap[1]**(math.sqrt(11) / gaussianity[1](pixFreq)) *
                pixEnergyMap[2]**(math.sqrt(11) / gaussianity[2](pixFreq))
                )
        elif len(tfmaps) == 3:
            likelihood = (
                (pixEnergyMap[1]**(math.sqrt(11) / gaussianity[1](pixFreq))) *
                (pixEnergyMap[2]**(math.sqrt(11) / gaussianity[2](pixFreq))) *
                (pixEnergyMap[3]**(math.sqrt(11) / gaussianity[3](pixFreq)))
              )

            wNull = cross(wFpDP / repmat(numpy.sqrt(Mpp),[1, 3]), wFcDP /
                       repmat(numpy.sqrt(Mcc),[1, 3]),2)
            N1 = numpy.abs(wNull[:,1])
            N2 = numpy.abs(wNull[:,2])
            N3 = numpy.abs(wNull[:,3])
            # ---- Compute total energy
            nullamplitude = zeros(size(pixFreq));
            for channelNumber in tfmaps:
                nullamplitude = (nullamplitude
                    + pixEnergyMap[channelNumber])

            # ---- Subtract standard likelihood energy.
            nullamplitude = nullamplitude - self.standard
            nullamplitude = numpy.sqrt(nullamplitude)
            likelihood = (likelihood *
                (1. / (1 + numpy.abs(nullamplitude / N1)**(math.sqrt(44) /
                 gaussianity[1](pixFreq)))
                    +
                 1. / (1 + numpy.abs(nullamplitude / N2)**(math.sqrt(44) /
                 gaussianity[2](pixFreq)))
                    +
                 1. / (1 + numpy.abs(nullamplitude / N3)**(math.sqrt(44) /
                 gaussianity[3](pixFreq)))
                ))
            likelihood = numpy.log(likelihood)

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
            self.wfptimefrequencymap_squared +
            self.wfctimefrequencymap_squared -
            2 * (self.wfptimefrequencymap.conj() *
                 self.wfctimefrequencymap).imag)

        likelihood1 = likelihood1 / (self.mpp + self.mcc)

        likelihood2 = (
            self.wfptimefrequencymap_squared +
            self.wfctimefrequencymap_squared +
            2 * (self.wfptimefrequencymap.conj() *
                 self.wfctimefrequencymap).imag)

        likelihood2 = likelihood2 / (self.mpp + self.mcc)

        likelihood = numpy.maximum(likelihood1, likelihood2)

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
        n_mcc = mcc / numpy.sqrt(mpp*mcc)
        n_mpp = mpp / numpy.sqrt(mpp*mcc)

        likelihood1 = (n_mcc * wfptimefrequencymap +
                       sqrt_of_neg_1 * n_mpp * wfctimefrequencymap).abs()**2/\
                      (mpp + mcc)

        likelihood2 = (n_mcc * wfptimefrequencymap -
                       sqrt_of_neg_1 * n_mpp * wfctimefrequencymap).abs()**2/\
                      (mpp + mcc)

        likelihood = numpy.minimum(likelihood1, likelihood2)

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
