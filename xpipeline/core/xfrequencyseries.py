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
import numpy
from collections import OrderedDict
from gwpy.frequencyseries import FrequencySeries
import copy

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['XFrequencySeriesDict',
           'XAntennaProjectedFrequencySeriesDict',
           'convert_to_dominant_polarization_frame']

SQRT_OF_NEG_1 = numpy.sqrt(numpy.array([-1],dtype=complex))

class XFrequencySeriesDict(OrderedDict):
    def project_onto_antenna_patterns(self, antenna_responses,
                                      to_dominant_polarization_frame=False,
                                      circular=True):
        """Shift timeseries by assosciated time delay

            Parameters
            ----------
            antenna_responses : `dict`
                key-wise pair of
                OrderedDict([('f_plus',
                              OrderedDict([('H1', array([-0.02424373])),
                                           ('L1', array([0.3089992]))])),
                             ('f_cross',
                              OrderedDict([('H1', array([-0.5677237])),
                                           ('L1', array([0.52872644]))])),
                             ('f_scalar',
                              OrderedDict([('H1', array([0.12427263])),
                                           ('L1', array([-0.30016348]))]))])

            to_dominant_polarization_frame " `bool`
                This boolean determines whether or not to calculate the
                relevant angle parameter that would project the data into
                the orthogonal cross plus polarization frame.
        """
        antenna_response_asds = XAntennaProjectedFrequencySeriesDict()
        for pattern, responses in antenna_responses.items():
            antenna_weighted_asds = XFrequencySeriesDict()
            for ifo, asd in self.items():
                abbr_ifo = ifo.split(':')[0]
                antenna_weighted_asds[ifo] = responses[abbr_ifo] / asd

            antenna_response_asds[pattern] = antenna_weighted_asds

        if to_dominant_polarization_frame:
            wfp = antenna_response_asds['f_plus'].to_array()
            wfc = antenna_response_asds['f_cross'].to_array()
            FpDP, FcDP, psi = convert_to_dominant_polarization_frame(wfp, wfc)
            tmp = copy.deepcopy(antenna_response_asds)
            for detidx, idet in enumerate(tmp['f_plus']):
                antenna_response_asds['f_plus'][idet] = (
                    numpy.cos(2*psi[:,detidx]) * tmp['f_plus'][idet] +
                    numpy.sin(2*psi[:,detidx]) * tmp['f_cross'][idet]
                    )
                antenna_response_asds['f_cross'][idet] = (
                    -numpy.sin(2*psi[:,detidx]) * tmp['f_plus'][idet] +
                    numpy.cos(2*psi[:,detidx]) * tmp['f_cross'][idet]
                    )

        if circular:
            antenna_response_asds.to_circular()

        return antenna_response_asds


    def to_array(self):
        """Convert to number of freq bins by number of detectors array
        """
        number_of_frequencies = list(self.values())[0].size
        number_of_detectors = len(self)
        array = numpy.zeros([number_of_frequencies, number_of_detectors])
        for idx, asd in enumerate(self.values()):
            array[:, idx] = asd

        return array


    def calculate_magnitude(self):
        """Matrix M_AB components.

           This is the dot product of the projected_asds, with
           itself, summed accross detectors.

           Returns:
               `gwpy.frequencyseries.FrequencySeries`
                   Units Hz
        """
        return numpy.sqrt(sum([v.real**2 + v.imag**2 for v in self.values()]))

    def slice_frequencies(self, frequencies):
        """select a subset of frequencies from XFrequencySeriesDict

           Parameters:
               indices (array):
                   an array of indexs to select from all elements
                   of `XFrequencySeriesDict`

           Returns:
               `XFrequencySeriesDict`
        """
        asd_subset = XFrequencySeriesDict()
        for det, asd in self.items():
            asd_subset[det] = self[det][numpy.in1d(self[det].xindex, frequencies)]

        return asd_subset

class XAntennaProjectedFrequencySeriesDict(OrderedDict):
    """Controls a dictionaries of projected asds
    """
    def calculate_magnitude(self):
        """Find unit vector of a series of asds

           Returns:
               `XAntennaProjectedFrequencySeriesDict`
        """
        asds_magnitude = XFrequencySeriesDict()
        for k, v in self.items():
            asds_magnitude[k] = v.calculate_magnitude()

        return asds_magnitude

    def to_unit(self):
        """Find unit vector of a series of asds

           Returns:
               `XAntennaProjectedFrequencySeriesDict`
        """
        for k, v in self.items():
            unit_dpf_asds = XFrequencySeriesDict()
            for k1, v1 in v.items():
                unit_dpf_asds[k1] = v1 / v.calculate_magnitude()
            self[k] = unit_dpf_asds

        return self

    def to_circular(self):
        """
        Find unit vector of a series of asds

           Returns:
               `XAntennaProjectedFrequencySeriesDict`
        """
        # Initialize projected circular keys
        self['f_left'] = XFrequencySeriesDict()
        self['f_right'] = XFrequencySeriesDict()
        self['f_left_null'] = XFrequencySeriesDict()
        self['f_right_null'] = XFrequencySeriesDict()

        # calculate the magntiude squared of the plus and cross projected asds
        f_plus_magnitude_squared = numpy.square(self['f_plus'].calculate_magnitude())
        f_cross_magnitude_squared = numpy.square(self['f_cross'].calculate_magnitude())

        for f_plus, f_cross in zip(self['f_plus'].values(), self['f_cross'].values()):
            self['f_right'][f_plus.name] = f_plus + SQRT_OF_NEG_1 * f_cross
            self['f_left'][f_plus.name] = f_plus - SQRT_OF_NEG_1 * f_cross
            self['f_right_null'][f_plus.name] = f_plus / f_plus_magnitude_squared - SQRT_OF_NEG_1 * f_cross / f_cross_magnitude_squared
            self['f_left_null'][f_plus.name] = f_plus / f_plus_magnitude_squared + SQRT_OF_NEG_1 * f_cross / f_cross_magnitude_squared

        return self

def convert_to_dominant_polarization_frame(Fp, Fc):
    """Take in stream of fplus and ≈f_cross and convert to DPF

       DPF is the Dominant Polarization Frame

        Parameters
        ----------
            Fp : `float`
            Fc :  `float`
    """
    # ---- Compute rotation needed to reach DP frame.
    psi = numpy.zeros([len(Fp), 1])
    psi[:, 0] = 1/4*numpy.arctan(2*(numpy.sum(Fp*Fc, 1))/(
                    numpy.sum(Fp*Fp, 1)-numpy.sum(Fc*Fc, 1))
                    )
    psi = psi.repeat(Fp.shape[1], 1)

    # ---- Rotate to DP frame.
    FpDP = numpy.cos(2*psi)*Fp + numpy.sin(2*psi)*Fc
    FcDP = -numpy.sin(2*psi)*Fp + numpy.cos(2*psi)*Fc

    # ---- Further rotate polarization by pi/4 if |Fp|<|Fc|.
    swapindex = (FpDP**2).sum(1) < (FcDP**2).sum(1)
    FpDP[swapindex, :] = FcDP[swapindex, :]
    FcDP[swapindex, :] = -FpDP[swapindex, :]
    psi[swapindex, :] = psi[swapindex, :] + numpy.pi/4

    return FpDP, FcDP, psi
