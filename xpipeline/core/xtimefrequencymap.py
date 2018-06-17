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
from collections import OrderedDict
from gwpy.signal.fft.ui import seconds_to_samples
from .xfrequencyseries import XFrequencySeriesDict
import numpy


__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['XTimeFrequencyMapDict', 'residual_time_shift', 'XTimeFrequencyMap']

class XTimeFrequencyMapDict(OrderedDict):
    def abs(self):
        """Take the absolute value of all maps in dict

           Returns:
               `XTimeFrequencyMapDict`:
                   power_map of all Fourier Grams in Dict
        """
        return XTimeFrequencyMapDict({k: v.abs() for k,v in self.items()})

    def to_coherent(self):
        """Sum all maps in the dict

           Returns:
               `XTimeFrequencyMap`:
                   A coherent TF-Map
        """
        return sum(self.values())

    def to_dominant_polarization_frame(self, projected_asds):
        """Project Tfmap to an antenna frame give a dict of asds

           Parameters:
               projected_asds : `dict`
                   key-wise dict of antenna pattern name
                   and values `XFrequencySeriesDict`

           Returns:
               `OrderedDict`:
                   A key-wise dict of antenna pattern name
                   and value `XTimeFrequencyMapDict`
        """
        projected_time_frequency_maps = OrderedDict()
        for pattern, asds in projected_asds.items():
            projected_time_frequency_maps[pattern] = XTimeFrequencyMapDict()
            for det, asd in asds.items():
                mask = numpy.in1d(asd.xindex, self[det].yindex)
                projected_time_frequency_maps[pattern][det] = self[det] * asd[mask]
        return projected_time_frequency_maps
        significant_pixels = OrderedDict()
        for key, pix_thres in blackpixel_percentile.items():
            significant_pixels[key] = {}
            significant_pixels[key]['pix_time'], significant_pixels[key]['pix_freq'] \
                = self[key].find_significant_pixels(blackpixel_percentile=pix_thres)

    def blackout_pixels(self, blackpixel_percentile):
        """Set pixels below certain energy level to zero

        Parameters:

            blackpixel_percentile : `dict`, `int`
                either a `dict` of (channel, `int`) pairs for key-wise
                significant pixel calc, or a single int to use as the
                threshold of all items.

        Returns:
            `dict`: key-wise pair of channel : freq,time indices
        """
        if not isinstance(blackpixel_percentile, dict):
            blackpixel_percentile = dict((c, blackpixel_percentile)
                                         for c in self)

        self_ = self
        for key, pix_thres in blackpixel_percentile.items():
            self_[key] = self_[key].blackout_pixels(
                                       blackpixel_percentile=pix_thres)

        return self_

    def plot(self, label='key', **kwargs):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        label : `str`, optional

            labelling system to use, or fixed label for all elements
            Special values include

            - ``'key'``: use the key of the `XTimeFrequencyMapDict`,
            - ``'name'``: use the :attr:`~XTimeSeries.name` of each element

            If anything else, that fixed label will be used for all lines.

        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        from gwpy.plotter import SpectrogramPlot
        figargs = dict()
        for key in ['figsize', 'dpi']:
            if key in kwargs:
                figargs[key] = kwargs.pop(key)
        plot_ = SpectrogramPlot(sep=True, **figargs)
        for lab, tfmap in self.items():
            if label.lower() == 'name':
                lab = tfmap.name
            elif label.lower() != 'key':
                lab = label
            plot_.add_spectrogram(tfmap, label=lab, newax=True, **kwargs)
        return plot_

    def gaussianity(self):
        """Calculate the gaussianity of all maps in this dict

           (.abs().percentile(99)/ .abs().median(0))**2
           Returns:
               `XFrequencySeriesDict`:
                   The gaussianity measure of each frequency bin
        """
        # make sure we are dealing with energy map not fftgram
        return XFrequencySeriesDict({k: v.gaussianity()
                                     for k,v in self.items()}
                                    )


class XTimeFrequencyMap(Spectrogram):
    def blackout_pixels(self, blackpixel_percentile):
        """Set pixels below certain energy level to zero

           Parameters:

           blackpixel_percentile : `int`
               the x-percentile loudest tile all pixels
               below which will be set to energy of 0
        """
        self_ = self
        energy_threshold  = numpy.percentile(self_, blackpixel_percentile,
                                             interpolation='midpoint')
        self_[self_.value <= energy_threshold] = 0
        return self_ 

    def gaussianity(self):
        """Calculate the gaussianity of this map

           (.abs().percentile(99)/ .abs().median(0))**2
           Returns:
               `gwpy.frequency.FrequencySeries`:
                   The gaussianity measure of each frequency bin
        """
        self_ = self.abs()
        return self_.percentile(99) / self_.median(0)

    def phaseshift(self, delta):
        """Phase shift this `spectrogram` by ``delta``

        This modifies the spectrogram in-place.

        Parameters
        ----------
        delta : `float`, `~astropy.units.Quantity`, `str`
            The amount by which to shift (in seconds if `float`), give
            a negative value to shift backwards in time
        """
        sqrt_of_neg_1 = numpy.sqrt(numpy.array([-1], dtype=complex))
        frequency_shift = residual_time_shift(delta,
                                              self.frequencies.to_value())
        return self * frequency_shift

    def circular_time_slide(self, seconds, sample_frequency, offset):
        """Slide the TF pixels of this map

        This should move the appropriate number of time bins
        such that slide represents a slide in seconds.

        Parameters
        ----------
        seconds : `int`,
            How many seconds we are sliding the map

        sample_frequency : `float`
            what is the sample frequency of the data

        offset : `float`
            what offset was used to make this spectrogram

        Returns:
            `XTimeFrequencyMap`:
                A time frequency map slide by the appropriate number
                of seconds
        """
        offsetlength = seconds_to_samples(offset, sample_frequency)
        ntimepixelshifted = int(seconds * sample_frequency//offsetlength)

        slided_map = numpy.roll(self.value, ntimepixelshifted)
        return XTimeFrequencyMap(slided_map,
                                 yindex=self.yindex, xindex=self.xindex)

    def to_dominant_polarization_frame(self, dpf_asd):
        return self * dpf_asd


def residual_time_shift(seconds, frequencies):
    # define sqrt of -1
    sqrt_of_neg_1 = numpy.sqrt(numpy.array([-1],dtype=complex))
    residual_time_shift_phases = numpy.exp(sqrt_of_neg_1 * 2 * numpy.pi *
                                           frequencies * seconds)

    return residual_time_shift_phases
