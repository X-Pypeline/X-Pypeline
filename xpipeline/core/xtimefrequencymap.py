#!/usr/bin/env python

# ---- Import standard modules to the python path.
from gwpy.spectrogram import Spectrogram
from collections import OrderedDict
import numpy as np


class XTimeFrequencyMapDict(OrderedDict):
    def to_coherent(self):
        """Sum all maps in the dict

           Returns:
               `XTimeFrequencyMap`:
                   A coherent TF-Map
        """
        coherent_map = 0
        for key, tfmap in self.items():
            coherent_map =+ tfmap

        return coherent_map


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
                projected_time_frequency_maps[pattern][det] = self[det] * asd
        return projected_time_frequency_maps


class XTimeFrequencyMap(Spectrogram):
    def find_significant_pixels(self, blackpixel_percentile=99):
        """Obtain the time-frequency indicies of the loudest pixels

           Parameters:

           blackpixel_percentile : `int`
               what 2D percentile value we will use as the threshold
        """
        energy_threshold  = np.percentile(self, blackpixel_percentile,
                                          interpolation='midpoint')
        pixel_time, pixel_freq = np.where(self.value > energy_threshold)

        return pixel_time, pixel_freq


    def phaseshift(self, delta):
        """Phase shift this `spectrogram` by ``delta``

        This modifies the spectrogram in-place.

        Parameters
        ----------
        delta : `float`, `~astropy.units.Quantity`, `str`
            The amount by which to shift (in seconds if `float`), give
            a negative value to shift backwards in time
        """
        sqrt_of_neg_1 = np.sqrt(np.array([-1], dtype=complex))
        frequency_shift = residual_time_shift(delta,
                                              self.frequencies.to_value())
        return self * frequency_shift


    def circular_time_slide(self, slide, sample_frequency, offset):
        """Phase shift this `spectrogram` by ``delta``

        This modifies the spectrogram in-place.

        Parameters
        ----------
        delta : `float`, `~astropy.units.Quantity`, `str`
            The amount by which to shift (in seconds if `float`), give
            a negative value to shift backwards in time
        """
        offsetlength = seconds_to_samples(offset, sample_frequency)
        ntimepixelshifted = slide * sample_frequency/offsetlength
        return np.roll(self.value, ntimepixelshifted)


    def to_dominant_polarization_frame(self, dpf_asd):
        return self * dpf_asd


def residual_time_shift(seconds, frequencies):
    # define sqrt of -1
    sqrt_of_neg_1 = np.sqrt(np.array([-1],dtype=complex))
    residual_time_shift_phases = np.exp(sqrt_of_neg_1 * 2 * np.pi * \
        frequencies * seconds)

    return residual_time_shift_phases
