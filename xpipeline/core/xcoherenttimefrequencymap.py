#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np
from .xtimefrequencymap import XTimeFrequencyMap


class XCoherentTimeFrequencyMap(XTimeFrequencyMap):
    def __init__(self, xtimeseries, fftlength, overlap=0, window='hann'):
        """This is an appropraitely shift tfmap

            Parameters
            ----------
            xtimeseries : `float`
                number of seconds in single FFT

            fftlength : `float`
                number of seconds in single FFT

            overlap : `float`
                number of seconds in single FFT

            window : `float`
                number of seconds in single FFT
        """

    def internal_time_slide(self, pixIdx, shift, tfMapLength):
        """This is an appropraitely shift tfmap

            Parameters
            ----------
            xtimeseries : `float`
                number of seconds in single FFT

            fftlength : `float`
                number of seconds in single FFT

            overlap : `float`
                number of seconds in single FFT

            window : `float`
                number of seconds in single FFT
        """
        # computes the new linear indices of pixels after they are shifted
        shiftedPixIdx = pixIdx  + shift;
        # wrap pixels that go beyond the left or right edge of the TF map
        shiftedPixIdx = shiftedPixIdx + \
            tfMapLength*((shiftedPixIdx<=0) - (shiftedPixIdx>tfMapLength))
        return shiftedPixIdx
