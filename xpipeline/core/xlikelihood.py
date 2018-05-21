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
