#!/usr/bin/env python
  
# ---- Import standard modules to the python path.
from gwpy.spectrogram import Spectrogram
import numpy as np

class XLikelihood(object):
    def standard(mpp, mcc, wfptimefrequencymap, wfctimefrequencymap):
        """Standard likelihood detection statistic.
            mpp :
                sum(asds_plus**2, 1)
            mcc :
                sum(asds_cross**2, 1)
            wfptimefrequencymap :
                plus weighted time-frequency-map
            wfctimefrequencymap :
                cross weighted time-frequency-map
        """
        likelihood = (
                     (mpp**(-1) * (wfptimefrequencymap.real**2 +
                                   wfptimefrequencymap.imag.^2))
                     +
                     (mcc**(-1) * (wfctimefrequencymap.real**2 + 
                                   wfctimefrequencymap.imag**2)
        return likelihood
