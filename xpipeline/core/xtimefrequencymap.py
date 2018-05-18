#!/usr/bin/env python

# ---- Import standard modules to the python path.
from .xtimeseries import XTimeSeries
from .xdetector import Detector

class XTimeFrequencyMap(XTimeSeries):
    def find_significant_pixels(self, blackpixel_percentile=99):
        """White this `TimeSeries` against its own ASD

            Parameters
            ----------
            fft_length : `float`
                number of seconds in single FFT
        """
        signficant_pixels = {}
        for (idet, itfmap) in self.tfmaps.items():
            signficant_pixels[idet] = {}
            black_pixels_rows, black_pixels_columns = np.where(
                                    itfmap.value >
                             np.percentile(itfmap.value, blackpixel_percentile,
                                    interpolation='midpoint')
                                    )
            signficant_pixels[idet]['rows'] = black_pixels_rows
            signficant_pixels[idet]['columns'] = black_pixels_columns
        self.signficant_pixels = signficant_pixels


    def spectrogram(self, fftlength, overlap=0, window='hann'):
        """White this `TimeSeries` against its own ASD

            Parameters
            ----------
            fft_length : `float`
                number of seconds in single FFT
        """
        tfmaps = TimeSeriesDict()
        for (idet, iseries) in self.event.whitened_timeseries.items():
            tfmap = iseries.spectrogram(stride=fftlength,
                                fftlength=fftlength, overlap=overlap,
                                window=window)
            tfmaps.append(
                          {idet : tfmap}
                         )
        self.tfmaps = tfmaps
        return tfmaps
