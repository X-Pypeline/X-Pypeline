#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np
from collections import OrderedDict
from gwpy.frequencyseries import FrequencySeries
import copy

class XFrequencySeriesDict(OrderedDict):
    def project_onto_antenna_patterns(self, antenna_responses,
                                      to_dominant_polarization_frame=False):
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
        antenna_response_asds = OrderedDict()
        for pattern, responses in antenna_responses.items():
            antenna_weighted_asds = XFrequencySeriesDict()
            for det, asd in self.items():
                abbr_det = det.split(':')[0]
                antenna_weighted_asds[det] = responses[abbr_det] / asd

            antenna_response_asds[pattern] = antenna_weighted_asds


        if to_dominant_polarization_frame:
            wfp = antenna_response_asds['f_plus'].to_array()
            wfc = antenna_response_asds['f_cross'].to_array()
            FpDP, FcDP, psi = convert_to_dominant_polarization_frame(wfp, wfc)
            tmp = copy.deepcopy(antenna_response_asds)
            for detidx, idet in enumerate(tmp['f_plus']):
                antenna_response_asds['f_plus'][idet] = (
                    np.cos(2*psi[:,detidx]) * tmp['f_plus'][idet] +
                    np.sin(2*psi[:,detidx]) * tmp['f_cross'][idet]
                    )
                antenna_response_asds['f_cross'][idet] = (
                    -np.sin(2*psi[:,detidx]) * tmp['f_plus'][idet] +
                    np.cos(2*psi[:,detidx]) * tmp['f_cross'][idet]
                    )
            

        return antenna_response_asds


    def to_array(self):
        """Convert to number of freq bins by number of detctors array
        """
        number_of_frequencies = list(self.values())[0].size
        number_of_detectors = len(self)
        array = np.zeros([number_of_frequencies, number_of_detectors])
        for idx, asd in enumerate(self.values()):
            array[:, idx] = asd

        return array


    def to_m_ab(self):
        """Matrix M_AB components.

           These are the dot products of wFp, with
           themselves and each other, for each frequency, computed in
           the DP frame.
        """
        return np.sum(self.to_array()**2, 1)


def convert_to_dominant_polarization_frame(Fp, Fc):
    """Take in stream of fplus and ≈f_cross and convert to DPF

       DPF is the Dominant Polarization Frame

        Parameters
        ----------
            Fp : `float`
            Fc :  `float`
    """
    # ---- Compute rotation needed to reach DP frame.
    psi = np.zeros([len(Fp), 1])
    psi[:, 0] = 1/4*np.arctan(2*(np.sum(Fp*Fc, 1))/(
                    np.sum(Fp*Fp, 1)-np.sum(Fc*Fc, 1))
                    )
    psi = psi.repeat(Fp.shape[1], 1)

    # ---- Rotate to DP frame.
    FpDP = np.cos(2*psi)*Fp + np.sin(2*psi)*Fc
    FcDP = -np.sin(2*psi)*Fp + np.cos(2*psi)*Fc

    # ---- Further rotate polarization by pi/4 if |Fp|<|Fc|.
    swapindex = (FpDP**2).sum(1) < (FcDP**2).sum(1)
    FpDP[swapindex, :] = FcDP[swapindex, :]
    FcDP[swapindex, :] = -FpDP[swapindex, :]
    psi[swapindex, :] = psi[swapindex, :] + np.pi/4

    return FpDP, FcDP, psi
