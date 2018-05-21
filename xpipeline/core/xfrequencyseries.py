#!/usr/bin/env python

# ---- Import standard modules to the python path.
import numpy as np
from collections import OrderedDict

class XFrequencySeriesDict(OrderedDict):
    def to_dominant_polarization_frame(self, phi, theta):
        """Shift timeseries by assosciated time delay

            Parameters
            ----------
            antenna_patterns : `dict`
                key-wise pair of antenna pattern name and value
                {'f_plus' : 0.045}
                available antenna patterns include
                'f_plus', 'f_cross', 'f_scalar',
                'f_long', 'f_one', 'f_two'

        """
        wFp = np.zeros([num_of_frequency_bins, num_of_detectors])
        wFc = np.zeros([num_of_frequency_bins, num_of_detectors])
        wFb = np.zeros([num_of_frequency_bins, num_of_detectors])
        wFL = np.zeros([num_of_frequency_bins, num_of_detectors])
        wF1 = np.zeros([num_of_frequency_bins, num_of_detectors])
        wF2 = np.zeros([num_of_frequency_bins, num_of_detectors])
        for idetidx, iasd in enumerate(asds.values()):
            detector = Detector(list(asds.keys())[idetidx])
            [Fp, Fc, Fb, FL, F1, F2] = detector.compute_antenna_response([phi], [theta])
            wFp[:, idetidx] = Fp / iasd
            wFc[:, idetidx] = Fc / iasd
            wFb[:, idetidx] = Fb / iasd
            wFL[:, idetidx] = FL / iasd
            wF1[:, idetidx] = F1 / iasd
            wF2[:, idetidx] = F2 / iasd

        wFp, wFc, psi = convert_to_dominant_polarization_frame(wFp, wFc)
        self['wFp'] = wFp
        self['wFc'] = wFc
        self['wFb'] = wFb
        self['wFL'] = wFL
        self['wF1'] = wF1
        self['wF2'] = wF2


def convert_to_dominant_polarization_frame(Fp, Fc):
    """Take in stream of fplus and f_cross and convert to DPF

       DPF is the Cominant Polarization Frame

        Parameters
        ----------
            Fp : `float`
            Fc :  `float`
    """
    # ---- Compute rotation needed to reach DP frame.
    psi = np.zeros([len(Fp), 1])
    psi[:, 0] = 1/4*np.arctan(2*(np.sum(Fp*Fc, 1))/(np.sum(Fp*Fp, 1)-np.sum(Fc*Fc, 1)))
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
