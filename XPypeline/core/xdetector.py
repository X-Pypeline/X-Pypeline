# Copyright (C) 2012  Alex Nitz
#
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""This module provides utilities for calculating detector responses.
This was grabbed from 
https://github.com/gwastro/pycbc/blob/master/pycbc/detector.py
"""
import lalsimulation
import numpy as np
import lal
from numpy import cos, sin
from gwpy.timeseries import TimeSeries


class Detector(object):
    """A gravitaional wave detector
    """
    def __init__(self, detector_name):
        self.name = str(detector_name)
        self.frDetector =  lalsimulation.DetectorPrefixToLALDetector(self.name)
        self.response = self.frDetector.response
        self.location = self.frDetector.location
        self.latitude = self.frDetector.frDetector.vertexLatitudeRadians
        self.longitude = self.frDetector.frDetector.vertexLongitudeRadians


    def light_travel_time_to_detector(self, det):
        """ Return the light travel time from this detector
        Parameters
        ----------
        detector: Detector
            The other detector to determine the light travel time to.
        Returns
        -------
        time: float
            The light travel time in seconds
        """
        return lal.LightTravelTime(self.frDetector, det.frDetector) * 1e-9


    def compute_antenna_response(self, phi, theta, psi=0):
        """Return the detector response.
        """
        d = self.response.flatten()
        if not psi:
            psi = np.zeros(len(phi))
        m1 = sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi)
        m2 = -cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi)
        m3 = sin(theta)*sin(psi)
        n1 = -sin(phi)*sin(psi)-cos(phi)*cos(theta)*cos(psi)
        n2 = cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi)
        n3 = sin(theta)*cos(psi)
        k1 = m2*n3 - m3*n2
        k2 = m3*n1 - m1*n3
        k3 = m1*n2 - m2*n1
        mm = np.array([m1*m1, m1*m2, m1*m3, m2*m1,
                      m2*m2, m2*m3, m3*m1, m3*m2, m3*m3])
        mn = np.array([m1*n1, m1*n2, m1*n3, m2*n1,
                      m2*n2, m2*n3, m3*n1, m3*n2, m3*n3])
        nm = np.array([n1*m1, n1*m2, n1*m3, n2*m1,
                      n2*m2, n2*m3, n3*m1, n3*m2, n3*m3])
        nn = np.array([n1*n1, n1*n2, n1*n3, n2*n1,
                      n2*n2, n2*n3, n3*n1, n3*n2, n3*n3])
        kk = np.array([k1*k1, k1*k2, k1*k3, k2*k1,
                      k2*k2, k2*k3, k3*k1, k3*k2, k3*k3])
        mk = np.array([m1*k1, m1*k2, m1*k3, m2*k1,
                      m2*k2, m2*k3, m3*k1, m3*k2, m3*k3])
        km = np.array([k1*m1, k1*m2, k1*m3, k2*m1,
                      k2*m2, k2*m3, k3*m1, k3*m2, k3*m3])
        nk = np.array([n1*k1, n1*k2, n1*k3, n2*k1,
                      n2*k2, n2*k3, n3*k1, n3*k2, n3*k3])
        kn = np.array([k1*n1, k1*n2, k1*n3, k2*n1,
                      k2*n2, k2*n3, k3*n1, k3*n2, k3*n3])
        e_plus = mm - nn
        e_cross = mn + nm
        e_breathing = mm + nn
        e_long = kk
        e_v1 = mk + km
        e_v2 = nk + kn

        # ----- Compute waveform projected onto antenna pattern.
        Fp = e_plus*d
        Fc = e_cross*d
        Fb = e_breathing*d
        FL = e_long*d
        F1 = e_v1*d
        F2 = e_v2*d
        return Fp, Fc, Fb, FL, F1, F2


    def time_delay_from_earth_center_phi_theta(self, phi, theta):
        """Return the time delay from the earth center
        """
        c =lal.C_SI
        # Construct the sky direction assumed from skyPosition
        omega = np.array(sin(theta)*cos(phi), # x
                 sin(theta)*sin(phi), # y
                 cos(theta))          # z
         
        # Calculate the time delay for each detector (in seconds)
        deltaT = -omega*self.location / c;

        return 


    def time_delay_from_earth_center(self, right_ascension, declination, t_gps):
        """Return the time delay from the earth center
        """
        return lal.TimeDelayFromEarthCenter(self.location,
                      float(right_ascension), float(declination), float(t_gps))


    def time_delay_from_detector(self, other_detector, right_ascension,
                                 declination, t_gps):
        """Return the time delay from the given to detector for a signal with
        the given sky location; i.e. return `t1 - t2` where `t1` is the
        arrival time in this detector and `t2` is the arrival time in the
        other detector. Note that this would return the same value as
        `time_delay_from_earth_center` if `other_detector` was geocentric.
        Parameters
        ----------
        other_detector : detector.Detector
            A detector instance.
        right_ascension : float
            The right ascension (in rad) of the signal.
        declination : float
            The declination (in rad) of the signal.
        t_gps : float
            The GPS time (in s) of the signal.
        Returns
        -------
        float
            The arrival time difference between the detectors.
        """
        return lal.ArrivalTimeDiff(self.location, other_detector.location,
                                   float(right_ascension), float(declination),
                                   float(t_gps))


    def project_wave(self, hp, hc, longitude, latitude, polarization):
        """Return the strain of a wave with given amplitudes and angles as
        measured by the detector.
        """
        h_lal = lalsimulation.SimDetectorStrainREAL8TimeSeries(
                hp.astype(np.float64).lal(), hc.astype(np.float64).lal(),
                longitude, latitude, polarization, self.frDetector)
        return TimeSeries(
                h_lal.data.data, delta_t=h_lal.deltaT, epoch=h_lal.epoch,
                dtype=np.float64)


    def optimal_orientation(self, t_gps):
        """Return the optimal orientation in right ascension and declination
           for a given GPS time.
        """
        ra = self.longitude + (lal.GreenwichMeanSiderealTime(t_gps) % (2*np.pi))
        dec = self.latitude
        return ra, dec


def overhead_antenna_pattern(right_ascension, declination, polarization):
    """Return the detector response where (0, 0) indicates an overhead source. 
    This functions uses coordinates such that the detector can be thought to
    be on the north pole.
    Parameters
    ----------
    right_ascention: float
    declination: float
    polarization: float
    Returns
    -------
    f_plus: float
    f_cros: float   
    """
    # convert from declination coordinate to polar (angle dropped from north axis)
    theta = np.pi / 2.0 - declination

    f_plus  = - (1.0/2.0) * (1.0 + cos(theta)*cos(theta)) * \
                cos (2.0 * right_ascension) * cos (2.0 * polarization) - \
                cos(theta) * sin(2.0*right_ascension) * sin (2.0 * polarization)

    f_cross =   (1.0/2.0) * (1.0 + cos(theta)*cos(theta)) * \
                cos (2.0 * right_ascension) * sin (2.0* polarization) - \
                cos(theta) * sin(2.0*right_ascension) * cos (2.0 * polarization)

    return f_plus, f_cross


def effective_distance(distance, inclination, f_plus, f_cross):
    return distance / np.sqrt( ( 1 + np.cos( inclination )**2 )**2 / 4 * f_plus**2 + np.cos( inclination )**2 * f_cross**2 )
