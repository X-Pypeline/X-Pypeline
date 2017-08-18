# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of XPypeline.
#
# GWpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy.  If not, see <http://www.gnu.org/licenses/>.

"""This module converts ra to dec and saves to a file if supplied
"""

import lal
import numpy as np

def radectoearth(ra, dec, gps, filename=''):
    gmst = lal.GreenwichMeanSiderealTime(gps)
    print(gmst)
    gmst_deg = gmst / 86400 * 360;
    # ---- Compute Earth-based coordinates of targeted sky location, in degrees.
    phi_deg = ra - gmst_deg;
    theta_deg  = 90 - dec;
    # ---- Convert to radians.
    phi = phi_deg * 2 * np.pi / 360;
    phi = np.mod(phi+np.pi,2*np.pi)-np.pi;
    theta = theta_deg * 2 * np.pi / 360;
    if filename:
        with open(filename, 'w') as f:
            f.write('{0} {1}'.format(phi, theta))
            f.close()
        return
    else:
        return phi, theta
