#!/usr/bin/env python

# ---- Import standard modules to the python path.

from __future__ import division

import sys
import os
import random
import string
import shutil
import ConfigParser
import optparse
import json
import rlcompleter
import pdb

import numpy as np

from scipy import signal
from scipy.interpolate import InterpolatedUnivariateSpline

from matplotlib import use
use('agg')
from matplotlib import (pyplot as plt, cm)
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from glue import datafind

from gwpy.timeseries import TimeSeries

pdb.Pdb.complete = rlcompleter.Completer(locals()).complete



if __name__ == '__main__':
    main()
