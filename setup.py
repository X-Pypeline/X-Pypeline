#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2018)
#
# This file is part of the xpipeline python package.
#
# xpipeline is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xpipeline is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with xpipeline.  If not, see <http://www.gnu.org/licenses/>.

"""Setup the xpipeline package
"""

from __future__ import print_function

import sys
if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")

import glob
import os.path

from setuptools import (setup, find_packages)
from distutils.extension import Extension

try:
    from Cython.Build import cythonize
except ImportError:
    raise ImportError("Cython needed for cpp extensions")

try:
   import numpy
except ImportError:
    raise ImportError("Building Cython extensions requires numpy.")

# set basic metadata
PACKAGENAME = 'xpipeline'
DISTNAME = 'xpipeline'
AUTHOR = 'Scott Coughlin'
AUTHOR_EMAIL = 'scott.coughlin@ligo.org'
LICENSE = 'GPLv3'

cmdclass = {}

# -- versioning ---------------------------------------------------------------

import versioneer
__version__ = versioneer.get_version()
cmdclass.update(versioneer.get_cmdclass())

# -- documentation ------------------------------------------------------------

# import sphinx commands
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    cmdclass['build_sphinx'] = BuildDoc

# -- dependencies -------------------------------------------------------------

setup_requires = [
    'setuptools',
    'pytest-runner',
]

install_requires = [
    'lalsuite>6.52',
    'cython>=0.29.5',
    'filelock>=3.0.10',
    'gwpy>=0.12',
    'lscsoft-glue>=2.0.0',
    'pandas >= 0.22 ; python_version >= \'3.5\'',
    'pandas < 0.21 ; python_version == \'3.4\'',
    'pandas >= 0.22 ; python_version == \'2.7\'',
]

tests_require = [
    'pytest'
]

extras_require = {
    'doc': [
        'ipython',
        'sphinx',
        'numpydoc',
        'sphinx_rtd_theme',
        'sphinxcontrib_programoutput',
    ],
}

extensions = [
    Extension("xpipeline.cluster.nearestneighbor",
              ["xpipeline/cluster/src/nearestneighbors.pyx"],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-std=c++11'],
              language='c++'
              ),

    Extension('xpipeline.cluster.clusterproperties',
              ['xpipeline/cluster/src/clusterproperties.pyx'],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-std=c++11'],
              language='c++'
              ),
    Extension('xpipeline.cluster.clustersum',
              ['xpipeline/cluster/src/clustersum.pyx'],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-std=c++11'],
              language='c++'
              ),
]

# -- run setup ----------------------------------------------------------------

packagenames = find_packages()
scripts = glob.glob(os.path.join('bin', '*'))

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version=__version__,
      description=None,
      long_description=None,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=packagenames,
      ext_modules=cythonize(extensions),
      include_package_data=True,
      cmdclass=cmdclass,
      scripts=scripts,
      setup_requires=setup_requires,
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require=extras_require,
      test_suite='xpipeline.tests',
      use_2to3=True,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
)
