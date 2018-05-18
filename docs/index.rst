.. xpipeline documentation master file, created by
   sphinx-quickstart on Thu Apr 21 14:05:08 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to xpipeline's documentation!
=====================================

`X Pipeline <https://trac.ligo.caltech.edu/xpipeline/wiki>`_ is a burst gravitational-wave search algorithm. It can perform SN, GRB, and ALLSKY searches. For full details see `Sutton et al. (IOP, 2010) <http://iopscience.iop.org/article/10.1088/1367-2630/12/5/053034/meta>`_ and `Was et al. (PRD, 2012) <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.86.022003>`_ 

Installing xpipeline
--------------------

The easiest method to install xpipeline is using `pip <https://pip.pypa.io/en/stable/>`_ directly from the `GitHub repository <https://github.com/xpipeline/xpipeline.git>`_:

.. code-block:: bash

   $ pip install git+https://github.com/xpipeline/xpipeline.git

How to run xpipeline
--------------------

The main product of this package is the command-line executable `setUpJobs`, which runs a script that sets up X-Pipeline gravitational-wave searches of SN GRB and ALLSKY/SPHRAD. 

To run an analysis:

.. code-block:: bash

   $ setUpJobs -p grb_full.ini --grb-name GRB160830 -g 1156610000 -r 307.65 -d 45.72 -i H1 L1 V1 -s grb --FAP=99.9 

where ``-g`` is the GPS time of the GRB event, and ``-p grb_full.ini`` is the path of your configuration file.

For a full list of command-line argument and options, run

.. code-block:: bash

   $ setUpJobs --help

For more details see :ref:`command-line`.

Package documentation
---------------------

Please consult these pages for more details on using xpipeline:

.. toctree::
   :maxdepth: 1

   command-line/index
   examples/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
