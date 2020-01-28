# X-Pypeline
An attempt to create a pythonic version of X-Pipeline: https://trac.ligo.caltech.edu/xpipeline/

# Installation

## CONDA
### Basic
`conda create -c conda-forge -n xpipeline-py37 gwpy cython pytables pandas python=3.7`
`pip install git+https://github.com/X-Pypeline/X-Pypeline.git`


# Reproduce Thesis or Paper Analysis on CIT
After cloning the repository locally to CIT copy the `./examples/snews` to somewhere you want to do the analysis on CIT
```
cp -r examples/snews/ .
cd snews
bash xpypeline-setup-on-cit.sh
condor_submit_dag grb_alljobs.dag
```


# Badges

[![Build Status](https://travis-ci.org/X-Pypeline/X-Pypeline.svg?branch=develop)](https://travis-ci.org/X-Pypeline/X-Pypeline)
[![Coverage Status](https://coveralls.io/repos/github/X-Pypeline/X-Pypeline/badge.svg?branch=develop)](https://coveralls.io/github/X-Pypeline/X-Pypeline?branch=develop)
