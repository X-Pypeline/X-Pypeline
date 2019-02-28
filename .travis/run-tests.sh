#!/bin/bash

coverage run ./setup.py test
coverage run --append `which xpipeline-workflow` --help
coverage run --append `which xpipeline-analysis` --help
