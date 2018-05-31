#!/bin/bash

coverage run ./setup.py test
coverage run --append `which setUpJobs` --help
