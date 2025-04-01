"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'compare_groups' is part of the 'phenotypes' subpackage
within the 'age_exercise' package.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This file is part of the project package directory 'age_exercise'
    (https://github.com/tcameronwaller/age_exercise/).

    Project 'age_exercise' supports data analysis with team in endocrinology.
    Copyright (C) 2025 Thomas Cameron Waller

    The code within project 'age_exercise' is free software: you can
    redistribute it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either version 3 of
    the GNU General Public License, or (at your option) any later version.

    The code within project 'age_exercise' is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'age_exercise'. If not, see <http://www.gnu.org/licenses/>.
"""


###############################################################################
# Notes


###############################################################################
# Installation and importation

# Standard
import sys
#print(sys.path)
import os
import os.path
import math
import statistics
import pickle
import copy
import random
import itertools
import time

# Relevant
import numpy
import scipy.stats
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"
pandas.set_option('future.no_silent_downcasting', True) # set option to suppress warnings
from warnings import simplefilter
simplefilter(action="ignore", category=pandas.errors.PerformanceWarning)

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
#import partner.regression as preg
import partner.plot as pplot
import partner.parallelization as prall
import age_exercise.phenotypes.organize_subject as aexph_sub
import age_exercise.phenotypes.organize_sample as extr_sample
import age_exercise.transcriptomics.organize_signal as extr_signal
import age_exercise.transcriptomics.select_gene_sets as extr_select

###############################################################################
# Functionality

# TODO; TCW; 31 March 2025
# Move to this module the portion of "organize_subject.py" about comparing
# quantitative features between groups of observations.
# Adapt functionality from module "transcriptomics.compare_sets_groups".

###############################################################################
# Procedure



###############################################################################
# End
