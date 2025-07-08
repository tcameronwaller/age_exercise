"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'merge_phenotype' is part of the 'transcriptomics' subpackage
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

# TCW; 2 July 2025
# I haven't implemented this yet. I decided instead to implement a more versatile
# method that can parallelize regressions over a large number of features such
# as genes or proteins.



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
import matplotlib

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.regression as preg
import partner.parallelization as prall
import age_exercise.phenotypes.organize_subject as aexph_sub
import age_exercise.phenotypes.organize_sample as exph_sample

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.


##########
# 2. Read source information from file.



###############################################################################
# Procedure


##########
# Manage main procedure.


def execute_procedure(
    path_directory_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

    raises:

    returns:

    """

    ##########
    # Parameters.
    project="age_exercise"
    routine="transcriptomics"
    procedure="regress_genes"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: regress_genes.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("routine: " + str(routine))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # 1. Initialize directories for read of source and write of product files.
    paths = aexph_sub.initialize_directories(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        initialize_routine=False,
        restore=True,
        report=report,
    )
    #paths["out_procedure_tables"]
    #paths["out_procedure_plot"]

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")


    ##########
    # 2. Read source information from file.
    if False:
        pail_source = read_source(
            paths=paths,
            report=report,
        )

    pass


###############################################################################
# End
