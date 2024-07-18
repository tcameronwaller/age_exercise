"""
Supply functionality for process and analysis of data from proteomics using
mass spectroscopy.

This module 'organization' is part of the 'proteomics' package within the
'exercise' package.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This file is part of the project package directory 'exercise'
    (https://github.com/tcameronwaller/exercise/).

    Project 'exercise' supports data analysis with team in endocrinology.
    Copyright (C) 2024 Thomas Cameron Waller

    The code within project 'exercise' is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the GNU
    General Public License, or (at your option) any later version.

    The code within project 'exercise' is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'exercise'. If not, see <http://www.gnu.org/licenses/>.
"""


###############################################################################
# Notes


###############################################################################
# Installation and importation

# Standard
import sys
#print(sys.path)
import os
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

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
#import partner.regression as preg
import partner.plot as pplot
import partner.parallelization as prall

###############################################################################
# Functionality


##########
# Initialization


def initialize_directories(
    project=None,
    technology=None,
    set=None,
    path_directory_dock=None,
    restore=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
        technology (str): name of technology, either 'transcriptomics' or
            'proteomics'
        set (str): name of set or step in process procedure
        path_directory_dock (str): path to dock directory for source and
            product directories and files
        restore (bool): whether to remove previous versions of data

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_directory_dock
    paths["in_data"] = os.path.join(
        paths["dock"], "in_data", str(project), str(technology),
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters", str(project), str(technology),
    )
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    paths["out_technology"] = os.path.join(
        paths["out_project"], str(technology),
    )
    paths["out_set"] = os.path.join(
        paths["out_technology"], str(set),
    )
    paths["out_test"] = os.path.join(
        paths["out_set"], "test",
    )
    paths["out_table"] = os.path.join(
        paths["out_set"], "table",
    )
    paths["out_plot"] = os.path.join(
        paths["out_set"], "plot",
    )
    paths_initialization = [
        paths["out_project"],
        paths["out_technology"],
        paths["out_set"],
        paths["out_test"],
        paths["out_table"],
        paths["out_plot"],
    ]
    # Remove previous files to avoid version or batch confusion.
    if restore:
        for path in paths_initialization:
            putly.remove_directory(path=path)
    # Initialize directories.
    for path in paths_initialization:
        putly.create_directories(
            path=path,
        )
    # Return information.
    return paths



###############################################################################
# Procedure


def execute_procedure(
    path_directory_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_dock (str): path to dock directory for source and
            product directories and files

    raises:

    returns:

    """

    # Report.
    print("system: local")
    print("project: exercise")
    print("technology: transcriptomics")
    print("set: organization")

    # Initialize directories.
    paths = initialize_directories(
        project="exercise",
        technology="transcriptomics",
        set="organization",
        path_directory_dock=path_directory_dock,
        restore=True,
    )


    pass


###############################################################################
# End
