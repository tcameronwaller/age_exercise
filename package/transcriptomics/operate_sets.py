"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'operate_sets' is part of the 'transcriptomics' package
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
    Copyright (C) 2024 Thomas Cameron Waller

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
import age_exercise.transcriptomics.organize_sample as extr_sample
import age_exercise.transcriptomics.organize_signal as extr_signal
import age_exercise.transcriptomics.select_gene_sets as extr_select
import age_exercise.proteomics.organize_subject as aexpr_sub

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.
# There is a hierarchy in these functions to initialize directories to manage
# the hierarchical tree structure of sub-procedures.


def initialize_directories_trunk(
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    restore=None,
    report=None,
):
    """
    Initialize directories for procedure's source and product files.

    arguments:
        project (str): name of project that normally corresponds to a single
            Python package
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics' that normally corresponds to a single Python package
            or subpackage
        procedure (str): name of procedure, a step in the routine process that
            normally corresponds to a single Python module within the package
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        restore (bool): whether to remove previous versions of data
        report (bool): whether to print reports

    raises:

    returns:
        (dict<str>): collection of paths to directories

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    # Broad.
    paths["dock"] = path_directory_dock
    paths["in_data"] = os.path.join(
        paths["dock"], "in_data",
    )
    paths["in_demonstration"] = os.path.join(
        paths["dock"], "in_demonstration",
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters",
    )
    paths["in_parameters_private"] = os.path.join(
        paths["dock"], "in_parameters_private",
    )
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    paths["out_routine"] = os.path.join(
        paths["out_project"], str(routine),
    )
    paths["out_procedure"] = os.path.join(
        paths["out_routine"], str(procedure),
    )
    # Specific.
    paths["in_sets_gene"] = os.path.join(
        paths["in_parameters_private"], project, routine,
        "sets_gene",
    )
    paths["out_procedure_data"] = os.path.join(
        paths["out_procedure"], "data",
    )
    paths["out_procedure_plot"] = os.path.join(
        paths["out_procedure"], "plot",
    )

    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"],
        paths["out_procedure_data"],
        paths["out_procedure_plot"],
    ]
    # Remove previous directories and files to avoid version or batch
    # confusion.
    if restore:
        for path in paths_initialization:
            putly.remove_directory(path=path) # caution
            pass
    # Create directories.
    for path in paths_initialization:
        putly.create_directories(
            path=path,
        )
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.compare_sets_groups.py")
        print("function: initialize_directories_trunk()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


##########
# 2. Read source information from file.



##########
# 3. Organize information from source.



###############################################################################
# Procedure


##########
# Execute main procedure at highest level of hierarchy.


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
    procedure="operate_sets"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.operate_sets.py")
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
    paths = initialize_directories_trunk(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=True,
        report=report,
    )

    ##########
    # 2.1. Read and count unique genes in sets.
    path_directory_sets_gene = os.path.join(
        paths["in_sets_gene"], "sets_gene_2025-01-09",
    )
    table_counts_sets_gene = (
        putly.read_child_files_text_list_count_unique_items(
            path_directory=path_directory_sets_gene,
            name_file_prefix="",
            name_file_suffix=".txt",
            name_file_not="",
            report=report,
    ))

    ##########
    # Subjective combination of sets of genes identified by threshold.

    # TODO: TCW; 9 January 2025
    # Organize operation below within a cleaner, more versatile function
    # arguments:
    # path_directory_read
    # names_sets <-- corresponding to names of files within directory
    # name_file_product
    # path_directory_write

    # Define paths to directories.
    path_directory_comparison_sets = os.path.join(
        paths["out_routine"], "compare_sets_groups", "data",
        "sets_subjective_match_directionality_younger_elder_omega3",
    )
    sets_match_younger_elder_omega3 = dict()
    sets_match_younger_elder_omega3["main"] = [
        "3_adipose_omega3_permissive",
        "4_adipose_omega3",
        "5_adipose_age_placebo_omega3",
    ]
    set_match_younger_elder_omega3 = (
        aexpr_sub.read_extract_combine_custom_feature_sets(
            names_sets=sets_match_younger_elder_omega3,
            features_available=None,
            path_directory=path_directory_comparison_sets,
            report=True,
    ))
    # Copy information.
    genes_match_younger_elder_omega3 = copy.deepcopy(
        set_match_younger_elder_omega3["main"]
    )
    # Collect unique names of genes in set.
    genes_match_younger_elder_omega3 = putly.collect_unique_elements(
        elements=genes_match_younger_elder_omega3,
    )
    # Collect information.
    # Collections of files.
    name_set = str(
        "genes_subjective_match_directionality_younger_elder_omega3"
    )
    pail_write_lists = dict()
    pail_write_lists[name_set] = (
        genes_match_younger_elder_omega3
    )
    # Define paths to directories.
    path_directory_data_sets = os.path.join(
        paths["out_procedure_data"], "sets",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_data_sets,
    )
    # Write product information to file.
    putly.write_lists_to_file_text(
        pail_write=pail_write_lists,
        path_directory=path_directory_data_sets,
        delimiter="\n",
    )


    ##########
    # Set operations and combinations.


    # Collections of files.
    pail_write_lists = dict()
    pail_write_lists["genes_adipose_1_up_25_up"] = (
        genes_1_up_25_up
    )
    pail_write_lists["genes_adipose_1_down_25_down"] = (
        genes_1_down_25_down
    )
    pail_write_lists["genes_adipose_same_1_25"] = (
        genes_same_1_25
    )
    pail_write_lists["genes_adipose_1_up_25_down"] = (
        genes_1_up_25_down
    )
    pail_write_lists["genes_adipose_1_down_25_down"] = (
        genes_up_1_up_25
    )
    pail_write_lists["genes_adipose_opposite_1_25"] = (
        genes_same_1_25
    )

    pass


###############################################################################
# End
