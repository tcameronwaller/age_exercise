"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'merge_phenotype' is part of the 'phenotypes' subpackage
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

# TODO: TCW; 2 May 2025




def read_organize_source_subject(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of source information read from file

    """

    # Define paths to parent directories.
    #paths["in_data"]
    #paths["in_demonstration"]
    #paths["in_parameters"]
    #paths["in_parameters_private"]

    # Define paths to child files.
    path_file_table_subject = os.path.join(
        paths["out_project"], "phenotypes", "organize_subject", "tables",
        "table_subject.tsv",
    )
    path_file_table_sample = os.path.join(
        paths["out_project"], "phenotypes", "organize_sample", "tables",
        "table_sample.tsv",
    )

    # Read information from file.

    # Table of properties and attributes for subjects and samples.
    types_columns = define_type_table_columns_subject_sample()
    table = pandas.read_csv(
        path_file_table_subject,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Table of parameters for organization of features.
    pail_feature = (
        aexph_sub.read_organize_source_table_subject_feature_organization(
            paths=paths,
            report=False,
    ))
    #pail_feature["types_columns"]
    #pail_feature["translations_feature_forward"]
    #pail_feature["translations_feature_reverse"]
    #pail_feature["columns_all"]
    #pail_feature["columns_quantitative"]
    #pail_feature["columns_olink_plasma"]
    #pail_feature["columns_olink_muscle"]
    #pail_feature["columns_olink_adipose"]

    # Organize information.
    # Determine names of columns for unique, relevant features.
    # The names from the keys of dictionary "types_columns" appear first in the
    # list, giving these priority in subsequent sorts of columns in the table.
    features_relevant = list(types_columns.keys())
    pail_quantitative = (
        define_type_table_columns_subject_sample_quantitative_continuous()
    )
    features_quantitative = pail_quantitative["names_features"]
    #features_relevant.extend(pail_feature["columns_quantitative"])
    #features_relevant.extend(pail_feature["columns_olink_plasma"])
    #features_relevant.extend(pail_feature["columns_olink_muscle"])
    #features_relevant.extend(pail_feature["columns_olink_adipose"])
    features_relevant = putly.collect_unique_elements(
        elements=features_relevant,
    )
    features_quantitative = putly.collect_unique_elements(
        elements=features_quantitative,
    )

    # Collect information.
    pail = dict()
    pail["table_subject_sample"] = table.copy(deep=True)
    pail["features_relevant"] = copy.deepcopy(features_relevant)
    pail["features_quantitative"] = copy.deepcopy(features_quantitative)
    pail["translations_feature_reverse"] = (
        pail_feature["translations_feature_reverse"]
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
        print("function: read_organize_source_subject_sample()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of source information read from file

    """

    # Read source information for subjects.
    pail_subject = read_organize_source_subject(
        paths=paths,
        report=False,
    )

    # Read source information for samples.
    pail_sample = read_organize_source_sample(
        paths=paths,
        report=False,
    )

    # Read source information for signals.
    pail_signal = read_organize_source_signal(
        paths=paths,
        report=False,
    )

    # Collect information.
    pail = dict()
    pail["table_subject"] = pail_subject["table"]
    pail["table_sample"] = pail_sample["table"]
    pail["table_signal"] = pail_signal["table"]

    #pail["features_relevant"] = pail_sample["features_relevant"]
    #pail["features_quantitative"] = pail_sample["features_quantitative"]
    #pail["instances_parameter"] = pail_parameter["records"]
    #pail["translations_feature_reverse"] = (
    #    pail_sample["translations_feature_reverse"]
    #)

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: merge_phenotype.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        #print("table of measurements in batch a: ")
        #print(pail["table_batch_a"].iloc[0:10, 0:])
        #putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail



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
    routine="phenotypes"
    procedure="compare_groups"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
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

    ##########
    # 2. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=report,
    )
    #pail_source["table_sample"]
    #pail_source["features_relevant"]
    #pail_source["instances_parameter"]

    ##########
    # 3. Merge phenotypes for subjects or samples to signals of genes.




    # ...



    if False:
        ##########
        # 5. Write information to file.
        # Collect information.
        # Collections of files.
        #pail_write_lists = dict()
        pail_write_tables = dict()
        #pail_write_tables[str("table_merge")] = pail_parts["table_group"]
        pail_write_objects = dict()
        #pail_write_objects[str("samples")]
        # Write product information to file.
        putly.write_tables_to_file(
            pail_write=pail_write_tables,
            path_directory=paths["out_procedure_tables"],
            reset_index_rows=False,
            write_index_rows=False,
            write_index_columns=True,
            type="text",
            delimiter="\t",
            suffix=".tsv",
        )
    pass


###############################################################################
# End
