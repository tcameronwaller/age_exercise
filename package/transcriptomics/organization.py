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
    paths["in_parameters_private"] = os.path.join(
        paths["dock"], "in_parameters_private", str(project), str(technology),
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


##########
# 1. Read information from file.


def define_table_column_types_quantity():
    """
    Defines the variable types of columns within table for values of intensity.

    Review: TCW; 24 June 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types in columns within table.
    types_columns = dict()
    types_columns["identifier_gene"] = "string"
    types_columns["gene_id"] = "string"
    types_columns["gene_name"] = "string"
    types_columns["exon_number"] = "string"
    types_columns["gene_type"] = "string"
    types_columns["chromosome"] = "string"
    # Return information.
    return types_columns


def define_table_column_types_sample():
    """
    Defines the variable types of columns within table for attributes of
    samples.

    Review: TCW; 24 June 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion"] = "int32"
    types_columns["identifier"] = "string"
    types_columns["path_file"] = "string"
    types_columns["sample_plate"] = "string"
    types_columns["plate"] = "string"
    types_columns["sample"] = "string"
    types_columns["sample_plate"] = "string"
    types_columns["condition_code"] = "string"
    types_columns["tissue"] = "string"
    types_columns["condition"] = "string"
    types_columns["subject"] = "string"
    types_columns["note_condition"] = "string"
    #types_columns["sex"] = "string"
    #types_columns["age"] = "int32"
    #types_columns["body_mass"] = "float32"
    #types_columns["body_mass_index"] = "float32"
    # Return information.
    return types_columns


def read_source(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of source information read from file

    """

    # Define paths to parent directories.
    #paths["in_data"]
    #paths["in_parameters"]

    # Define paths to child files.
    path_file_table_sample = os.path.join(
        paths["in_parameters_private"], "quantification_2024-07-14",
        "attributes_samples", "table_samples_attributes.tsv",
    )
    path_file_table_quantity_adipose = os.path.join(
        paths["in_data"], "quantification_2024-07-14",
        "organization", "quantification_rna_reads_gene_adipose.tsv",
    )
    path_file_table_quantity_muscle = os.path.join(
        paths["in_data"], "quantification_2024-07-14",
        "organization", "quantification_rna_reads_gene_muscle.tsv",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of samples and their attributes.
    types_columns = define_table_column_types_sample()
    pail["table_sample"] = pandas.read_csv(
        path_file_table_sample,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )

    # Table of values of intensity across samples and proteins.
    types_columns = define_table_column_types_quantity()
    pail["table_quantity_adipose"] = pandas.read_csv(
        path_file_table_quantity_adipose,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    pail["table_quantity_muscle"] = pandas.read_csv(
        path_file_table_quantity_muscle,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("module: exercise.transcriptomics.organization.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=4)
        print("sample table: ")
        print(pail["table_sample"])
        putly.print_terminal_partition(level=4)
        count_rows = (pail["table_quantity_adipose"].shape[0])
        count_columns = (pail["table_quantity_adipose"].shape[1])
        print("quantity adipose table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_quantity_adipose"])
        putly.print_terminal_partition(level=4)
        count_rows = (pail["table_quantity_muscle"].shape[0])
        count_columns = (pail["table_quantity_muscle"].shape[1])
        print("quantity muscle table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_quantity_muscle"])
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail


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

    ##########
    # Report.
    print("system: local")
    print("project: exercise")
    print("technology: transcriptomics")
    print("procedure: 2_organization")
    print("set: organization")

    ##########
    # Initialize directories.
    paths = initialize_directories(
        project="exercise",
        technology="transcriptomics",
        set="organization",
        path_directory_dock=path_directory_dock,
        restore=True,
    )

    ##########
    # 1. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=True,
    )




    pass


###############################################################################
# End
