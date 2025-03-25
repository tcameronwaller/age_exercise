"""
Supply functionality for process and analysis of data from proteomics.

This module 'organize_spectroscopy' is part of the 'proteomics' package within
the 'age_exercise' package.

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
import math
import statistics
import pickle
import copy
import random
import itertools
import time
from datetime import datetime
#import dateutil # requires explicit installation

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
import partner.decomposition as pdecomp
import partner.plot as pplot
import partner.parallelization as prall
import age_exercise.proteomics.organize_subject as aexpr_sub

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.


def initialize_directories(
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
    #paths["out_procedure_lists"] = os.path.join(
    #    paths["out_procedure"], "lists",
    #)
    paths["out_procedure_tables"] = os.path.join(
        paths["out_procedure"], "tables",
    )
    paths["out_procedure_plot"] = os.path.join(
        paths["out_procedure"], "plot",
    )

    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"],
        paths["out_procedure_lists"],
        paths["out_procedure_tables"],
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
        print("module: age_exercise.proteomics.organize_spectroscopy.py")
        print("function: initialize_directories()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


##########
# 2. Read source information from file.


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

    # Define paths to parent directories.
    #paths["in_data"]
    #paths["in_parameters"]

    # Define paths to child files.
    path_file_table_batch_a = os.path.join(
        paths["in_data"], "study_age_exercise", "proteomics",
        "spectroscopy_ms2_global_zcr_2024-06-11_clean",
        "zr_061124_MS2_TMT_HK_TMTA_Proteins.tsv",
    )
    path_file_table_batch_b = os.path.join(
        paths["in_data"], "study_age_exercise", "proteomics",
        "spectroscopy_ms2_global_zcr_2024-06-11_clean",
        "zr_061424_MS2_TMT_HK_TMTB_Proteins.tsv",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of measurements.

    # Table of parameters for organization of the table of attributes for
    # subjects and samples.
    #types_columns = (
    #    aexpr_sub.define_type_columns_table_subject_feature_organization()
    #)
    pail["table_batch_a"] = pandas.read_csv(
        path_file_table_batch_a,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_batch_b"] = pandas.read_csv(
        path_file_table_batch_a=b,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.organize_spectroscopy.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of measurements in batch a: ")
        print(pail["table_batch_a"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Organize information.

# TCW; 25 March 2025
# manual clean-up of the output spectroscopy tables
# 1. filter to identity confidence "High"
# 2. delete columns:
# "Protein FDR Confidence: Combined"
# "Master"
# "Description"
# "Abundance Ratio: (Post, Sample) / (Pre, Sample)"
# everything after, except for: "Entrez Gene ID", "Ensembl Gene ID"
# 3. rename columns:
# "Accession": "identifier_protein_uniprot"
# "Gene Symbol": "identifier_gene_symbol"
# "Entrez Gene ID": "identifier_gene_entrez"
# "Ensembl Gene ID": "identifier_gene_ensembl"
# "": ""
# "": ""

# 4. organize sets of columns (samples)
# samples in batch "a"
# samples in batch "b"

# programmatic organization




def define_samples_by_batch():
    """
    Defines parameters.

    arguments:

    raises:

    returns:
        (dict): collection of information

    """

    # Collect information.
    pail = dict()

    # Samples in batch A.
    pail["names_samples_batch_a"] = [
        "BOH319_pre",
        "BOH319_post",
        "BUR210_pre",
        "BUR210_post",
        "COR775_pre",
        "COR775_post",
        "DIE776_pre",
        "DIE776_post",
        "HAN172_pre",
        "HAN172_post",
    ]
    # Samples in batch B.
    pail["names_samples_batch_b"] = [
        "ARM897_pre",
        "ARM897_post",
        "BEC65_pre",
        "BEC65_post",
        "BUC505_pre",
        "BUC505_post",
        "CLO122_pre",
        "CLO122_post",
        "HAL405_pre",
        "HAL405_post",
    ]

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
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

    raises:

    returns:

    """

    ##########
    # Parameters.
    project="age_exercise"
    routine="proteomics"
    procedure="organize_spectroscopy"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_olink.py")
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
    paths = initialize_directories(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=True,
        report=report,
    )

    ##########
    # 2. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=report,
    )
    print(pail_source["table_batch_a"])
    print(pail_source["table_batch_b"])

    ##########
    # 3.

    pail_batch_samples = define_samples_by_batch()
    #pail_batch_samples["names_samples_batch_a"]
    #pail_batch_samples["names_samples_batch_b"]


    pass


###############################################################################
# End
