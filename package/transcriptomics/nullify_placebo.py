"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

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


# adipose_14 nullify from adipose_15

# adipose_27 nullify from adipose_28

# for each gene in placebo list...
# find the gene's row in the rank table
# set the rank factor to 0



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


##########
# 1. Initialize directories for read of source and write of product files.


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
    path_directory_sets_gene = os.path.join(
        paths["out_project"], "transcriptomics", "select_gene_sets", "data",
        "sets_gene",
    )
    path_directory_tables = os.path.join(
        paths["out_project"], "transcriptomics", "select_gene_sets", "data",
        "pickle",
    )

    # Define paths to child files.
    path_file_table_adipose_15 = os.path.join(
        path_directory_tables,
        "table_adipose_15_age-elder-intervention-omega3_visit.pickle",
    )
    path_file_table_adipose_28 = os.path.join(
        path_directory_tables,
        "table_adipose_28_age-elder-intervention-any_placebo-omega3.pickle",
    )
    path_file_table_adipose_32 = os.path.join(
        path_directory_tables,
        str(
            "table_adipose_32_age-elder-intervention-any_placebo-omega3_no-" +
            "subject.pickle"
        ),
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Sets of genes.
    pail["genes_adipose_14"] = aexph_sub.read_extract_set_features(
        name_set="adipose_14_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    pail["genes_adipose_27"] = aexph_sub.read_extract_set_features(
        name_set="adipose_27_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    pail["genes_adipose_31"] = aexph_sub.read_extract_set_features(
        name_set="adipose_31_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )


    # Table of properties for subjects.
    pail["table_adipose_15"] = pandas.read_pickle(
        path_file_table_adipose_15,
    )
    pail["table_adipose_28"] = pandas.read_pickle(
        path_file_table_adipose_28,
    )
    pail["table_adipose_32"] = pandas.read_pickle(
        path_file_table_adipose_32,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.merge_phenotype.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table: ")
        print(pail["table_adipose_28"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    procedure="nullify_placebo"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.nullify_placebo.py")
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

    ##########
    # 2. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=report,
    )

    # define null value
    value_null = float(0.0)

    ##########
    # adipose_15
    table = pail_source["table_adipose_15"].copy(deep=True)
    genes_placebo = copy.deepcopy(pail_source["genes_adipose_14"])
    table["rank_fold_p_placebo_null"] = table.apply(
        lambda row:
            (value_null)
            if (
                row["gene_identifier_base"] in genes_placebo
            ) else (row["rank_fold_p"]),
        axis="columns", # apply function to each row
    )
    # Sort rows within table.
    table.sort_values(
        by=[
            "rank_fold_p_placebo_null",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Write information to file.
    # Collect information.
    # Collections of files.
    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_adipose_15")] = table
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

    ##########
    # adipose_28
    table = pail_source["table_adipose_28"].copy(deep=True)
    genes_placebo = copy.deepcopy(pail_source["genes_adipose_27"])
    table["rank_fold_p_placebo_null"] = table.apply(
        lambda row:
            (value_null)
            if (
                row["gene_identifier_base"] in genes_placebo
            ) else (row["rank_fold_p"]),
        axis="columns", # apply function to each row
    )
    # Sort rows within table.
    table.sort_values(
        by=[
            "rank_fold_p_placebo_null",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Write information to file.
    # Collect information.
    # Collections of files.
    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_adipose_28")] = table
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

    ##########
    # adipose_32
    table = pail_source["table_adipose_32"].copy(deep=True)
    genes_placebo = copy.deepcopy(pail_source["genes_adipose_31"])
    table["rank_fold_p_placebo_null"] = table.apply(
        lambda row:
            (value_null)
            if (
                row["gene_identifier_base"] in genes_placebo
            ) else (row["rank_fold_p"]),
        axis="columns", # apply function to each row
    )
    # Sort rows within table.
    table.sort_values(
        by=[
            "rank_fold_p_placebo_null",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Write information to file.
    # Collect information.
    # Collections of files.
    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_adipose_32")] = table
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
