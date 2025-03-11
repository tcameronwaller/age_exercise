"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'scratch' is part of the 'transcriptomics' package within
the 'exercise' package.

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

# consider changing the name of this module to "select_changes"


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
import age_exercise.transcriptomics.organize_signal as exrosig
import age_exercise.transcriptomics.select_gene_sets as exrosel
import age_exercise.proteomics.organize_subject as aexpr_sub

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.
# There is a hierarchy in these functions to initialize directories to manage
# the hierarchical tree structure of sub-procedures.



##########
# 2. Read source information from file.


##########
# Plot charts.


###############################################################################
# Procedure


##########
# Execute main procedure.


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
    routine="main"
    procedure="scratch"
    report = True

    # Define paths to directories.
    # Broad.
    paths = dict()
    paths["dock"] = path_directory_dock
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    # Define paths to child files.
    path_file_table_signal = os.path.join(
        paths["out_project"], "transcriptomics", "organize_signal", "whole",
        "preparation", str("table_signal_scale_" + "adipose" + ".pickle"),
    )
    # Read information from file.
    table_signal = pandas.read_pickle(
        path_file_table_signal,
    )
    print(table_signal)

    ##########
    # Copy information.
    table = table_signal.copy(deep=True)
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table.set_index(
        "identifier_gene",
        append=False,
        drop=True,
        inplace=True
    )
    table.columns.rename(
        "sample",
        inplace=True,
    ) # single-dimensional index
    # Transpose table.
    table_transpose = table.transpose(copy=True)
    # Organize indices in table.
    table_transpose.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_transpose.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index
    # Copy information.
    table_signal_transpose = table_transpose.copy(deep=True)

    print(table_signal_transpose)

    ##########
    # Copy information.
    table = table_signal_transpose.copy(deep=True)
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Convert variables to lowest float type.
    columns = copy.deepcopy(table.columns.to_list())
    putly.convert_table_columns_variables_types_float(
        columns=columns,
        table=table,
    )
    # Calculate pairwise correlations between columns.
    table_correlation = table.corr(
        method="spearman",
        min_periods=10,
        numeric_only=True,
    )
    # Organize indices in table.
    table_correlation.rename_axis(
        "identifier_gene",
        axis="index",
        inplace=True,
    )
    table_transpose.columns.rename(
        "identifier_gene",
        inplace=True,
    ) # single-dimensional index
    # Copy information.
    table_signal_correlation = table_correlation.copy(deep=True)
    print(table_signal_correlation)

    ##########
    # Cluster.
    # Copy information.
    table = table_signal_correlation.copy(deep=True)
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster rows in table.
    table = porg.cluster_table_rows(
        table=table,
    )
    # Cluster columns in table.
    table = porg.cluster_table_columns(
        table=table,
    )
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Copy information.
    table_signal_correlation_cluster = table.copy(deep=True)
    print(table_signal_correlation_cluster)


    if False:
        heatmap = (
            aexpr_sub.plot_heatmap_features_observations_labels(
                table=table_signal_correlation_cluster,
                index_columns="group_observations",
                index_rows=index_features,
                report=False,
        ))


    # Path to out directory.
    path_scratch_out = os.path.join(
        paths["dock"], str("out_" + project), "scratch",
    )
    # Create directories.
    putly.create_directories(
        path=path_scratch_out,
    )


    ##########
    # Collect information.
    # Collections of files.
    pail_write_data = dict()
    pail_write_data[str("table_signal_correlation_cluster")] = (
        table_signal_correlation_cluster
    )
    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_data,
        path_directory=path_scratch_out,
        reset_index_rows=False,
        write_index_rows=True,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_data,
        path_directory=path_scratch_out,
        reset_index_rows=None,
        write_index_rows=None,
        write_index_columns=None,
        type="pickle",
        delimiter=None,
        suffix=".pickle",
    )


    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.scratch.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("routine: " + str(routine))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass



    pass


###############################################################################
# End
