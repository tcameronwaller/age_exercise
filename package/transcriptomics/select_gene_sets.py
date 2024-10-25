"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'select_gene_sets' is part of the 'transcriptomics' package within
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
import exercise.transcriptomics.organize_signal as exrosig

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.


def preinitialize_directories(
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    restore=None,
    report=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        restore (bool): whether to remove previous versions of data
        report (bool): whether to print reports

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_directory_dock
    paths["in_data"] = os.path.join(
        paths["dock"], "in_data", str(project), str(routine),
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters", str(project), str(routine),
    )
    paths["in_parameters_private"] = os.path.join(
        paths["dock"], "in_parameters_private", str(project), str(routine),
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
    paths["out_data_overall"] = os.path.join(
        paths["out_procedure"], "data",
    )
    paths["out_data_overall_text"] = os.path.join(
        paths["out_data_overall"], "text",
    )
    paths["out_data_overall_pickle"] = os.path.join(
        paths["out_data_overall"], "pickle",
    )
    paths["out_plot_overall"] = os.path.join(
        paths["out_procedure"], "plot",
    )
    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"],
        paths["out_data_overall"],
        paths["out_data_overall_text"],
        paths["out_data_overall_pickle"],
        paths["out_plot_overall"],
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
        print("module: exercise.transcriptomics.select_gene_sets.py")
        print("function: preinitialize_directories()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


def initialize_directories(
    project=None,
    routine=None,
    procedure=None,
    tissue=None,
    group=None,
    name_instance=None,
    path_directory_dock=None,
    restore=None,
    report=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for
            selection of samples in cohort and definition of analysis
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        restore (bool): whether to remove previous versions of data
        report (bool): whether to print reports

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    ##########
    # Preinitialize directories before parallel branches.
    paths = preinitialize_directories(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=False,
        report=report,
    )

    # Define paths to directories.
    paths["out_tissue"] = os.path.join(
        paths["out_procedure"], str(tissue),
    )
    paths["out_group"] = os.path.join(
        paths["out_tissue"], str(group),
    )
    paths["out_instance"] = os.path.join(
        paths["out_group"], str(name_instance),
    )
    #paths["out_test"] = os.path.join(
    #    paths["out_instance"], "test",
    #)
    paths["out_data"] = os.path.join(
        paths["out_instance"], "data",
    )
    paths["out_plot"] = os.path.join(
        paths["out_instance"], "plot",
    )
    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        #paths["out_procedure"], # omit to avoid conflict in parallel branches
        paths["out_instance"],
        paths["out_data"],
        paths["out_plot"],
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
        print("module: exercise.transcriptomics.select_gene_sets.py")
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


def define_column_types_table_deseq2():
    """
    Defines the names and types of variable columns within table from DESeq2.

    Review: TCW; 7 August 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["identifier_gene"] = "string"
    types_columns["baseMean"] = "float32"
    types_columns["log2FoldChange"] = "float32"
    types_columns["lfcSE"] = "float32"
    types_columns["stat"] = "float32"
    types_columns["pvalue"] = "float64"
    types_columns["padj"] = "float64"
    types_columns["gene_identifier"] = "string"
    types_columns["gene_identifier_base"] = "string"
    types_columns["gene_name"] = "string"
    types_columns["gene_type"] = "string"
    types_columns["gene_chromosome"] = "string"
    # Return information.
    return types_columns


def define_column_types_table_gene_emphasis():
    """
    Defines the names and types of variable columns within table from DESeq2.

    Review: TCW; 7 August 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["identifier"] = "string"
    types_columns["gene_identifier"] = "string"
    types_columns["gene_name"] = "string"
    types_columns["gene_type"] = "string"
    types_columns["gene_chromosome"] = "string"
    types_columns["aerobic_acute"] = "float32"
    types_columns["aerobic_chronic"] = "float32"
    types_columns["resistance_acute"] = "float32"
    types_columns["resistance_chronic"] = "float32"
    types_columns["inactivity"] = "float32"
    types_columns["WP_INSULIN_SIGNALING_IN_ADIPOCYTES_NORMAL_CONDITION"] = (
        "float32"
    )
    # Return information.
    return types_columns


def read_source(
    tissue=None,
    group=None,
    name_instance=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for
            selection of samples in cohort and definition of analysis
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
    path_file_table_deseq2 = os.path.join(
        paths["out_routine"], "deseq2", tissue, group,
        str("table_result_deseq2_" + name_instance + ".tsv"),
    )
    if False:
        path_file_table_gene_emphasis = os.path.join(
            paths["in_parameters_private"], "gene_sets_emphasis",
            "table_gene_reference_adipose.tsv",
        )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of fold changes for genes with differential expression from DESeq2.
    types_columns = define_column_types_table_deseq2()
    pail["table_deseq2"] = pandas.read_csv(
        path_file_table_deseq2,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    if False:
        # Table of information about genes for special emphasis.
        types_columns = define_column_types_table_gene_emphasis()
        table_gene_emphasis = pandas.read_csv(
            path_file_table_gene_emphasis,
            sep="\t",
            header=0,
            dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
        if False:
            pail["table_gene_emphasis"] = table_gene_emphasis.loc[
                (
                    (table_gene_emphasis["aerobic_acute"] == 1) |
                    (table_gene_emphasis["aerobic_chronic"] == 1) |
                    (table_gene_emphasis["resistance_acute"] == 1) |
                    (table_gene_emphasis["resistance_chronic"] == 1) |
                    (table_gene_emphasis["inactivity"] == 1)
                ), :
            ]
        if True:
            set = "inclusion"
            pail["table_gene_emphasis"] = table_gene_emphasis.loc[
                (
                    (table_gene_emphasis[set] == 1)
                ), :
            ]

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.select_gene_sets.py")
        print("function: read_source()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=4)
        count_rows = (pail["table_deseq2"].shape[0])
        count_columns = (pail["table_deseq2"].shape[1])
        print("deseq2 table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_deseq2"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail





##########
# 3. Organize from source the information about differential expression of
#    genes.


def define_column_sequence_table_change_deseq2():
    """
    Defines the columns in sequence within table.

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "identifier_gene",
        #"baseMean",
        "fold_change_log2",
        "fold_change_log2_standard_error",
        #"stat",
        "p_value",
        "p_value_fill",
        "p_value_negative_log10",
        "q_value",
        "q_value_fill",
        "q_value_negative_log10",
        "q_value_check",
        "q_significance_check",
        "gene_identifier_base",
        "rank_fold_p",
        "gene_identifier",
        #"gene_identifier_base",
        "gene_name",
        "gene_type",
        "gene_chromosome",
    ]
    # Return information.
    return columns_sequence


def organize_table_change_deseq2(
    table=None,
    columns_sequence=None,
    tissue=None,
    name_instance=None,
    report=None,
):
    """
    Organizes information in table about differential expression of genes.

    arguments:
        table (object): Pandas data-frame table of information about genes
            that demonstrate differential expression
        columns_sequence (list<str>): names of columns in sequence by which to
            filter and sort columns in table
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        name_instance (str): name for set of samples and parameters in the
            analysis of differential expression
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_change = table.copy(deep=True)
    # Copy other information.
    columns_sequence = copy.deepcopy(columns_sequence)

    # Translate names of columns.
    translations = dict()
    translations["log2FoldChange"] = "fold_change_log2"
    translations["lfcSE"] = "fold_change_log2_standard_error"
    translations["pvalue"] = "p_value"
    translations["padj"] = "q_value"
    table_change.rename(
        columns=translations,
        inplace=True,
    )

    # Replace values of zero for p-value and q-value with reasonable
    # approximation of precision.
    # Set threshold on very small values of p-value and q-value.
    # DESeq2 returns values of zero for p-value and q-value when the actual
    # value is smaller than a threshold of approximately 1E-250.
    #table_change["p_value_fill"] = table_change["p_value"].replace(
    #    to_replace=0.0,
    #    value=1E-250,
    #)
    #table_change["q_value_fill"] = table_change["q_value"].replace(
    #    to_replace=0.0,
    #    value=1E-250,
    #)
    table_change = table_change.loc[
        (
            (pandas.notna(table_change["p_value"])) &
            (pandas.notna(table_change["q_value"])) &
            (table_change["p_value"] >= float(0)) &
            (table_change["q_value"] >= float(0))
        ), :
    ]
    table_change["p_value_fill"] = table_change.apply(
        lambda row:
            2.23E-308 if (float(row["p_value"]) < 2.23E-308) else row["p_value"],
        axis="columns", # apply function to each row
    )
    table_change["q_value_fill"] = table_change.apply(
        lambda row:
            2.23E-308 if (float(row["q_value"]) < 2.23E-308) else row["q_value"],
        axis="columns", # apply function to each row
    )
    # Filter rows in table for selection of non-missing values for fold change.
    table_change = table_change.loc[
        (
            (pandas.notna(table_change["fold_change_log2"])) &
            (pandas.notna(table_change["fold_change_log2_standard_error"])) &
            (pandas.notna(table_change["p_value_fill"])) &
            (pandas.notna(table_change["q_value_fill"])) &
            (table_change["p_value_fill"] > float(0)) &
            (table_change["q_value_fill"] > float(0))
        ), :
    ]

    # Calculate the negative, base-ten logarithm of the p-value for
    # differential expression of each gene.
    #table_change["p_value_negative_log10"] = numpy.log10(
    #    table_change["p_value"]
    #)
    table_change["p_value_negative_log10"] = table_change.apply(
        lambda row: (-1*math.log(row["p_value_fill"], 10)),
        axis="columns", # apply function to each row
    )
    table_change["q_value_negative_log10"] = table_change.apply(
        lambda row: (-1*math.log(row["q_value_fill"], 10)),
        axis="columns", # apply function to each row
    )

    ##########
    # Calculate Benjamini-Hochberg q-values for False-Discovery Rate (FDR).
    # Calculate q-values across all comparisons in table.
    # FDR 5% (q <= 0.05).
    # The goal is to check the adjusted p-values from DESeq2.
    table_change = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=0.05, # alpha; family-wise error rate
        name_column_p_value="p_value_fill",
        name_column_q_value="q_value_check",
        name_column_significance="q_significance_check",
        table=table_change,
    )

    # Calculate a metric by which to rank genes with differential expression
    # for subsequent analysis by gene set enrichment (GSEA).
    # Use the p-value rather than the q-value since the q-value tends to have
    # ties between groups of comparisons.
    table_change["rank_fold_p"] = table_change.apply(
        lambda row: (
            (row["fold_change_log2"])*(row["p_value_negative_log10"])
        ),
        axis="columns", # apply function to each row
    )

    # Sort rows within table.
    table_change.sort_values(
        by=[
            "rank_fold_p",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Filter and sort columns within table.
    table_change = porg.filter_sort_table_columns(
        table=table_change,
        columns_sequence=columns_sequence,
        report=report,
    )
    # Filter rows within table.
    table_significance = table_change.loc[
        (
            (table_change["q_value_fill"] < 0.05)
        ), :
    ].copy(deep=True)

    # Organize indices in table.
    table_change.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Collect information.
    pail = dict()
    pail["table"] = table_change
    pail["table_significance"] = table_significance
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: organize_table_change_deseq2()")
        putly.print_terminal_partition(level=5)
        print("tissue: " + tissue)
        print("name_instance: " + name_instance)
        putly.print_terminal_partition(level=4)
        count_rows = (pail["table_significance"].shape[0])
        count_columns = (pail["table_significance"].shape[1])
        print("table of genes with significant differential expression: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_significance"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 4. Select sets of genes with differential expression.


def select_sets_differential_expression_gene(
    table=None,
    column_identifier=None,
    column_name=None,
    column_fold_change=None,
    column_significance=None,
    threshold_fold_change=None,
    threshold_significance=None,
    tissue=None,
    name_instance=None,
    report=None,
):
    """
    Selects sets of genes applying filters at specific thresholds on their
    statistics for differential expression.

    arguments:
        table (object): Pandas data-frame table of information about genes
            that demonstrate differential expression
        column_identifier (str): name of column in table for the unique
            identifier corresponding to the fold change
        column_name (str): name of column in table for the name corresponding
            to the fold change
        column_fold_change (str): name of column in table on which to apply the
            threshold for the fold change
        column_significance (str): name of column in table on which to apply
            the threshold for significance, corresponding to the p-value or
            q-value for the estimate of fold change
        threshold_fold_change (float): value for threshold on fold change
            (fold change > threshold) that is on the same scale, such as
            base-two logarithm, as the actual values themselves
        threshold_significance (float): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        name_instance (str): name for set of samples and parameters in the
            analysis of differential expression
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Filter rows in table for selection of sets of genes that demonstrate
    # specific characteristics of differential expression.
    pail_threshold = porg.segregate_fold_change_values_by_thresholds(
        table=table,
        column_fold_change=column_fold_change,
        column_significance=column_significance,
        threshold_fold_change=threshold_fold_change,
        threshold_significance=threshold_significance,
        report=False,
    )
    # Extract identifiers genes with differential expression beyond thresholds.
    genes_change = copy.deepcopy(
        pail_threshold["table_pass_any"][column_identifier].to_list()
    )
    genes_up = copy.deepcopy(
        pail_threshold["table_pass_up"][column_identifier].to_list()
    )
    genes_down = copy.deepcopy(
        pail_threshold["table_pass_down"][column_identifier].to_list()
    )

    # Collect information.
    pail = dict()
    pail["genes_change"] = genes_change
    pail["genes_up"] = genes_up
    pail["genes_down"] = genes_down
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: select_sets_differential_expression_gene()")
        putly.print_terminal_partition(level=5)
        print("tissue: " + tissue)
        print("name_instance: " + name_instance)
        putly.print_terminal_partition(level=4)
        count_both = (len(pail["genes_change"]))
        count_up = (len(pail["genes_up"]))
        count_down = (len(pail["genes_down"]))
        print("count of genes beyond thresholds: " + str(count_both))
        print("count of those with positive change: " + str(count_up))
        print("count of those with negative change: " + str(count_down))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 6. Create chart to represent fold changes and write to file.


def plot_write_chart_fold_change_volcano(
    table=None,
    column_identifier=None,
    column_name=None,
    column_plot_fold=None,
    column_plot_p=None,
    column_fold_change=None,
    column_significance=None,
    threshold_fold_change=None,
    threshold_significance=None,
    identifiers_emphasis=None,
    tissue=None,
    name_instance=None,
    paths=None,
    report=None,
):
    """
    Create chart representation of fold change and write to file.

    Selection of genes come from the study below.
    Pillon et al, Nature Communications, 2020; (PubMed:31980607)

    arguments:
        table (object): Pandas data-frame table of information about genes
            that demonstrate differential expression
        column_identifier (str): name of column in table for the unique
            identifier corresponding to the fold change
        column_name (str): name of column in table for the name corresponding
            to the fold change
        column_plot_fold (str): name of column in table for values of fold
            change for representation on the plot
        column_plot_p (str): name of column in table for values of p-value or
            representation of significance for representation on the plot
        column_fold_change (str): name of column in table on which to apply the
            threshold for the fold change
        column_significance (str): name of column in table on which to apply
            the threshold for significance, corresponding to the p-value or
            q-value corresponding to the estimate of fold change
        threshold_fold_change (float): value for threshold on fold change
            (fold change > threshold) that is on the same scale, such as
            base-two logarithm, as the actual values themselves
        threshold_significance (float): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        identifiers_emphasis (list<str>): identifiers corresponding to a
            special selection of fold changes for which to emphasize points on
            chart and for which to create text labels adjacent to the points
            on the chart
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        name_instance (str): name for set of samples and parameters in the
            analysis of differential expression
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    #column_plot_fold="fold_change_log2",
    #column_plot_p="p_value_negative_log10",
    #column_significance_q="q_value_fill",
    #threshold_fold=math.log(float(1.1), 2), # base two logarithm
    #threshold_q=(-1 * math.log(float(0.05), 10)), # negative base ten logarithm


    # Organize parameters.
    name_figure = name_instance
    #path_directory = paths["out_plot"]
    path_directory = paths["out_plot_overall"]

    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_scatter_fold_change_volcano(
        table=table,
        column_identifier=column_identifier,
        column_name=column_name,
        column_plot_fold=column_plot_fold,
        column_plot_p=column_plot_p,
        column_fold_change=column_fold_change,
        column_significance=column_significance,
        threshold_fold_change=threshold_fold_change,
        threshold_significance=threshold_significance,
        identifiers_emphasis=identifiers_emphasis,
        lines_thresholds=True,
        emphasis_label=True,
        count_label=True,
        minimum_abscissa=(
            (numpy.nanmin(table[column_plot_fold].to_numpy())) - 0.5
        ),
        maximum_abscissa=(
            (numpy.nanmax(table[column_plot_fold].to_numpy())) + 0.5
        ),
        minimum_ordinate=-0.1,
        maximum_ordinate=(
            (numpy.nanmax(table[column_plot_p].to_numpy())) + 1
        ),
        title_abscissa="log2(fold change)",
        title_ordinate="-1*log10(p-value)", # "-1*log10(B.-H. FDR q-value)"
        size_title_abscissa="eight", # ten
        size_title_ordinate="eight", # ten
        size_label_abscissa="twelve", # multi-panel: ten; individual: twelve
        size_label_ordinate="twelve", # multi-panel: ten; individual: twelve
        size_label_emphasis="seventeen",
        size_label_count="ten",
        aspect="landscape", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Write figure to file.
    pplot.write_product_plot_figure(
        figure=figure,
        format="jpg", # jpg, png, svg
        resolution=150,
        name_file=name_figure,
        path_directory=path_directory,
    )
    # Return information.
    return figure


# TODO: TCW; 8 August 2024
# NEXT STEPS
# 1. filter genes by fold change and by p-value
#    - potentially consider the 95%CI of fold change for this threshold filter
#    - this would look funny on the volcano plot...
#    - maybe the volcano plot could look normal, except that the genes that pass
#    - filter by 95% CI could be the ones that I actually highlight in orange
# 2.


##########
# Downstream visualizations

# volcano plots of fold changes
# Venn diagrams of sets of differentially expressed genes (up or down)
# heatmaps of gene signals across persons
# pairwise correlation matrices with hierarchical clustering
#



###############################################################################
# Procedure


##########
# Control procedure within branch for parallelization.


def control_procedure_part_branch(
    tissue=None,
    group=None,
    name_instance=None,
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for
            selection of samples in cohort and definition of analysis
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # 1. Initialize directories for read of source and write of product files.
    paths = initialize_directories(
        project=project,
        routine=routine,
        procedure=procedure,
        tissue=tissue,
        group=group,
        name_instance=name_instance,
        path_directory_dock=path_directory_dock,
        restore=True,
        report=report,
    )

    ##########
    # 2. Read source information from file.
    pail_source = read_source(
        tissue=tissue,
        group=group,
        name_instance=name_instance,
        paths=paths,
        report=report,
    )

    ##########
    # 3. Organize from source the information about differential expression of
    #    genes.
    columns_sequence = define_column_sequence_table_change_deseq2()
    pail_organization = organize_table_change_deseq2(
        table=pail_source["table_deseq2"],
        columns_sequence=columns_sequence,
        tissue=tissue,
        name_instance=name_instance,
        report=report,
    )
    #pail_organization["table_change"]

    ##########
    # 4. Select sets of genes with differential expression.
    #threshold_p=(-1 * math.log(float(0.05), 10)), # negative base ten logarithm
    pail_selection = select_sets_differential_expression_gene(
        table=pail_organization["table"],
        column_identifier="gene_identifier_base",
        column_name="gene_name",
        column_fold_change="fold_change_log2",
        column_significance="q_value_fill",
        threshold_fold_change=math.log(float(1), 2), # base two logarithm
        threshold_significance=float(0.05),
        tissue=tissue,
        name_instance=name_instance,
        report=report,
    )

    ##########
    # 5. Prepare rank lists of genes for analysis by gene set enrichment
    #    analysis.
    if False:
        pail_rank = rank_list_gene(
            table=pail_organization["table"],
            column_identifier="gene_identifier",
            column_name="gene_name",
            column_rank="rank_fold_p",
            tissue=tissue,
            name_instance=name_instance,
            report=report,
        )

    ##########
    # 6. Create chart to represent fold changes and write to file.
    identifiers_emphasis = []
    plot_write_chart_fold_change_volcano(
        table=pail_organization["table"],
        column_identifier="gene_identifier_base",
        column_name="gene_name",
        column_plot_fold="fold_change_log2",
        column_plot_p="p_value_negative_log10",
        column_fold_change="fold_change_log2",
        column_significance="q_value_fill",
        threshold_fold_change=math.log(float(1), 2), # base two logarithm
        threshold_significance=float(0.05),
        identifiers_emphasis=identifiers_emphasis,
        tissue=tissue,
        name_instance=name_instance,
        paths=paths,
        report=report,
    )


    ##########
    # Collect information.
    # Collections of files.
    pail_write_lists = dict()
    pail_write_lists[str("genes_change")] = (
        pail_selection["genes_change"]
    )
    pail_write_lists[str("genes_up")] = (
        pail_selection["genes_up"]
    )
    pail_write_lists[str("genes_down")] = (
        pail_selection["genes_down"]
    )
    pail_write_tables = dict()
    pail_write_tables[str("table_" + name_instance)] = (
        pail_organization["table"]
    )

    ##########
    # Write product information to file.
    putly.write_lists_to_file_text(
        pail_write=pail_write_lists,
        path_directory=paths["out_data"],
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_data_overall_text"],
        reset_index=False,
        write_index=True,
        type="text",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_data_overall_pickle"],
        reset_index=False,
        write_index=True,
        type="pickle",
    )
    pass


##########
# Manage parallelization.


def control_parallel_instance(
    instance=None,
    parameters=None,
):
    """
    Control procedure to organize within tables the information about genetic
    correlations from LDSC.

    arguments:
        instance (dict): parameters specific to current instance
            tissue (str): name of tissue that distinguishes study design and
                set of relevant samples, either 'adipose' or 'muscle'
            group (str): name of a group of analyses
            name_instance (str): name of instance set of parameters for
                selection of samples in cohort and definition of analysis
        parameters (dict): parameters common to all instances
            project (str): name of project
            routine (str): name of routine, either 'transcriptomics' or
                'proteomics'
            procedure (str): name of procedure, a set or step in the routine
                process
            path_directory_dock (str): path to dock directory for procedure's
                source and product directories and files
            report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Extract parameters.
    # Extract parameters specific to each instance.
    tissue = instance["tissue"]
    group = instance["group"]
    name_instance = instance["name_instance"]
    # Extract parameters common across all instances.
    project = parameters["project"]
    routine = parameters["routine"]
    procedure = parameters["procedure"]
    path_directory_dock = parameters["path_directory_dock"]
    report = parameters["report"]

    ##########
    # Control procedure with split for parallelization.
    control_procedure_part_branch(
        tissue=tissue,
        group=group,
        name_instance=name_instance,
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    pass


def control_parallel_instances(
    instances=None,
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control procedure for parallel instances.

    arguments:
        instances (list<dict>): parameters to control individual instances in
            parallel
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Collect parameters common across all instances.
    parameters = dict()
    parameters["project"] = project
    parameters["routine"] = routine
    parameters["procedure"] = procedure
    parameters["path_directory_dock"] = path_directory_dock
    parameters["report"] = report

    # Execute procedure iteratively with parallelization across instances.
    if True:
        prall.drive_procedure_parallel(
            function_control=(
                control_parallel_instance
            ),
            instances=instances,
            parameters=parameters,
            cores=5,
            report=True,
        )
    else:
        # Execute procedure directly for testing.
        control_parallel_instance(
            instance=instances[1],
            parameters=parameters,
        )
    pass


##########
# Call main procedure.


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
    project="exercise"
    routine="transcriptomics"
    procedure="select_gene_sets"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.select_gene_sets.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("routine: " + str(routine))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Preinitialize directories before parallel branches.
    paths = preinitialize_directories(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=True,
        report=report,
    )

    ##########
    # Read and organize parameters for parallel instances.
    instances = exrosig.read_organize_source_parameter_instances(
        paths=paths,
        report=report,
    )

    ##########
    # Control procedure for parallel instances.
    control_parallel_instances(
        instances=instances,
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass


###############################################################################
# End
