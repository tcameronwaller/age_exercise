"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'compare_sets_groups' is part of the 'transcriptomics' package
within the 'exercise' package.

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

# TODO: TCW; 15 October 2024
# TODO: simplify this process of organizing the signal information
# TODO: It might help to separate the preparation of the overall signal table
# TODO: from the stratified signal tables for individual analyses.



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
import exercise.transcriptomics.organize_sample as extr_sample
import exercise.transcriptomics.organize_signal as extr_signal
import exercise.transcriptomics.select_gene_sets as extr_select


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
        "sets_gene_2024-11-21",
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
        print("module: exercise.transcriptomics.compare_sets_groups.py")
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


# 2.2. Read and organize main source information from file.


def define_column_types_table_demonstration():
    """
    Defines the types of variables for columns in a table.

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types in columns within table.
    types_columns = dict()
    types_columns["observation"] = "string"
    types_columns["group"] = "string"
    types_columns["feature_1"] = "float"
    types_columns["feature_2"] = "float"
    types_columns["feature_3"] = "float"
    types_columns["feature_4"] = "float"
    types_columns["feature_5"] = "float"
    # Return information.
    return types_columns


def read_source(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (str): name of tissue, either 'muscle' or 'adipose', that
            distinguishes the study design and its relevant data
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
    path_file_table_demonstration = os.path.join(
        paths["in_demonstration"], "partner",
        "table_columns_features_rows_observations_groups.tsv",
    )
    path_file_table_sample = os.path.join(
        paths["out_routine"], "organize_sample", "data",
        "table_sample.pickle",
    )
    path_file_table_gene = os.path.join(
        paths["out_routine"], "organize_signal", "whole", "preparation",
        str("table_gene_" + tissue + ".pickle"),
    )
    path_file_table_signal = os.path.join(
        paths["out_project"], "transcriptomics", "organize_signal", "whole",
        "preparation", str("table_signal_scale_" + tissue + ".pickle"),
    )

    # Collect information.
    pail = dict()
    # Read information from file.
    #types_columns = define_column_types_table_demonstration()
    pail["table_demonstration"] = pandas.read_csv(
        path_file_table_demonstration,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_sample"] = pandas.read_pickle(
        path_file_table_sample,
    )
    pail["table_gene"] = pandas.read_pickle(
        path_file_table_gene,
    )
    pail["table_signal"] = pandas.read_pickle(
        path_file_table_signal,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.compare_sets_groups.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("tissue: " + str(tissue))
        putly.print_terminal_partition(level=5)
        count_rows = (pail["table_demonstration"].shape[0])
        count_columns = (pail["table_demonstration"].shape[1])
        print("demonstration table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_demonstration"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# 2.3. Read and organize information about parameters for instances.


def define_column_types_table_parameter_instances():
    """
    Defines the variable types of columns within table for attributes of
    samples.

    Review: TCW; 14 August 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion"] = "string" # "int32"
    types_columns["tissue"] = "string"
    types_columns["sort"] = "int32" # "int32"
    types_columns["sort_group"] = "int32" # "int32"
    types_columns["group"] = "string"
    types_columns["instance"] = "string"
    types_columns["selection_samples_primary"] = "string"
    types_columns["selection_samples_secondary"] = "string"
    types_columns["name_set_gene_inclusion"] = "string"
    types_columns["names_sets_gene_cluster"] = "string"
    types_columns["names_sets_gene_allocation"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


def read_organize_source_parameter_instances(
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
        (list<dict>): collections of source information about instances for
           selection of information that is specific to separate analyses

    """

    # Define paths to parent directories.
    #paths["in_data"]
    #paths["in_parameters"]
    #paths["in_parameters_private"]

    # Define paths to child files.
    path_file_table_parameter = os.path.join(
        paths["in_parameters_private"], "exercise", "transcriptomics",
        "table_comparisons_gene_sets_between_sample_groups.tsv",
    )

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_column_types_table_parameter_instances()
    table = pandas.read_csv(
        path_file_table_parameter,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Collect information.
    instances = list()
    for index, row in table.iterrows():
        if (int(row["inclusion"]) == 1):
            # Collect information.
            pail = dict()
            pail["sort"] = str(row["sort"])
            pail["sort_group"] = str(row["sort_group"])
            pail["group"] = str(row["group"])
            pail["name_group"] = "_".join([
                str(row["sort_group"]),
                str(row["group"])
            ])
            pail["tissue"] = str(row["tissue"])
            pail["instance"] = str(row["instance"])
            pail["name_instance"] = "_".join([
                str(row["tissue"]),
                str(row["sort"]),
                str(row["instance"])
            ])
            # set: selection_samples_primary
            pail["selection_samples_primary"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["selection_samples_primary"],
                )
            )["features_values"]
            # set: selection_samples_secondary
            pail["selection_samples_secondary"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["selection_samples_secondary"],
                )
            )["features_values"]
            pail["names_sets_gene_inclusion"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["names_sets_gene_inclusion"],
                )
            )["features_values"]
            pail["names_sets_gene_cluster"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["names_sets_gene_cluster"],
                )
            )["features_values"]
            pail["names_sets_gene_allocation"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["names_sets_gene_allocation"],
                )
            )["features_values"]
            pail["sequence_sets_gene_inclusion"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["names_sets_gene_inclusion"],
                )
            )["features"]
            pail["sequence_sets_gene_cluster"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["names_sets_gene_cluster"],
                )
            )["features"]
            pail["sequence_sets_gene_allocation"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["names_sets_gene_allocation"],
                )
            )["features"]
            # Extract names of columns corresponding to feature variables for
            # which to calculate tertiles.
            #columns_tertile = extract_source_columns_for_tertiles(
            #    selection_samples_set=pail["selection_samples_secondary"],
            #    report=report,
            #)
            # Collect names of unique features or columns relevant to current
            # instance.
            features = list()
            dictionaries = [
                "selection_samples_primary",
                "selection_samples_secondary",
            ]
            for dictionary in dictionaries:
                if pail[dictionary] is not None:
                    features.extend(list(pail[dictionary].keys()))
                    pass
                pass
            # Collect unique names of features.
            features_unique = putly.collect_unique_elements(
                elements=features,
            )
            pail["features"] = features_unique
            instances.append(pail)
            pass
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.compare_sets_groups.py")
        print("function: read_organize_source_parameter_instances()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        #print("instances:")
        #print(instances)
        print("instance[0]:")
        print(instances[0])
        pass
    # Return information.
    return instances


##########
# 3. Organize information from source.


def prepare_sets_gene(
    names_sets_gene_inclusion=None,
    names_sets_gene_cluster=None,
    names_sets_gene_allocation=None,
    genes_available=None,
    paths=None,
    report=None,
):
    """
    Prepare sets of genes for analysis and description.

    arguments:
        names_sets_gene_inclusion (dict<list<str>>): names of sets of genes to
            include within a single custom combination set with custom name;
            this is the set of total genes for inclusion in description and
            visual representation
        names_sets_gene_cluster (dict<list<str>>): names of sets of genes to
            include within custom combination sets with custom names; these are
            the sets or groups within which to cluster signals for individual
            genes across samples; the sets must include all of the total
            inclusion genes, and each set must be mutually exclusive
        names_sets_gene_allocation (dict<list<str>>): names of sets of genes to
            include within custom combination sets with custom names; these are
            the sets or groups within which to allocate the total inclusion
            genes for visual representation
        genes_available (list<str>): identifiers of genes with available
            signals
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    ##########
    # Copy information.
    names_sets_gene_inclusion = copy.deepcopy(names_sets_gene_inclusion)
    names_sets_gene_cluster = copy.deepcopy(names_sets_gene_cluster)
    names_sets_gene_allocation = copy.deepcopy(names_sets_gene_allocation)
    genes_available = copy.deepcopy(genes_available)

    ##########
    # Genes for inclusion.
    sets_gene_inclusion = extr_select.read_extract_combine_custom_sets_genes(
        names_sets=names_sets_gene_inclusion,
        genes_available=genes_available,
        path_directory=paths["in_sets_gene"],
        report=report,
    )
    # From the sets of genes for inclusion, only use the combination custom set
    # with name "main".

    ##########
    # Genes for clustering.
    sets_gene_cluster = extr_select.read_extract_combine_custom_sets_genes(
        names_sets=names_sets_gene_cluster,
        genes_available=genes_available,
        path_directory=paths["in_sets_gene"],
        report=report,
    )

    ##########
    # Genes for allocation.
    sets_gene_allocation = extr_select.read_extract_combine_custom_sets_genes(
        names_sets=names_sets_gene_allocation,
        genes_available=genes_available,
        path_directory=paths["in_sets_gene"],
        report=report,
    )

    # Collect information.
    pail = dict()
    pail["sets_gene_inclusion"] = sets_gene_inclusion
    pail["sets_gene_cluster"] = sets_gene_cluster
    pail["sets_gene_allocation"] = sets_gene_allocation
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.compare_sets_groups.py")
        print("function: prepare_sets_gene()")
        putly.print_terminal_partition(level=5)
        print("only use custom combination set 'main' for total inclusion.")
        count_inclusion = (len(pail["sets_gene_inclusion"]["main"]))
        print(
            "count of genes for total inclusion: " + str(count_inclusion)
        )
        putly.print_terminal_partition(level=5)

        pass
    # Return information.
    return pail


def determine_gene_sets_allocation(
    genes_inclusion=None,
    sets_gene_allocation=None,
    index_genes=None,
    column_other=None,
    columns_sequence=None,
    report=None,
):
    """
    Determine the sets to which each relevant gene belongs.

    arguments:
        genes_inclusion (list<str>): identifiers of total genes for inclusion
            in description and visual representation
        sets_gene_allocation (dict<list<str>>): names of sets of genes to
            include within custom combination sets with custom names; these are
            the sets or groups within which to allocate the total inclusion
            genes for visual representation
        index_genes (str): name for index corresponding to genes, which
            is a column in the original source table of signals
        column_other (str): name for column indication of allocation to none of
            the query sets
        columns_sequence (list<str>): names of columns in sequence by which to
            filter and sort columns in table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information.
    genes_inclusion = copy.deepcopy(genes_inclusion)
    sets_gene_allocation = copy.deepcopy(sets_gene_allocation)
    columns_sequence = copy.deepcopy(columns_sequence)

    # Determine whether to include column for other sets.
    if (
        (column_other is not None) and
        (len(column_other) > 0)
    ):
        columns_sequence.append(column_other)
        pass

    ##########
    # Determine set allocations for each relevant gene.

    # Iterate on relevant genes.
    # Collect information about set allocations for each gene.
    records = list()
    for gene in genes_inclusion:
        # Collect information.
        record = dict()
        record[index_genes] = gene
        if (
            (column_other is not None) and
            (len(column_other) > 0)
        ):
            record[column_other] = 1
        # Iterate on sets for allocation.
        for name_set in sets_gene_allocation.keys():
            if (gene in sets_gene_allocation[name_set]):
                record[name_set] = 1
                if (
                    (column_other is not None) and
                    (len(column_other) > 0)
                ):
                    record[column_other] = 0
            else:
                record[name_set] = 0
            pass
        # Collect information.
        records.append(record)
        pass
    # Create pandas data-frame table.
    table = pandas.DataFrame(data=records)

    # Filter and sort columns within table.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=False,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.compare_sets_groups.py")
        print("function: determine_gene_sets_allocation()")
        putly.print_terminal_partition(level=5)
        print(table)
        # Determine whether to describe column for other sets.
        if (
            (column_other is not None) and
            (len(column_other) > 0)
        ):
            # Filter rows within table.
            table_nonother = table.loc[
                (table[column_other] == 0), :
            ].copy(deep=True)
            putly.print_terminal_partition(level=5)
            print(table_nonother)
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def sort_table_rows_gene_sets_allocation(
    table_gene_sets=None,
    table_signal=None,
    index_genes=None,
    indices_rows_signal=None,
    report=None,
):
    """
    Sort genes to match those in table of signals.

    arguments:
        table_gene_sets (object): Pandas data-frame table
        table_signal (object): Pandas data-frame table
        index_genes (str): name for index corresponding to genes, which
            is a column in the original source table of signals
        indices_rows_signal (list<str>): names for columns in index across rows
            of signal table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_gene_sets = table_gene_sets.copy(deep=True)
    table_signal = table_signal.copy(deep=True)

    # Extract identifiers of genes.
    identifiers_gene = copy.deepcopy(table_signal.columns.to_list())
    for index in indices_rows_signal:
        identifiers_gene.remove(index)
        pass

    # Prepare indices for sort.
    reference_sort = dict(zip(identifiers_gene, range(len(identifiers_gene))))

    # Sort rows in table reference indices.
    table = porg.sort_table_rows_by_single_column_reference(
        table=table_gene_sets,
        index_rows=index_genes,
        column_reference=index_genes,
        column_sort_temporary="sort_temporary",
        reference_sort=reference_sort,
    )

    # Report.
    if report:
        # Extract index across rows for comparison.
        values_index_rows = copy.deepcopy(
            table[index_genes].to_list()
        )
        # Confirm that sets of indices from both sources are inclusive.
        check_inclusion = putly.compare_lists_by_mutual_inclusion(
            list_primary=identifiers_gene,
            list_secondary=values_index_rows,
        )
        # Confirm that sets of indices from both sources are identical across
        # their respective sequences.
        check_identity = putly.compare_lists_by_elemental_identity(
            list_primary=identifiers_gene,
            list_secondary=values_index_rows,
        )
        # Confirm that sets of indices from both sources are equal.
        check_equality = (
            identifiers_gene == values_index_rows
        )

        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.compare_sets_groups.py")
        print("function: sort_table_rows_gene_sets_allocation()")
        putly.print_terminal_partition(level=4)
        print(
            "confirm that gene identifiers across rows are identical to " +
            "reference"
        )
        print(
            "check coherence of gene identifiers"
        )
        print("check inclusion: " + str(check_inclusion))
        print("check identity: " + str(check_identity))
        print("check equality: " + str(check_equality))
        putly.print_terminal_partition(level=5)
        print("product table after sort")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# 4. Plot information to represent visually in charts.


def create_plot_chart_box(
    table=None,
    column_feature=None,
    column_group=None,
    report=None,
):
    """
    Create and plot a chart of the box type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 16 October 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement
        column_feature (str): name of column for values corresponding to a
            specific feature
        column_group (str): name of column in table to use for groups
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract names and indices of groups to preserve their original sequence.
    groups_sequence = copy.deepcopy(table[column_group].unique().tolist())
    sequence_groups = dict()
    index = 0
    for name in groups_sequence:
        sequence_groups[name] = index
        index += 1
        pass
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.set_index(
        [column_group],
        append=False,
        drop=True,
        inplace=True
    )
    # Split rows within table by factor columns.
    groups = table.groupby(
        level=column_group,
    )
    # Collect information.
    names_groups = list()
    values_groups = list()
    # Iterate on groups, apply operations, and collect information from each.
    for name_group, table_group in groups:
        # Copy information in table.
        table_group = table_group.copy(deep=True)
        # Extract information.
        values = table_group[column_feature].dropna().to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )
        # Collect information.
        names_groups.append(name_group)
        values_groups.append(values)
        pass
    # Sort information for groups and values.
    sequence_sort = list()
    for name in names_groups:
        sequence_sort.append(sequence_groups[name])
        pass
    names_groups_sort = [
        y for _,y in sorted(zip(sequence_sort, names_groups))
    ]
    values_groups_sort = [
        y for _,y in sorted(zip(sequence_sort, values_groups))
    ]

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_boxes_groups(
        values_groups=values_groups_sort,
        title_ordinate="scale-normal(gene signal)",
        title_abscissa="",
        titles_abscissa_groups=names_groups_sort,
        colors_groups=None,
        label_top_center="",
        label_top_left="",
        label_top_right="",
        aspect="landscape",
        orientation_box="vertical",
        axis_linear_minimum=0.0,
        fonts=fonts,
        colors=colors,
        report=report,
    )

    # Return information.
    return figure


def create_plot_chart_heatmap_individual(
    table_signal=None,
    table_feature_sets=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 17 October 2024

    arguments:
        table_signal (object): Pandas data-frame table of floating-point values
            of a signal corresponding features in columns across observations
            in rows
        table_feature_sets (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table of signals
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    table_signal_extract = table_signal.copy(deep=True)
    table_feature_sets = table_feature_sets.copy(deep=True)
    # Organize indices in table.
    table_signal_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_signal_extract.set_index(
        [index_rows, column_group,],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_signal_extract.to_numpy())
    value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    #features_sets_in_observations_groups

    # plot_heatmap_signal_label_features_groups_of_observations

    figure = pplot.plot_heatmap_signal_features_sets_in_observations_groups(
        table_signal=table_signal,
        table_feature_sets=table_feature_sets,
        format_table_signal=2, # 2: features in columns; observations, groups in rows
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=True,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        show_labels_ordinate=False,
        show_labels_abscissa=False,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        show_scale_bar=True, # whether to show scale bar on individual figures
        title_ordinate="",
        title_abscissa="",
        title_bar="gene signal (z-score)",
        size_title_ordinate="eight", # ten
        size_title_abscissa="eight", # ten
        size_label_ordinate="seventeen", # multi-panel: ten; individual: twelve
        size_label_abscissa="eleven", # multi-panel: ten; individual: twelve
        size_title_bar="twelve", # twelve
        size_label_bar="thirteen", # thirteen for whole; five for bar itself
        aspect="portrait", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )

    # Return information.
    return figure


def create_plot_chart_heatmap_mean(
    table=None,
    index_columns=None,
    index_rows=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 17 October 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values of a
            single, specific type of descriptive statistics (usually either
            mean or median) corresponding to groups of observations across
            columns and features across rows
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)
    table_extract = table.copy(deep=True)
    # Organize indices in table.
    table_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_extract.set_index(
        [index_rows],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_extract.to_numpy())
    value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heatmap_signal_label_features_observations(
        table=table,
        format_table=1, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        show_labels_ordinate=True,
        #show_labels_abscissa=True,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        show_scale_bar=True, # whether to show scale bar on individual figures
        title_ordinate="",
        title_abscissa="",
        title_bar="gene signal (z-score)",
        size_title_ordinate="eight", # ten
        size_title_abscissa="eight", # ten
        size_label_ordinate="eleven", # multi-panel: ten; individual: twelve
        size_label_abscissa="eleven", # multi-panel: ten; individual: twelve
        size_title_bar="twelve", # twelve
        size_label_bar="thirteen", # thirteen for whole; five for bar itself
        aspect="square", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )

    # Return information.
    return figure


def manage_plot_charts(
    table_box=None,
    table_heatmap_individual_1=None,
    table_heatmap_individual_2=None,
    table_heatmap_mean=None,
    table_allocation_1=None,
    table_allocation_2=None,
    box_features=None,
    index_features=None,
    index_observations=None,
    box=None,
    heatmap_individual=None,
    heatmap_mean=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        table_box (object): Pandas data-frame table of values of signal
            intensity for features across columns and sample observations
            in groups across rows
        table_heatmap_individual_1 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_individual_2 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_mean (object): Pandas data-frame table of descriptive
            statistics for values of signal intensity for features across
            rows and groups of sample observations across columns
        table_allocation_1 (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table
            'table_heatmap_individual_1' and the sequence of genes across rows
            in table 'table_heatmap_mean'
        table_allocation_2 (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table
            'table_heatmap_individual_2'
        box_features (list<str>): identifiers of features for which to create
            box charts
        index_features (str): name for index corresponding to features
        index_observations (str): name for index corresponding to observations
        box (bool): whether to create box charts
        heatmap_individual (bool): whether to create heatmap chart for
            individual values of signal intensity
        heatmap_mean (bool): whether to create heatmap chart for means of
            signal intensity
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of figure objects from MatPlotLib

    """

    # Copy information in table.
    table_box = table_box.copy(deep=True)
    table_heatmap_individual_1 = table_heatmap_individual_1.copy(deep=True)
    table_heatmap_individual_2 = table_heatmap_individual_2.copy(deep=True)
    table_heatmap_mean = table_heatmap_mean.copy(deep=True)
    table_allocation_1 = table_allocation_1.copy(deep=True)
    table_allocation_2 = table_allocation_2.copy(deep=True)
    # Copy other information.
    box_features = copy.deepcopy(box_features)

    ##########
    # Box
    if (
        box and (len(box_features) < 50)
    ):
        figures_box = list()
        for feature in box_features:
            record_box = dict()
            record_box["feature"] = feature
            record_box["name"] = str("box_" + feature)
            record_box["figure"] = create_plot_chart_box(
                table=table_box,
                column_feature=feature,
                column_group="group",
                report=False,
            )
            figures_box.append(record_box)
            pass
    else:
        figures_box = None
        pass


    ##########
    # Heatmap Individual
    if heatmap_individual:
        figure_heatmap_individual_1 = create_plot_chart_heatmap_individual(
            table_signal=table_heatmap_individual_1,
            table_feature_sets=table_allocation_1,
            index_columns=index_features,
            index_rows=index_observations,
            column_group="group",
            report=report,
        )
        figure_heatmap_individual_2 = create_plot_chart_heatmap_individual(
            table_signal=table_heatmap_individual_2,
            table_feature_sets=table_allocation_2,
            index_columns=index_features,
            index_rows=index_observations,
            column_group="group",
            report=report,
        )
    else:
        figure_heatmap_individual_1 = None
        figure_heatmap_individual_2 = None
        pass


    ##########
    # Heatmap Mean
    if heatmap_mean:
        figure_heatmap_mean = create_plot_chart_heatmap_mean(
            table=table_heatmap_mean,
            index_columns="group_observations",
            index_rows="feature",
            report=False,
        )
    else:
        figure_heatmap_mean = None
        pass

    ##########
    # Collect information.
    pail = dict()
    pail["box"] = figures_box
    pail["heatmap_individual_1"] = figure_heatmap_individual_1
    pail["heatmap_individual_2"] = figure_heatmap_individual_2
    pail["heatmap_mean"] = figure_heatmap_mean

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.description.py")
        function = str(
            "manage_plot_charts" +
            "()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)

    # Return information.
    return pail


###############################################################################
# Procedure


##########
# Control procedure within branch for iteration.


def control_procedure_part_branch(
    name_group_instances=None,
    instances_parameter=None,
    index_genes=None,
    index_samples=None,
    paths=None,
    report=None,
):
    """
    Control branch of procedure.

    table_sample (object): Pandas data-frame table of information about
        samples
    table_gene (object): Pandas data-frame table of information about genes
    table_signal (object): Pandas data-frame table of values of signal
        intensity corresponding to genes across rows and samples across
        columns

    arguments:
        name_group_instances (str): name for a group of instances in the
            collection of parameters
        instances_parameter (list<dict>): collection of instances of parameters
        index_genes (str): name for index corresponding to genes, which
            is a column in the original source table of signals
        index_samples (str): name for index corresponding to samples,
            which is not a column in the original source table of signals but
            will be in a novel product table
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Organize information.
    # Copy other information.
    instances_parameter = copy.deepcopy(instances_parameter)

    ##########
    # Extract parameters for instances in current group.
    # Filter instances by group.
    instances_group = list(filter(
        lambda record: (str(record["name_group"]) == name_group_instances),
        instances_parameter
    ))

    ##########
    # Read and organize source information from file.

    pail_source = read_source(
        tissue=instances_group[0]["tissue"],
        paths=paths,
        report=report,
    )
    # Organize table of information about sample observations.
    table_sample = pail_source["table_sample"]
    table_sample["inclusion"] = table_sample["inclusion"].astype("str")
    # Organize table of information about signals.
    table_signal = pail_source["table_signal"]
    # Organize indices in table.
    table_signal.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Organize table of information about genes.
    table_gene = pail_source["table_gene"]
    # Organize indices in table.
    table_gene.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    ##########
    # Prepare information about genes.
    # Notice that the parameters only use sets of genes from the first instance
    # in each group of instances.
    # Determine total set of genes with available signals.
    genes_signal = copy.deepcopy(
        table_signal[index_genes].unique().tolist()
    )
    # Prepare information about genes in sets.
    pail_gene_sets = prepare_sets_gene(
        names_sets_gene_inclusion=(
            instances_group[0]["names_sets_gene_inclusion"]
        ),
        names_sets_gene_cluster=instances_group[0]["names_sets_gene_cluster"],
        names_sets_gene_allocation=(
            instances_group[0]["names_sets_gene_allocation"]
        ),
        genes_available=genes_signal,
        paths=paths,
        report=report,
    )
    # Prepare table of allocation of total inclusion genes to sets.
    genes_inclusion = copy.deepcopy(
        pail_gene_sets["sets_gene_inclusion"]["main"]
    )
    columns_sequence = copy.deepcopy(
        instances_group[0]["sequence_sets_gene_allocation"]
    )
    columns_sequence.insert(0, index_genes)
    #columns_sequence.append("other")
    table_gene_sets_allocation = determine_gene_sets_allocation(
        genes_inclusion=genes_inclusion,
        sets_gene_allocation=pail_gene_sets["sets_gene_allocation"],
        index_genes=index_genes,
        column_other=None, # None, "", or "other"
        columns_sequence=columns_sequence,
        report=report,
    )
    # Prepare translations for genes
    table_gene_selection = table_gene.loc[(
        table_gene["gene_identifier_base"].isin(genes_signal)
    ), :].copy(deep=True)
    # Extract information for translation of names of columns.
    table_translations = table_gene_selection.filter(
        items=["gene_identifier_base", "gene_name",],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_translations["gene_name"].to_list(),
        index=table_translations["gene_identifier_base"],
    )
    translations_gene = series_translations.to_dict()
    # Translate identifiers of genes in set.
    genes_inclusion_translation = copy.deepcopy(genes_inclusion)
    if (translations_gene is not None):
        genes_inclusion_translation = list(map(
            lambda feature: (translations_gene[feature]),
            genes_inclusion_translation
        ))
        pass

    ##########
    # Prepare information about samples in groups.
    # Collect information.
    names_groups_samples_sequence = list()
    groups_samples = dict()
    samples_inclusion = list()
    # Iterate on instances of parameters in current group.
    for instance in instances_group:
        # Filter and extract identifiers of cohort sample observations
        # corresponding to selection criteria for current instance.
        samples = (
            porg.filter_extract_table_row_identifiers_by_columns_categories(
                table=table_sample,
                column_identifier="identifier_signal",
                name=instance["instance"], # or "name_instance"
                columns_categories=instance["selection_samples_primary"],
                report=report,
        ))
        # Collect information.
        names_groups_samples_sequence.append(instance["instance"])
        groups_samples[instance["instance"]] = samples
        samples_inclusion.extend(samples)
        pass
    # Collect unique names of sample observations.
    samples_inclusion = putly.collect_unique_elements(
        elements=samples_inclusion,
    )

    ##########
    # Prepare basic tables.
    # Prepare tables for signals.
    pail_tables = (
        pdesc.describe_signals_for_features_sets_in_observations_groups(
            table=table_signal,
            index_features=index_genes,
            index_observations=index_samples, # assigned in new tables
            features_inclusion=genes_inclusion,
            observations_inclusion=samples_inclusion,
            groups_features=pail_gene_sets["sets_gene_cluster"],
            groups_observations=groups_samples,
            names_groups_features_sequence=(
                instances_group[0]["sequence_sets_gene_cluster"]
            ),
            names_groups_observations_sequence=(
                names_groups_samples_sequence
            ),
            translations_features=translations_gene,
            translations_observations=None,
            report=report,
    ))

    # Prepare table for gene set allocation.
    # Translate names of features and observations.
    table_allocation_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_gene_sets_allocation,
            index_rows=index_genes,
            translations_columns=None,
            translations_rows=translations_gene,
            remove_redundancy=True,
            report=report,
    ))

    # Sort sequence of genes in tables of allocation to sets to match the
    # sequences in tables of signals.
    table_allocation_3 = sort_table_rows_gene_sets_allocation(
        table_gene_sets=table_allocation_translation,
        table_signal=pail_tables["table_3"],
        index_genes=index_genes,
        indices_rows_signal=[index_samples, "group",],
        report=report,
    )
    table_allocation_4 = sort_table_rows_gene_sets_allocation(
        table_gene_sets=table_allocation_translation,
        table_signal=pail_tables["table_4"],
        index_genes=index_genes,
        indices_rows_signal=[index_samples, "group",],
        report=report,
    )

    ##########
    # Prepare plot charts.
    pail_plot = manage_plot_charts(
        table_box=pail_tables["table_2"],
        table_heatmap_individual_1=pail_tables["table_3"],
        table_heatmap_individual_2=pail_tables["table_4"],
        table_heatmap_mean=pail_tables["table_7"],
        table_allocation_1=table_allocation_3,
        table_allocation_2=table_allocation_4,
        box_features=genes_inclusion_translation,
        index_features=index_genes,
        index_observations=index_samples,
        box=True,
        heatmap_individual=True,
        heatmap_mean=True,
        report=report,
    )

    ##########
    # Collect information.
    pail_write_plot = dict()
    if (pail_plot["box"] is not None) and (len(pail_plot["box"]) > 0):
        for record_box in pail_plot["box"]:
            pail_write_plot[record_box["name"]] = record_box["figure"]
            pass
        pass
    pail_write_plot["heatmap_individual_1"] = pail_plot["heatmap_individual_1"]
    pail_write_plot["heatmap_individual_2"] = pail_plot["heatmap_individual_2"]
    pail_write_plot["heatmap_mean"] = pail_plot["heatmap_mean"]

    ##########
    # _. Write product information to file.
    #paths["out_data"]
    #paths["out_plot"]

    # Define paths to directories.
    path_directory_plot_group = os.path.join(
        paths["out_procedure_plot"], str(name_group_instances),
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_plot_group,
    )
    # Write figures to file.
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot,
        format="jpg", # jpg, png, svg
        resolution=300,
        path_directory=path_directory_plot_group,
    )
    pass


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
    project="exercise"
    routine="transcriptomics"
    procedure="compare_sets_groups"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.compare_sets_groups.py")
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
        paths["in_sets_gene"],
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
    # 2.2. Read and organize information about parameters for instances.
    instances = read_organize_source_parameter_instances(
        paths=paths,
        report=report,
    )
    # Organize information for groups of instances.
    groups_instances = list()
    for instance in instances:
        groups_instances.append(instance["name_group"])
    # Collect unique names of features.
    groups_instances_unique = putly.collect_unique_elements(
        elements=groups_instances,
    )

    ##########
    # Execute procedure iteratively across groups of instances of parameters.
    if True:
        for group_instances in groups_instances_unique:
            control_procedure_part_branch(
                name_group_instances=group_instances,
                instances_parameter=instances,
                index_genes="identifier_gene",
                index_samples="identifier_sample",
                paths=paths,
                report=report,
            )
            pass
    else:
        control_procedure_part_branch(
            name_group_instances="2_muscle_exercise_interaction_age_exercise_genes_interaction",
            instances_parameter=instances,
            index_genes="identifier_gene",
            index_samples="identifier_sample",
            paths=paths,
            report=report,
        )
    pass


###############################################################################
# End
