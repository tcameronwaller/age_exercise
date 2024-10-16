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
import partner.plot_cluster as pcplot
import partner.parallelization as prall
import exercise.transcriptomics.organize_signal as exrosig

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.


def initialize_directories(
    project=None,
    procedure=None,
    path_directory_dock=None,
    restore=None,
    report=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
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
        paths["dock"], "in_data",
    )
    paths["in_demonstration"] = os.path.join(
        paths["dock"], "in_demonstration",
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters",
    )
    paths["in_parameters_private"] = os.path.join(
        paths["dock"], "in_parameters_private", str(project),
    )
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    paths["out_procedure"] = os.path.join(
        paths["out_project"], str(procedure),
    )
    paths["out_data"] = os.path.join(
        paths["out_procedure"], "data",
    )
    paths["out_plot"] = os.path.join(
        paths["out_procedure"], "plot",
    )

    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"],
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
        print("module: exercise.transcriptomics.interaction.py")
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
    paths["in_demonstration_partner"] = os.path.join(
        paths["in_demonstration"], "partner",
    )
    paths["in_transcriptomics_signal"] = os.path.join(
        paths["out_project"], "transcriptomics", "organize_signal", "whole",
        "preparation",
    )

    # Define paths to child files.
    path_file_table_demonstration = os.path.join(
        paths["in_demonstration_partner"],
        "table_plot_heatmap_label_features_groups_cluster_observations.tsv",
    )
    path_file_table_sample = os.path.join(
        paths["out_project"], "transcriptomics", "organize_sample", "data",
        "table_sample.pickle",
    )
    path_file_table_gene = os.path.join(
        paths["out_project"], "transcriptomics", "organize_signal", "whole",
        "preparation", "table_gene_adipose.pickle",
    )
    path_file_table_signal = os.path.join(
        paths["out_project"], "transcriptomics", "organize_signal", "whole",
        "preparation", "table_signal_scale_adipose.pickle",
    )

    # Collect information.
    pail = dict()
    # Read information from file.
    types_columns = define_column_types_table_demonstration()
    pail["table_demonstration"] = pandas.read_csv(
        path_file_table_demonstration,
        sep="\t",
        header=0,
        dtype=types_columns,
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
        print("module: exercise.transcriptomics.scratch.py")
        print("function: read_source()")
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
    types_columns["group"] = "string"
    types_columns["instance"] = "string"
    types_columns["cohort_selection_primary"] = "string"
    types_columns["cohort_selection_secondary"] = "string"
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
        paths["in_parameters_private"], "transcriptomics",
        "table_cohort_group_comparison.tsv",
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
            pail["group"] = str(row["group"])
            pail["tissue"] = str(row["tissue"])
            pail["instance"] = str(row["instance"])
            pail["name_instance"] = "_".join([
                str(row["tissue"]),
                str(row["sort"]),
                str(row["instance"])
            ])
            pail["cohort_selection_primary"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["cohort_selection_primary"],
                )
            )
            pail["cohort_selection_secondary"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["cohort_selection_secondary"],
                )
            )
            # Extract names of columns corresponding to feature variables for
            # which to calculate tertiles.
            #columns_tertile = extract_source_columns_for_tertiles(
            #    cohort_selection=pail["cohort_selection_secondary"],
            #    report=report,
            #)
            # Collect names of unique features or columns relevant to current
            # instance.
            features = list()
            dictionaries = [
                "cohort_selection_primary",
                "cohort_selection_secondary",
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
        print("module: exercise.transcriptomics.comparison.py")
        print("function: read_organize_source_parameter_instances()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        print("instances:")
        print(instances)
        print("instance[0]:")
        print(instances[0])
        pass
    # Return information.
    return instances



###############################################################################
# Procedure


##########
# Execute main procedure.

# TODO: TCW; 16 October 2024
# NEXT!!!
# 1. Organize new high-level driver function
# 2. organize a few new mid-level functions to organize necessary operations
#    - reading in parameters for groups of observations
#    - determining gene translations
#    - etc
#3. write tables to file
#_. write plots to file

# 2. new low-level driver function
#    use information from the tables above to prepare plot charts
#    1. box plots
#    2. heat map for the summary information (table 6, I think)
#    3. heat map with group bar for clustered individual signals (clustered table 3)
#  arguments:
#    bar_chart: True/False
#    heatmap_mean: True/False
#    heatmap_median: True/False
#    heatmap_groups_individuals: True/False
#    report: True/False




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
    procedure="scratch"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.scratch.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # 1. Initialize directories for read of source and write of product files.
    paths = initialize_directories(
        project=project,
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



    ################################
    # Parameters for future function...
    index_genes = "identifier_gene"
    identifiers_genes = [
        "ENSG00000101405",
        "ENSG00000100344",
        "ENSG00000146674",
        "ENSG00000164530",
        "ENSG00000183785",
        "ENSG00000196169",
        "ENSG00000183671",
        "ENSG00000176387",
        "ENSG00000196208",
        "ENSG00000099337",
        "ENSG00000114737",
    ]
    #collections_groups_observations[0] <-- argument to new function... determine separately
    # translations_gene <-- argument to new function... determine separately

    ############################
    # Prepare translations for genes
    # Ensure that all genes are in the table of signals.
    genes_all = copy.deepcopy(
        table_signal[index_genes].unique().tolist()
    )
    identifiers_genes_available = list(filter(
        lambda gene: (gene in genes_all),
        identifiers_genes
    ))
    # Prepare information about genes.
    table_gene_selection = table_gene.loc[(
        table_gene["gene_identifier_base"].isin(identifiers_genes_available)
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


    # Read and organize information about parameters for instances.
    instances = read_organize_source_parameter_instances(
        paths=paths,
        report=report,
    )
    collections = list()
    for instance in instances:
        collections.append(instance["group"])
    # Collect unique names of features.
    collections_unique = putly.collect_unique_elements(
        elements=collections,
    )
    # Collect information.
    collections_groups_observations = list()
    # Iterate on groups of cohort instances.
    for collection in collections_unique:
        # Collect information.
        groups_observations = dict()
        # Filter instances by group.
        instances_collection = list(filter(
            lambda record: (str(record["group"]) == collection),
            instances
        ))
        # Iterate on cohort instances in current group.
        for instance in instances_collection:
            # Filter and extract identifiers of cohort sample observations
            # corresponding to selection criteria for current instance.
            observations = (
                porg.filter_extract_table_row_identifiers_by_columns_categories(
                    table=table_sample,
                    column_identifier="identifier_signal",
                    name=instance["instance"], # or "name_instance"
                    columns_categories=instance["cohort_selection_primary"],
                    report=report,
            ))
            # Collect information.
            groups_observations[instance["instance"]] = observations
            pass
        # Collect information.
        collections_groups_observations.append(groups_observations)
        pass

    ##########
    # 3. Do stuff.


    ##################
    # Prepare basic tables.
    pail = pdesc.extract_describe_signals_for_features_in_observations_groups(
        table=pail_source["table_signal"],
        index_features="identifier_gene",
        index_observations="identifier_sample", # assigned in new tables
        features=identifiers_genes_available,
        groups_observations=collections_groups_observations[0],
        translations_features=translations_gene,
        translations_observations=None,
        report=report,
    )
    #pail["table_3"] <-- clustered heatmap of individual signals
    #pail["table_6"] <-- simple heatmap of means

    ###################
    # Cluster table for heatmap of individual signals

    # Cluster rows in table within groups.
    table_3_cluster = porg.cluster_table_rows_by_group(
        table=pail["table_3"],
        index_rows="identifier_sample",
        column_group="group",
    )
    # Organize indices in table.
    table_3_cluster.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_3_cluster.set_index(
        ["identifier_sample", "group"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster columns in table.
    table_3_cluster = porg.cluster_table_columns(
        table=table_3_cluster,
    )
    table_3_cluster.index = pandas.MultiIndex.from_tuples(
        table_3_cluster.index,
        names=["identifier_sample", "group"]
    )
    # Organize indices in table.
    table_3_cluster.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    print("!!!!!!!!!!!!!!!!!!!!!!!!! table after cluster !!!!!!!!!!!!!!!!!")
    print(table_3_cluster)


    ##########
    # _. Write product information to file.
    #paths["out_data"]
    #paths["out_plot"]


    pass


###############################################################################
# End
