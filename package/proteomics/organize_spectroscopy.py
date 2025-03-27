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
    paths["out_procedure_lists"] = os.path.join(
        paths["out_procedure"], "lists",
    )
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
        path_file_table_batch_b,
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

# eventually define a parameter table for each batch of proteomics data in a table
# parameters:
# translations for names of columns (similar to "table_subject_sample_feature_organization")
# indications of columns in relevant subsets
# example: names of samples in each batch (need a special column to indicate the batch)
# example: names of columns corresponding to attributes of the analytes themselves (identifiers, names, etc)



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
# 1. separate columns for analyte properties to separate table from columns for analyte signals in samples
# 2. filter to remove analytes with inadequate signal (missing or not within threshold range) coverage across samples (start with 100%)
# 3.



def define_names_analyte_attributes():
    """
    Defines parameters.

    arguments:

    raises:

    returns:
        (dict): collection of information

    """

    # Collect information.
    pail = dict()

    # Main identifier of analytes.
    pail["identifier_analyte"] = [
        "identifier_protein_uniprot",
    ]

    # Attributes of analytes.
    pail["attributes_analyte"] = [
        "identifier_protein_uniprot",
        "identifier_gene_symbol",
        "identifier_gene_entrez",
        "identifier_gene_ensembl",
    ]

    # Return information.
    return pail


def define_names_samples_by_batch():
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
    pail["samples_batch_a"] = [
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
    pail["samples_batch_b"] = [
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


def separate_organize_table_measurement_analyte_signal(
    table_measurement=None,
    name_identifier_analyte=None,
    names_analyte=None,
    names_signal=None,
    explicate_indices=None,
    report=None,
):
    """
    Separates information from table for measurements.

    ----------
    Format of table for measurements from instrument and preprocessing (name:
    "table_measurement")
    Format of source table is in wide format with features across rows and
    observations across columns. Multiple special columns give identifiers and
    other information about attributes of analytes. Another special column
    gives signals from measurement of each analyte in a special bridge sample.
    For versatility, this table does not have explicitly defined indices
    across columns or rows. Values for observations of features are on a
    quantitative, continuous, interval or ratio scale of measurement.
    ----------
    samples     attribute_1 attribute_2 bridge sample_1 sample_2 sample_3 ...
    analyte
    analyte_1   stuff       stuff       0.001  0.001    0.001    0.001    ...
    analyte_2   stuff       stuff       0.001  0.001    0.001    0.001    ...
    analyte_3   stuff       stuff       0.001  0.001    0.001    0.001    ...
    analyte_4   stuff       stuff       0.001  0.001    0.001    0.001    ...
    analyte_5   stuff       stuff       0.001  0.001    0.001    0.001    ...
    ----------

    ----------
    Format of table for attributes of analytes (name: "table_analyte")
    Depending on parameters to this function, this table optionally has
    explicitly defined indices across columns and rows.
    ----------
    attributes   attribute_1 attribute_2 attribute_3 attribute_4 ...
    analyte
    analyte_1    stuff       stuff       stuff       stuff       ...
    analyte_2    stuff       stuff       stuff       stuff       ...
    analyte_3    stuff       stuff       stuff       stuff       ...
    analyte_4    stuff       stuff       stuff       stuff       ...
    analyte_5    stuff       stuff       stuff       stuff       ...
    ----------

    ----------
    Format of table for signals of analytes across samples (name:
    "table_signal")
    The names of columns in the table for signals indicate the sample
    corresponding to measurements of the analyte; however, use terminology of
    "identifier_signal" to distinguish these identifiers, which might differ
    from those for samples in the table for attributes of samples. Depending on
    parameters to this function, this table optionally has explicit definitions
    of indices across columns and rows.
    ----------
    signal      bridge sample_1  sample_2  sample_3  sample_4  sample_5 ...
    analyte
    analyte_1   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_2   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_3   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_4   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_5   0.001  0.001     0.001     0.001     0.001     0.001    ...
    ----------

    arguments:
        table_measurement (object): Pandas data-frame table of attributes of
            analytes and signals from their measurements across samples
        name_identifier_analyte (str): name of column in table corresponding to
            identifiers of analytes, which must also be in the 'names_analyte'
            list
        names_analyte (list<str>): names of columns in table corresponding to
            attributes of analytes
        names_signal (list<str>): names of columns in table corresponding to
            signals from measurements of analytes across samples
        explicate_indices (bool): whether to explicate, define, or specify
            explicit indices across columns and rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table_measurement = table_measurement.copy(deep=True)
    names_analyte = copy.deepcopy(names_analyte)
    names_signal = copy.deepcopy(names_signal)

    # Organize identifiers for analytes.
    table_measurement["identifier_analyte"] = table_measurement[
        name_identifier_analyte
    ]

    # Separate columns into separate tables.
    #names_analyte.remove("identifier_analyte")
    names_signal.insert(0, "identifier_analyte")
    table_analyte = table_measurement.loc[
        :, table_measurement.columns.isin(names_analyte)
    ].copy(deep=True)
    table_signal = table_measurement.loc[
        :, table_measurement.columns.isin(names_signal)
    ].copy(deep=True)

    # Organize identifiers for analytes.
    table_analyte["identifier_analyte"] = table_analyte[
        name_identifier_analyte
    ]

    # Organize columns in tables.
    names_analyte.insert(0, "identifier_analyte")
    table_analyte = porg.filter_sort_table_columns(
        table=table_analyte,
        columns_sequence=names_analyte,
        report=report,
    )
    table_signal = porg.filter_sort_table_columns(
        table=table_signal,
        columns_sequence=names_signal,
        report=report,
    )

    # Translate names of columns.
    #translations = dict()
    #translations["identifier_analyte"] = "identifier"
    #table_analyte.rename(
    #    columns=translations,
    #    inplace=True,
    #)

    # Organize indices in table.
    table_analyte = porg.explicate_table_indices_columns_rows_single_level(
        table=table_analyte,
        index_columns="attributes",
        index_rows="identifier_analyte",
        explicate_indices=explicate_indices,
        report=False,
    )
    table_signal = porg.explicate_table_indices_columns_rows_single_level(
        table=table_signal,
        index_columns="identifier_signal",
        index_rows="identifier_analyte",
        explicate_indices=explicate_indices,
        report=False,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.proteomics.organize_spectroscopy.py")
        name_function = str(
            "separate_organize_table_measurement_analyte_signal()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("source table of measurements:")
        count_rows = (table_measurement.shape[0])
        count_columns = (table_measurement.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_measurement.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("table of information about analytes:")
        count_rows = (table_analyte.shape[0])
        count_columns = (table_analyte.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_analyte.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("table of information about signals:")
        count_rows = (table_signal.shape[0])
        count_columns = (table_signal.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_signal.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Collect information.
    pail = dict()
    pail["table_analyte"] = table_analyte
    pail["table_signal"] = table_signal
    # Return information.
    return pail


def filter_fill_intensity_table_signal(
    table_signal=None,
    name_index_columns=None,
    name_index_rows=None,
    name_bridge=None,
    names_samples=None,
    explicate_indices=None,
    report=None,
):
    """
    Filters analytes by proportion of signal intensities that are nonmissing
    and valid, then fills missing values by imputation across samples.

    ----------
    Format of table for signals of analytes across samples (name:
    "table_signal")
    The names of columns in the table for signals indicate the sample
    corresponding to measurements of the analyte; however, use terminology of
    "identifier_signal" to distinguish these identifiers, which might differ
    from those for samples in the table for attributes of samples. For
    versatility, this table does not have explicit defininitions of indices
    across columns or rows.
    ----------
    signal      bridge sample_1  sample_2  sample_3  sample_4  sample_5 ...
    analyte
    analyte_1   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_2   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_3   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_4   0.001  0.001     0.001     0.001     0.001     0.001    ...
    analyte_5   0.001  0.001     0.001     0.001     0.001     0.001    ...
    ----------

    arguments:
        table_signal (object): Pandas data-frame table of signals for analytes
            from their measurements across samples
        name_index_columns (str): name of index across columns in table
        name_index_rows (str): name of index across rows in table
        name_bridge (str): name of column in table corresponding to signals
            from measurements of analytes in the pooled bridge sample
        names_samples (list<str>): names of columns in table corresponding to
            signals from measurements of analytes across samples
        explicate_indices (bool): whether to explicate, define, or specify
            explicit indices across columns and rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table_signal = table_signal.copy(deep=True)
    table_filter = table_signal.copy(deep=True)
    names_samples = copy.deepcopy(names_samples)
    names_samples_all = copy.deepcopy(names_samples)

    # Extract names or identifiers of rows for analytes.
    #names_analytes = copy.deepcopy(list(table_filter.index.values))
    names_analytes = copy.deepcopy(
        table_filter[name_index_rows].unique().tolist()
    )

    # Organize names of columns for signals of analytes in samples.
    names_samples_all.insert(0, name_bridge)

    # Replace values of zero for signal intensity with missing values.
    # Only replace values within table's columns for samples.
    # This implementation is more concise than iteration across specific
    # columns.
    table_filter[names_samples_all] = table_filter[names_samples_all].replace(
        to_replace=0,
        value=pandas.NA,
    )
    # Replace values less than zero with missing values.
    table_filter[names_samples_all][table_filter[names_samples_all] < 0] = (
        pandas.NA
    )

    # Filter table's rows corresponding to analytes by their proportions of
    # nonmissing and valid signal intensities in measurements across the bridge
    # sample or samples.
    table_filter.dropna(
        axis="index",
        how="all",
        subset=[name_bridge],
        inplace=True,
    )
    names_analytes = copy.deepcopy(
        table_filter[name_index_rows].unique().tolist()
    )
    table_filter = (
        porg.filter_table_rows_by_proportion_nonmissing_threshold(
            table=table_filter,
            index_columns=name_index_columns,
            index_rows=name_index_rows,
            columns_selection=[name_bridge],
            rows_selection=names_analytes,
            threshold_low=0.0,
            threshold_high=None,
            proportion=1.0,
            report=report,
    ))

    # Filter table's rows corresponding to analytes by their proportions of
    # nonmissing and valid signal intensities in measurements across samples.
    table_filter.dropna(
        axis="index",
        how="all",
        subset=names_samples,
        inplace=True,
    )
    names_analytes = copy.deepcopy(
        table_filter[name_index_rows].unique().tolist()
    )
    table_filter = (
        porg.filter_table_rows_by_proportion_nonmissing_threshold(
            table=table_filter,
            index_columns=name_index_columns,
            index_rows=name_index_rows,
            columns_selection=names_samples,
            rows_selection=names_analytes,
            threshold_low=0.0,
            threshold_high=None,
            proportion=0.9,
            report=report,
    ))

    # Copy information.
    table_fill = table_filter.copy(deep=True)

    # Fill missing values of signal intensity within table's rows corresponding
    # to analytes with measurements across samples.
    # Notice that this fill procedure occurs after the filters on proportions
    # of missing signal intensities across features and observations. Hence,
    # those previous filters regulate the extent of fills on missing values.
    # fill missing values for each analyte (across rows)
    names_analytes = copy.deepcopy(
        table_fill[name_index_rows].unique().tolist()
    )
    table_fill = porg.fill_missing_values_table_by_row(
        table=table_fill,
        index_columns=name_index_columns,
        index_rows=name_index_rows,
        columns_selection=names_samples, # samples???
        rows_selection=names_analytes,
        method="half_minimum",
        report=report,
    )

    # Organize indices in table.
    table_fill = porg.explicate_table_indices_columns_rows_single_level(
        table=table_fill,
        index_columns=name_index_columns,
        index_rows=name_index_rows,
        explicate_indices=explicate_indices,
        report=False,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: age_exercise")
        print("subpackage: proteomics")
        print("module: organize_spectroscopy.py")
        name_function = str(
            "filter_fill_intensity_table_signal()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("source table of signals:")
        count_rows = (table_signal.shape[0])
        count_columns = (table_signal.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_signal)
        putly.print_terminal_partition(level=5)
        print("table of signals after filters:")
        count_rows = (table_filter.shape[0])
        count_columns = (table_filter.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_filter)
        putly.print_terminal_partition(level=5)
        print("table of signals after filters and fills:")
        count_rows = (table_fill.shape[0])
        count_columns = (table_fill.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_fill)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_fill





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

    # Copy information.
    table_batch_a = pail_source["table_batch_a"].copy(deep=True)
    table_batch_b = pail_source["table_batch_b"].copy(deep=True)

    # Determine names of columns in table for attributes of analytes, including
    # identifiers.
    pail_names_attributes = define_names_analyte_attributes()

    # Determine names of columns in table for measurements of analytes across
    # samples.
    pail_names_samples = define_names_samples_by_batch()
    #pail_names_samples["names_samples_batch_a"]
    #pail_names_samples["names_samples_batch_b"]

    # Filter analytes in the measurement table by their attributes.
    # TODO: TCW; 26 March 2025
    # 1. mainly filter by the confidence of identification by Proteome Discoverer: "High"

    # From measurement table, separate columns for analyte attributes from
    # columns for signals from measurements of analytes across samples.
    # Copy information.
    names_attributes = copy.deepcopy(pail_names_attributes["attributes_analyte"])
    names_signal_a = copy.deepcopy(pail_names_samples["samples_batch_a"])
    names_signal_b = copy.deepcopy(pail_names_samples["samples_batch_b"])
    names_signal_a.insert(0, "bridge")
    names_signal_b.insert(0, "bridge")
    # Separate and organize columns.
    pail_split_a = separate_organize_table_measurement_analyte_signal(
        table_measurement=table_batch_a,
        name_identifier_analyte=pail_names_attributes["identifier_analyte"],
        names_analyte=names_attributes,
        names_signal=names_signal_a,
        explicate_indices=False,
        report=report,
    )
    pail_split_b = separate_organize_table_measurement_analyte_signal(
        table_measurement=table_batch_b,
        name_identifier_analyte=pail_names_attributes["identifier_analyte"],
        names_analyte=names_attributes,
        names_signal=names_signal_b,
        explicate_indices=False,
        report=report,
    )

    #################
    # New Stuff in progress

    # In signal table, filter analytes to remove those with inadequate values
    # of signal that are nonmissing and withing threshold range.
    # In signal table, fill missing values by imputation.

    # Copy information.
    names_samples_a = copy.deepcopy(pail_names_samples["samples_batch_a"])
    names_samples_b = copy.deepcopy(pail_names_samples["samples_batch_b"])
    # Filter and fill signals.
    table_signal_a = filter_fill_intensity_table_signal(
        table_signal=pail_split_a["table_signal"],
        name_index_columns="identifier_signal",
        name_index_rows="identifier_analyte",
        name_bridge="bridge",
        names_samples=names_samples_a,
        explicate_indices=True,
        report=report,
    )
    table_signal_b = filter_fill_intensity_table_signal(
        table_signal=pail_split_b["table_signal"],
        name_index_columns="identifier_signal",
        name_index_rows="identifier_analyte",
        name_bridge="bridge",
        names_samples=names_samples_b,
        explicate_indices=True,
        report=report,
    )

    # Correct for artifactual batch effects.


    # Scale values of signal intensity for measurements of analytes across
    # samples.
    # The goal of this scaling is to make the individual samples more
    # comparable to each other.
    # This scaling can decrease the variance or noise in measurements between
    # samples.


    pass


###############################################################################
# End
