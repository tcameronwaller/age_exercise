"""
Supply functionality for process and analysis of data from proteomics using
mass spectroscopy.

This module 'proteomics_mass_spectroscopy' is part of the 'exercise' package.

This module is not directly executable.

This subpackage 'exercise' provides executable functionality under the
management of a higher level package. Importation paths require this hierarchy.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This initialization file is part of the project package directory 'exercise'
    (https://github.com/tcameronwaller/exercise/).

    Project 'exercise' supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

    The code within project 'exercise' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the GNU
    General Public License, or (at your option) any later version.

    The code within project 'exercise' is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'exercise'. If not, see <http://www.gnu.org/licenses/>.
"""


################################################################################
# Notes


################################################################################
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

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.description as pdesc
#import partner.regression as preg
#import partner.plot as pplot
import partner.parallelization as prall

################################################################################
# Functionality


##########
# Initialization


def initialize_directories(
    project=None,
    set=None,
    path_directory_dock=None,
    restore=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project for parent directory
        set (str): name of set of process and analysis for child directory
        path_directory_dock (str): path to dock directory for source and product
            directories and files
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
        paths["dock"], "in_data", "scratch", str(project), str(set),
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters", "scratch", str(project), str(set),
    )
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    paths["out_set"] = os.path.join(
        paths["out_project"], str(set),
    )
    paths["out_test"] = os.path.join(
        paths["out_set"], "test",
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        putly.remove_directory(path=paths["out_project"])
        putly.remove_directory(path=paths["out_set"])
        putly.remove_directory(path=paths["out_test"])
    # Initialize directories.
    putly.create_directories(
        path=paths["out_project"]
    )
    putly.create_directories(
        path=paths["out_set"]
    )
    putly.create_directories(
        path=paths["out_test"]
    )
    # Return information.
    return paths


##########
# Read from source


def define_table_column_types_main():
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
    types_columns["Protein FDR Confidence: Combined"] = "string"
    types_columns["Accession"] = "string"
    types_columns["Description"] = "string"
    types_columns["Intensity.1"] = "float32"
    types_columns["Intensity.2"] = "float32"
    types_columns["Intensity.3"] = "float32"
    types_columns["Intensity.4"] = "float32"
    types_columns["Intensity.5"] = "float32"
    types_columns["Intensity.6"] = "float32"
    types_columns["Intensity.7"] = "float32"
    types_columns["Intensity.8"] = "float32"
    types_columns["Intensity.9"] = "float32"
    types_columns["Intensity.10"] = "float32"
    types_columns["# Peptides"] = "int32"
    types_columns["# Unique Peptides"] = "int32"
    types_columns["Gene Symbol"] = "string"
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
    types_columns["sort"] = "int32"
    types_columns["group"] = "string"
    types_columns["sort_group"] = "int32"
    types_columns["identifier_original"] = "string"
    types_columns["identifier_novel"] = "string"
    types_columns["description"] = "string"
    types_columns["pair"] = "string"
    types_columns["abbreviation"] = "string"
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
    path_file_table_main = os.path.join(
        paths["in_data"], "table_intensity_mouse_zcr_2024-06-21.tsv",
    )
    path_file_table_sample = os.path.join(
        paths["in_parameters"], "table_samples_attributes_filter_sort.tsv",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of values of intensity across samples and proteins.
    types_columns = define_table_column_types_main()
    pail["table_main"] = pandas.read_csv(
        path_file_table_main,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",],
    )
    # Table of attributes of samples.
    types_columns = define_table_column_types_sample()
    pail["table_sample"] = pandas.read_csv(
        path_file_table_sample,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",],
    )

    # Return information.
    return pail


##########
# Organization


def organize_table_main(
    table_sample=None,
    table_main=None,
    report=None,
):
    """
    Organizes information in table.

    arguments:
        table_sample (object): Pandas data-frame table of information about
            samples represented in the main table
        table_main (object): Pandas data-frame table of values of intensity
            across samples in columns and proteins in rows
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)
    table_main = table_main.copy(deep=True)

    # Organize information in table.

    # Determine original and novel names for columns representing samples.
    # Filter information about samples.
    table_sample_inclusion = table_sample.loc[
        (table_sample["inclusion"] == 1), :
    ]
    samples = copy.deepcopy(
        table_sample_inclusion["identifier_novel"].to_list()
    )

    # remove all columns except for "identifier_original" and "identifier_novel"
    # Filter and sort table's columns.
    table_sample_translation = porg.filter_sort_table_columns(
        table=table_sample_inclusion,
        columns_sequence=[
            "identifier_original",
            "identifier_novel",
        ],
        report=False,
    )
    # Organize information.
    table_sample_translation.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_sample_translation.set_index(
        ["identifier_original"],
        append=False,
        drop=True,
        inplace=True,
    )
    translations_sample = table_sample_translation["identifier_novel"].to_dict()
    # Translate names of columns.
    translations = dict()
    translations["Protein FDR Confidence: Combined"] = "confidence"
    translations["Accession"] = "identifier_protein_uniprot"
    translations["Description"] = "description"
    translations["# Peptides"] = "coverage_peptides"
    translations["# Unique Peptides"] = "coverage_peptides_unique"
    translations["Gene Symbol"] = "name_gene"
    translations.update(translations_sample)
    table_main.rename(
        columns=translations,
        inplace=True,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        columns = copy.deepcopy(table_main.columns.to_list())
        print("Names of columns: ")
        print(columns)
        putly.print_terminal_partition(level=4)
        print("Table: ")
        print(table_main)
        putly.print_terminal_partition(level=4)
        pass

    # Collect information.
    pail = dict()
    pail["samples"] = samples
    pail["table_main"] = table_main
    # Return information.
    return pail


##########
# Filter


def define_column_sequence_table_protein():
    """
    Defines the columns in sequence within table.

    arguments:

    raises:

    returns:
        (list<str>): variable types of columns within table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "identifier_protein_uniprot",
        "name_gene",
        "description",
        "confidence",
        "coverage_peptides",
        "coverage_peptides_unique",
    ]
    # Return information.
    return columns_sequence

def define_row_test_filter():
    """
    Defines a row to test filter operations.

    values: 0.4526, 0.3325, 0.2953, 0.3015, 0.4379, 0.3127, 0.1578, 0.2351
    mean: 0.3157
    standard deviation: 0.0971
    (mean - (3 * standard deviation)): 0.0243
    minimum: 0.1578
    (minimum / 2): 0.0789

    arguments:

    raises:

    returns:
        (object): Pandas series of values of intensity across samples for a
            single protein

    """

    row = pandas.Series(
        data={
            "identifier_protein_uniprot": "KF0KFU",
            "name_gene": "MPC1",
            "description": "mitochondrial pyruvate carrier 1",
            "confidence": "High",
            "coverage_peptides": 2,
            "coverage_peptides_unique": 2,
            "control_1": 0.4526,
            "control_2": 0.3325,
            "control_3": float("nan"),
            "control_4": 0.2953,
            "control_5": 0.3015,
            "intervention_1": float("nan"),
            "intervention_2": 0.4379,
            "intervention_3": 0.3127,
            "intervention_4": 0.1578,
            "intervention_5": 0.2351,
        },
        index=[
            "identifier_protein_uniprot",
            "name_gene",
            "description",
            "confidence",
            "coverage_peptides",
            "coverage_peptides_unique",
            "control_1",
            "control_2",
            "control_3",
            "control_4",
            "control_5",
            "intervention_1",
            "intervention_2",
            "intervention_3",
            "intervention_4",
            "intervention_5",
        ],
    )
    # Return information.
    return row

def match_keep_table_row_identification(
    confidences_keep=None,
    peptides_keep=None,
    peptides_unique_keep=None,
    row_identifier=None,
    row_confidence=None,
    row_peptides=None,
    row_peptides_unique=None,
    report=None,
):
    """
    Determines whether to keep a row from a table.

    arguments:
        confidences_keep (list<str>): discrete values of confidence to keep
        peptides_keep (int): count of peptides to keep
        peptides_unique_keep (int): count of unique peptides to keep
        row_identifier (str): current row's identifier
        row_confidence (str): current row's indication of confidence
        row_peptides (int): current row's coverage count of peptides
        row_peptides_unique (int): current row's coverage count of unique
            peptides
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Determine whether to keep current row from table.
    # if (any(row_confidence in item for item in confidences_keep)):
    if (
        (pandas.notna(row_identifier)) and
        (len(str(row_identifier)) > 0) and
        (str(row_confidence).lower() in confidences_keep) and
        (row_peptides >= peptides_keep) and
        (row_peptides_unique >= peptides_unique_keep)
    ):
        indicator = 1
    else:
        indicator = 0
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Match function for filter by identification.")
        putly.print_terminal_partition(level=4)
        print("Indicator: " + str(indicator))
        putly.print_terminal_partition(level=4)
    # Return information.
    return indicator

def match_keep_table_row_quantification(
    row=None,
    columns_intensity=None,
    proportion_validity=None,
    report=None,
):
    """
    Determines whether to keep a row from a table.

    arguments:
        row (object): Pandas series of values of intensity across samples for
            a single protein
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        proportion_validity (float): proportion of values of intensity that must
            have valid, non-missing values to keep the row for each protein
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Copy information in row.
    row = row.copy(deep=True)
    # Extract information for current row from table.
    count_samples = len(columns_intensity)
    array_intensities = row[columns_intensity].to_numpy(
        dtype="float",
        na_value=numpy.nan,
        copy=True,
    )
    #array_intensities = array_intensities.astype(numpy.float32)
    count_validity = numpy.count_nonzero(~numpy.isnan(array_intensities))
    proportion_actual = float(count_validity / count_samples)
    # Determine whether to keep current row from table.
    if (
        (proportion_actual >= proportion_validity)
    ):
        indicator = 1
    else:
        indicator = 0
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Match function for filter by quantification.")
        putly.print_terminal_partition(level=4)
        print("Array of intensity values accross samples:")
        print(array_intensities.shape)
        print("Proportion Valid, Actual: " + str(proportion_actual))
        print("Proportion Valid, Threshold: " + str(proportion_validity))
        print("Indicator: " + str(indicator))
        putly.print_terminal_partition(level=4)
    # Return information.
    return indicator

def filter_table_main(
    table=None,
    columns_intensity=None,
    proportion_validity=None,
    report=None,
):
    """
    Filters information in table.

    arguments:
        table (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        proportion_validity (float): proportion of values of intensity that must
            have valid, non-missing values to keep the row for each protein
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table_filter = table.copy(deep=True)
    # Copy other information.
    columns_intensity = copy.deepcopy(columns_intensity)

    # Filter information in table.

    # Filter rows in table.

    ##########
    # Filter rows in table on basis of protein identification.
    table_filter["match_keep"] = table_filter.apply(
        lambda row:
            match_keep_table_row_identification(
                confidences_keep=["high", "medium",],
                peptides_keep=1,
                peptides_unique_keep=1,
                row_identifier=row["identifier_protein_uniprot"],
                row_confidence=row["confidence"],
                row_peptides=row["coverage_peptides"],
                row_peptides_unique=row["coverage_peptides_unique"],
                report=False,
            ),
        axis="columns", # apply function to each row
    )
    table_filter = table_filter.loc[
        (table_filter["match_keep"] == 1), :
    ]
    # Remove unnecessary columns.
    table_filter.drop(
        labels=["match_keep",],
        axis="columns",
        inplace=True
    )

    ##########
    # Filter rows in table on basis of protein quantification.
    table_filter["match_keep"] = table_filter.apply(
        lambda row:
            match_keep_table_row_quantification(
                row=row,
                columns_intensity=columns_intensity,
                proportion_validity=proportion_validity,
                report=False,
            ),
        axis="columns", # apply function to each row
    )
    table_filter = table_filter.loc[
        (table_filter["match_keep"] == 1), :
    ]
    # Remove unnecessary columns.
    table_filter.drop(
        labels=["match_keep",],
        axis="columns",
        inplace=True
    )

    # Filter and sort columns in table.
    columns_sequence_protein = define_column_sequence_table_protein()
    columns_sequence = copy.deepcopy(columns_sequence_protein)
    columns_sequence.extend(columns_intensity)
    table_filter = porg.filter_sort_table_columns(
        table=table_filter,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_rows = (table_filter.shape[0])
        count_columns = (table_filter.shape[1])
        columns = copy.deepcopy(table_filter.columns.to_list())
        print("Count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=4)
        print("Names of columns: ")
        print(columns)
        putly.print_terminal_partition(level=4)
        print("Table: ")
        print(table_filter)
        putly.print_terminal_partition(level=4)

    # Return information.
    return table_filter


##########
# Split


def split_table_main_columns(
    table=None,
    columns_protein=None,
    columns_intensity=None,
    report=None,
):
    """
    Splits or separates information in table by columns.

    arguments:
        table (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows
        columns_protein (list<str>): names of columns corresponding to
            information about proteins
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table_split = table.copy(deep=True)
    # Copy other information.
    columns_intensity = copy.deepcopy(columns_intensity)
    # Organize information.
    table_split.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_split.set_index(
        ["identifier_protein_uniprot"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Separate information into separate tables.
    table_protein = table_split.loc[
        :, table_split.columns.isin(columns_protein)
    ]
    table_intensity = table_split.loc[
        :, table_split.columns.isin(columns_intensity)
    ]
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        columns_protein = copy.deepcopy(table_protein.columns.to_list())
        print("Table-Protein")
        print("Names of columns:")
        print(columns_protein)
        print("Table:")
        print(table_protein)
        putly.print_terminal_partition(level=4)
        columns_intensity = copy.deepcopy(table_intensity.columns.to_list())
        print("Table-Intensity")
        print("Names of columns:")
        print(columns_intensity)
        print("Table:")
        print(table_intensity)
        pass
    # Collect information.
    pail = dict()
    pail["table_protein"] = table_protein
    pail["table_intensity"] = table_intensity
    # Return information.
    return pail


##########
# Fill


def match_table_row_fill_missing(
    row=None,
    columns_intensity=None,
    report=None,
):
    """
    Determines whether the current row from a table needs fill of missing
    values.

    arguments:
        row (object): Pandas series of values of intensity across samples for
            a single protein
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Copy information in row.
    row = row.copy(deep=True)
    # Extract information for current row from table.
    array_intensities = row[columns_intensity].to_numpy(
        dtype="float",
        na_value=numpy.nan,
        copy=True,
    )
    #array_intensities = array_intensities.astype(numpy.float32)
    count_missing = numpy.count_nonzero(numpy.isnan(array_intensities))
    # Determine whether the current row has missing values.
    if (count_missing > 0):
        indicator = 1
    else:
        indicator = 0
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Match missing value fill.")
        print("Indicator: " + str(indicator))
        putly.print_terminal_partition(level=4)
    # Return information.
    return indicator

def fill_missing_values_intensity_row(
    row=None,
    column=None,
    columns_intensity=None,
    method=None,
    report=None,
):
    """
    Determines whether the current row from a table needs fill of missing
    values.

    Within a normal distribution, 99.7% of values occur within 3 standard
    deviations of the mean, either above or below.
    3 standard deviations: 99.73%
    2 standard deviations: 95.45%
    1 standard deviations: 68.27%

    Reference:
    https://en.wikipedia.org/wiki/68-95-99.7_rule

    arguments:
        row (object): Pandas series of values of intensity across samples for
            a single protein
        column (str): name of current column within row
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        method (str): name of method to use for filling missing values, either
            'triple_standard_deviation_below', 'half_minimum', or 'zero'
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Copy information in row.
    row = row.copy(deep=True)
    # Determine whether the current row has a missing value in the current
    # column.
    check_fill = 0
    if (
        (row["match_missing"] == 1) and
        (pandas.isna(row[column]))
    ):
        # Fill missing value.
        check_fill = 1
        # Extract information for current row from table.
        array_intensities = row[columns_intensity].to_numpy(
            dtype="float",
            na_value=numpy.nan,
            copy=True,
        )
        #array_intensities = array_intensities.astype(numpy.float32)
        # Determine which method to use for fill.
        if (method == "triple_standard_deviation_below"):
            mean = numpy.nanmean(array_intensities)
            standard_deviation = numpy.nanstd(
                array_intensities,
                ddof=1, # divisor is (n - 1) for sample standard deviation
            )
            value = float(mean - (3 * standard_deviation))
            pass
        elif (method == "half_minimum"):
            minimum = numpy.nanmin(array_intensities)
            value = float(minimum / 2)
            pass
        elif (method == "zero"):
            value = float(0)
            pass
    else:
        value = row[column]
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Fill missing value check: " + str(check_fill))
        print("Fill value: " + str(value))
        putly.print_terminal_partition(level=4)
    # Return information.
    return value

def fill_missing_values_intensity_table(
    table=None,
    columns_intensity=None,
    method=None,
    report=None,
):
    """
    Filters information in table.

    Within a normal distribution, 99.7% of values occur within 3 standard
    deviations of the mean, either above or below.
    3 standard deviations: 99.73%
    2 standard deviations: 95.45%
    1 standard deviations: 68.27%

    Reference:
    https://en.wikipedia.org/wiki/68-95-99.7_rule

    arguments:
        table (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        method (str): name of method to use for filling missing values, either
            'triple_standard_deviation_below', 'half_minimum', or 'zero'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table_fill = table.copy(deep=True)
    # Copy other information.
    columns_intensity = copy.deepcopy(columns_intensity)

    # Fill missing values across all rows for sample corresponding to
    # current column.
    table_fill["match_missing"] = table_fill.apply(
        lambda row:
            match_table_row_fill_missing(
                row=row,
                columns_intensity=columns_intensity,
                report=False,
            ),
        axis="columns", # apply function to each row
    )
    table_fill_actual = table_fill.loc[
        (table_fill["match_missing"] == 1), :
    ]

    # Iterate across columns for values of intensity in each sample.
    for column in columns_intensity:
        # Fill missing values across all rows for sample corresponding to
        # current column.
        table_fill[column] = table_fill.apply(
            lambda row:
                fill_missing_values_intensity_row(
                    row=row,
                    column=column,
                    columns_intensity=columns_intensity,
                    method=method,
                    report=False,
                ),
            axis="columns", # apply function to each row
        )
        pass

    # Remove unnecessary columns.
    table_fill.drop(
        labels=["match_missing",],
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_rows_missing = (table_fill_actual.shape[0])
        count_rows_total = (table_fill.shape[0])
        proportion_missing = round(
            float(count_rows_missing / count_rows_total), 2
        )
        print(
            "Count of rows requiring missing fill: " + str(count_rows_missing)
        )
        print(
            "Proportion of rows requiring missing fill: " +
            str(proportion_missing)
        )
        putly.print_terminal_partition(level=4)
        print("Table: ")
        print(table_fill)
        putly.print_terminal_partition(level=4)

    # Return information.
    return table_fill




################################################################################
# Procedure


def execute_procedure(
    path_directory_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Initialize directories.
    paths = initialize_directories(
        project="2024_exercise_age",
        set="wrangler_ms",
        path_directory_dock=path_directory_dock,
        restore=True,
    )

    ##########
    # 1. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=True,
    )

    ##########
    # TODO: TCW; 25 June 2024
    # Steps 2-6 will need to be performed separately on the raw table from each
    # separate run before combination in step 7.
    ##########

    ##########
    # 2. Organize information from source.
    pail_organization = organize_table_main(
        table_sample=pail_source["table_sample"],
        table_main=pail_source["table_main"],
        report=True,
    )

    ##########
    # 3. Filter columns and rows in table.
    if False:
        # Test.
        row = define_row_test_filter()
        indicator = match_keep_table_row_identification(
            confidences_keep=["high", "medium",],
            peptides_keep=3,
            peptides_unique_keep=3,
            row_identifier=row["identifier_protein_uniprot"],
            row_confidence=row["confidence"],
            row_peptides=row["coverage_peptides"],
            row_peptides_unique=row["coverage_peptides_unique"],
            report=True,
        )
        pass
    if False:
        # Test.
        row = define_row_test_filter()
        indicator = match_keep_table_row_quantification(
            row=row,
            columns_intensity=pail_organization["samples"],
            proportion_validity=0.80,
            report=True,
        )
        pass
    table_filter = filter_table_main(
        table=pail_organization["table_main"],
        columns_intensity=pail_organization["samples"],
        proportion_validity=0.90,
        report=True,
    )

    ##########
    # 4. Separate tables for information of discrete types.
    # table_sample
    # table_protein
    # table_intensity
    columns_protein = define_column_sequence_table_protein()
    columns_protein.remove("identifier_protein_uniprot")
    pail_split = split_table_main_columns(
        table=table_filter,
        columns_protein=columns_protein,
        columns_intensity=pail_organization["samples"],
        report=True,
    )

    ##########
    # 5. Fill missing values of intensity.
    # - performed row-by-row (protein-by-protein) across columns (samples)
    # - fill missing values with 3 SD below minimum
    # methods:
    # "triple_standard_deviation_below"
    # "half_minimum"
    # "zero"
    if False:
        # Test.
        # values: 0.4526, 0.3325, 0.2953, 0.3015, 0.4379, 0.3127, 0.1578, 0.2351
        # mean: 0.3157
        # standard deviation: 0.0971
        # (mean - (3 * standard deviation)): 0.0243
        # minimum: 0.1578
        # (minimum / 2): 0.0789
        row = define_row_test_filter()
        row["match_missing"] = 1
        value = fill_missing_values_intensity_row(
            row=row,
            column="control_3", # "control_3", "intervention_1"
            columns_intensity=pail_organization["samples"],
            method="triple_standard_deviation_below",
            report=True,
        )
        pass
    table_intensity = fill_missing_values_intensity_table(
        table=pail_split["table_intensity"],
        columns_intensity=pail_organization["samples"],
        method="triple_standard_deviation_below",
        report=True,
    )

    ##########
    # 6. Scale values of intensity.
    # - median scaling
    # - performed row-by-row (protein-by-protein) across columns (samples)

    ##########
    # 7. Combine tables from multiple runs.
    #   This combination must happen after standardizing the scales of each run.




    ##########
    # 7. Normalize values of intensity.



    ##########
    # 8. Transform table of intensities to wide format.

    # TODO: TCW; 24 June 2024
    # TODO: AFTER filling missing values
    # TODO: AFTER scaling values (median normalization or scaling).
    # TODO: AFTER normalizing scaled intensity values
    #    The median scaling will be performed row-by-row across proteins.
    #    It will be more convenient for the proteins to be oriented across rows.
    # TODO: I think that I'll eventually need to have the sample intensities
    # across rows.
    # TODO: WHY?????
    # TODO: 1. I need to transfer more information from the Table-Sample.
    # TODO:     Specifically "group" and "pair" need to be transfered according to sample identifier (identifier_novel)
    # TODO: 2. The regression analyses will be more intuitive to design with
    #          samples across rows.

    ##########
    # 9. Transfer information about experimental groups and pairs.



    ##########
    # Draft material for new "inquisitor_ms" module...
    # 10. Analyze values across samples with respect for pairs of observations.
    #     - Run this with parallelization across the columns for proteins.
    #     - Pass the list of columns for proteins to the parallel function along with full table


    pass


################################################################################
# End
