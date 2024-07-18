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
        paths["dock"], "in_data", "exercise", str(project), str(set),
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
    paths["out_table"] = os.path.join(
        paths["out_set"], "table",
    )
    paths["out_plot"] = os.path.join(
        paths["out_set"], "plot",
    )
    paths["out_plot_raw"] = os.path.join(
        paths["out_plot"], "raw",
    )
    paths["out_plot_scale_sample"] = os.path.join(
        paths["out_plot"], "scale_sample",
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        putly.remove_directory(path=paths["out_project"])
        putly.remove_directory(path=paths["out_set"])
        putly.remove_directory(path=paths["out_test"])
        putly.remove_directory(path=paths["out_table"])
        putly.remove_directory(path=paths["out_plot"])
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
    putly.create_directories(
        path=paths["out_table"]
    )
    putly.create_directories(
        path=paths["out_plot"]
    )
    putly.create_directories(
        path=paths["out_plot_raw"]
    )
    putly.create_directories(
        path=paths["out_plot_scale_sample"]
    )
    # Return information.
    return paths


##########
# 1. Read information from file.


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
    types_columns["Intensity.1"] = "float64"
    types_columns["Intensity.2"] = "float64"
    types_columns["Intensity.3"] = "float64"
    types_columns["Intensity.4"] = "float64"
    types_columns["Intensity.5"] = "float64"
    types_columns["Intensity.6"] = "float64"
    types_columns["Intensity.7"] = "float64"
    types_columns["Intensity.8"] = "float64"
    types_columns["Intensity.9"] = "float64"
    types_columns["Intensity.10"] = "float64"
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
    types_columns["batch"] = "string"
    types_columns["sort_batch"] = "int32"
    types_columns["group"] = "string"
    types_columns["sort_group"] = "int32"
    types_columns["pair"] = "string"
    types_columns["sort_pair"] = "int32"
    types_columns["identifier_original"] = "string"
    types_columns["identifier_novel"] = "string"
    types_columns["description"] = "string"
    types_columns["abbreviation"] = "string"
    types_columns["quantity"] = "float32"
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
        paths["in_parameters"], "table_samples_attributes.tsv",
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
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Table of attributes of samples.
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

    # Return information.
    return pail


##########
# 2. Organize information.


def extract_organize_values_from_series(
    series=None,
    report=None,
):
    """
    Extract and organize values from a Pandas series, such as a single column
    or row from a Pandas data-frame table.

    arguments:
        series (object): Pandas series of values of intensity
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in series.
    series = series.copy(deep=True)
    # Extract and organize information from series.
    #values_raw = values_raw.astype(numpy.float32)
    values_raw = series.to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    values_positive = values_raw
    values_positive[values_positive < 0] = numpy.nan
    values_valid = values_positive[~numpy.isnan(values_positive)]
    values_log = numpy.log(values_valid)
    # Collect information.
    pail = dict()
    pail["values_raw"] = values_raw
    pail["values_positive"] = values_positive
    pail["values_valid"] = values_valid
    pail["values_log"] = values_log
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("Report:")
        print("exercise")
        print("proteomics_mass_spectroscopy")
        print("extract_organize_values_from_series()")
        putly.print_terminal_partition(level=4)
        print("Original source series:")
        print(series)
        putly.print_terminal_partition(level=4)
        print("Raw values:")
        print(values_raw)
        putly.print_terminal_partition(level=4)
        print("Positive values:")
        print(values_positive)
        putly.print_terminal_partition(level=4)
        print("Valid values:")
        print(values_valid)
        putly.print_terminal_partition(level=4)
        print("Logarithmic values:")
        print(values_log)
    # Return information.
    return pail


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
        table_main (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)
    table_main = table_main.copy(deep=True)

    # Filter information about samples.
    table_sample_inclusion = table_sample.loc[
        (table_sample["inclusion"] == 1), :
    ]
    table_sample_control = table_sample_inclusion.loc[
        (table_sample_inclusion["group"] == "control"), :
    ]
    table_sample_intervention = table_sample_inclusion.loc[
        (table_sample_inclusion["group"] == "intervention"), :
    ]
    # Extract information about samples.
    samples = copy.deepcopy(
        table_sample_inclusion["identifier_novel"].to_list()
    )
    samples_control = copy.deepcopy(
        table_sample_control["identifier_novel"].to_list()
    )
    samples_intervention = copy.deepcopy(
        table_sample_intervention["identifier_novel"].to_list()
    )

    # Determine original and novel names for columns representing samples.
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

    # Replace values of zero with missing values.
    # Only replace values within table's columns for samples.
    # This implementation is more concise than iteration across specific
    # columns.
    table_main[samples] = table_main[samples].replace(
        to_replace=0.0,
        value=pandas.NA,
    )
    # Replace values less than zero with missing values.
    table_main[samples][table_main[samples] < 0] = pandas.NA

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
    pail["samples_control"] = samples_control
    pail["samples_intervention"] = samples_intervention
    pail["table_main"] = table_main
    # Return information.
    return pail


##########
# 3. Filter rows and columns in table.


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
        proportion_validity (float): proportion of values of intensity that
            must have non-missing, valid values in order to keep the row for
            each protein
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Copy information in row.
    row = row.copy(deep=True)
    # Extract and organize information from series.
    pail_values = extract_organize_values_from_series(
        series=row[columns_intensity],
        report=False,
    )
    values_nonmissing = pail_values["values_raw"][~numpy.isnan(
        pail_values["values_raw"]
    )]
    #count_nonmissing = numpy.count_nonzero(~numpy.isnan(
    #    pail_values["values_raw"]
    #))
    count_nonmissing = int(values_nonmissing.size)
    count_validity = numpy.count_nonzero(values_nonmissing > 0)
    count_samples = len(columns_intensity)
    count_invalidity = (count_samples - count_validity)
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
        print(
            "Match function for filter by quantification of intensity " +
            "values across samples for a single protein."
        )
        putly.print_terminal_partition(level=4)
        print("Array of intensity values accross samples:")
        print("Count of samples: " + str(count_samples))
        print("Count of nonmissing: " + str(count_nonmissing))
        print("Count of valid: " + str(count_validity))
        print("Count of invalid: " + str(count_invalidity))
        putly.print_terminal_partition(level=4)
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
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        proportion_validity (float): proportion of values of intensity that
            must have non-missing, valid values in order to keep the row for
            each protein
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
# 4. Separate tables for information of distinct types.


def split_table_main_columns(
    table=None,
    columns_protein=None,
    columns_intensity=None,
    report=None,
):
    """
    Splits or separates information in table by columns.

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
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
    # Organize indices in table.
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
    # Organize indices in table.
    table_protein.columns.rename(
        "attribute",
        inplace=True,
    ) # single-dimensional index
    table_intensity.columns.rename(
        "sample",
        inplace=True,
    ) # single-dimensional index

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
# 5. Fill missing values of intensity.


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
    # Extract and organize information from series.
    pail_values = extract_organize_values_from_series(
        series=row[columns_intensity],
        report=False,
    )
    #count_missing = numpy.count_nonzero(
    #    numpy.isnan(pail_values["values_raw"])
    #)
    # Determine whether the current row has missing values.
    if (pail_values["values_valid"].size < len(columns_intensity)):
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

    For the 'triple_standard_deviation_below' method, this function transforms
    raw values via natural logarithm before determining their mean and
    standard deviation so that these statistics are representative of a more
    normal distribution. This function calculates the intermediate fill value
    and then inverts the transformation by exponentiation.

    Reference:
    https://en.wikipedia.org/wiki/68-95-99.7_rule

    Review:
    On 27 June 2024, TCW confirmed that the procedure was actually recognizing
    and filling missing values appropriately. TCW also used a separate
    spreadsheet to check the math that determines the fill value via the
    'triple_standard_deviation_below' method. As expected, the fill value
    differed by whether the method transformed by logarithm before calculation
    of mean and standard deviation.

    arguments:
        row (object): Pandas series of values of intensity across samples for
            a single protein
        columns_intensity (list<str>): names of columns corresponding to values
            of intensity
        method (str): name of method to use for filling missing values, either
            'triple_standard_deviation_below', 'half_minimum', or 'zero'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas series of values of intensity across samples for
            a single protein

    """

    # Copy information in row.
    row = row.copy(deep=True)

    # Determine whether the current row has a missing value in the current
    # column.
    check_fill = 0
    value_fill = float("nan")
    if (row["match_missing"] == 1):
        # Fill missing value.
        check_fill = 1
        # Extract and organize information from series.
        pail_values = extract_organize_values_from_series(
            series=row[columns_intensity],
            report=report,
        )
        # Determine which method to use for fill.
        if (method == "triple_standard_deviation_below"):
            mean = numpy.nanmean(pail_values["values_log"])
            standard_deviation = numpy.nanstd(
                pail_values["values_log"],
                ddof=1, # divisor is (n - 1) for sample standard deviation
            )
            value_fill_log = float(mean - (3 * standard_deviation))
            value_fill = numpy.exp(value_fill_log)
            row_fill = row.replace(
                to_replace=pandas.NA,
                value=value_fill,
            )
        elif (method == "half_minimum"):
            minimum = numpy.nanmin(array_intensities)
            value_fill = float(minimum / 2)
            row_fill = row.replace(
                to_replace=pandas.NA,
                value=value_fill,
            )
        elif (method == "zero"):
            value_fill = 0.0
            row_fill = row.replace(
                to_replace=pandas.NA,
                value=value_fill,
            )
        # Report.
        if report:
            if (value_fill <= 0):
                putly.print_terminal_partition(level=1)
                print("**********")
                print("WARNING!")
                print("**********")
                print(
                    "Method selection to fill missing values returned a value " +
                    "less than or equal to zero!"
                )
                print("**********")
                putly.print_terminal_partition(level=1)
            putly.print_terminal_partition(level=4)
            print("Fill missing value check: " + str(check_fill))
            print("Fill value: " + str(value_fill))
            putly.print_terminal_partition(level=4)
            print("Row with fill values: " + str(row_fill))
            putly.print_terminal_partition(level=4)
            pass
        pass
    else:
        row_fill = row
        pass

    # Report.
    if report:
        if (value_fill <= 0):
            putly.print_terminal_partition(level=1)
            print("**********")
            print("WARNING!")
            print("**********")
            print(
                "Method selection to fill missing values returned a value " +
                "less than or equal to zero!"
            )
            print("**********")
            putly.print_terminal_partition(level=1)
        putly.print_terminal_partition(level=4)
        print("Fill missing value check: " + str(check_fill))
        print("Fill value: " + str(value_fill))
        putly.print_terminal_partition(level=4)
        print("Row with fill values: " + str(row_fill))
        putly.print_terminal_partition(level=4)
    # Return information.
    return row_fill


def fill_missing_values_intensity_table(
    table=None,
    columns_intensity=None,
    method=None,
    report=None,
):
    """
    Filters information in table.

    Table's format and orientation

    Table has values of intensity for each protein oriented across rows with
    samples oriented across columns.

                sample_1 sample_2 sample_3 sample_4 sample_5
    protein_1   ...      ...      ...      ...      ...
    protein_2   ...      ...      ...      ...      ...
    protein_3   ...      ...      ...      ...      ...
    protein_4   ...      ...      ...      ...      ...
    protein_5   ...      ...      ...      ...      ...

    Refer to function 'calculate_standard_score_gene_signal_by_gene' from module
    'tissue' within package 'bimodality'.

    Within a normal distribution, 99.7% of values occur within 3 standard
    deviations of the mean, either above or below.
    3 standard deviations: 99.73%
    2 standard deviations: 95.45%
    1 standard deviations: 68.27%

    Reference:
    https://en.wikipedia.org/wiki/68-95-99.7_rule

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
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

    # Fill missing values for all columns and across all rows.
    # Apply the function to each row.
    table_fill = table_fill.apply(
        lambda row:
            fill_missing_values_intensity_row(
                row=row,
                columns_intensity=columns_intensity,
                method=method,
                report=False,
            ),
        axis="columns", # apply function to each row
    )
    # Remove unnecessary columns.
    table_fill.drop(
        labels=["match_missing",],
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Report: Fill missing values")
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


##########
# 6. Scale intensity by quantity.



##########
# 7. Scale intensity by batch.



##########
# 8. Combine intensities for samples from separate batches.



##########
# 9. Scale intensity by sample.


def shift_values_greater_zero_row(
    row=None,
):
    """
    Shifts values by minimum if necessary to ensure that all values are greater
    that zero.

    This arithmetic shift would distort the data. Currently this function is
    not in use, but TCW kept the function for reference.

    Review: 27 June 2024

    arguments:
        row (object): Pandas series of values for observations of a single
            feature

    raises:

    returns:
        (object): Pandas series of values for observations of a single feature

    """

    # Copy information in row.
    row = row.copy(deep=True)
    # Extract information for current row from table.
    values = row.to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    #array_intensities = array_intensities.astype(numpy.float32)
    minimum = numpy.nanmin(values)
    minimum_absolute = math.fabs(minimum)
    # Determine whether it is necessary to shift values.
    if (minimum < 0):
        row_shift = row.add(minimum_absolute)
    else:
        row_shift = row
        pass
    # Return information.
    return row_shift


def scale_values_intensity_table(
    table=None,
    method=None,
    report=None,
):
    """
    Scales values in table.

    Table's format and orientation

    Table has values of intensity for each protein oriented across rows with
    samples oriented across columns.

                sample_1 sample_2 sample_3 sample_4 sample_5
    protein_1   ...      ...      ...      ...      ...
    protein_2   ...      ...      ...      ...      ...
    protein_3   ...      ...      ...      ...      ...
    protein_4   ...      ...      ...      ...      ...
    protein_5   ...      ...      ...      ...      ...

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
        method (str): name of method to use for scaling values, currently only
            'median_ratio'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Organize indices in table.
    table = porg.change_names_table_indices_columns_rows(
        table=table,
        name_columns_novel="observations",
        name_rows_original="identifier_protein_uniprot",
        name_rows_novel="features",
        report=False,
    )
    # Determine method for scaling.
    if (method == "median_ratio"):
        table_scale = (
            pscl.scale_feature_values_between_observations_by_median_ratio(
                table=table,
                name_columns="observations",
                name_rows="features",
                report=report,
        ))
    # Organize indices in table.
    table_scale = porg.change_names_table_indices_columns_rows(
        table=table_scale,
        name_columns_novel="sample",
        name_rows_original="features",
        name_rows_novel="identifier_protein_uniprot",
        report=False,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Report:")
        print("exercise")
        print("proteomics_mass_spectroscopy")
        print("scale_values_intensity_table()")
        putly.print_terminal_partition(level=4)
        print(table_scale)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_scale


def compare_least_change_proteins_by_experimental_groups(
    table=None,
    name_columns=None,
    name_rows=None,
    columns_first=None,
    columns_second=None,
    count_quantile=None,
    report=None,
):
    """
    Compares the sets of proteins that demonstrate the least change between
    samples corresponding to two experiomental groups, such as control and
    intervention.

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        columns_first (list<str>): names of columns corresponding to samples in
            first experimental group, such as control or intervention
        columns_second (list<str>): names of columns corresponding to samples
            in second experimental group, such as control or intervention
        count_quantile (int): odd count of quantiles, such as three for
            tertiles
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Separate table by columns corresponding to samples in each group.
    table_first = table[columns_first].copy(deep=True)
    table_second = table[columns_second].copy(deep=True)
    # Compare sets of proteins that demonstrate least change between samples in
    # each experimental group.
    pail = (
        pscl.compare_middle_quantile_feature_sets_by_ratio_to_geometric_mean(
            table_first=table_first,
            table_second=table_second,
            name_columns=name_columns,
            name_rows=name_rows,
            count_quantile=count_quantile,
            report=report,
    ))
    # For the sets of proteins that demonstrate Describe the variance of the ratios of change across all samples.
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("Report:")
        print("exercise")
        print("proteomics_mass_spectroscopy")
        print("compare_constant_proteins_by_experimental_groups()")
        putly.print_terminal_partition(level=4)
        print("Columns for samples in first group:")
        print(columns_first)
        print("Columns for samples in second group:")
        print(columns_second)
    # Return information.
    pass


##########
# 10. Normalize distribution in values of intensity.


def normalize_values_intensity_table(
    table=None,
    columns=None,
    logarithm=None,
    z_score=None,
    report=None,
):
    """
    Normalize values of intensity in table.

    Table's format and orientation

    Table has values of intensity for each protein oriented across rows with
    samples oriented across columns.

                sample_1 sample_2 sample_3 sample_4 sample_5
    protein_1   ...      ...      ...      ...      ...
    protein_2   ...      ...      ...      ...      ...
    protein_3   ...      ...      ...      ...      ...
    protein_4   ...      ...      ...      ...      ...
    protein_5   ...      ...      ...      ...      ...

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
        columns (list<str>): names of columns corresponding to values of
            intensity of each protein across samples
        logarithm (bool): whether to transform values by logarithm
        z_score (bool): whether to transform values by z-score standardization
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table_normal = table.copy(deep=True)
    # Determine method for scaling.
    if (logarithm):
        for column in columns:
            # math.log() # optimal for scalar values
            # numpy.log() # optimal for array values
            table_normal[column] = table_normal.apply(
                lambda row:
                    math.log(row[column])
                    if row[column] != 0.0 else pandas.NA,
                axis="columns", # apply function to each row
            )
            pass
        pass
    if (z_score):
        table_normal = table_normal.apply(
            lambda row:
                scipy.stats.zscore(
                    row[columns],
                    axis=0,
                    ddof=1, # divisor is (n - 1) for sample standard deviation
                    nan_policy="omit",
                ),
            axis="columns", # apply function to each row
        )
        pass

    # Calculate mean values for each protein across samples.
    means = table_normal.aggregate(
        lambda row: numpy.nanmean(
            row.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
        ),
        axis="columns", # apply function to each row
    )
    # Calculate mean values for each protein across samples.
    standard_deviations = table_normal.aggregate(
        lambda row: numpy.nanstd(
            row.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            ),
            ddof=1, # divisor is (n - 1) for sample standard deviation
        ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Report:")
        print("exercise")
        print("proteomics_mass_spectroscopy")
        print("normalize_values_intensity_table()")
        putly.print_terminal_partition(level=5)
        print(table_normal)
        putly.print_terminal_partition(level=5)
        print("Means across proteins")
        print(means)
        putly.print_terminal_partition(level=5)
        print("Standard deviations across proteins")
        print(standard_deviations)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_normal


##########
# 11. Prepare information for analysis.


def prepare_table_for_analysis(
    table_sample=None,
    table_protein=None,
    table_intensity=None,
    report=None,
):
    """
    Prepares table of information for analysis.

    arguments:
        table_sample (object): Pandas data-frame table of information about
            samples
        table_protein (object): Pandas data-frame table of information about
            proteins
        table_intensity (object): Pandas data-frame table of values of
            intensity for samples across columns and for proteins across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)
    table_protein = table_protein.copy(deep=True)
    table_intensity = table_intensity.copy(deep=True)
    # Extract information.
    proteins = copy.deepcopy(
        table_protein.index.get_level_values(
            "identifier_protein_uniprot"
        ).to_list()
    )
    # Transpose table.
    table_main = table_intensity.transpose(copy=True)
    # Organize indices in table.
    table_sample.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_main.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Transfer extra attributes of samples.
    attributes = [
        "sort", "batch", "group", "pair", "description", "abbreviation",
    ]
    table_main = porg.transfer_table_rows_attributes_reference(
        table_main=table_main,
        column_main_key="sample",
        table_reference=table_sample,
        column_reference_key="identifier_novel",
        columns_reference_transfer=attributes,
        prefix_reference_main="",
        suffix_reference_main="",
        report=report,
    )
    # Sort rows within table.
    table_main.sort_values(
        by=["sort",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Sort columns within table.
    attributes.remove("sort")
    attributes.insert(0, "sample")
    columns = copy.deepcopy(attributes)
    columns.extend(proteins)
    table_main = porg.filter_sort_table_columns(
        table=table_main,
        columns_sequence=columns,
        report=report,
    )
    # Remove unnecessary columns.
    if False:
        table_main.drop(
            labels=["sort"],
            axis="columns",
            inplace=True
        )
    # Organize indices in table.
    table_main.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_main.set_index(
        ["sample"],
        append=False,
        drop=True,
        inplace=True,
    )
    table_main.columns.rename(
        "protein_attribute",
        inplace=True,
    ) # single-dimensional index
    table_main = porg.change_names_table_indices_columns_rows(
        table=table_main,
        name_columns_novel="feature",
        name_rows_original="sample",
        name_rows_novel="observation",
        report=False,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Report:")
        print("exercise")
        print("proteomics_mass_spectroscopy")
        print("prepare_table_for_analysis()")
        putly.print_terminal_partition(level=5)
        print("Table-Sample")
        print(table_sample)
        putly.print_terminal_partition(level=5)
        print("Table-Protein")
        print(table_protein)
        putly.print_terminal_partition(level=5)
        print("Table-Intensity")
        print(table_intensity)
        putly.print_terminal_partition(level=5)
        print("Table-Main")
        print(table_main)
        pass
    # Collect information.
    pail = dict()
    pail["table_main"] = table_main
    pail["columns"] = columns
    pail["attributes"] = attributes
    pail["proteins"] = proteins
    # Return information.
    return pail




##########
# Plot values of intensity in histogram
# Values of intensity across all proteins for each sample


def create_write_figure_histogram(
    values=None,
    name_figure=None,
    path_directory=None,
    report=None,
):
    """
    Creation dot plot and write to file.

    Assume that rows for records in the source table are already in proper sort
    order corresponding to the categorical labels on the abscissa (horizontal
    axis).

    arguments:
        values (object): NumPy array of values
        name_figure (str): name of figure to use in file name
        path_directory (str): path to parent directory within which to write a
            file for figure
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Organize information for figure.
    mean = numpy.nanmean(values)

    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_distribution_histogram(
        array=values,
        title=name_figure,
        bin_method="count",
        bin_count=30,
        bar_width=0.5,
        label_bins="Natural Log (Intensity)",
        label_counts="Protein Instances per bin",
        fonts=fonts,
        colors=colors,
        line=True,
        line_position=mean,
        label_title=name_figure,
        label_report=True,
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


def plot_histogram_values_across_proteins_each_sample(
    table=None,
    name_columns=None,
    name_rows=None,
    logarithm=None,
    path_directory_parent=None,
    report=None,
):
    """
    Plots in a histogram values of intensity for all proteins in each sample.

    Table's format and orientation

    Table has values of intensity for each protein oriented across rows with
    samples oriented across columns.

                sample_1 sample_2 sample_3 sample_4 sample_5
    protein_1   ...      ...      ...      ...      ...
    protein_2   ...      ...      ...      ...      ...
    protein_3   ...      ...      ...      ...      ...
    protein_4   ...      ...      ...      ...      ...
    protein_5   ...      ...      ...      ...      ...

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for proteins across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        logarithm (bool): whether to transform values logarithmically to
            normalize their distribution
        path_directory_parent (str): path to parent directory within which to
            write files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy names of columns in original table.
    names_columns = copy.deepcopy(
        table.columns.get_level_values(name_columns).to_list()
    )
    names_rows = copy.deepcopy(
        table.index.get_level_values(name_rows).to_list()
    )
    # Create separate plot for each sample.
    for column in names_columns:
        # Extract and organize information from series.
        pail_values = extract_organize_values_from_series(
            series=table[column],
            report=False,
        )
        if logarithm:
            values_plot = pail_values["values_log"]
        else:
            values_plot = pail_values["values_valid"]
        # Create figure and write to file.
        create_write_figure_histogram(
            values=values_plot,
            name_figure=column,
            path_directory=path_directory_parent,
            report=True,
        )
        pass
    # Return information.
    pass




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

    # Initialize directories.
    paths = initialize_directories(
        project="2024_exercise_age",
        set="preparation",
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
    # separate batch run before combination in step 7.
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
    table_filter = filter_table_main(
        table=pail_organization["table_main"],
        columns_intensity=pail_organization["samples"],
        proportion_validity=0.9, # keep only rows without any missing values
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
    putly.print_terminal_partition(level=3)
    print("Step 5: Fill missing values of intensity.")
    putly.print_terminal_partition(level=3)
    if False:
        # Test.
        # values: 0.4526, 0.3325, 0.2953, 0.3015, 0.4379, 0.3127, 0.1578, 0.2351
        # mean: 0.3157
        # standard deviation: 0.0971
        # (mean - (3 * standard deviation)): 0.0243
        # minimum: 0.1578
        # (minimum / 2): 0.0789
        # Samples with missing values: "control_3", "intervention_1"
        row = define_row_test_filter()
        row["match_missing"] = 1
        row_fill = fill_missing_values_intensity_row(
            row=row,
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
    if True:
        # Describe variance of ratios to geometric mean across the proteins that
        # demonstrate the least change.
        pail = pscl.describe_variance_across_features_with_least_change(
            table=table_intensity,
            name_columns="sample",
            name_rows="identifier_protein_uniprot",
            count_quantile=7,
            report=True,
        )
        # Compare sets of proteins that demonstrate least change.
        compare_least_change_proteins_by_experimental_groups(
            table=table_intensity,
            name_columns="sample",
            name_rows="identifier_protein_uniprot",
            columns_first=pail_organization["samples_control"],
            columns_second=pail_organization["samples_intervention"],
            count_quantile=7,
            report=True,
        )
        plot_histogram_values_across_proteins_each_sample(
            table=table_intensity,
            name_columns="sample",
            name_rows="identifier_protein_uniprot",
            logarithm=True,
            path_directory_parent=paths["out_plot_raw"],
            report=True,
        )
        pass


    ##########
    # 6. Scale values of intensity by batch.
    # - performed row-by-row (protein-by-protein) across columns (samples)

    # TODO: TCW; 25 June 2024
    # "Batch Scaling"
    # The goal of this scaling is to make it reasonable to combine samples from separate
    # batches or "runs" on the mass spec
    # I think that calculating ratios of each protein in each sample to the corresponding
    # protein in the pooled bridge sample from the same batch will work.

    ##########
    # 7. Combine tables from multiple batches.
    #   This combination must happen after standardizing the scales of each run.

    # TODO: TCW; 1 July 2024
    # Use the protein identifiers as index keys
    # Combine the two, batch-specific tables by addition columns
    # partner.organization.merge_columns_two_tables()


    ##########
    # 8. Scale overall values of intensity by sample.
    # The goal of this scaling is to make the individual samples more
    # comparable to each other.
    # There is potential for drift in instrument sensitivity even between runs
    # of individual samples.
    # This scaling can decrease the variance or noise in measurements between
    # samples.
    putly.print_terminal_partition(level=3)
    print("Step 8: Scale overall values of intensity by sample.")
    putly.print_terminal_partition(level=3)
    table_scale = scale_values_intensity_table(
        table=table_intensity,
        method="median_ratio",
        report=True,
    )
    if True:
        # Describe variance of ratios to geometric mean across the proteins that
        # demonstrate the least change.
        pail = pscl.describe_variance_across_features_with_least_change(
            table=table_scale,
            name_columns="sample",
            name_rows="identifier_protein_uniprot",
            count_quantile=7,
            report=True,
        )
        # Compare sets of proteins that demonstrate least change.
        compare_least_change_proteins_by_experimental_groups(
            table=table_scale,
            name_columns="sample",
            name_rows="identifier_protein_uniprot",
            columns_first=pail_organization["samples_control"],
            columns_second=pail_organization["samples_intervention"],
            count_quantile=7,
            report=True,
        )
        plot_histogram_values_across_proteins_each_sample(
            table=table_scale,
            name_columns="sample",
            name_rows="identifier_protein_uniprot",
            logarithm=True,
            path_directory_parent=paths["out_plot_scale_sample"],
            report=True,
        )
        pass


    ##########
    # 9. Normalize distribution in values of intensity.
    # The goal of this normalization is different than the scaling above.
    # The goal of this normalization is to make the distributions of values of
    # intensity more usable in subsequent analyses.
    # Other methods of normalization.
    # robust scaling ((x - median(x) / interquartile_range(x)))
    # logarithmic transformation
    #   mitigate skewness, especially left
    # square root transformation
    #   mitigate skewness in presence of zeros
    # square transformation
    #   mitigate skewness, especially left
    # exponential transformation
    #   mitigate skewness, especially when relationship is exponential
    putly.print_terminal_partition(level=3)
    print("Step 9: Normalize distribution in values of intensity.")
    putly.print_terminal_partition(level=3)
    table_normal = normalize_values_intensity_table(
        table=table_scale,
        columns=pail_organization["samples"],
        logarithm=True,
        z_score=True,
        report=True,
    )


    ##########
    # 10. Prepare information for analysis.
    # Transpose the table of intensities to orient proteins across columns and
    # samples across rows.
    # Transfer information about samples to the table.
    putly.print_terminal_partition(level=3)
    print("Step 10: Prepare information for analysis.")
    putly.print_terminal_partition(level=3)
    pail_preparation = prepare_table_for_analysis(
        table_sample=pail_source["table_sample"],
        table_protein=pail_split["table_protein"],
        table_intensity=table_normal,
        report=True,
    )



    ##########
    # Draft material for new "inquisitor_ms" module...
    # 10. Analyze values across samples with respect for pairs of observations.
    #     - Run this with parallelization across the columns for proteins.
    #     - Pass the list of columns for proteins to the parallel function along with full table

    ##########
    # Collect information.
    # Collections of files.
    pail_write_tables = dict()
    pail_write_tables[str("table_main")] = pail_preparation["table_main"]
    pail_write_objects = dict()
    pail_write_objects[str("columns")] = pail_preparation["columns"]
    pail_write_objects[str("attributes")] = pail_preparation["attributes"]
    pail_write_objects[str("proteins")] = pail_preparation["proteins"]
    pail_write_objects[str("samples")] = pail_organization["samples"]

    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_table"],
        type="text",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_table"],
        type="pickle",
    )
    putly.write_objects_to_file_pickle(
        pail_write=pail_write_objects,
        path_directory=paths["out_table"],
    )
    pass


################################################################################
# End
