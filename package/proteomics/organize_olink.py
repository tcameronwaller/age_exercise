"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'organize_olink' is part of the 'proteomics' package within the
'age_exercise' package.

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
import age_exercise.phenotypes.organize_subject as aexph_sub

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
    paths["out_procedure_data"] = os.path.join(
        paths["out_procedure"], "data",
    )
    paths["out_procedure_data_lists"] = os.path.join(
        paths["out_procedure_data"], "lists",
    )
    paths["out_procedure_data_tables"] = os.path.join(
        paths["out_procedure_data"], "tables",
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
        paths["out_procedure_data_lists"],
        paths["out_procedure_data_tables"],
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
        print("module: age_exercise.proteomics.organize_olink.py")
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
    path_file_table_feature_organization = os.path.join(
        paths["in_data"], "study_age_exercise", "subject_sample",
        "table_subject_sample_feature_organization.tsv",
    )
    path_file_table_subject = os.path.join(
        paths["out_project"], "phenotypes", "organize_subject", "data",
        "tables", "table_subject.pickle",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of parameters for organization of the table of attributes for
    # subjects and samples.
    types_columns = (
        aexph_sub.define_type_columns_table_subject_feature_organization()
    )
    pail["table_feature_organization"] = pandas.read_csv(
        path_file_table_feature_organization,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail_parse = aexph_sub.parse_extract_table_sample_feature_organization(
        table=pail["table_feature_organization"],
        inclusion="inclusion_proteomics",
        report=report,
    )
    pail["columns_quantitative"] = pail_parse["columns_quantitative"]
    pail["columns_olink_plasma"] = pail_parse["columns_olink_plasma"]
    pail["columns_olink_muscle"] = pail_parse["columns_olink_muscle"]
    pail["columns_olink_adipose"] = pail_parse["columns_olink_adipose"]

    # Table of properties for subjects.
    pail["table_subject"] = pandas.read_pickle(
        path_file_table_subject,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.organize_sample.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of properties for subjects: ")
        print(pail["table_subject"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Filter and fill information to prepare table for analysis.


def filter_fill_table_subject_olink(
    table=None,
    index_columns=None,
    index_rows=None,
    columns_quantitative=None,
    columns_olink_plasma=None,
    columns_olink_muscle=None,
    columns_olink_adipose=None,
    fill_missing=None,
    report=None,
):
    """
    Filter columns and rows in table and optionally fill missing values. This
    procedure focuses on signal intensities from measurements of peptides by
    Olink technology.

    By intentional design, this function does not apply any transformation of
    scale or normalization of distribution to the values of signal intensity.

    Format of source table (name: "table_source")
    Format of source table is in wide format with floating-point values of
    signal intensities or other measurement values corresponding to features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    observation     feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    Review: 11 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        columns_quantitative (list<str>): names of columns in original source
            table for a selection of features on a quantitative, continuous,
            interval or ratio scale of measurement
        columns_olink_plasma (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in plasma
        columns_olink_muscle (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in muscle
        columns_olink_adipose (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in adipose
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    table_source = table.copy(deep=True)
    columns_quantitative = copy.deepcopy(columns_quantitative)
    columns_olink_plasma = copy.deepcopy(columns_olink_plasma)
    columns_olink_muscle = copy.deepcopy(columns_olink_muscle)
    columns_olink_adipose = copy.deepcopy(columns_olink_adipose)

    ##########
    # Filter rows in table by specific categorical values in columns for
    # specific features.
    # Olink measurements are only available for "Pre" or "first" clinical
    # visits of the study.
    columns_categories = dict()
    columns_categories["age_cohort_text"] = ["younger", "elder",]
    #columns_categories["age_cohort_text"] = ["younger",]
    #columns_categories["age_cohort_text"] = ["elder",]
    columns_categories["sex_text"] = ["female", "male",]
    #columns_categories["sex_text"] = ["female",]
    #columns_categories["sex_text"] = ["male",]
    #columns_categories["visit_text"] = ["first", "second",]
    columns_categories["visit_text"] = ["first",]
    #columns_categories["visit_text"] = ["second",]
    table_category = porg.filter_select_table_rows_by_columns_categories(
        table=table_source,
        columns_categories=columns_categories,
        report=report,
    )

    ##########
    # Organize selections of features and observations.
    # Extract names or identifiers of all columns and rows in the table.
    columns_all = copy.deepcopy(table_category.columns.to_list())
    columns_all.remove(index_rows)
    rows_selection = copy.deepcopy(
        table_category[index_rows].unique().tolist()
    )
    # Filter sets of features for those that are in the table.
    columns_quantitative_available = list(filter(
        lambda feature: (feature in columns_all),
        columns_quantitative
    ))
    columns_olink_plasma_available = list(filter(
        lambda feature: (feature in columns_all),
        columns_olink_plasma
    ))
    columns_olink_muscle_available = list(filter(
        lambda feature: (feature in columns_all),
        columns_olink_muscle
    ))
    columns_olink_adipose_available = list(filter(
        lambda feature: (feature in columns_all),
        columns_olink_adipose
    ))
    # Combine sets of features.
    columns_selection = list()
    columns_selection.extend(columns_quantitative_available)
    columns_selection.extend(columns_olink_plasma_available)
    columns_selection.extend(columns_olink_muscle_available)
    columns_selection.extend(columns_olink_adipose_available)

    ##########
    # Filter table's rows and columns by proportion of nonmissing signal
    # intensities.
    table_filter = (
        porg.filter_table_rows_columns_by_proportion_nonmissing_threshold(
            table=table_category,
            index_columns=index_columns,
            index_rows=index_rows,
            columns_selection=columns_selection,
            rows_selection=rows_selection,
            threshold_low=None,
            threshold_high=None,
            proportion_columns=0.9,
            proportion_rows=0.9,
            report=report,
        )
    )

    ##########
    # Organize selections of features and observations.
    # Extract names or identifiers of all columns and rows in the table that
    # remain after filters.
    columns_all = copy.deepcopy(table_filter.columns.to_list())
    columns_all.remove(index_rows)
    rows_selection = copy.deepcopy(
        table_filter[index_rows].unique().tolist()
    )
    # Copy information.
    columns_quantitative_filter = copy.deepcopy(columns_quantitative_available)
    columns_olink_plasma_filter = copy.deepcopy(columns_olink_plasma_available)
    columns_olink_muscle_filter = copy.deepcopy(columns_olink_muscle_available)
    columns_olink_adipose_filter = copy.deepcopy(
        columns_olink_adipose_available
    )
    # Filter sets of features for those that remain in the table after filters.
    columns_quantitative_filter = list(filter(
        lambda feature: (feature in columns_all),
        columns_quantitative_filter
    ))
    columns_olink_plasma_filter = list(filter(
        lambda feature: (feature in columns_all),
        columns_olink_plasma_filter
    ))
    columns_olink_muscle_filter = list(filter(
        lambda feature: (feature in columns_all),
        columns_olink_muscle_filter
    ))
    columns_olink_adipose_filter = list(filter(
        lambda feature: (feature in columns_all),
        columns_olink_adipose_filter
    ))
    # Combine sets of features.
    columns_selection = list()
    columns_selection.extend(columns_quantitative_filter)
    columns_selection.extend(columns_olink_plasma_filter)
    columns_selection.extend(columns_olink_muscle_filter)
    columns_selection.extend(columns_olink_adipose_filter)

    ##########
    # Fill missing values of signal intensities.
    # Notice that this fill procedure occurs after the filters on proportions
    # of missing signal intensities across features and observations. Hence,
    # those previous filters regulate the extent of fills on missing values.
    table_fill = porg.fill_missing_values_table_by_column(
        table=table_filter,
        index_columns=index_columns,
        index_rows=index_rows,
        columns_selection=columns_selection,
        rows_selection=rows_selection,
        method="mean",
        report=report,
    )

    # Collect information.
    pail = dict()
    pail["table_filter"] = table_filter
    pail["table_fill"] = table_fill
    pail["columns_quantitative"] = columns_quantitative_filter
    pail["columns_olink_plasma"] = columns_olink_plasma_filter
    pail["columns_olink_muscle"] = columns_olink_muscle_filter
    pail["columns_olink_adipose"] = columns_olink_adipose_filter

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_olink.py")
        print("function: filter_fill_table_subject_olink()")
        putly.print_terminal_partition(level=5)
        print("table after filter: ")
        print(pail["table_filter"].iloc[0:10, 0:])
        print(pail["table_filter"])
        putly.print_terminal_partition(level=5)
        print("table after filter and fill: ")
        print(pail["table_fill"].iloc[0:10, 0:])
        print(pail["table_fill"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 4. Calculate principal components on features for Olink measurements.


def drive_manage_calculate_principal_components(
    table=None,
    index_rows=None,
    sets_columns=None,
    report=None,
):
    """
    Organizes the calculation and integration of principal components across
    measurements of proteomics by O-Link technology in specific tissues.

    arguments:
        table (object): Pandas data-frame table of information about samples
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        sets_columns (dict<list<str>>): names of sets of features and names of
            columns corresponding to each feature in each set
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    sets_columns = copy.deepcopy(sets_columns)

    # Collect information.
    pail_collection_total = dict()
    pail_component_columns = dict()
    # Iterate on sets of features.
    for name_set in sets_columns.keys():
        pail_reduction = (
            pdecomp.calculate_principal_components_table_columns_selection(
                table=table,
                index_rows=index_rows,
                columns_selection=sets_columns[name_set],
                prefix=name_set,
                separator="_",
                threshold_proportion=0.005, # float or None to keep all
                report=report,
        ))
        # Copy information in table.
        table = pail_reduction["table"].copy(deep=True)
        # Collect information.
        pail_collection_total[name_set] = pail_reduction
        pail_component_columns[name_set] = (
            pail_reduction["columns_component_scores"]
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("subpackage: proteomics")
        print("module: organize_olink.py")
        print("function: drive_manage_calculate_principal_components()")
        putly.print_terminal_partition(level=5)
        print("table with principal components: ")
        #print(table.iloc[0:10, 0:])
        print(table)
        putly.print_terminal_partition(level=5)
        print("columns for principal components corresponding to each set: ")
        print(pail_component_columns)
        putly.print_terminal_partition(level=5)
        pass
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["sets_columns"] = pail_component_columns
    # Return information.
    return pail


##########
# 5. Prepare derivative, deliverable, product tables for signals.
# Tables represent features for individual Olink targets with measurement
# signal intensities across observations for separate individual subjects.


def prepare_derivative_deliverable_product_tables_signal(
    table=None,
    index_columns=None,
    index_rows=None,
    columns_quantitative=None,
    columns_olink_plasma=None,
    columns_olink_muscle=None,
    columns_olink_adipose=None,
    report=None,
):
    """
    Prepare derivative, deliverable, product table 1. Table 1 represents
    features for individual measurement signal intensities of Olink targets
    across observations for separate individual subjects.

    Format of source table (name: "table_source")
    Format of source table is in wide format with floating-point values of
    signal intensities for measurements of individual Olink targets
    corresponding to features across columns and distinct individual
    observations across rows. A special header row gives identifiers or names
    corresponding to each feature across columns, and a special column gives
    identifiers or names corresponding to each observation across rows. For
    versatility, this table does not have explicitly defined indices across
    rows or columns.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    Review: 12 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        columns_quantitative (list<str>): names of columns in original source
            table for a selection of features on a quantitative, continuous,
            interval or ratio scale of measurement
        columns_olink_plasma (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in plasma
        columns_olink_muscle (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in muscle
        columns_olink_adipose (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in adipose
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    ##########
    # Copy information.
    table_source = table.copy(deep=True)
    columns_quantitative = copy.deepcopy(columns_quantitative)
    columns_olink_plasma = copy.deepcopy(columns_olink_plasma)
    columns_olink_muscle = copy.deepcopy(columns_olink_muscle)
    columns_olink_adipose = copy.deepcopy(columns_olink_adipose)

    ##########
    # Prepare parameters for the preparation of derivative tables.
    names_groups_features_sequence = list()
    names_groups_features_sequence.append("olink_plasma")
    names_groups_features_sequence.append("olink_muscle")
    names_groups_features_sequence.append("olink_adipose")

    names_groups_observations_sequence = list()
    names_groups_observations_sequence.append("younger_female")
    names_groups_observations_sequence.append("younger_male")
    names_groups_observations_sequence.append("elder_female")
    names_groups_observations_sequence.append("elder_male")

    groups_features = dict()
    #groups_features["other"] = columns_quantitative
    groups_features["olink_plasma"] = columns_olink_plasma
    groups_features["olink_muscle"] = columns_olink_muscle
    groups_features["olink_adipose"] = columns_olink_adipose

    features_selection = list()
    for group_features in groups_features.keys():
        features_selection.extend(groups_features[group_features])

    groups_observations = dict()
    groups_observations["younger_female"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier="subject_visit",
            name="younger_female", # or "name_instance"
            columns_categories={
                "visit_text": ["first",],
                "age_cohort_text": ["younger",],
                "sex_text": ["female",],
            },
            report=report,
    ))
    groups_observations["younger_male"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier="subject_visit",
            name="younger_male", # or "name_instance"
            columns_categories={
                "visit_text": ["first",],
                "age_cohort_text": ["younger",],
                "sex_text": ["male",],
            },
            report=report,
    ))
    groups_observations["elder_female"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier="subject_visit",
            name="elder_female", # or "name_instance"
            columns_categories={
                "visit_text": ["first",],
                "age_cohort_text": ["elder",],
                "sex_text": ["female",],
            },
            report=report,
    ))
    groups_observations["elder_male"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier="subject_visit",
            name="elder_male", # or "name_instance"
            columns_categories={
                "visit_text": ["first",],
                "age_cohort_text": ["elder",],
                "sex_text": ["male",],
            },
            report=report,
    ))
    observations_selection = list()
    for group_observations in groups_observations.keys():
        observations_selection.extend(groups_observations[group_observations])

    ##########
    # Prepare tables for signals.
    pail_signal = (
        aexph_sub.prepare_tables_signals_features_sets_observations_groups(
            table=table_source,
            transpose_source_table=False,
            index_features=index_columns,
            index_observations=index_rows,
            features_selection=features_selection,
            observations_selection=observations_selection,
            groups_features=groups_features,
            groups_observations=groups_observations,
            names_groups_features_sequence=(
                names_groups_features_sequence
            ),
            names_groups_observations_sequence=(
                names_groups_observations_sequence
            ),
            translations_features=None,
            translations_observations=None,
            report=report,
        )
    )

    ##########
    # Prepare tables for allocation of features to sets.
    # Notice that "groups_features" here is the same collection of groups or
    # sets of features for use in clustering the table's columns within those
    # groups. Alternatively, it would be possible to define a separate
    # collection of groups or sets of features for use in the allocation.
    table_allocation_1 = (
        aexph_sub.prepare_table_features_sets_allocation_match_table_signal(
            table_signal=pail_signal["table_3"],
            index_features=index_columns,
            indices_observations=[index_rows, "group",],
            groups_features=groups_features,
            names_groups_features_sequence=(
                names_groups_features_sequence
            ),
            translations_features=None,
            report=report,
        )
    )
    table_allocation_2 = (
        aexph_sub.prepare_table_features_sets_allocation_match_table_signal(
            table_signal=pail_signal["table_4"],
            index_features=index_columns,
            indices_observations=[index_rows, "group",],
            groups_features=groups_features,
            names_groups_features_sequence=(
                names_groups_features_sequence
            ),
            translations_features=None,
            report=report,
        )
    )

    ##########
    # Collect information.
    pail_return = dict()
    pail_return["table_signal_1"] = pail_signal["table_3"]
    pail_return["table_signal_2"] = pail_signal["table_4"]
    pail_return["table_feature_1"] = table_allocation_1
    pail_return["table_feature_2"] = table_allocation_2
    pail_return["table_mean"] = pail_signal["table_7"]

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_olink.py")
        function = str(
            "prepare_derivative_deliverable_product_tables_signal()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_return


##########
# 6. Prepare derivative, deliverable, product tables for correlations between
# features across observations.
# Tables represent correlations between features for individual Olink targets
# with measurement signal intensities across observations for separate
# individual subjects.


def prepare_derivative_deliverable_product_tables_correlation(
    table=None,
    index_columns=None,
    index_rows=None,
    columns_quantitative=None,
    columns_olink_plasma=None,
    columns_olink_muscle=None,
    columns_olink_adipose=None,
    report=None,
):
    """
    Prepare derivative, deliverable, product table 1. Table 1 represents
    features for individual measurement signal intensities of Olink targets
    across observations for separate individual subjects.

    Format of source table (name: "table_source")
    Format of source table is in wide format with floating-point values of
    signal intensities for measurements of individual Olink targets
    corresponding to features across columns and distinct individual
    observations across rows. A special header row gives identifiers or names
    corresponding to each feature across columns, and a special column gives
    identifiers or names corresponding to each observation across rows. For
    versatility, this table does not have explicitly defined indices across
    rows or columns.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    Review: 16 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        columns_quantitative (list<str>): names of columns in original source
            table for a selection of features on a quantitative, continuous,
            interval or ratio scale of measurement
        columns_olink_plasma (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in plasma
        columns_olink_muscle (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in muscle
        columns_olink_adipose (list<str>): names of columns in original source
            table for features corresponding to measurements by Olink
            technology of peptides in adipose
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    ##########
    # Copy information.
    table_source = table.copy(deep=True)
    columns_quantitative = copy.deepcopy(columns_quantitative)
    columns_olink_plasma = copy.deepcopy(columns_olink_plasma)
    columns_olink_muscle = copy.deepcopy(columns_olink_muscle)
    columns_olink_adipose = copy.deepcopy(columns_olink_adipose)

    ##########
    # Prepare parameters for the preparation of derivative tables.
    columns_olink_all = list()
    columns_olink_all.extend(columns_olink_plasma)
    columns_olink_all.extend(columns_olink_muscle)
    columns_olink_all.extend(columns_olink_adipose)

    names_groups_features_sequence = list()
    names_groups_features_sequence.append("olink_plasma")
    names_groups_features_sequence.append("olink_muscle")
    names_groups_features_sequence.append("olink_adipose")

    groups_features = dict()
    #groups_features["other"] = columns_quantitative
    groups_features["olink_plasma"] = columns_olink_plasma
    groups_features["olink_muscle"] = columns_olink_muscle
    groups_features["olink_adipose"] = columns_olink_adipose

    features_selection = list()
    for group_features in groups_features.keys():
        features_selection.extend(groups_features[group_features])
    features_selection.extend(columns_quantitative)

    observations_selection = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier="subject_visit",
            name="observations", # or "name_instance"
            columns_categories={
                "visit_text": ["first",],
                "age_cohort_text": ["younger", "elder",],
                "sex_text": ["female", "male",],
            },
            report=report,
    ))

    ##########
    # Prepare tables for correlations.
    pail_correlation = (
        aexph_sub.prepare_tables_correlations_of_features_across_observations(
            table=table_source,
            index_features=index_columns,
            index_observations=index_rows,
            features_selection=features_selection,
            observations_selection=observations_selection,
            groups_features=groups_features,
            names_groups_features_sequence=(
                names_groups_features_sequence
            ),
            features_primary=columns_quantitative,
            features_secondary=columns_olink_all,
            translations_features=None,
            method_priority="spearman", # "pearson", "spearman", "kendall",
            report=None,
        )
    )

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!table_3")
    print(pail_correlation["table_3"])

    if False:

        ##########
        # Prepare tables for allocation of features to sets.
        # Notice that "groups_features" here is the same collection of groups or
        # sets of features for use in clustering the table's columns within those
        # groups. Alternatively, it would be possible to define a separate
        # collection of groups or sets of features for use in the allocation.
        table_allocation_1 = (
            aexph_sub.prepare_table_features_sets_allocation_match_table_signal(
                table_signal=pail_correlation["table_2"],
                index_columns=index_columns,
                indices_observations=[index_rows, "group",],
                groups_features=groups_features,
                names_groups_features_sequence=(
                    names_groups_features_sequence
                ),
                translations_features=None,
                report=report,
            )
        )


    if False:

        table_allocation_2 = (
            aexph_sub.prepare_table_features_sets_allocation_match_table_signal(
                table_signal=pail_signal["table_4"],
                index_features=index_columns,
                indices_observations=[index_rows, "group",],
                groups_features=groups_features,
                names_groups_features_sequence=(
                    names_groups_features_sequence
                ),
                translations_features=None,
                report=report,
            )
        )
        pass

    ##########
    # Collect information.
    pail_return = dict()
    pail_return["table_correlation_1"] = pail_correlation["table_1"]
    pail_return["table_correlation_2"] = pail_correlation["table_2"]
    pail_return["table_correlation_3"] = pail_correlation["table_3"]
    #pail_return["table_feature_1"] = table_allocation_1
    #pail_return["table_feature_2"] = table_allocation_2

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_olink.py")
        function = str(
            "prepare_derivative_deliverable_product_tables_correlation()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_return

##########
# _. Plot charts


# Histograms for Olink


def create_plot_chart_histogram(
    table=None,
    column_feature=None,
    report=None,
):
    """
    Create and plot a chart of the histogram type.

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement
        column_feature (str): name of column for values corresponding to a
            specific feature
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)

    # Extract nonmissing values from column in table.
    # partner.organization.extract_filter_array_values_from_series()
    values_raw = table[column_feature].to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    values_nonmissing = numpy.copy(values_raw[~numpy.isnan(values_raw)])
    mean_values = round(numpy.nanmean(values_nonmissing), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_distribution_histogram(
        array=values_nonmissing,
        title="",
        bin_method="count",
        bin_count=20,
        bar_width=0.5,
        label_bins="Value Distribution",
        label_counts="Counts",
        fonts=fonts,
        colors=colors,
        line=True,
        line_position=mean_values,
        label_title=column_feature,
        label_report=True,
    )

    # Return information.
    return figure


def drive_manage_plot_write_histograms_charts(
    table=None,
    sets_columns=None,
    paths=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        table (object): Pandas data-frame table of values of features across
            columns and sample observations across rows
        sets_columns (dict<list<str>>): names of sets of features and names of
            columns corresponding to each feature in each set
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of figure objects from MatPlotLib

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    sets_columns = copy.deepcopy(sets_columns)

    # Iterate on sets of features.
    for name_set in sets_columns.keys():
        # Initialize child directory in which to write charts to file.
        path_directory = os.path.join(
            paths["out_procedure_plot"], "histogram", str(name_set),
        )
        putly.create_directories(
            path=path_directory,
        )
        # Iterate on features in set.
        for column_feature in sets_columns[name_set]:
            # Create histogram for values of feature across observations.
            figure = create_plot_chart_histogram(
                table=table,
                column_feature=column_feature,
                report=report,
            )
            # Write figure to file.
            pplot.write_product_plot_figure(
                figure=figure,
                format="jpg", # jpg, png, svg
                resolution=300,
                name_file=column_feature,
                path_directory=path_directory,
            )
            pass
        pass
    pass


# Scatters for principal components


def define_response_features_principal_components_olink():
    """
    Defines names of columns in a table for features of interest as response,
    outcome, or dependent variables.

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Define records for features of interest and parameters for their
    # handling.
    records = list()
    if True:
        records.append({
            "name": "sex_text",
            "type": "category",
        })
        records.append({
            "name": "age_cohort_text",
            "type": "category",
        })
    if True:
        records.append({
            "name": "body_mass_index",
            "type": "continuity",
        })
        records.append({
            "name": "insulin_sensitivity",
            "type": "continuity",
        })
    # Return information.
    return records


def create_plot_chart_scatter_point_response(
    table=None,
    column_response=None,
    column_abscissa=None,
    column_ordinate=None,
    type_response=None,
    report=None,
):
    """
    Create and plot a chart of the scatter point type.

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement
        column_response (str): name of column for values corresponding to a
            specific feature
        column_abscissa (str): name of column with values for representation on
            abscissa horizontal axis
        column_ordinate (str): name of column with values for representation on
            ordinate vertical axis
        type_response (bool): type of the response feature, either 'continuity'
            or 'category'
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_scatter(
        data=table,
        abscissa=column_abscissa,
        ordinate=column_ordinate,
        title_abscissa=column_abscissa,
        title_ordinate=column_ordinate,
        fonts=fonts,
        colors=colors,
        size=10,
    )
    figure = pplot.plot_scatter_point_color_response_discrete_or_continuous(
        table=table,
        column_identifier="subject_visit",
        column_name="subject_visit",
        column_response=column_response,
        column_abscissa=column_abscissa,
        column_ordinate=column_ordinate,
        type_response=type_response,
        minimum_abscissa=None,
        maximum_abscissa=None,
        minimum_ordinate=None,
        maximum_ordinate=None,
        set_axis_limits=False, # helpful if using minima and maxima to filter values
        title_response=column_response,
        title_abscissa=column_abscissa,
        title_ordinate=column_ordinate,
        identifiers_emphasis=None,
        emphasis_label=None,
        line_diagonal=True, # diagonal is not proportional to respective ranges of axes
        size_title_abscissa="nine",
        size_title_ordinate="nine",
        size_title_legend_bar="ten",
        size_label_abscissa="thirteen",
        size_label_ordinate="thirteen",
        size_label_emphasis="thirteen",
        size_label_legend_bar="thirteen",
        aspect="landscape",
        fonts=fonts,
        colors=colors,
        report=None,
    )
    # Return information.
    return figure


def drive_manage_plot_write_principal_component_scatter_charts(
    table=None,
    features_response=None,
    sets_columns=None,
    paths=None,
    report=None,
):
    """
    Plot scatter point chart representations of values of principal components
    on selections of features across observations.

    arguments:
        table (object): Pandas data-frame table of values of features across
            columns and sample observations across rows
        features_response (list<dict<str>>): names of columns and type
            description corresponding to features of interest as response,
            outcome, or dependent variables
        sets_columns (dict<list<str>>): names of sets of principal components
            corresponding to specific features and names of columns for
            each principal component in each set
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of figure objects from MatPlotLib

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    features_response = copy.deepcopy(features_response)
    sets_columns = copy.deepcopy(sets_columns)

    # Iterate on selection of columns for extra features.
    for record in features_response:
        # Extract information.
        name_feature = record["name"]
        type_feature = record["type"]
        # Iterate on sets of principal components corresponding to features.
        for name_set in sets_columns.keys():
            # Copy information in table.
            table_excerpt = table.copy(deep=True)
            # Separate information in table for a specific selection of columns.
            columns_excerpt = copy.deepcopy(sets_columns[name_set])
            columns_excerpt.insert(0, name_feature)
            # Filter columns in table.
            #table_excerpt = table_excerpt.loc[
            #    :, table_excerpt.columns.isin(columns_excerpt)
            #].copy(deep=True)
            table_excerpt = table_excerpt.filter(
                items=columns_excerpt,
                axis="columns",
            )
            # Filter rows in table for non-missing values across relevant columns.
            table_excerpt.dropna(
                axis="index",
                how="any",
                subset=columns_excerpt,
                inplace=True,
            )
            # Initialize child directory in which to write charts to file.
            path_directory = os.path.join(
                paths["out_procedure_plot"], "scatter_pca",
                str(name_feature), str(name_set),
            )
            putly.create_directories(
                path=path_directory,
            )
            # Define comparisons between principal components.
            comparisons = [
                ["_1", "_2",],
                ["_1", "_3",],
                ["_1", "_4",],
                ["_2", "_3",],
                ["_2", "_4",],
                ["_3", "_4",],
            ]
            # Iterate on comparisons.
            for comparison in comparisons:
                # Determine name for chart.
                name_chart = str(
                    name_feature + "_" + name_set +
                    "_pc" + comparison[0] + comparison[1]
                )
                # Plot chart of type scatter point as representation of principal
                # components.
                figure = create_plot_chart_scatter_point_response(
                    table=table,
                    column_response=name_feature,
                    column_abscissa=str(name_set + comparison[0]),
                    column_ordinate=str(name_set + comparison[1]),
                    type_response=type_feature,
                    report=report,
                )
                # Write figure to file.
                pplot.write_product_plot_figure(
                    figure=figure,
                    format="jpg", # jpg, png, svg
                    resolution=300,
                    name_file=name_chart,
                    path_directory=path_directory,
                )
                pass
            pass
        pass
    pass


# Heatmaps.


def plot_write_heatmap_chart_signal(
    table_signal=None,
    table_feature=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    name_file=None,
    paths=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        table_signal (object): Pandas data-frame table of floating-point values
            of a signal corresponding features in columns across observations
            in rows
        table_feature (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table of signals
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        name_file (str): name for file
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of figure objects from MatPlotLib

    """

    # Copy information.
    table_signal = table_signal.copy(deep=True)
    table_feature = table_feature.copy(deep=True)

    # Initialize child directory in which to write charts to file.
    path_directory = os.path.join(
        paths["out_procedure_plot"], "heatmap_signal",
    )
    putly.create_directories(
        path=path_directory,
    )
    # Create heatmap.
    figure = (
        aexph_sub.plot_heatmap_signal_features_sets_observations_groups(
            table_signal=table_signal,
            table_feature=table_feature,
            index_columns=index_columns,
            index_rows=index_rows,
            column_group="group",
            report=report,
    ))
    # Write figure to file.
    pplot.write_product_plot_figure(
        figure=figure,
        format="jpg", # jpg, png, svg
        resolution=300,
        name_file=name_file,
        path_directory=path_directory,
    )
    pass


def plot_write_heatmap_chart_mean(
    table_mean=None,
    index_columns=None,
    index_rows=None,
    name_file=None,
    paths=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        table_mean (object): Pandas data-frame table of floating-point values
            of mean signals
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        name_file (str): name for file
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of figure objects from MatPlotLib

    """

    # Copy information.
    table_mean = table_mean.copy(deep=True)

    # Initialize child directory in which to write charts to file.
    path_directory = os.path.join(
        paths["out_procedure_plot"], "heatmap_signal",
    )
    putly.create_directories(
        path=path_directory,
    )
    # Create heatmap.
    figure = (
        aexph_sub.plot_heatmap_signal_mean(
            table=table_mean,
            index_columns=index_columns,
            index_rows=index_rows,
            report=False,
    ))
    # Write figure to file.
    pplot.write_product_plot_figure(
        figure=figure,
        format="jpg", # jpg, png, svg
        resolution=300,
        name_file=name_file,
        path_directory=path_directory,
    )
    pass




###############################################################################
# Procedure


# TODO: TCW; 6 December 2024
# 1. calculate correlations between...
#   - subject continuous traits and Olink signals
#   - pairwise Olink signals
#   - subject continuous traits and PCs on Olink signals
#   -

# IDEA: TCW; 6 December 2024
# - implement a main parameter table for this procedure as in "transcriptomics.compare_sets_groups"
# - it would not be reasonable to compare PCs side-by-side for different groups of observations
# - hence filters on observations should be more or less central before calculating the PCs

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
    procedure="organize_olink"
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

    ##########
    # 3. Filter and fill information in table for analysis.
    pail_filter_fill = filter_fill_table_subject_olink(
        table=pail_source["table_subject"],
        index_columns="features",
        index_rows="subject_visit",
        columns_quantitative=pail_source["columns_quantitative"],
        columns_olink_plasma=pail_source["columns_olink_plasma"],
        columns_olink_muscle=pail_source["columns_olink_muscle"],
        columns_olink_adipose=pail_source["columns_olink_adipose"],
        fill_missing=True,
        report=report,
    )
    #pail_filter_fill["table_filter"]
    #pail_filter_fill["table_fill"]
    #pail_filter_fill["columns_quantitative"]
    #pail_filter_fill["columns_olink_plasma"]
    #pail_filter_fill["columns_olink_muscle"]
    #pail_filter_fill["columns_olink_adipose"]

    ##########
    # 4. Transform signals for Olink targets to different scale.
    if False:
        # Transform signals to their original, non-logarithmic scale and then
        # back to logarithmic scale at a different base.
        columns_transform = list()
        columns_transform.extend(
            copy.deepcopy(pail_filter_fill["columns_olink_plasma"])
        )
        columns_transform.extend(
            copy.deepcopy(pail_filter_fill["columns_olink_muscle"])
        )
        columns_transform.extend(
            copy.deepcopy(pail_filter_fill["columns_olink_adipose"])
        )
        table_transform = pscl.transform_exponent_by_table_columns(
            table=pail_filter_fill["table_fill"],
            columns=columns_transform,
            base=2.0,
            report=report,
        )
        table_transform = pscl.transform_logarithm_by_table_columns(
            table=table_transform,
            columns=columns_transform,
            base=10.0,
            report=report,
        )
    else:
        # Copy information.
        table_transform = pail_filter_fill["table_fill"].copy(deep=True)

    ##########
    # 5. Calculate principal components on features for Olink measurements.
    # Organize sets of features.
    columns_olink_plasma = copy.deepcopy(
        pail_filter_fill["columns_olink_plasma"]
    )
    columns_olink_muscle = copy.deepcopy(
        pail_filter_fill["columns_olink_muscle"]
    )
    columns_olink_adipose = copy.deepcopy(
        pail_filter_fill["columns_olink_adipose"]
    )
    columns_olink_all = list()
    columns_olink_all.extend(columns_olink_plasma)
    columns_olink_all.extend(columns_olink_muscle)
    columns_olink_all.extend(columns_olink_adipose)
    pail_columns_decomposition = dict()
    pail_columns_decomposition["olink_all"] = columns_olink_all
    pail_columns_decomposition["olink_plasma"] = columns_olink_plasma
    pail_columns_decomposition["olink_muscle"] = columns_olink_muscle
    pail_columns_decomposition["olink_adipose"] = columns_olink_adipose
    # Calculate principal components.
    pail_reduction = drive_manage_calculate_principal_components(
            table=table_transform,
            index_rows="subject_visit",
            sets_columns=pail_columns_decomposition,
            report=report,
    )
    if False:
        pail_pca_test = (
            pdecomp.calculate_principal_components_table_columns_selection(
                table=table_transform,
                index_rows="subject_visit",
                columns_selection=columns_olink_plasma,
                prefix="olink_plasma_test",
                separator="_",
                threshold_proportion=0.02, # float or None to keep all
                report=report,
        ))
    # Notice that the names of columns corresponding to principal components
    # of features in each original set are in ascending sort order.
    #pail_reduction["sets_columns"]


    ##########
    # 6. Prepare derivative, deliverable product tables for actual signals.
    # This set of tables represents features and sets of features for
    # individual Olink targets with measurement signal intensities across
    # observations and groups of observations for separate individual subjects.
    pail_signal = prepare_derivative_deliverable_product_tables_signal(
        table=pail_reduction["table"],
        index_columns="features",
        index_rows="subject_visit",
        columns_quantitative=pail_filter_fill["columns_quantitative"],
        columns_olink_plasma=pail_filter_fill["columns_olink_plasma"],
        columns_olink_muscle=pail_filter_fill["columns_olink_muscle"],
        columns_olink_adipose=pail_filter_fill["columns_olink_adipose"],
        report=report,
    )
    if False:
        # Continue with scores of principal components.
        pail_signal = prepare_derivative_deliverable_product_tables_signal(
            table=pail_reduction["table"],
            index_columns="features",
            index_rows="subject_visit",
            columns_quantitative=list(),
            columns_olink_plasma=pail_reduction["sets_columns"]["olink_plasma"],
            columns_olink_muscle=pail_reduction["sets_columns"]["olink_muscle"],
            columns_olink_adipose=pail_reduction["sets_columns"]["olink_adipose"],
            report=report,
        )
        pass

    ##########
    # Note about eventual parameter table
    # 1. information about selection of observations
    # 2. information about sets of features for correlations (primary and secondary)
    # 3. sets or groups of features for both primary and secondary contexts
    #    - reference library of feature sets in private parameters directory
    # 4. for primary features on horizontal axis of heatmap, whether to
    #    - represent a few features by explicit name
    #      or
    #    - represent many features by sets, similar to vertical axis

    # TODO: TCW; 17 December 2024
    # Integrate optional "features_groups" for both primary and secondary sets
    # of features in the correlation procedure.
    # Keep 1 single "translations_features"
    # Clean up accross-the-board translation of feature names before proceeding with procedure




    ##########
    # 7. Prepare derivative, deliverable product tables for correlations
    # between features across observations.
    pail_correlation = (
        prepare_derivative_deliverable_product_tables_correlation(
            table=pail_reduction["table"],
            index_columns="features",
            index_rows="subject_visit",
            columns_quantitative=pail_filter_fill["columns_quantitative"],
            columns_olink_plasma=pail_filter_fill["columns_olink_plasma"],
            columns_olink_muscle=pail_filter_fill["columns_olink_muscle"],
            columns_olink_adipose=pail_filter_fill["columns_olink_adipose"],
            report=report,
    ))
    # table_signal=pail_correlation["table_correlation_2"],
    # table_feature=pail_correlation["table_feature"],

    #pail_correlation["table_correlation_1"]
    #pail_correlation["table_correlation_2"]
    #pail_correlation["table_correlation_3"]


    # Correlations...
    # Olink plasma, muscle, and adipose with Olink plasma, muscle, and adipose
    # Olink plasma, muscle, and adipose with selection of continuous clinical variables (BMI, SI, etc)
    # Olink PC's with Olink PC's
    # Olink PC's with selection of continuous clinical variables (BMI, SI, etc)

    # TODO: Plot correlations
    # need separate functions to plot correlations either in symmetrical heatmaps
    # or in assymmetrical heatmaps.
    # Symmetrical heatmaps have sets of features on both horizontal and vertical axes
    # Assymmetrical heatmaps have sets of features only on vertical axis and names on horizontal axis


    ###########################################################3

    ##########
    # Collect information.
    # Collections of files.

    pail_write_lists = dict()
    pail_write_lists["olink_all"] = columns_olink_all
    pail_write_lists["olink_plasma"] = columns_olink_plasma
    pail_write_lists["olink_muscle"] = columns_olink_muscle
    pail_write_lists["olink_adipose"] = columns_olink_adipose
    pail_write_tables = dict()
    pail_write_tables[str("table_subject")] = pail_reduction["table"]
    pail_write_tables[str("table_correlation_1")] = pail_correlation["table_correlation_1"]
    pail_write_tables[str("table_correlation_2")] = pail_correlation["table_correlation_2"]
    pail_write_tables[str("table_correlation_3")] = pail_correlation["table_correlation_3"]
    pail_write_objects = dict()

    #pail_write_objects[str("samples")]
    ##########
    # Write product information to file.
    putly.write_lists_to_file_text(
        pail_write=pail_write_lists,
        path_directory=paths["out_procedure_data_lists"],
        delimiter="\n",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_procedure_data_tables"],
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_procedure_data_tables"],
        reset_index_rows=None,
        write_index_rows=None,
        write_index_columns=None,
        type="pickle",
        delimiter=None,
        suffix=".pickle",
    )
    if False:
        putly.write_objects_to_file_pickle(
            pail_write=pail_write_objects,
            path_directory=paths["out_procedure_data"],
        )


    ##########
    # Plot histogram charts for features on quantitative continuous interval
    # or ratio scales of measurement.
    if False:
        pail_columns_histogram = dict()
        pail_columns_histogram["other"] = copy.deepcopy(
            pail_parse["columns_quantitative"]
        )
        pail_columns_histogram["olink_plasma"] = copy.deepcopy(
            pail_parse["columns_olink_plasma"]
        )
        pail_columns_histogram["olink_muscle"] = copy.deepcopy(
            pail_parse["columns_olink_muscle"]
        )
        pail_columns_histogram["olink_adipose"] = copy.deepcopy(
            pail_parse["columns_olink_adipose"]
        )
        drive_manage_plot_write_histograms_charts(
            table=table_subject,
            sets_columns=pail_columns_histogram,
            paths=paths,
            report=report,
        )
        pass

    if False:

        ##########
        # Plot heatmaps of signals.
        plot_write_heatmap_chart_signal(
            table_signal=pail_signal["table_signal_1"],
            table_feature=pail_signal["table_feature_1"],
            index_columns="features",
            index_rows="subject_visit",
            column_group="group",
            name_file="signal_1",
            paths=paths,
            report=report,
        )
        plot_write_heatmap_chart_signal(
            table_signal=pail_signal["table_signal_2"],
            table_feature=pail_signal["table_feature_2"],
            index_columns="features",
            index_rows="subject_visit",
            column_group="group",
            name_file="signal_2",
            paths=paths,
            report=report,
        )
        plot_write_heatmap_chart_mean(
            table_mean=pail_signal["table_mean"],
            index_columns="groups",
            index_rows="feature",
            name_file="mean",
            paths=paths,
            report=report,
        )

    ##########
    # Plot scatter charts to represent principal components on sets of features
    # for Olink measurements.
    if False:
        features_response = define_response_features_principal_components_olink()
        drive_manage_plot_write_principal_component_scatter_charts(
            table=pail_reduction["table"],
            features_response=features_response,
            sets_columns=pail_reduction["sets_columns"],
            paths=paths,
            report=report,
        )

    pass


###############################################################################
# End
