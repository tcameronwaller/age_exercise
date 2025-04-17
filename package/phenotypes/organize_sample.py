"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'organize_sample' is part of the 'phenotypes' subpackage within
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
import partner.plot as pplot
import partner.parallelization as prall
import age_exercise.phenotypes.organize_subject as aexph_sub

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.



##########
# 2. Read source information from file.


def define_type_columns_table_sample_file():
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
    types_columns["identifier_signal"] = "string"
    types_columns["path_file"] = "string"
    types_columns["sample_plate"] = "string"
    types_columns["plate"] = "string"
    types_columns["identifier_sample"] = "string"
    types_columns["identifier_subject"] = "string"
    types_columns["tissue"] = "string"
    types_columns["condition_code"] = "string"
    types_columns["condition_correction"] = "string"
    types_columns["condition_interpretation"] = "string"
    types_columns["note_condition"] = "string"
    #types_columns["sex"] = "string"
    #types_columns["age"] = "int32"
    #types_columns["body_mass"] = "float32"
    #types_columns["body_mass_index"] = "float32"
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

    # Define paths to child files.
    path_file_table_feature_organization = os.path.join(
        paths["in_data"], "study_age_exercise", "subject_sample",
        "table_subject_sample_feature_organization.tsv",
    )
    path_file_table_subject = os.path.join(
        paths["out_project"], "phenotypes", "organize_subject", "tables",
        "table_subject.pickle",
    )
    path_file_table_sample_file = os.path.join(
        paths["in_data"], "study_age_exercise", "subject_sample",
        "table_sample_file_rnaseq.tsv",
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
    #pail_parse = parse_extract_table_sample_feature_organization(
    #    table=pail["table_feature_organization"],
    #    inclusion="inclusion",
    #    report=report,
    #)

    # Table of properties for subjects.
    pail["table_subject"] = pandas.read_pickle(
        path_file_table_subject,
    )

    # Table of matches between samples and files.
    types_columns = define_type_columns_table_sample_file()
    pail["table_sample_file"] = pandas.read_csv(
        path_file_table_sample_file,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of properties for subjects: ")
        print(pail["table_subject"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("table of matches between samples and files: ")
        print(pail["table_sample_file"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Organize table of properties for study subjects.


def define_translation_columns_table_subject_property():
    """
    Defines translations for the names of columns in a table.

    arguments:

    raises:

    returns:
        (dict<str>): translations for names of columns in a table

    """


    # Translate names of columns.
    translations = dict()
    translations["identifier_subject"] = "identifier_subject_study"
    translations["visit_text"] = "visit_text_subject"
    # Return information.
    return translations


def determine_match_subject_sample_file_forward(
    subject=None,
    visit_text=None,
):
    """
    Determines a designator to match samples from their files of signals with
    their attributes.

    arguments:
        subject (str): identifier of study participant subject
        visit_text (str): indicator of clinical visit in the study at
            which collection of a sample occurred, either 'first' or 'second'

    raises:

    returns:
        (str): common designator to match samples from their files of signals
            to their attributes

    """

    # Determine designator.
    if (
        (pandas.notna(subject)) and
        (len(str(subject).strip()) > 0) and
        (pandas.notna(visit_text)) and
        (len(str(visit_text).strip()) > 0)
    ):
        # There is adequate information.
        subject = str(subject).strip()
        visit_text = str(visit_text).strip().lower()
        designator = str(subject + "_" + visit_text)
    else:
        designator = ""
        pass
    # Return information.
    return designator


def organize_table_subject_property(
    table=None,
    translations_column=None,
    columns_original=None,
    columns_novel=None,
    report=None,
):
    """
    Organizes information in table that provides attributes of samples.

    This function prepares the table of sample attributes for merge with the
    table of matches between samples and files.

    arguments:
        table (object): Pandas data-frame table of subjects, samples, and their
            attribute features
        translations_column (dict<str>): translations for names of columns in a
            table
        columns_original (list<str>): names of original columns in sequence by
            which to filter and sort columns in table
        columns_novel (list<str>): names of original columns in sequence by
            which to filter and sort columns in table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    translations_column = copy.deepcopy(translations_column)
    columns_original = copy.deepcopy(columns_original)
    columns_novel = copy.deepcopy(columns_novel)

    # Translate names of columns.
    table.rename(
        columns=translations_column,
        inplace=True,
    )

    # Filter rows in table.
    table = table.loc[
        (
            (table["identifier_subject_study"].str.len() > 0)
        ), :
    ].copy(deep=True)
    table.dropna(
        how="all",
        axis="index",
    )

    # Determine designation to match sample to attribute.
    table["match_subject_sample_file_transcriptomics"] = table.apply(
        lambda row:
            determine_match_subject_sample_file_forward(
                subject=row["identifier_subject_study"],
                visit_text=row["visit_text_subject"],
            ),
        axis="columns", # apply function to each row
    )

    # Sort rows within table.
    table.sort_values(
        by=[
            "age_cohort_elder",
            "intervention_omega3",
            "identifier_subject_study",
            "visit_text_subject",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Filter and sort columns within table.
    #columns_sequence.insert(0, column_index)
    columns_sequence = copy.deepcopy(columns_novel)
    columns_sequence.extend(columns_original)
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Collect information.
    pail = dict()
    pail["columns_sequence"] = columns_sequence
    pail["table"] = table

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: organize_table_subject_property()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table"].iloc[0:10, 0:])
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 4. Organize table of matches between samples and files.


def define_translation_columns_table_sample_file():
    """
    Defines translations for the names of columns in a table.

    arguments:

    raises:

    returns:
        (dict<str>): translations for names of columns in a table

    """


    # Translate names of columns.
    translations = dict()
    #translations["identifier"] = "identifier_signal"
    translations["condition_interpretation"] = "condition_obsolete"
    # Return information.
    return translations


def define_sequence_columns_table_sample_file():
    """
    Defines names of columns in sequence by which to filter and sort columns in
    a table.

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "inclusion",
        "identifier_subject",
        "identifier_sample",
        "identifier_signal",
        "tissue",
        "condition_correction",
        "visit_text",
        "exercise_time_point",
        "match_subject_sample_file_transcriptomics",
        #"path_file",
        #"sample_plate",
        #"plate",
        #"condition_code",
        #"condition_interpretation",
        #"note_condition",
    ]
    # Return information.
    return columns_sequence


def determine_sample_visit_text(
    tissue=None,
    instance=None,
):
    """
    Determines the clinical visit of the study at which collection of a
    sample occurred.

    arguments:
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        instance (str): designation of study instance in terms of clinical
            visit for sample collection

    raises:

    returns:
        (str): indicator of clinical visit in the study at which collection of
            a sample occurred, either 'first' or 'second'

    """

    # Determine indicator.
    if (
        (pandas.notna(tissue)) and
        (len(str(tissue).strip()) > 0) and
        (pandas.notna(instance)) and
        (len(str(instance).strip()) > 0)
    ):
        # There is adequate information.
        if (
            (str(tissue).strip().lower() == "muscle") and
            (str(instance).strip() in ["1B", "2B", "3B"])
        ):
            indicator = "first"
        elif (
            (str(tissue).strip().lower() == "adipose") and
            (str(instance).strip() == "B")
        ):
            indicator = "first"
        elif (
            (str(tissue).strip().lower() == "adipose") and
            (str(instance).strip() == "PI")
        ):
            indicator = "second"
        else:
            indicator = ""
    else:
        indicator = ""
        pass
    # Return information.
    return indicator


def determine_muscle_exercise_time_point(
    tissue=None,
    instance=None,
):
    """
    Determines the approximate categorical or ordinal duration of time after
    exercise at which sample collection of muscle tissue occurred.

    arguments:
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        instance (str): designation of study instance in terms of clinical
            visit for sample collection

    raises:

    returns:
        (str): approximate categorical or ordinal duration of time after
            exercise, either '0_hour', '3_hour', or '48_hour'

    """

    # Determine indicator.
    if (
        (pandas.notna(tissue)) and
        (len(str(tissue).strip()) > 0) and
        (pandas.notna(instance)) and
        (len(str(instance).strip()) > 0)
    ):
        # There is adequate information.
        if (
            (str(tissue).strip().lower() == "muscle") and
            (str(instance).strip() == "1B")
        ):
            indicator = "0_hour"
        elif (
            (str(tissue).strip().lower() == "muscle") and
            (str(instance).strip() == "2B")
        ):
            indicator = "3_hour"
        elif (
            (str(tissue).strip().lower() == "muscle") and
            (str(instance).strip() == "3B")
        ):
            indicator = "48_hour"
        else:
            indicator = ""
    else:
        indicator = ""
        pass
    # Return information.
    return indicator


def determine_match_subject_sample_file_reverse(
    subject=None,
    visit_text=None,
):
    """
    Determines a designator to match samples from their files of signals with
    their attributes.

    arguments:
        subject (str): identifier of study participant subject
        visit_text (str): indicator of clinical visit in the study at
            which collection of a sample occurred, either 'first' or 'second'

    raises:

    returns:
        (str): common designator to match samples from their files of signals
            to their attributes

    """

    # Determine designator.
    if (
        (pandas.notna(subject)) and
        (len(str(subject).strip()) > 0) and
        (pandas.notna(visit_text)) and
        (len(str(visit_text).strip()) > 0)
    ):
        # There is adequate information.
        subject = str(subject).strip()
        visit_text = str(visit_text).strip()
        designator = str(subject + "_" + visit_text)
    else:
        designator = ""
        pass
    # Return information.
    return designator


def organize_table_sample_file(
    table=None,
    translations_column=None,
    columns_sequence=None,
    report=None,
):
    """
    Organizes information in table that designates matches between samples and
    their corresponding files of data.

    This function prepares the table of matches between samples and files for
    merge with table of attributes for samples.

    arguments:
        table (object): Pandas data-frame table of information about samples
        translations_column (dict<str>): translations for names of columns in a
            table
        columns_sequence (list<str>): names of columns in sequence by which to
            filter and sort columns in table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    translations_column = copy.deepcopy(translations_column)
    columns_sequence = copy.deepcopy(columns_sequence)

    # Translate names of columns.
    table.rename(
        columns=translations_column,
        inplace=True,
    )
    # Sort rows within table.
    table.sort_values(
        by=[
            "tissue",
            "identifier_subject",
            "condition_correction",
        ],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Determine the clinical visit of the study at which collection of the
    # sample occurred.
    table["visit_text"] = table.apply(
        lambda row:
            determine_sample_visit_text(
                tissue=row["tissue"],
                instance=row["condition_correction"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine designation to match sample to attribute.
    table["match_subject_sample_file_transcriptomics"] = table.apply(
        lambda row:
            determine_match_subject_sample_file_reverse(
                subject=row["identifier_subject"],
                visit_text=row["visit_text"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine designation of time point from the study of exercise in muscle.
    table["exercise_time_point"] = table.apply(
        lambda row:
            determine_muscle_exercise_time_point(
                tissue=row["tissue"],
                instance=row["condition_correction"],
            ),
        axis="columns", # apply function to each row
    )
    # Filter and sort columns within table.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=report,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: organize_table_sample_file()")
        putly.print_terminal_partition(level=5)
        print("table of matches between samples and files: ")
        print(table.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# 5. Combine within the same table the matches between samples and files
# along with their further attributes.


def combine_table_subject_sample_file_property(
    table_sample_file=None,
    table_subject=None,
    columns_transfer=None,
    report=None,
):
    """
    Combines in the same table information about samples.

    arguments:
        table_sample_file (object): Pandas data-frame table of information
            about samples at the level of individual files of data
        table_subject (object): Pandas data-frame table of information
            about samples at the level of individual subjects and their
            clinical visits for the study
        columns_transfer (list<str>): names of columns for attributes to
            transfer
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_sample_file = table_sample_file.copy(deep=True)
    table_subject = table_subject.copy(deep=True)

    # Transfer attributes.
    table = porg.transfer_table_rows_attributes_reference(
        table_main=table_sample_file,
        column_main_key="match_subject_sample_file_transcriptomics",
        table_reference=table_subject,
        column_reference_key="match_subject_sample_file_transcriptomics",
        columns_reference_transfer=columns_transfer,
        prefix_reference_main="",
        suffix_reference_main="",
        report=report,
    )
    # Copy information in table.
    table = table.copy(deep=True)

    # Sort rows within table.
    table.sort_values(
        by=[
            "tissue",
            "age_cohort_elder",
            "intervention_omega3",
            "sex_y",
            "identifier_subject",
            "visit_text",
            "exercise_time_point",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: combine_table_subject_sample_file_property()")
        putly.print_terminal_partition(level=5)
        print("table of files and attributes for samples: ")
        print(table.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table



##########
# 6. Prepare combinations of categorical factor variables for analyses of
# interaction.


def define_interaction_combination_categorical_factor():
    """
    Defines names of columns for interaction combinations of categorical factor
    variables.

    arguments:

    raises:

    returns:
        (dict<str>): names of columns for interaction combinations of
            categorical factor variables and their specific single values for
            designation as not 'other'

    """

    # Specify sequence of columns within table.
    columns_sequence = dict()
    columns_sequence["age_cohort_text_by_sex_text"] = "elder_by_male"
    columns_sequence["age_cohort_text_by_exercise_time_point"] = (
        "elder_by_3_hour"
    )
    columns_sequence["sex_text_by_exercise_time_point"] = "male_by_3_hour"
    columns_sequence["intervention_text_by_visit_text"] = (
        "omega3_by_second"
    )
    columns_sequence["sex_text_by_visit_text"] = "male_by_second"

    columns_sequence["age_cohort_text_by_intervention_text"] = "elder_by_omega3" # other: elder_by_placebo, younger_by_none
    # Return information.
    return columns_sequence


# TODO: TCW; 18 October 2024
# This function could additionally prepare binary "dummies" for categories
# and then calculate the product combinations for interaction effects.


def organize_table_sample_interaction_combinations(
    table=None,
    columns_interaction=None,
    report=None,
):
    """
    Organizes information in table that provides attributes of samples.

    This function prepares combinations of categorical features for analysis of
    interaction effects.

    arguments:
        table (object): Pandas data-frame table of subjects, samples, and their
            attribute features
        columns_interaction (dict<str>): names of columns for interaction
            combinations of categorical factor variables and their specific
            single values for designation as not 'other'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    columns_interaction = copy.deepcopy(columns_interaction)

    # Define interaction combinations of categorical factor variables.
    #interactions = define_interaction_combination_categorical_factor()
    for interaction in columns_interaction.keys():
        parts = list(str(interaction).strip().split("_by_"))
        part_first = parts[0]
        part_second = parts[1]
        # Create interaction combination of categorical factor variables.
        table[interaction] = table.apply(
            lambda row: "_by_".join([
                str(row[part_first]),
                str(row[part_second]),
            ]),
            axis="columns", # apply function to each row
        )
        # Desginate all values of interaction combination as 'other' that do
        # not match the specific single effect value.
        match = str(columns_interaction[interaction]).strip()
        table[interaction] = table.apply(
            lambda row:
                row[interaction] if (row[interaction] == match) else "other",
            axis="columns", # apply function to each row
        )
        pass

    # Filter and sort columns within table.
    # TODO: TCW; 18 October 2024
    # Consider whether to filter columns again at this point.

    # Collect information.
    pail = dict()
    pail["columns_interaction"] = columns_interaction
    pail["table"] = table

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: organize_table_sample_interaction_combinations()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table"].iloc[0:10, 0:])
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 7. Describe factors in table of samples.


def describe_table_sample_factors(
    table_sample=None,
    report=None,
):
    """
    Describes factors in table of information about samples.

    arguments:
        table_sample (object): Pandas data-frame table of information about
            samples that correspond to signals within accompanying main table
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)

    # Filter rows in table for selection of relevant samples.
    table_inclusion = table_sample.loc[
        (table_sample["inclusion"] == 1), :
    ].copy(deep=True)
    table_elder = table_inclusion.loc[
        (table_inclusion["age_cohort_text"] == "elder"), :
    ].copy(deep=True)
    table_tissue = table_elder.loc[
        (table_elder["tissue"] == "adipose"), :
    ].copy(deep=True)

    # Create cross tabulation.
    cross_tabulation = pandas.crosstab(
        table_tissue["visit_text"].astype("category"),
        table_tissue["intervention_text"].astype("category"),
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: describe_table_sample_factors()")
        putly.print_terminal_partition(level=4)
        print("tissue: adipose")
        print(cross_tabulation)
        putly.print_terminal_partition(level=4)
        pass


    # Return information.
    pass


##########
# 7. Describe sets of samples for specific analyses.


# TODO: TCW; 9 September 2024
# TODO: obsolete
def define_selections_sample_set():
    """
    Defines selection criteria for sets of samples in specific analyses.

    arguments:

    raises:

    returns:
        (list<dict<str>>): names and values of features for selection of
            samples in sets for specific analyses

    """

    # Define instances that determine sets of samples.

    # format for 'name_set': {tissue}_{cohort-details}_{condition-details}

    instances = [
        {
            "name_set": "muscle_elder_exercise-3hr",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "age_cohort_text": ["elder",],
            },
            "factor_availability": {
                "exercise_time_point": ["0_hour", "3_hour",],
            },
        },
        {
            "name_set": "muscle_elder_exercise-48hr",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "age_cohort_text": ["elder",],
            },
            "factor_availability": {
                "exercise_time_point": ["0_hour", "48_hour",],
            },
        },
        {
            "name_set": "muscle_younger_exercise-3hr",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "age_cohort_text": ["younger",],
            },
            "factor_availability": {
                "exercise_time_point": ["0_hour", "3_hour",],
            },
        },
        {
            "name_set": "muscle_younger_exercise-48hr",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "age_cohort_text": ["younger",],
            },
            "factor_availability": {
                "exercise_time_point": ["0_hour", "48_hour",],
            },
        },
        {
            "name_set": "muscle_exercise-0hr_age",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "exercise_time_point": ["0_hour",],
            },
            "factor_availability": {
                "age_cohort_text": ["younger", "elder",],
            },
        },
        {
            "name_set": "muscle_exercise-3hr_age",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "exercise_time_point": ["3_hour",],
            },
            "factor_availability": {
                "age_cohort_text": ["younger", "elder",],
            },
        },
        {
            "name_set": "muscle_exercise-48hr_age",
            "tissue": "muscle",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["muscle",],
                "exercise_time_point": ["48_hour",],
            },
            "factor_availability": {
                "age_cohort_text": ["younger", "elder",],
            },
        },
        {
            "name_set": "adipose_visit-first_age",
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "visit_text": ["first",],
            },
            "factor_availability": {
                "age_cohort_text": ["younger", "elder",],
            },
        },
        {
            "name_set": str(
                "adipose_elder_visit-intervention"
            ),
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "age_cohort_text": ["elder",],
            },
            "factor_availability": {
                "visit_text": ["first", "second",],
                "intervention_text": ["placebo", "omega3",],
            },
        },


        {
            "name_set": str(
                "adipose_elder-visit-first_intervention"
            ),
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "visit_text": ["first",],
            },
            "factor_availability": {
                "intervention_text": ["placebo", "omega3",],
            },
        },
        {
            "name_set": str(
                "adipose_elder-visit-second_intervention"
            ),
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "visit_text": ["second",],
            },
            "factor_availability": {
                "intervention_text": ["placebo", "omega3",],
            },
        },
        {
            "name_set": str(
                "adipose_elder-placebo_visit"
            ),
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "intervention_text": ["placebo",],
            },
            "factor_availability": {
                "visit_text": ["first", "second",],
            },
        },
        {
            "name_set": str(
                "adipose_elder-omega3_visit"
            ),
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "intervention_text": ["omega3",],
            },
            "factor_availability": {
                "visit_text": ["first", "second",],
            },
        },
    ]
    # Return information.
    return instances


# TODO: TCW; 9 September 2024
# TODO: obsolete
def describe_table_sample_sets(
    table_sample=None,
    selections=None,
    report=None,
):
    """
    Describes samples in sets for specific analyses.

    arguments:
        table_sample (object): Pandas data-frame table of information about
            samples that correspond to signals within accompanying main table
        selections (list<dict<str>>): names and values of features for
            selection of samples in sets for specific analyses
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: describe_table_sample_sets()")
        putly.print_terminal_partition(level=5)
        print("count of selections: " + str(len(selections)))
        putly.print_terminal_partition(level=5)
        pass


    # Copy information in table.
    table_sample = table_sample.copy(deep=True)
    # Copy other information.
    selections = copy.deepcopy(selections)

    # Iterate on instances of selection criteria for sets of samples.
    for selection in selections:
        # Report.
        if report:
            putly.print_terminal_partition(level=3)
            print("name_set: " + str(selection["name_set"]))
            print("tissue: " + str(selection["tissue"]))
            putly.print_terminal_partition(level=4)
            pass
        # Iterate on features and values for selection of samples in cohort.
        table_cohort = table_sample.copy(deep=True)
        for feature in selection["cohort_selection"].keys():
            table_cohort = table_cohort.loc[(
                table_cohort[feature].isin(
                    selection["cohort_selection"][feature])
            ), :].copy(deep=True)
        # Iterate on factors and values for selection of samples on the basis
        # of availability.
        table_factor = table_cohort.copy(deep=True)
        for factor in selection["factor_availability"].keys():
            table_factor = table_factor.loc[(
                table_factor[factor].isin(
                    selection["factor_availability"][factor])
            ), :].copy(deep=True)
            # Report.
            if report:
                putly.print_terminal_partition(level=5)
                print("factor: " + str(factor))
                print(
                    "counts of samples with each unique categorical value of "
                    + "factor:")
                print(table_factor[factor].value_counts(dropna=False))
                putly.print_terminal_partition(level=5)
                pass
            pass
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("module: age_exercise.phenotypes.organize_sample.py")
        print("function: describe_table_sample_sets()")
        putly.print_terminal_partition(level=5)
        print("description complete")
        putly.print_terminal_partition(level=5)
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
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

    raises:

    returns:

    """

    ##########
    # Parameters.
    project="age_exercise"
    routine="phenotypes"
    procedure="organize_sample"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: organize_sample.py")
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
    # pail_source["table_subject"]

    ##########
    # 3. Parse table of parameters describing properties for subjects.
    pail_parse = aexph_sub.parse_extract_table_sample_feature_organization(
        table=pail_source["table_feature_organization"],
        inclusion="inclusion",
        report=report,
    )

    ##########
    # 4. Organize table of properties for study subjects.
    translations_subject = define_translation_columns_table_subject_property()
    columns_subject_original = pail_parse["columns_all"]
    columns_subject_original.remove("identifier_subject")
    #columns_subject_original.append("identifier_subject_study")
    columns_subject_original.insert(5, "identifier_subject_study")
    columns_subject_novel = (
        aexph_sub.define_sequence_columns_novel_subject_feature()
    )
    columns_subject_novel.remove("visit_text")
    #columns_subject_novel.append("visit_text_subject")
    columns_subject_novel.insert(5, "visit_text_subject")
    columns_subject_novel.insert(
        5, "match_subject_sample_file_transcriptomics"
    )
    pail_subject = organize_table_subject_property(
        table=pail_source["table_subject"],
        translations_column=translations_subject,
        columns_original=columns_subject_original,
        columns_novel=columns_subject_novel,
        report=report,
    )
    table_subject = pail_subject["table"]

    ##########
    # 5. Organize table of matches between samples and files.
    translations_sample_file = define_translation_columns_table_sample_file()
    columns_sample_file = define_sequence_columns_table_sample_file()
    table_sample_file = organize_table_sample_file(
        table=pail_source["table_sample_file"],
        translations_column=translations_sample_file,
        columns_sequence=columns_sample_file,
        report=report,
    )

    ##########
    # 6. Combine within the same table the matches between samples and files
    # along with their further attributes.
    # The main reason that this combination is necessary and special is that
    # there are records in the table of data for subjects that match multiple
    # samples for RNAseq transcriptomics, especially in the study on muscle.
    columns_transfer = copy.deepcopy(columns_subject_original)
    columns_transfer.extend(columns_subject_novel)
    columns_transfer.remove("match_subject_sample_file_transcriptomics")
    table_sample_merge = combine_table_subject_sample_file_property(
        table_sample_file=table_sample_file,
        table_subject=table_subject,
        columns_transfer=columns_transfer,
        report=report,
    )

    ##########
    # 7. Prepare combinations of categorical factor variables for analyses of
    # interaction.
    columns_interaction = (
        define_interaction_combination_categorical_factor()
    )
    pail_interaction = (
        organize_table_sample_interaction_combinations(
            table=table_sample_merge,
            columns_interaction=columns_interaction,
            report=report,
    ))

    ##########
    # 8. Describe factors in table of samples.
    describe_table_sample_factors(
        table_sample=pail_interaction["table"],
        report=report,
    )

    ##########
    # 7. Describe sets of samples for specific analyses.
    if False:
        selections = define_selections_sample_set()
        describe_table_sample_sets(
            table_sample=pail_interaction["table"],
            selections=selections,
            report=report,
        )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_tables = dict()
    pail_write_tables[str("table_sample")] = pail_interaction["table"]
    pail_write_objects = dict()
    #pail_write_objects[str("samples")]

    ##########
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
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_procedure_tables"],
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
            path_directory=paths["out_data"],
        )

    pass


###############################################################################
# End
