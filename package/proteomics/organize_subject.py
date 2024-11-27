"""
Supply functionality for process and analysis of data from proteomics.

This module 'organize_subject' is part of the 'proteomics' package within
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
import exercise.transcriptomics.organize_sample as extr_sample

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
        print("module: exercise.proteomics.organize_subject.py")
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


def define_type_columns_table_subject_feature_organization():
    """
    Defines the variable types of columns within table for organization of
    attributes of samples.

    Review: TCW; 22 August 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion_transcriptomics"] = "int32"
    types_columns["inclusion_proteomics"] = "int32"
    types_columns["category_raw"] = "string"
    types_columns["category"] = "string"
    types_columns["name_source"] = "string"
    types_columns["name_intermediate"] = "string"
    types_columns["name_product"] = "string"
    types_columns["type"] = "string"
    types_columns["description"] = "string"
    # ...
    # Return information.
    return types_columns


def parse_extract_table_sample_feature_organization(
    table=None,
    inclusion=None,
    report=None,
):
    """
    Determines the variable types of columns within table for attributes of
    samples.

    Review: TCW; 22 August 2024

    arguments:
        table (object): Pandas data-frame table of information about samples
        inclusion (str): name of column in table to use as a filter on the
            rows from table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Filter rtable_feature_organizationn of relevant samples.
    table_inclusion = table.loc[
        (table[inclusion] == 1), :
    ].copy(deep=True)

    # Extract information for types of columns.
    table_types_columns = table_inclusion.filter(
        items=["name_source", "type",],
        axis="columns",
    )
    series_types_columns = pandas.Series(
        table_types_columns["type"].to_list(),
        index=table_types_columns["name_source"],
    )
    types_columns = series_types_columns.to_dict()

    # Extract information for translation of names of columns.
    table_name_product = table_inclusion.loc[
        (table_inclusion["name_product"].str.len() > 0), :
    ].copy(deep=True)
    table_translations = table_name_product.filter(
        items=["name_intermediate", "name_product",],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_translations["name_product"].to_list(),
        index=table_translations["name_intermediate"],
    )
    translations_column = series_translations.to_dict()

    # Extract product names of all columns.
    columns_all = copy.deepcopy(
        table_inclusion["name_product"].to_list()
    )

    # Extract names of columns corresponding to O-Link measurements in each
    # type of tissue.
    columns_olink_plasma = copy.deepcopy(table_inclusion.loc[
        (table_inclusion["category"] == "olink_plasma"), :
    ].copy(deep=True)["name_product"].to_list())
    columns_olink_muscle = copy.deepcopy(table_inclusion.loc[
        (table_inclusion["category"] == "olink_muscle"), :
    ].copy(deep=True)["name_product"].to_list())
    columns_olink_adipose = copy.deepcopy(table_inclusion.loc[
        (table_inclusion["category"] == "olink_adipose"), :
    ].copy(deep=True)["name_product"].to_list())

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_sample.py")
        print("function: parse_extract_table_sample_feature_organization()")
        putly.print_terminal_partition(level=5)
        print("types for columns upon read: ")
        print(types_columns)
        putly.print_terminal_partition(level=5)
        print(
            "translations for columns after initial, intermediate " +
            "translation:"
        )
        print(translations_column)
        putly.print_terminal_partition(level=5)
        print("Counts of each category of columns")
        print(table_inclusion["category"].value_counts(dropna=False))
        putly.print_terminal_partition(level=5)
        print(
            "Counts of all columns."
        )
        count_all = len(columns_all)
        print(str("all: " + str(count_all)))
        putly.print_terminal_partition(level=5)
        print(
            "Counts of columns corresponding to O-Link measurements in each " +
            "type of tissue."
        )
        count_plasma = len(columns_olink_plasma)
        count_muscle = len(columns_olink_muscle)
        count_adipose = len(columns_olink_adipose)
        print(str("plasma: " + str(count_plasma)))
        print(str("muscle: " + str(count_muscle)))
        print(str("adipose: " + str(count_adipose)))
        putly.print_terminal_partition(level=5)
        pass
    # Collect information.
    pail = dict()
    pail["types_columns"] = types_columns
    pail["translations_column"] = translations_column
    pail["columns_all"] = columns_all
    pail["columns_olink_plasma"] = columns_olink_plasma
    pail["columns_olink_muscle"] = columns_olink_muscle
    pail["columns_olink_adipose"] = columns_olink_adipose
    # Return information.
    return pail


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
    path_file_table_subject_property = os.path.join(
        paths["in_data"], "study_age_exercise", "subject_sample",
        "table_subject_sample_feature_olink_hek_2024-11-22.csv",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of parameters for organization of the table of attributes for
    # subjects and samples.
    types_columns = define_type_columns_table_subject_feature_organization()
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
    #    inclusion="inclusion_proteomics",
    #    report=report,
    #)

    # Table of attributes for samples.
    pail["table_subject_property"] = pandas.read_csv(
        path_file_table_subject_property,
        sep=",",
        header=0,
        dtype=pail_parse["types_columns"],
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Fill information about intervention in experimental condition.
    #pail["table_subject_property"]["Intervention"] = "Placebo"

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.proteomics.organize_subject.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for subjects and samples: ")
        print(pail["table_subject_property"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Organize table of attributes for samples.


def define_sequence_columns_novel_sample_feature():
    """
    Defines names of columns in sequence by which to filter and sort columns in
    a table.

    This list represents the columns that are novel derivations of the original
    columns.

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "match_sample_attribute_file_transcriptomics",
        "cohort_age",
        "cohort_age_text",
        #"cohort_age_letter",
        "intervention",
        "intervention_text",
        #"identifier_subject_attribute",
        #"study_clinic_visit_relative",
        "date_visit_text",
        #"date_visit_text_raw",
        "sex_y",
        #"sex_letter",
        "sex_text",
        #"age",
        #"body_mass_index",
        #"body_fat_percent",
        #"body_fat_mass",
        #"body_lean_mass",
        #"oxygen_consumption",
        #"tertiles_body_mass_index",
        #"tertiles_body_skeletal_muscle_index",
        #"tertiles_body_fat_percent",
        #"tertiles_insulin_sensitivity",
        #"tertiles_activity_steps",
    ]
    # Return information.
    return columns_sequence


def determine_match_sample_file_forward(
    subject=None,
    study_clinic_visit=None,
):
    """
    Determines a designator to match samples from their files of signals with
    their attributes.

    arguments:
        subject (str): identifier of study participant subject
        study_clinic_visit (str): indicator of clinical visit in the study at
            which collection of a sample occurred, either 'Pre' or 'Post'

    raises:

    returns:
        (str): common designator to match samples from their files of signals
            to their attributes

    """

    # Determine designator.
    if (
        (pandas.notna(subject)) and
        (len(str(subject).strip()) > 0) and
        (pandas.notna(study_clinic_visit)) and
        (len(str(study_clinic_visit).strip()) > 0)
    ):
        # There is adequate information.
        subject = str(subject).strip()
        study_clinic_visit = str(study_clinic_visit).strip().lower()
        if (study_clinic_visit == "pre"):
            visit = str("first")
            designator = str(subject + "_" + visit)
        elif (study_clinic_visit == "post"):
            visit = str("second")
            designator = str(subject + "_" + visit)
        else:
            designator = ""
    else:
        designator = ""
        pass
    # Return information.
    return designator


def determine_cohort_age_text(
    cohort_age_letter=None,
):
    """
    Determines a designator for cohort by age.

    arguments:
        cohort_age_letter (str): designator of cohort by age

    raises:

    returns:
        (str): designator of cohort by age

    """

    # Determine designator.
    if (
        (pandas.notna(cohort_age_letter)) and
        (len(str(cohort_age_letter).strip()) > 0)
    ):
        # There is adequate information.
        cohort_age_letter = str(cohort_age_letter).strip().upper()
        if (cohort_age_letter == "E"):
            designator = str("elder")
        elif (cohort_age_letter == "Y"):
            designator = str("younger")
        else:
            designator = ""
    else:
        designator = ""
        pass
    # Return information.
    return designator


def determine_cohort_age(
    cohort_age_text=None,
):
    """
    Determines a designator for cohort by age.

    arguments:
        cohort_age_text (str): designator of cohort by age

    raises:

    returns:
        (float): designator of cohort by age

    """

    # Determine designator.
    if (
        (pandas.notna(cohort_age_text)) and
        (len(str(cohort_age_text).strip()) > 0)
    ):
        # There is adequate information.
        cohort_age_text = str(cohort_age_text).strip().lower()
        if (cohort_age_text == "elder"):
            designator = 1
        elif (cohort_age_text == "younger"):
            designator = 0
        else:
            designator = float("nan")
    else:
        designator = float("nan")
        pass
    # Return information.
    return designator


def determine_intervention(
    intervention_text=None,
):
    """
    Determines a designator for intervention.

    arguments:
        intervention_text (str): designator of intervention

    raises:

    returns:
        (float): designator of intervention

    """

    # Determine designator.
    if (
        (pandas.notna(intervention_text)) and
        (len(str(intervention_text).strip()) > 0)
    ):
        # There is adequate information.
        intervention_text = str(intervention_text).strip().lower()
        if (intervention_text == "active"):
            designator = 1
        elif (intervention_text == "placebo"):
            designator = 0
        else:
            designator = float("nan")
    else:
        designator = float("nan")
        pass
    # Return information.
    return designator


def determine_sex_text(
    sex_letter=None,
):
    """
    Determines a designator for sex.

    arguments:
        sex_letter (str): designator of sex

    raises:

    returns:
        (str): designator of sex

    """

    # Determine designator.
    if (
        (pandas.notna(sex_letter)) and
        (len(str(sex_letter).strip()) > 0)
    ):
        # There is adequate information.
        sex_letter = str(sex_letter).strip().upper()
        if (sex_letter == "F"):
            designator = str("female")
        elif (sex_letter == "M"):
            designator = str("male")
        else:
            designator = ""
    else:
        designator = ""
        pass
    # Return information.
    return designator


def determine_sex_y(
    sex_text=None,
):
    """
    Determines a designator for sex.

    arguments:
        sex_text (str): designator of sex

    raises:

    returns:
        (float): designator of sex

    """

    # Determine designator.
    if (
        (pandas.notna(sex_text)) and
        (len(str(sex_text).strip()) > 0)
    ):
        # There is adequate information.
        sex_text = str(sex_text).strip().lower()
        if (sex_text == "female"):
            designator = 0
        elif (sex_text == "male"):
            designator = 1
        else:
            designator = float("nan")
    else:
        designator = float("nan")
        pass
    # Return information.
    return designator


def determine_date_visit_text(
    date_visit_text_raw=None,
):
    """
    Determines a designator for date of visit to the clinic for study.

    arguments:
        date_visit_text (str): designator of date

    raises:

    returns:
        (object): designator of date

    """

    # Determine designator.
    if (
        (pandas.notna(date_visit_text_raw)) and
        (len(str(date_visit_text_raw).strip()) > 0)
    ):
        # There is adequate information.
        date_visit_text_raw = str(date_visit_text_raw).strip()
        date_visit = datetime.strptime(date_visit_text_raw, '%m/%d/%Y').date()
        designator = date_visit.strftime('%Y-%m-%d')
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

    # Translate names of columns to remove white space.
    #table.columns = [column.replace("\n", "_") for column in table.columns]
    table.columns = table.columns.str.strip()
    table.columns = table.columns.str.replace(" ", "_")
    table.columns = table.columns.str.replace("\n", "_")

    # Translate names of columns.
    table.rename(
        columns=translations_column,
        inplace=True,
    )

    # Filter rows in table.
    table = table.loc[
        (
            (table["identifier_subject_attribute"].str.len() > 0)
        ), :
    ].copy(deep=True)
    table.dropna(
        how="all",
        axis="index",
    )

    # Determine designation to match sample to attribute.
    table["match_sample_attribute_file_transcriptomics"] = table.apply(
        lambda row:
            determine_match_sample_file_forward(
                subject=row["identifier_subject_attribute"],
                study_clinic_visit=row["study_clinic_visit_relative"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine designations of cohort by age.
    table["cohort_age_text"] = table.apply(
        lambda row:
            determine_cohort_age_text(
                cohort_age_letter=row["cohort_age_letter"],
            ),
        axis="columns", # apply function to each row
    )
    table["cohort_age"] = table.apply(
        lambda row:
            determine_cohort_age(
                cohort_age_text=row["cohort_age_text"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine designations of intervention.
    table["intervention_text"] = table.apply(
        lambda row: str(row["intervention_text_raw"]).strip().lower(),
        axis="columns", # apply function to each row
    )
    table["intervention"] = table.apply(
        lambda row:
            determine_intervention(
                intervention_text=row["intervention_text"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine designations of sex.
    table["sex_text"] = table.apply(
        lambda row:
            determine_sex_text(
                sex_letter=row["sex_letter"],
            ),
        axis="columns", # apply function to each row
    )
    table["sex_y"] = table.apply(
        lambda row:
            determine_sex_y(
                sex_text=row["sex_text"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine date of visit to the clinic for study.
    table["date_visit_text"] = table.apply(
        lambda row:
            determine_date_visit_text(
                date_visit_text_raw=row["date_visit_text_raw"],
            ),
        axis="columns", # apply function to each row
    )
    # Replace values of zero for oxygen consumption with missing values.
    # Use the values of maximal oxygen consumption adjusted to lean body mass.
    table["oxygen_consumption"] = table["oxygen_consumption"].replace(
        to_replace=0,
        value=pandas.NA,
    )

    # Clean values for C-reactive protein variable.
    table["c_react_protein"] = table.apply(
        lambda row:
            str(row["c_react_protein"]).replace(
                "<.2", "0.1"
            ).replace("<0.2", "0.1"),
        axis="columns", # apply function to each row
    )
    table["c_react_protein"] = table["c_react_protein"].astype("float32")

    # Sort rows within table.
    table.sort_values(
        by=[
            "cohort_age",
            "intervention",
            "identifier_subject_attribute",
            "study_clinic_visit_relative",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Filter and sort columns within table.
    #columns_sequence.insert(0, column_index)
    columns_sequence = copy.deepcopy(columns_original)
    columns_sequence.extend(columns_novel)
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
        print("module: exercise.proteomics.organize_sample_olink.py")
        print("function: organize_table_subject_property()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table"].iloc[0:10, 0:])
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail





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
    columns_sequence["cohort_age_text_by_sex_text"] = "elder_by_male"
    columns_sequence["cohort_age_text_by_exercise_time_point"] = (
        "elder_by_3_hour"
    )
    columns_sequence["sex_text_by_exercise_time_point"] = "male_by_3_hour"
    columns_sequence["intervention_text_by_study_clinic_visit"] = (
        "active_by_second"
    )
    columns_sequence["sex_text_by_study_clinic_visit"] = "male_by_second"
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
        print("module: exercise.proteomics.organize_sample_olink.py")
        print("function: organize_table_sample_interaction_combinations()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table"].iloc[0:10, 0:])
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# TODO: TCW; 2024-08-27
# Ideally, the calculation of the PCs should happen AFTER stratifying the
# samples for the cohort that's relevant to the analysis.


##########
# 4. Organize proteomics from measurements by O-Link technology.


def organize_olink_principal_components_tissue(
    table=None,
    column_index=None,
    columns_olink=None,
    prefix=None,
    report=None,
):
    """
    Organizes the calculation and integration of principal components across
    measurements of proteomics by O-Link technology in a specific tissue.

    arguments:
        table (object): Pandas data-frame table of information about samples
        column_index (str): name of column for index of unique values across
            rows
        columns_olink (list<str>): names of columns corresponding to O-Link
            measurements in a single tissue
        prefix (str): prefix for names of columns corresponding to scores for
            principal components
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_main = table.copy(deep=True)
    table_excerpt = table.copy(deep=True)
    # Copy other information.
    columns_olink = copy.deepcopy(columns_olink)
    # Organize indices in table.
    table_main.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_main.set_index(
        [column_index],
        append=False,
        drop=True,
        inplace=True,
    )
    table_excerpt.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Separate information in table for proteomics from measurements by O-Link
    # technology.
    columns_excerpt = copy.deepcopy(columns_olink)
    columns_excerpt.insert(0, column_index)
    #table_excerpt = table_excerpt.loc[
    #    :, table_excerpt.columns.isin(columns_excerpt)
    #].copy(deep=True)
    table_excerpt = table_excerpt.filter(
        items=columns_excerpt,
        axis="columns",
    )
    table_excerpt.set_index(
        [column_index],
        append=False,
        drop=True,
        inplace=True,
    )

    # Calculate principal components to represent variance with a reduction of
    # dimensionality.
    pail_reduction = (
        pdecomp.organize_principal_components_by_singular_value_decomposition(
            table=table_excerpt,
            index_name=column_index,
            prefix=prefix,
            separator="_",
            report=False, # report is extensive
        )
    )
    table_component_scores = (
        pail_reduction["table_component_scores"].copy(deep=True)
    )
    # Organize indices in table.
    table_component_scores.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_component_scores.set_index(
        column_index,
        append=False,
        drop=True,
        inplace=True
    )
    # Extract names of columns corresponding to scores for principal
    # components.
    columns_component_scores = (
        copy.deepcopy(table_component_scores.columns.to_list())
    )
    # Merge scores for principal components with the main table.
    table_merge = porg.merge_columns_two_tables(
        identifier_first=column_index,
        identifier_second=column_index,
        table_first=table_main,
        table_second=table_component_scores,
        preserve_index=True,
        report=report,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.proteomics.organize_sample_olink.py")
        print("function: organize_olink_principal_components_tissue()")
        putly.print_terminal_partition(level=5)
        print("table of excerpt information for decomposition: ")
        print(table_excerpt)
        putly.print_terminal_partition(level=5)
        print("table of scores for principal components: ")
        print(table_component_scores.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("columns of scores for principal components: ")
        count_columns_component_scores = len(columns_component_scores)
        print(columns_component_scores)
        print(str(
            "count of scores for principal components: " +
            str(count_columns_component_scores)
        ))
        putly.print_terminal_partition(level=5)
        # Compare methods for calculation of Pincipal Component Analysis.
        pdecomp.compare_principal_components_methods(
            table=table_excerpt,
            index_name=column_index,
            prefix=prefix,
            separator="_",
            report=True,
        )
        pass
    # Collect information.
    pail = dict()
    pail["table"] = table_merge
    pail["columns_component_scores"] = columns_component_scores
    # Return information.
    return pail


def organize_olink_principal_components_tissues(
    table=None,
    column_index=None,
    columns_olink_plasma=None,
    columns_olink_muscle=None,
    columns_olink_adipose=None,
    report=None,
):
    """
    Organizes the calculation and integration of principal components across
    measurements of proteomics by O-Link technology in specific tissues.

    arguments:
        table (object): Pandas data-frame table of information about samples
        column_index (str): name of column for index of unique values across
            rows
        columns_olink_plasma (list<str>): names of columns corresponding to
            O-Link measurements in plasma tissue
        columns_olink_muscle (list<str>): names of columns corresponding to
            O-Link measurements in muscle tissue
        columns_olink_adipose (list<str>): names of columns corresponding to
            O-Link measurements in adipose tissue
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    columns_olink_plasma = copy.deepcopy(columns_olink_plasma)
    columns_olink_muscle = copy.deepcopy(columns_olink_muscle)
    columns_olink_adipose = copy.deepcopy(columns_olink_adipose)


    pail_plasma = organize_olink_principal_components_tissue(
        table=table,
        column_index=column_index,
        columns_olink=columns_olink_plasma,
        prefix="olink_plasma_component",
        report=report,
    )


    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.proteomics.organize_sample_olink.py")
        print("function: organize_olink_principal_components_tissues()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(table.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    #return table
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
    project="exercise"
    routine="proteomics"
    procedure="organize_subject"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.proteomics.organize_subject.py")
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

    # source["table_subject_property"]

    ##########
    # 3. Parse table of parameters describing properties for subjects.
    pail_parse = parse_extract_table_sample_feature_organization(
        table=pail_source["table_feature_organization"],
        inclusion="inclusion_proteomics",
        report=report,
    )

    ##########
    # 4. Organize table of attributes for samples.
    columns_original = pail_parse["columns_all"]
    columns_novel = define_sequence_columns_novel_sample_feature() # <-- this might move to "transcriptomics"
    pail_organization = organize_table_subject_property(
        table=pail_source["table_subject_property"],
        translations_column=pail_parse["translations_column"],
        columns_original=columns_original,
        columns_novel=columns_novel,
        report=report,
    )
    table_sample_attribute = pail_sample_attribute["table"]



    ##########
    # 3. Organize table of properties for subjects.
    columns_interaction = define_interaction_combination_categorical_factor()
    columns_novel = define_sequence_columns_novel_sample_feature()
    table_subject_property = organize_table_subject_property(
        table=pail_source["table_subject_property"],
        translations_column=pail_parse["translations_column"],
        columns_original=pail_parse["columns_all"],
        columns_novel=columns_novel,
        report=report,
    )

    ##########
    # 4. Organize proteomics from measurements by O-Link technology.
    table_sample_olink_components = (
        organize_olink_principal_components_tissues(
            table=table_subject_property,
            column_index="identifier_subject_attribute",
            columns_olink_plasma=pail_parse["columns_olink_plasma"],
            columns_olink_muscle=pail_parse["columns_olink_muscle"],
            columns_olink_adipose=pail_parse["columns_olink_adipose"],
            report=report,
    ))

    # 4.1. extract from the parameter table the columns of O-Link targets in
    # each tissue



    # "table_subject.tsv"
    # "table_subject.pickle"



    pass


###############################################################################
# End
