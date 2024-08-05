"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'organize_sample' is part of the 'transcriptomics' package within
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
import partner.plot as pplot
import partner.parallelization as prall

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
    #paths["out_test"] = os.path.join(
    #    paths["out_procedure"], "test",
    #)
    paths["out_data"] = os.path.join(
        paths["out_procedure"], "data",
    )
    #paths["out_plot"] = os.path.join(
    #    paths["out_procedure"], "plot",
    #)
    # Initialize directories in main branch.
    paths_initialization = [
        paths["out_project"],
        paths["out_routine"],
        paths["out_procedure"],
        paths["out_data"],
    ]
    # Remove previous directories and files to avoid version or batch
    # confusion.
    if restore:
        for path in paths_initialization:
            putly.remove_directory(path=path)
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
        print("module: exercise.transcriptomics.organize_sample.py")
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
    types_columns["identifier"] = "string"
    types_columns["path_file"] = "string"
    types_columns["sample_plate"] = "string"
    types_columns["plate"] = "string"
    types_columns["sample"] = "string"
    types_columns["subject"] = "string"
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


def define_type_columns_table_sample_attribute():
    """
    Defines the variable types of columns within table for attributes of
    samples.

    Review: TCW; 30 July 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["Name"] = "string"
    types_columns["Visit"] = "string"
    types_columns["Date"] = "string"
    types_columns["Age Group"] = "string"
    types_columns["Sex"] = "string"
    types_columns["Intervention"] = "string"
    types_columns["Age"] = "int32"
    types_columns["BMI (kg/m2)"] = "float32"
    types_columns["Total % Fat"] = "float32"
    types_columns["Total Fat Mass"] = "float32"
    types_columns["Lean Mass"] = "float32"
    # ...
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
    path_file_table_sample_file = os.path.join(
        paths["in_parameters_private"], "table_sample_file_rnaseq.tsv",
    )
    path_file_table_sample_attribute = os.path.join(
        paths["in_parameters_private"], "table_sample_attribute.csv",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

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

    # Table of attributes for samples.
    types_columns = define_type_columns_table_sample_attribute()
    pail["table_sample_attribute"] = pandas.read_csv(
        path_file_table_sample_attribute,
        sep=",",
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
        print("module: exercise.transcriptomics.organize_sample.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of matches between samples and files: ")
        print(pail["table_sample_file"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table_sample_attribute"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Organize table of matches between samples and files.


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
    translations["identifier"] = "identifier_signal"
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
        "identifier_signal",
        "tissue",
        "sample",
        "subject",
        "condition_correction",
        "study_clinic_visit",
        "exercise_time_point",
        "match_sample_file_attribute",
        #"path_file",
        #"sample_plate",
        #"plate",
        #"condition_code",
        #"condition_interpretation",
        #"note_condition",
    ]
    # Return information.
    return columns_sequence


def determine_sample_study_clinic_visit(
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


def determine_match_sample_file_attribute(
    subject=None,
    study_clinic_visit=None,
):
    """
    Determines a designator to match samples from their files of signals with
    their attributes.

    arguments:
        subject (str): identifier of study participant subject
        study_clinic_visit (str): indicator of clinical visit in the study at
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
        (pandas.notna(study_clinic_visit)) and
        (len(str(study_clinic_visit).strip()) > 0)
    ):
        # There is adequate information.
        subject = str(subject).strip()
        study_clinic_visit = str(study_clinic_visit).strip()
        designator = str(subject + "_" + study_clinic_visit)
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
            "subject",
            "condition_correction",
        ],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Determine the clinical visit of the study at which collection of the
    # sample occurred.
    table["study_clinic_visit"] = table.apply(
        lambda row:
            determine_sample_study_clinic_visit(
                tissue=row["tissue"],
                instance=row["condition_correction"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine designation to match sample to attribute.
    table["match_sample_file_attribute"] = table.apply(
        lambda row:
            determine_match_sample_file_attribute(
                subject=row["subject"],
                study_clinic_visit=row["study_clinic_visit"],
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
        print("module: exercise.transcriptomics.organize_sample.py")
        print("function: organize_table_sample_file()")
        putly.print_terminal_partition(level=5)
        print("table of matches between samples and files: ")
        print(table.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# 4. Organize table of attributes for samples.


def define_translation_columns_table_sample_attribute():
    """
    Defines translations for the names of columns in a table.

    arguments:

    raises:

    returns:
        (dict<str>): translations for names of columns in a table

    """

    # Translate names of columns.
    translations = dict()
    translations["Name"] = "subject_attribute"
    translations["Visit"] = "study_clinic_visit_relative"
    translations["Date"] = "date_visit_text_raw"
    translations["Age Group"] = "cohort_age_letter"
    translations["Sex"] = "sex_letter"
    translations["Intervention"] = "intervention_text"
    translations["Age"] = "age"
    translations["BMI (kg/m2)"] = "body_mass_index"
    translations["Total % Fat"] = "body_fat_percent"
    translations["Total Fat Mass"] = "body_fat_mass"
    translations["Lean Mass"] = "body_lean_mass"
    #translations["Weight (cm)"] = "weight_kg"
    # Return information.
    return translations


def define_sequence_columns_table_sample_attribute():
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
        "match_sample_file_attribute",
        "cohort_age",
        "cohort_age_text",
        "cohort_age_letter",
        "intervention",
        "intervention_text",
        "subject_attribute",
        "study_clinic_visit_relative",
        "date_visit_text",
        "date_visit_text_raw",
        "sex_y",
        "sex_letter",
        "sex_text",
        "age",
        "body_mass_index",
        "body_fat_percent",
        "body_fat_mass",
        "body_lean_mass",
    ]
    # Return information.
    return columns_sequence


def determine_match_sample_attribute_file(
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


def organize_table_sample_attribute(
    table=None,
    translations_column=None,
    columns_sequence=None,
    report=None,
):
    """
    Organizes information in table that provides attributes of samples.

    This function prepares the table of sample attributes for merge with the
    table of matches between samples and files.

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
    # Determine designation to match sample to attribute.
    table["match_sample_file_attribute"] = table.apply(
        lambda row:
            determine_match_sample_attribute_file(
                subject=row["subject_attribute"],
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
        lambda row: str(row["intervention_text"]).strip().lower(),
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

    # Sort rows within table.
    table.sort_values(
        by=[
            "cohort_age",
            "intervention",
            "subject_attribute",
            "study_clinic_visit_relative",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
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
        print("module: exercise.transcriptomics.organize_sample.py")
        print("function: organize_table_sample_attribute()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(table.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# 4. Combine within the same table the matches between samples and files
# along with their further attributes.


# "match_sample_file_attribute"
    table_sample = combine_table_sample_file_attribute(
        table_sample_file=table_sample_file,
        table_sample_attribute=table_sample_attribute,
        columns_transfer=columns_transfer,
        report=True,
    )


##########
# 5. Combine within the same table the matches between samples and files
# along with their further attributes.


def combine_table_sample_file_attribute(
    table_sample_file=None,
    table_sample_attribute=None,
    columns_transfer=None,
    report=None,
):
    """
    Combines in the same table information about samples.

    arguments:
        table_sample_file (object): Pandas data-frame table of information
            about samples at the level of individual files of data
        table_sample_attribute (object): Pandas data-frame table of information
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
    table_sample_attribute = table_sample_attribute.copy(deep=True)

    # Transfer attributes.
    table = porg.transfer_table_rows_attributes_reference(
        table_main=table_sample_file,
        column_main_key="match_sample_file_attribute",
        table_reference=table_sample_attribute,
        column_reference_key="match_sample_file_attribute",
        columns_reference_transfer=columns_transfer,
        prefix_reference_main="",
        suffix_reference_main="",
        report=report,
    )

    # Sort rows within table.
    table.sort_values(
        by=[
            "tissue",
            "cohort_age",
            "intervention",
            "sex_y",
            "subject",
            "study_clinic_visit",
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
        print("module: exercise.transcriptomics.organize_sample.py")
        print("function: combine_table_sample_file_attribute()")
        putly.print_terminal_partition(level=5)
        print("table of files and attributes for samples: ")
        print(table.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# 6. Describe factors in table of samples.


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

    # Filter rows for selection of relevant samples.
    table_inclusion = table_sample.loc[
        (table_sample["inclusion"] == 1), :
    ].copy(deep=True)
    table_elder = table_inclusion.loc[
        (table_inclusion["cohort_age_text"] == "elder"), :
    ].copy(deep=True)
    table_tissue = table_elder.loc[
        (table_elder["tissue"] == "adipose"), :
    ].copy(deep=True)

    # Create cross tabulation.
    cross_tabulation = pandas.crosstab(
        table_tissue["study_clinic_visit"].astype("category"),
        table_tissue["intervention_text"].astype("category"),
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("module: exercise.transcriptomics.organize_sample.py")
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
                "cohort_age_text": ["elder",],
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
                "cohort_age_text": ["elder",],
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
                "cohort_age_text": ["younger",],
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
                "cohort_age_text": ["younger",],
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
                "cohort_age_text": ["younger", "elder",],
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
                "cohort_age_text": ["younger", "elder",],
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
                "cohort_age_text": ["younger", "elder",],
            },
        },
        {
            "name_set": "adipose_visit-first_age",
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "study_clinic_visit": ["first",],
            },
            "factor_availability": {
                "cohort_age_text": ["younger", "elder",],
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
                "cohort_age_text": ["elder",],
            },
            "factor_availability": {
                "study_clinic_visit": ["first", "second",],
                "intervention_text": ["placebo", "active",],
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
                "cohort_age_text": ["elder",],
                "study_clinic_visit": ["first",],
            },
            "factor_availability": {
                "intervention_text": ["placebo", "active",],
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
                "cohort_age_text": ["elder",],
                "study_clinic_visit": ["second",],
            },
            "factor_availability": {
                "intervention_text": ["placebo", "active",],
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
                "cohort_age_text": ["elder",],
                "intervention_text": ["placebo",],
            },
            "factor_availability": {
                "study_clinic_visit": ["first", "second",],
            },
        },
        {
            "name_set": str(
                "adipose_elder-active_visit"
            ),
            "tissue": "adipose",
            "cohort_selection": {
                "inclusion": [1,],
                "tissue": ["adipose",],
                "cohort_age_text": ["elder",],
                "intervention_text": ["active",],
            },
            "factor_availability": {
                "study_clinic_visit": ["first", "second",],
            },
        },
    ]
    # Return information.
    return instances


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
        print("module: exercise.transcriptomics.organize_sample.py")
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
        print("module: exercise.transcriptomics.organize_sample.py")
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
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_sample.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: exercise")
        print("routine: transcriptomics")
        print("procedure: organize_sample")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # 1. Initialize directories for read of source and write of product files.
    paths = initialize_directories(
        project="exercise",
        routine="transcriptomics",
        procedure="organize_sample",
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
    # 3. Organize table of matches between samples and files.
    translations_sample_file = define_translation_columns_table_sample_file()
    columns_sample_file = define_sequence_columns_table_sample_file()
    table_sample_file = organize_table_sample_file(
        table=pail_source["table_sample_file"],
        translations_column=translations_sample_file,
        columns_sequence=columns_sample_file,
        report=report,
    )

    ##########
    # 4. Organize table of attributes for samples.
    translations_sample_attribute = (
        define_translation_columns_table_sample_attribute()
    )
    columns_sample_attribute = define_sequence_columns_table_sample_attribute()
    table_sample_attribute = organize_table_sample_attribute(
        table=pail_source["table_sample_attribute"],
        translations_column=translations_sample_attribute,
        columns_sequence=columns_sample_attribute,
        report=report,
    )

    ##########
    # 5. Combine within the same table the matches between samples and files
    # along with their further attributes.
    columns_transfer = copy.deepcopy(columns_sample_attribute)
    columns_transfer.remove("match_sample_file_attribute")
    table_sample = combine_table_sample_file_attribute(
        table_sample_file=table_sample_file,
        table_sample_attribute=table_sample_attribute,
        columns_transfer=columns_transfer,
        report=report,
    )

    ##########
    # 6. Describe factors in table of samples.
    describe_table_sample_factors(
        table_sample=table_sample,
        report=report,
    )

    ##########
    # 7. Describe sets of samples for specific analyses.
    selections = define_selections_sample_set()
    describe_table_sample_sets(
        table_sample=table_sample,
        selections=selections,
        report=report,
    )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_tables = dict()
    pail_write_tables[str("table_sample")] = table_sample
    pail_write_objects = dict()
    #pail_write_objects[str("samples")]

    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_data"],
        reset_index=False,
        write_index=True,
        type="text",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_data"],
        reset_index=False,
        write_index=True,
        type="pickle",
    )
    if False:
        putly.write_objects_to_file_pickle(
            pail_write=pail_write_objects,
            path_directory=paths["out_data"],
        )

    pass


###############################################################################
# End
