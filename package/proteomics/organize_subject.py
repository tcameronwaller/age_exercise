"""
Supply functionality for process and analysis of data from proteomics.

This module 'organize_subject' is part of the 'proteomics' package within
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
        print("module: age_exercise.proteomics.organize_subject.py")
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
    types_columns["selection_continuous"] = "int32"
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
        table_inclusion["name_product"].unique().tolist()
    )

    # Determine names of columns for feature variables on quantitative
    # continuous interval or ratio scales of measurement.
    table_continuous = table.loc[
        (
            (table[inclusion] == 1) &
            (table["selection_continuous"] == 1) &
            (~table["category"].str.contains("olink_"))
        ), :
    ].copy(deep=True)
    columns_continuous = copy.deepcopy(
        table_continuous["name_product"].to_list()
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
        print("module: age_exercise.transcriptomics.organize_sample.py")
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
    pail["columns_continuous"] = columns_continuous
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
    pail_parse = parse_extract_table_sample_feature_organization(
        table=pail["table_feature_organization"],
        inclusion="inclusion_proteomics",
        report=report,
    )
    pail["columns_all"] = pail_parse["columns_all"]
    pail["translations_column"] = pail_parse["translations_column"]

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
        print("module: age_exercise.proteomics.organize_subject.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for subjects and samples: ")
        print(pail["table_subject_property"].iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Organize table of properties for study subjects.


def define_sequence_columns_novel_subject_feature():
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
        "cohort_age",
        "cohort_age_text",
        #"cohort_age_letter",
        "intervention",
        "intervention_text",
        #"identifier_subject",
        #"study_clinic_visit_relative",
        "study_clinic_visit",
        "subject_visit",
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


def determine_subject_study_clinic_visit(
    visit_relative=None,
):
    """
    Determines the clinical visit of the study at which collection of a
    sample occurred.

    arguments:
        visit_relative (str): indicator of clinical visit in the study at which
            collection of a sample occurred, either 'Pre' or 'Post'

    raises:

    returns:
        (str): indicator of clinical visit in the study at which collection of
            a sample occurred, either 'first' or 'second'

    """

    # Determine designator.
    if (
        (pandas.notna(visit_relative)) and
        (len(str(visit_relative).strip()) > 0)
    ):
        # There is adequate information.
        visit_relative = str(visit_relative).strip().lower()
        if (visit_relative == "pre"):
            visit = str("first")
        elif (visit_relative == "post"):
            visit = str("second")
        else:
            visit = ""
    else:
        visit = ""
        pass
    # Return information.
    return visit


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

# TODO: TCW; 27 November 2024
# TODO: include another "selection" column in the parameter table to
# designate features for logarithmic scale. Then will need to include
# those is a new column list in the "parse" function.


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
            (table["identifier_subject"].str.len() > 0)
        ), :
    ].copy(deep=True)
    table.dropna(
        how="all",
        axis="index",
    )

    # Determine designation for subject's first or second clinical visit of the
    # study.
    table["study_clinic_visit"] = table.apply(
        lambda row:
            determine_subject_study_clinic_visit(
                visit_relative=row["study_clinic_visit_relative"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine combination designation for subject and visit.
    table["subject_visit"] = table.apply(
        lambda row: str(
            row["identifier_subject"] + "_" + row["study_clinic_visit"]
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

    # Transform values to logarithmic scale for a selection of features on a
    # quantitative continuous interval or ratio scale of measurement.

    # TODO: TCW; 27 November 2024
    # TODO: include another "selection" column in the parameter table to
    # designate features for logarithmic scale. Then will need to include
    # those is a new column list in the "parse" function.

    # Sort rows within table.
    table.sort_values(
        by=[
            "cohort_age",
            "intervention",
            "identifier_subject",
            "study_clinic_visit",
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
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        print("function: organize_table_subject_property()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table"].iloc[0:10, 0:])
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    procedure="organize_subject"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
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
    # 3. Organize table of properties for study subjects.
    columns_original = pail_source["columns_all"]
    columns_novel = define_sequence_columns_novel_subject_feature()
    pail_organization = organize_table_subject_property(
        table=pail_source["table_subject_property"],
        translations_column=pail_source["translations_column"],
        columns_original=columns_original,
        columns_novel=columns_novel,
        report=report,
    )

    # TODO: TCW; 27 November 2024
    # TODO: include another "selection" column in the parameter table to
    # designate features for logarithmic scale. Then will need to include
    # those is a new column list in the "parse" function.


    ##########
    # 4. Collect information.
    # Collections of files.

    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_subject")] = pail_organization["table"]
    pail_write_objects = dict()
    #pail_write_objects[str("samples")]

    ##########
    # 5. Write product information to file.
    if False:
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

    pass


###############################################################################
# End
