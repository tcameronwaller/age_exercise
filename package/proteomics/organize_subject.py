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
        print("module: exercise.proteomics.organize_subject.py")
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
                threshold_proportion=0.02, # float or None to keep all
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
        print("module: organize_subject.py")
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
# 5. Prepare derivative tables of information about features and sets of
# features across observations and groups of observations.


def organize_preliminary_information_to_prepare_tables(
    index_features=None,
    index_observations=None,
    features_all=None,
    observations_all=None,
    features_inclusion=None,
    observations_inclusion=None,
    groups_features=None,
    groups_observations=None,
    names_groups_features_sequence=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    report=None,
):
    """
    Organize preliminary information for subsequent use in preparation of
    derivative tables.

    Features could correspond to columns in an original source table, or they
    could correspond to rows. Likewise, observations could correspond either
    to rows or columns. This function handles features and observations
    equivalently, so the difference is semantic.

    Review: ____

    arguments:
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        features_all (list<str>): identifiers of all features in the original
            source table
        observations_all (list<str>): identifiers of all observations in the
            original source table
        features_inclusion (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_inclusion (list<str>): identifiers of observations for
            which to include and describe values of signal intensity across
            features
        groups_features (dict<list<str>>): names of groups and identifiers
            of features that belong to each of these groups; clustering of
            features will have constraint within these groups
        groups_observations (dict<list<str>>): names of groups and identifiers
            of observations that belong to each of these groups; clustering of
            observations will have constraint within these groups
        names_groups_features_sequence (list<str>): names of groups for
            features in specific sequence
        names_groups_observations_sequence (list<str>): names of groups for
            observations in specific sequence
        translations_features (dict<str>): translations for names or
            identifiers of features
        translations_observations (dict<str>): translations for names or
            identifiers of observations
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Organize and collect information.
    pail = dict()

    # Copy other information.
    pail["index_features"] = index_features
    pail["index_observations"] = index_observations
    pail["features_all"] = copy.deepcopy(features_all)
    pail["observations_all"] = copy.deepcopy(observations_all)
    pail["features_inclusion"] = copy.deepcopy(features_inclusion)
    pail["observations_inclusion"] = copy.deepcopy(observations_inclusion)
    pail["groups_features"] = copy.deepcopy(groups_features)
    pail["groups_observations"] = copy.deepcopy(groups_observations)
    pail["names_groups_features_sequence"] = copy.deepcopy(
        names_groups_features_sequence
    )
    pail["names_groups_observations_sequence"] = copy.deepcopy(
        names_groups_observations_sequence
    )
    pail["translations_features"] = copy.deepcopy(translations_features)
    pail["translations_observations"] = copy.deepcopy(
        translations_observations
    )

    ##########
    # Prepare source information.
    # Prepare reference to sort rows or columns in a subsequent table.
    # Option 1.
    if True:
        # Features.
        sequence_groups_features = dict()
        index = 0
        for name in pail["names_groups_features_sequence"]:
            sequence_groups_features[name] = index
            index += 1
            pass
        # Observations.
        sequence_groups_observations = dict()
        index = 0
        for name in pail["names_groups_observations_sequence"]:
            sequence_groups_observations[name] = index
            index += 1
            pass
    # Option 2.
    if False:
        sequence_groups_observations = dict(zip(
            pail["names_groups_observations_sequence"],
            range(len(names_groups_observations_sequence))
        ))
    # Option 3.
    if False:
        sequence_groups_observations = {
            key: value for value, key in enumerate(
                pail["names_groups_observations_sequence"]
            )
        }
    pail["sequence_groups_features"] = sequence_groups_features
    pail["sequence_groups_observations"] = sequence_groups_observations
    # Ensure that features are unique.
    features_inclusion_unique = putly.collect_unique_elements(
        elements=pail["features_inclusion"],
    )
    pail["features_inclusion_unique"] = features_inclusion_unique
    # Ensure that all features are in the source table.
    features_available = list(filter(
        lambda feature: (feature in pail["features_all"]),
        pail["features_inclusion_unique"]
    ))
    pail["features_available"] = features_available
    # Translate identifiers of features.
    features_available_translation = copy.deepcopy(
        features_available
    )
    if (pail["translations_features"] is not None):
        features_available_translation = list(map(
            lambda feature: (pail["translations_features"][feature]),
            features_available_translation
        ))
        pass
    # Ensure that features are unique after translation.
    features_available_translation = putly.collect_unique_elements(
        elements=features_available_translation,
    )
    pail["features_available_translation"] = features_available_translation
    # Ensure that observations are unique.
    observations_inclusion_unique = putly.collect_unique_elements(
        elements=pail["observations_inclusion"],
    )
    pail["observations_inclusion_unique"] = observations_inclusion_unique
    # Ensure that all observations are in the source table.
    observations_available = list(filter(
        lambda observation: (observation in pail["observations_all"]),
        observations_inclusion_unique
    ))
    pail["observations_available"] = observations_available
    # Translate identifiers of features.
    observations_available_translation = copy.deepcopy(
        observations_available
    )
    if (pail["translations_observations"] is not None):
        observations_available_translation = list(map(
            lambda observation: (
                pail["translations_observations"][observation]
            ),
            observations_available_translation
        ))
        pass
    # Ensure that features are unique after translation.
    observations_available_translation = putly.collect_unique_elements(
        elements=observations_available_translation,
    )
    pail["observations_available_translation"] = (
        observations_available_translation
    )

    # Return information.
    return pail


def prepare_tables_signals_for_features_sets_in_observations_groups(
    table=None,
    index_features=None,
    index_observations=None,
    features_inclusion=None,
    observations_inclusion=None,
    groups_features=None,
    groups_observations=None,
    names_groups_features_sequence=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    report=None,
):
    """
    Prepare derivative tables of information about signals, measurements, or
    other values corresponding to features and sets of features across
    observations and groups of observations.

    Any translations of identifiers or names of features and observations occur
    after the selection of features and observations. Hence identifiers or
    names for selection of features and observations must match those in the
    original source table.

    Each observation must only belong to a single group. That is, the groups of
    observations must be exclusive.

    Each feature can belong to multiple sets. That is, the sets of features are
    not exclusive.

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

    Format of product table 1 (name: "table_product_1")
    Format of product table 1 is in wide format with floating-point values of
    signal intensities or other measurement values corresponding to features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns. This novel product table
    shares its format with the original source table. The difference between
    them is that this novel product table only includes a specific selection of
    features and observations from the original source table.
    ----------
    observation     feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    Review: ____

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across rows and observations across
            columns
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        features_inclusion (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_inclusion (list<str>): identifiers of observations for
            which to include and describe values of signal intensity across
            features
        groups_features (dict<list<str>>): names of groups and identifiers
            of features that belong to each of these groups; clustering of
            features will have constraint within these groups
        groups_observations (dict<list<str>>): names of groups and identifiers
            of observations that belong to each of these groups; clustering of
            observations will have constraint within these groups
        names_groups_features_sequence (list<str>): names of groups for
            features in specific sequence
        names_groups_observations_sequence (list<str>): names of groups for
            observations in specific sequence
        translations_features (dict<str>): translations for names or
            identifiers of features
        translations_observations (dict<str>): translations for names or
            identifiers of observations
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information in table.
    table_source = table.copy(deep=True)

    ##########
    # Organize preliminary information.
    features_all = copy.deepcopy(table_source.columns.to_list())
    features_all.remove(index_observations)
    observations_all = copy.deepcopy(
        table_source[index_observations].unique().tolist()
    )
    pail = organize_preliminary_information_to_prepare_tables(
        index_features=index_features,
        index_observations=index_observations, # assigned in new tables
        features_all=features_all,
        observations_all=observations_all,
        features_inclusion=features_inclusion,
        observations_inclusion=observations_inclusion,
        groups_features=groups_features,
        groups_observations=groups_observations,
        names_groups_features_sequence=names_groups_features_sequence,
        names_groups_observations_sequence=(
            names_groups_observations_sequence
        ),
        translations_features=translations_features,
        translations_observations=translations_observations,
        report=report,
    )


    pass




#############################
# Below here can move to "plot_subject.py"

##########
# 4. Histograms


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
    # partner.organization.extract_organize_values_from_series()
    values_raw = table[column_feature].to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    values_nonmissing = numpy.copy(values_raw[~numpy.isnan(values_raw)])
    mean_values = round(numpy.nanmean(pail_values["values_nonmissing"]), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_distribution_histogram(
        array=pail_values["values_nonmissing"],
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


##########
# 6. Plot charts of scatter points for principal components


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
            "name": "cohort_age_text",
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


###############################################################################
# Procedure


# TODO: TCW; 5 December 2024
# DONE
# 1. write out lists of the IDs for Olink signals from Adipose, Plasma, and Muscle
# 2. place these file lists of Olink IDs in parameter directory "sets_olink"

# TODO
# 3. calculate correlations between...
#   - subject continuous traits and Olink signals
#   - pairwise Olink signals
#   - subject continuous traits and PCs on Olink signals
#   -


# table = table_subject
# index_features = "features" # across columns... function will assign
# index_observations = "subject_visit"

#    pail_tables = (
#        pdesc.describe_signals_for_features_sets_in_observations_groups(
#            table=table_signal,
#            index_features=index_genes,
#            index_observations=index_samples, # assigned in new tables
#            features_inclusion=genes_inclusion,
#            observations_inclusion=samples_inclusion,
#            groups_features=pail_gene_sets["sets_gene_cluster"],
#            groups_observations=groups_samples,
#            names_groups_features_sequence=(
#                instances_group[0]["sequence_sets_gene_cluster"]
#            ),
#            names_groups_observations_sequence=(
#                names_groups_samples_sequence
#            ),
#            translations_features=translations_gene,
#            translations_observations=None,
#            report=report,
#    ))

# separate function <-- prepare correlations...



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
    # 4. Organize table of properties for study subjects.
    columns_original = pail_parse["columns_all"]
    columns_novel = define_sequence_columns_novel_subject_feature()
    pail_organization = organize_table_subject_property(
        table=pail_source["table_subject_property"],
        translations_column=pail_parse["translations_column"],
        columns_original=columns_original,
        columns_novel=columns_novel,
        report=report,
    )
    table_subject = pail_organization["table"]

    ##########
    # 5. Plot histogram charts for features on quantitative continuous interval
    # or ratio scales of measurement.

    # TODO: TCW; 27 November 2024
    # TODO: include another "selection" column in the parameter table to
    # designate features for logarithmic scale. Then will need to include
    # those is a new column list in the "parse" function.
    if False:
        pail_columns_histogram = dict()
        pail_columns_histogram["other"] = copy.deepcopy(
            pail_parse["columns_continuous"]
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


    ######################################
    # TODO: TCW; 30 November 2024
    # TODO: stratify the subjects by sex and age group before calculating
    # principal components
    #######################################
    # Filter rows in table by specific categorical values of specific columns.
    columns_categories = dict()
    columns_categories["cohort_age_text"] = ["younger", "elder",]
    #columns_categories["cohort_age_text"] = ["younger",]
    #columns_categories["cohort_age_text"] = ["elder",]
    columns_categories["sex_text"] = ["female", "male",]
    #columns_categories["sex_text"] = ["female",]
    #columns_categories["sex_text"] = ["male",]
    table_subject = porg.filter_select_table_rows_by_columns_categories(
        table=table_subject,
        columns_categories=columns_categories,
        report=report,
    )

    ##########
    # 6. Calculate principal components on features for Olink measurements.
    # Olink measurements are only available for "Pre" or "first" clinical
    # visits of the study.
    columns_olink_plasma = copy.deepcopy(pail_parse["columns_olink_plasma"])
    columns_olink_muscle = copy.deepcopy(pail_parse["columns_olink_muscle"])
    columns_olink_adipose = copy.deepcopy(pail_parse["columns_olink_adipose"])
    columns_olink_all = copy.deepcopy(columns_olink_plasma)
    columns_olink_all.extend(columns_olink_muscle)
    columns_olink_all.extend(columns_olink_adipose)
    pail_columns_decomposition = dict()
    pail_columns_decomposition["olink_all"] = columns_olink_all
    pail_columns_decomposition["olink_plasma"] = columns_olink_plasma
    pail_columns_decomposition["olink_muscle"] = columns_olink_muscle
    pail_columns_decomposition["olink_adipose"] = columns_olink_adipose
    if True:
        pail_reduction = drive_manage_calculate_principal_components(
                table=table_subject,
                index_rows="subject_visit",
                sets_columns=pail_columns_decomposition,
                report=report,
        )
    else:
        pail_pca_test = (
            pdecomp.calculate_principal_components_table_columns_selection(
                table=table_subject,
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


    if False:
        ##########
        # 7. Plot charts of scatter points to represent principal components on
        # sets of features for Olink measurements.
        features_response = define_response_features_principal_components_olink()
        drive_manage_plot_write_principal_component_scatter_charts(
            table=pail_reduction["table"],
            features_response=features_response,
            sets_columns=pail_reduction["sets_columns"],
            paths=paths,
            report=report,
        )

    ##########
    # 8. Prepare derivative tables of information about features and sets of
    # features across observations and groups of observations.
    pail_tables = (
        prepare_tables_signals_for_features_sets_in_observations_groups(
            table=table_subject,
            index_features=None,
            index_observations=None,
            features_inclusion=None,
            observations_inclusion=None,
            groups_features=None,
            groups_observations=None,
            names_groups_features_sequence=None,
            names_groups_observations_sequence=None,
            translations_features=None,
            translations_observations=None,
            report=None,
    ))


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

    pass


###############################################################################
# End
