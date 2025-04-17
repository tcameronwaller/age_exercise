"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'organize_subject' is part of the 'phenotypes' subpackage within
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
    initialize_routine=None,
    restore=None,
    report=None,
):
    """
    Initialize directories for procedure's source and product files.

    arguments:
        project (str): name of project that normally corresponds to a single
            Python package
        routine (str): name of routine that normally corresponds to a single
            Python package or subpackage, either 'phenotypes',
            'transcriptomics', or 'proteomics'
        procedure (str): name of procedure, a step in the routine process that
            normally corresponds to a single Python module within the package
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        initialize_routine (bool): whether to initialize directory for routine
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
    paths["out_procedure_object"] = os.path.join(
        paths["out_procedure"], "object",
    )

    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"],
        paths["out_procedure_lists"],
        paths["out_procedure_tables"],
        paths["out_procedure_plot"],
        paths["out_procedure_object"],
    ]
    if (initialize_routine):
        paths_initialization.insert(0, paths["out_routine"],)
        pass
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
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: organize_subject.py")
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


def read_source_demonstration(
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
    path_file_table_demonstration_signal = os.path.join(
        paths["in_demonstration"], "partner",
        "table_columns_features_rows_observations_groups.tsv",
    )
    path_file_table_demonstration_set = os.path.join(
        paths["in_demonstration"], "partner",
        "table_allocation_features_sets.tsv",
    )

    # Collect information.
    pail = dict()
    # Read information from file.
    #types_columns = define_column_types_table_demonstration()
    pail["table_demonstration_signal"] = pandas.read_csv(
        path_file_table_demonstration_signal,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_demonstration_set"] = pandas.read_csv(
        path_file_table_demonstration_set,
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
        print("module: age_exercise.transcriptomics.compare_sets_groups.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        #count_rows = (pail["table_demonstration"].shape[0])
        #count_columns = (pail["table_demonstration"].shape[1])
        #print("demonstration table: ")
        #print("count of rows in table: " + str(count_rows))
        #print("Count of columns in table: " + str(count_columns))
        #print(pail["table_demonstration"])
        #putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    types_columns["inclusion"] = "int32"
    types_columns["selection_quantitative"] = "int32"
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

    # Extract information for forward translation of names of columns.
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
    translations_feature_forward = series_translations.to_dict()

    # Extract information for reverse translation of names of columns.
    table_name_product = table_inclusion.loc[
        (table_inclusion["name_product"].str.len() > 0), :
    ].copy(deep=True)
    table_translations = table_name_product.filter(
        items=["name_source", "name_product",],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_translations["name_source"].to_list(),
        index=table_translations["name_product"],
    )
    translations_feature_reverse = series_translations.to_dict()

    # Extract product names of all columns.
    columns_all = copy.deepcopy(
        table_inclusion["name_product"].unique().tolist()
    )

    # Determine names of columns for feature variables on quantitative
    # continuous interval or ratio scales of measurement.
    table_quantitative = table.loc[
        (
            (table[inclusion] == 1) &
            (table["selection_quantitative"] == 1) &
            (~table["category"].str.contains("olink_"))
        ), :
    ].copy(deep=True)
    columns_quantitative = copy.deepcopy(
        table_quantitative["name_product"].to_list()
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
            "forward translations for columns after initial, intermediate " +
            "translation:"
        )
        print(translations_feature_forward)
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
    pail["translations_feature_forward"] = translations_feature_forward
    pail["translations_feature_reverse"] = translations_feature_reverse
    pail["columns_all"] = columns_all
    pail["columns_quantitative"] = columns_quantitative
    pail["columns_olink_plasma"] = columns_olink_plasma
    pail["columns_olink_muscle"] = columns_olink_muscle
    pail["columns_olink_adipose"] = columns_olink_adipose
    # Return information.
    return pail


def read_organize_source_table_subject_feature_organization(
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

    # Read information from file.

    # Table of parameters for organization of the table of attributes for
    # subjects and samples.
    types_columns = define_type_columns_table_subject_feature_organization()
    table_feature_organization = pandas.read_csv(
        path_file_table_feature_organization,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail = parse_extract_table_sample_feature_organization(
        table=table_feature_organization,
        inclusion="inclusion",
        report=report,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: organize_subject.py")
        name_function = str(
            "read_organize_source_table_subject_feature_organization()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass
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
    path_file_table_subject_property = os.path.join(
        paths["in_data"], "study_age_exercise", "subject_sample",
        "table_subject_sample_feature_olink_hek_2024-11-22.csv",
    )

    # Read information from file.

    # Table of parameters for organization of the table of attributes for
    # subjects and samples.
    pail_parse = read_organize_source_table_subject_feature_organization(
        paths=paths,
        report=paths,
    )

    # Table of properties and attributes for subjects and samples.
    table_subject_property = pandas.read_csv(
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

    # Collect information.
    pail = dict()
    pail["table_subject_property"] = table_subject_property
    pail["translations_feature_forward"] = (
        pail_parse["translations_feature_forward"]
    )
    pail["translations_feature_reverse"] = (
        pail_parse["translations_feature_reverse"]
    )
    pail["columns_all"] = pail_parse["columns_all"]
    pail["columns_quantitative"] = pail_parse["columns_quantitative"]

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
# Exercise caution and pay attention to the definitions of the logical binary
# features.
# Some definitions are 'lazy', such that they have values of zero (0) when they
# perhaps ought to have missing values.
# For this reason, for clarity, filter and stratify observations by the
# corresponding features with categorical text values.


def define_sequence_columns_novel_subject_feature():
    """
    Defines names of columns in sequence by which to filter and sort columns in
    a table.

    This list represents the columns that are novel derivations of the original
    columns. These features are not in the source table
    "table_subject_sample_feature_organization".

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "identifier_observation_trial",
        #"identifier_subject",
        "age_cohort_text",
        "age_cohort_elder",
        "sex_text",
        "sex_female",
        "sex_y",
        "intervention_text",
        "intervention_placebo",
        "intervention_omega3",
        "intervention_text_placebo_versus_other",
        "intervention_text_omega3_versus_other",
        "intervention_placebo_other",
        "intervention_omega3_other",
        "intervention_text_after_placebo_versus_other",
        "intervention_text_after_omega3_versus_other",
        "intervention_after_placebo",
        "intervention_after_omega3",
        "visit_text",
        "visit_second",
        "visit_text_fill_second_for_age_younger",
        #"subject_tissue_visit",
        "subject_visit",
        "date_visit_text",
        #"date_visit_text_raw",
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


def determine_sex_text(
    sex_text_raw=None,
):
    """
    Determines a designator for sex.

    arguments:
        sex_text_raw (str): designator of sex

    raises:

    returns:
        (str): designator of sex

    """

    # Determine designator.
    if (
        (pandas.notna(sex_text_raw)) and
        (len(str(sex_text_raw).strip()) > 0)
    ):
        # There is adequate information.
        sex_text_raw = str(sex_text_raw).strip().upper()
        if (sex_text_raw == "F"):
            designator = str("female")
        elif (sex_text_raw == "M"):
            designator = str("male")
        else:
            designator = ""
    else:
        designator = ""
        pass
    # Return information.
    return designator


def determine_subject_visit_text(
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


def determine_age_cohort_text(
    age_cohort_text_raw=None,
):
    """
    Determines a designator for cohort by age.

    arguments:
        age_cohort_text_raw (str): designator of cohort by age

    raises:

    returns:
        (str): designator of cohort by age

    """

    # Determine designator.
    if (
        (pandas.notna(age_cohort_text_raw)) and
        (len(str(age_cohort_text_raw).strip()) > 0)
    ):
        # There is adequate information.
        age_cohort_text_raw = str(age_cohort_text_raw).strip().upper()
        if (age_cohort_text_raw == "E"):
            designator = str("elder")
        elif (age_cohort_text_raw == "Y"):
            designator = str("younger")
        else:
            designator = ""
    else:
        designator = ""
        pass
    # Return information.
    return designator


def determine_intervention_text(
    intervention_text_raw=None,
    age_cohort_text=None,
):
    """
    Determines a designator for intervention.

    arguments:
        intervention_text_raw (str): designator of intervention
        age_cohort_text (str): designator of cohort by age

    raises:

    returns:
        (float): designator of intervention

    """

    # Determine designator.
    if (
        (pandas.notna(intervention_text_raw)) and
        (len(str(intervention_text_raw).strip()) > 0)
    ):
        # Determine designator.
        designator = str(intervention_text_raw).strip().lower()
        # Translate designator.
        if (designator == "active"):
            designator = "omega3"
    else:
        # Fill with parsable term for not applicable.
        designator = "none"
        pass
    # Interpret designator and fill or change as necessary.
    if (
        (pandas.notna(age_cohort_text)) and
        (len(str(age_cohort_text).strip()) > 0) and
        (
            (age_cohort_text == "younger") or
            (designator not in ["placebo", "omega3",])
        )
    ):
        # Fill with parsable term for not applicable.
        designator = "none"
        pass
    # Return information.
    return designator


def determine_intervention_text_after_intervention_versus_other(
    type=None,
    intervention_text=None,
    visit_text=None,
):
    """
    Determines a designator.

    arguments:
        type (str): type of intervention to designate, either 'placebo' or
            'omega3'
        intervention_text (str): designator of intervention
        visit_text (str): designator of clinical visit in the study at which
            collection of a sample occurred, either 'first' or 'second'

    raises:

    returns:
        (float): designator of intervention

    """

    #table["intervention_text_after_placebo_versus_other"] = table.apply(
    #    lambda series_row: (
    #        str("other")
    #        if (series_row["visit_text"] in ["first",])
    #        else series_row["intervention_text_placebo_versus_other"]
    #    ),
    #    axis="columns", # apply function to each row
    #)
    #table["intervention_text_after_omega3_versus_other"] = table.apply(
    #    lambda series_row: (
    #        str("other")
    #        if (series_row["visit_text"] in ["first",])
    #        else series_row["intervention_text_omega3_versus_other"]
    #    ),
    #    axis="columns", # apply function to each row
    #)

    # Determine designator.
    if (
        (visit_text == "second") and
        (intervention_text == type)
    ):
        # Determine designator.
        designator = type
    else:
        # Fill with parsable term for not applicable.
        designator = "other"
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


def wrangle_messy_features(
    table=None,
    columns_original=None,
    columns_novel=None,
    report=None,
):
    """
    Wrangle, parse, clean, or organize information for features that are messy
    or otherwise need extra help.

    arguments:
        table (object): Pandas data-frame table of subjects, samples, and their
            attribute features
        columns_original (list<str>): names of original columns in sequence by
            which to filter and sort columns in table
        columns_novel (list<str>): names of original columns in sequence by
            which to filter and sort columns in table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information.
    table = table.copy(deep=True)
    columns_original = copy.deepcopy(columns_original)
    columns_novel = copy.deepcopy(columns_novel)

    # Combine original and novel columns.
    columns_all = copy.deepcopy(columns_original)
    columns_all.extend(columns_novel)

    # Clean values for counts of white blood cells.
    if ("white_blood_cells" in columns_all):
        table["white_blood_cells"] = pandas.to_numeric(
            table["white_blood_cells"],
            errors="coerce", # force any parse error values to missing "NaN"
            downcast="float", # cast type to smallest float type
        )
        pass

    # Clean values for counts of eosinophils.
    # Use the half minimum method of imputation.
    # Imputation seems justifiable since there is more information than just a
    # missing value.
    if ("eosinophils" in columns_all):
        table["eosinophils"] = table.apply(
            lambda row:
                str(row["eosinophils"]).replace(
                    "<.03", "0.015"
                ).replace("<0.03", "0.015"),
            axis="columns", # apply function to each row
        )
        table["eosinophils"] = table["eosinophils"].astype("float32")
        pass

    # Clean values for counts of basophils.
    # Use the half minimum method of imputation.
    # Imputation seems justifiable since there is more information than just a
    # missing value.
    if ("basophils" in columns_all):
        table["basophils"] = table.apply(
            lambda row:
                str(row["basophils"]).replace(
                    "<.03", "0.015"
                ).replace("<0.03", "0.015"),
            axis="columns", # apply function to each row
        )
        table["basophils"] = table["basophils"].astype("float32")
        table["basophils"].where(
            table["basophils"] <= 100,
            other=pandas.NA,
            inplace=True,
        )
        pass

    # Clean values for C-reactive protein variable.
    # Use the half minimum method of imputation.
    # Imputation seems justifiable since there is more information than just a
    # missing value.
    if ("c_react_protein" in columns_all):
        table["c_react_protein"] = table.apply(
            lambda row:
                str(row["c_react_protein"]).replace(
                    "<.2", "0.1"
                ).replace("<0.2", "0.1"),
            axis="columns", # apply function to each row
        )
        table["c_react_protein"] = table["c_react_protein"].astype("float32")
        table["c_react_protein"].where(
            table["c_react_protein"] <= 30,
            other=pandas.NA,
            inplace=True,
        )
        pass

    # Replace values of zero for oxygen consumption with missing values.
    # Use the values of maximal oxygen consumption adjusted to lean body mass.
    if ("oxygen_consumption" in columns_all):
        table["oxygen_consumption"] = table["oxygen_consumption"].replace(
            to_replace=0,
            value=pandas.NA,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        print("function: wrangle_messy_features()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail["table"].iloc[0:10, 0:])
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


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

    # Wrangle, parse, clean, or organize information for features that are
    # messy or otherwise need extra help.
    table = wrangle_messy_features(
        table=table,
        columns_original=columns_original,
        columns_novel=columns_novel,
        report=False,
    )

    # Identifiers of observations.
    # Name generic index in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.index.set_names("identifier_observation_trial", inplace=True)
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    # Determine designations of sex.
    table["sex_text"] = table.apply(
        lambda row:
            determine_sex_text(
                sex_text_raw=row["sex_text_raw"],
            ),
        axis="columns", # apply function to each row
    )
    table["sex_female"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row["sex_text"],
                values_1=["female",],
                values_0=["male",],
            ),
        axis="columns", # apply function to each row
    )
    table["sex_y"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row["sex_text"],
                values_1=["male",],
                values_0=["female",],
            ),
        axis="columns", # apply function to each row
    )

    # Determine designations of cohort by age.
    table["age_cohort_text"] = table.apply(
        lambda row:
            determine_age_cohort_text(
                age_cohort_text_raw=row["age_cohort_text_raw"],
            ),
        axis="columns", # apply function to each row
    )

    table["age_cohort_elder"] = table.apply(
        lambda row: putly.determine_logical_binary_designator(
                category_text=row["age_cohort_text"],
                values_1=["elder",],
                values_0=["younger",],
            ),
        axis="columns", # apply function to each row
    )

    # Determine designation for subject's first or second clinical visit of the
    # study.
    table["visit_text"] = table.apply(
        lambda row:
            determine_subject_visit_text(
                visit_relative=row["study_clinic_visit_text_raw"],
            ),
        axis="columns", # apply function to each row
    )
    table["visit_second"] = table.apply(
        lambda row:
        putly.determine_logical_binary_designator(
            category_text=row["visit_text"],
            values_1=["second",],
            values_0=["first",],
        ),
        axis="columns", # apply function to each row
    )
    table["visit_text_fill_second_for_age_younger"] = table.apply(
        lambda series_row: (
            str("second")
            if (series_row["age_cohort_text"] == "younger")
            else series_row["visit_text"]
        ),
        axis="columns", # apply function to each row
    )

    # Determine combination designation for subject and visit.
    table["subject_visit"] = table.apply(
        lambda row: str(
            row["identifier_subject"] + "_" + row["visit_text"]
        ),
        axis="columns", # apply function to each row
    )
    #table["subject_tissue_visit"] = table.apply(
    #    lambda row: str(
    #        row["identifier_subject"] + "_" +
    #        row["tissue"] + "_" +
    #        row["visit_text"]
    #    ),
    #    axis="columns", # apply function to each row
    #)

    # Determine designations of intervention, either placebo or active.
    table["intervention_text"] = table.apply(
        lambda row:
        determine_intervention_text(
            intervention_text_raw=row["intervention_text_raw"],
            age_cohort_text=row["age_cohort_text"],
        ),
        axis="columns", # apply function to each row
    )
    table["intervention_placebo"] = table.apply(
        lambda row:
        putly.determine_logical_binary_designator(
            category_text=row["intervention_text"],
            values_1=["placebo",],
            values_0=["omega3",],
        ),
        axis="columns", # apply function to each row
    )
    table["intervention_omega3"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row["intervention_text"],
                values_1=["omega3",],
                values_0=["placebo",],
            ),
        axis="columns", # apply function to each row
    )
    table["intervention_text_placebo_versus_other"] = table.apply(
        lambda series_row: (
            str("other")
            if (series_row["intervention_text"] in ["omega3", "none",])
            else series_row["intervention_text"]
        ),
        axis="columns", # apply function to each row
    )
    table["intervention_text_omega3_versus_other"] = table.apply(
        lambda series_row: (
            str("other")
            if (series_row["intervention_text"] in ["placebo", "none",])
            else series_row["intervention_text"]
        ),
        axis="columns", # apply function to each row
    )
    abbreviation = "intervention_text_placebo_versus_other"
    table["intervention_placebo_other"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row[abbreviation],
                values_1=["placebo",],
                values_0=["other",],
            ),
        axis="columns", # apply function to each row
    )
    abbreviation = "intervention_text_omega3_versus_other"
    table["intervention_omega3_other"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row[abbreviation],
                values_1=["omega3",],
                values_0=["other",],
            ),
        axis="columns", # apply function to each row
    )
    table["intervention_text_after_placebo_versus_other"] = table.apply(
        lambda row:
            determine_intervention_text_after_intervention_versus_other(
                type="placebo",
                intervention_text=row["intervention_text"],
                visit_text=row["visit_text"],
            ),
        axis="columns", # apply function to each row
    )
    table["intervention_text_after_omega3_versus_other"] = table.apply(
        lambda row:
            determine_intervention_text_after_intervention_versus_other(
                type="omega3",
                intervention_text=row["intervention_text"],
                visit_text=row["visit_text"],
            ),
        axis="columns", # apply function to each row
    )
    abbreviation = "intervention_text_after_placebo_versus_other"
    table["intervention_after_placebo"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row[abbreviation],
                values_1=["placebo",],
                values_0=["other",],
            ),
        axis="columns", # apply function to each row
    )
    abbreviation = "intervention_text_after_omega3_versus_other"
    table["intervention_after_omega3"] = table.apply(
        lambda row:
            putly.determine_logical_binary_designator(
                category_text=row[abbreviation],
                values_1=["omega3",],
                values_0=["other",],
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

    # Transform values to logarithmic scale for a selection of features on a
    # quantitative continuous interval or ratio scale of measurement.

    # TODO: TCW; 27 November 2024
    # TODO: include another "selection" column in the parameter table to
    # designate features for logarithmic scale. Then will need to include
    # those is a new column list in the "parse" function.

    # Sort rows within table.
    table.sort_values(
        by=[
            "age_cohort_elder",
            "intervention_omega3",
            "identifier_subject",
            "visit_text",
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


##########
# Functionality of use in multiple modules of the "age_exercise" package.


# Prepare derivative, deliverable, product tables for individual signals and
# their means corresponding to features in sets across observations in groups.


def organize_preliminary_information_to_prepare_tables_signal(
    index_features=None,
    index_observations=None,
    features_all=None,
    observations_all=None,
    features_selection=None,
    observations_selection=None,
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

    Review: 6 December 2024

    arguments:
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        features_all (list<str>): identifiers of all features in the original
            source table
        observations_all (list<str>): identifiers of all observations in the
            original source table
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_selection (list<str>): identifiers of observations for
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

    # Copy information.
    pail["index_features"] = index_features
    pail["index_observations"] = index_observations
    pail["features_all"] = copy.deepcopy(features_all)
    pail["observations_all"] = copy.deepcopy(observations_all)
    pail["features_selection"] = copy.deepcopy(features_selection)
    pail["observations_selection"] = copy.deepcopy(observations_selection)
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
    features_selection_unique = putly.collect_unique_elements(
        elements=pail["features_selection"],
    )
    pail["features_selection_unique"] = features_selection_unique
    # Ensure that all features are in the source table.
    features_available = list(filter(
        lambda feature: (feature in pail["features_all"]),
        pail["features_selection_unique"]
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
    observations_selection_unique = putly.collect_unique_elements(
        elements=pail["observations_selection"],
    )
    pail["observations_selection_unique"] = observations_selection_unique
    # Ensure that all observations are in the source table.
    observations_available = list(filter(
        lambda observation: (observation in pail["observations_all"]),
        observations_selection_unique
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


# TODO: TCW; 29 January 2025
# prepare_table_signal_summaries_features_observation_groups()
# this function belongs in the "description" module of package "partner"
def prepare_table_signal_summaries_features_observation_groups(
    table=None,
    index_columns=None,
    column_group=None,
    columns_features=None,
    index_rows=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    summary=None,
    report=None,
):
    """
    Prepare derivative tables of information about statistical summaries of
    measurement signal intensities corresponding to individual features across
    individual observations in groups.

    Any translations of identifiers or names of features and observations occur
    after the selection of features and observations. Hence identifiers or
    names for selection of features and observations must match those in the
    original source table.

    ----------
    Format of source table (name: "table_source")
    ----------
    Format of source table is in wide format with features across columns and
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations, with multiple observations in each
    group. For versatility, this table does not have explicitly defined indices
    across columns or rows. Values for observations of features are on a
    quantitative, continuous, interval or ratio scale of measurement.
    ----------
    observation     group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_1   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_2   0.001     0.001     0.001     0.001     0.001
    observation_6   group_2   0.001     0.001     0.001     0.001     0.001
    ----------

    ----------
    Format of product table (name: "table_product")
    ----------
    Format of product table is in wide format with floating-point values of
    a single, specific type of descriptive statistics (usually either mean or
    median) corresponding to features across rows and groups of observations
    across columns.
    ----------
    group     group_1 group_2 group_3 group_4
    feature
    feature_1 0.01    0.001   0.001   0.015
    feature_2 0.01    0.001   0.001   0.015
    feature_3 0.01    0.001   0.001   0.015
    feature_4 0.01    0.001   0.001   0.015
    feature_5 0.01    0.001   0.001   0.015
    ----------

    Review: 30 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        index_columns (str): name for header row index corresponding to
            features across columns in the original source table
        column_group (str): name of column corresponding to groups of
            observations across rows in the original source table
        columns_features (list<str>): names of columns in source table for
            features
        index_rows (str): name of a column for index corresponding to
            observations across rows in the original source table
        names_groups_observations_sequence (list<str>): names of groups for
            observations in specific sequence
        translations_features (dict<str>): translations for names or
            identifiers of features
        summary (str): name of statistical measure to select for summaries of
            signal intensities within groups of observations
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    table_source = table.copy(deep=True)
    columns_features = copy.deepcopy(columns_features)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_features = copy.deepcopy(translations_features)

    # Prepare summary of descriptive statistics for features across
    # observations in groups.
    table_long = pdesc.describe_table_features_by_groups(
        table=table_source,
        column_group=column_group,
        columns_features=columns_features,
        index_feature=index_columns,
        translations_feature=None,
        threshold_observations=5,
        digits_round=3,
        report=False,
    )

    # Organize indices in table.
    table_long.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_wide = table_long.pivot(
        columns=[column_group,],
        index=index_columns,
        values=summary,
    )
    # Organize indices in table.
    table_wide.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_wide.set_index(
        [index_columns,],
        append=False,
        drop=True,
        inplace=True,
    )
    table_wide.columns.rename(
        "groups_observations",
        inplace=True,
    ) # single-dimensional index
    # Cluster rows in table.
    table_wide = porg.cluster_table_rows(
        table=table_wide,
    )
    # Organize indices in table.
    table_wide.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_wide.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index
    # Filter and sort columns in table.
    columns_sequence = copy.deepcopy(names_groups_observations_sequence)
    columns_sequence.insert(0, index_columns)
    table_product = porg.filter_sort_table_columns(
        table=table_wide,
        columns_sequence=columns_sequence,
        report=False,
    )
    # Translate names of features.
    table_product_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product,
            index_rows=index_columns,
            translations_columns=None,
            translations_rows=translations_features,
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Collect information.
    pail_return = dict()
    pail_return["table_summary"] = table_product
    pail_return["table_summary_translation"] = table_product_translation

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "prepare_table_signal_summaries_features_observation_groups()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_return


def prepare_tables_signals_features_sets_observations_groups(
    table=None,
    transpose_source_table=None,
    index_features=None,
    index_observations=None,
    features_selection=None,
    observations_selection=None,
    groups_features=None,
    groups_observations=None,
    names_groups_features_sequence=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    report=None,
):
    """
    Prepare derivative tables of information about measurement signal
    intensities corresponding to individual features in sets across individual
    observations in groups.

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

    ----------
    Format of source table (name: "table_source")
    ----------
    Format of source table is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------
    Notice that it is also possible for the original source table to have a
    format corresponding to the simple transposition of the format in the
    description above. In this case, it is optional to transpose the table
    before further modification.

    ----------
    Format of product table 1 (name: "table_product_1")
    ----------
    Format of product table 1 is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns. This novel product table
    shares its format with the original source table. The difference between
    them is that this novel product table only includes only a specific
    selection of features and observations from the original source table.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product table 2 (name: "table_product_2")
    ----------
    Format of product table 2 is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns. A special column gives identifiers or names corresponding to each
    observation across rows, and another special column provides names of
    categorical groups for these observations. For versatility, this table does
    not have explicitly defined indices across rows or columns.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    ----------
    Format of product tables 3, 4, and 5 (name: "table_product_3",
    "table_product_4", "table_product_5")
    ----------
    Format of product tables 3, 4, and 5 is in wide format with floating-point
    values of signal intensities for measurements corresponding to individual
    features across columns and distinct individual observations across rows. A
    special header row gives identifiers or names corresponding to each feature
    across columns. A special column gives identifiers or names corresponding
    to each observation across rows, and another special column provides names
    of categorical groups for these observations. For versatility, these tables
    do not have explicitly defined indices across rows or columns.

    Product tables 3, 4, and 5 are derivations from product table 2. This
    derivation includes standardization of values of signal intensities and
    clustering of columns for features and rows for observations by their
    respective values of signal intensities. The standardization transforms by
    z score the values of signal intensity for each feature such that these
    values have a mean of zero and standard deviation of one across all
    observations. This standardization simplifies the scale and distribution of
    values for subsequent visual representation on charts, especially heatmaps.

    The difference between product tables 3, 4, and 5 is whether the clustering
    of columns for features occurs with constrait by sets and whether the
    cluster of rows for observations occurs with constraint by groups.
    Clustering for table 3 has constraint by sets of features and by groups of
    observations. Clustering for table 4 has constraint by groups of
    observations only. Clustering for table 5 does not have any constraints.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------


    ----------
    Format of product table 6 (name: "table_product_6")
    ----------
    Format of product table 6 is in partial long format with floating-point
    values of statistics corresponding to type of descriptive statistic across
    columns and features and groups across rows.
    Product table 6 is a derivation from product table 2.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    ----------
    Format of product table 7 (name: "table_product_7")
    ----------
    Format of product table 7 is in partial long format with floating-point
    values of statistics corresponding to type of descriptive statistic across
    columns and features and groups across rows.
    Product table 7 is a derivation from product table 3.
    Product table 7 shares its format with product table 6. The difference
    between them is that product table 7 is a derivation after z score
    standardization of signals for each feature across all observations.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    ----------
    Format of product table 8 (name: "table_product_8")
    ----------
    Format of product table 8 is in wide format with floating-point values of
    a single, specific type of descriptive statistics (usually either mean or
    median) corresponding to features across rows and groups of observations
    across columns.
    Product table 8 is a derivation from product table 7.
    ----------
    group     group_1 group_2 group_3 group_4
    feature
    feature_1 0.01    0.001   0.001   0.015
    feature_2 0.01    0.001   0.001   0.015
    feature_3 0.01    0.001   0.001   0.015
    feature_4 0.01    0.001   0.001   0.015
    feature_5 0.01    0.001   0.001   0.015
    ----------

    Review: 24 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        transpose_source_table (bool): whether to transpose the original source
            table before further transformation
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_selection (list<str>): identifiers of observations for
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
    # Copy information.
    table_source = table.copy(deep=True)

    ##########
    # Optional preliminary transposition.
    # Copy information.
    table_source_format = table_source.copy(deep=True)
    # Determine whether to apply optional transposition.
    if (transpose_source_table):
        # Organize indices in table.
        table_source_format.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_source_format.set_index(
            index_features,
            append=False,
            drop=True,
            inplace=True
        )
        table_source_format.columns.rename(
            index_observations,
            inplace=True,
        ) # single-dimensional index
        # Transpose table.
        table_source_format = table_source_format.transpose(copy=True)
        # Organize indices in table.
        table_source_format.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        table_source_format.columns.rename(
            None,
            inplace=True,
        ) # single-dimensional index
        pass

    ##########
    # Organize preliminary information.
    features_all = copy.deepcopy(table_source_format.columns.to_list())
    features_all.remove(index_observations)
    observations_all = copy.deepcopy(
        table_source_format[index_observations].unique().tolist()
    )
    pail = organize_preliminary_information_to_prepare_tables_signal(
        index_features=index_features,
        index_observations=index_observations, # assigned in new tables
        features_all=features_all,
        observations_all=observations_all,
        features_selection=features_selection,
        observations_selection=observations_selection,
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

    ##########
    # Tables 1, 2, 3, 4, and 5 for individual observations.

    ##########
    # Prepare product table 1.
    # Filter specific features and observations from table.
    table_selection = porg.filter_select_table_columns_rows_by_identifiers(
        table=table_source_format,
        index_rows=pail["index_observations"],
        identifiers_columns=pail["features_available"],
        identifiers_rows=pail["observations_available"],
        report=False,
    )
    # Copy information.
    table_product_1 = table_selection.copy(deep=True)
    # Translate names of features and observations.
    table_product_1_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_1,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 2.
    # Determine and fill groups of observations.
    table_group = porg.determine_fill_table_groups_rows(
        table=table_selection,
        column_group="group",
        index_rows=pail["index_observations"],
        groups_rows=pail["groups_observations"],
        report=False,
    )
    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )
    # Copy information.
    table_product_2 = table_group.copy(deep=True)
    # Translate names of features and observations.
    table_product_2_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_2,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 3.
    # Calculate z scores of values for each feature across observations to
    # standardize their scales and distributions.
    table_scale = pscl.transform_standard_z_score_by_table_columns(
        table=table_group,
        columns=pail["features_available"],
        report=False,
    )
    # Cluster rows in table by groups of observations.
    table_product_3 = porg.cluster_table_rows_by_group(
        table=table_scale,
        index_rows=pail["index_observations"],
        column_group="group",
    )
    # Sort rows in table by groups.
    table_product_3 = porg.sort_table_rows_by_single_column_reference(
        table=table_product_3,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )
    # Cluster columns in table by groups of features.
    table_product_3 = porg.cluster_table_columns_by_external_group(
        table=table_product_3,
        indices_rows=[pail["index_observations"], "group",],
        groups_columns=pail["groups_features"],
        names_groups_sequence=pail["names_groups_features_sequence"],
        report=False,
    )
    # Translate names of features and observations.
    table_product_3_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_3,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 4.
    # Copy information.
    table_product_4 = table_scale.copy(deep=True)
    # Cluster rows in table by groups of observations.
    table_product_4 = porg.cluster_table_rows_by_group(
        table=table_product_4,
        index_rows=pail["index_observations"],
        column_group="group",
    )
    # Sort rows in table by groups.
    table_product_4 = porg.sort_table_rows_by_single_column_reference(
        table=table_product_4,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )
    # Organize indices in table.
    table_product_4.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_product_4.set_index(
        [pail["index_observations"], "group"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster columns in table.
    table_product_4 = porg.cluster_table_columns(
        table=table_product_4,
    )
    table_product_4.index = pandas.MultiIndex.from_tuples(
        table_product_4.index,
        names=[pail["index_observations"], "group"]
    )
    # Organize indices in table.
    table_product_4.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Translate names of features and observations.
    table_product_4_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_4,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 5.
    # Copy information.
    table_product_5 = table_scale.copy(deep=True)
    # Organize indices in table.
    table_product_5.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_product_5.set_index(
        [pail["index_observations"], "group"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster rows in table.
    table_product_5 = porg.cluster_table_rows(
        table=table_product_5,
    )
    table_product_5.index = pandas.MultiIndex.from_tuples(
        table_product_5.index,
        names=[pail["index_observations"], "group"]
    )
    # Cluster columns in table.
    table_product_5 = porg.cluster_table_columns(
        table=table_product_5,
    )
    table_product_5.index = pandas.MultiIndex.from_tuples(
        table_product_5.index,
        names=[pail["index_observations"], "group"]
    )
    # Organize indices in table.
    table_product_5.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Translate names of features and observations.
    table_product_5_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_5,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Tables 6, 7, and 8 for summaries using descriptive statistics on groups
    # of observations.

    ##########
    # Prepare product table 6.
    # Calculate descriptive statistics for each feature across observations.
    table_product_6 = pdesc.describe_table_features_by_groups(
        table=table_product_2_translation,
        column_group="group",
        columns_features=pail["features_available_translation"],
        index_feature=pail["index_features"],
        translations_feature=None,
        threshold_observations=5,
        digits_round=3,
        report=False,
    )

    ##########
    # Prepare product table 7.
    table_product_7 = pdesc.describe_table_features_by_groups(
        table=table_product_3_translation,
        column_group="group",
        columns_features=pail["features_available_translation"],
        index_feature=pail["index_features"],
        translations_feature=None,
        threshold_observations=5,
        digits_round=3,
        report=False,
    )

    ##########
    # Prepare product table 8.
    pail_table_8 = prepare_table_signal_summaries_features_observation_groups(
        table=table_scale,
        index_columns=pail["index_features"],
        column_group="group",
        columns_features=pail["features_available"],
        index_rows=pail["index_observations"],
        names_groups_observations_sequence=(
            pail["names_groups_observations_sequence"]
        ),
        summary="mean", # "mean", "median",
        translations_features=pail["translations_features"],
        report=report,
    )

    ##########
    # Collect information.
    pail_return = dict()
    pail_return["table_1"] = table_product_1
    pail_return["table_1_translation"] = table_product_1_translation
    pail_return["table_2"] = table_product_2
    pail_return["table_2_translation"] = table_product_2_translation
    pail_return["table_3"] = table_product_3
    pail_return["table_3_translation"] = table_product_3_translation
    pail_return["table_4"] = table_product_4
    pail_return["table_4_translation"] = table_product_4_translation
    pail_return["table_5"] = table_product_5
    pail_return["table_5_translation"] = table_product_5_translation
    pail_return["table_6"] = table_product_6
    pail_return["table_7"] = table_product_7
    pail_return["table_8"] = pail_table_8["table_summary"]
    pail_return["table_8_translation"] = (
       pail_table_8["table_summary_translation"]
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "prepare_tables_signals_individual_features_sets_observations_" +
            "groups()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_return


# Prepare derivative, deliverable product tables for correlations between
# features across observations.


# sequence_groups_features_primary <-- STILL NEED
# sequence_groups_features_secondary <-- STILL NEED


def organize_preliminary_information_to_prepare_tables_correlation(
    index_features=None,
    index_observations=None,
    features_all=None,
    observations_all=None,
    features_selection=None,
    observations_selection=None,
    groups_features=None,
    names_groups_features_sequence=None,
    features_primary=None,
    features_secondary=None,
    translations_features=None,
    report=None,
):
    """
    Organize preliminary information for subsequent use in preparation of
    derivative tables.

    Review: 16 December 2024

    arguments:
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        features_all (list<str>): identifiers of all features in the original
            source table
        observations_all (list<str>): identifiers of all observations in the
            original source table
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_selection (list<str>): identifiers of observations for
            which to include and describe values of signal intensity across
            features
        groups_features (dict<list<str>>): names of groups and identifiers
            of features that belong to each of these groups; clustering of
            features will have constraint within these groups
        names_groups_features_sequence (list<str>): names of groups for
            features in specific sequence
        features_primary (list<str>): identifiers of features to become the
            primary in pairs for correlation
        features_secondary (list<str>): identifiers of features to become the
            secondary in pairs for correlation
        translations_features (dict<str>): translations for names or
            identifiers of features
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Organize and collect information.
    pail = organize_preliminary_information_to_prepare_tables_signal(
        index_features=index_features,
        index_observations=index_observations, # assigned in new tables
        features_all=features_all,
        observations_all=observations_all,
        features_selection=features_selection,
        observations_selection=observations_selection,
        groups_features=groups_features,
        groups_observations=dict(),
        names_groups_features_sequence=names_groups_features_sequence,
        names_groups_observations_sequence=list(),
        translations_features=translations_features,
        translations_observations=None,
        report=report,
    )

    # Translate identifiers of features in sets.
    if (pail["translations_features"] is not None):
        groups_features_translation = dict()
        for group_features in pail["groups_features"].keys():
            features_group = pail["groups_features"][group_features]
            features_group_translation = list(map(
                lambda feature: (pail["translations_features"][feature]),
                features_group
            ))
            groups_features_translation[group_features] = (
                features_group_translation
            )
            pass
    else:
        groups_features_translation = copy.deepcopy(pail["groups_features"])
        pass
    pail["groups_features_translation"] = groups_features_translation

    # Copy information.
    features_primary = copy.deepcopy(features_primary)
    features_secondary = copy.deepcopy(features_secondary)

    # Prepare source information.
    # Ensure that features are unique.
    features_primary_unique = putly.collect_unique_elements(
        elements=features_primary,
    )
    features_secondary_unique = putly.collect_unique_elements(
        elements=features_secondary,
    )
    # Ensure that all features are in the source table.
    # Consider the selection of features of interest.
    pail["features_primary_available"] = list(filter(
        lambda feature: (feature in pail["features_available"]),
        features_primary_unique
    ))
    pail["features_secondary_available"] = list(filter(
        lambda feature: (feature in pail["features_available"]),
        features_secondary_unique
    ))
    # Translate identifiers of features.
    pail["features_primary_translation"] = copy.deepcopy(
        pail["features_primary_available"]
    )
    pail["features_secondary_translation"] = copy.deepcopy(
        pail["features_secondary_available"]
    )
    if (pail["translations_features"] is not None):
        pail["features_primary_translation"] = list(map(
            lambda feature: (pail["translations_features"][feature]),
            pail["features_primary_translation"]
        ))
        pail["features_secondary_translation"] = list(map(
            lambda feature: (pail["translations_features"][feature]),
            pail["features_secondary_translation"]
        ))
        # Ensure that features are unique after translation.
        pail["features_primary_translation"] = putly.collect_unique_elements(
            elements=pail["features_primary_translation"],
        )
        pail["features_secondary_translation"] = putly.collect_unique_elements(
            elements=pail["features_secondary_translation"],
        )
        pass

    # Return information.
    return pail


def prepare_tables_correlations_of_features_across_observations(
    table=None,
    index_features=None,
    index_observations=None,
    translations_features=None,
    features_selection=None,
    observations_selection=None,
    features_primary=None,
    features_secondary=None,
    groups_features_primary=None,
    groups_features_secondary=None,
    names_groups_sequence_features_primary=None,
    names_groups_sequence_features_secondary=None,
    method_priority=None,
    report=None,
):
    """
    Prepare derivative tables of information about correlations between
    features across observations.

    Any translation of identifiers or names of features occurs before the
    selection of features and observations. This sequence of actions avoids
    potential problems from translations to overlapping, redundant identifiers
    or names.

    By intentional design, this function does not apply any transformation of
    scale or normalization of distribution to the values of signal intensity.

    ----------
    Format of source table (name: "table_source")
    ----------
    Format of source table is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product table 1 (name: "table_product_1")
    ----------
    Format of product table 1 is in long format with information about pairs of
    features and calculations of their correlations across observations.
    ----------
    feature_1 feature_2   observations correlation_pearson probability_pears...

    name_1    name_2      100          0.1                 0.01             ...
    name_1    name_2      200          0.1                 0.01             ...
    name_1    name_2      300          0.1                 0.01             ...
    name_1    name_2      400          0.1                 0.01             ...
    name_1    name_2      500          0.1                 0.01             ...
    ----------

    ----------
    Format of product table 2 (name: "table_product_2")
    ----------
    Format of product table 2 is in wide format with information about a single
    type of correlation coefficient between features across columns and
    features across rows.
    ----------
    feature_2   name_2 name_2 name_2 name_2 name_2
    feature_1
    name_1      0.1    0.1    0.1    0.1    0.1
    name_1      0.1    0.1    0.1    0.1    0.1
    name_1      0.1    0.1    0.1    0.1    0.1
    name_1      0.1    0.1    0.1    0.1    0.1
    name_1      0.1    0.1    0.1    0.1    0.1
    ----------

    Review: 17 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table

        translations_features (dict<str>): translations for names or
            identifiers of features


        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_selection (list<str>): identifiers of observations for
            which to include and describe values of signal intensity across
            features
        groups_features (dict<list<str>>): names of groups and identifiers
            of features that belong to each of these groups; clustering of
            features will have constraint within these groups


        names_groups_sequence_features_primary (list<str>): names of groups for
            features in specific sequence
        features_primary (list<str>): identifiers of features to become the
            primary in pairs for correlation
        features_secondary (list<str>): identifiers of features to become the
            secondary in pairs for correlation
        method_priority (str): name of method of correlation to prioritize in
            preparation of the product table in wide format, either 'pearson',
            'spearman', or 'kendall'
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    table_source = table.copy(deep=True)

    ##########
    # Organize preliminary information.
    features_all = copy.deepcopy(table_source.columns.to_list())
    features_all.remove(index_observations)
    observations_all = copy.deepcopy(
        table_source[index_observations].unique().tolist()
    )
    pail = organize_preliminary_information_to_prepare_tables_correlation(
        index_features=index_features,
        index_observations=index_observations, # assigned in new tables
        features_all=features_all,
        observations_all=observations_all,
        features_selection=features_selection,
        observations_selection=observations_selection,
        groups_features=groups_features,
        #groups_observations=dict(),
        names_groups_features_sequence=names_groups_features_sequence,
        #names_groups_observations_sequence=list(),
        features_primary=features_primary,
        features_secondary=features_secondary,
        translations_features=translations_features,
        #translations_observations=translations_observations,
        report=report,
    )

    # Translate names of features.
    table_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_source,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=None,
            remove_redundancy=True,
            report=False,
    ))
    # Filter specific features and observations from table.
    table_selection = porg.filter_select_table_columns_rows_by_identifiers(
        table=table_translation,
        index_rows=pail["index_observations"],
        identifiers_columns=pail["features_available_translation"],
        identifiers_rows=pail["observations_available_translation"],
        report=False,
    )
    # Convert values in columns to float type.
    features_correlation_all = list()
    features_correlation_all.extend(copy.deepcopy(
        pail["features_primary_translation"]
    ))
    features_correlation_all.extend(copy.deepcopy(
        pail["features_secondary_translation"]
    ))
    table_source = putly.convert_table_columns_variables_types_float(
        columns=features_correlation_all,
        table=table_selection,
    )

    ##########
    # Prepare product table 1.
    # Correlations between pairs of features in long format.
    # Collect information.
    records = list()
    columns_sequence = None
    # Iterate on pairwise combinations of primary and secondary features.
    for feature_primary in pail["features_primary_translation"]:
        for feature_secondary in pail["features_secondary_translation"]:
            # Calculate correlations.
            pail_correlation = pdesc.calculate_correlations_table_columns_pair(
                table=table_selection,
                column_primary=feature_primary,
                column_secondary=feature_secondary,
                count_minimum_observations=10,
                report=False,
            )
            # Collect information.
            record = dict()
            record["feature_primary"] = feature_primary
            record["feature_secondary"] = feature_secondary
            if (columns_sequence is None):
                columns_sequence = copy.deepcopy(pail_correlation["names"])
                pass
            del pail_correlation["names"]
            record.update(pail_correlation)
            records.append(record)
            pass
        pass
    # Organize table.
    table_correlation_long = pandas.DataFrame(data=records)
    # Define sequence of columns in table.
    columns_sequence.insert(0, "feature_secondary",)
    columns_sequence.insert(0, "feature_primary",)
    # Filter and sort columns in table.
    table_correlation_long = porg.filter_sort_table_columns(
        table=table_correlation_long,
        columns_sequence=columns_sequence,
        report=False,
    )
    # Copy information.
    table_product_1 = table_correlation_long.copy(deep=True)

    ##########
    # Prepare product table 2.
    # Correlations between pairs of features in wide format.
    # Copy information.
    table_correlation_long_temporary = table_correlation_long.copy(deep=True)
    # Arrange table from full long format to wide format.
    if (method_priority == "pearson"):
        table_correlation_wide = table_correlation_long_temporary.pivot(
            columns=["feature_primary",],
            index="feature_secondary",
            values="correlation_pearson",
        )
    elif (method_priority == "spearman"):
        table_correlation_wide = table_correlation_long_temporary.pivot(
            columns=["feature_primary",],
            index="feature_secondary",
            values="correlation_spearman",
        )
        pass
    # Organize indices in table.
    table_correlation_wide.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Determine and fill groups of features across rows.
    table_correlation_wide = porg.determine_fill_table_groups_rows(
        table=table_correlation_wide,
        column_group="group",
        index_rows="feature_secondary",
        groups_rows=pail["groups_features_secondary"],
        report=False,
    )
    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )

    # Copy information.
    table_wide_cluster_group = table_correlation_wide.copy(deep=True)
    # Add a group column to satisfy requirement for an index with multiple
    # levels across rows.
    table_wide_cluster_group["group"] = "group"
    # Cluster columns in table by groups of features.
    table_wide_cluster_group = porg.cluster_table_columns_by_external_group(
        table=table_wide_cluster_group,
        indices_rows=["feature_secondary", "group",],
        groups_columns=pail["groups_features_translation"],
        names_groups_sequence=pail["names_groups_features_sequence"],
        report=False,
    )
    # Remove unnecessary columns.
    table_wide_cluster_group.drop(
        labels=["group",],
        axis="columns",
        inplace=True
    )
    # Copy information.
    table_product_2 = table_wide_cluster_group.copy(deep=True)

    ##########
    # Prepare product table 3.
    # Correlations between pairs of features in wide format.
    # Copy information.
    table_wide_cluster = table_correlation_wide.copy(deep=True)
    # Organize indices in table.
    table_wide_cluster.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_wide_cluster.set_index(
        ["feature_secondary",],
        append=False,
        drop=True,
        inplace=True,
    )
    table_wide_cluster.columns.rename(
        "feature_primary",
        inplace=True,
    ) # single-dimensional index
    # Cluster columns in table.
    table_wide_cluster = porg.cluster_table_columns(
        table=table_wide_cluster,
    )
    # Organize indices in table.
    table_wide_cluster.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_wide_cluster.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index
    # Copy information.
    table_product_3 = table_wide_cluster.copy(deep=True)

    print("!!!!!!!!!!!!!!!!!!!!!table_correlation_wide")
    print(table_correlation_wide)


    ##########
    # Collect information.
    pail_return = dict()
    pail_return["table_1"] = table_product_1
    pail_return["table_2"] = table_product_2
    pail_return["table_3"] = table_product_3

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "prepare_tables_correlations_of_features_across_observations()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_return


# Read and organize sets of features.


def read_extract_set_features(
    name_set=None,
    suffix_file=None,
    path_directory=None,
    report=None,
):
    """
    Read and extract from a source file the identifiers of features in a set.

    arguments:
        name_set (str): name for a set of features that corresponds to the name
            of a file in the source directory
        suffix_file (str): suffix in name of file
        path_directory (str): path to directory within which to find files of
            sets of features
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of features in a set from file

    """

    # Determine whether parameter is valid.
    if (
        (name_set is not None) and
        (str(name_set).strip().lower() != "none")
    ):
        # Define paths to files.
        path_file = os.path.join(
            path_directory, str(name_set + suffix_file),
        )
        # Read information from file.
        features_set = putly.read_file_text_list(
            delimiter="\n",
            path_file=path_file,
        )
    else:
        features_set = list()
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = "read_extract_set_features"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=4)
        count_items = len(features_set)
        print("count of items in set or list: " + str(count_items))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return features_set


def read_extract_combine_custom_feature_sets(
    names_sets=None,
    features_available=None,
    path_directory=None,
    report=None,
):
    """
    Reads and extracts from a source file the identifiers of genes in a set.

    Review: TCW; 23 December 2024

    arguments:
        names_sets (dict<list<str>>): names of sets of genes to include within
            custom combination sets with custom names
        features_available (list<str>): identifiers of features with available
            signals
        path_directory (str): path to directory from which to find and read
            files for sets of features
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): identifiers of genes within custom combination sets

    """

    # Copy information.
    names_sets = copy.deepcopy(names_sets)
    features_available = copy.deepcopy(features_available)
    # Iterate on features and values for sets of genes.
    if (
        (names_sets is not None) and
        (str(names_sets).strip().lower() != "none")
    ):
        sets_features = dict()
        for name_merge in names_sets.keys():
            features_set_collection = list()
            for name_set in names_sets[name_merge]:
                if (
                    (name_set is not None) and
                    (str(name_set).strip().lower() != "none")
                ):
                    # Read and extract identifiers of features in set.
                    features_set_temporary = read_extract_set_features(
                        name_set=name_set,
                        suffix_file=".txt",
                        path_directory=path_directory,
                        report=False,
                    )
                    # Collect information.
                    features_set_collection.extend(features_set_temporary)
                    pass
                pass
            # Organize custom combination set of genes.
            # Collect unique names of genes in set.
            features_set_collection = putly.collect_unique_elements(
                elements=features_set_collection,
            )
            # Ensure that all genes in set are in the table of signals.
            if (features_available is not None):
                features_set_collection = list(filter(
                    lambda feature: (feature in features_available),
                    features_set_collection
                ))
                pass
            # Collect genes in custom set.
            sets_features[name_merge] = copy.deepcopy(features_set_collection)
            pass
        pass
    else:
        sets_features = None
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = "read_extract_combine_custom_feature_sets"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=4)
        # Iterate on custom sets of genes.
        if (sets_features is not None):
            for name_set in sets_features.keys():
                count_set = len(sets_features[name_set])
                print(
                    "count of items in set " +
                    str(name_set) + ": " +
                    str(count_set)
                )
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return sets_features


def read_organize_feature_sets(
    names_sets_inclusion=None,
    names_sets_cluster=None,
    names_sets_allocation=None,
    features_available=None,
    directory_child=None,
    path_directory_parent=None,
    report=None,
):
    """
    Prepare sets of genes for analysis and description.

    arguments:
        names_sets_inclusion (dict<list<str>>): names of sets of genes to
            include within a single custom combination set with custom name;
            this is the set of total genes for inclusion in description and
            visual representation
        names_sets_cluster (dict<list<str>>): names of sets of genes to
            include within custom combination sets with custom names; these are
            the sets or groups within which to cluster signals for individual
            genes across samples; the sets must include all of the total
            inclusion genes, and each set must be mutually exclusive
        names_sets_allocation (dict<list<str>>): names of sets of genes to
            include within custom combination sets with custom names; these are
            the sets or groups within which to allocate the total inclusion
            genes for visual representation
        features_available (list<str>): identifiers of features with available
            signals
        directory_child (str): name of child directory from which to read
            sets of features
        path_directory_parent : (str): path to parent directory within which to
            find child directory from which to read files for sets of features
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    ##########
    # Copy information.
    names_sets_inclusion = copy.deepcopy(names_sets_inclusion)
    names_sets_cluster = copy.deepcopy(names_sets_cluster)
    names_sets_allocation = copy.deepcopy(names_sets_allocation)
    features_available = copy.deepcopy(features_available)

    # Define paths to child files.
    path_directory_sets_gene = os.path.join(
        path_directory_parent, directory_child,
    )

    ##########
    # Genes for inclusion.
    sets_inclusion = read_extract_combine_custom_feature_sets(
        names_sets=names_sets_inclusion,
        features_available=features_available,
        path_directory=path_directory_sets_gene,
        report=report,
    )
    # Determine whether set "main" already exists in the sets for inclusion.
    # If set "main" does not already exist, then define set "main" as the
    # unique union of all other sets.
    if ("main" not in sets_inclusion.keys()):
        sets_inclusion["main"] = list()
        for name_set in sets_inclusion.keys():
            genes_set = copy.deepcopy(sets_inclusion[name_set])
            sets_inclusion["main"].extend(genes_set)
            pass
        # Collect unique names of genes in set.
        sets_inclusion["main"] = putly.collect_unique_elements(
            elements=sets_inclusion["main"],
        )
        # Ensure that all genes in set are in the table of signals.
        sets_inclusion["main"] = list(filter(
            lambda gene: (gene in features_available),
            sets_inclusion["main"]
        ))
        pass
    # From the sets of genes for inclusion, only use the combination custom set
    # with name "main".

    ##########
    # Genes for clustering.
    sets_cluster = read_extract_combine_custom_feature_sets(
        names_sets=names_sets_cluster,
        features_available=features_available,
        path_directory=path_directory_sets_gene,
        report=report,
    )
    # Determine whether there were any valid sets of genes for constraint on
    # the clustering operation.
    # If there were not any sets of genes for clustering constraint, then use
    # set "main" from the sets for inclusion.
    # Notice that any sets for constraint on clustering operation must be
    # mutually exclusive while being all inclusive of the available genes.
    if (sets_cluster is None):
        sets_cluster = dict()
        sets_cluster["main"] = copy.deepcopy(sets_inclusion["main"])
        pass

    ##########
    # Genes for allocation.
    sets_allocation = read_extract_combine_custom_feature_sets(
        names_sets=names_sets_allocation,
        features_available=features_available,
        path_directory=path_directory_sets_gene,
        report=report,
    )

    # Collect information.
    pail = dict()
    pail["sets_inclusion"] = sets_inclusion
    pail["sets_cluster"] = sets_cluster
    pail["sets_allocation"] = sets_allocation
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "read_organize_feature_sets()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("use special combination set 'main' for total inclusion")
        count_inclusion = (len(pail["sets_inclusion"]["main"]))
        print(
            "count of features for total inclusion: " + str(count_inclusion)
        )
        putly.print_terminal_partition(level=5)

        pass
    # Return information.
    return pail


def combine_sets_items_union_unique(
    sets_items=None,
    report=None,
):
    """
    Combines items from multiple sets by taking the union and then collecting
    unique items.

    Review: TCW; 11 February 2025

    arguments:
        sets_items (list<list<str>>): lists of items in distinct sets
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): items in combination set

    """

    # Copy information.
    sets_items = copy.deepcopy(sets_items)
    # Collect.
    items_combination = list()
    # Iterate.
    for items_set in sets_items:
        # Combine by union.
        items_combination.extend(items_set)
        pass
    # Collect unique.
    items_combination_unique = putly.collect_unique_elements(
        elements=items_combination,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = "combine_sets_items_union_unique"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=4)
        # Summarize.
        count_items = len(items_combination_unique)
        print(
            "count of items in combination set: " + str(count_items)
        )
        putly.print_terminal_partition(level=5)
        pass
    # Return.
    return items_combination_unique


def combine_sets_items_difference_unique(
    items_inclusion=None,
    items_exclusion=None,
    report=None,
):
    """
    Combines items from multiple sets by taking the difference and then
    collecting unique items.

    Review: TCW; 11 February 2025

    arguments:
        items_inclusion (list<str>): items to include, main set
        items_exclusion (list<str>): items to exclude from main, inclusion set
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): items in combination set

    """

    # Copy information.
    items_inclusion = copy.deepcopy(items_inclusion)
    items_exclusion = copy.deepcopy(items_exclusion)
    # Filter.
    items_combination = list(filter(
        lambda item: (item not in items_exclusion),
        items_inclusion
    ))
    # Collect unique.
    items_combination_unique = putly.collect_unique_elements(
        elements=items_combination,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = "combine_sets_items_difference_unique"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=4)
        # Summarize.
        count_items = len(items_combination_unique)
        print(
            "count of items in combination set: " + str(count_items)
        )
        putly.print_terminal_partition(level=5)
        pass
    # Return.
    return items_combination_unique


def combine_sets_items_intersection_unique(
    items_first=None,
    items_second=None,
    report=None,
):
    """
    Combines items from multiple sets by taking the intersection and then
    collecting unique items.

    Review: TCW; 11 February 2025

    arguments:
        items_first (list<str>): items to include in intersection
        items_second (list<str>): items to include in intersection
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): items in combination set

    """

    # Copy information.
    items_first = copy.deepcopy(items_first)
    items_second = copy.deepcopy(items_second)
    # Combine items by union.
    items_union = combine_sets_items_union_unique(
        sets_items=[
            items_first,
            items_second,
        ],
        report=False,
    )
    # Filter.
    items_combination = list(filter(
        lambda item: ((item in items_first) and (item in items_second)),
        items_union
    ))
    # Collect unique.
    items_combination_unique = putly.collect_unique_elements(
        elements=items_combination,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = "combine_sets_items_intersection_unique"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=4)
        # Summarize.
        count_items = len(items_combination_unique)
        print(
            "count of items in combination set: " + str(count_items)
        )
        putly.print_terminal_partition(level=5)
        pass
    # Return.
    return items_combination_unique



# Prepare table for allocation of features to sets.


def assemble_table_features_sets_allocation(
    features_selection=None,
    groups_features=None,
    names_groups_features_sequence=None,
    group_other=None,
    index_features=None,
    report=None,
):
    """
    Assemble table with logical binary designations of allocation for a
    selection of features to sets.

    ----------
    Format of product table
    ----------
    sets        set_1 set_2 set_3
    features
    feature_1   1     0     0
    feature_2   0     1     0
    feature_3   0     0     1
    feature_4   0     1     0
    feature_5   1     0     0
    ----------

    Review: 12 December 2024

    arguments:
        features_selection (list<str>): identifiers of features to allocate to
            sets
        groups_features (dict<list<str>>): names of groups and identifiers
            of features that belong to each of these groups
        names_groups_features_sequence (list<str>): names of groups for
            features in specific sequence
        group_other (str): name for group to which to assign all features, even
            without allocation to any other groups
        index_features (str): name for index corresponding to features across
            rows in the allocation table
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information.
    features_selection = copy.deepcopy(features_selection)
    groups_features = copy.deepcopy(groups_features)
    names_groups_features_sequence = copy.deepcopy(
        names_groups_features_sequence
    )

    # Determine whether to include column for other sets.
    if (
        (group_other is not None) and
        (len(group_other) > 0)
    ):
        names_groups_features_sequence.append(group_other)
        pass

    ##########
    # Determine set allocations for each relevant gene.
    # Iterate on relevant genes.
    # Collect information about set allocations for each gene.
    records = list()
    for feature in features_selection:
        # Collect information.
        record = dict()
        record[index_features] = feature
        if (
            (group_other is not None) and
            (len(group_other) > 0)
        ):
            record[group_other] = 1
        # Iterate on sets for allocation.
        for name_set in groups_features.keys():
            if (feature in groups_features[name_set]):
                record[name_set] = 1
                if (
                    (group_other is not None) and
                    (len(group_other) > 0)
                ):
                    record[group_other] = 0
            else:
                record[name_set] = 0
            pass
        # Collect information.
        records.append(record)
        pass
    # Create pandas data-frame table.
    table = pandas.DataFrame(data=records)

    # Filter and sort columns within table.
    columns_sequence = copy.deepcopy(names_groups_features_sequence)
    columns_sequence.insert(0, index_features)
    if (
        (group_other is not None) and
        (len(group_other) > 0)
    ):
        columns_sequence.append("other")
        pass
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=False,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "assemble_table_features_sets_allocation()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("table with allocations of features to sets")
        print(table)
        # Determine whether to describe column for other sets.
        if (
            (group_other is not None) and
            (len(group_other) > 0)
        ):
            # Filter rows within table.
            table_nonother = table.loc[
                (table[group_other] == 0), :
            ].copy(deep=True)
            putly.print_terminal_partition(level=5)
            print(table_nonother)
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def sort_check_table_rows_sequence_custom(
    table=None,
    index_rows=None,
    values_sort_sequence=None,
    report=None,
):
    """
    Sort rows in a table by a custom sequence of values corresponding to the
    table's index across rows. Check that the table's rows after sort match the
    custom sequence.

    arguments:
        table (object): Pandas data-frame table
        index_rows (str): name for index across rows by which to sort rows
        values_sort_sequence (list<str>): values corresponding to index across
            rows in custom sequence by which to sort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    values_sort_sequence = copy.deepcopy(values_sort_sequence)

    # Prepare indices for sort.
    reference_sort = dict(zip(
        values_sort_sequence, range(len(values_sort_sequence))
    ))

    # Sort rows in table reference indices.
    table_sort = porg.sort_table_rows_by_single_column_reference(
        table=table,
        index_rows=index_rows,
        column_reference=index_rows,
        column_sort_temporary="sort_temporary_982416",
        reference_sort=reference_sort,
    )

    # Report.
    if report:
        # Extract index across rows for comparison.
        values_index_rows = copy.deepcopy(
            table_sort[index_rows].to_list()
        )
        # Count identifiers of features.
        count_sequence = len(values_sort_sequence)
        count_sort = len(values_index_rows)
        # Confirm that sets of indices from both sources are inclusive.
        check_inclusion = putly.compare_lists_by_mutual_inclusion(
            list_primary=values_sort_sequence,
            list_secondary=values_index_rows,
        )
        # Confirm that sets of indices from both sources are identical across
        # their respective sequences.
        check_identity = putly.compare_lists_by_elemental_identity(
            list_primary=values_sort_sequence,
            list_secondary=values_index_rows,
        )
        # Confirm that sets of indices from both sources are equal.
        check_equality = (
            values_sort_sequence == values_index_rows
        )

        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "sort_check_table_rows_sequence_custom()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)
        print(
            "count of values in sort sequence: " +
            str(count_sequence)
        )
        print(
            "count of values of index in table: " +
            str(count_sort)
        )
        putly.print_terminal_partition(level=5)
        print(
            "confirm that values of index across rows are identical to " +
            "reference"
        )
        print("check inclusion: " + str(check_inclusion))
        print("check identity: " + str(check_identity))
        print("check equality: " + str(check_equality))
        putly.print_terminal_partition(level=5)
        print("product table after sort")
        print(table_sort)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_sort


def prepare_table_features_sets_allocation_match_table_signal(
    table_signal=None,
    index_features=None,
    indices_observations=None,
    groups_features=None,
    names_groups_features_sequence=None,
    translations_features=None,
    report=None,
):
    """
    Prepare tables for allocation of features to sets that match with a table
    of signals.

    Optional translation of the identifiers or names of features occurs before
    assembly of the allocation table and before comparison and sort to match
    the signal table.

    ----------
    Format of source table for signals (name: "table_signal")
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity corresponding to features across columns and observations
            across rows
        index_features (str): name for index corresponding to features across
            columns in the signal table and across rows in the allocation table
        indices_observations (str): names for indices corresponding to
            observations across rows in the signal table
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        groups_features (dict<list<str>>): names of groups and identifiers
            of features that belong to each of these groups
        names_groups_features_sequence (list<str>): names of groups for
            features in specific sequence
        translations_features (dict<str>): translations for names or
            identifiers of features
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    table_signal = table_signal.copy(deep=True)
    groups_features = copy.deepcopy(groups_features)
    names_groups_features_sequence = copy.deepcopy(
        names_groups_features_sequence
    )

    # Extract identifiers of features from the columns in the table of signals.
    features_signal = copy.deepcopy(table_signal.columns.to_list())
    for index in indices_observations:
        features_signal.remove(index)
        pass
    features_signal_unique = putly.collect_unique_elements(
        elements=features_signal,
    )

    # Determine whether to translate identifiers of features in sets to match
    # the identifiers of features from the table of signals.
    if (translations_features is not None):
        groups_features_translation = dict()
        for group_features in groups_features.keys():
            features_group = groups_features[group_features]
            features_group_translation = list(map(
                lambda feature: (translations_features[feature]),
                features_group
            ))
            groups_features_translation[group_features] = (
                features_group_translation
            )
            pass
    else:
        groups_features_translation = copy.deepcopy(groups_features)
        pass

    # Assemble table with logical binary designations of allocation for a
    # selection of features to sets.
    table_allocation = assemble_table_features_sets_allocation(
        features_selection=features_signal_unique,
        groups_features=groups_features_translation,
        names_groups_features_sequence=names_groups_features_sequence,
        group_other=None, # None, "", or "other"
        index_features=index_features,
        report=report,
    )
    # Sort the sequence of features across rows in the set allocation table to
    # match the sequence of features across columns in the signal table.
    # Extract identifiers of features in original sequence from the columns in
    # the table of signals.
    features_sequence = copy.deepcopy(table_signal.columns.to_list())
    for index in indices_observations:
        features_sequence.remove(index)
        pass
    table_allocation_sort = sort_check_table_rows_sequence_custom(
        table=table_allocation,
        index_rows=index_features,
        values_sort_sequence=features_sequence,
        report=report,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "prepare_table_features_sets_allocation_match_table_signal()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_allocation_sort


# Plot charts.


def plot_heatmap_features_sets_observations_groups(
    table_signal=None,
    table_feature=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 27 December 2024

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
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    table_signal_extract = table_signal.copy(deep=True)
    table_feature = table_feature.copy(deep=True)
    # Organize indices in table.
    table_signal_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_signal_extract.set_index(
        [index_rows, column_group,],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_signal_extract.to_numpy())
    #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
    round_offset = (numpy.nanmin(matrix) * 0.10)
    value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
    value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heatmap_signal_features_sets_observations_groups(
        table_signal=table_signal,
        table_feature_sets=table_feature,
        format_table_signal=2, # 2: features in columns; observations, groups in rows
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=True,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar="individual signal (z-score)",
        size_title_ordinate="ten",
        size_title_abscissa="ten",
        size_title_bar="thirteen",
        size_label_feature_set="fifteen",
        size_label_legend_observation_group="fourteen",
        size_label_bar="fifteen",
        show_scale_bar=True,
        aspect="portrait",
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Return information.
    return figure


def plot_heatmap_features_sets_observations_labels(
    table_signal=None,
    table_feature=None,
    index_columns=None,
    index_rows=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 26 December 2024

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
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    table_signal_extract = table_signal.copy(deep=True)
    table_feature = table_feature.copy(deep=True)
    # Organize indices in table.
    table_signal_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_signal_extract.set_index(
        [index_rows,],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_signal_extract.to_numpy())
    #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
    round_offset = (numpy.nanmin(matrix) * 0.10)
    value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
    value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heatmap_signal_features_sets_observations_labels(
        table_signal=table_signal,
        table_feature_sets=table_feature,
        format_table_signal=1, # 1: features in rows, observations or groups in columns
        index_columns=index_columns,
        index_rows=index_rows,
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar="Mean Signal (z-score)",
        labels_abscissa_categories=None,
        size_title_ordinate="eight",
        size_title_abscissa="eight",
        size_title_bar="nine",
        size_label_feature_set="eleven",
        size_label_abscissa="nine",
        size_label_bar="eleven",
        show_labels_abscissa=True,
        show_scale_bar=True, # whether to show scale bar on individual figures
        aspect="square", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Return information.
    return figure


def plot_heatmap_features_observations_labels(
    table=None,
    index_columns=None,
    index_rows=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 27 December 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values of a
            single, specific type of descriptive statistics (usually either
            mean or median) corresponding to groups of observations across
            columns and features across rows
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)
    table_extract = table.copy(deep=True)
    # Organize indices in table.
    table_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_extract.set_index(
        [index_rows],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_extract.to_numpy())
    #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
    round_offset = (numpy.nanmin(matrix) * 0.10)
    value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
    value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heatmap_signal_features_observations_labels(
        table=table,
        format_table=1, # 1: features in rows, observations or groups in columns
        index_columns=index_columns,
        index_rows=index_rows,
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar="Mean Signal (z-score)",
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        size_title_ordinate="eight",
        size_title_abscissa="eight",
        size_title_bar="nine",
        size_label_ordinate="fifteen", # determine automatically if "None"; "fifteen"
        size_label_abscissa="nine", # determine automatically if "None"
        size_label_bar="eleven",
        show_labels_ordinate=True,
        show_labels_abscissa=True,
        show_scale_bar=True,
        aspect="portrait", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Return information.
    return figure


def manage_plot_charts_demonstration(
    paths=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:


    """

    ##########
    # 1. Read source information from file.
    pail_source = read_source_demonstration(
        paths=paths,
        report=report,
    )
    # Copy information.
    table_signal = pail_source["table_demonstration_signal"].copy(deep=True)
    table_set_signal = pail_source["table_demonstration_set"].copy(deep=True)
    table_set_summary = pail_source["table_demonstration_set"].copy(deep=True)

    ##########
    # 2. Heatmap Individual
    figure_heatmap_individual = (
        plot_heatmap_features_sets_observations_groups(
            table_signal=table_signal,
            table_feature=table_set_signal,
            index_columns="feature",
            index_rows="observation",
            column_group="group",
            report=report,
    ))

    ##########
    # Heatmap Mean
    # Prepare table.
    columns_features = list()
    for index in range(30):
        column = str("feature_" + str(index + 1))
        columns_features.append(column)
        pass
    # Calculate z scores of values for each feature across observations to
    # standardize their scales and distributions.
    table_scale = pscl.transform_standard_z_score_by_table_columns(
        table=table_signal,
        columns=columns_features,
        report=False,
    )
    # Prepare table of summary descriptive statistics.
    pail_table = prepare_table_signal_summaries_features_observation_groups(
        table=table_scale,
        index_columns="feature",
        column_group="group",
        columns_features=columns_features,
        index_rows="observation",
        names_groups_observations_sequence=[
            "group_1", "group_2", "group_3", "group_4",
        ],
        summary="mean", # "mean", "median",
        translations_features=None,
        report=report,
    )
    # Prepare table for allocation of features to sets.
    # Sort the sequence of features across rows in the set allocation table to
    # match the sequence of features across columns in the signal table.
    features_sequence = copy.deepcopy(
        pail_table["table_summary"]["feature"].tolist()
    )
    table_set_summary = sort_check_table_rows_sequence_custom(
        table=table_set_summary,
        index_rows="feature",
        values_sort_sequence=features_sequence,
        report=report,
    )

    # Create heatmaps.
    figure_heatmap_mean_set = plot_heatmap_features_sets_observations_labels(
        table_signal=pail_table["table_summary"],
        table_feature=table_set_summary,
        index_columns="group_observations",
        index_rows="feature",
        report=False,
    )
    figure_heatmap_mean_label = plot_heatmap_features_observations_labels(
        table=pail_table["table_summary_translation"],
        index_columns="group_observations",
        index_rows="feature",
        report=False,
    )

    ##########
    # Collect information.
    pail_write_plot = dict()
    pail_write_plot["heatmap_individual"] = figure_heatmap_individual
    pail_write_plot["heatmap_mean_set"] = figure_heatmap_mean_set
    pail_write_plot["heatmap_mean_label"] = figure_heatmap_mean_label

    ##########
    # _. Write product information to file.
    #paths["out_data"]
    #paths["out_plot"]

    # Define paths to directories.
    path_directory_plot = os.path.join(
        paths["out_procedure_plot"], "demonstration",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_plot,
    )
    # Write figures to file.
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot,
        format="jpg", # jpg, png, svg
        resolution=300,
        path_directory=path_directory_plot,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "manage_plot_charts_demonstration()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)

    # Return information.
    return pail_write_plot




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
    procedure="organize_subject"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
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
        initialize_routine=True,
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
    # 3. Demonstrate and test functions to create plot charts.
    if False:
        manage_plot_charts_demonstration(
            paths=paths,
            report=report,
        )

    ##########
    # 3. Organize table of properties for study subjects.
    columns_original = pail_source["columns_all"]
    columns_novel = define_sequence_columns_novel_subject_feature()
    pail_organization = organize_table_subject_property(
        table=pail_source["table_subject_property"],
        translations_column=pail_source["translations_feature_forward"],
        columns_original=columns_original,
        columns_novel=columns_novel,
        report=report,
    )

    ##########
    # 5. Calculate principal components of features quantitative, continuous,
    # ratio or interval scales of measurement.

    #pail_source["columns_quantitative"]
    #pail_source["columns_olink_plasma"] # only available for second visit
    #pail_source["columns_olink_muscle"] # only available for second visit
    #pail_source["columns_olink_adipose"] # only available for second visit



    ##########
    # 6. Collect information.
    # Collections of files.

    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_subject")] = pail_organization["table"]
    pail_write_objects = dict()
    #pail_write_objects[str("samples")]

    ##########
    # 7. Write product information to file.
    if False:
        putly.write_lists_to_file_text(
            pail_write=pail_write_lists,
            path_directory=paths["out_procedure_lists"],
            delimiter="\n",
        )
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
            path_directory=paths["out_procedure_data"],
        )

    pass


###############################################################################
# End
