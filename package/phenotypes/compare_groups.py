"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'compare_groups' is part of the 'phenotypes' subpackage
within the 'age_exercise' package.

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
import os.path
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
from warnings import simplefilter
simplefilter(action="ignore", category=pandas.errors.PerformanceWarning)

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.regression as preg
import partner.plot as pplot
import partner.parallelization as prall
import age_exercise.phenotypes.organize_subject as aexph_sub
import age_exercise.phenotypes.organize_sample as exph_sample

###############################################################################
# Functionality

# TODO; TCW; 31 March 2025
# Move to this module the portion of "organize_subject.py" about comparing
# quantitative features between groups of observations.
# Adapt functionality from module "transcriptomics.compare_sets_groups".


##########
# 1. Initialize directories for read of source and write of product files.


##########
# 2. Read source information from file.


def define_type_table_columns_parameter():
    """
    Defines the types of variables for columns in table.

    Review: TCW; 10 April 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    types_columns["execution"] = "string" # "int32"
    types_columns["sequence"] = "string" # "int32"
    types_columns["group"] = "string"
    types_columns["name"] = "string"
    #types_columns["name_combination"] = "string"
    types_columns["selection_observations"] = "string"
    types_columns["review"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


def read_organize_source_parameter(
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
    #paths["in_demonstration"]
    #paths["in_parameters"]
    #paths["in_parameters_private"]

    # Define paths to child files.
    path_file_table_parameter = os.path.join(
        paths["in_parameters_private"], "age_exercise", "phenotypes",
        "table_groups_observations.tsv",
    )

    # Read information from file.

    # Table of parameters for instances.
    types_columns = define_type_table_columns_parameter()
    table = pandas.read_csv(
        path_file_table_parameter,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    table["execution"] = pandas.to_numeric(
        table["execution"],
        downcast="integer",
        errors="coerce",
    )
    table["sequence"] = pandas.to_numeric(
        table["sequence"],
        downcast="integer",
        errors="coerce",
    )

    # Collect information.
    records = list()
    for index, row in table.iterrows():
        # Collect information and parameters from current row in table.
        record = dict()
        record["execution"] = int(row["execution"])
        record["sequence"] = row["sequence"]
        record["group"] = str(row["group"]).strip()
        record["name"] = str(row["name"]).strip() # name for instance of parameters
        record["name_combination"] = "_".join([
            str(row["group"]).strip(),
            str(row["sequence"]).strip(),
            str(row["name"]).strip(),
        ])
        record["selection_observations"] = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=row["selection_observations"],
            )
        )["features_values"]
        record["review"] = str(row["review"]).strip()
        record["note"] = str(row["note"]).strip()

        # Collect unique names of columns relevant to instance of parameters
        # from current row in table.
        features_relevant = list()
        dictionaries = [
            "selection_observations",
        ]
        for dictionary in dictionaries:
            if record[dictionary] is not None:
                features_relevant.extend(list(record[dictionary].keys()))
                pass
            pass
        #features_relevant.extend(record["any_others"])
        features_relevant = putly.collect_unique_elements(
            elements=features_relevant,
        )
        record["features_relevant"] = copy.deepcopy(features_relevant)
        # Collect information and parameters for current row in table.
        records.append(record)
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table.copy(deep=True)
    pail["records"] = copy.deepcopy(records)

    # Report.
    if report:
        # Organize information.
        count_records = len(records)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
        print("function: read_organize_source_parameter()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        print("count of records or instances: " + str(count_records))
        putly.print_terminal_partition(level=5)
        print("instance[0]:")
        print(records[0])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def define_type_table_columns_subject_sample():
    """
    Defines the types of variables for columns in table.

    Review: TCW; 10 April 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()

    types_columns["subject_visit"] = "string"
    types_columns["identifier_subject"] = "string"
    types_columns["identifier_sample"] = "string"
    types_columns["identifier_signal"] = "string"

    types_columns["sex_text"] = "string"
    types_columns["sex_female"] = "float32"
    types_columns["sex_y"] = "float32"
    types_columns["age_cohort_text"] = "string"
    types_columns["age_cohort"] = "float32"
    types_columns["visit_text"] = "string"
    types_columns["visit_second"] = "float32"
    types_columns["tissue"] = "string"
    types_columns["intervention_text"] = "string"
    types_columns["intervention_placebo"] = "float32"
    types_columns["intervention_omega3"] = "float32"
    types_columns["intervention_placebo_other"] = "float32"
    types_columns["intervention_omega3_other"] = "float32"
    types_columns["intervention_after_placebo"] = "float32"
    types_columns["intervention_after_omega3"] = "float32"
    types_columns["exercise_time_point"] = "string"

    types_columns["age"] = "float32"
    types_columns["body_mass_index"] = "float32"

    types_columns["triglyceride"] = "float32"
    types_columns["omega3_eicosapentaenoate"] = "float32"
    types_columns["omega3_docosahexaenoate"] = "float32"
    types_columns["cholesterol"] = "float32"
    types_columns["lipoprotein_hdl"] = "float32"
    types_columns["lipoprotein_nonhdl"] = "float32"
    types_columns["lipoprotein_ldl"] = "float32"

    types_columns["thyroid_stimulate_hormone"] = "float32"
    types_columns["insulin"] = "float32"
    types_columns["insulin_sensitivity"] = "float32"
    types_columns["homa_insulin_resist"] = "float32"

    types_columns["c_react_protein"] = "float32"
    types_columns["adipocyte_diameter"] = "float32"
    types_columns["adipocyte_lipid_content"] = "float32"
    types_columns["mitochondrial_respiration_maximum"] = "float32"
    types_columns["oxygen_consumption"] = "float32"

    types_columns["white_blood_cells"] = "float32"
    types_columns["neutrophils"] = "float32"
    types_columns["lymphocytes"] = "float32"
    types_columns["monocytes"] = "float32"
    types_columns["eosinophils"] = "float32"
    types_columns["basophils"] = "float32"
    types_columns["cd68_adipose_percent"] = "float32"
    types_columns["p16_adipose_percent"] = "float32"
    types_columns["cd14_adipose_percent"] = "float32"
    types_columns["cd206_adipose_percent"] = "float32"

    # Return information.
    return types_columns


def read_organize_source_subject_sample(
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
    #paths["in_demonstration"]
    #paths["in_parameters"]
    #paths["in_parameters_private"]

    # Define paths to child files.
    path_file_table_subject = os.path.join(
        paths["out_project"], "phenotypes", "organize_subject", "tables",
        "table_subject.tsv",
    )
    path_file_table_sample = os.path.join(
        paths["out_project"], "phenotypes", "organize_sample", "tables",
        "table_sample.tsv",
    )

    # Read information from file.

    # Table of properties and attributes for subjects and samples.
    types_columns = define_type_table_columns_subject_sample()
    table = pandas.read_csv(
        path_file_table_sample,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Table of parameters for organization of features.
    pail_feature = (
        aexph_sub.read_organize_source_table_subject_feature_organization(
            paths=paths,
            report=False,
    ))
    #pail_feature["types_columns"]
    #pail_feature["translations_feature_forward"]
    #pail_feature["translations_feature_reverse"]
    #pail_feature["columns_all"]
    #pail_feature["columns_quantitative"]
    #pail_feature["columns_olink_plasma"]
    #pail_feature["columns_olink_muscle"]
    #pail_feature["columns_olink_adipose"]

    # Organize information.
    # Determine names of columns for unique, relevant features.
    # The names from the keys of dictionary "types_columns" appear first in the
    # list, giving these priority in subsequent sorts of columns in the table.
    features_relevant = list(types_columns.keys())
    features_relevant.extend(pail_feature["columns_quantitative"])
    features_relevant.extend(pail_feature["columns_olink_plasma"])
    features_relevant.extend(pail_feature["columns_olink_muscle"])
    features_relevant.extend(pail_feature["columns_olink_adipose"])
    features_relevant = putly.collect_unique_elements(
        elements=features_relevant,
    )

    # Collect information.
    pail = dict()
    pail["table_subject_sample"] = table.copy(deep=True)
    pail["features_relevant"] = copy.deepcopy(features_relevant)

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
        print("function: read_organize_source_subject_sample()")
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

    # Read source information for parameters.
    pail_parameter = read_organize_source_parameter(
        paths=paths,
        report=report,
    )
    #pail["table"]
    #pail["records"]

    # Read source information for parameters.
    pail_sample = read_organize_source_subject_sample(
        paths=paths,
        report=False,
    )
    #pail["features_relevant"]
    #pail["table_subject_sample"]

    # Collect information.
    pail = dict()
    pail["table_sample"] = pail_sample["table_subject_sample"]
    pail["features_relevant"] = pail_sample["features_relevant"]
    pail["instances_parameter"] = pail_parameter["records"]

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        #print("table of measurements in batch a: ")
        #print(pail["table_batch_a"].iloc[0:10, 0:])
        #putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Organize information in table for features and observations.






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
    procedure="compare_groups"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
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
    #paths["out_procedure_tables"]
    #paths["out_procedure_plot"]

    ##########
    # 2. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=report,
    )
    #pail_source["instances_parameter"]
    #pail_source["table_sample"]

    table = pail_source["table_sample"].copy(deep=True)
    table["observations"] = table["subject_visit"]


    selection_observations = (
        pail_source["instances_parameter"][0]["selection_observations"]
    )

    table = porg.prepare_table_features_observations_for_analysis(
        table=table,
        selection_observations=selection_observations,
        features_relevant=pail_source["features_relevant"],
        features_essential=list(),
        features_continuity_scale=list(),
        index_columns_source="features",
        index_columns_product="features",
        index_rows_source="observations",
        index_rows_product="observations",
        remove_missing=False,
        remove_redundancy=True,
        adjust_scale=False,
        method_scale=None, # "z_score" or "unit_range"
        explicate_indices=True,
        report=False,
    )
    print(table)




    pass


###############################################################################
# End
