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

# TODO: TCW; 14 April 2025
# I'm still trying to figure out the introduction of multiple groups for an observations.
# I've been encountering difficulty due to the presence of both "adipose" and "muscle"
# records for subjects.
# Rather than using "table_sample", use the simpler, more original "table_subject".
# The tissue information only became available with RNA Seq on adipose or muscle samples.



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
            str(row["sequence"]).strip(),
            str(row["group"]).strip(),
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


def define_type_table_columns_subject_sample_quantitative_continuous():
    """
    Defines the types of variables for columns in table.

    Review: TCW; 14 April 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Define names of columns for features on quantitative, continuous scales
    # of measurement.
    names_features = [
        "age",
        "body_mass_index",

        "triglyceride",
        "omega3_eicosapentaenoate",
        "omega3_docosahexaenoate",
        "cholesterol",
        "lipoprotein_hdl",
        "lipoprotein_nonhdl",
        "lipoprotein_ldl",

        "thyroid_stimulate_hormone",
        "insulin",
        "insulin_sensitivity",
        "homa_insulin_resist",

        "c_react_protein",
        "adipocyte_diameter",
        "adipocyte_lipid_content",
        "mitochondrial_respiration_maximum",
        "oxygen_consumption",

        "white_blood_cells",
        "neutrophils",
        "lymphocytes",
        "monocytes",
        "eosinophils",
        "basophils",
        "cd68_adipose_percent",
        "p16_adipose_percent",
        "cd14_adipose_percent",
        "cd206_adipose_percent",
    ]

    # Specify types of variables in columns of table.
    types_columns = dict()
    for name in names_features:
        types_columns[name] = "float32"
        pass

    # Collect information.
    pail = dict()
    pail["names_features"] = names_features
    pail["types"] = types_columns
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

    #types_columns["subject_tissue_visit"] = "string"
    types_columns["subject_visit"] = "string"
    types_columns["identifier_subject"] = "string"
    #types_columns["identifier_sample"] = "string"
    #types_columns["identifier_signal"] = "string"

    types_columns["sex_text"] = "string"
    types_columns["sex_female"] = "float32"
    types_columns["sex_y"] = "float32"
    types_columns["age_cohort_text"] = "string"
    types_columns["age_cohort"] = "float32"
    types_columns["visit_text"] = "string"
    types_columns["visit_second"] = "float32"
    types_columns["intervention_text"] = "string"
    types_columns["intervention_placebo"] = "float32"
    types_columns["intervention_omega3"] = "float32"
    types_columns["intervention_placebo_other"] = "float32"
    types_columns["intervention_omega3_other"] = "float32"
    types_columns["intervention_after_placebo"] = "float32"
    types_columns["intervention_after_omega3"] = "float32"
    #types_columns["tissue"] = "string"
    #types_columns["exercise_time_point"] = "string"

    pail_quantitative = (
        define_type_table_columns_subject_sample_quantitative_continuous()
    )
    types_columns_quantitative = pail_quantitative["types"]
    types_columns.update(types_columns_quantitative)

    # Return information.
    return types_columns


# In this function it is convenient to designate additional features as
# relevant to analyses.
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
        path_file_table_subject,
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
    pail_quantitative = (
        define_type_table_columns_subject_sample_quantitative_continuous()
    )
    features_quantitative = pail_quantitative["names_features"]
    #features_relevant.extend(pail_feature["columns_quantitative"])
    #features_relevant.extend(pail_feature["columns_olink_plasma"])
    #features_relevant.extend(pail_feature["columns_olink_muscle"])
    #features_relevant.extend(pail_feature["columns_olink_adipose"])
    features_relevant = putly.collect_unique_elements(
        elements=features_relevant,
    )
    features_quantitative = putly.collect_unique_elements(
        elements=features_quantitative,
    )

    # Collect information.
    pail = dict()
    pail["table_subject_sample"] = table.copy(deep=True)
    pail["features_relevant"] = copy.deepcopy(features_relevant)
    pail["features_quantitative"] = copy.deepcopy(features_quantitative)
    pail["translations_feature_reverse"] = (
        pail_feature["translations_feature_reverse"]
    )

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
    pail["features_quantitative"] = pail_sample["features_quantitative"]
    pail["instances_parameter"] = pail_parameter["records"]
    pail["translations_feature_reverse"] = (
        pail_sample["translations_feature_reverse"]
    )

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
# 3. Select columns and rows in table for features and observations that are
#    relevant to each analysis.
#    Organize information in table for features and observations.


# Alternative definition of custom groups of observations.
def define_custom_selection_groups_observations(
    table=None,
    index_rows=None,
    report=None,
):
    """
    Defines groups of observations.

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    table_source = table.copy(deep=True)

    # Collect information.
    pail = dict()

    # Names of groups in sequence for sort order.
    pail["names_groups_observations_sequence"] = list()
    pail["names_groups_observations_sequence"].append("younger_older")
    pail["names_groups_observations_sequence"].append("younger")
    pail["names_groups_observations_sequence"].append("older")
    pail["names_groups_observations_sequence"].append("older_placebo_omega3")

    # Prepare reference to sort rows by group in a subsequent table.
    pail["sequence_groups_observations"] = dict(zip(
        pail["names_groups_observations_sequence"],
        range(len(pail["names_groups_observations_sequence"]))
    ))

    # Selections of observations in groups.
    pail["groups_observations"] = dict()
    pail["groups_observations"]["younger_older"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier=index_rows,
            name="younger_control", # or "name_instance"
            columns_categories={
                #"tissue": ["adipose",], # column "tissue" not in original subject table
                "age_cohort_text": ["younger",],
                "sex_text": ["female","male",],
                "visit_text": ["first",],
            },
            report=report,
    ))
    pail["groups_observations"]["older_control"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier=index_rows,
            name="older_control", # or "name_instance"
            columns_categories={
                #"tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "sex_text": ["female","male",],
                "visit_text": ["first",],
            },
            report=report,
    ))
    pail["groups_observations"]["older_placebo"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier=index_rows,
            name="older_placebo", # or "name_instance"
            columns_categories={
                #"tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "sex_text": ["female","male",],
                "visit_text": ["second",],
                "intervention_text": ["placebo",],
            },
            report=report,
    ))
    pail["groups_observations"]["older_omega3"] = (
        porg.filter_extract_table_row_identifiers_by_columns_categories(
            table=table_source,
            column_identifier=index_rows,
            name="older_omega3", # or "name_instance"
            columns_categories={
                #"tissue": ["adipose",],
                "age_cohort_text": ["elder",],
                "sex_text": ["female","male",],
                "visit_text": ["second",],
                "intervention_text": ["omega3",],
            },
            report=report,
    ))

    # Strategy to filter rows in table by collection of all that appear in the
    # various relevant groups.

    # Total selection of observations.
    pail["observations_selection"] = list()
    for group_observations in pail["groups_observations"].keys():
        pail["observations_selection"].extend(
            pail["groups_observations"][group_observations]
        )
        pass

    # Filter rows in table.
    #table_selection = table_source.loc[
    #    table_source[index_rows].isin(
    #        pail["observations_selection"]
    #    ), :
    #].copy(deep=True)

    # Return information.
    return pail


def collect_entries_tables_identifiers_groups_observations(
    instances=None,
    parameters_group=None,
    replicate_groups=None,
    table=None,
    index_columns=None,
    index_rows=None,
    features_relevant=None,
    features_essential=None,
    paths=None,
    report=None,
):
    """
    Collect clean tables for features across groups of observations, as defined
    by instances of parameters.

    This procedure selects from a source table the columns corresponding to
    relevant features and the rows corresponding to relevant observations
    within specific groups. The procedure also applies some general operations
    to prepare information in the table for analyses. The procedure organizes
    the product information both as partial, separate, stratified tables within
    a dictionary and as a whole table with a new column that designates groups
    of observations.

    For the collection of stratified tables within the dictionary
    'entries_tables', any individual observation can belong to multiple groups
    and can appear in multiple stratified tables accordingly.

    For the table 'table_group', it is optional to replicate rows for
    individual observations that occur in multiple groups. It is important to
    be aware and cautious when using this option of replication of records for
    observations. This option can be convenient for describing features in
    groups of observations that overlap.

    ----------
    Format of source data table (name: 'table')
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    identifiers     feature_1 feature_2 feature_3 feature_4 feature_5 ...

    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product data tables, partial (name: 'table')
    ----------
    Format of product, partial data tables organized within entries of a
    dictionary is in wide format with features across columns and values
    corresponding to their observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. The table has explicitly named indices across
    columns and rows.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observations
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product data table, whole (name: 'table_group')
    ----------
    Format of product, whole data table is in wide format with features across
    columns and values corresponding to their observations across rows. A
    special header row gives identifiers or names corresponding to each feature
    across columns, and a special column gives identifiers or names
    corresponding to each observation across rows. Another special column
    provides names of categorical groups for these observations. The table has
    explicitly named indices across columns and rows.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    Review: TCW; 14 April 2025

    arguments:
        instances (list<dict>): multiple instances, each with parameters for
            the preparation and analysis of a selection of data from a main
            table with information about features and their measurements across
            observations
            execution (int): logical binary indicator of whether to execute and
                handle the parameters for the current instance
            sequence (int): sequential index for instance's name and sort order
            group (str): categorical group of instances
            name (str): name or designator for instance of parameters
            name_combination (str): compound name for instance of parameters
            selection_observations (dict<list<str>>): names of columns in data
                table for feature variables and their categorical values by
                which to filter rows for observations in data table
            review (str):
            note (str):
        parameters_group (list<str>): parameters to combine as names of groups
            of instances from the parameters
        replicate_groups (bool): whether to replicate records or rows in table
            for individual observations that belong to multiple groups
        table (object): Pandas data-frame table of data with features
            and observations for analysis
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant, for which to keep columns and
            to remove rows of observations with redundancy across all
        features_essential (list<str>): names of columns in data table for
            feature variables that are essential, for which to remove rows of
            observations with any missing values
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
            groups_sequence (list<str>): unique names of groups of observations
                corresponding to the separate tables in the original sequence
                of instances of parameters
            entries_tables (dict<object>): entries with names as keys and as
                values, Pandas data-frame tables of data with features and
                observations for analysis
            table_group (object): Pandas data-frame table of data with features
                and observations for analysis

    """

    # Define subordinate functions for internal use.
    def alias_filter_extract(
        table=None,
        column_identifier=None,
        name=None,
        columns_categories=None,
        report=None,
    ):
        identifiers = (
            porg.filter_extract_table_row_identifiers_by_columns_categories(
                table=table,
                column_identifier=column_identifier,
                name=name,
                columns_categories=columns_categories,
                report=report,
        ))
        return identifiers

    # Copy information.
    table_source = table.copy(deep=True)
    instances = copy.deepcopy(instances)
    features_relevant_first = copy.deepcopy(features_relevant)

    # Organize information in table.
    # To avoid errors in the management of indices, the source name of index
    # across rows needs to be different than the product name.
    table_source["observations_preliminary"] = table_source[index_rows]
    features_relevant_first.insert(0, "observations_preliminary")

    # Organize information in table as a whole before interating on instances
    # of parameters for selection of features and observations in groups.
    table = porg.prepare_table_features_observations_for_analysis(
        table=table_source,
        selection_observations=dict(),
        features_relevant=features_relevant_first,
        features_essential=features_essential,
        features_continuity_scale=list(),
        index_columns_source=index_columns,
        index_columns_product="features",
        index_rows_source="observations_preliminary", # not same as product
        index_rows_product="observations",
        remove_missing=False,
        remove_redundancy=True,
        adjust_scale=False,
        method_scale=None, # "z_score" or "unit_range"
        explicate_indices=False,
        report=False,
    )

    # Copy information.
    features_relevant_second = copy.deepcopy(features_relevant)
    # Organize information.
    features_relevant_second.insert(0, "observations")
    # Collect information.
    pail = dict()
    pail["groups_sequence"] = list()
    pail["sequence_groups_observations"] = dict()
    pail["identifiers_rows"] = dict()
    pail["entries_tables"] = dict()
    # Iterate across instances of parameters.
    #instances = [instances[0]]
    for instance in instances:
        if (int(instance["execution"]) == 1):
            # Determine name of group.
            #name_combination = str(instance["name_combination"])
            #name_group = str(instance["group"] + "_-_" + instance["name"])
            items = list()
            for parameter in parameters_group:
                items.append(instance[parameter])
                pass
            name_group = "_-_".join(items)
            # Collect identifiers of rows in table.
            identifiers_rows = (
                alias_filter_extract(
                    table=table,
                    column_identifier="observations",
                    name=instance["name"],
                    columns_categories=instance["selection_observations"],
                    report=report,
            ))
            # Organize information in table.
            table_part = porg.prepare_table_features_observations_for_analysis(
                table=table,
                selection_observations=instance["selection_observations"],
                features_relevant=features_relevant_second,
                features_essential=features_essential,
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
            # Collect information.
            pail["groups_sequence"].append(name_group)
            pail["identifiers_rows"][name_group] = copy.deepcopy(
                identifiers_rows
            )
            pail["entries_tables"][name_group] = table_part
            pail["sequence_groups_observations"][name_group] = (
                int(instance["sequence"])
            )
            pass
        pass
    # Total selection of observations.
    observations_selection = list()
    for group_observations in pail["identifiers_rows"].keys():
        observations_selection.extend(
            pail["identifiers_rows"][group_observations]
        )
        pass
    # Filter rows in table.
    table_selection = table.loc[
        table["observations"].isin(
            observations_selection
        ), :
    ].copy(deep=True)

    # Determine and fill groups of observations.
    # This function creates a new column named "group" and assigns categorical
    # names corresponding to specific sets of observations.
    # Each observation only belongs to a single group.
    if (replicate_groups):
        table_group = porg.determine_fill_table_groups_rows_with_replicates(
            table=table_selection,
            index_rows="observations",
            column_group="group",
            groups_rows=pail["identifiers_rows"],
            report=False,
        )
    else:
        table_group = porg.determine_fill_table_groups_rows(
            table=table_selection,
            index_rows="observations",
            column_group="group",
            groups_rows=pail["identifiers_rows"],
            report=False,
        )
        pass
    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=index_rows,
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )
    # Separate information.
    table_group["group_group"] = table_group.apply(
        lambda row: str(row["group"]).split(sep="_-_")[0],
        axis="columns", # apply function to each row
    )
    table_group["group_name"] = table_group.apply(
        lambda row: str(row["group"]).split(sep="_-_")[1],
        axis="columns", # apply function to each row
    )

    # Filter and sort columns within table.
    columns_sequence = copy.deepcopy(features_relevant)
    columns_sequence.insert(0, "group")
    columns_sequence.insert(0, "group_name")
    columns_sequence.insert(0, "group_group")
    columns_sequence.insert(0, "observations")
    # Organize information in table as a whole before interating on instances
    # of parameters for selection of features and observations in groups.
    table_group = porg.prepare_table_features_observations_for_analysis(
        table=table_group,
        selection_observations=dict(),
        features_relevant=columns_sequence,
        features_essential=features_essential,
        features_continuity_scale=list(),
        index_columns_source="features",
        index_columns_product="features",
        index_rows_source=None, # not same as product
        index_rows_product="observations_novel",
        remove_missing=False,
        remove_redundancy=True,
        adjust_scale=False,
        method_scale=None, # "z_score" or "unit_range"
        explicate_indices=True,
        report=False,
    )

    # Collect information.
    pail["table_group"] = table_group

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: phenotypes")
        print("module: compare_groups.py")
        print("function: collect_entries_tables_groups_observations()")
        putly.print_terminal_partition(level=5)
        print(str(
            "whole table after basic preparation and introduction of groups:"
        ))
        print(table_group)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 4. Create plot charts for quantitative features in groups of observations.


def manage_plot_write_groups_observations_box(
    table=None,
    column_feature=None,
    column_group=None,
    translations_feature=None,
    path_directory_parent=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            corresponding to features across columns and observations across
            rows
        column_feature (str): name of column in original source table for a
            feature on a quantitative, continuous, interval, or ratio scale of
            measurement
        column_group (str): name of column in table to use for groups
        translations_feature (dict<str>): translations for names of features
        path_directory_parent ((str): path to parent directory within which to
            write files
        report (bool): whether to print reports

    raises:

    returns:


    """

    ##########
    # Organize information for plot.
    pail_extract = porg.extract_array_values_from_table_column_by_groups_rows(
        table=table,
        column_group=column_group,
        column_feature=column_feature,
        report=False,
    )
    # Determine title.
    if (
        (translations_feature is not None) and
        (column_feature in translations_feature.keys())
    ):
        title = translations_feature[column_feature]
    else:
        title = column_feature
        pass

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Extract parameters for colors.
    colors_names_groups = [
        "purple_violet",
        "blue_sky",
        "green_mint",
        "yellow_sunshine",
    ]
    if (
        (colors_names_groups is not None) and
        (len(colors_names_groups) > 0)
    ):
        colors_groups = list(map(
            lambda color_name: copy.deepcopy(colors[color_name]),
            colors_names_groups
        ))
    else:
        colors_groups = None
        pass

    # Create figure.
    figure = pplot.plot_boxes_groups(
        values_groups=pail_extract["values_nonmissing_groups"],
        names_groups=pail_extract["names_groups"],
        title_ordinate=title,
        title_abscissa="",
        title_chart_top_center="",
        colors_groups=colors_groups,
        label_top_center="",
        label_top_left="",
        label_top_right="",
        aspect="landscape",
        orientation_box="vertical",
        axis_linear_minimum=0.0,
        fonts=fonts,
        colors=colors,
        report=report,
    )

    ##########
    # Collect information.
    pail_write_plot_box = dict()
    pail_write_plot_box[column_feature] = figure_box
    pail_write_plot_vioin = dict()
    pail_write_plot_vioin[column_feature] = figure_violin

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_box = os.path.join(
        path_directory_parent, "box",
    )
    path_directory_violin = os.path.join(
        path_directory_parent, "violin",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_box,
    )
    putly.create_directories(
        path=path_directory_violin,
    )
    # Write figures to file.
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot_box,
        format="svg", # jpg, png, svg
        resolution=300,
        path_directory=path_directory_box,
    )
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot_violin,
        format="svg", # jpg, png, svg
        resolution=300,
        path_directory=path_directory_violin,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "manage_plot_write_groups_observations_box()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)

    # Return information.
    pass


##########
# 5. Describe quantitative features in groups of observations.


def describe_quantitative_features_by_observations_groups(
    entries_tables=None,
    groups_sequence=None,
    index_columns=None,
    index_rows=None,
    index_features=None,
    columns_features=None,
    translations_feature=None,
    table_group=None,
    column_group=None,
    ttest_one=None,
    ttest_two=None,
    report=None,
):
    """
    Describe features on quantitative, continuous, interval or ratio scale of
    measurement in terms of their values within groups of observations. For
    each feature, this function prepares a table of summary, descriptive,
    statistical measurements.

    ----------
    Format of source table
    ----------
    Format of source table is in wide format with features across columns and
    values corresponding to their observations across rows. A special header
    row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. Another special column provides names of
    categorical groups for these observations. The table has explicitly named
    indices across columns and rows.
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
    Format of product table
    ----------
    Format of product table is in partial long format with summary statistics
    and measures across columns and features across rows. A special column
    gives identifiers corresponding to each feature across rows. Another
    special column provides names of categorical groups of observations. For
    versatility, this table does not have explicity defined indices across
    columns or rows.
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
    An alternative format for a product table.
    ----------
    A disadvantage of this format is that values within columns would have
    different variable types.
    ----------
    group                group_1   group_2   group_3   group_4   ...
    measure
    mean                 0.001     0.001     0.001     0.001     ...
    standard_error       0.001     0.001     0.001     0.001     ...
    standard_deviation   0.001     0.001     0.001     0.001     ...
    95_confidence_low    0.001     0.001     0.001     0.001     ...
    95_confidence_high   0.001     0.001     0.001     0.001     ...
    minimum              0.001     0.001     0.001     0.001     ...
    maximum              0.001     0.001     0.001     0.001     ...
    median               0.001     0.001     0.001     0.001     ...
    count_total          30        30        30        30        ...
    count_valid          25        25        25        25        ...
    ----------

    Review: TCW; 14 April 2025

    arguments:
        entries_tables (dict<object>): entries with names as keys and as
            values, Pandas data-frame tables with features across columns and
            observations across rows
        groups_sequence (list<str>): unique names of groups of observations
            corresponding to the separate tables in the original sequence
            of instances of parameters
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        index_features (str): name for index corresponding to features across
            rows in the novel, product table
        columns_features (list<str>): names of columns in original source
            table for a selection of features on a quantitative, continuous,
            interval or ratio scale of measurement
        translations_feature (dict<str>): translations for names of features
        table_group (object): Pandas data-frame table of features across
            columns and observations across rows with values on a quantitative,
            continuous, interval or ratio scale of measurement
        column_group (str): name of column in original source table that
            designates groups of observations across rows
        ttest_one (dict): collection of parameters for T-test
            name (str): name for a column in the summary table that reports the
                p-value of the T-test
            groups (list<str>): names of groups of observations between which
                to perform a T-test, only use names of two groups
            equal_variances (bool): whether to assume that values in both groups
                have equal variances
            independent_groups (bool): whether the groups are independent
                ('True'), or whether the groups are paired, related, or
                otherwise dependent ('False')
            hypothesis_alternative (str): description of the alternative
                hypothesis, either a 'one-sided' or 'one-tailed' analysis in
                which only 'less' or 'greater' is relevant or a 'two-sided' or
                'two-tailed' analysis in which both lesser and greater are
                relevant
        ttest_two (dict): collection of parameters for T-test
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information.
    table_group = table_group.copy(deep=True)
    columns_features = copy.deepcopy(columns_features)
    translations_feature = copy.deepcopy(translations_feature)
    entries_tables = copy.deepcopy(entries_tables)
    groups_sequence = copy.deepcopy(groups_sequence)
    ttest_one = copy.deepcopy(ttest_one)
    ttest_two = copy.deepcopy(ttest_two)

    ##########
    # Describe features in groups of observations.
    table_description_priority = (
        pdesc.describe_features_from_columns_by_separate_tables_rows(
            entries_tables=entries_tables,
            groups_sequence=groups_sequence,
            index_columns=index_columns,
            index_rows=index_rows,
            index_features=index_features,
            columns_features=columns_features,
            translations_feature=translations_feature,
            key_group=column_group,
            threshold_observations=5,
            digits_round=3,
            ttest_one=ttest_one,
            ttest_two=ttest_two,
            report=report,
    ))
    table_description_check = (
        pdesc.describe_features_from_table_columns_by_groups_rows(
            table_group=table_group,
            index_columns=index_columns,
            index_rows=index_rows,
            index_features=index_features,
            column_group=column_group,
            groups_sequence=groups_sequence,
            columns_features=columns_features,
            translations_feature=translations_feature,
            key_group=column_group,
            threshold_observations=5,
            digits_round=5,
            ttest_one=ttest_one,
            ttest_two=ttest_two,
            report=report,
    ))

    # Collect information.
    pail = dict()
    pail["table_priority"] = table_description_priority
    pail["table_check"] = table_description_check

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "describe_quantitative_features_by_observations_groups()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("Table of descriptive statistics and T-tests")
        print(pail["table_priority"])
        print(pail["table_check"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


###############################################################################
# Procedure


##########
# Manage main procedure.


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
    #pail_source["table_sample"]
    #pail_source["features_relevant"]
    #pail_source["instances_parameter"]

    ##########
    # 3. Collect clean tables for features across groups of observations, as
    # defined by instances of parameters.
    pail_parts = collect_entries_tables_identifiers_groups_observations(
        instances=pail_source["instances_parameter"],
        parameters_group=["group", "name",],
        replicate_groups=True, # caution
        table=pail_source["table_sample"],
        index_columns="features",
        index_rows="subject_visit",
        features_relevant=pail_source["features_relevant"],
        features_essential=[
            "identifier_subject",
            "sex_text",
            "age_cohort_text",
            "visit_text",
            #"intervention_text",
            "age",
        ],
        paths=paths,
        report=report,
    )
    #pail_parts["table_group"]
    #pail_parts["entries_tables"]
    #pail_parts["identifiers_rows_super"]
    #pail_parts["identifiers_rows"]
    #pail_parts["sequence_groups_observations"]
    print(pail_parts["table_group"])
    #print(pail_parts["entries_tables"]["younger_older"])

    ##########
    # 4. Create plot charts for quantitative features in groups of observations.
    # Plot.
    if False:
        for column_quantitative in pail_source["features_quantitative"]:
            manage_plot_write_groups_observations_box_violin(
                table=pail_parts["table_group"],
                column_feature=column_quantitative,
                column_group="group",
                column_directory="group_group",
                translations_feature=pail["translations_feature_reverse"],
                path_directory_parent=paths["out_procedure_plot"],
                report=False,
            )
            pass
        pass

    ##########
    # 5. Describe quantitative features in groups of observations.
    pail_description = describe_quantitative_features_by_observations_groups(
        entries_tables=pail_parts["entries_tables"],
        groups_sequence=pail_parts["groups_sequence"],
        index_columns="features",
        index_rows="observations",
        index_features="features",
        columns_features=pail_source["features_quantitative"],
        translations_feature=pail_source["translations_feature_reverse"],
        table_group=pail_parts["table_group"],
        column_group="group",
        ttest_one={
            "name": "p_ttest_age",
            "groups": [
                "adipose_age_-_younger", "adipose_age_-_older",
            ],
            "equal_variances": True,
            "independent_groups": True,
            "hypothesis_alternative": "two-sided",
        }, # or None
        ttest_two={
            "name": "p_ttest_sex",
            "groups": [
                "adipose_age_-_older_female", "adipose_age_-_older_male",
            ],
            "equal_variances": True,
            "independent_groups": True,
            "hypothesis_alternative": "two-sided",
        }, # or None
        report=report,
    )

    ##########
    # Perform T-test, ANOVA, and regression analyses.

    # T-tests and regressions for features against age cohorts.
    #entries_tables["1_adipose_age_younger_older"]
    # T-tests and regressions for features against omega-3 intervention.
    #entries_tables["8_adipose_diet_older_placebo_omega3"]


    # Collect information.
    # Collections of files.
    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_group")] = pail_parts["table_group"]
    pail_write_tables[str("table_description_priority")] = (
        pail_description["table_priority"]
    )
    pail_write_tables[str("table_description_check")] = (
        pail_description["table_check"]
    )
    pail_write_objects = dict()
    #pail_write_objects[str("samples")]

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

    pass


###############################################################################
# End
