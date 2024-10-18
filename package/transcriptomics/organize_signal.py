"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'organize_signal' is part of the 'transcriptomics' package within
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

# TODO: TCW; 15 October 2024
# TODO: simplify this process of organizing the signal information
# TODO: It might help to separate the preparation of the overall signal table
# TODO: from the stratified signal tables for individual analyses.



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
#import partner.regression as preg
import partner.plot as pplot
import partner.parallelization as prall
import exercise.transcriptomics.organize_sample as extr_sample

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.
# There is a hierarchy in these functions to initialize directories to manage
# the hierarchical tree structure of sub-procedures.


def initialize_directories_trunk(
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
        paths["dock"], "in_data",
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
    paths["out_whole"] = os.path.join(
        paths["out_procedure"], "whole",
    )
    paths["out_whole_preparation"] = os.path.join(
        paths["out_whole"], "preparation",
    )
    paths["out_whole_description"] = os.path.join(
        paths["out_whole"], "description",
    )
    paths["out_whole_plot"] = os.path.join(
        paths["out_whole"], "plot",
    )
    paths["out_parts"] = os.path.join(
        paths["out_procedure"], "parts",
    )
    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"],
        paths["out_whole"],
        paths["out_whole_preparation"],
        paths["out_whole_description"],
        paths["out_whole_plot"],
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
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: initialize_directories_trunk()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


def initialize_directories_before_branch(
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

    ##########
    # Initialize directories for trunk procedure.
    paths = initialize_directories_trunk(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=False,
        report=report,
    )

    ##########
    # Define paths to directories.

    # Initialize directories in main branch.
    paths_initialization = [
        paths["out_parts"],
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
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: initialize_directories_before_branch()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


def initialize_directories_branch_tissue(
    project=None,
    routine=None,
    procedure=None,
    tissue=None,
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
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        restore (bool): whether to remove previous versions of data
        report (bool): whether to print reports

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    ##########
    # Initialize directories for trunk procedure.
    paths = initialize_directories_trunk(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=False,
        report=report,
    )

    ##########
    # Initialize directories for instance-specific parallel branch procedure.
    # Define paths to directories.
    paths["out_tissue"] = os.path.join(
        paths["out_parts"], str(tissue),
    )
    paths["out_summary"] = os.path.join(
        paths["out_tissue"], "summary",
    )
    paths["out_summary_instances"] = os.path.join(
        paths["out_summary"], "instances",
    )
    # Initialize directories in main branch.
    paths_initialization = [
        paths["out_tissue"],
        paths["out_summary"],
        paths["out_summary_instances"],
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
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: initialize_directories_branch_tissue()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


def initialize_directories_branch_instance(
    project=None,
    routine=None,
    procedure=None,
    tissue=None,
    name_instance=None,
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
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        restore (bool): whether to remove previous versions of data
        report (bool): whether to print reports

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    ##########
    # Initialize directories for trunk procedure.
    paths = initialize_directories_branch_tissue(
        project=project,
        routine=routine,
        procedure=procedure,
        tissue=tissue,
        path_directory_dock=path_directory_dock,
        restore=False,
        report=report,
    )

    ##########
    # Initialize directories for instance-specific parallel branch procedure.
    # Define paths to directories.
    paths["out_instance"] = os.path.join(
        paths["out_tissue"], str(name_instance),
    )
    #paths["out_test"] = os.path.join(
    #    paths["out_instance"], "test",
    #)
    paths["out_data"] = os.path.join(
        paths["out_instance"], "data",
    )
    #paths["out_plot"] = os.path.join(
    #    paths["out_instance"], "plot",
    #)
    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_tissue"], # omit to avoid conflict in parallel branches
        paths["out_instance"],
        paths["out_data"],
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
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: initialize_directories_branch()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


##########
# 2. Read source information from file.


def define_column_types_table_parameter_instances():
    """
    Defines the variable types of columns within table for attributes of
    samples.

    Review: TCW; 14 August 2024

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion"] = "string" # "int32"
    types_columns["tissue"] = "string"
    types_columns["sort"] = "int32" # "int32"
    types_columns["group"] = "string"
    types_columns["instance"] = "string"
    types_columns["cohort_selection_primary"] = "string"
    types_columns["factor_availability"] = "string"
    types_columns["cohort_selection_secondary"] = "string"
    types_columns["continuity_scale"] = "string"
    types_columns["formula_text"] = "string"
    types_columns["condition"] = "string"
    types_columns["levels_condition"] = "string"
    types_columns["supplement_1"] = "string"
    types_columns["levels_supplement_1"] = "string"
    types_columns["supplement_2"] = "string"
    types_columns["levels_supplement_2"] = "string"
    types_columns["supplement_3"] = "string"
    types_columns["levels_supplement_3"] = "string"
    types_columns["subject"] = "string"
    types_columns["threshold_significance"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns

# TODO: TCW; 16 October 2024
# use putly.parse_extract_text_keys_values_semicolon_colon_comma()

def read_organize_source_parameter_instances(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): collections of source information about instances for
           selection of information that is specific to separate analyses

    """

    # Define paths to parent directories.
    #paths["in_data"]
    #paths["in_parameters"]
    #paths["in_parameters_private"]

    # Define paths to child files.
    path_file_table_parameter = os.path.join(
        paths["in_parameters_private"],
        "table_set_differential_expression.tsv",
    )

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_column_types_table_parameter_instances()
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

    # TODO: TCW; 24 September 2024
    # write a generic function for parsing information from each of the
    # columns:
    # cohort_selection_primary
    # factor_availability
    # cohort_selection_secondary

    # Collect information.
    instances = list()
    for index, row in table.iterrows():
        if (int(row["inclusion"]) == 1):
            # Collect parameters for current instance.
            pail = dict()
            # Collect names of unique columns relevant to current instance.
            columns = list()
            pail["sort"] = str(row["sort"])
            pail["group"] = str(row["group"])
            pail["tissue"] = str(row["tissue"])
            pail["name_instance"] = "_".join([
                str(row["tissue"]),
                str(row["sort"]),
                str(row["instance"])
            ])
            if (str(row["cohort_selection_primary"]).strip() != "none"):
                pail["cohort_selection_primary"] = dict()
                for part in row["cohort_selection_primary"].strip().split(";"):
                    part_split = part.split(":")
                    pail["cohort_selection_primary"][part_split[0]] = (
                        part_split[1].split(",")
                    )
                    pass
                columns.extend(list(pail["cohort_selection_primary"].keys()))
                pass
            else:
                pail["cohort_selection_primary"] = None
                pass
            if (str(row["factor_availability"]).strip() != "none"):
                pail["factor_availability"] = dict()
                for part in row["factor_availability"].strip().split(";"):
                    part_split = part.split(":")
                    pail["factor_availability"][part_split[0]] = (
                        part_split[1].split(",")
                    )
                    pass
                columns.extend(list(pail["factor_availability"].keys()))
                pass
            else:
                pail["factor_availability"] = None
                pass
            if (str(row["cohort_selection_secondary"]).strip() != "none"):
                pail["cohort_selection_secondary"] = dict()
                for part in row["cohort_selection_secondary"].strip(
                ).split(";"):
                    part_split = part.split(":")
                    pail["cohort_selection_secondary"][part_split[0]] = (
                        part_split[1].split(",")
                    )
                    pass
                columns.extend(list(pail["cohort_selection_secondary"].keys()))
                # Extract names of columns corresponding to feature variables for
                # which to calculate tertiles.
                columns_tertile = extract_source_columns_for_tertiles(
                    cohort_selection=pail["cohort_selection_secondary"],
                    report=report,
                )
                columns.extend(columns_tertile)
                pass
            else:
                pail["cohort_selection_secondary"] = None
                pass
            if (str(row["continuity_scale"]) != "none"):
                pail["continuity_scale"] = (
                    row["continuity_scale"].strip().split(",")
                )
                columns.extend(pail["continuity_scale"])
            else:
                pail["continuity_scale"] = None
                pass
            # Collect names of unique columns relevant to current instance.
            #columns.append("identifier_signal")
            #columns.append("inclusion")
            #columns.append("tissue")
            # Extract and include names of other columns.
            columns_formula = row["formula_text"].strip().split(",")
            # Remove any interaction terms, since these are not columns.
            columns_formula = list(filter(
                lambda column: (":" not in str(column)),
                columns_formula
            ))
            columns.extend(columns_formula)
            columns = putly.collect_unique_elements(
                elements=columns,
            )
            columns.remove("inclusion")
            pail["columns_set"] = columns
            instances.append(pail)
            pass
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: read_organize_source_parameter_instances()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        print("instances:")
        print(instances)
        print("instance[0]:")
        print(instances[0])
        pass
    # Return information.
    return instances


def define_column_types_table_main():
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
    types_columns["identifier_gene"] = "string"
    types_columns["gene_id"] = "string"
    types_columns["gene_name"] = "string"
    types_columns["exon_number"] = "string"
    types_columns["gene_type"] = "string"
    types_columns["chromosome"] = "string"
    # Return information.
    return types_columns


def read_source_sample(
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
    path_file_table_sample_file = os.path.join(
        paths["in_data"], "study_exercise_age", "subject_sample",
        "table_sample_file_rnaseq.tsv",
    )
    path_file_table_sample = os.path.join(
        paths["out_routine"], "organize_sample", "data",
        "table_sample.pickle",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of matches between samples and files.
    types_columns = extr_sample.define_type_columns_table_sample_file()
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
    # Table of samples and their attributes.
    pail["table_sample"] = pandas.read_pickle(
        path_file_table_sample,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: read_source_sample()")
        putly.print_terminal_partition(level=5)
        count_rows = (pail["table_sample"].shape[0])
        count_columns = (pail["table_sample"].shape[1])
        print("sample table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_sample"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_main(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
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
    if (tissue == "adipose"):
        path_file_table_main = os.path.join(
            paths["in_data"], "study_exercise_age", "transcriptomics",
            "quantification_2024-07-14",
            "organization", "quantification_rna_reads_gene_adipose.tsv",
        )
    elif (tissue == "muscle"):
        path_file_table_main = os.path.join(
            paths["in_data"], "study_exercise_age", "transcriptomics",
            "quantification_2024-07-14",
            "organization", "quantification_rna_reads_gene_muscle.tsv",
        )
        pass

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of values of intensity across samples and proteins.
    types_columns = define_column_types_table_main()
    pail["table_main"] = pandas.read_csv(
        path_file_table_main,
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
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: read_source()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        count_rows = (pail["table_main"].shape[0])
        count_columns = (pail["table_main"].shape[1])
        print("main table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_main"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_signal_for_description(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of source information read from file

    """

    # Define paths to parent directories.
    # paths["out_whole_preparation"]

    # Define paths to child files.
    path_file_table_gene = os.path.join(
        paths["out_whole_preparation"], str(
            "table_gene_" + str(tissue) + ".pickle"
        ),
    )
    path_file_table_signal = os.path.join(
        paths["out_whole_preparation"], str(
            "table_signal_" + str(tissue) + ".pickle"
        ),
    )
    path_file_table_signal_scale = os.path.join(
        paths["out_whole_preparation"], str(
            "table_signal_scale_" + str(tissue) + ".pickle"
        ),
    )
    # Collect information.
    pail = dict()
    # Read information from file.
    pail["table_gene"] = pandas.read_pickle(
        path_file_table_gene,
    )
    pail["table_signal"] = pandas.read_pickle(
        path_file_table_signal,
    )
    pail["table_signal_scale"] = pandas.read_pickle(
        path_file_table_signal_scale,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: read_source_signal_for_description()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        count_rows = (pail["table_signal_scale"].shape[0])
        count_columns = (pail["table_signal_scale"].shape[1])
        print("signal scale table: ")
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(pail["table_signal_scale"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_organize_write_summary_instances_tissue(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): collections of source information about instances for
           selection of information that is specific to separate analyses

    """

    # Extract and filter complete paths to child files within parent directory.
    paths_file = putly.extract_filter_child_file_names_paths(
        path_directory=paths["out_summary_instances"],
        name_file_prefix="record_",
        name_file_suffix=".pickle",
        name_file_not="_blarg_blarg_blarg_",
        report=report,
    )
    # Read information from files.
    # Iterate on paths to files.
    # Collect information from read of source files.
    records = list()
    for path_file in paths_file:
        record = putly.read_object_from_file_pickle(
            path_file=path_file,
        )
        records.append(record)
        pass
    # Create pandas data-frame table.
    table = pandas.DataFrame(data=records)
    # Sort rows within table.
    #table["sort"] = pandas.to_numeric(
    #    table["sort"],
    #    errors="coerce", # force any parse error values to missing "NaN"
    #    downcast="float", # cast type to smallest float type
    #)
    table["sort"] = table["sort"].astype("int32")
    table.sort_values(
        by=[
            "sort",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: read_organize_write_summary_instances()")
        putly.print_terminal_partition(level=5)
        print("table summary of counts of samples in each set:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Write product information to file.
    putly.write_table_to_file(
        table=table,
        name_file="table_counts_sets_sample",
        path_directory=paths["out_summary"],
        reset_index=False,
        write_index=False,
        type="text",
    )
    # Return information.
    pass


##########
# 3. Select set of samples for specific analyses.

# TODO: TCW; 16 October 2024
# use porg.filter_extract_table_row_identifiers_by_columns_categories()
def select_sets_identifier_table_sample(
    table_sample=None,
    name_instance=None,
    tissue=None,
    cohort_selection=None,
    factor_availability=None,
    report=None,
):
    """
    Selects sets of samples relevant to specific analyses.

    This function preserves the original sequence of samples from the source
    table in its extraction of identifiers. It is important to preserve the
    definitive sequence of samples from the table of their attributes as this
    sequence will determine the sort sequence of values in the table of
    signals, which must correspond exactly to the table of samples.

    arguments:
        table_sample (object): Pandas data-frame table of information about
            samples that correspond to signals within accompanying main table
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        cohort_selection (dict<list<str>>): filters on rows in table for
            selection of samples relevant to cohort for analysis
        factor_availability (dict<list<str>>): features and their values
            corresponding to factors of interest in the analysis
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)
    # Copy other information.
    cohort_selection = copy.deepcopy(cohort_selection)
    factor_availability = copy.deepcopy(factor_availability)

    # Organize indices in table.
    #table_sample.reset_index(
    #    level=None,
    #    inplace=True,
    #    drop=False, # remove index; do not move to regular columns
    #)

    # Translate names of columns.
    #translations = dict()
    #translations["identifier_signal"] = "identifier"
    #table_sample.rename(
    #    columns=translations,
    #    inplace=True,
    #)

    # Filter rows in table for selection of relevant samples.
    # Filter by inclusion indicator.
    table_sample["inclusion"] = table_sample["inclusion"].astype("str")
    table_inclusion = table_sample.loc[
        (table_sample["inclusion"] == "1"), :
    ].copy(deep=True)
    # Separate information about sets of samples for difference experimental
    # conditions or groups.
    table_tissue = table_inclusion.loc[
        (table_inclusion["tissue"] == tissue), :
    ].copy(deep=True)
    #table_factor = table_tissue.loc[
    #    (
    #        (table_tissue["condition"] == "control") |
    #        (table_tissue["condition"] == "intervention_1")
    #    ), :
    #].copy(deep=True)

    # Filter rows in table by rules for selection of a specific set.
    # Iterate on features and values for selection of samples in cohort.
    table_cohort = table_inclusion.copy(deep=True)
    if (cohort_selection is not None):
        for feature in cohort_selection.keys():
            table_cohort = table_cohort.loc[(
                table_cohort[feature].isin(cohort_selection[feature])
            ), :].copy(deep=True)
            pass
        pass

    # Iterate on factors and values for selection of samples on the basis
    # of availability.
    table_factor = table_cohort.copy(deep=True)
    if (factor_availability is not None):
        for factor in factor_availability.keys():
            table_factor = table_factor.loc[(
                table_factor[factor].isin(factor_availability[factor])
            ), :].copy(deep=True)
            pass

    # Copy information in table.
    table_selection = table_factor.copy(deep=True)

    # Separate samples for unique values of factor.
    # This operation would necessarily assume that there were two and only two
    # unique values of a single factor.
    # table_control =
    # table_case =

    # Extract identifiers of samples in separate groups.
    samples_inclusion = copy.deepcopy(
        table_inclusion["identifier_signal"].to_list()
    )
    samples_tissue = copy.deepcopy(
        table_tissue["identifier_signal"].to_list()
    )
    samples_selection = copy.deepcopy(
        table_selection["identifier_signal"].to_list()
    )

    # Organize indices in tables.
    if False:
        tables = [
            table_inclusion,
            table_tissue,
            table_selection,
        ]
        for table in tables:
            table.reset_index(
                level=None,
                inplace=True,
                drop=True, # remove index; do not move to regular columns
            )
            table.set_index(
                ["identifier_signal"],
                append=False,
                drop=True,
                inplace=True,
            )
            pass
        pass

    # Collect information.
    pail = dict()
    pail["table_inclusion"] = table_inclusion
    pail["table_tissue"] = table_tissue
    pail["table_selection"] = table_selection
    pail["samples_inclusion"] = samples_inclusion
    pail["samples_tissue"] = samples_tissue
    pail["samples_selection"] = samples_selection

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: select_sets_identifier_table_sample_primary()")
        print("name_instance: " + str(name_instance))
        print("tissue: " + tissue)
        #print("features for cohort selection:")
        #print(str(cohort_selection.keys()))
        #print("factors for availability of specific values:")
        #print(str(factor_availability.keys()))
        putly.print_terminal_partition(level=5)
        print("sample table, filtered by inclusion and tissue:")
        print(table_tissue.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("sample table, filtered by set selection rules:")
        print(table_selection.iloc[0:10, 0:])
        #putly.print_terminal_partition(level=5)
        #print("description of first categorical factor:")
        #print("first factor: " + str(factor_availability.keys()[0]))
        #table_selection[list(factor_availability.keys())[0]].describe(
        #    include=["category",]
        #)
        putly.print_terminal_partition(level=4)
        #print(
        #    "counts of samples with each unique categorical value of each " +
        #    "factor:"
        #)
        #for factor in factor_availability.keys():
        #    print("factor: " + factor)
        #    print(
        #        table_selection[factor].value_counts(
        #            dropna=False,
        #        )
        #    )
        #    putly.print_terminal_partition(level=4)
        #    pass
        #putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def extract_source_columns_for_tertiles(
    cohort_selection=None,
    report=None,
):
    """
    Extract from the parameters the names of columns for which to determine
    tertiles.

    arguments:
        cohort_selection (dict<list<str>>): filters on rows in table for
            selection of samples relevant to cohort for analysis
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): names of columns corresponding to feature variables for
            which to determine tertiles

    """

    # Copy other information.
    cohort_selection = copy.deepcopy(cohort_selection)
    # Determine whether to calculate tertiles for any feature variables with
    # continuous values.
    columns_tertile = list()
    if (
        (cohort_selection is not None) and
        any("tertiles_" in item for item in list(cohort_selection.keys()))
    ):
        # Extract names of columns corresponding to feature variables for which
        # to calculate tertiles.
        features_tertile = list()
        for feature in cohort_selection.keys():
            if ("tertiles_" in str(feature)):
                features_tertile.append(feature)
        # Extract names of columns.
        columns_tertile = list(map(
            lambda feature: feature.replace("tertiles_", ""),
            features_tertile
        ))
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: extract_source_columns_for_tertiles()")
        putly.print_terminal_partition(level=5)
        print("columns of feature variables for tertiles:")
        print(columns_tertile)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return columns_tertile


def organize_describe_summarize_table_sample_tertiles(
    table=None,
    cohort_selection=None,
    group=None,
    name_instance=None,
    tissue=None,
    paths=None,
    report=None,
):
    """
    Organize the determination of tertiles on for multiple feature variables
    with values on a continuous interval or ratio measurement scale.

    arguments:
        table (object): Pandas data-frame table of subjects, samples, and their
            attribute features
        cohort_selection (dict<list<str>>): filters on rows in table for
            selection of samples relevant to cohort for analysis
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    cohort_selection = copy.deepcopy(cohort_selection)

    # Determine whether to calculate tertiles for any feature variables with
    # continuous values.
    if (
        (cohort_selection is not None) and
        any("tertiles_" in item for item in list(cohort_selection.keys()))
    ):
        # Extract names of columns corresponding to feature variables for which
        # to calculate tertiles.
        columns_tertile = extract_source_columns_for_tertiles(
            cohort_selection=cohort_selection,
            report=report,
        )
        # Determine tertiles for stratification of sample cohorts.
        for column_source in columns_tertile:
            column_product = str("tertiles_" + column_source)
            pail_tertile = pdesc.determine_describe_quantiles_ordinal(
                table=table,
                column_source=column_source,
                column_product=column_product,
                count=3,
                text_string=True,
                name_prefix="tertile_",
                report=True,
            )
            pdesc.describe_quantiles_ordinal(
                table=pail_tertile["table"],
                column_source=column_source,
                column_product=column_product,
                columns_category=["sex_text",],
                report=True,
            )
            path_file = os.path.join(
                paths["out_summary"], str(column_product + ".txt"),
            )
            if not os.path.isfile(path_file):
                bins_text = list(map(
                    lambda value: str(round(value, 3)),
                    pail_tertile["bins"].tolist()
                ))
                putly.write_list_to_file_text(
                    elements=bins_text,
                    delimiter=", ",
                    name_file=column_product,
                    path_directory=paths["out_summary"],
                )
            pass
    else:
        pail_tertile = dict()
        pail_tertile["table"] = table
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: organize_table_sample_tertiles()")
        putly.print_terminal_partition(level=5)
        print("table of attributes for samples: ")
        print(pail_tertile["table"].iloc[0:10, 0:])
        print(pail_tertile["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_tertile["table"]


def select_sets_final_identifier_table_sample(
    table_sample=None,
    name_instance=None,
    tissue=None,
    continuity_scale=None,
    columns_set=None,
    report=None,
):
    """
    Selects sets of samples relevant to specific analyses.

    This function preserves the original sequence of samples from the source
    table in its extraction of identifiers. It is important to preserve the
    definitive sequence of samples from the table of their attributes as this
    sequence will determine the sort sequence of values in the table of
    signals, which must correspond exactly to the table of samples.

    arguments:
        table_sample (object): Pandas data-frame table of information about
            samples that correspond to signals within accompanying main table
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        continuity_scale (list<str>): names of columns for covariates with
            values on continuous scale of measurement, interval or ratio,
            for which to standardize the scale by z score
        columns_set (list<str>): names of columns for feature variables that
            are relevant to the current set or instance of parameters
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_selection = table_sample.copy(deep=True)

    # Copy other information.
    continuity_scale = copy.deepcopy(continuity_scale)
    columns_set = copy.deepcopy(columns_set)

    # Adjust the scale of specific variables by transformation to z score.
    if (continuity_scale is not None):
        # Calculate standard z scores.
        table_selection = (
            pscl.drive_transform_variables_distribution_scale_z_score(
                table=table_selection,
                columns=continuity_scale,
                report=report,
        ))
        pass

    # Filter rows in table for non-missing values across relevant columns.
    table_selection.dropna(
        axis="index",
        how="any",
        subset=columns_set,
        inplace=True,
    )

    # Filter and sort columns within table.
    # Here it is necessary to filter to unique names of columns because
    # otherwise the column "tissue" has a redundant occurrence due to the
    # extraction from the parameters for selection of samples in the cohort.
    columns_sequence = copy.deepcopy(columns_set)
    columns_sequence.insert(0, "tissue")
    columns_sequence.insert(0, "inclusion")
    columns_sequence.insert(0, "identifier_signal")
    columns_sequence = putly.collect_unique_elements(
        elements=columns_sequence,
    )
    if True:
        table_selection = porg.filter_sort_table_columns(
            table=table_selection,
            columns_sequence=columns_sequence,
            report=report,
        )

    # Organize indices in table.
    table_selection.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_selection.set_index(
        ["identifier_signal"],
        append=False,
        drop=True,
        inplace=True,
    )

    # Extract identifiers of samples in separate groups.
    #samples_selection = copy.deepcopy(
    #    table_selection["identifier_signal"].to_list()
    #)
    samples_selection = copy.deepcopy(
        table_selection.index.get_level_values(
            "identifier_signal"
        ).unique().to_list()
    )

    # Collect information.
    pail = dict()
    pail["table_selection"] = table_selection
    pail["samples_selection"] = samples_selection
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: select_sets_final_identifier_table_sample()")
        print("name_instance: " + str(name_instance))
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("sample table, filtered by set selection rules:")
        print(table_selection.iloc[0:10, 0:])
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail


def report_write_count_samples(
    samples=None,
    sort=None,
    group=None,
    name_instance=None,
    tissue=None,
    paths=None,
    report=None,
):
    """
    Writes to file information for a report summary.

    arguments:
        samples (list<str>): identifiers of samples corresponding to names of
            columns for measurement values of signal intensity across features
        sort (int): sequential index for sort order
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:


    """

    # Copy other information.
    samples = copy.deepcopy(samples)

    # Collect information.
    record = dict()
    record["sort"] = sort
    record["group"] = group
    record["name_instance"] = name_instance
    record["tissue"] = tissue
    record["count_samples"] = int(len(samples))
    # Write product information to file.
    putly.write_object_to_file_pickle(
        object=record,
        name_file=str("record_" + name_instance),
        path_directory=paths["out_summary_instances"],
    )
    # Return information.
    pass


##########
# 4. Organize information from source.


def define_column_sequence_table_main_gene():
    """
    Defines the columns in sequence within table.

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "identifier_gene",
        "gene_identifier",
        "gene_identifier_base",
        "gene_name",
        #"gene_exon_number",
        "gene_type",
        "gene_chromosome",
        #"Unnamed: 160",
    ]
    # Return information.
    return columns_sequence


def organize_table_main(
    table_main=None,
    columns_gene=None,
    samples=None,
    tissue=None,
    report=None,
):
    """
    Organizes information in tables about samples and measurement signals.

    arguments:
        table_main (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows, with a few additional columns for attributes
            of gene features
        columns_gene (list<str>): names of columns corresponding to
            information about genes
        samples (list<str>): identifiers of samples corresponding to names of
            columns for measurement values of signal intensity across features
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_main = table_main.copy(deep=True)
    # Copy other information.
    columns_gene = copy.deepcopy(columns_gene)
    samples = copy.deepcopy(samples)

    # Translate names of columns.
    translations = dict()
    translations["gene_id"] = "gene_identifier"
    translations["exon_number"] = "gene_exon_number"
    translations["chromosome"] = "gene_chromosome"
    table_main.rename(
        columns=translations,
        inplace=True,
    )
    # Determine base Ensembl identifier for genes.
    # reference: https://useast.ensembl.org/Help/Faq?id=488
    table_main["gene_identifier_base"] = table_main.apply(
        lambda row: str(row["gene_identifier"]).strip().split(".")[0],
        axis="columns", # apply function to each row
    )
    table_main["identifier_gene"] = table_main["gene_identifier_base"]
    # Replace values of zero for signal intensity with missing values.
    # Only replace values within table's columns for samples.
    # This implementation is more concise than iteration across specific
    # columns.
    # This operation is inaccurate for quantification of RNA sequence read
    # data. In this case, it might even be more reasonable to replace missing
    # values with values of zero.
    #table_main[samples] = table_main[samples].replace(
    #    to_replace=0,
    #    value=pandas.NA,
    #)
    # Replace values less than zero with missing values.
    table_main[samples][table_main[samples] < 0] = pandas.NA

    # Filter and sort columns within table.
    columns_sequence = copy.deepcopy(columns_gene)
    columns_sequence.extend(samples)
    table_main = porg.filter_sort_table_columns(
        table=table_main,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Collect information.
    pail = dict()
    pail["table_main"] = table_main
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: organize_table_main()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("main table: ")
        print(table_main.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        print("description of categorical gene type:")
        print(table_main["gene_type"].describe(include=["category",]))
        putly.print_terminal_partition(level=5)
        print(
            "counts of genes with each unique categorical value of "
            + "gene type:")
        print(table_main["gene_type"].value_counts(dropna=False))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 5. Filter columns and rows in main table.


def define_keep_types_gene_broad():
    """
    Defines the categorical types of genes for which to keep rows in the main
    table.

    Reference:
    https://www.gencodegenes.org/pages/biotypes.html

    arguments:

    raises:

    returns:
        (list<str>): types of gene to keep

    """

    # Specify categorical types.
    types_gene = [
        "protein_coding",
        "lncRNA",
        #"processed_pseudogene",
        #"unprocessed_pseudogene",
        "misc_RNA",
        "snRNA",
        "miRNA",
        "TEC",
        #"transcribed_unprocessed_pseudogene"
        "snoRNA",
        #"transcribed_processed_pseudogene",
        #"rRNA_pseudogene",
        #"IG_V_pseudogene",
        #"transcribed_unitary_pseudogene",
        "IG_V_gene",
        "TR_V_gene",
        #"unitary_pseudogene",
        "TR_J_gene",
        "rRNA",
        "scaRNA",
        "IG_D_gene",
        #"TR_V_pseudogene",
        "Mt_tRNA",
        #"artifact",
        "IG_J_gene",
        "IG_C_gene",
        #"IG_C_pseudogene",
        "ribozyme",
        "TR_C_gene",
        "sRNA",
        "TR_D_gene",
        #"pseudogene",
        "vault_RNA",
        #"TR_J_pseudogene",
        #"IG_J_pseudogene",
        "Mt_rRNA",
        #"translated_processed_pseudogene",
        "scRNA",
        #"IG_pseudogene",
    ]
    # Return information.
    return types_gene


def define_keep_types_gene_narrow():
    """
    Defines the categorical types of genes for which to keep rows in the main
    table.

    Reference:
    https://www.gencodegenes.org/pages/biotypes.html

    arguments:

    raises:

    returns:
        (list<str>): types of gene to keep

    """

    # Specify categorical types.
    types_gene = [
        "protein_coding",
        #"IG_C_gene",
        #"IG_D_gene",
        #"IG_J_gene",
        #"IG_LV_gene",
        #"IG_V_gene",
        #"TR_C_gene",
        #"TR_J_gene",
        #"TR_V_gene",
        #"TR_D_gene",
    ]
    # Return information.
    return types_gene


def define_keep_gene_chromosomes(
    remove_sex_chromosomes=None,
):
    """
    Defines the categorical types of genes for which to keep rows in the main
    table.

    Reference:
    https://www.gencodegenes.org/pages/biotypes.html

    arguments:
        remove_sex_chromosomes (bool): whether to remove all genes on sex
            chromosomes

    raises:

    returns:
        (list<str>): identifiers of chromosomes for which to keep genes

    """

    # Specify categorical types.
    if remove_sex_chromosomes:
        chromosomes = [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            #"chrX",
            #"chrY",
            "chrM",
        ]
    else:
        chromosomes = [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
            "chrM",
        ]
    # Return information.
    return chromosomes


def determine_keep_series_by_identity(
    row_identifier=None,
    row_type=None,
    row_chromosome=None,
    identifier_prefix=None,
    types_gene=None,
    chromosomes=None,
):
    """
    Determines whether to keep a row from a table.

    arguments:
        row_identifier (str): current row's identifier of gene
        row_type (str): current row's type of gene
        row_chromosome (str): identifier of chromosome of gene
        identifier_prefix (str): prefix in identifier of gene to keep
        types_gene (list<str>): types of gene to keep
        chromosomes (list<str>): identifiers of chromosomes for which to keep
            genes

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Determine whether to keep current row from table.
    #(any(row_type == type_gene for type_gene in types_gene))
    if (
        (pandas.notna(row_identifier)) and
        (len(str(row_identifier)) > 0) and
        (str(identifier_prefix) in str(row_identifier)) and
        (str(row_type) in types_gene) and
        (str(row_chromosome) in chromosomes)
    ):
        indicator = 1
    else:
        indicator = 0
        pass
    # Return information.
    return indicator


def filter_table_main_rows_signal(
    table_main=None,
    samples_all=None,
    samples_control=None,
    samples_intervention=None,
    filter_rows_signal_by_condition=None,
    threshold_signal_low=None,
    threshold_signal_high=None,
    proportion_signal_all=None,
    proportion_signal_control=None,
    proportion_signal_intervention=None,
    tissue=None,
    report=None,
):
    """
    Filters information in table for gene features across rows.

    This function filters rows for gene features by different thresholds on the
    proportion of a gene's signal intensities that must have non-missing,
    within-threshold, valid values across sets of samples in either the control
    or intervention experimental conditions or both. The reason for using these
    different thresholds is to allow genes to have signals that are detectable
    in only one experimental condition or the other. An example would be an
    experimental intervention that activates or deactivates transcriptional
    expression of specific genes.

    arguments:
        table_main (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows, with a few additional columns for attributes
            of gene features
        samples_all (list<str>): identifiers of samples in both control and
            intervention experimental conditions corresponding to names of
            columns for measurement values of signal intensity across features
        samples_control (list<str>): identifiers of samples in control
            experimental condition corresponding to names of columns for
            measurement values of signal intensity across features
        samples_intervention (list<str>): identifiers of samples in
            intervention experimental condition corresponding to names of
            columns for measurement values of signal intensity across features
        filter_rows_signal_by_condition (bool): whether to filter rows by
            signal with separate consideration of columns for different
            experimental conditions
        threshold_signal_low (float): threshold below which
            (value <= threshold) all values are considered invalid and missing
        threshold_signal_high (float): threshold above which
            (value > threshold) all values are considered invalid and missing
        proportion_signal_all (float): proportion of signal intensities
            across all samples in any experimental condition that must have
            non-missing, within-threshold, valid values in order to keep the
            row for each gene feature
        proportion_signal_control (float): proportion of signal intensities
            across samples in control experimental condition that must have
            non-missing, within-threshold, valid values in order to keep the
            row for each gene feature
        proportion_signal_intervention (float): proportion of values of signal
            intensities across samples in intervention experimental condition
            that must have non-missing, within-threshold, valid values in
            order to keep the row for each gene feature
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table_filter = table_main.copy(deep=True)
    # Copy other information.
    samples_all = copy.deepcopy(samples_all)
    samples_control = copy.deepcopy(samples_control)
    samples_intervention = copy.deepcopy(samples_intervention)

    # Collect names of temporary columns.
    names_columns_temporary = list()

    ##########
    # Filter rows in table on basis of signal validity.

    # Filter rows in table by signal validity across all samples.
    # porg.test_extract_organize_values_from_series()
    table_filter["match_keep_signal_all"] = table_filter.apply(
        lambda row:
            porg.determine_series_signal_validity_threshold(
                series=row,
                keys_signal=samples_all,
                threshold_low=threshold_signal_low,
                threshold_high=threshold_signal_high,
                proportion=proportion_signal_all,
                report=False,
            ),
        axis="columns", # apply function to each row
    )
    table_filter = table_filter.loc[
        (table_filter["match_keep_signal_all"] == 1), :
    ]
    names_columns_temporary.append("match_keep_signal_all")

    # Filter rows in table by signal validity across samples in control or
    # intervention experimental conditions, respectively.
    if filter_rows_signal_by_condition:
        table_filter["match_keep_signal_control"] = table_filter.apply(
            lambda row:
                porg.determine_series_signal_validity_threshold(
                    series=row,
                    keys_signal=samples_control,
                    threshold_low=threshold_signal_low,
                    threshold_high=threshold_signal_high,
                    proportion=proportion_signal_control,
                    report=False,
                ),
            axis="columns", # apply function to each row
        )
        table_filter["match_keep_signal_intervention"] = table_filter.apply(
            lambda row:
                porg.determine_series_signal_validity_threshold(
                    series=row,
                    keys_signal=samples_intervention,
                    threshold_low=threshold_signal_low,
                    threshold_high=threshold_signal_high,
                    proportion=proportion_signal_intervention,
                    report=False,
                ),
            axis="columns", # apply function to each row
        )
        table_filter = table_filter.loc[
            (
                (table_filter["match_keep_signal_control"] == 1) |
                (table_filter["match_keep_signal_intervention"] == 1)
            ), :
        ]
        names_columns_temporary.append("match_keep_signal_control")
        names_columns_temporary.append("match_keep_signal_intervention")
        pass

    # Remove unnecessary columns.
    table_filter.drop(
        labels=names_columns_temporary,
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: filter_table_main_rows_signal()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_filter


def filter_table_main(
    table_main=None,
    columns_gene=None,
    samples_all=None,
    samples_control=None,
    samples_intervention=None,
    filter_rows_identity=None,
    remove_sex_chromosomes=None,
    filter_rows_signal=None,
    filter_rows_signal_by_condition=None,
    threshold_signal_low=None,
    threshold_signal_high=None,
    proportion_signal_all=None,
    proportion_signal_control=None,
    proportion_signal_intervention=None,
    tissue=None,
    report=None,
):
    """
    Filters information in table.

    This function filters rows for gene features by different thresholds on the
    proportion of a gene's signal intensities that must have non-missing,
    within-threshold, valid values across sets of samples in either the control
    or intervention experimental conditions or both. The reason for using these
    different thresholds is to allow genes to have signals that are detectable
    in only one experimental condition or the other. An example would be an
    experimental intervention that activates or deactivates transcriptional
    expression of specific genes.

    arguments:
        table_main (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows, with a few additional columns for attributes
            of gene features
        columns_gene (list<str>): names of columns corresponding to
            information about genes
        samples_all (list<str>): identifiers of samples in both control and
            intervention experimental conditions corresponding to names of
            columns for measurement values of signal intensity across features
        samples_control (list<str>): identifiers of samples in control
            experimental condition corresponding to names of columns for
            measurement values of signal intensity across features
        samples_intervention (list<str>): identifiers of samples in
            intervention experimental condition corresponding to names of
            columns for measurement values of signal intensity across features
        filter_rows_identity (bool): whether to filter rows by identity
        remove_sex_chromosomes (bool): whether to remove all genes on sex
            chromosomes
        filter_rows_signal (bool): whether to filter rows by signal
        filter_rows_signal_by_condition (bool): whether to filter rows by
            signal with separate consideration of columns for different
            experimental conditions
        threshold_signal_low (float): threshold below which
            (value <= threshold) all values are considered invalid and missing
        threshold_signal_high (float): threshold above which
            (value > threshold) all values are considered invalid and missing
        proportion_signal_all (float): proportion of signal intensities
            across all samples in any experimental condition that must have
            non-missing, within-threshold, valid values in order to keep the
            row for each gene feature
        proportion_signal_control (float): proportion of signal intensities
            across samples in control experimental condition that must have
            non-missing, within-threshold, valid values in order to keep the
            row for each gene feature
        proportion_signal_intervention (float): proportion of values of signal
            intensities across samples in intervention experimental condition
            that must have non-missing, within-threshold, valid values in
            order to keep the row for each gene feature
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_filter = table_main.copy(deep=True)
    # Copy other information.
    columns_gene = copy.deepcopy(columns_gene)
    samples_all = copy.deepcopy(samples_all)
    samples_control = copy.deepcopy(samples_control)
    samples_intervention = copy.deepcopy(samples_intervention)

    # Filter information in table.

    ##########
    # Filter rows within table on basis of gene feature identity.
    if filter_rows_identity:
        types_gene = define_keep_types_gene_narrow()
        chromosomes = define_keep_gene_chromosomes(
            remove_sex_chromosomes=remove_sex_chromosomes,
        )
        table_filter["match_keep_identity"] = table_filter.apply(
            lambda row:
                determine_keep_series_by_identity(
                    row_identifier=row["identifier_gene"],
                    row_type=row["gene_type"],
                    row_chromosome=row["gene_chromosome"],
                    identifier_prefix="ENSG",
                    types_gene=types_gene,
                    chromosomes=chromosomes,
                ),
            axis="columns", # apply function to each row
        )
        table_filter = table_filter.loc[
            (table_filter["match_keep_identity"] == 1), :
        ]
        # Remove unnecessary columns.
        table_filter.drop(
            labels=["match_keep_identity",],
            axis="columns",
            inplace=True
        )
        pass

    ##########
    # Filter rows within table on basis of signal validity.
    if filter_rows_signal:
        table_filter = filter_table_main_rows_signal(
            table_main=table_filter,
            samples_all=samples_all,
            samples_control=samples_control,
            samples_intervention=samples_intervention,
            filter_rows_signal_by_condition=filter_rows_signal_by_condition,
            threshold_signal_low=threshold_signal_low,
            threshold_signal_high=threshold_signal_high,
            proportion_signal_all=proportion_signal_all,
            proportion_signal_control=proportion_signal_control,
            proportion_signal_intervention=proportion_signal_intervention,
            tissue=tissue,
            report=report,
        )
        pass

    ##########
    # Filter and sort columns within table.
    columns_sequence = copy.deepcopy(columns_gene)
    #columns_sequence.extend(samples_control)
    #columns_sequence.extend(samples_intervention)
    columns_sequence.extend(samples_all)
    table_filter = porg.filter_sort_table_columns(
        table=table_filter,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: filter_table_main()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("source main table:")
        count_rows = (table_main.shape[0])
        count_columns = (table_main.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        print(table_main)
        putly.print_terminal_partition(level=5)
        print("product main table:")
        print(table_filter)
        putly.print_terminal_partition(level=5)
        count_rows = (table_filter.shape[0])
        count_columns = (table_filter.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("Count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=5)
        print("description of categorical gene type:")
        print(table_filter["gene_type"].describe(include=["category",]))
        putly.print_terminal_partition(level=5)
        print(
            "counts of genes with each unique categorical value of "
            + "gene type:")
        print(table_filter["gene_type"].value_counts(dropna=False))
        putly.print_terminal_partition(level=5)
        print(
            "counts of genes with each unique categorical value of "
            + "gene chromosome:")
        print(table_filter["gene_chromosome"].value_counts(dropna=False))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_filter


##########
# 6. Separate tables for information of distinct types.


def separate_table_main_columns(
    table_main=None,
    columns_gene=None,
    columns_signal=None,
    tissue=None,
    report=None,
):
    """
    Splits or separates information in table by columns.

    arguments:
        table_main (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows, with a few additional columns for attributes
            of gene features
        columns_gene (list<str>): names of columns corresponding to
            information about genes
        columns_signal (list<str>): names of columns corresponding to
            identifiers of samples in both control and intervention
            experimental conditions for measurement values of signal intensity
            across features
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_split = table_main.copy(deep=True)
    # Copy other information.
    columns_gene = copy.deepcopy(columns_gene)
    columns_signal = copy.deepcopy(columns_signal)
    # Separate information into separate tables.
    #columns_gene.remove("identifier_gene")
    columns_signal.insert(0, "identifier_gene")
    table_gene = table_split.loc[
        :, table_split.columns.isin(columns_gene)
    ].copy(deep=True)
    table_signal = table_split.loc[
        :, table_split.columns.isin(columns_signal)
    ].copy(deep=True)
    # Translate names of columns.
    #translations = dict()
    #translations["identifier_gene"] = "identifier"
    #table_gene.rename(
    #    columns=translations,
    #    inplace=True,
    #)
    # Organize indices in table.
    table_gene.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_gene.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    table_signal.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Organize indices in table.
    table_gene.columns.rename(
        "attribute",
        inplace=True,
    ) # single-dimensional index
    table_signal.columns.rename(
        "identifier_signal",
        inplace=True,
    ) # single-dimensional index

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: separate_table_main_columns()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("source main table:")
        count_rows = (table_split.shape[0])
        count_columns = (table_split.shape[1])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        print(table_main.iloc[0:10, 0:])
        putly.print_terminal_partition(level=5)
        count_rows = (table_gene.shape[0])
        count_columns = (table_gene.shape[1])
        print("table of information about genes:")
        print(table_gene.iloc[0:10, 0:])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=5)
        count_rows = (table_signal.shape[0])
        count_columns = (table_signal.shape[1])
        print("table of information about signals:")
        print(table_signal.iloc[0:10, 0:])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=5)
        pass
    # Collect information.
    pail = dict()
    pail["table_gene"] = table_gene
    pail["table_signal"] = table_signal
    # Return information.
    return pail


##########
# 7. Fill missing values of signal intensity.
# Functionality for this operation is in the "partner" package.


##########
# 8. Scale and normalize values of signal intensity across genes within each
# sample.


def scale_normalize_values_intensity_signal_table(
    table=None,
    method=None,
    report=None,
):
    """
    Adjust the scale and normalize the distribution of values of signal
    intensity across genes within each sample.

    Table's format and orientation

    Table has values of signal intensity for each gene oriented across rows
    with their samples oriented across columns.

    The table must have a single-level index (potentially with name
    "observations") across columns and a single-level index (potentially with
    name "features") across rows.

    identifier_signal   sample_1 sample_2 sample_3 sample_4 sample_5
    identifier_gene
    gene_1              ...      ...      ...      ...      ...
    gene_2              ...      ...      ...      ...      ...
    gene_3              ...      ...      ...      ...      ...
    gene_4              ...      ...      ...      ...      ...
    gene_5              ...      ...      ...      ...      ...

    arguments:
        table (object): Pandas data-frame table of values of intensity for
            samples across columns and for genes across rows
        method (str): name of method to use for scaling values, currently only
            'deseq'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values of intensity across
            samples in columns and proteins in rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Organize indices in table.
    table = porg.translate_names_table_indices_columns_rows(
        table=table,
        index_columns_product="observations",
        index_rows_source="identifier_gene",
        index_rows_product="features",
        report=False,
    )
    # Determine method for scaling.
    if (method == "deseq"):
        table_scale = (
            pscl.scale_feature_values_between_observations_by_deseq(
                table=table,
                name_columns="observations",
                name_rows="features",
                report=report,
        ))
    # Organize indices in table.
    table_scale = porg.translate_names_table_indices_columns_rows(
        table=table_scale,
        index_columns_product="identifier_signal",
        index_rows_source="features",
        index_rows_product="identifier_gene",
        report=False,
    )
    # Calculate the natural logarithm of signal intensity values.
    # math.log() # optimal for scalar values
    # numpy.log() # optimal for array values
    #table_normal = table_scale.map(
    #    lambda value: math.log(value),
    #    na_action="ignore", # ignore missing values in calculation
    #)
    table_normal = table_scale.transform(
        lambda row: numpy.log(row.to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )),
        axis="columns", # apply function to each row
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: scale_values_intensity_signal_table()")
        putly.print_terminal_partition(level=4)
        print(table_normal)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_normal


##########
# 9. Combine signals across samples and genes with supplemental information
# about genes.


def combine_organize_table_signal_genes_samples(
    table_signal=None,
    table_gene=None,
    table_sample=None,
    columns_signal=None,
    columns_gene=None,
    tissue=None,
    report=None,
):
    """
    Splits or separates information in table by columns.

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows
        table_gene (object): Pandas data-frame table of supplemental or
            contextual information about genes
        table_sample (object): Pandas data-frame table of information about
            samples
        columns_signal (list<str>): names of columns corresponding to
            identifiers of samples in both control and intervention
            experimental conditions for measurement values of signal intensity
            across features
        columns_gene (list<str>): names of columns corresponding to
            information about genes
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    table_gene = table_gene.copy(deep=True)
    table_sample = table_sample.copy(deep=True)
    # Copy other information.
    columns_signal = copy.deepcopy(columns_signal)
    columns_gene = copy.deepcopy(columns_gene)

    # Merge tables.
    table_merge = porg.merge_columns_two_tables(
        identifier_first="identifier_gene",
        identifier_second="identifier_gene",
        table_first=table_signal,
        table_second=table_gene,
        preserve_index=True,
        report=report,
    )

    # Filter rows in table for non-missing values across relevant columns.
    table_merge.dropna(
        axis="index",
        how="all",
        subset=columns_signal,
        inplace=True,
    )
    # Extract identifiers of relevant samples.
    table_sample = table_sample.loc[(
        table_sample["identifier_signal"].isin(columns_signal)
    ), :].copy(deep=True)
    columns_sample = copy.deepcopy(
        table_sample["identifier_sample"].to_list()
    )
    # Translate names of columns corresponding to samples.
    table_sample = table_sample.filter(
        items=["identifier_signal", "identifier_sample",],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_sample["identifier_sample"].to_list(),
        index=table_sample["identifier_signal"],
    )
    translations_column = series_translations.to_dict()
    table_merge.rename(
        columns=translations_column,
        inplace=True,
    )

    # Organize indices in table.
    table_merge.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_merge.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Organize indices in table.
    table_merge.columns.rename(
        "features",
        inplace=True,
    ) # single-dimensional index

    # Determine sort order from designations of chromosomes.
    table_merge["gene_chromosome_sort"] = table_merge.apply(
        lambda row: (
            "chr23" if ("X" in str(row["gene_chromosome"]))
            else row["gene_chromosome"]
        ),
        axis="columns", # apply function to each row
    )
    table_merge["gene_chromosome_sort"] = table_merge.apply(
        lambda row: (
            "chr24" if ("Y" in str(row["gene_chromosome_sort"]))
            else row["gene_chromosome_sort"]
        ),
        axis="columns", # apply function to each row
    )
    table_merge["gene_chromosome_sort"] = table_merge.apply(
        lambda row: (
            "chr25" if ("M" in str(row["gene_chromosome_sort"]))
            else row["gene_chromosome_sort"]
        ),
        axis="columns", # apply function to each row
    )
    table_merge["gene_chromosome_sort"] = table_merge.apply(
        lambda row: int(
            str(row["gene_chromosome_sort"]).strip().replace("chr", "")
        ),
        axis="columns", # apply function to each row
    )

    # Sort rows within table.
    table_merge.sort_values(
        by=[
            "gene_chromosome_sort",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Filter and sort columns within table.
    columns_gene.insert(len(columns_gene), "gene_chromosome_sort")
    columns_sequence = copy.deepcopy(columns_gene)
    columns_sequence.extend(columns_sample)
    columns_sequence.remove("identifier_gene")
    table_merge = porg.filter_sort_table_columns(
        table=table_merge,
        columns_sequence=columns_sequence,
        report=report,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: combine_organize_table_signal_genes_samples()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        count_rows = (table_merge.shape[0])
        count_columns = (table_merge.shape[1])
        print("table of information about signals across genes and samples:")
        print(table_merge.iloc[0:10, 0:])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_merge


def rank_genes_by_mean_signal(
    identifiers_gene=None,
    table_signal=None,
    report=None,
):
    """
    Ranks a list of genes by the magnitude of their mean signals.

    arguments:
        identifiers_gene (list<str): identifiers of genes
        table_signal (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    # Organize indices in table.
    table_signal.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Filter rows by identifiers of stable genes.
    table_signal = table_signal.loc[(
        table_signal["identifier_gene"].isin(identifiers_gene)
    ), :].copy(deep=True)
    # Organize indices in table.
    table_signal.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Calculate mean of values for each feature across all observations.
    table_signal["mean"] = table_signal.apply(
        lambda row:
            numpy.nanmean(
                row.to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                )
            ),
        axis="columns", # apply function to each row
    )
    # Calculate geometric mean of values for each feature across all
    # observations.
    table_signal["mean_geometric"] = table_signal.apply(
        lambda row:
            scipy.stats.mstats.gmean(
                row.to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                ),
                nan_policy="omit",
            ),
        axis="columns", # apply function to each row
    )
    # Sort rows within table.
    table_signal.sort_values(
        by=[
            "mean_geometric",
        ],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Organize indices in table.
    table_signal.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Extract identifiers of genes.
    identifiers_gene_rank = copy.deepcopy(
        table_signal["identifier_gene"].to_list()
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: rank_genes_by_mean_signal()")
        putly.print_terminal_partition(level=5)
        print("table in rank order")
        print(table_signal)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return identifiers_gene_rank



##########
# 10. Check the coherence of separate tables for analysis.


def check_coherence_table_sample_table_signal(
    table_sample=None,
    table_signal=None,
    tissue=None,
    name_instance=None,
    report=None,
):
    """
    Checks the coherence, specifically identity and sequence of samples, of
    information in separate tables for sample attributes and gene signals
    across samples.

    arguments:
        table_sample (object): Pandas data-frame table of attributes for
            samples
        table_signal (object): Pandas data-frame table of values of signal
            intensity for sample observations across columns and for gene
            features across rows
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)
    table_sample_extract = table_sample.copy(deep=True)
    table_signal = table_signal.copy(deep=True)
    # Organize indices in table.
    table_sample_extract.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Extract identifiers of samples from each separate table.
    samples_sample = copy.deepcopy(
        table_sample_extract["identifier_signal"].to_list()
    )
    samples_signal = copy.deepcopy(
        table_signal.columns.to_list()
    )
    # Confirm that both sets of samples are inclusive.
    inclusion = putly.compare_lists_by_mutual_inclusion(
        list_primary=samples_sample,
        list_secondary=samples_signal,
    )
    # Confirm that both sets of samples are identical across sequence.
    identity = putly.compare_lists_by_elemental_identity(
        list_primary=samples_sample,
        list_secondary=samples_signal,
    )
    # Confirm that both sets of samples are equal.
    equality = (samples_sample == samples_signal)
    # Perform tests of comparisons.
    test_primary = ["a", "b", "c", "d", "e", "f", "g",]
    test_secondary = ["a", "b", "c", "d", "g", "f", "e",]
    test_inclusion = putly.compare_lists_by_mutual_inclusion(
        list_primary=test_primary,
        list_secondary=test_secondary,
    )
    test_identity = putly.compare_lists_by_elemental_identity(
        list_primary=test_primary,
        list_secondary=test_secondary,
    )
    test_equality = (test_primary == test_secondary)
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: check_coherence_table_sample_table_signal()")
        print("tissue: " + tissue)
        print("name_instance: " + name_instance)
        putly.print_terminal_partition(level=5)
        count_rows = (table_sample.shape[0])
        count_columns = (table_sample.shape[1])
        print("table of information about samples:")
        print(table_sample.iloc[0:10, 0:])
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=5)
        count_rows = (table_signal.shape[0])
        count_columns = (table_signal.shape[1])
        print("table of information about signals between genes and samples:")
        print(table_signal)
        print("count of rows in table: " + str(count_rows))
        print("count of columns in table: " + str(count_columns))
        putly.print_terminal_partition(level=5)
        print("sample-signal identifiers from rows of sample table:")
        print(samples_sample)
        putly.print_terminal_partition(level=5)
        print("sample-signal identifiers from columns of signal table:")
        print(samples_signal)
        putly.print_terminal_partition(level=5)
        print("real comparisons of sample-signal identifiers:")
        print("inclusion: " + str(inclusion))
        print("identity: " + str(identity))
        print("equality: " + str(equality))
        putly.print_terminal_partition(level=5)
        print("fake comparisons to check the tests themselves:")
        print("inclusion: " + str(test_inclusion))
        print("identity: " + str(test_identity))
        print("equality: " + str(test_equality))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    pass


##########
# 11. Create charts to represent distribution of signal values for individual
# gene features across sample observations both with and without adjustment
# of scale and normalization between sample observations.


def create_write_chart_feature_signal_observations_distribution(
    identifiers_gene=None,
    colors_names_groups=None,
    tissue=None,
    table_gene=None,
    table_raw=None,
    table_scale=None,
    column_identifier=None,
    column_name=None,
    paths=None,
    report=None,
):
    """
    Create chart representation of distribution of signal values.

    arguments:
        identifiers_gene (list<str>): identifier of single gene feature
        colors_names_groups (list<str>): names of colors for groups
        tissue (list<str>): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        table_gene (object): Pandas data-frame table of information about genes
        table_raw (object): Pandas data-frame table of signals without
            adjustment of scale or normalization
        table_scale (object): Pandas data-frame table of signals with
            adjustment of scale or normalization
        column_identifier (str): name of column in table for the unique
            identifier corresponding to the fold change
        column_name (str): name of column in table for the name corresponding
            to the fold change
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Organize parameters.
    name_figure = str("chart_gene_signal_scale_normal_" + tissue)
    #path_directory = paths["out_plot"]
    path_directory = paths["out_whole_plot"]

    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()

    # Copy information in table.
    table_gene = table_gene.copy(deep=True)
    table_raw = table_raw.copy(deep=True)
    table_scale = table_scale.copy(deep=True)

    # Extract parameters for colors.
    colors_groups = list(map(
        lambda color_name: copy.deepcopy(colors[color_name]),
        colors_names_groups
    ))
    # Collect information.
    names_groups = list()
    values_groups = list()
    # Iterate on features.
    for identifier_gene in identifiers_gene:
        # Extract name of gene.
        name_gene = str(table_gene.at[identifier_gene, "gene_name"])
        # Iterate on tables of signal values.
        types = ["raw", "scale",]
        for type in types:
            # Prepare name for group.
            name_group = str(name_gene + "_" + type)
            # Extract values.
            if (type == "raw"):
                series = table_raw.loc[identifier_gene]
                values = numpy.log(series.dropna().to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                ))
            elif (type == "scale"):
                series = table_scale.loc[identifier_gene]
                values = series.dropna().to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                )
            # Collect information.
            names_groups.append(name_group)
            values_groups.append(values)
        pass
    # Specify colors.

    # Create figure.
    figure = pplot.plot_boxes_groups(
        values_groups=values_groups,
        title_ordinate="log(gene signal)",
        title_abscissa="",
        titles_abscissa_groups=names_groups,
        colors_groups=colors_groups,
        label_top_center="",
        label_top_left="",
        label_top_right="",
        aspect="landscape",
        orientation_box="vertical",
        axis_linear_minimum=0.0,
        fonts=fonts,
        colors=colors,
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




##########
# TODO: TCW; 18 July 2024
# For analysis (by regression, for example) it will make most sense to transpose the "signal" table
# to orient samples across rows (convenience in merging in sample attributes)
# and genes across columns (features, along with sample attributes as covariates).



##########
# Downstream analyses

# differential expression in DESeq2
# regression analyses in custom models

##########
# Downstream visualizations

# volcano plots of fold changes
# Venn diagrams of sets of differentially expressed genes (up or down)
# heatmaps of gene signals across persons
# pairwise correlation matrices with hierarchical clustering
#



###############################################################################
# Procedure


##########
# Control procedure for whole signal data for each tissue type.


def control_procedure_whole_trunk_preparation(
    tissue=None,
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # 1. Initialize directories for read of source and write of product files.
    # Initialize directories for trunk procedure.
    paths = initialize_directories_trunk(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=False,
        report=report,
    )

    ##########
    # 2. Read source information from file.
    pail_source_sample = read_source_sample(
        paths=paths,
        report=report,
    )
    # Extract identifiers of relevant samples.
    table_sample = pail_source_sample["table_sample"].loc[(
        pail_source_sample["table_sample"]["tissue"] == tissue
    ), :].copy(deep=True)
    samples = copy.deepcopy(
        table_sample["identifier_signal"].to_list()
    )

    ##########
    # 3. Prepare data for signals across genes and samples with stratification
    # for a specific instance of analysis.
    pail_signal = control_procedure_part_branch_signal(
        samples=samples,
        tissue=tissue, # adipose, muscle
        scale_combination=True,
        paths=paths,
        report=report,
    )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_table = dict()
    pail_write_table[str("table_gene_" + tissue)] = (
        pail_signal["table_gene"]
    )
    pail_write_table[str("table_signal_" + tissue)] = (
        pail_signal["table_signal"]
    )
    pail_write_table[str("table_signal_scale_" + tissue)] = (
        pail_signal["table_signal_scale"]
    )
    pail_write_table[str("table_signal_gene_sample_" + tissue)] = (
        pail_signal["table_combination"]
    )

    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_table,
        path_directory=paths["out_whole_preparation"],
        reset_index=False,
        write_index=True,
        type="text",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_table,
        path_directory=paths["out_whole_preparation"],
        reset_index=False,
        write_index=True,
        type="pickle",
    )
    pass


def control_procedure_whole_trunk_description(
    tissue=None,
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # 1. Initialize directories for read of source and write of product files.
    # Initialize directories for trunk procedure.
    paths = initialize_directories_trunk(
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        restore=False,
        report=report,
    )

    ##########
    # 2. Read source information from file.
    pail_source = read_source_signal_for_description(
        tissue=tissue,
        paths=paths,
        report=report,
    )

    ##########
    # 3. Extract identifiers of common gene features that have the most stable
    # signal values across all sample observations both with and without
    # adjustment of scale and normalization between sample observations.
    pail_stable = porg.compare_stable_feature_sets_before_after_scale(
        table_raw=pail_source["table_signal"],
        table_scale=pail_source["table_signal_scale"],
        name_columns="identifier_signal",
        name_rows="identifier_gene",
        count_quantile=21,
        report=report,
    )
    #identifiers_genes_stable = list(pail_stable["union_middle"])
    identifiers_gene_stable_rank = rank_genes_by_mean_signal(
        identifiers_gene=list(pail_stable["union_stable"]),
        table_signal=pail_source["table_signal_scale"],
        report=report,
    )
    #count_half = int(len(identifiers_gene_stable_rank) // 2)
    #genes_selection = identifiers_gene_stable_rank[count_half:(count_half+5)]
    genes_selection = identifiers_gene_stable_rank[3:8]
    print(genes_selection)

    ##########
    # 4. For a selection of gene features, describe their distribution of
    # signal values across sample observations both with and without adjustment
    # of scale and normalization between sample observations.
    identifiers_gene = genes_selection
    colors_names_groups = [
        "red_burgundy",
        "red_crimson",
        "yellow_sunshine",
        "yellow_sunflower",
        "green_forest",
        "green_kelly",
        "blue_navy",
        "blue_sky",
        "purple",
        "purple_lavender",
    ]
    # Create and write to file charts to represent distribution of
    # signals.
    create_write_chart_feature_signal_observations_distribution(
        identifiers_gene=identifiers_gene,
        colors_names_groups=colors_names_groups,
        tissue=tissue,
        table_gene=pail_source["table_gene"],
        table_raw=pail_source["table_signal"],
        table_scale=pail_source["table_signal_scale"],
        column_identifier="identifier_gene",
        column_name="gene_name",
        paths=paths,
        report=report,
    )
    pass


##########
# Control procedure within branch for parallelization.


def control_procedure_part_branch_sample(
    sort=None,
    group=None,
    name_instance=None,
    tissue=None,
    cohort_selection_primary=None,
    factor_availability=None,
    cohort_selection_secondary=None,
    continuity_scale=None,
    columns_set=None,
    project=None,
    routine=None,
    procedure=None,
    paths=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        sort (int): sequential index for sort order
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        cohort_selection_primary (dict<list<str>>): filters on rows in table
            for selection of samples relevant to cohort for analysis
        factor_availability (dict<list<str>>): features and their values
            corresponding to factors of interest in the analysis
        cohort_selection_secondary (dict<list<str>>): filters on rows in table
            for selection of samples relevant to cohort for analysis
        continuity_scale (list<str>): names of columns for covariates with
            values on continuous scale of measurement, interval or ratio,
            for which to standardize the scale by z score
        columns_set (list<str>): names of columns for feature variables that
            are relevant to the current set or instance of parameters
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Read source information from file.
    pail_source_sample = read_source_sample(
        paths=paths,
        report=report,
    )
    # Select set of samples relevant for analysis.
    pail_sample_primary = select_sets_identifier_table_sample(
        table_sample=pail_source_sample["table_sample"],
        name_instance=name_instance,
        tissue=tissue,
        cohort_selection=cohort_selection_primary,
        factor_availability=factor_availability,
        report=report,
    )
    # Organize tertiles for feature variables with continuous values that the
    # parameters specify.
    table_tertile = organize_describe_summarize_table_sample_tertiles(
        table=pail_sample_primary["table_selection"],
        cohort_selection=cohort_selection_secondary,
        group=group,
        name_instance=name_instance,
        tissue=tissue,
        paths=paths,
        report=report,
    )
    # Select sets of samples on the basis of tertiles.
    pail_sample_secondary = select_sets_identifier_table_sample(
        table_sample=table_tertile,
        name_instance=name_instance,
        tissue=tissue,
        cohort_selection=cohort_selection_secondary,
        factor_availability=factor_availability,
        report=report,
    )
    # Filter rows in table for non-missing values across relevant columns.
    pail_sample = select_sets_final_identifier_table_sample(
        table_sample=pail_sample_secondary["table_selection"],
        name_instance=name_instance,
        tissue=tissue,
        continuity_scale=continuity_scale,
        columns_set=columns_set,
        report=report,
    )
    # Summarize the counts of samples corresponding to each set of parameters.
    report_write_count_samples(
        samples=pail_sample["samples_selection"],
        sort=sort,
        group=group,
        name_instance=name_instance,
        tissue=tissue,
        paths=paths,
        report=report,
    )
    # Return information.
    return pail_sample


def control_procedure_part_branch_signal(
    samples=None,
    tissue=None,
    scale_combination=None,
    paths=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        samples (list<str>): identifiers of samples corresponding to
            measurement values of signal intensity that are relevant
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        scale_combination (bool): whether to adjust the scale of signals and
            combine with supplemental information about genes
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Read source information from file.
    pail_source_sample = read_source_sample(
        paths=paths,
        report=report,
    )
    pail_source_main = read_source_main(
        tissue=tissue,
        paths=paths,
        report=report,
    )

    # Organize from source the information about signals.
    columns_gene = define_column_sequence_table_main_gene()
    pail_organization_main = organize_table_main(
        table_main=pail_source_main["table_main"],
        columns_gene=columns_gene,
        samples=samples,
        tissue=tissue,
        report=report,
    )

    # Filter columns and rows in main table.
    table_filter = filter_table_main(
        table_main=pail_organization_main["table_main"],
        columns_gene=columns_gene,
        samples_all=samples,
        samples_control=[],
        samples_intervention=[],
        filter_rows_identity=True,
        remove_sex_chromosomes=True,
        filter_rows_signal=True,
        filter_rows_signal_by_condition=False, # separate cases from controls
        threshold_signal_low=10, # 10 is DESeq2 recommendation for bulk RNAseq
        threshold_signal_high=None,
        proportion_signal_all=0.33, # proportion of smaller condition to total
        proportion_signal_control=0.5, # inactive
        proportion_signal_intervention=0.5, # inactive
        tissue=tissue,
        report=report,
    )

    # Separate tables for information of distinct types.
    pail_separate = separate_table_main_columns(
        table_main=table_filter,
        columns_gene=columns_gene,
        columns_signal=samples,
        tissue=tissue,
        report=report,
    )

    # Fill missing values of signal intensity.
    table_signal = porg.fill_missing_values_table_by_row(
        table=pail_separate["table_signal"],
        columns=samples,
        method="zero",
        report=report,
    )

    # Determine whether to adjust the scale of signals and combine with
    # supplemental information about genes.
    if scale_combination:
        # Scale values of signal intensity across genes within each sample.
        # The goal of this scaling is to make the individual samples more
        # comparable to each other.
        # This scaling can decrease the variance or noise in measurements between
        # samples.
        # Normalize distribution values of signal intensity by taking the
        # natural logarithm.
        table_signal_scale = scale_normalize_values_intensity_signal_table(
            table=table_signal,
            method="deseq",
            report=report,
        )
        # Combine and organize signals across samples and genes with supplemental
        # information about genes.
        table_combination = combine_organize_table_signal_genes_samples(
            table_signal=table_signal_scale,
            table_gene=pail_separate["table_gene"],
            table_sample=pail_source_sample["table_sample_file"],
            columns_signal=samples,
            columns_gene=columns_gene,
            tissue=tissue,
            report=report,
        )
    else:
        table_signal_scale = pandas.DataFrame()
        table_combination = pandas.DataFrame()
        pass

    ##########
    # Collect information.
    pail = dict()
    pail["table_gene"] = pail_separate["table_gene"]
    pail["table_signal"] = table_signal
    pail["table_signal_scale"] = table_signal_scale
    pail["table_combination"] = table_combination
    # Return information.
    return pail


def control_procedure_part_branch(
    sort=None,
    group=None,
    name_instance=None,
    tissue=None,
    cohort_selection_primary=None,
    factor_availability=None,
    cohort_selection_secondary=None,
    continuity_scale=None,
    columns_set=None,
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        sort (int): sequential index for sort order
        group (str): name of a group of analyses
        name_instance (str): name of instance set of parameters for selection
            of samples in cohort and definition of analysis
        tissue (str): name of tissue that distinguishes study design and
            set of relevant samples, either 'adipose' or 'muscle'
        cohort_selection_primary (dict<list<str>>): filters on rows in table
            for selection of samples relevant to cohort for analysis
        factor_availability (dict<list<str>>): features and their values
            corresponding to factors of interest in the analysis
        cohort_selection_secondary (dict<list<str>>): filters on rows in table
            for selection of samples relevant to cohort for analysis
        continuity_scale (list<str>): names of columns for covariates with
            values on continuous scale of measurement, interval or ratio,
            for which to standardize the scale by z score
        columns_set (list<str>): names of columns for feature variables that
            are relevant to the current set or instance of parameters
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # 1. Initialize directories for read of source and write of product files.
    # Initialize directories for instance-specific parallel branch procedure.
    paths = initialize_directories_branch_instance(
        project=project,
        routine=routine,
        procedure=procedure,
        tissue=tissue,
        name_instance=name_instance,
        path_directory_dock=path_directory_dock,
        restore=True,
        report=report,
    )

    ##########
    # 2. Prepare data for subjects and samples with stratification for a
    # specific instance of analysis.
    pail_sample = control_procedure_part_branch_sample(
        sort=sort,
        group=group,
        name_instance=name_instance,
        tissue=tissue, # adipose, muscle
        cohort_selection_primary=cohort_selection_primary,
        factor_availability=factor_availability,
        cohort_selection_secondary=cohort_selection_secondary,
        continuity_scale=continuity_scale,
        columns_set=columns_set,
        project=project,
        routine=routine,
        procedure=procedure,
        paths=paths,
        report=report,
    )

    ##########
    # 3. Prepare data for signals across genes and samples with stratification
    # for a specific instance of analysis.
    pail_signal = control_procedure_part_branch_signal(
        samples=pail_sample["samples_selection"],
        tissue=tissue, # adipose, muscle
        scale_combination=False,
        paths=paths,
        report=report,
    )

    ##########
    # 4. Check the coherence of separate tables for analysis.
    check_coherence_table_sample_table_signal(
        table_sample=pail_sample["table_selection"],
        table_signal=pail_signal["table_signal"],
        tissue=tissue,
        name_instance=name_instance,
        report=report,
    )


    # TODO: TCW; 24 July 2024
    # TODO:
    # 1. - In future, add a pseudo-count of 1 to any zeros or missing values.
    #    - Silverman et al, 2020 (PubMed:33101615) is a good review of methods for handling zeros.
    # 2. - In future, apply a variety of scale-adjustment and normalization methods.
    #    - DESeq2 algorithm (PubMed:25516281), edgeR algorithm (PubMed:20196867), MRN algorithm (PubMed:26442135)
    #    - Evans, 2018 (PubMed:28334202)
    # 3. - In future, apply logarithmic transformation and evaluate distributions (histograms).
    # 4. - In future, evaluate coefficient of variance of each gene across control samples.

    ##########
    # Collect information.
    # Collections of files.
    pail_write_data = dict()
    pail_write_data[str("table_sample")] = (
        pail_sample["table_selection"]
    )
    pail_write_data[str("table_gene")] = (
        pail_signal["table_gene"]
    )
    pail_write_data[str("table_signal")] = (
        pail_signal["table_signal"]
    )

    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_data,
        path_directory=paths["out_data"],
        reset_index=False,
        write_index=True,
        type="text",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_data,
        path_directory=paths["out_data"],
        reset_index=False,
        write_index=True,
        type="pickle",
    )
    pass


##########
# Manage parallelization.


def control_parallel_instance(
    instance=None,
    parameters=None,
):
    """
    Control procedure to organize within tables the information about genetic
    correlations from LDSC.

    arguments:
        instance (dict): parameters specific to current instance
            sort (int): sequential index for sort order
            group (str): name of a group of analyses
            name_instance (str): name of instance set of parameters for
                selection of samples in cohort and definition of analysis
            tissue (str): name of tissue that distinguishes study design and
                set of relevant samples, either 'adipose' or 'muscle'
            cohort_selection_primary (dict<list<str>>): filters on rows in
                table for selection of samples relevant to cohort for analysis
            factor_availability (dict<list<str>>): features and their values
                corresponding to factors of interest in the analysis
            cohort_selection_secondary (dict<list<str>>): filters on rows in
                table for selection of samples relevant to cohort for analysis
            continuity_scale (list<str>): names of columns for covariates with
                values on continuous scale of measurement, interval or ratio,
                for which to standardize the scale by z score
            columns_set (list<str>): names of columns for feature variables
                that are relevant to the current set or instance of parameters
        parameters (dict): parameters common to all instances
            project (str): name of project
            routine (str): name of routine, either 'transcriptomics' or
                'proteomics'
            procedure (str): name of procedure, a set or step in the routine
                process
            path_directory_dock (str): path to dock directory for procedure's
                source and product directories and files
            report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Extract parameters.
    # Extract parameters specific to each instance.
    sort = instance["sort"]
    group = instance["group"]
    name_instance = instance["name_instance"]
    tissue = instance["tissue"]
    cohort_selection_primary = instance["cohort_selection_primary"]
    factor_availability = instance["factor_availability"]
    cohort_selection_secondary = instance["cohort_selection_secondary"]
    continuity_scale = instance["continuity_scale"]
    columns_set = instance["columns_set"]
    # Extract parameters common across all instances.
    project = parameters["project"]
    routine = parameters["routine"]
    procedure = parameters["procedure"]
    path_directory_dock = parameters["path_directory_dock"]
    report = parameters["report"]

    ##########
    # Control procedure with split branch for parallelization.
    control_procedure_part_branch(
        sort=sort,
        group=group,
        name_instance=name_instance,
        tissue=tissue, # adipose, muscle
        cohort_selection_primary=cohort_selection_primary,
        factor_availability=factor_availability,
        cohort_selection_secondary=cohort_selection_secondary,
        continuity_scale=continuity_scale,
        columns_set=columns_set,
        project=project,
        routine=routine,
        procedure=procedure,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    pass


def collect_scrap_parallel_instances_for_analysis_sets(
):
    """
    Collect scrap parallel instances for analysis sets.

    arguments:

    raises:

    returns:

    """

    # Collect parameters specific to each instance.
    # tissue: adipose
    instances = [
        {
            "name_instance": str(
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
            "name_instance": str(
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
    # tissue: muscle
    instances = [
        {
            "name_instance": str(
                "muscle_exercise-0hr_age"
            ),
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
            "name_instance": str(
                "muscle_younger_exercise"
            ),
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
            "name_instance": str(
                "muscle_elder_exercise"
            ),
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
    ]
    pass


def control_parallel_instances(
    instances=None,
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control procedure for parallel instances.

    arguments:
        instances (list<dict>): parameters to control individual instances in
            parallel
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Collect parameters common across all instances.
    parameters = dict()
    parameters["project"] = project
    parameters["routine"] = routine
    parameters["procedure"] = procedure
    parameters["path_directory_dock"] = path_directory_dock
    parameters["report"] = report

    # Execute procedure iteratively with parallelization across instances.
    if True:
        prall.drive_procedure_parallel(
            function_control=(
                control_parallel_instance
            ),
            instances=instances,
            parameters=parameters,
            cores=5,
            report=True,
        )
    else:
        # Execute procedure directly for testing.
        control_parallel_instance(
            instance=instances[0],
            parameters=parameters,
        )
    pass


##########
# Call main procedure.


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
    routine="transcriptomics"
    procedure="organize_signal"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("routine: " + str(routine))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Trunk procedure to prepare tables of signals with adjustment of scale
    # and normalization.
    if True:
        # Initialize directories.
        paths = initialize_directories_trunk(
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            restore=True,
            report=report,
        )
        # Control procedure for preparation of signal data as a whole.
        control_procedure_whole_trunk_preparation(
            tissue="muscle",
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            report=report,
        )
        control_procedure_whole_trunk_preparation(
            tissue="adipose",
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            report=report,
        )

    ##########
    # Trunk procedure to describe tables of signals with adjustment of scale
    # and normalization.
    if True:
        # Control procedure for description of signal data as a whole.
        control_procedure_whole_trunk_description(
            tissue="muscle",
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            report=report,
        )
        control_procedure_whole_trunk_description(
            tissue="adipose",
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            report=report,
        )
        pass

    ##########
    # Branch procedure to prepare tables of signals without adjustment of scale
    # or normalization for selections of samples in specific instances of
    # analysis.

    # The current implementation requires manual switch on or off according to
    # the tissues that are included in the batch.
    if True:
        # Initialize directories before branch procedure.
        paths = initialize_directories_before_branch(
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            restore=True,
            report=report,
        )
        paths_muscle = initialize_directories_branch_tissue(
            project=project,
            routine=routine,
            procedure=procedure,
            tissue="muscle",
            path_directory_dock=path_directory_dock,
            restore=True,
            report=report,
        )
        paths_adipose = initialize_directories_branch_tissue(
            project=project,
            routine=routine,
            procedure=procedure,
            tissue="adipose",
            path_directory_dock=path_directory_dock,
            restore=True,
            report=report,
        )

        ##########
        # Read and organize parameters for parallel instances.
        instances = read_organize_source_parameter_instances(
            paths=paths,
            report=report,
        )
        ##########
        # Control procedure for parallel instances.
        control_parallel_instances(
            instances=instances,
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            report=report,
        )
        ##########
        # Organize summary information about all instances overall.
        if False:
            read_organize_write_summary_instances_tissue(
                tissue="muscle",
                paths=paths_muscle,
                report=report,
            )
        if True:
            read_organize_write_summary_instances_tissue(
                tissue="adipose",
                paths=paths_adipose,
                report=report,
            )

    pass


###############################################################################
# End
