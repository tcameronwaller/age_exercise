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
# Initialization


def initialize_directories(
    project=None,
    routine=None,
    procedure=None,
    tissues=None,
    path_directory_dock=None,
    restore=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        tissues (list<str>): names of tissue that distinguish study design and
            sets of samples
        path_directory_dock (str): path to dock directory for source and
            product directories and files
        restore (bool): whether to remove previous versions of data

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
    # Initialize directories in main branch.
    paths_initialization = [
        paths["out_project"],
        paths["out_routine"],
        paths["out_procedure"],
    ]
    # Remove previous directories and files to avoid version or batch
    # confusion.
    if restore:
        for path in paths_initialization:
            putly.remove_directory(path=path)
    # Create directories.
    for path in paths_initialization:
        putly.create_directories(
            path=path,
        )
    # Initialize directories in branch forks.
    for tissue in tissues:
        paths["out_tissue"] = os.path.join(
            paths["out_procedure"], str(tissue),
        )
        #paths[str(str(tissue) + "_out_test")] = os.path.join(
        #    paths["out_tissue"], "test",
        #)
        paths[str(str(tissue) + "_out_data")] = os.path.join(
            paths["out_tissue"], "data",
        )
        #paths[str(str(tissue) + "_out_plot")] = os.path.join(
        #    paths["out_tissue"], "plot",
        #)
        # Initialize directories.
        paths_initialization = [
            paths["out_tissue"],
            #paths[str(str(tissue) + "_out_test")],
            paths[str(str(tissue) + "_out_data")],
            #paths[str(str(tissue) + "_out_plot")],
        ]
        # Remove previous directories and files to avoid version or batch
        # confusion.
        if restore:
            for path in paths_initialization:
                putly.remove_directory(path=path)
        # Create directories.
        for path in paths_initialization:
            putly.create_directories(
                path=path,
            )
    # Return information.
    return paths


##########
# 1. Read source information from file.


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
# 2. Organize table of matches between samples and files.


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
        "identifier",
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
        (dict<object>): collection of information

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
# 3. Organize table of attributes for samples.


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
    translations["Date"] = "date_visit"
    translations["Age Group"] = "cohort_age_letter"
    translations["Sex"] = "sex_letter"
    translations["Intervention"] = "intervention"
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
        "cohort_age_letter",
        "intervention",
        "subject_attribute",
        "study_clinic_visit_relative",
        "date_visit",
        "sex_letter",
        #"sex_text",
        #"sex_y",
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
        (dict<object>): collection of information

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
            "cohort_age_letter",
            "intervention",
            "subject_attribute",
            "study_clinic_visit_relative",
        ],
        axis="index",
        ascending=True,
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
# 2. Organize information from source.


def define_column_sequence_table_main_gene():
    """
    Defines the columns in sequence within table.

    arguments:

    raises:

    returns:
        (list<str>): variable types of columns within table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        "identifier_gene",
        "gene_identifier",
        "gene_name",
        #"gene_exon_number",
        "gene_type",
        "gene_chromosome",
        #"Unnamed: 160",
    ]
    # Return information.
    return columns_sequence


# dateutil.parser.parser


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
        print("module: exercise.transcriptomics.organization.py")
        print("function: organize_table_main()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("main table: ")
        print(table_main)
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





###############################################################################
# Procedure


##########
# Control procedure with split for parallelization.


def control_branch_procedure(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Control branch of procedure.

    arguments:
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:


    """

    ##########
    # 2. Organize information from source.
    columns_sample = define_column_sequence_table_sample()
    columns_gene = define_column_sequence_table_main_gene()

    pail_organization_supplement = organize_table_sample_supplement(
        table_sample=pail_source["table_sample"],
        columns_sample=columns_sample,
        tissue=tissue,
        report=report,
    )



    pail_organization_sample = organize_table_sample(
        table_sample=pail_source["table_sample"],
        columns_sample=columns_sample,
        tissue=tissue,
        report=report,
    )
    pail_organization_main = organize_table_main(
        table_main=pail_source["table_main"],
        columns_gene=columns_gene,
        samples=pail_organization_sample["samples_selection"],
        tissue=tissue,
        report=report,
    )

    ##########
    # 3. Filter columns and rows in main table.
    table_filter = filter_table_main(
        table_main=pail_organization_main["table_main"],
        columns_gene=columns_gene,
        samples_all=pail_organization_sample["samples_selection"],
        samples_control=pail_organization_sample["samples_control"],
        samples_intervention=(
            pail_organization_sample["samples_intervention_1"]
        ),
        filter_rows_identity=True,
        filter_rows_signal=True,
        filter_rows_signal_by_condition=True,
        threshold_signal_low=10, # DESeq2 recommendation for bulk RNAseq
        threshold_signal_high=None,
        proportion_signal_all=0.10, # proportion smaller condition to total
        proportion_signal_control=0.5,
        proportion_signal_intervention=0.5,
        tissue=tissue,
        report=report,
    )

    ##########
    # 4. Separate tables for information of distinct types.
    pail_separate = separate_table_main_columns(
        table_main=table_filter,
        columns_gene=columns_gene,
        columns_signal=pail_organization_sample["samples_selection"],
        tissue=tissue,
        report=report,
    )

    ##########
    # 5. Fill missing values of signal intensity.
    table_signal = porg.fill_missing_values_table_by_row(
        table=pail_separate["table_signal"],
        columns=pail_organization_sample["samples_selection"],
        method="zero",
        report=report,
    )


    ##########
    # 6. Check the coherence of separate tables for analysis.
    check_coherence_table_sample_table_signal(
        table_sample=pail_organization_sample["table_sample_selection"],
        table_signal=table_signal,
        tissue=tissue,
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
        pail_organization_sample["table_sample_selection"]
    )
    pail_write_data[str("table_gene")] = (
        pail_separate["table_gene"]
    )
    pail_write_data[str("table_signal")] = (
        pail_separate["table_signal"]
    )

    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_data,
        path_directory=paths[str(str(tissue) + "_out_data")],
        reset_index=False,
        write_index=True,
        type="text",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_data,
        path_directory=paths[str(str(tissue) + "_out_data")],
        reset_index=False,
        write_index=True,
        type="pickle",
    )
    pass


def execute_procedure(
    path_directory_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_dock (str): path to dock directory for source and
            product directories and files

    raises:

    returns:

    """

    ##########
    # Report.
    print("system: local")
    print("project: exercise")
    print("routine: transcriptomics")
    print("procedure: organize_sample")

    ##########
    # Initialize directories.
    paths = initialize_directories(
        project="exercise",
        routine="transcriptomics",
        procedure="organize_sample",
        tissues=["muscle", "adipose",],
        path_directory_dock=path_directory_dock,
        restore=True,
    )

    ##########
    # 1. Read source information from file.
    pail_source = read_source(
        paths=paths,
        report=True,
    )

    ##########
    # 2. Organize table of matches between samples and files.
    translations_sample_file = define_translation_columns_table_sample_file()
    columns_sample_file = define_sequence_columns_table_sample_file()
    table_sample_file = organize_table_sample_file(
        table=pail_source["table_sample_file"],
        translations_column=translations_sample_file,
        columns_sequence=columns_sample_file,
        report=True,
    )

    ##########
    # 3. Organize table of attributes for samples.
    translations_sample_attribute = (
        define_translation_columns_table_sample_attribute()
    )
    columns_sample_attribute = define_sequence_columns_table_sample_attribute()
    table_sample_attribute = organize_table_sample_attribute(
        table=pail_source["table_sample_attribute"],
        translations_column=translations_sample_attribute,
        columns_sequence=columns_sample_attribute,
        report=True,
    )





    if False:
        ##########
        # Control procedure branches with split by tissue type.
        control_branch_procedure(
            tissue="adipose", # adipose, muscle
            paths=paths,
            report=True,
        )
        control_branch_procedure(
            tissue="muscle", # adipose, muscle
            paths=paths,
            report=True,
        )

    pass


###############################################################################
# End
