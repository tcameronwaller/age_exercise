"""
Supply functionality for process and analysis of data from proteomics using
mass spectroscopy.

This module 'organization' is part of the 'proteomics' package within the
'exercise' package.

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
    technology=None,
    set=None,
    path_directory_dock=None,
    restore=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
        technology (str): name of technology, either 'transcriptomics' or
            'proteomics'
        set (str): name of set or step in process procedure
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
        paths["dock"], "in_data", str(project), str(technology),
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters", str(project), str(technology),
    )
    paths["in_parameters_private"] = os.path.join(
        paths["dock"], "in_parameters_private", str(project), str(technology),
    )
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    paths["out_technology"] = os.path.join(
        paths["out_project"], str(technology),
    )
    paths["out_set"] = os.path.join(
        paths["out_technology"], str(set),
    )
    paths["out_test"] = os.path.join(
        paths["out_set"], "test",
    )
    paths["out_table"] = os.path.join(
        paths["out_set"], "table",
    )
    paths["out_plot"] = os.path.join(
        paths["out_set"], "plot",
    )
    paths_initialization = [
        paths["out_project"],
        paths["out_technology"],
        paths["out_set"],
        paths["out_test"],
        paths["out_table"],
        paths["out_plot"],
    ]
    # Remove previous files to avoid version or batch confusion.
    if restore:
        for path in paths_initialization:
            putly.remove_directory(path=path)
    # Initialize directories.
    for path in paths_initialization:
        putly.create_directories(
            path=path,
        )
    # Return information.
    return paths


##########
# 1. Read information from file.



def define_table_column_types_sample():
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
    types_columns["sample_plate"] = "string"
    types_columns["condition_code"] = "string"
    types_columns["tissue"] = "string"
    types_columns["condition"] = "string"
    types_columns["subject"] = "string"
    types_columns["note_condition"] = "string"
    #types_columns["sex"] = "string"
    #types_columns["age"] = "int32"
    #types_columns["body_mass"] = "float32"
    #types_columns["body_mass_index"] = "float32"
    # Return information.
    return types_columns


def define_table_column_types_main():
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


def read_source(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
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
    path_file_table_sample = os.path.join(
        paths["in_parameters_private"], "quantification_2024-07-14",
        "attributes_samples", "table_attributes_samples.tsv",
    )
    if (tissue == "adipose"):
        path_file_table_main = os.path.join(
            paths["in_data"], "quantification_2024-07-14",
            "organization", "quantification_rna_reads_gene_adipose.tsv",
        )
    elif (tissue == "muscle"):
        path_file_table_main = os.path.join(
            paths["in_data"], "quantification_2024-07-14",
            "organization", "quantification_rna_reads_gene_muscle.tsv",
        )
        pass

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of samples and their attributes.
    types_columns = define_table_column_types_sample()
    pail["table_sample"] = pandas.read_csv(
        path_file_table_sample,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )

    # Table of values of intensity across samples and proteins.
    types_columns = define_table_column_types_main()
    pail["table_main"] = pandas.read_csv(
        path_file_table_main,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organization.py")
        print("function: read_source()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("sample table: ")
        print(pail["table_sample"])
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


##########
# 2. Organize information.


def organize_table_sample(
    tissue=None,
    table_sample=None,
    report=None,
):
    """
    Organizes information in tables about samples and measurement signals.

    arguments:
        tissue (str): name of tissue, either 'adipose' or 'muscle', which
            distinguishes study design and sets of samples
        table_sample (object): Pandas data-frame table of information about
            samples represented in the main table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_sample = table_sample.copy(deep=True)

    # Organize indices in table.
    table_sample.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_sample.set_index(
        ["identifier"],
        append=False,
        drop=True,
        inplace=True,
    )
    table_sample.columns.rename(
        "attribute",
        inplace=True,
    ) # single-dimensional index

    # Filter and separate information about samples.
    table_sample_tissue = table_sample.loc[
        (table_sample["tissue"] == tissue), :
    ]
    table_sample_inclusion = table_sample_tissue.loc[
        (table_sample_tissue["inclusion"] == 1), :
    ]
    table_sample_control = table_sample_inclusion.loc[
        (table_sample_inclusion["condition"] == "control"), :
    ]
    table_sample_intervention_1 = table_sample_inclusion.loc[
        (table_sample_inclusion["condition"] == "intervention_1"), :
    ]
    table_sample_intervention_2 = table_sample_inclusion.loc[
        (table_sample_inclusion["condition"] == "intervention_2"), :
    ]
    # Extract identifiers of samples in separate groups.
    samples_tissue = copy.deepcopy(
        table_sample_tissue["identifier"].to_list()
    )
    samples_inclusion = copy.deepcopy(
        table_sample_inclusion["identifier"].to_list()
    )
    samples_control = copy.deepcopy(
        table_sample_control["identifier"].to_list()
    )
    samples_intervention_1 = copy.deepcopy(
        table_sample_intervention_1["identifier"].to_list()
    )
    samples_intervention_2 = copy.deepcopy(
        table_sample_intervention_2["identifier"].to_list()
    )
    # Collect information.
    pail = dict()
    pail["table_sample"] = table_sample
    pail["samples"] = samples_tissue
    pail["samples_inclusion"] = samples_inclusion
    pail["samples_control"] = samples_control
    pail["samples_intervention_1"] = samples_intervention_1
    pail["samples_intervention_2"] = samples_intervention_2
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organization.py")
        print("function: organize_table_sample()")
        print("tissue: " + tissue)
        putly.print_terminal_partition(level=5)
        print("sample table: ")
        print(table_sample)
        putly.print_terminal_partition(level=5)
        print("description of categorical experimental conditions:")
        print(table_sample_tissue["condition"].describe(include=["category",]))
        print(
            "counts of samples with each unique categorical value of "
            + "experimental condition:")
        print(table_sample_tissue["condition"].value_counts(dropna=False))
        putly.print_terminal_partition(level=5)
        for name_list in pail.keys():
            count_list = len(pail[name_list])
            print("count sample list " + name_list + ": " + str(count_list))
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_table_main(
    samples=None,
    table_main=None,
    report=None,
):
    """
    Organizes information in tables about samples and measurement signals.

    arguments:
        samples (list<str>): identifiers of samples corresponding to names of
            columns for measurement signals across features
        table_main (object): Pandas data-frame table of values of signal
            intensity for samples across columns and for genes across rows,
            with a few additional columns of attributes for genes
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table_main = table_main.copy(deep=True)

    # Translate names of columns.
    translations = dict()
    translations["gene_id"] = "gene_identifier"
    translations["exon_number"] = "gene_exon_number"
    translations["chromosome"] = "gene_chromosome"
    translations.update(translations_sample)
    table_main.rename(
        columns=translations,
        inplace=True,
    )
    # Organize indices in table.
    table_main.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_main.set_index(
        ["identifier_gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Replace values of zero for signal intensity with missing values.
    # Only replace values within table's columns for samples.
    # This implementation is more concise than iteration across specific
    # columns.
    table_main[samples] = table_main[samples].replace(
        to_replace=0,
        value=pandas.NA,
    )
    # Replace values less than zero with missing values.
    table_main[samples][table_main[samples] < 0] = pandas.NA

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
        print(
            "counts of genes with each unique categorical value of "
            + "gene type:")
        print(table_main["gene_type"].value_counts(dropna=False))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail



# TODO: TCW; 18 July 2024
# Filter rows:
# "ENSG" not it "identifier_gene"


# TODO: TCW; 18 July 2024
# For analysis (by regression, for example) it will make most sense to transpose the "signal" table
# to orient samples across rows (convenience in merging in sample attributes)
# and genes across columns (features, along with sample attributes as covariates).




###############################################################################
# Procedure


##########
# Control procedure with split for parallelization.

def control_split_procedure(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Organizes information in tables about samples, features, and measurement
    signals.

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
    # 1. Read source information from file.
    pail_source = read_source(
        tissue=tissue,
        paths=paths,
        report=report,
    )

    ##########
    # 2. Organize information from source.
    pail_organization_sample = organize_table_sample(
        tissue=tissue,
        table_sample=pail_source["table_sample"],
        report=report,
    )
    pail_organization_main = organize_table_main(
        table_main=pail_source["table_main"],
        samples=pail_organization_sample["samples"],
        report=report,
    )


    ##########
    # Write product information to file.
    ##########
    # Return information.
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
    print("technology: transcriptomics")
    print("procedure: 2_organization")
    print("set: organization")

    ##########
    # Initialize directories.
    paths = initialize_directories(
        project="exercise",
        technology="transcriptomics",
        set="organization",
        path_directory_dock=path_directory_dock,
        restore=True,
    )

    ##########
    # Control procedure with split for parallelization.
    control_split_procedure(
        tissue="adipose",
        paths=paths,
        report=True,
    )


    pass


###############################################################################
# End
