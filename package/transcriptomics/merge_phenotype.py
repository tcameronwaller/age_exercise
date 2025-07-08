"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'merge_phenotype' is part of the 'transcriptomics' subpackage
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

# For subjects in a study cohort, this procedure merges together phenotype
# features with a small number of gene features. The procedure then organizes
# and prepares the table for analysis. The advantage of this approach is that
# the product table is accessible for external review and for the design of
# specific, custom analyses. The disadvantage of this approach is that it would
# be less efficient to merge together phenotype features with the signals of
# all 15,000 or so protein-coding genes. For an alternative approach, refer to
# the "regress_genes.py" module within the "transcriptomics" subpackage of the
# "age_exercise" package.

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
import matplotlib

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.regression as preg
import partner.parallelization as prall
import age_exercise.phenotypes.organize_subject as aexph_sub
import age_exercise.phenotypes.organize_sample as exph_sample

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.


##########
# 2. Read source information from file.


def read_source_subject(
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
    path_file_table_subject = os.path.join(
        paths["out_project"], "phenotypes", "organize_subject", "tables",
        "table_subject.pickle",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of properties for subjects.
    pail["table_subject"] = pandas.read_pickle(
        path_file_table_subject,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.merge_phenotype.py")
        print("function: read_source_sample()")
        putly.print_terminal_partition(level=5)
        print("subject table: ")
        print(pail["table_subject"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    path_file_table_sample = os.path.join(
        paths["out_project"], "phenotypes", "organize_sample", "tables",
        "table_sample.pickle",
    )

    # Collect information.
    pail = dict()
    # Read information from file.

    # Table of samples and their attributes.
    pail["table_sample"] = pandas.read_pickle(
        path_file_table_sample,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.merge_phenotype.py")
        print("function: read_source_sample()")
        putly.print_terminal_partition(level=5)
        print("sample table: ")
        print(pail["table_sample"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_gene(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (str): name of tissue, either 'muscle' or 'adipose', that
            distinguishes the study design and its relevant data
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
    path_file_table_gene = os.path.join(
        paths["out_routine"], "organize_signal", "whole", "preparation",
        str("table_gene_" + tissue + ".pickle"),
    )

    # Read information from file.

    # Table of properties and attributes of genes.
    table = pandas.read_pickle(
        path_file_table_gene,
    )
    # Organize indices of table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index

    # Collect information.
    pail = dict()
    pail["table_gene"] = table

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.merge_phenotype.py")
        print("function: read_source_gene()")
        putly.print_terminal_partition(level=5)
        print("gene table: ")
        print(pail["table_gene"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_signal(
    tissue=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        tissue (str): name of tissue, either 'muscle' or 'adipose', that
            distinguishes the study design and its relevant data
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
    path_file_table_signal = os.path.join(
        paths["out_project"], "transcriptomics", "organize_signal", "whole",
        "preparation", str("table_signal_scale_normal_" + tissue + ".pickle"),
    )

    # Read information from file.

    # Table of signals for genes across samples.
    table = pandas.read_pickle(
        path_file_table_signal,
    )
    # Organize indices of table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index

    # Collect information.
    pail = dict()
    pail["table_signal"] = table

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.merge_phenotype.py")
        print("function: read_source_signal()")
        putly.print_terminal_partition(level=5)
        print("signal table: ")
        print(pail["table_signal"])
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

    # Read source information for subjects.
    pail_subject = read_source_subject(
        paths=paths,
        report=False,
    )

    # Read source information for samples.
    pail_sample = read_source_sample(
        paths=paths,
        report=False,
    )

    # Read source information for genes.
    pail_gene_adipose = read_source_gene(
        tissue="adipose",
        paths=paths,
        report=False,
    )
    pail_gene_muscle = read_source_gene(
        tissue="muscle",
        paths=paths,
        report=False,
    )

    # Read source information for signals.
    pail_signal_adipose = read_source_signal(
        tissue="adipose",
        paths=paths,
        report=False,
    )
    pail_signal_muscle = read_source_signal(
        tissue="muscle",
        paths=paths,
        report=False,
    )

    # Collect information.
    pail = dict()
    pail["table_subject"] = pail_subject["table_subject"]
    pail["table_sample"] = pail_sample["table_sample"]
    pail["table_gene_adipose"] = pail_gene_adipose["table_gene"]
    pail["table_gene_muscle"] = pail_gene_muscle["table_gene"]
    pail["table_signal_adipose"] = pail_signal_adipose["table_signal"]
    pail["table_signal_muscle"] = pail_signal_muscle["table_signal"]

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: merge_phenotype.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 3. Merge information for phenotypes of samples with signals of genes.


def define_sequence_columns_priority():
    """
    Define the sequence of columns in table for a priority selection of
    features.

    Review: TCW; 8 May 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Define sequence of columns in table.
    columns_sequence = [
        "inclusion",
        "identifier_subject",
        "identifier_sample",
        "identifier_signal",
        "age",
        "age_cohort_text",
        "age_cohort_elder",
        "sex_text",
        "sex_female",
        "sex_y",
        "visit_text",
        "visit_second",
        "tissue",
        "intervention_text",
        "intervention_omega3",
        "exercise_duration_text",
    ]
    # Return information.
    return columns_sequence


def merge_organize_table_signal_tissues(
    table_signal_adipose=None,
    table_signal_muscle=None,
    table_gene_adipose=None,
    table_gene_muscle=None,
    selection_genes=None,
    prefix_name_gene=None,
    report=None,
):
    """
    Merge together tables of information for phenotypes of samples and signals
    oif genes.

    arguments:
        table_signal_adipose (object): Pandas data-frame table
        table_signal_muscle (object): Pandas data-frame table
        table_gene_adipose (object): Pandas data-frame table
        table_gene_muscle (object): Pandas data-frame table
        selection_genes (list<str>): identifiers of genes for which to include
            signals across samples
        prefix_name_gene (str): prefix for name of columns in table for genes
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table_signal_adipose = table_signal_adipose.copy(deep=True)
    table_signal_muscle = table_signal_muscle.copy(deep=True)
    table_gene_adipose = table_gene_adipose.copy(deep=True)
    table_gene_muscle = table_gene_muscle.copy(deep=True)
    table_gene = table_gene_muscle.copy(deep=True)
    selection_genes = copy.deepcopy(selection_genes)

    # Filter genes according to selection.
    if (len(selection_genes) > 0):
        table_gene = table_gene.loc[
            table_gene["identifier_gene"].isin(selection_genes), :
        ].copy(deep=True)
        pass
    table_gene = table_gene.loc[
        (table_gene["gene_name"].str.len() > 0), :
    ].copy(deep=True)
    table_gene = table_gene.loc[
        (table_gene["gene_identifier_base"].str.len() > 0), :
    ].copy(deep=True)
    # Extract information for translation of names of columns.
    table_translations = table_gene.filter(
        items=["gene_identifier_base", "gene_name",],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_translations["gene_name"].to_list(),
        index=table_translations["gene_identifier_base"],
    )
    translations_gene = series_translations.to_dict()

    # Append prefix to names of columns in table for genes.
    for key in translations_gene.keys():
        name = str(translations_gene[key])
        name_prefix = str(prefix_name_gene + name)
        translations_gene[key] = name_prefix
        pass
    # Filter signals for selection of genes.
    if (len(selection_genes) > 0):
        table_signal_adipose = table_signal_adipose.loc[
            table_signal_adipose["identifier_gene"].isin(selection_genes), :
        ].copy(deep=True)
        table_signal_muscle = table_signal_muscle.loc[
            table_signal_muscle["identifier_gene"].isin(selection_genes), :
        ].copy(deep=True)
        pass

    # Remove unnecessary columns.
    if ("index" in table_signal_adipose.columns.to_list()):
        table_signal_adipose.drop(
            labels=["index",],
            axis="columns",
            inplace=True
        )
    if ("index" in table_signal_muscle.columns.to_list()):
        table_signal_muscle.drop(
            labels=["index",],
            axis="columns",
            inplace=True
        )
        pass

    # Merge phenotypes of subjects together with signals of genes.
    table_merge = porg.merge_columns_two_tables(
        identifier_first="identifier_gene",
        identifier_second="identifier_gene",
        table_first=table_signal_adipose,
        table_second=table_signal_muscle,
        preserve_index=False,
        report=report,
    )

    # Organize indices of table.
    table_merge = (
        porg.explicate_table_indices_columns_rows_single_level(
            table=table_merge,
            index_columns="identifier_signal",
            index_rows="identifier_gene",
            explicate_indices=True,
            report=report,
    ))
    # Transpose table.
    table_merge = table_merge.transpose(copy=True)
    # Organize indices of table.
    table_merge.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_merge.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index
    # Translate names of columns in table.
    table_merge.rename(
        columns=translations_gene,
        inplace=True,
    )

    # Collect information.
    pail = dict()
    pail["table_merge"] = table_merge
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: merge_phenotype.py")
        print("function: merge_table_signal_tissues()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def merge_organize_tables_phenotype_signal(
    table_sample=None,
    table_signal=None,
    selection_genes=None,
    columns_sequence_priority=None,
    report=None,
):
    """
    Merge together tables of information for phenotypes of samples and signals
    of genes.

    arguments:
        table_sample (object): Pandas data-frame table
        table_signal (object): Pandas data-frame table
        selection_genes (list<str>): identifiers of genes for which to include
            signals across samples
        columns_sequence_priority (list<str>): sequence of columns in table for
            a selection of priority features
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table_sample = table_sample.copy(deep=True)
    table_signal = table_signal.copy(deep=True)
    selection_genes = copy.deepcopy(selection_genes)

    # Filter signals for selection of genes.
    if (len(selection_genes) > 0):
        columns = list()
        columns.extend(selection_genes)
        columns.insert(0, "identifier_signal")
        table_signal = table_signal.loc[
            :, table_signal.columns.isin(columns)
        ].copy(deep=True)
        pass

    # Merge phenotypes of subjects together with signals of genes.
    table_merge = porg.merge_columns_two_tables(
        identifier_first="identifier_signal",
        identifier_second="identifier_signal",
        table_first=table_sample,
        table_second=table_signal,
        preserve_index=False,
        report=report,
    )

    # Sort sequence of columns in table.
    table_merge = porg.sort_table_columns_explicit_other(
        table=table_merge,
        columns_sequence=columns_sequence_priority,
        report=report,
    )

    # Sort sequence of rows in table.
    table_merge.sort_values(
        by=["age_cohort_text", "sex_text", "visit_text", "tissue",],
        axis="index",
        ascending=True,
        inplace=True,
    )

    # Collect information.
    pail = dict()
    pail["table_merge"] = table_merge
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: merge_phenotype.py")
        print("function: merge_tables_sample_gene_signal()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


##########
# 4. Calculate product interaction.


def calculate_product_terms_interaction_effect(
    table=None,
    features_first=None,
    features_second=None,
    delimiter_name=None,
    report=None,
):
    """
    Calculate product of quantitative values between pairs of features for
    evaluation of interaction effects in regression.

    arguments:
        table (object): Pandas data-frame table
        features_first (list<str>): identifiers or names of features for which
            to calculate products
        features_second (list<str>): identifiers or names of features for which
            to calculate products
        delimiter_name (str): delimiter to place between text items when
            combining identifiers or names of features
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    features_first = copy.deepcopy(features_first)
    features_second = copy.deepcopy(features_second)

    # Iterate across combinations of features.
    for feature_first in features_first:
        for feature_second in features_second:
            name = str(feature_first + delimiter_name + feature_second)
            table[name] = table.apply(
                lambda series_row: float(
                    series_row[feature_first] * series_row[feature_second]
                ),
                axis="columns", # apply function to each row
            )
            pass
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: merge_phenotype.py")
        print("function: calculate_product_terms_interaction_effect()")
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
    routine="transcriptomics"
    procedure="merge_phenotype"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise")
        print("subpackage: transcriptomics")
        print("module: merge_phenotype.py")
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

    # Read and count unique genes in sets.
    path_directory_sets_gene = os.path.join(
        paths["out_routine"], "operate_sets", "lists",
    )
    genes_age_omega3 = aexph_sub.read_extract_set_features(
        name_set="genes_age_and_28_omega3_not_placebo",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )

    ##########
    # 3. Merge signals of genes from samples in tissue of adipose and muscle.
    selection_genes = [
        "ENSG00000198431", # TXNRD1
        "ENSG00000265972", # TXNIP
        "ENSG00000066336", # SPI1
        "ENSG00000291237", # SOD2
        "ENSG00000142168", # SOD1
        "ENSG00000073169", # SELENOO
        "ENSG00000198832", # SELENOM
        "ENSG00000113811", # SELENOK
        "ENSG00000117592", # PRDX6
        "ENSG00000181019", # NQO1
        "ENSG00000086991", # NOX4
        "ENSG00000100365", # NCF4
        "ENSG00000116701", # NCF2
        "ENSG00000158517", # NCF1
        "ENSG00000167468", # GPX4
        "ENSG00000233276", # GPX1
        "ENSG00000001084", # GCLC
        "ENSG00000165168", # CYBB
        "ENSG00000051523", # CYBA
        "ENSG00000121691", # CAT
    ]
    pail_merge_signal = merge_organize_table_signal_tissues(
        table_signal_adipose=pail_source["table_signal_adipose"],
        table_signal_muscle=pail_source["table_signal_muscle"],
        table_gene_adipose=pail_source["table_gene_adipose"],
        table_gene_muscle=pail_source["table_gene_muscle"],
        selection_genes=selection_genes,
        prefix_name_gene="rnaseq_",
        report=report,
    )
    #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    #print(pail_merge_signal["table_merge"])
    #print(pail_merge_signal["table_merge"].columns.to_list())

    ##########
    # 4. Merge phenotypes for subjects or samples to signals of genes.
    selection_genes = [
        "rnaseq_TXNRD1",
        "rnaseq_TXNIP",
        "rnaseq_SPI1",
        "rnaseq_SOD2",
        "rnaseq_SOD1",
        "rnaseq_SELENOO",
        "rnaseq_SELENOM",
        "rnaseq_SELENOK",
        "rnaseq_PRDX6",
        "rnaseq_NQO1",
        "rnaseq_NOX4",
        "rnaseq_NCF4",
        "rnaseq_NCF2",
        "rnaseq_NCF1",
        "rnaseq_GPX4",
        "rnaseq_GPX1",
        "rnaseq_GCLC",
        "rnaseq_CYBB",
        "rnaseq_CYBA",
        "rnaseq_CAT",
    ]
    columns_sequence_priority = define_sequence_columns_priority()
    pail_merge = merge_organize_tables_phenotype_signal(
        table_sample=pail_source["table_sample"],
        table_signal=pail_merge_signal["table_merge"],
        selection_genes=selection_genes,
        columns_sequence_priority=columns_sequence_priority,
        report=report,
    )

    ##########
    # 5. Copy information.
    table_sample = pail_source["table_sample"].copy(deep=True)
    table_signal = pail_merge_signal["table_merge"].copy(deep=True)
    table_merge = pail_merge["table_merge"].copy(deep=True)
    table_merge["age_scale"] = table_merge["age"].copy(deep=True)

    ##########
    # 6. Adjust scale of signals for genes.
    features_continuity_scale = [
        "age_scale",
        "rnaseq_TXNRD1",
        "rnaseq_TXNIP",
        "rnaseq_SPI1",
        "rnaseq_SOD2",
        "rnaseq_SOD1",
        "rnaseq_SELENOO",
        "rnaseq_SELENOM",
        "rnaseq_SELENOK",
        "rnaseq_PRDX6",
        "rnaseq_NQO1",
        "rnaseq_NOX4",
        "rnaseq_NCF4",
        "rnaseq_NCF2",
        "rnaseq_NCF1",
        "rnaseq_GPX4",
        "rnaseq_GPX1",
        "rnaseq_GCLC",
        "rnaseq_CYBB",
        "rnaseq_CYBA",
        "rnaseq_CAT",
    ]
    table_merge = pscl.manage_transform_scale_feature_by_table_columns(
        table=table_merge,
        features_continuity_scale=features_continuity_scale,
        adjust_scale=True,
        method_scale="z_score",
        report=report,
    )

    ##########
    # 6. Calculate predictor terms for interaction effects.
    selection_genes = [
        "rnaseq_TXNRD1",
        "rnaseq_TXNIP",
        "rnaseq_SPI1",
        "rnaseq_SOD2",
        "rnaseq_SOD1",
        "rnaseq_SELENOO",
        "rnaseq_SELENOM",
        "rnaseq_SELENOK",
        "rnaseq_PRDX6",
        "rnaseq_NQO1",
        "rnaseq_NOX4",
        "rnaseq_NCF4",
        "rnaseq_NCF2",
        "rnaseq_NCF1",
        "rnaseq_GPX4",
        "rnaseq_GPX1",
        "rnaseq_GCLC",
        "rnaseq_CYBB",
        "rnaseq_CYBA",
        "rnaseq_CAT",
    ]
    pail_interaction = calculate_product_terms_interaction_effect(
        table=table_merge,
        features_first=selection_genes,
        features_second=["age_cohort_elder",],
        delimiter_name="_-_",
        report=report,
    )
    pail_interaction = calculate_product_terms_interaction_effect(
        table=pail_interaction["table"],
        features_first=selection_genes,
        features_second=["age_scale",],
        delimiter_name="_-_",
        report=report,
    )

    ##########
    # 7. Write information to file.
    # Collect information.
    # Collections of files.
    #pail_write_lists = dict()
    pail_write_tables = dict()
    pail_write_tables[str("table_sample")] = table_sample
    pail_write_tables[str("table_signal")] = table_signal
    pail_write_tables[str("table_merge")] = pail_interaction["table"]
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
