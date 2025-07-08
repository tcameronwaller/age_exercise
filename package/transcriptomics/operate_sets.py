"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'operate_sets' is part of the 'transcriptomics' package
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
#import partner.regression as preg
import partner.plot as pplot
import partner.parallelization as prall
import age_exercise.phenotypes.organize_subject as aexph_sub
import age_exercise.phenotypes.organize_sample as extr_sample
import age_exercise.transcriptomics.organize_signal as extr_signal
import age_exercise.transcriptomics.select_gene_sets as extr_select


###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.
# There is a hierarchy in these functions to initialize directories to manage
# the hierarchical tree structure of sub-procedures.


##########
# 2. Read source information from file.



##########
# 3. Organize information from source.



###############################################################################
# Procedure


##########
# Execute main procedure at highest level of hierarchy.


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
    procedure="operate_sets"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: age_exercise.transcriptomics.operate_sets.py")
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

    ##########
    # Summary of available sets.
    # 2.1. Read and count unique genes in sets.
    #path_directory_sets_gene = os.path.join(
    #    paths["in_parameters_private"], project, routine, "sets_gene",
    #    "sets_gene_deseq2_2025-06-05",
    #)
    path_directory_sets_gene = os.path.join(
        paths["out_project"], "transcriptomics", "select_gene_sets", "data",
        "sets_gene",
    )
    table_counts_sets_gene = (
        putly.read_child_files_text_list_count_unique_items(
            path_directory=path_directory_sets_gene,
            name_file_prefix="",
            name_file_suffix=".txt",
            name_file_not="",
            report=report,
    ))

    ##########
    # Read and extract identifiers of genes in sets.

    # age
    genes_adipose_1 = aexph_sub.read_extract_set_features(
        name_set="adipose_1_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_1_down = aexph_sub.read_extract_set_features(
        name_set="adipose_1_genes_down",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_1_up = aexph_sub.read_extract_set_features(
        name_set="adipose_1_genes_up",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )

    # placebo
    genes_adipose_14 = aexph_sub.read_extract_set_features(
        name_set="adipose_14_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_27 = aexph_sub.read_extract_set_features(
        name_set="adipose_27_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_31 = aexph_sub.read_extract_set_features(
        name_set="adipose_31_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )

    genes_adipose_15 = aexph_sub.read_extract_set_features(
        name_set="adipose_15_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_15_down = aexph_sub.read_extract_set_features(
        name_set="adipose_15_genes_down",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_15_up = aexph_sub.read_extract_set_features(
        name_set="adipose_15_genes_up",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_28 = aexph_sub.read_extract_set_features(
        name_set="adipose_28_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_28_down = aexph_sub.read_extract_set_features(
        name_set="adipose_28_genes_down",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_28_up = aexph_sub.read_extract_set_features(
        name_set="adipose_28_genes_up",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_32 = aexph_sub.read_extract_set_features(
        name_set="adipose_32_genes_change",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_32_down = aexph_sub.read_extract_set_features(
        name_set="adipose_32_genes_down",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )
    genes_adipose_32_up = aexph_sub.read_extract_set_features(
        name_set="adipose_32_genes_up",
        suffix_file=".txt",
        path_directory=path_directory_sets_gene,
        report=report,
    )

    ##########
    # Combine the identifiers of genes in sets.
    # Age.
    print("age union")
    genes_age = aexph_sub.combine_sets_items_union_unique(
        sets_items=[
            genes_adipose_1,
        ],
        report=True,
    )
    # Notice that the genes from "adipose_25" and "adipose_32" that associated
    # with omega-3 already had integral control for placebo; however, the genes
    # from "adipose_12" need exclusion of genes from "adipose_11" that
    # associated with placebo.
    # Placebo.
    print("placebo union")
    genes_placebo_union = aexph_sub.combine_sets_items_union_unique(
        sets_items=[
            genes_adipose_14,
            genes_adipose_27,
            genes_adipose_31,
        ],
        report=True,
    )

    # Omega-3.
    print("omega-3 union triple: 15, 28, 32")
    genes_omega3_union_triple = aexph_sub.combine_sets_items_union_unique(
        sets_items=[
            genes_adipose_15,
            genes_adipose_28,
            genes_adipose_32,
        ],
        report=True,
    )
    genes_omega3_union_triple_down = aexph_sub.combine_sets_items_union_unique(
        sets_items=[
            genes_adipose_15_down,
            genes_adipose_28_down,
            genes_adipose_32_down,
        ],
        report=True,
    )
    genes_omega3_union_triple_up = aexph_sub.combine_sets_items_union_unique(
        sets_items=[
            genes_adipose_15_up,
            genes_adipose_28_up,
            genes_adipose_32_up,
        ],
        report=True,
    )

    print("omega-3 union double: 28, 32")
    genes_omega3_union_double = aexph_sub.combine_sets_items_union_unique(
        sets_items=[
            genes_adipose_28,
            genes_adipose_32,
        ],
        report=True,
    )

    # Omega-3 correction for placebo.
    print("adipose_28 difference adipose_27")
    genes_adipose_28_not_27 = aexph_sub.combine_sets_items_difference_unique(
        items_inclusion=genes_adipose_28,
        items_exclusion=genes_adipose_27,
        report=True,
    )
    genes_adipose_28_down_not_27 = aexph_sub.combine_sets_items_difference_unique(
        items_inclusion=genes_adipose_28_down,
        items_exclusion=genes_adipose_27,
        report=True,
    )
    genes_adipose_28_up_not_27 = aexph_sub.combine_sets_items_difference_unique(
        items_inclusion=genes_adipose_28_up,
        items_exclusion=genes_adipose_27,
        report=True,
    )

    ##########
    # Combine by difference of the identifiers of genes in sets.
    # This approach is unnecessarily conservative.
    # Some genes might change in one direction in placebo and in the other
    # direction in omega-3 groups.
    print("omega-3 difference placebo")
    genes_omega3_not_placebo = aexph_sub.combine_sets_items_difference_unique(
        items_inclusion=genes_omega3_union_triple,
        items_exclusion=genes_placebo_union,
        report=True,
    )
    genes_omega3_down_not_placebo = aexph_sub.combine_sets_items_difference_unique(
        items_inclusion=genes_omega3_union_triple_down,
        items_exclusion=genes_placebo_union,
        report=True,
    )
    genes_omega3_up_not_placebo = aexph_sub.combine_sets_items_difference_unique(
        items_inclusion=genes_omega3_union_triple_up,
        items_exclusion=genes_placebo_union,
        report=True,
    )

    ##########
    # Combine by intersection of the identifiers of genes in sets.
    print("age intersection omega-3 triple")
    genes_age_and_omega3_not_placebo = (
        aexph_sub.combine_sets_items_intersection_unique(
            items_first=genes_age,
            items_second=genes_omega3_not_placebo,
            report=True,
    ))
    genes_age_and_28_omega3_not_placebo = (
        aexph_sub.combine_sets_items_intersection_unique(
            items_first=genes_age,
            items_second=genes_adipose_28_not_27,
            report=True,
    ))
    # adipose_1 (age) and adipose_28 (omega-3) without adipose_27 (placebo)
    genes_adipose_1_down_28_down_not_27 = (
        aexph_sub.combine_sets_items_intersection_unique(
            items_first=genes_adipose_1_down,
            items_second=genes_adipose_28_down_not_27,
            report=True,
    ))
    genes_adipose_1_down_28_up_not_27 = (
        aexph_sub.combine_sets_items_intersection_unique(
            items_first=genes_adipose_1_down,
            items_second=genes_adipose_28_up_not_27,
            report=True,
    ))
    genes_adipose_1_up_28_down_not_27 = (
        aexph_sub.combine_sets_items_intersection_unique(
            items_first=genes_adipose_1_up,
            items_second=genes_adipose_28_down_not_27,
            report=True,
    ))
    genes_adipose_1_up_28_up_not_27 = (
        aexph_sub.combine_sets_items_intersection_unique(
            items_first=genes_adipose_1_up,
            items_second=genes_adipose_28_up_not_27,
            report=True,
    ))

    ##########
    # Write information to file.
    # Collect information.
    # Collections of files.
    pail_write_lists = dict()
    pail_write_lists["genes_adipose_28_not_27"] = (
        genes_adipose_28_not_27
    )
    pail_write_lists["genes_adipose_28_down_not_27"] = (
        genes_adipose_28_down_not_27
    )
    pail_write_lists["genes_adipose_28_up_not_27"] = (
        genes_adipose_28_up_not_27
    )
    pail_write_lists["genes_omega3_not_placebo"] = (
        genes_omega3_not_placebo
    )
    pail_write_lists["genes_omega3_down_not_placebo"] = (
        genes_omega3_down_not_placebo
    )
    pail_write_lists["genes_omega3_up_not_placebo"] = (
        genes_omega3_up_not_placebo
    )
    pail_write_lists["genes_age_and_omega3_not_placebo"] = (
        genes_age_and_omega3_not_placebo
    )
    pail_write_lists["genes_age_and_28_omega3_not_placebo"] = (
        genes_age_and_28_omega3_not_placebo
    )

    # adipose_1 (age) and adipose_28 (omega-3) without adipose_27 (placebo)
    pail_write_lists["genes_adipose_1_down_28_down_not_27"] = (
        genes_adipose_1_down_28_down_not_27
    )
    pail_write_lists["genes_adipose_1_down_28_up_not_27"] = (
        genes_adipose_1_down_28_up_not_27
    )
    pail_write_lists["genes_adipose_1_up_28_down_not_27"] = (
        genes_adipose_1_up_28_down_not_27
    )
    pail_write_lists["genes_adipose_1_up_28_up_not_27"] = (
        genes_adipose_1_up_28_up_not_27
    )


    # Define paths to directories.
    # dock/out_age_exercise/transcriptomics/operate_sets/lists
    path_directory_write = os.path.join(
        paths["out_procedure_lists"],
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_write,
    )
    # Write product information to file.
    putly.write_lists_to_file_text(
        pail_write=pail_write_lists,
        path_directory=path_directory_write,
        delimiter="\n",
    )


    ###################
    # Scrap
    ##################


    if False:
        ##########
        # Subjective combination of sets of genes identified by threshold.

        # TODO: TCW; 9 January 2025
        # Organize operation below within a cleaner, more versatile function
        # arguments:
        # path_directory_read
        # names_sets <-- corresponding to names of files within directory
        # name_file_product
        # path_directory_write

        # Define paths to directories.
        path_directory_comparison_sets = os.path.join(
            paths["out_routine"], "compare_sets_groups", "data",
            "sets_subjective_match_directionality_younger_elder_omega3",
        )
        sets_match_younger_elder_omega3 = dict()
        sets_match_younger_elder_omega3["main"] = [
            "3_adipose_omega3_permissive",
            "4_adipose_omega3",
            "5_adipose_age_placebo_omega3",
        ]
        set_match_younger_elder_omega3 = (
            aexph_sub.read_extract_combine_custom_feature_sets(
                names_sets=sets_match_younger_elder_omega3,
                features_available=None,
                path_directory=path_directory_comparison_sets,
                report=True,
        ))
        # Copy information.
        genes_match_younger_elder_omega3 = copy.deepcopy(
            set_match_younger_elder_omega3["main"]
        )
        # Collect unique names of genes in set.
        genes_match_younger_elder_omega3 = putly.collect_unique_elements(
            elements=genes_match_younger_elder_omega3,
        )
        # Collect information.
        # Collections of files.
        name_set = str(
            "genes_subjective_match_directionality_younger_elder_omega3"
        )
        pail_write_lists = dict()
        pail_write_lists[name_set] = (
            genes_match_younger_elder_omega3
        )

    pass


###############################################################################
# End
