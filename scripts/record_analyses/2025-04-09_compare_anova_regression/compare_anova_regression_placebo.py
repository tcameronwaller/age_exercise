"""
Compare Analysis of Variance (ANOVA) against Linear Regression in the analysis
of published data from a placebo-controlled clinical trial.

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

    The code within project 'partner' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the GNU
    General Public License, or (at your option) any later version.

    The code within project 'partner' is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'partner'. If not, see <http://www.gnu.org/licenses/>.
"""

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, first execution: 7 April 2025
# Date, last execution: 9 April 2025
# Review: TCW; 9 April 2025
################################################################################
# Note

# title: 'Comparing ANOVA to Mixed Effects Linear Regression for analysis of
#         placebo-controlled repeat measures'

# The purpose of this Python script is to demonstrate a comparison between
# the methods, analysis of variance (ANOVA) with mixed effects and repeated
# measures against linear regression with mixed effects for the analysis of
# data from a placebo-controlled clinical trial with repeated measures before
# and after intervention.

# References:

# 1. Scientific literature article with peer review comparing analysis of
#     variance (ANOVA) with repeat measures to linear regression with mixed
#     effects to analyze the effect of treatment between groups that receive
#     intervention or placebo control
#  PubMed: 15081686
#  title: 'Mixed-model regression analysis and dealing with interindividual
#          differences'
#  journal: 'Methods in Enzymology'
#  year: 2004
#  authors: Hans P.A. Van Dongen; ... Greg Maislin
#  note:
#  - The study design included measurement at baseline before treatment and at
#    3 time points after treatment by either placebo or intervention.
#  - TCW transferred the raw data from Table 1 to file "table_data.tsv" within
#    the directory 'partner/repository/demonstration/15081686_vandongen_2004'.
#  - TCW encoded the placebo group as value '0' for feature 'condition'.
#  - TCW encoded the intervention group as value '1' for feature 'condition'.
#  - TCW endoced the 'day 0', 'day 1', 'day 2', 'day 3' time points of
#    measurement as values '0', '1', '2', and '3' respectively for feature
#    'time_point'.
#  - In the main analysis, TCW only considered values 'day 0' and 'day 1' of
#    'time_point'.

################################################################################
# Installation and importation

# Standard
import sys
import os
import copy
import textwrap

# Relevant
import pandas
import scipy
import numpy
#import statsmodels.api
#import pingouin

# Custom
import partner.utility as putly
import partner.parallelization as prall
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.regression as preg

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Read source information from file.


def define_type_table_columns_vandongen():
    """
    Defines the types of variables for columns in table.

    Review: TCW; 7 April 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    types_columns["identifier_subject"] = "string"
    types_columns["condition"] = "int32"
    types_columns["time_point"] = "int32"
    types_columns["measurement"] = "float32"
    # Return information.
    return types_columns


def read_source_table_data_vandongen(
    name_file_table_data=None,
    path_directory_source=None,
    report=None,
):
    """
    Read and organize source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 7 April 2025

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of data with features and
            observations for regression

    """

    # Define paths to child files.
    path_file_table_data = os.path.join(
        path_directory_source, name_file_table_data,
    )

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_type_table_columns_vandongen()
    table = pandas.read_csv(
        path_file_table_data,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    table["condition"] = pandas.to_numeric(
        table["condition"],
        downcast="integer",
        errors="coerce",
    )
    table["time_point"] = pandas.to_numeric(
        table["time_point"],
        downcast="integer",
        errors="coerce",
    )
    table["measurement"] = pandas.to_numeric(
        table["measurement"],
        downcast="float",
        errors="coerce",
    )

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "read_source_table_data_vandongen()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("table:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table



def define_type_table_columns_waller():
    """
    Defines the types of variables for columns in table.

    Review: TCW; 7 April 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()

    types_columns["tissue"] = "string"
    types_columns["subject_visit"] = "string"
    types_columns["identifier_subject"] = "string"
    types_columns["sex_text"] = "string"
    types_columns["sex_y"] = "float32"
    types_columns["age_cohort_text"] = "string"
    types_columns["age_cohort"] = "float32"
    types_columns["intervention_text"] = "string"
    types_columns["intervention_placebo"] = "float32"
    types_columns["intervention_omega3"] = "float32"
    types_columns["intervention_placebo_other"] = "float32"
    types_columns["intervention_omega3_other"] = "float32"
    types_columns["intervention_after_placebo"] = "float32"
    types_columns["intervention_after_omega3"] = "float32"
    types_columns["visit_text"] = "string"
    types_columns["visit_second"] = "float32"

    types_columns["age"] = "float32"
    types_columns["body_mass_index"] = "float32"

    types_columns["triglyceride"] = "float32"
    types_columns["omega3_eicosapentaenoate"] = "float32"
    types_columns["omega3_docosahexaenoate"] = "float32"
    types_columns["c_react_protein"] = "float32"
    types_columns["adipocyte_diameter"] = "float32"
    types_columns["adipocyte_lipid_content"] = "float32"
    types_columns["insulin_sensitivity"] = "float32"
    types_columns["homa_insulin_resist"] = "float32"
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


def read_source_table_data_waller(
    name_file_table_data=None,
    path_directory_source=None,
    report=None,
):
    """
    Read and organize source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 7 April 2025

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of data with features and
            observations for regression

    """

    # Define paths to child files.
    path_file_table_data = os.path.join(
        path_directory_source, name_file_table_data,
    )

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_type_table_columns_waller()
    table = pandas.read_csv(
        path_file_table_data,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    #table["condition"] = pandas.to_numeric(
    #    table["condition"],
    #    downcast="integer",
    #    errors="coerce",
    #)

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "read_source_table_data_waller()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("table:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


################################################################################
# Procedure


##########
# Manage procedure for demonstration data from Van Dongen et al.


def manage_procedure_vandongen(
    name_file_table_data=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Manage a major branch in the module's main behavior.

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "manage_procedure_vandongen()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("name_file_table_data: " + str(name_file_table_data))
        print("path_directory_source: " + str(path_directory_source))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    table_source = read_source_table_data_vandongen(
        name_file_table_data=name_file_table_data,
        path_directory_source=path_directory_source,
        report=report,
    )

    ##########
    # Organize information in table.
    table_source["observations"] = table_source["identifier_subject"]
    table_source["identifier_subject_raw"] = table_source["identifier_subject"]
    table_source["identifier_subject"] = table_source.apply(
        lambda row: str(
            "subject_" + row["identifier_subject_raw"]
        ),
        axis="columns", # apply function to each row
    )
    selection_observations = dict()
    selection_observations["condition"] = [0, 1,]
    selection_observations["time_point"] = [0, 1,]
    features_relevant = [
        "observations",
        "identifier_subject",
        "condition",
        "time_point",
        "time_point_0",
        "time_point_1",
        "condition_by_time_point_1",
        "caffeine_at_time_point_1",
        "placebo_at_time_point_1",
        "measurement",
    ]
    table = preg.organize_table_data(
        table=table_source,
        selection_observations=selection_observations,
        features_relevant=features_relevant,
        features_regression=features_relevant,
        features_continuity_scale=list(),
        index_columns_source="features",
        index_columns_product="features",
        index_rows_source="observations",
        index_rows_product="observations",
        adjust_scale=False,
        method_scale=None, # 'z_score' or 'unit_range'
        explicate_indices=True,
        report=report,
    )

    # Two-Way ANOVA with Repeated Measures.
    # This version of ANOVA does not fit the study design and data, because the
    # predictor feature for 'condition' distinguishes between subjects and not
    # within subjects.
    #pail_anova = preg.determine_type_analysis_variance_anova(
    #    table=table,
    #    index_columns="features",
    #    index_rows="observations",
    #    type_anova="2way_repeat",
    #    formula_text="(measurement) ~ (condition) + (time_point)",
    #    feature_response="measurement",
    #    features_predictor_between=["condition",],
    #    features_predictor_within=["time_point_1"],
    #    subject="identifier_subject",
    #    report=report,
    #)

    # Two-Way ANOVA with Repeated Measures and Mixed Effects.
    description_analysis = str(
        "Two-Way Analysis of Variance (ANOVA) with Mixed Effects and " +
        "Repeated Measures"
    )
    formula_text = str(
        "(measurement) ~ (condition) + (time_point_1) + " +
        "(condition_by_time_point_1) + (subject)"
    )
    description_response = str(
        "Score on psychomotor vigilance task (PVT)"
    )
    description_groups_random = str(
        "Identifiers of individual subjects, people in study"
    )
    description_predictor = textwrap.dedent("""\
        effects between subjects:
           condition (0: placebo; 1: caffeine)
        effects within subjects:
           time_point_1 (0: before intervention; 1: after 1 day of sleep
           deprivation)
    """)
    pail_anova_mix = preg.determine_type_analysis_variance_anova(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_anova="2way_repeat_mixed",
        formula_text=formula_text,
        feature_response="measurement",
        features_predictor_between=["condition",],
        features_predictor_within=["time_point_1"],
        subject="identifier_subject",
        report=report,
    )
    summary_text_anova = preg.prepare_text_summary_regression_anova(
        title="Van Dongen et al, Methods Enzymology, 2024, PubMed:15081686",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_anova_mix["anova"].round(3)),
        summary_2=str(pail_anova_mix["t_test_within"]),
        summary_3=str(pail_anova_mix["t_test_between"]),
        report=report,
    )

    # Linear Regression with Mixed Effects.
    # This is the standard design of analysis for placebo-controlled clinical
    # trial with paired, repeated measures in subjects before and after
    # treatment.
    description_analysis = str(
        "Linear Regression with Mixed Effects: Fixed Slopes and Random " +
        " Slopes and Intercepts; Standard Model"
    )
    formula_text = str(
        "(measurement) ~ (condition) + (time_point_1) + " +
        "(condition_by_time_point_1) + (subject)"
    )
    description_response = str(
        "Score on psychomotor vigilance task (PVT)"
    )
    description_groups_random = str(
        "Identifiers of individual subjects, people in study"
    )
    description_predictor = textwrap.dedent("""\
        fixed effects:
           condition (0: placebo; 1: caffeine);
           time_point_1 (0: before intervention; 1: after 1 day of sleep
           deprivation);
           condition * time_point_1
        random effects, intercepts:
           for each subject
        random effects, slope coefficients:
           time_point_1
           - Allow each subject to respond differently to sleep deprivation,
             regardless of placebo or caffeine.
           - In this study, the sleep deprivation was an additional
             intervention.
    """)
    pail_regression = preg.determine_type_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_regression="continuous_mixed",
        formula_text=formula_text,
        feature_response="measurement",
        features_predictor_fixed=[
            "condition", "time_point_1", "condition_by_time_point_1",
        ],
        features_predictor_random=["time_point_1",], # list()
        groups_random="identifier_subject",
        method_fit=None, # use default method to fit model
        report=report,
    )
    summary_text_regression_1 = preg.prepare_text_summary_regression_anova(
        title="Van Dongen et al, Methods Enzymology, 2024, PubMed:15081686",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_regression["table_summary"]),
        summary_2="",
        summary_3="",
        report=report,
    )

    # Linear Regression with Mixed Effects.
    # This is an alternative design of analysis for placebo-controlled clinical
    # trial with paired, repeated measures in subjects before and after
    # treatment.
    # group 1: after real treatment intervention
    # group 2: after placebo mock treatment intervention
    # group 3: before real or mock treatment intervention
    # Represent these three experimental groups as categories using two (n - 1)
    # dummy indicator variables.
    # The dummy indicator variable for "group 1" (real treatment) is
    # effectively the same as the interaction variable in the standard design
    # above, "condition_by_time_point_1".
    # While the basic features "condition" and "time_point_1" can be included
    # in the model, these features are themselves less relevant to the
    # hypothesis.
    # Including the dummy indicator variable for "group 2" (placebo mock
    # treatment) completes the representation of all three experimental groups
    # and also asks directly about the effect of placebo mock treatment.
    description_analysis = str(
        "Linear Regression with Mixed Effects: Fixed Slopes and Random " +
        " Slopes and Intercepts; Alternative Model"
    )
    formula_text = str(
        "(measurement) ~ (caffeine_at_time_point_1) + " +
        "(placebo_at_time_point_1) + (subject)"
    )
    description_response = str(
        "Score on psychomotor vigilance task (PVT)"
    )
    description_groups_random = str(
        "Identifiers of individual subjects, people in study"
    )
    description_predictor = textwrap.dedent("""\
        fixed effects:
           caffeine_at_time_point_1 (0: other; 1: caffeine at time_point 1);
           placebo_at_time_point_1 (0: other; 1: placebo at time_point 1);
        random effects, intercepts:
           for each subject
        random effects, slope coefficients:
           caffeine_at_time_point_1
           - Allow each subject to respond differently to sleep deprivation
             with caffeine.
    """)
    pail_regression = preg.determine_type_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_regression="continuous_mixed",
        formula_text=formula_text,
        feature_response="measurement",
        features_predictor_fixed=[
            "caffeine_at_time_point_1", "placebo_at_time_point_1",
        ],
        features_predictor_random=["caffeine_at_time_point_1",], # list()
        groups_random="identifier_subject",
        method_fit=None, # use default method to fit model
        report=report,
    )
    summary_text_regression_2 = preg.prepare_text_summary_regression_anova(
        title="Van Dongen et al, Methods Enzymology, 2024, PubMed:15081686",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_regression["table_summary"]),
        summary_2="",
        summary_3="",
        report=report,
    )

    ##########
    # Write information to file.
    # Write for each parallel instance of regression.
    # A subsequent procedure will read the information from file and collect it
    # within a summary for all instances of regression.

    # Collect information.
    # Collections of files.
    pail_write_text = dict()
    pail_write_text[str("vandongen_anova")] = summary_text_anova
    pail_write_text[str("vandongen_regression_1")] = summary_text_regression_1
    pail_write_text[str("vandongen_regression_2")] = summary_text_regression_2
    # Write product information to file.
    putly.write_character_strings_to_file_text(
        pail_write=pail_write_text,
        path_directory=path_directory_product,
    )

    pass


##########
# Manage procedure for demonstration data from Waller et al.


def manage_procedure_waller(
    name_file_table_data=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Manage a major branch in the module's main behavior.

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "manage_procedure_waller()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("name_file_table_data: " + str(name_file_table_data))
        print("path_directory_source: " + str(path_directory_source))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    table_source = read_source_table_data_waller(
        name_file_table_data=name_file_table_data,
        path_directory_source=path_directory_source,
        report=report,
    )

    #intervention_text_placebo_other_after_only:placebo,other;intervention_text_active_other_after_only:active,other

    ##########
    # Organize information in table.
    #table_source["observations"] = table_source["subject_visit"]
    selection_observations = dict()
    selection_observations["tissue"] = ["adipose",]
    selection_observations["age_cohort_text"] = ["elder",]
    selection_observations["sex_text"] = ["female", "male",]
    selection_observations["visit_text"] = ["first", "second",]
    selection_observations["intervention_text"] = ["placebo", "omega3",]
    #features_relevant = [
    #    "observations",
    #    "identifier_subject",
    #    "condition",
    #    "time_point",
    #    "measurement",
    #]
    types_columns = define_type_table_columns_waller()
    features_relevant = list(types_columns.keys())
    #features_relevant.insert(0, "observations",)
    features_continuity_scale = [
        #"age",
        #"body_mass_index",
        "triglyceride",
        "omega3_eicosapentaenoate",
        "omega3_docosahexaenoate",
    ]
    table = preg.organize_table_data(
        table=table_source,
        selection_observations=selection_observations,
        features_relevant=features_relevant,
        features_regression=features_relevant,
        features_continuity_scale=features_continuity_scale,
        index_columns_source="features",
        index_columns_product="features",
        index_rows_source="subject_visit",
        index_rows_product="observations",
        adjust_scale=False,
        method_scale=None, # 'z_score' or 'unit_range'
        explicate_indices=True,
        report=report,
    )

    # Two-Way ANOVA with Repeated Measures.
    # This version of ANOVA does not fit the study design and data, because the
    # predictor feature for 'condition' distinguishes between subjects and not
    # within subjects.
    #pail_anova = preg.determine_type_analysis_variance_anova(
    #    table=table,
    #    index_columns="features",
    #    index_rows="observations",
    #    type_anova="2way_repeat",
    #    formula_text="(measurement) ~ (condition) + (time_point)",
    #    feature_response="measurement",
    #    features_predictor_between=["condition",],
    #    features_predictor_within=["time_point_1"],
    #    subject="identifier_subject",
    #    report=report,
    #)

    # Two-Way ANOVA with Repeated Measures and Mixed Effects.
    description_analysis = str(
        "Two-Way Analysis of Variance (ANOVA) with Mixed Effects and " +
        "Repeated Measures"
    )
    formula_text = str(
        "(triglyceride) ~ (intervention_omega3) + (visit_second) + " +
        "(intervention_after_omega3) + (identifier_subject)"
    )
    description_response = str(
        "Levels of triglycerides in blood (mg/dL)"
    )
    description_groups_random = str(
        "Identifiers of individual subjects, people in study"
    )
    description_predictor = textwrap.dedent("""\
        effects between subjects:
           intervention_omega3 (0: placebo; 1: omega3)
        effects within subjects:
           visit_second (0: before intervention; 1: after intervention)
    """)
    pail_anova_mix = preg.determine_type_analysis_variance_anova(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_anova="2way_repeat_mixed",
        formula_text=formula_text,
        feature_response="triglyceride",
        features_predictor_between=["intervention_omega3",],
        features_predictor_within=["visit_second"],
        subject="identifier_subject",
        report=report,
    )
    summary_text_anova = preg.prepare_text_summary_regression_anova(
        title="Waller et al, clinical study on omega-3 fatty acids",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_anova_mix["anova"].round(3)),
        summary_2=str(pail_anova_mix["t_test_within"]),
        summary_3=str(pail_anova_mix["t_test_between"]),
        report=report,
    )

    # Linear Regression with Mixed Effects.
    # This is the standard design of analysis for placebo-controlled clinical
    # trial with paired, repeated measures in subjects before and after
    # treatment.
    description_analysis = str(
        "Linear Regression with Mixed Effects: Fixed Slopes and Random " +
        " Slopes and Intercepts; Standard Model"
    )
    formula_text = str(
        "(triglyceride) ~ (intervention_omega3) + (visit_second) + " +
        "(intervention_after_omega3) + (identifier_subject)"
    )
    description_response = str(
        "Levels of triglycerides in blood (mg/dL)"
    )
    description_groups_random = str(
        "Identifiers of individual subjects, people in study"
    )
    description_predictor = textwrap.dedent("""\
        fixed effects:
           intervention_omega3 (0: placebo; 1: omega3);
           visit_second (0: before intervention; 1: after intervention);
           intervention_after_omega3 (0: placebo or before; 1: after w3);
        random effects, intercepts:
           for each subject
        random effects, slope coefficients:
           none
           - Here, do not allow random slopes to simplify.
           - Random slopes for visit_second would allow each subject to respond
             differently to the duration of the study, independently of whether
             the intervention was placebo or omega3.
    """)
    pail_regression = preg.determine_type_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_regression="continuous_mixed",
        formula_text=formula_text,
        feature_response="triglyceride",
        features_predictor_fixed=[
            "intervention_omega3", "visit_second",
            "intervention_after_omega3",
        ],
        features_predictor_random=list(), # ["visit_second",]
        groups_random="identifier_subject",
        method_fit=None, # use default method to fit model
        report=report,
    )
    summary_text_regression_1 = preg.prepare_text_summary_regression_anova(
        title="Waller et al, clinical study on omega-3 fatty acids",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_regression["table_summary"]),
        summary_2="",
        summary_3="",
        report=report,
    )

    # Linear Regression with Mixed Effects.
    # This is an alternative design of analysis for placebo-controlled clinical
    # trial with paired, repeated measures in subjects before and after
    # treatment.
    # group 1: after real treatment intervention
    # group 2: after placebo mock treatment intervention
    # group 3: before real or mock treatment intervention
    # Represent these three experimental groups as categories using two (n - 1)
    # dummy indicator variables.
    # The dummy indicator variable for "group 1" (real treatment) is
    # effectively the same as the interaction variable in the standard design
    # above, "condition_by_time_point_1".
    # While the basic features "condition" and "time_point_1" can be included
    # in the model, these features are themselves less relevant to the
    # hypothesis.
    # Including the dummy indicator variable for "group 2" (placebo mock
    # treatment) completes the representation of all three experimental groups
    # and also asks directly about the effect of placebo mock treatment.
    description_analysis = str(
        "Linear Regression with Mixed Effects: Fixed Slopes and Random " +
        " Slopes and Intercepts; Alternative Model"
    )
    formula_text = str(
        "(triglyceride) ~ (intervention_after_omega3) + " +
        "(intervention_after_placebo) + (identifier_subject)"
    )
    description_response = str(
        "Levels of triglycerides in blood (mg/dL)"
    )
    description_groups_random = str(
        "Identifiers of individual subjects, people in study"
    )
    description_predictor = textwrap.dedent("""\
        fixed effects:
           intervention_after_omega3 (0: placebo or before; 1: after w3);
           intervention_after_placebo (0: omega3 or before; 1: after placebo);
        random effects, intercepts:
           for each subject
        random effects, slope coefficients:
           none
           - Here, do not allow random slopes, to simplify.
    """)
    pail_regression = preg.determine_type_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_regression="continuous_mixed",
        formula_text=formula_text,
        feature_response="triglyceride",
        features_predictor_fixed=[
            "intervention_after_omega3", "intervention_after_placebo",
        ],
        features_predictor_random=list(),
        groups_random="identifier_subject",
        method_fit=None,
        report=report,
    )
    summary_text_regression_2 = preg.prepare_text_summary_regression_anova(
        title="Waller et al, clinical study on omega-3 fatty acids",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_regression["table_summary"]),
        summary_2="",
        summary_3="",
        report=report,
    )

    ##########
    # Write information to file.
    # Write for each parallel instance of regression.
    # A subsequent procedure will read the information from file and collect it
    # within a summary for all instances of regression.

    # Collect information.
    # Collections of files.
    pail_write_text = dict()
    pail_write_text[str("waller_anova")] = summary_text_anova
    pail_write_text[str("waller_regression_1")] = summary_text_regression_1
    pail_write_text[str("waller_regression_2")] = summary_text_regression_2
    # Write product information to file.
    putly.write_character_strings_to_file_text(
        pail_write=pail_write_text,
        path_directory=path_directory_product,
    )

    pass


##########
# Execute main procedure.


def execute_procedure(
    name_file_table_data_vandongen=None,
    path_directory_source_vandongen=None,
    path_directory_product_vandongen=None,
    name_file_table_data_waller=None,
    path_directory_source_waller=None,
    path_directory_product_waller=None,
    path_directory_dock=None,
    report=None,
):
    """
    Execute module's main behavior.

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parameters.
    if (
        (report is not None) and
        (str(report) != "") and
        (str(report) == "true")
    ):
        report = True
    else:
        report = False
        pass

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "execute_procedure()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("system: local")
        #print("name_file_table_data: " + str(name_file_table_data))
        #print("path_directory_source: " + str(path_directory_source))
        #print("path_directory_product: " + str(path_directory_product))
        #print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass


    ##########
    # Manage procedure for data from Van Dongen et al.
    if True:
        manage_procedure_vandongen(
            name_file_table_data=name_file_table_data_vandongen,
            path_directory_source=path_directory_source_vandongen,
            path_directory_product=path_directory_product_vandongen,
            path_directory_dock=path_directory_dock,
            report=report,
        )

    ##########
    # Manage procedure for data from Waller et al.
    if True:
        manage_procedure_waller(
            name_file_table_data=name_file_table_data_waller,
            path_directory_source=path_directory_source_waller,
            path_directory_product=path_directory_product_waller,
            path_directory_dock=path_directory_dock,
            report=report,
        )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    name_file_table_data_vandongen = sys.argv[1]
    path_directory_source_vandongen = sys.argv[2]
    path_directory_product_vandongen = sys.argv[3]
    name_file_table_data_waller = sys.argv[4]
    path_directory_source_waller = sys.argv[5]
    path_directory_product_waller = sys.argv[6]
    path_directory_dock = sys.argv[7]
    report = sys.argv[8]

    # Call function for procedure.
    execute_procedure(
        name_file_table_data_vandongen=name_file_table_data_vandongen,
        path_directory_source_vandongen=path_directory_source_vandongen,
        path_directory_product_vandongen=path_directory_product_vandongen,
        name_file_table_data_waller=name_file_table_data_waller,
        path_directory_source_waller=path_directory_source_waller,
        path_directory_product_waller=path_directory_product_waller,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass


#
