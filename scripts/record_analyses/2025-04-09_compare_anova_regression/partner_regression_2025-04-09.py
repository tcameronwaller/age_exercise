"""
Supply functionality for regression analysis.

This module 'regression' is part of the 'partner' package.

This module is not directly executable.

This subpackage 'partner' provides executable functionality under the
management of a higher level package. Importation paths require this hierarchy.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This module file is part of the project package directory 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
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

###############################################################################
# Notes

# TODO: TCW; 26 June 2024
# Support more regression models
# Ordinary Least Squares
#   - https://www.statsmodels.org/dev/examples/notebooks/generated/ols.html
# Generalized Linear Models
#   - https://www.statsmodels.org/stable/glm.html
# Mixed effects models
# - support pairs of samples from the same subject
# ANOVA
# ANCOVA


###############################################################################
# Installation and importation

# Standard
import os
import csv
import copy
import textwrap
import string
import gzip
import shutil
import textwrap
import itertools
import math

# Relevant

import pandas
import scipy
import numpy
import statsmodels.api
import statsmodels.stats.outliers_influence
import pingouin

# Custom
import partner.utility as putly
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Prepare table of parameters for many different response features, all with
# identical predictor features and other parameters.

# Refer to Bash script "template_create_parameter_table.sh".



##########
# Organize information in table for features and observations.


def organize_table_data(
    table=None,
    selection_observations=None,
    features_relevant=None,
    features_regression=None,
    features_continuity_scale=None,
    index_columns_source=None,
    index_columns_product=None,
    index_rows_source=None,
    index_rows_product=None,
    adjust_scale=None,
    method_scale=None,
    explicate_indices=None,
    report=None,
):
    """
    Organize table of data with features and observations for regression.

    1. filter columns in table for relevant features
    2. filter rows in table for relevant observations
    3. standardize scale of features with measurement values on quantitative,
       continuous scale of measurement, ratio or interval

    ----------
    Format of source data table (name: "table")
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
    Format of product data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        selection_observations (dict<list<str>>): names of columns in data
            table for feature variables and their categorical values by
            which to filter rows for observations in data table
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        features_regression (list<str>): names of columns in data table for
            feature variables that are relevant to the actual regression model
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        index_columns_source (str): name of single-level index across columns
            in table
        index_columns_product (str): name of single-level index across columns
            in table
        index_rows_source (str): name of single-level index across rows in
            table
        index_rows_product (str): name of single-level index across rows in
            table
        adjust_scale (bool): whether to adjust or standardize the scale of
            values for features across observations
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        explicate_indices (bool): whether to explicate, define, or specify
            explicit indices across columns and rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    selection_observations = copy.deepcopy(selection_observations)
    features_relevant = copy.deepcopy(features_relevant)
    features_regression = copy.deepcopy(features_regression)
    features_continuity_scale = copy.deepcopy(features_continuity_scale)

    # Restore or reset indices to generic default.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index

    # Filter rows in table for observations.
    table = porg.filter_select_table_rows_by_columns_categories(
        table=table,
        columns_categories=selection_observations,
        report=report,
    )
    # Filter columns in table for features.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=features_relevant,
        report=report,
    )

    # Remove rows from table for observations with any missing values for
    # relevant features that are part of the model for regression.
    table.dropna(
        axis="index",
        how="any",
        subset=features_regression,
        inplace=True,
    )

    # Standardize scale of values for observations of features.
    if (
        (adjust_scale) and
        (len(features_continuity_scale) > 0) and
        (method_scale is not None) and
        (str(method_scale).strip().lower() != "none")
    ):
        if (
            (str(method_scale).strip().lower() == "z_score")
        ):
            table = pscl.transform_standard_z_score_by_table_columns(
                    table=table,
                    columns=features_continuity_scale,
                    report=report,
            )
        elif (
            (str(method_scale).strip().lower() == "unit_range")
        ):
            table = pscl.transform_unit_range_by_table_columns(
                    table=table,
                    columns=features_continuity_scale,
                    report=report,
            )
            pass
        pass

    # Organize indices in table.
    # Determine whether parameters specified a column in the table for
    # identifiers of observations that will become index across rows.
    if (
        (index_rows_source is not None) and
        (index_rows_source in (table.columns.tolist()))
    ):
        # Create index across rows from column that already exists in table.
        table = porg.explicate_table_indices_columns_rows_single_level(
            table=table,
            index_columns=index_columns_source,
            index_rows=index_rows_source,
            explicate_indices=explicate_indices,
            report=False,
        )
        # Standardize names of indices.
        table = porg.translate_names_table_indices_columns_rows(
            table=table,
            index_columns_product=index_columns_product,
            index_rows_source=index_rows_source,
            index_rows_product=index_rows_product,
            report=None,
        )
    else:
        # Name generic index in table.
        table.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table.index.set_names(index_rows_product, inplace=True)
        # Standardize names of indices.
        table = porg.translate_names_table_indices_columns_rows(
            table=table,
            index_columns_product=index_columns_product,
            index_rows_source=index_rows_product,
            index_rows_product=index_rows_product,
            report=None,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "organize_table_data()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


# Evaluate information in table for features and observations.


def evaluate_table_data(
    table=None,
    selection_observations=None,
    features_relevant=None,
    features_regression=None,
    features_continuity_scale=None,
    index_columns_source=None,
    index_columns_product=None,
    index_rows_source=None,
    index_rows_product=None,
    adjust_scale=None,
    method_scale=None,
    explicate_indices=None,
    report=None,
):
    """
    Evaluate table of data with features and observations for regression.

    This function is useful to evaluate the variance of values for features
    across observations after adjusting their scale to the common unit range.
    This evaluation enables a manual, supervised approach to a method of
    feature selection that is called 'variance thresholding'. The rationale is
    that features with very little variance are unlikely to contribute to the
    regression model, at least not with a favorable ratio of signal to noise.

    From basic testing, it seems that a random number with a uniform
    distribution between zero (0) and one (1) has a variance of about 10% or
    0.1. After adjusting the scales of features to have unit range, a
    reasonable threshold on the variance might be approximately 1% or 0.01.

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        selection_observations (dict<list<str>>): names of columns in data
            table for feature variables and their categorical values by
            which to filter rows for observations in data table
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        features_regression (list<str>): names of columns in data table for
            feature variables that are relevant to the actual regression model
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        index_columns_source (str): name of single-level index across columns
            in table
        index_columns_product (str): name of single-level index across columns
            in table
        index_rows_source (str): name of single-level index across rows in
            table
        index_rows_product (str): name of single-level index across rows in
            table
        adjust_scale (bool): whether to adjust or standardize the scale of
            values for features across observations
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        explicate_indices (bool): whether to explicate, define, or specify
            explicit indices across columns and rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Organize information in table.
    # For all features in the regression model (including the response), adjust
    # the scale of values to have common unit range.
    features_scale = copy.deepcopy(features_regression)
    features_category = copy.deepcopy(list(selection_observations.keys()))
    features_category.insert(0, index_rows_source)
    features_scale = list(filter(
        lambda feature: (feature not in features_category),
        features_scale
    ))
    table = preg.organize_table_data(
        table=table,
        selection_observations=selection_observations,
        features_relevant=features_relevant,
        features_regression=features_regression,
        features_continuity_scale=features_scale,
        index_columns_source=index_columns_source,
        index_columns_product=index_columns_product,
        index_rows_source=index_rows_source,
        index_rows_product=index_rows_product,
        adjust_scale=True,
        method_scale="unit_range", # "z_score" or "unit_range"
        explicate_indices=True,
        report=False,
    )
    # Filter columns in table for features.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=features_scale,
        report=report,
    )

    # Calculate descriptive statistical measures for features.
    table_mean = table.aggregate(
        lambda series: series.mean(),
        axis="index", # apply function to each column
    )
    table_variance = table.aggregate(
        lambda series: series.var(),
        axis="index", # apply function to each column
    )
    table_deviation = table.aggregate(
        lambda series: series.std(),
        axis="index", # apply function to each column
    )
    table_minimum = table.aggregate(
        lambda series: series.min(),
        axis="index", # apply function to each column
    )
    table_maximum = table.aggregate(
        lambda series: series.max(),
        axis="index", # apply function to each column
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "evaluate_table_data()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        print("table before aggregation:")
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=4)
        print("series after aggregation:")
        putly.print_terminal_partition(level=5)
        print("mean")
        putly.print_terminal_partition(level=6)
        print(table_mean.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("variance")
        putly.print_terminal_partition(level=6)
        print(table_variance.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("standard deviation")
        putly.print_terminal_partition(level=6)
        print(table_deviation.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("minimum")
        putly.print_terminal_partition(level=6)
        print(table_minimum.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("maximum")
        putly.print_terminal_partition(level=6)
        print(table_maximum.iloc[0:10])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


# Check parameters and table of data for performing regression analysis.


def check_parameters_table_data_regression(
    table=None,
    name_instance=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    groups_random=None,
    threshold_features_variance=None,
    threshold_observations_count=None,
    measure_variance=None,
    report=None,
):
    """
    Check parameters for regression, including the table of data with features
    and observations.

    For the check on the variance of features across observations, it is most
    accurate and informative to compare the variances after transforming the
    values of features to have unit range between zero (0) and one (1). Do not
    compare variances after transforming the values of features to standard
    z scores, as this transformation forces values of all features to have the
    same variance, unit standard deviation.

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        name_instance (str): compound name for instance of parameters
        type_regression (str): name of type of regression model to use
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        groups_random (str): name of column in data table for identifiers or
            designations of groups of observations for which to allow random
            effects in the regression model
        threshold_features_variance (float): threshold minimal variance for
            values of features across observations
        threshold_observations_count (int): threshold minimal count of
            observations that must have nonmissing values of all features
        measure_variance (str): name of measure to use for variance, either
            'variance', 'standard_deviation', or 'coefficient_variation'
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_fixed = copy.deepcopy(features_predictor_fixed)
    features_predictor_random = copy.deepcopy(features_predictor_random)

    # Check that functionality currently supports the type of regression.
    types = [
        "continuous_ols",
        "continuous_mixed",
        "discrete_logit",
        #"linear_r",
        #"logistic_r",
    ]
    if (
        (type_regression is not None) and
        (len(str(type_regression)) > 0) and
        (str(type_regression).strip().lower() != "none") and
        (str(type_regression) in types)
    ):
        check_type = True
    else:
        check_type = False
        pass

    # Check for features in the regression model that features exist in the
    # table of data.
    columns_available = copy.deepcopy(table.columns.to_list())
    check_fixed = putly.compare_lists_by_inclusion(
        items_dominant=columns_available,
        items_subordinate=features_predictor_fixed,
    )
    check_random = putly.compare_lists_by_inclusion(
        items_dominant=columns_available,
        items_subordinate=features_predictor_random,
    )
    if (
        (len(columns_available) > 1) and
        (feature_response in columns_available) and
        (check_fixed) and
        (check_random)
    ):
        check_features_exist = True
    else:
        check_features_exist = False
        pass

    # Check for features in the regression model that adequate observations
    # with nonmissing values exist in the table of data.
    def check_features_nonmissing_observations(
        table=None,
        features=None,
        threshold_observations_count=None,
    ):
        """
        Define subordinate function for internal use within dominant function.
        """
        for feature in features:
            count = table[feature].dropna().count()
            if (count >= threshold_observations_count):
                return True
            else:
                return False
            pass
        pass
    # Copy information.
    table_count = table.copy(deep=True)
    # Evaluate features separately.
    #check_fixed = check_features_nonmissing_observations(
    #    table=table_count,
    #    features=features_predictor_fixed,
    #    threshold_observations_count=threshold_observations_count,
    #)
    # Evaluate features together.
    features_regression = list()
    features_regression.append(feature_response)
    features_regression.extend(features_predictor_fixed)
    features_regression.extend(features_predictor_random)
    table_count.dropna(
        axis="index",
        how="any",
        subset=features_regression,
        inplace=True,
    )
    count_rows = int(table_count.shape[0])
    if (count_rows >= threshold_observations_count):
        check_observations_count = True
    else:
        check_observations_count = False
        pass

    # Check for features in the regression model that values across
    # observations of each feature have adequate variance.
    def check_features_variance(
        table=None,
        features=None,
        measure_variance=None,
        threshold_features_variance=None,
    ):
        """
        Define subordinate function for internal use within dominant function.
        """
        checks = list()
        for feature in features:
            values_raw = table[feature].to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
            pail = pdesc.calculate_variance_measures(
                array=values_raw,
            )
            if (pail[measure_variance] >= threshold_features_variance):
                checks.append(True)
            else:
                checks.append(False)
                pass
            pass
        return all(checks)
    # Evaluate individual features separately.
    check_fixed = check_features_variance(
        table=table,
        features=features_predictor_fixed,
        measure_variance=measure_variance,
        threshold_features_variance=threshold_features_variance,
    )
    # Features with random effects will generall be categorical.
    if False:
        check_random = check_features_variance(
            table=table,
            features=features_predictor_random,
            measure_variance=measure_variance,
            threshold_features_variance=threshold_features_variance,
        )
    else:
        check_random = True
        pass
    if (check_fixed and check_random):
        check_variance = True
    else:
        check_variance = False
        pass

    # Check whether there was failure of any checks overall.
    if (
        (check_type) and
        (check_features_exist) and
        (check_observations_count) and
        (check_variance)
    ):
        check_overall = True
    else:
        check_overall = False
        pass

    # Collect information.
    pail = dict()
    pail["check_type"] = check_type
    pail["check_features_exist"] = check_features_exist
    pail["check_observations_count"] = check_observations_count
    pail["check_variance"] = check_variance
    pail["check_overall"] = check_overall

    # Warn.
    if (not check_overall):
        # Print.
        putly.print_terminal_partition(level=3)
        putly.print_terminal_warning()
        putly.print_terminal_partition(level=5)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "check_parameters_table_data_regression()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(
            "Parameters and or data for a regression analysis failed checks!"
        )
        print("name_instance: " + name_instance)
        putly.print_terminal_partition(level=5)

        pass

    # Report.
    if report:
        # Organize.
        #count_records = len(records)
        # Print.
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "check_parameters_table_data_regression()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("check_type: " + str(check_type))
        print("check_features_exist: " + str(check_features_exist))
        print("check_observations_count: " + str(check_observations_count))
        print("check_variance: " + str(check_variance))
        print("check_overall: " + str(check_overall))
        pass
    # Return information.
    return pail


##########
# Organize information from regresion.


def define_sequence_features(
    features=None,
):
    """
    Define specific sequence of features.

    Review: TCW; 2 April 2025

    arguments:
        features (list<str>): names of features for which to define specific
            sequence

    raises:

    returns:
        (dict<str>): collection of information

    """

    # Copy information.
    features = copy.deepcopy(features)

    # Option 1.
    if True:
        sequence_features = dict()
        index = 0
        for name in features:
            sequence_features[index] = name
            index += 1
            pass
    # Option 2.
    if False:
        sequence_features = dict(zip(
            range(len(features)),
            features
        ))
    # Option 3.
    if False:
        sequence_features = {
            key: value for key, value in enumerate(features)
        }

    # Return information.
    return sequence_features


def define_features_sequence(
    features=None,
):
    """
    Define specific sequence of features.

    Review: TCW; 2 April 2025

    arguments:
        features (list<str>): names of features for which to define specific
            sequence

    raises:

    returns:
        (dict<int>): collection of information

    """

    # Copy information.
    features = copy.deepcopy(features)

    # Option 1.
    if True:
        features_sequence = dict()
        index = 0
        for name in features:
            features_sequence[name] = index
            index += 1
            pass
    # Option 2.
    if False:
        features_sequence = dict(zip(
            features,
            range(len(features))
        ))
    # Option 3.
    if False:
        features_sequence = {
            key: value for value, key in enumerate(features)
        }

    # Return information.
    return features_sequence


def parse_regression_model_predictor(
    table_predictor_intercept=None,
    indicator=None,
    name_source=None,
    name_product=None,
    parameters=None,
    standard_errors=None,
    pvalues=None,
    variance_inflation=None,
):
    """
    Parse raw information about predictor independent variables from a model
    of a linear or logistic regression.

    Review: TCW; 3 April 2025

    arguments:
        table_predictor_intercept (object): Pandas data-frame table of data
            with features and observations for regression
        indicator (int): numeric sequence and predictor feature
        name_source (str): name of predictor feature
        name_product (str): name of predictor feature
        parameters (dict): parameters for predictors from regression model
        standard_errors (dict): standard errors for predictors from regression
            model
        pvalues (dict): pvalues for predictors from regression model
        variance_inflation (bool): whether to determine variance inflation
            factor for the predictor feature

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Determine code.
    code = str("predictor_" + str(indicator))
    # Collect information.
    pail = dict()
    pail["predictor"] = dict()
    pail["entries"] = [
        str(code + "_name"),
        str(code + "_parameter"),
        str(code + "_error"),
        str(code + "_ci99_low"),
        str(code + "_ci99_high"),
        str(code + "_pvalue"),
        str(code + "_variance_inflation"),
    ]
    pail["predictor"][str(code + "_name")] = name_product
    pail["predictor"][str(code + "_parameter")] = float(
        parameters[name_source]
    )
    pail["predictor"][str(code + "_error")] = float(
        standard_errors[name_source]
    )
    pail_confidence = pdesc.determine_95_99_confidence_intervals_ranges(
        estimate=parameters[name_source],
        standard_error=standard_errors[name_source],
    )
    pail["predictor"][str(code + "_ci99_low")] = float(
        pail_confidence["range_99_low"]
    )
    pail["predictor"][str(code + "_ci99_high")] = float(
        pail_confidence["range_99_high"]
    )
    pail["predictor"][str(code + "_pvalue")] = float(pvalues[name_source])
    if (variance_inflation):
        inflation_value = float(
            statsmodels.stats.outliers_influence.variance_inflation_factor(
                table_predictor_intercept.to_numpy(),
                indicator
            )
        )
        pail["predictor"][str(code + "_variance_inflation")] = float(
            inflation_value
        )
    else:
        pail["predictor"][str(code + "_variance_inflation")] = float("nan")
        pass

    # Return information.
    return pail


def parse_regression_intercept_predictors(
    table_predictor_intercept=None,
    features_predictor=None,
    features_sequence=None,
    name_intercept_source=None,
    name_intercept_product=None,
    parameters=None,
    standard_errors=None,
    pvalues=None,
):
    """
    Parse raw information about predictor independent variables from a model
    of a linear or logistic regression.

    Review: TCW; 3 April 2025

    arguments:
        table_predictor_intercept (object): Pandas data-frame table of data
            with features and observations for regression
        features_predictor (list<str>): names of features to include in
            regression model as predictor independent variables
        features_sequence (dict<int>): numeric sequence and names of features
        name_intercept_source (str): name of parameters for intercept
        name_intercept_product (str): name of parameters for intercept
        parameters (dict): parameters for predictors from regression model
        standard_errors (dict): standard errors for predictors from regression
            model
        pvalues (dict): pvalues for predictors from regression model

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table_predictor_intercept = table_predictor_intercept.copy(deep=True)
    features_predictor = copy.deepcopy(features_predictor)
    features_sequence = copy.deepcopy(features_sequence)
    parameters = copy.deepcopy(parameters)
    standard_errors = copy.deepcopy(standard_errors)
    pvalues = copy.deepcopy(pvalues)

    # Collect information.
    pail = dict()
    pail["intercept_predictors"] = dict()
    pail["entries"] = list()

    # Intercept.
    if (
        (name_intercept_source in parameters.index) and
        (name_intercept_source in standard_errors.index) and
        (name_intercept_source in pvalues.index)
    ):
        # Invert dictionary for convenience.
        #sequence_features = {v : k for k, v in features_sequence.items()}
        indicator = features_sequence[name_intercept_product]
        # Collect information.
        pail_predictor = parse_regression_model_predictor(
            table_predictor_intercept=table_predictor_intercept,
            indicator=indicator,
            name_source=name_intercept_source,
            name_product=name_intercept_product,
            parameters=parameters,
            standard_errors=standard_errors,
            pvalues=pvalues,
            variance_inflation=False,
        )
        pail["intercept_predictors"].update(pail_predictor["predictor"])
        pail["entries"].extend(pail_predictor["entries"])
        pass

    # Main predictor features.
    for feature in features_predictor:
        indicator = features_sequence[feature]
        # Collect information.
        pail_predictor = parse_regression_model_predictor(
            table_predictor_intercept=table_predictor_intercept,
            indicator=indicator,
            name_source=feature,
            name_product=feature,
            parameters=parameters,
            standard_errors=standard_errors,
            pvalues=pvalues,
            variance_inflation=True,
        )
        pail["intercept_predictors"].update(pail_predictor["predictor"])
        pail["entries"].extend(pail_predictor["entries"])
        pass

    # Return information.
    return pail


def parse_organize_regression_intercept_predictors(
    table_predictor_intercept=None,
    features_predictor=None,
    parameters=None,
    standard_errors=None,
    pvalues=None,
    report=None,
):
    """
    Parse and organize information about predictor independent variables from a
    model of a linear or logistic regression.

    Review: TCW; 2 April 2025

    arguments:
        table_predictor_intercept (object): Pandas data-frame table of data
            with features and observations for regression
        features_predictor (list<str>): names of features to include in
            regression model as predictor independent variables
        parameters (dict): parameters for predictors from regression model
        standard_errors (dict): standard errors for predictors from regression
            model
        pvalues (dict): pvalues for predictors from regression model
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table_predictor_intercept = table_predictor_intercept.copy(deep=True)
    features_predictor = copy.deepcopy(features_predictor)
    features_predictor_intercept = copy.deepcopy(features_predictor)

    # Determine count and sequence of predictors.
    if ("const" in table_predictor_intercept.columns.to_list()):
        features_predictor_intercept.insert(0, "intercept")
        pass
    count_intercept_predictors = len(features_predictor_intercept)
    features_sequence = define_features_sequence(
        features=features_predictor_intercept,
    )

    # Parse information for predictor features from regression model.
    pail = parse_regression_intercept_predictors(
        table_predictor_intercept=table_predictor_intercept,
        features_predictor=features_predictor,
        features_sequence=features_sequence,
        name_intercept_source="const",
        name_intercept_product="intercept",
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
    )

    # Collect information.
    pail["intercept_predictors"]["count_intercept_predictors"] = (
        count_intercept_predictors
    )
    pail["entries"].insert(0, "count_intercept_predictors")

    # Return information.
    return pail


##########
# Regression models


def regress_continuous_linear_ordinary_least_squares(
    table=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Regress a response feature on a quantitative, continuous, ratio or interval
    scale of measurement against either a single or multiple predictor features
    by a linear ordinary least squares (OLS) model, applying an implementation
    from the StatsModels package in Python.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Format of information for response dependent variable is a one-dimensional
    vector of scalar values on a quantitative, continuous, ratio or interval
    scale of measurement.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Implementation: 'statsmodels.regression.linear_model.OLS()'

    References:

    1. Documentation for implementation in StatsModels
       - title: 'statsmodels.regression.linear_model.OLS'
       - site: https://www.statsmodels.org/stable/generated/statsmodels.
                regression.linear_model.OLS.html

    Review: TCW; 3 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Extract values of response and predictor variables.
    # It is possible to pass values of response and predictor features either
    # as NumPy arrays or as Pandas Series and data-frame table respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable
    # names.
    #values_response = table[feature_response].to_numpy()

    # Keep information for predictor features in Pandas data-frame table to
    # preserve names of variables.
    #values_predictor = table.loc[ :, features_predictor_fixed].to_numpy()
    table_predictor = table.loc[ :, features_predictor_fixed].copy(deep=True)

    # Introduce constant intercept to predictor features.
    # If any column in the table of information for predictor features already
    # has constant values across observations, then the function skips it by
    # default. It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_predictor_intercept = statsmodels.api.add_constant(
        table_predictor,
        prepend=True, # insert intercept constant before other features
        has_constant="add", # force introduction of new intercept constant
    )
    columns_predictor = copy.deepcopy(
        table_predictor_intercept.columns.to_list()
    )
    #matrix_predictor = table.to_numpy()
    # Define model.
    model = statsmodels.api.OLS(
        table[feature_response],
        table_predictor_intercept,
        missing="drop",
    )
    # Fit the model.
    handle_model = model.fit()

    # Collect information.
    pail = dict()
    pail["model"] = handle_model
    pail["table_predictor_intercept"] = table_predictor_intercept

    # Return information.
    return pail


def regress_continuous_linear_mixed_effects_random_intercepts_slopes(
    table=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    groups_random=None,
    method_fit=None,
    report=None,
):
    """
    Regress a response feature on a quantitative, continuous, ratio or interval
    scale of measurement against either a single or multiple predictor features
    by a model with mixed fixed and random effects, applying an implementation
    from the StatsModels package in Python.

    The model accommodates a random effect of different intercepts between
    groups of observations. It also accommodates random effects of slope
    coefficients that differ between groups of observations, requiring these
    random slopes to be correlated, I think with the predictor features for
    random effects. It is possible to force these random slopes to be
    correlated or to force these random slopes to be independent and
    uncorrelated. Refer to the documentation of the 'MixedLM' method in the
    StatsModels package in Python.

    The implementation in StatsModels allows additional parameters through the
    'fit' method (https://www.statsmodels.org/dev/generated/statsmodels.
    regression.mixed_linear_model.MixedLM.fit.html).

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Format of information for response dependent variable is a one-dimensional
    vector of scalar values on a quantitative, continuous, ratio or interval
    scale of measurement.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Implementation: 'statsmodels.regression.mixed_linear_model.MixedLM()'

    Fit methods
    ’newton’: Newton-Raphson # unavailable
    ‘nm’: Nelder-Mead
    ’bfgs’: Broyden-Fletcher-Goldfarb-Shanno (BFGS) # faulty error
    ’lbfgs’: limited-memory BFGS with optional box constraints # faulty error
    ’powell’: modified Powell’s method # seems robust to colinearity error
    ’cg’: conjugate gradient
    ’ncg’: Newton-conjugate gradient # unavailable
    ’basinhopping’: global basin-hopping solver
    ’minimize’: generic wrapper of scipy minimize (BFGS by default)

    References:

    1. Documentation for implementation in StatsModels
       - title: 'statsmodels.regression.mixed_linear_model.MixedLM'
       - site: https://www.statsmodels.org/stable/generated/statsmodels.
                regression.mixed_linear_model.MixedLM.html

    2. Documentation for model fit method in StatsModels
       - title: 'statsmodels.regression.mixed_linear_model.MixedLM.fit'
       - site: https://www.statsmodels.org/dev/generated/statsmodels.
                regression.mixed_linear_model.MixedLM.fit.html

    3. Report of issue in implementation with failure to converge in model fit
       - title: 'Singular matrix error on mixedlm fit #7051'
       - site: https://github.com/statsmodels/statsmodels/issues/7051
       - site: https://github.com/statsmodels/statsmodels/pull/7052
       - note:
          - 'I get LinAlgError: Singular matrix when I try to fit a mixedlm
             model to my data, even though there are no colinear components.'

    Review: TCW; 9 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        groups_random (str): name of column in data table for identifiers or
            designations of groups of observations for which to allow random
            effects in the regression model
        method_fit (list<str>): names of methods to apply for fitting the
            model, especially with random slopes some models can fail to fit
            the model to convergence
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Extract values of response and predictor variables.
    # It is possible to pass values of response and predictor features either
    # as NumPy arrays or as Pandas Series and data-frame table respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable
    # names.
    #values_response = table[feature_response].to_numpy()

    # Keep information for predictor features in Pandas data-frame table to
    # preserve names of variables.
    #values_predictor = table.loc[ :, features_predictor_fixed].to_numpy()
    #matrix_predictor = table.to_numpy()
    table_predictor_fixed = (
        table.loc[ :, features_predictor_fixed].copy(deep=True)
    )
    if (
        (features_predictor_random is not None) and
        (len(features_predictor_random) > 0) and
        (str(features_predictor_random[0]).strip().lower() != "none")
    ):
        table_predictor_random = (
            table.loc[ :, features_predictor_random].copy(deep=True)
        )
    else:
        table_predictor_random = None

    # Do not introduce constant intercept to predictor features.
    # For part of the random effects, the model will create separate intercepts
    # for each group of observations.

    # Define model.
    # Because it is the default to include random intercepts for the groups of
    # observations, it is unnecessary to specify this explicitly using the
    # argument "exog_re=table_predictor_intercept["const"],".
    if (table_predictor_random is None):
        # Define model with random intercepts but not random slopes.
        model = statsmodels.api.MixedLM(
            table[feature_response],
            table_predictor_fixed,
            groups=table[groups_random],
            exog_re=None, # default to random intercept for each group
            missing="drop",
        )
        # Fit the model.
        handle_model = model.fit(
            reml=True, # Restricted Maximum Likelihood (REML) or ML
            method=method_fit,
        )
    else:
        # Define model with random intercepts and random slopes.
        model = statsmodels.api.MixedLM(
            table[feature_response],
            table_predictor_fixed,
            groups=table[groups_random],
            exog_re=table_predictor_random,
            missing="drop",
        )
        # Fit the model.
        handle_model = model.fit(
            reml=True, # Restricted Maximum Likelihood (REML) or ML
            method=method_fit,
            #maxiter=10000,
        )

        pass

    # Collect information.
    pail = dict()
    pail["model"] = handle_model
    pail["table_predictor_fixed"] = table_predictor_fixed

    # Return information.
    return pail


def regress_discrete_logistic_logit(
    table=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Regress a response feature on a discrete, binary scale of measurement
    against either a single or multiple predictor features by a generalized
    linear model (GLM) with the logit link function, applying an implementation
    from the StatsModels package in Python.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Format of information for response dependent variable is a one-dimensional
    vector of scalar values on a discrete, binary scale of measurement.
    [0, 1, 1, 0, 1, 1, 0, 0, 0, 1,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Complete or Quasi-Complete Separation of discrete (or binary) values of any
    independent variable between discrete values of the dependent variable can
    prevent convergence of the logistic regression either by maximum likelihood
    or regularized maximum likelihood methods. Evaluate the relationships
    between dependent and independent variables to avoid any Complete or
    Quasi-Complete Separation.

    Implementation: 'statsmodels.discrete.discrete_model.Logit()'

    References:

    1. Documentation for implementation in StatsModels
       - title: 'statsmodels.discrete.discrete_model.Logit'
       - site: https://www.statsmodels.org/stable/generated/statsmodels.
                discrete.discrete_model.Logit.html

    Review: TCW; 3 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Extract values of response and predictor variables.
    # It is possible to pass values of response and predictor features either
    # as NumPy arrays or as Pandas Series and data-frame table respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable
    # names.
    #values_response = table[feature_response].to_numpy()

    # Keep information for predictor features in Pandas data-frame table to
    # preserve names of variables.
    #values_predictor = table.loc[ :, features_predictor_fixed].to_numpy()
    table_predictor = table.loc[ :, features_predictor_fixed].copy(deep=True)

    # Introduce constant intercept to predictor features.
    # If any column in the table of information for predictor features already
    # has constant values across observations, then the function skips it by
    # default. It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_predictor_intercept = statsmodels.api.add_constant(
        table_predictor,
        prepend=True, # insert intercept constant before other features
        has_constant="add", # force introduction of new intercept constant
    )
    columns_predictor = copy.deepcopy(
        table_predictor_intercept.columns.to_list()
    )
    #matrix_predictor = table.to_numpy()
    # Define model.
    model = statsmodels.api.Logit(
        table[feature_response],
        table_predictor_intercept,
        missing="drop",
    )
    # Fit the model.
    handle_model = model.fit(
        maxiter=100, # maximal count of iterations to convergence
    )

    # Regularization of the fit method can help to avoid overfitting and promote
    # convergence.
    # Fit the model using regularized maximum likelihood.
    #handle_model = model.fit_regularized(
    #    method="l1", # solver method
    #    maxiter=100,
    #)

    # Collect information.
    pail = dict()
    pail["model"] = handle_model
    pail["table_predictor_intercept"] = table_predictor_intercept

    # Return information.
    return pail


##########
# Manage regression analysis.


def manage_regression_continuous_linear_ordinary_least_squares(
    table=None,
    index_columns=None,
    index_rows=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Manage regression analysis according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 3 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_regression (str): name of type of regression model to use, either
            'continuous_ols', 'continuous_mixed_random_intercepts',
            'continuous_mixed_random_intercepts_slopes',
            'continuous_gls', 'discrete_logit', 'discrete_probit',
            'discrete_poisson', 'discrete_negbinom'
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_fixed = copy.deepcopy(features_predictor_fixed)

    # Regress.
    pail_regression = regress_continuous_linear_ordinary_least_squares(
        table=table,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        report=report,
    )
    # Report attributes.
    #print(dir(pail_regression))
    # Extract residuals.
    residuals = pail_regression["model"].resid

    # Organize information from regression for predictor features.
    parameters = pandas.Series(data=pail_regression["model"].params)
    standard_errors = pandas.Series(data=pail_regression["model"].bse)
    pvalues = pandas.Series(data=pail_regression["model"].pvalues)
    pail_predictors = parse_organize_regression_intercept_predictors(
        table_predictor_intercept=pail_regression["table_predictor_intercept"],
        features_predictor=features_predictor_fixed,
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
        report=report,
    )

    # Collect information.
    pail = dict()
    # Regression model as a whole.
    pail["entries_model"] = [
        "implementation",
        "count_observations_table",
        "count_observations_model",
        "degrees_freedom_model",
        "degrees_freedom_residual",
        "r_square",
        "r_square_adjust",
        "r_square_pseudo",
        "log_likelihood",
        "akaike",
        "bayes",
        "condition",
    ]
    pail["implementation"] = str(
        "statsmodels.regression.linear_model.OLS()"
    )
    pail["count_observations_table"] = int(table.shape[0])
    pail["count_observations_model"] = int(pail_regression["model"].nobs)
    pail["degrees_freedom_model"] = int(
        pail_regression["model"].df_model
    )
    pail["degrees_freedom_residual"] = int(
        pail_regression["model"].df_resid
    )
    pail["r_square"] = pail_regression["model"].rsquared
    pail["r_square_adjust"] = pail_regression["model"].rsquared_adj
    pail["r_square_pseudo"] = float("nan") # hold place for logistic regression
    pail["log_likelihood"] = pail_regression["model"].llf
    pail["akaike"] = pail_regression["model"].aic
    pail["bayes"] = pail_regression["model"].bic
    pail["condition"] = pail_regression["model"].condition_number
    # Intercept and predictors as parts of regression model.
    pail["entries_intercept_predictors"] = pail_predictors["entries"]
    pail.update(
        pail_predictors["intercept_predictors"]
    )
    pail["table_summary"] = pail_regression["model"].summary()

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_regression_continuous_linear_ordinary_least_squares()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for regression model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("summary from regression model:")
        print(pail_regression["model"].summary())
        #print(dir(pail_regression))
        #print(pail_regression.params)
        #print(pail_regression.pvalues)
        pass

    # Return information.
    return pail


# TODO: TCW; 6 April 2025
# enhance the information extracted from the model?
def manage_regression_continuous_linear_mixed_effects(
    table=None,
    index_columns=None,
    index_rows=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    groups_random=None,
    method_fit=None,
    report=None,
):
    """
    Manage regression analysis according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 6 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_regression (str): name of type of regression model to use, either
            'continuous_ols', 'continuous_mixed_random_intercepts',
            'continuous_mixed_random_intercepts_slopes',
            'continuous_gls', 'discrete_logit', 'discrete_probit',
            'discrete_poisson', 'discrete_negbinom'
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        groups_random (str): name of column in data table for identifiers or
            designations of groups of observations for which to allow random
            effects in the regression model
        method_fit (list<str>): names of methods to apply for fitting the
            model, especially with random slopes some models can fail to fit
            the model to convergence
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_fixed = copy.deepcopy(features_predictor_fixed)
    features_predictor_random = copy.deepcopy(features_predictor_random)

    # Regress.
    pail_regression = (
        regress_continuous_linear_mixed_effects_random_intercepts_slopes(
            table=table,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            groups_random=groups_random,
            method_fit=method_fit,
            report=report,
    ))
    # Report attributes.
    #print(dir(pail_regression["model"]))
    # Extract residuals.
    residuals = pail_regression["model"].resid

    print(pail_regression["model"].summary())

    # Organize information from regression for predictor features.
    parameters = pandas.Series(data=pail_regression["model"].params)
    standard_errors = pandas.Series(data=pail_regression["model"].bse)
    pvalues = pandas.Series(data=pail_regression["model"].pvalues)
    pail_predictors = parse_organize_regression_intercept_predictors(
        table_predictor_intercept=pail_regression["table_predictor_fixed"],
        features_predictor=features_predictor_fixed,
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
        report=report,
    )

    # Collect information.
    pail = dict()
    # Regression model as a whole.
    pail["entries_model"] = [
        "implementation",
        "count_observations_table",
        "count_observations_model",
        "degrees_freedom_model",
        "degrees_freedom_residual",
        "r_square",
        "r_square_adjust",
        "r_square_pseudo",
        "log_likelihood",
        "akaike",
        "bayes",
        "condition",
    ]
    pail["implementation"] = str(
        "statsmodels.regression.linear_model.OLS()"
    )
    pail["count_observations_table"] = int(table.shape[0])
    pail["count_observations_model"] = int(pail_regression["model"].nobs)
    pail["degrees_freedom_model"] = int(
        pail_regression["model"].df_modelwc # what does this mean?
    )
    pail["degrees_freedom_residual"] = int(
        pail_regression["model"].df_resid
    )
    pail["degrees_freedom_residual"] = float("nan")
    #pail["r_square"] = pail_regression["model"].rsquared
    pail["r_square"] = float("nan")
    #pail["r_square_adjust"] = pail_regression["model"].rsquared_adj
    pail["r_square_adjust"] = float("nan")
    pail["r_square_pseudo"] = float("nan") # hold place for logistic regression
    pail["log_likelihood"] = pail_regression["model"].llf
    pail["akaike"] = pail_regression["model"].aic
    pail["bayes"] = pail_regression["model"].bic
    #pail["condition"] = pail_regression["model"].condition_number
    pail["condition"] = float("nan")
    # Intercept and predictors as parts of regression model.
    pail["entries_intercept_predictors"] = pail_predictors["entries"]
    pail.update(
        pail_predictors["intercept_predictors"]
    )
    pail["table_summary"] = pail_regression["model"].summary()

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_regression_continuous_linear_mixed_effects()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for regression model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("summary from regression model:")
        print(pail_regression["model"].summary())
        #print(dir(pail_regression))
        #print(pail_regression.params)
        #print(pail_regression.pvalues)
        pass

    # Return information.
    return pail


def manage_regression_discrete_logistic_logit(
    table=None,
    index_columns=None,
    index_rows=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Manage regression analysis according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 3 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_regression (str): name of type of regression model to use, either
            'continuous_ols', 'continuous_mixed_random_intercepts',
            'continuous_mixed_random_intercepts_slopes',
            'continuous_gls', 'discrete_logit', 'discrete_probit',
            'discrete_poisson', 'discrete_negbinom'
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_fixed = copy.deepcopy(features_predictor_fixed)

    # Regress.
    pail_regression = regress_discrete_logistic_logit(
        table=table,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        report=report,
    )
    # Extract residuals.
    residuals = pail_regression["model"].resid_generalized

    # Organize information from regression for predictor features.
    parameters = pandas.Series(data=pail_regression["model"].params)
    standard_errors = pandas.Series(data=pail_regression["model"].bse)
    pvalues = pandas.Series(data=pail_regression["model"].pvalues)
    pail_predictors = parse_organize_regression_intercept_predictors(
        table_predictor_intercept=pail_regression["table_predictor_intercept"],
        features_predictor=features_predictor_fixed,
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
        report=report,
    )

    # Collect information.
    pail = dict()
    # Regression model as a whole.
    pail["entries_model"] = [
        "implementation",
        "count_observations_table",
        "count_observations_model",
        "degrees_freedom_model",
        "degrees_freedom_residual",
        "r_square",
        "r_square_adjust",
        "r_square_pseudo",
        "log_likelihood",
        "akaike",
        "bayes",
        "condition",
    ]
    pail["implementation"] = str(
        "statsmodels.discrete.discrete_model.Logit()"
    )
    pail["count_observations_table"] = int(table.shape[0])
    pail["count_observations_model"] = int(pail_regression["model"].nobs)
    pail["degrees_freedom_model"] = int(
        pail_regression["model"].df_model
    )
    pail["degrees_freedom_residual"] = int(
        pail_regression["model"].df_resid
    )
    pail["r_square"] = float("nan") # hold place for linear regression
    pail["r_square_adjust"] = float("nan") # hold place for linear regression
    pail["r_square_pseudo"] = pail_regression["model"].prsquared
    pail["log_likelihood"] = pail_regression["model"].llf
    pail["akaike"] = pail_regression["model"].aic
    pail["bayes"] = pail_regression["model"].bic
    pail["condition"] = float("nan") # hold place for linear regression
    # Intercept and predictors as parts of regression model.
    pail["entries_intercept_predictors"] = pail_predictors["entries"]
    pail.update(
        pail_predictors["intercept_predictors"]
    )
    pail["table_summary"] = pail_regression["model"].summary()

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_regression_discrete_logistic_logit()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for regression model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("summary from regression model:")
        print(pail_regression["model"].summary())
        #print(dir(pail_regression))
        #print(pail_regression.params)
        #print(pail_regression.pvalues)
        pass

    # Return information.
    return pail


# See document string for informative notes about regressions.
def determine_type_regression_analysis(
    table=None,
    index_columns=None,
    index_rows=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    groups_random=None,
    method_fit=None,
    report=None,
):
    """
    Determine the type of regression analysis according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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
    Notes

    Scale transformation of predictor independent variables
    - Normally do not transform the scale of the dependent variable response
       feature. In some cases, this can make it difficult for the fit of the
       model to converge (Linear Regression with Mixed Effects).
    - In multiple linear regression, it is common to standardize the scale of
       of predictor independent variables that are on a quantitative,
       continuous scale of measurement, such as ratio or interval.
    - Standardizing the scale of these predictor independent variables makes it
       more reasonable and convenient to compare the magnitudes of estimates
       for parameter coefficients between these different predictor independent
       variables.
    - When the regression model also includes additional predictor independent
       variables that are on a binary scale of measurement, such as dummy
       variables for categories, it is most reasonable to leave these on the
       binary scale.
    - The rough interpretation of the slope parameter coefficient for a binary
       predictor independent variable is the difference in means of the
       response dependent variable between the two groups of observations
       corresponding to the two values of the binary predictor.
    - In regression with multiple predictor independent variables, the rough
       interpretation of the slope for a binary predictor is the same, while
       holding the other predictor covariates constant, perhaps at their
       respective means.

    References:

    1. Forum discussion about standardizing the scales of dummy variables
       title: 'Standardizing dummy variable in multiple linear regression?'
       site: https://stats.stackexchange.com/questions/414053/standardizing-
              dummy-variable-in-multiple-linear-regression

    2. Forum discussion about interpretation of regression for binary variable
       title: 'Binary variable in multiple linear regression?'
       site: https://www.reddit.com/r/AskStatistics/comments/kaiwk2/
              binary_variable_in_multiple_linear_regression/

    ----------
    Notes

    Interpretation
    - The parameter coefficient for the intercept is the model's prediction of
       the response when all predictors have values of zero (0). If the model
       includes dummy predictors to represent features with categorical values,
       then the value of the intercept corresponds to the reference category
       that corresponds to values of zero (0) for dummy predictors that
       represent all the other categories. To represent "n" categories", "n-1"
       dummy predictors are necessary.
    - The parameter coefficient for a continuous predictor indicates that, when
       holding all other covariate predictors constant at any of their
       respective values, a one unit increase in the predictor of interest
       corresponds on average to the coefficient's unit change in the response.

    ----------
    Notes

    Types of Regression Models
    - Continuous response feature
     - Linear Ordinary Least Squares (OLS)
      - requires homoscedasticity of variance in response across ranges of all
         predictors
     - Linear with Mixed Effects
      - types of random effects in the regression model:
       - 1. random intercept different between each group of observations
        - this is the simplest degree of random effects
       - 2. random slope coefficient between each group of observations, with
             respect to one of the predictor features with fixed effects, such
             as time
        - 2.1. random slope coefficients can be correlated
         - this is the intermediate degree of random effects
        - 2.2. random slope coefficients can be uncorrelated and independent
         - this is the extreme, complex degree of random effects
     - Linear Generalized Least Squares (GLS)
      - allows heteroscedasticity of variance in response across ranges of all
         predictors
    - Linear
    - Discrete response feature
     - Generalized Linear Model with Logit link function
      - response is on discrete binary scale of measurement
      - link function is the logarithm of the odds ratio
     - Generalized Linear Model with Probit (probability unit) link function
      - response is on discrete binary scale of measurement
      - link function bases on the cumulative density function of the standard
         normal distribution
      - less accurate and less interpretable for odds ratios
     - Generalized Linear Model with Poisson link function
      - response is on discrete ordinal scale of measurement
     - Generalized Linear Model with Negative Binomial link function
      - response is on discrete ordinal scale of measurement
     - General Linear Model with Mixed Effects
      -

    References:

    1. Types of linear regression models in StatsModels
       title: 'Linear Regression'
       site: https://www.statsmodels.org/stable/regression.html

    2. Linear Models with Mixed Effects in StatsModels
       title: 'Linear Mixed Effects Models'
       site: https://www.statsmodels.org/stable/mixed_linear.html
       site: https://www.statsmodels.org/stable/examples/notebooks/generated/
              mixed_lm_example.html

    2. Types of discrete response regression models in StatsModels
       title: 'Regression with Discrete Dependent Variable'
       site: https://www.statsmodels.org/stable/discretemod.html

    3. Types of regression with mixed effects in StatsModels
       title: 'Linear Mixed Effects Models'
       site: https://www.statsmodels.org/stable/mixed_linear.html

    4. Forum discussion of differences between logit, probit, and other link
       functions for regression on discrete responses
       title: 'Difference between logit and probit models'
       site: https://stats.stackexchange.com/questions/20523/difference-
              between-logit-and-probit-models/30909#30909

    5. StatsModels Generalized Linear Models with Mixed Effects
       title: 'Generalized Linear Mixed Effects Models'
       site: https://www.statsmodels.org/stable/mixed_glm.html


    ----------
    Notes

    Linear regression against predictors with mixed fixed and random effects
    - ...

    References:

    1. Article about linear regression analysis with mixed fixed and random
       effects
       title: 'Linear Mixed Effect Models in python using mtcars'
       site: https://medium.com/@josef.waples/linear-mixed-effect-models-in-
              python-using-mtcars-064c6c3e5b32
       note:
       - This article also provides a clear and insightful interpretation of
          the results of regression analysis.

    ----------
    Review: TCW; 3 April 2025
    - On 3 April 2025, TCW confirmed that extracted values from linear OLS
       and discrete generalized Logit regression models, intercept, and
       predictors in the summary table match those reported directly in the
       summary from the implementations of the respective regression models on
       demonstration data.

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_regression (str): name of type of regression model to use, either
            'continuous_ols', 'continuous_mixed_random_intercepts',
            'continuous_mixed_random_intercepts_slopes',
            'continuous_gls', 'discrete_logit', 'discrete_probit',
            'discrete_poisson', 'discrete_negbinom'
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        groups_random (str): name of column in data table for identifiers or
            designations of groups of observations for which to allow random
            effects in the regression model
        method_fit (list<str>): names of methods to apply for fitting the
            model, especially with random slopes some models can fail to fit
            the model to convergence
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Determine type of regression.
    if (type_regression == "continuous_ols"):
        # Perform continuous linear ordinary least squares regression.
        pail = manage_regression_continuous_linear_ordinary_least_squares(
            table=table,
            index_columns=index_columns,
            index_rows=index_rows,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            report=report,
        )
        pass
    elif (
        type_regression in [
            "continuous_mixed",
        ]
    ):
        # Perform continuous linear regression with mixed effects.
        pail = manage_regression_continuous_linear_mixed_effects(
            table=table,
            index_columns=index_columns,
            index_rows=index_rows,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            groups_random=groups_random,
            method_fit=method_fit,
            report=report,
        )
        pass
    elif (type_regression == "discrete_logit"):
        # Perform logistic regression.
        pail = manage_regression_discrete_logistic_logit(
            table=table,
            index_columns=index_columns,
            index_rows=index_rows,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            report=report,
        )
        pass
    else:
        pail = dict()
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "determine_type_regression_analysis()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail


##########
# Analysis of Variance (ANOVA)


def analyze_variance_anova_two_way_repeat_measures(
    table=None,
    feature_response=None,
    features_predictor_within=None,
    subject=None,
    report=None,
):
    """
    Perform Analysis of Variance (ANOVA) for a response feature on a
    quantitative, continuous, ratio or interval scale of measurement against
    either one or two predictor features with repeated measures on groups of
    samples or observations, applying an implementation from the 'pingouin'
    package in Python.

    Notice that this analysis only allows predictor features that correspond to
    the repeated measurements of the same subject.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Format of information for response dependent variable is a one-dimensional
    vector of scalar values on a quantitative, continuous, ratio or interval
    scale of measurement.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Implementation: 'pingouin.rm_anova()'

    References:

    1. Documentation for implementation in 'pingouin'
       - title: 'pingouin.rm_anova()'
       - site: https://pingouin-stats.org/build/html/generated/pingouin.
               rm_anova.html

    Review: TCW; 7 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations
        feature_response (str): name of column in data table for feature
            variable to include in ANOVA model as response dependent
            variable
        features_predictor_within (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements within separate
            subjects
        subject (str): name of column in data table for identifier of subject
            or individual source of pairs of samples or observations
            corresponding to repeat measures
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from ANOVA

    """

    # Analysis of Variance (ANOVA).
    handle_anova = pingouin.rm_anova(
        data=table,
        dv=feature_response,
        within=features_predictor_within,
        subject=subject,
        effsize="ng2",
        detailed=True, # whether to return full summary table
    )
    # T-tests.
    handle_t = pingouin.pairwise_tests(
        data=table,
        dv=feature_response,
        between=None,
        within=features_predictor_within,
        subject=subject,
        within_first=False,
    )
    # Collect information.
    pail = dict()
    pail["anova"] = handle_anova
    pail["t_test"] = handle_t

    # Return information.
    return pail


def analyze_variance_anova_two_way_repeat_measures_mixed(
    table=None,
    feature_response=None,
    features_predictor_between=None,
    features_predictor_within=None,
    subject=None,
    report=None,
):
    """
    Perform Analysis of Variance (ANOVA) for a response feature on a
    quantitative, continuous, ratio or interval scale of measurement against
    either one or two predictor features with repeated measures on groups of
    samples or observations, applying an implementation from the 'pingouin'
    package in Python.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Format of information for response dependent variable is a one-dimensional
    vector of scalar values on a quantitative, continuous, ratio or interval
    scale of measurement.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Implementation: 'pingouin.mixed_anova()'

    References:

    1. Documentation for implementation in 'pingouin'
       - title: 'pingouin.mixed_anova()'
       - site: https://pingouin-stats.org/build/html/generated/pingouin.
               mixed_anova.html

    Review: TCW; 7 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations
        feature_response (str): name of column in data table for feature
            variable to include in ANOVA model as response dependent
            variable
        features_predictor_between (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements between
            separate groups of observations
        features_predictor_within (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements within separate
            groups of observations
        subject (str): name of column in data table for identifier of subject
            or individual source of pairs of samples or observations
            corresponding to repeat measures
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from ANOVA

    """

    # Analysis of Variance (ANOVA).
    handle_anova = pingouin.mixed_anova(
        data=table,
        dv=feature_response,
        between=features_predictor_between[0],
        within=features_predictor_within[0],
        subject=subject,
        effsize="ng2",
    )
    # T-tests.
    handle_t_within_first = pingouin.pairwise_tests(
        data=table,
        dv=feature_response,
        between=features_predictor_between[0],
        within=features_predictor_within[0],
        subject=subject,
        within_first=True,
    )
    handle_t_between_first = pingouin.pairwise_tests(
        data=table,
        dv=feature_response,
        between=features_predictor_between[0],
        within=features_predictor_within[0],
        subject=subject,
        within_first=False,
    )
    # Collect information.
    pail = dict()
    pail["anova"] = handle_anova
    pail["t_test_within"] = handle_t_within_first
    pail["t_test_between"] = handle_t_between_first

    # Return information.
    return pail


##########
# Manage analysis of variance (ANOVA).


def manage_anova_two_way_repeat_measures(
    table=None,
    index_columns=None,
    index_rows=None,
    type_anova=None,
    formula_text=None,
    feature_response=None,
    features_predictor_between=None,
    features_predictor_within=None,
    subject=None,
    report=None,
):
    """
    Manage Analysis of Variance (ANOVA) according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 7 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_anova (str): name of type of analysis of variance (ANOVA) to use,
            either '2way_repeat', or '2way_repeat_mixed'
        formula_text (str): human readable formula for ANOVA model, treated as
            a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in ANOVA model as response dependent
            variable
        features_predictor_between (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements between
            separate groups of observations
        features_predictor_within (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements within separate
            groups of observations
        subject (str): name of column in data table for identifier of subject
            or individual source of pairs of samples or observations
            corresponding to repeat measures
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from ANOVA

    """

    # Copy information.
    table = table.copy(deep=True)
    #features_predictor_between = copy.deepcopy(features_predictor_between)
    features_predictor_within = copy.deepcopy(features_predictor_within)
    # Combine predictor features.
    #features_predictor_within.extend(features_predictor_between)

    # Analysis of Variance (ANOVA).
    pail = (
        analyze_variance_anova_two_way_repeat_measures(
            table=table,
            feature_response=feature_response,
            features_predictor_within=features_predictor_within,
            subject=subject,
            report=report,
    ))
    # Report attributes.
    #print(dir(pail["anova"]))
    #print(dir(pail["t_test"]))
    # Extract attributes.
    #degrees_freedom_numerator = pail["anova"].ddof1
    #degrees_freedom_denominator = pail["anova"].ddof2

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_anova_two_way_repeat_measures()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for ANOVA model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("type ANOVA: " + type_anova)
        putly.print_terminal_partition(level=5)
        print("summary from ANOVA:")
        print(pail["anova"].round(3))
        putly.print_terminal_partition(level=5)
        print("summary from T-tests:")
        print(pail["t_test"])
        pass

    # Return information.
    return pail


def manage_anova_two_way_repeat_measures_mixed_effects(
    table=None,
    index_columns=None,
    index_rows=None,
    type_anova=None,
    formula_text=None,
    feature_response=None,
    features_predictor_between=None,
    features_predictor_within=None,
    subject=None,
    report=None,
):
    """
    Manage Analysis of Variance (ANOVA) according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 7 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_anova (str): name of type of analysis of variance (ANOVA) to use,
            either '2way_repeat', or '2way_repeat_mixed'
        formula_text (str): human readable formula for ANOVA model, treated as
            a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in ANOVA model as response dependent
            variable
        features_predictor_between (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements between
            separate groups of observations
        features_predictor_within (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements within separate
            groups of observations
        subject (str): name of column in data table for identifier of subject
            or individual source of pairs of samples or observations
            corresponding to repeat measures
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from ANOVA

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_between = copy.deepcopy(features_predictor_between)
    features_predictor_within = copy.deepcopy(features_predictor_within)

    # Analysis of Variance (ANOVA).
    pail = (
        analyze_variance_anova_two_way_repeat_measures_mixed(
            table=table,
            feature_response=feature_response,
            features_predictor_between=features_predictor_between,
            features_predictor_within=features_predictor_within,
            subject=subject,
            report=report,
    ))
    # Report attributes.
    #print(dir(pail["anova"]))
    #print(dir(pail["t_test"]))
    # Extract attributes.
    #degrees_freedom_numerator = pail["anova"].ddof1
    #degrees_freedom_denominator = pail["anova"].ddof2

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_anova_two_way_repeat_measures_mixed_effects()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for ANOVA model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("type ANOVA: " + type_anova)
        putly.print_terminal_partition(level=5)
        print("summary from ANOVA:")
        print(pail["anova"].round(3))
        putly.print_terminal_partition(level=5)
        print("summary from T-tests, within first:")
        print(pail["t_test_within"])
        putly.print_terminal_partition(level=5)
        print("summary from T-tests, between first:")
        print(pail["t_test_between"])
        pass

    # Return information.
    return pail


# See document string for informative notes about ANOVA.
def determine_type_analysis_variance_anova(
    table=None,
    index_columns=None,
    index_rows=None,
    type_anova=None,
    formula_text=None,
    feature_response=None,
    features_predictor_between=None,
    features_predictor_within=None,
    subject=None,
    report=None,
):
    """
    Determine the type of Analysis of Variance (ANOVA) according to parameters.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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
    Review: TCW; 7 April 2025

    arguments:

        table (object): Pandas data-frame table of data with features
            and observations
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_anova (str): name of type of analysis of variance (ANOVA) to use,
            either '2way_repeat', or '2way_repeat_mixed'
        formula_text (str): human readable formula for ANOVA model, treated as
            a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in ANOVA model as response dependent
            variable
        features_predictor_between (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements between
            separate groups of observations
        features_predictor_within (list<str>): names of columns in data table
            for feature variables to include in ANOVA model as predictor
            independent variables that distinguish measurements within separate
            groups of observations
        subject (str): name of column in data table for identifier of subject
            or individual source of pairs of samples or observations
            corresponding to repeat measures
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Determine type of regression.
    if (type_anova == "2way_repeat"):
        # Perform continuous linear ordinary least squares regression.
        pail = manage_anova_two_way_repeat_measures(
            table=table,
            index_columns=index_columns,
            index_rows=index_rows,
            type_anova=type_anova,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_between=features_predictor_between,
            features_predictor_within=features_predictor_within,
            subject=subject,
            report=report,
        )
    elif (type_anova == "2way_repeat_mixed"):
        # Perform continuous linear regression with mixed effects.
        pail = manage_anova_two_way_repeat_measures_mixed_effects(
            table=table,
            index_columns=index_columns,
            index_rows=index_rows,
            type_anova=type_anova,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_between=features_predictor_between,
            features_predictor_within=features_predictor_within,
            subject=subject,
            report=report,
        )
    else:
        pail = dict()
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "determine_type_analysis_variance_anova()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail


# Prepare information for write to file.


def prepare_text_summary_regression_anova(
    title=None,
    description_analysis=None,
    formula_text=None,
    description_response=None,
    description_groups_random=None,
    description_predictor=None,
    summary_1=None,
    summary_2=None,
    summary_3=None,
    report=None,
):
    """
    Prepare text information to print to terminal or write to file.

    Review: TCW; 8 April 2025

    arguments:
        title (str): title for summary text
        description_analysis (str): description
        formula_text (str): human readable formula, treated as a note for
            clarification
        description_response (str): description_response
        description_groups_random (str): description_groups_random
        description_predictor (str): description_predictor
        summary_1 (str): summary from analysis
        summary_2 (str): summary from analysis
        summary_3 (str): summary from analysis
        report (bool): whether to print reports

    raises:

    returns:
        (str): character string representation of information

    """

    summary_text = str(
        textwrap.dedent("""\

            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



        """) +
        str(title) +
        textwrap.dedent("""\

            --------------------------------------------------



        """) +
        str("Description analysis: \n" + description_analysis) +
        textwrap.dedent("""\

            ----------

        """) +
        str("Formula: \n" + formula_text) +
        textwrap.dedent("""\

            ----------

        """) +
        str(
            "Description of response feature: \n" +
            description_response
        ) +
        textwrap.dedent("""\

            ----------

        """) +
        str(
            "Description of groups of observations for random effects: \n" +
            description_groups_random
        ) +
        textwrap.dedent("""\

            ----------

        """) +
        str(
            "Description of predictor features: \n" +
            description_predictor
        ) +
        textwrap.dedent("""\

            ----------

        """) +
        textwrap.dedent("""\

            ----------

        """) +
        str("Summary from analysis: \n") +
        textwrap.dedent("""\



        """) +
        str(summary_1) +
        textwrap.dedent("""\



            ----------



        """) +
        str("Summary from analysis: \n") +
        textwrap.dedent("""\



        """) +
        str(summary_2) +
        textwrap.dedent("""\



            ----------



        """) +
        str("Summary from analysis: \n") +
        textwrap.dedent("""\



        """) +
        str(summary_3) +
        textwrap.dedent("""\



            ----------



        """) +
        textwrap.dedent("""\
            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------

        """)
    )

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "regression.py"
        )
        print("module: " + name_module)
        name_function = str(
            "prepare_text_summary_regression_anova()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("summary_text:")
        print(summary_text)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return summary_text


##########
# Collect record of information for regression analysis.


def collect_organize_record_regression_analysis(
    table=None,
    index_columns=None,
    index_rows=None,
    sequence=None,
    group=None,
    instance=None,
    name_instance=None,
    check_overall=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    groups_random=None,
    record_extra=None,
    delimiter_list_items=None,
    report=None,
):
    """
    Collect information within a record to describe a regression analysis.

    For an example on the application of this function in the context of
    extensive ancillary functionality to organize parameters and data and to
    drive multiple instances of regression concurrently in parallel, refer to
    the Python script 'script_drive_regressions_from_table_parameters.py'
    within the 'partner' repository. Within its subdirectory 'demonstration',
    the 'partner' repository also includes an example of the table of
    parameters that controls this Python script.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        sequence (int): sequential index for instance's name and sort order
        group (str): categorical group of instances
        instance (str): name or designator of instance
        name_instance (str): compound name for instance of parameters
        check_overall (bool): whether to execute the regression or create
            missing values
        type_regression (str): name of type of regression model to use
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        groups_random (str): name of column in data table for identifiers or
            designations of groups of observations for which to allow random
            effects in the regression model
        record_extra (dict): collection of secondary information to join to the
            record of primary information from the regression
        delimiter_list_items (str): delimiter to place between items when
            collapsing or flattening lists
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # pail (dict): collection of information about a regression analysis
    #     entries_parameter_instance (list<str>): names and sequence of entries
    #         with information about the instance of parameters for the
    #         regression
    #     entries_model (list<str>): names and sequence of entries with
    #         information about the regression model as a whole
    #     entries_intercept_predictors (list<str>): names and sequence of
    #         entries with with information about model's intercept and
    #         predictors
    #     record (dict<str>): record for collection and assembly as rows within
    #         a Pandas data-frame table, using the lists of entries to sort the
    #         columns in the table

    # Collect information.
    pail = dict()
    pail["entries_parameter_instance"] = [
        "sequence",
        "group",
        "instance",
        "name_instance",
        "check_overall",
        "type_regression",
        "formula_text",
        "feature_response",
        "features_predictor_fixed",
        "features_predictor_random",
        "groups_random",
    ]
    pail["entries_parameter_instance"].extend(list(record_extra.keys()))
    pail["record"] = dict()
    pail["record"]["sequence"] = sequence
    pail["record"]["group"] = group
    pail["record"]["instance"] = instance
    pail["record"]["name_instance"] = name_instance
    pail["record"]["check_overall"] = check_overall
    pail["record"]["type_regression"] = type_regression
    pail["record"]["formula_text"] = formula_text
    pail["record"]["feature_response"] = feature_response
    pail["record"]["features_predictor_fixed"] = delimiter_list_items.join(
        features_predictor_fixed
    )
    pail["record"]["features_predictor_random"] = delimiter_list_items.join(
        features_predictor_random
    )
    pail["record"]["groups_random"] = groups_random
    pail["record"].update(record_extra)
    # Determine whether to perform regression.
    if (check_overall):
        # Perform regression.
        pail_regression = determine_type_regression_analysis(
            table=table,
            index_columns=index_columns,
            index_rows=index_rows,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            groups_random=groups_random,
            report=report,
        )
        # Collect information.
        pail["entries_model"] = copy.deepcopy(pail_regression["entries_model"])
        pail["entries_intercept_predictors"] = copy.deepcopy(
            pail_regression["entries_intercept_predictors"]
        )
        pail["record"].update(pail_regression)
        del pail["record"]["entries_model"]
        del pail["record"]["entries_intercept_predictors"]
        del pail["record"]["table_summary"]
        pail["table_summary"] = pail_regression["table_summary"]
    else:
        # Collect information.
        pail["entries_model"] = None
        pail["entries_intercept_predictors"] = None
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "collect_organize_record_regression_analysis()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail




###############################################################################
# End
