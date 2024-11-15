"""
Supply functionality for process and analysis of data from transcriptomics.

This module 'scratch' is part of the 'transcriptomics' package within
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

# consider changing the name of this module to "select_changes"


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
import exercise.transcriptomics.organize_signal as exrosig
import exercise.transcriptomics.select_gene_sets as exrosel

###############################################################################
# Functionality


##########
# 1. Initialize directories for read of source and write of product files.
# There is a hierarchy in these functions to initialize directories to manage
# the hierarchical tree structure of sub-procedures.



##########
# 2. Read source information from file.


##########
# Plot charts.


def create_plot_chart_box(
    table=None,
    column_feature=None,
    column_group=None,
    report=None,
):
    """
    Create and plot a chart of the box type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 16 October 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement
        column_feature (str): name of column for values corresponding to a
            specific feature
        column_group (str): name of column in table to use for groups
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract names and indices of groups to preserve their original sequence.
    groups_sequence = copy.deepcopy(table[column_group].unique().tolist())
    sequence_groups = dict()
    index = 0
    for name in groups_sequence:
        sequence_groups[name] = index
        index += 1
        pass
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.set_index(
        [column_group],
        append=False,
        drop=True,
        inplace=True
    )
    # Split rows within table by factor columns.
    groups = table.groupby(
        level=column_group,
    )
    # Collect information.
    names_groups = list()
    values_groups = list()
    # Iterate on groups, apply operations, and collect information from each.
    for name_group, table_group in groups:
        # Copy information in table.
        table_group = table_group.copy(deep=True)
        # Extract information.
        values = table_group[column_feature].dropna().to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )
        # Collect information.
        names_groups.append(name_group)
        values_groups.append(values)
        pass
    # Sort information for groups and values.
    sequence_sort = list()
    for name in names_groups:
        sequence_sort.append(sequence_groups[name])
        pass
    names_groups_sort = [
        y for _,y in sorted(zip(sequence_sort, names_groups))
    ]
    values_groups_sort = [
        y for _,y in sorted(zip(sequence_sort, values_groups))
    ]

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_boxes_groups(
        values_groups=values_groups_sort,
        title_ordinate="scale-normal(gene signal)",
        title_abscissa="",
        titles_abscissa_groups=names_groups_sort,
        colors_groups=None,
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

    # Return information.
    return figure


def create_plot_chart_heatmap_individual(
    table=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 17 October 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values of a
            single, specific type of descriptive statistics (usually either
            mean or median) corresponding to groups of observations across
            columns and features across rows
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
        [index_rows, column_group,],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_extract.to_numpy())
    value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heatmap_signal_label_features_groups_of_observations(
        table=table,
        format_table=2, # 2: features in columns; observations, groups in rows
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=True,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        show_labels_ordinate=True,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        show_scale_bar=True, # whether to show scale bar on individual figures
        title_ordinate="",
        title_abscissa="",
        title_bar="gene signal (z-score)",
        size_title_ordinate="eight", # ten
        size_title_abscissa="eight", # ten
        size_label_ordinate="seventeen", # multi-panel: ten; individual: twelve
        size_label_abscissa="eleven", # multi-panel: ten; individual: twelve
        size_title_bar="twelve", # twelve
        size_label_bar="thirteen", # thirteen for whole; five for bar itself
        aspect="portrait", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )

    # Return information.
    return figure


def create_plot_chart_heatmap_mean(
    table=None,
    index_columns=None,
    index_rows=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 17 October 2024

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
    value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heatmap_signal_label_features_observations(
        table=table,
        format_table=1, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        show_labels_ordinate=True,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        show_scale_bar=True, # whether to show scale bar on individual figures
        title_ordinate="",
        title_abscissa="",
        title_bar="gene signal (z-score)",
        size_title_ordinate="eight", # ten
        size_title_abscissa="eight", # ten
        size_label_ordinate="eleven", # multi-panel: ten; individual: twelve
        size_label_abscissa="eleven", # multi-panel: ten; individual: twelve
        size_title_bar="twelve", # twelve
        size_label_bar="thirteen", # thirteen for whole; five for bar itself
        aspect="square", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )

    # Return information.
    return figure


def manage_plot_charts(
    index_features=None,
    index_observations=None,
    table_box=None,
    table_heatmap_individual_1=None,
    table_heatmap_individual_2=None,
    table_heatmap_mean=None,
    box_features=None,
    box=None,
    heatmap_individual=None,
    heatmap_mean=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        index_features (str): name for index corresponding to features
        index_observations (str): name for index corresponding to observations
        table_box (object): Pandas data-frame table of values of signal
            intensity for features across columns and sample observations
            in groups across rows
        table_heatmap_individual_1 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_individual_2 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_mean (object): Pandas data-frame table of descriptive
            statistics for values of signal intensity for features across
            rows and groups of sample observations across columns
        box_features (list<str>): identifiers of features for which to create
            box charts
        box (bool): whether to create box charts
        heatmap_individual (bool): whether to create heatmap chart for
            individual values of signal intensity
        heatmap_mean (bool): whether to create heatmap chart for means of
            signal intensity
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of figure objects from MatPlotLib

    """

    # Copy information in table.
    table_box = table_box.copy(deep=True)
    table_heatmap_individual_1 = table_heatmap_individual_1.copy(deep=True)
    table_heatmap_individual_2 = table_heatmap_individual_2.copy(deep=True)
    table_heatmap_mean = table_heatmap_mean.copy(deep=True)
    # Copy other information.
    box_features = copy.deepcopy(box_features)

    ##########
    # Box
    if (
        box and (len(box_features) < 50)
    ):
        figures_box = list()
        for feature in box_features:
            record_box = dict()
            record_box["feature"] = feature
            record_box["name"] = str("box_" + feature)
            record_box["figure"] = create_plot_chart_box(
                table=table_box,
                column_feature=feature,
                column_group="group",
                report=False,
            )
            figures_box.append(record_box)
            pass
    else:
        figures_box = None
        pass


    ##########
    # Heatmap Individual
    if heatmap_individual:
        figure_heatmap_individual_1 = create_plot_chart_heatmap_individual(
            table=table_heatmap_individual_1,
            index_columns=index_features,
            index_rows=index_observations,
            column_group="group",
            report=report,
        )
        figure_heatmap_individual_2 = create_plot_chart_heatmap_individual(
            table=table_heatmap_individual_2,
            index_columns=index_features,
            index_rows=index_observations,
            column_group="group",
            report=report,
        )
    else:
        figure_heatmap_individual_1 = None
        figure_heatmap_individual_2 = None
        pass


    ##########
    # Heatmap Mean
    if heatmap_mean:
        figure_heatmap_mean = create_plot_chart_heatmap_mean(
            table=table_heatmap_mean,
            index_columns="group_observations",
            index_rows="feature",
            report=False,
        )
    else:
        figure_heatmap_mean = None
        pass

    ##########
    # Collect information.
    pail = dict()
    pail["box"] = figures_box
    pail["heatmap_individual_1"] = figure_heatmap_individual_1
    pail["heatmap_individual_2"] = figure_heatmap_individual_2
    pail["heatmap_mean"] = figure_heatmap_mean

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.description.py")
        function = str(
            "manage_plot_charts" +
            "()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)

    # Return information.
    return pail




###############################################################################
# Procedure


##########
# Execute main procedure.


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
    procedure="scratch"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.scratch.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass



    pass


###############################################################################
# End
