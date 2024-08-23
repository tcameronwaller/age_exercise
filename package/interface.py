"""
Manage execution of procedures.

This module 'interface' is part of the 'exercise' package.

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

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.

import exercise.transcriptomics.organize_sample
import exercise.transcriptomics.organize_signal
import exercise.transcriptomics.select_gene_sets
import exercise.proteomics.organize_sample_olink

#dir()
#importlib.reload()

###############################################################################
# Functionality


def define_interface_parsers():
    """
    Defines and parses arguments from terminal's interface.

    arguments:

    raises:

    returns:
        (object): arguments from terminal

    """

    # Define description.
    description = define_main_description()
    # Define epilog.
    epilog = define_main_epilog()
    # Define arguments.
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(title="procedures")
    #parser_main = define_main_subparser(subparsers=subparsers)
    parser_transcriptomics = define_subparser_transcriptomics(
        project="exercise",
        routine="transcriptomics",
        subparsers=subparsers,
    )
    parser_proteomics = define_subparser_proteomics(
        project="exercise",
        routine="proteomics",
        subparsers=subparsers,
    )
    # Parse arguments.
    return parser.parse_args()


def define_main_description():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
        Description of main project package

        Welcome to the interface of the "exercise" package in Python.
        This package supports project-specific wrangling, processing,
        and analysis of data from transcriptomics and proteomics
        technologies.

        --------------------------------------------------
    """)
    return description


def define_main_epilog():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\
        --------------------------------------------------
        Epilog of main project package

        Well, I guess we did some fun stuff with data.
        Sure was nice visiting. Hope y'all come back soon!
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_subparser_transcriptomics(
    project=None,
    routine=None,
    subparsers=None,
):
    """
    Defines subparser for a specific routine of procedures.

    arguments:
        project (str): name of project
        routine (str): name of routine
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    # Define description.
    description = define_routine_description()
    # Define epilog.
    epilog = define_routine_epilog()
    # Define parser.
    parser_routine = subparsers.add_parser(
        name="transcriptomics",
        description=description,
        epilog=epilog,
        help="Help for transcriptomics routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser_routine.add_argument(
        "-path_directory_dock", "--path_directory_dock",
        dest="path_directory_dock", type=str, required=True,
        help=(
            "Path to dock directory for source and product " +
            "directories and files."
        )
    )
    parser_routine.add_argument(
        "-organize_sample",
        "--organize_sample",
        dest="organize_sample",
        action="store_true",
        help=(
            "Organize information about samples."
        )
    )
    parser_routine.add_argument(
        "-organize_signal",
        "--organize_signal",
        dest="organize_signal",
        action="store_true",
        help=(
            "Organize information about signals."
        )
    )
    parser_routine.add_argument(
        "-select_gene_sets",
        "--select_gene_sets",
        dest="select_gene_sets",
        action="store_true",
        help=(
            "Analyze information about genes with differential expression."
        )
    )

    # Define behavior.
    parser_routine.set_defaults(func=evaluate_parameters_transcriptomics)
    # Return parser.
    return parser_routine


def define_subparser_proteomics(
    project=None,
    routine=None,
    subparsers=None,
):
    """
    Defines subparser for a specific routine of procedures.

    arguments:
        project (str): name of project
        routine (str): name of routine
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    # Define description.
    description = define_routine_description()
    # Define epilog.
    epilog = define_routine_epilog()
    # Define parser.
    parser_routine = subparsers.add_parser(
        name="proteomics",
        description=description,
        epilog=epilog,
        help="Help for proteomics routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser_routine.add_argument(
        "-path_directory_dock", "--path_directory_dock",
        dest="path_directory_dock", type=str, required=True,
        help=(
            "Path to dock directory for source and product " +
            "directories and files."
        )
    )
    parser_routine.add_argument(
        "-organize_sample_olink",
        "--organize_sample_olink",
        dest="organize_sample",
        action="store_true",
        help=(
            "Organize information about samples, including measurements by " +
            "O-Link technology."
        )
    )

    # Define behavior.
    parser_routine.set_defaults(func=evaluate_parameters_proteomics)
    # Return parser.
    return parser_routine


def define_routine_description(
    project=None,
    routine=None,
):
    """
    Defines description for terminal interface.

    arguments:
        project (str): name of project
        routine (str): name of routine

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
        Description of routine within main project package

        project: {project}
        routine: {routine}

        --------------------------------------------------
    """).format(
        project=project,
        routine=routine,
    )
    return description


def define_routine_epilog(
    project=None,
    routine=None,
):
    """
    Defines epilog for terminal interface.

    arguments:
        project (str): name of project
        routine (str): name of routine

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\
        --------------------------------------------------
        Epilog of routine within main project package

        project: {project}
        routine: {routine}

        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """).format(
        project=project,
        routine=routine,
    )
    return epilog


def evaluate_parameters_transcriptomics(arguments):
    """
    Evaluates parameters for routine.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to transcriptomics routine ...")
    # Execute procedure.
    if arguments.organize_sample:
        # Report status.
        print(
           "... executing exercise.transcriptomics.organize_sample " +
           "procedure ..."
        )
        # Execute procedure.
        exercise.transcriptomics.organize_sample.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.organize_signal:
        # Report status.
        print(
           "... executing exercise.transcriptomics.organize_signal " +
           "procedure ..."
        )
        # Execute procedure.
        exercise.transcriptomics.organize_signal.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.select_gene_sets:
        # Report status.
        print(
           "... executing exercise.transcriptomics.select_gene_sets " +
           "procedure ..."
        )
        # Execute procedure.
        exercise.transcriptomics.select_gene_sets.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )

    pass


def evaluate_parameters_proteomics(arguments):
    """
    Evaluates parameters for routine.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to transcriptomics routine ...")
    # Execute procedure.
    if arguments.organize_sample_olink:
        # Report status.
        print(
           "... executing exercise.proteomics.organize_sample_olink " +
           "procedure ..."
        )
        # Execute procedure.
        exercise.proteomics.organize_sample_olink.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )

    pass




###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    # Parse arguments from terminal.
    arguments = define_interface_parsers()
    # Call the appropriate function.
    arguments.func(arguments)


if (__name__ == "__main__"):
    execute_procedure()
