"""
Studies of age, exercise, and dietary omega-3 in skeletal muscle and
subcutaneous adipose of healthy adults.

This module 'interface' is part of the 'age_exercise' package.

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

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.

import age_exercise.scratch
import age_exercise.phenotypes.organize_subject
import age_exercise.phenotypes.organize_sample
import age_exercise.phenotypes.compare_groups
import age_exercise.proteomics.organize_olink
import age_exercise.proteomics.organize_spectroscopy
import age_exercise.transcriptomics.organize_signal
import age_exercise.transcriptomics.select_gene_sets
import age_exercise.transcriptomics.compare_sets_groups
import age_exercise.transcriptomics.operate_sets

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
    parser_main = define_subparser_main(
        project="age_exercise",
        routine="main",
        subparsers=subparsers,
    )
    parser_phenotypes = define_subparser_phenotypes(
        project="age_exercise",
        routine="phenotypes",
        subparsers=subparsers,
    )
    parser_proteomics = define_subparser_proteomics(
        project="age_exercise",
        routine="proteomics",
        subparsers=subparsers,
    )
    parser_transcriptomics = define_subparser_transcriptomics(
        project="age_exercise",
        routine="transcriptomics",
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

        Welcome to the interface of the 'age_exercise' package in Python.

        This package supports studies of age, exercise, and dietary omega-3 in
        skeletal muscle and subcutaneous adipose of healthy adults. Operations
        include the wrangling, processing, organization, visual representation,
        and analysis of data from measurements of phenotypes in the clinic and
        laboratory and measurements from transcriptomics and proteomics
        technologies.

        Due to dependencies, here is the typical sequence of routines and
        procedures.

        subpackage: phenotypes
        1. phenotypes.organize_subject
        2. phenotypes.organize_sample
        3. phenotypes.compare_groups

        subpackage.transcriptomics
        _. transcriptomics.organize_signal
        _. transcriptomics.select_gene_sets
        _. ...

        subpackage: proteomics
        _. proteomics.organize_olink
        _. proteomics.organize_spectroscopy
        _. ...

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


def define_subparser_main(
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
        name=routine,
        description=description,
        epilog=epilog,
        help="Help for main routine.",
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
        "-scratch",
        "--scratch",
        dest="scratch",
        action="store_true",
        help=(
            "Test implementation of new functionality."
        )
    )

    # Define behavior.
    parser_routine.set_defaults(func=evaluate_parameters_main)
    # Return parser.
    return parser_routine


def define_subparser_phenotypes(
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
        name=routine,
        description=description,
        epilog=epilog,
        help="Help for phenotypes routine.",
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
        "-organize_subject",
        "--organize_subject",
        dest="organize_subject",
        action="store_true",
        help=(
            "Organize information about phenotypes at the level of study "
            "subjects, including measurements from the clinic and laboratory."
        )
    )
    parser_routine.add_argument(
        "-organize_sample",
        "--organize_sample",
        dest="organize_sample",
        action="store_true",
        help=(
            "Organize information about phenotypes at the level of study "
            "samples, including measurements from the clinic and laboratory."
        )
    )
    parser_routine.add_argument(
        "-compare_groups",
        "--compare_groups",
        dest="compare_groups",
        action="store_true",
        help=(
            "Compare phenotype features between groups of observations."
        )
    )

    # Define behavior.
    parser_routine.set_defaults(func=evaluate_parameters_phenotypes)
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
        name=routine,
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
        "-organize_olink",
        "--organize_olink",
        dest="organize_olink",
        action="store_true",
        help=(
            "Organize information about samples at the level of study "
            "subjects, with special attention to measurements from Olink " +
            "technology."
        )
    )
    parser_routine.add_argument(
        "-organize_spectroscopy",
        "--organize_spectroscopy",
        dest="organize_spectroscopy",
        action="store_true",
        help=(
            "Organize information about analytes (genes) and their signals " +
            "(RNA sequence read counts mapped to genes)."
        )
    )

    # Define behavior.
    parser_routine.set_defaults(func=evaluate_parameters_proteomics)
    # Return parser.
    return parser_routine


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
        name=routine,
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
        "-organize_signal",
        "--organize_signal",
        dest="organize_signal",
        action="store_true",
        help=(
            "Organize information about analytes (genes) and their signals " +
            "(RNA sequence read counts mapped to genes)."
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
    parser_routine.add_argument(
        "-compare_sets_groups",
        "--compare_sets_groups",
        dest="compare_sets_groups",
        action="store_true",
        help=(str(
            "Compare signals across genes within sets between groups of " +
            "samples."
        )),
    )
    parser_routine.add_argument(
        "-operate_sets",
        "--operate_sets",
        dest="operate_sets",
        action="store_true",
        help=(str(
            "Combine and filter sets of features by inclusion or exclusion."
        )),
    )

    # Define behavior.
    parser_routine.set_defaults(func=evaluate_parameters_transcriptomics)
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


def evaluate_parameters_main(arguments):
    """
    Evaluates parameters for routine.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to main routine ...")
    # Execute procedure.
    if arguments.scratch:
        # Report status.
        print(
           "... executing age_exercise.scratch " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.scratch.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )

    pass


def evaluate_parameters_phenotypes(arguments):
    """
    Evaluates parameters for routine.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to proteomics routine ...")
    # Execute procedure.
    if arguments.organize_subject:
        # Report status.
        print(
           "... executing age_exercise.phenotypes.organize_subject " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.phenotypes.organize_subject.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.organize_sample:
        # Report status.
        print(
           "... executing age_exercise.phenotypes.organize_sample " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.phenotypes.organize_sample.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )

    if arguments.compare_groups:
        # Report status.
        print(
           "... executing age_exercise.phenotypes.compare_groups " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.phenotypes.compare_groups.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
        pass
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
    print("... call to proteomics routine ...")
    # Execute procedure.
    if arguments.organize_olink:
        # Report status.
        print(
           "... executing age_exercise.proteomics.organize_olink " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.proteomics.organize_olink.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.organize_spectroscopy:
        # Report status.
        print(
           "... executing age_exercise.proteomics.organize_spectroscopy " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.proteomics.organize_spectroscopy.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )

    pass


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
    if arguments.organize_signal:
        # Report status.
        print(
           "... executing age_exercise.transcriptomics.organize_signal " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.transcriptomics.organize_signal.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.select_gene_sets:
        # Report status.
        print(
           "... executing age_exercise.transcriptomics.select_gene_sets " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.transcriptomics.select_gene_sets.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.compare_sets_groups:
        # Report status.
        print(
           "... executing age_exercise.transcriptomics.compare_sets_groups " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.transcriptomics.compare_sets_groups.execute_procedure(
            path_directory_dock=arguments.path_directory_dock
        )
    if arguments.operate_sets:
        # Report status.
        print(
           "... executing age_exercise.transcriptomics.operate_sets " +
           "procedure ..."
        )
        # Execute procedure.
        age_exercise.transcriptomics.operate_sets.execute_procedure(
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
