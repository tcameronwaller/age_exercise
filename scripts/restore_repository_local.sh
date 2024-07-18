#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 July 2024
# Date, last execution or modification: 17 July 2024
# Review: TCW; 17 July 2024
###############################################################################
# Note

# This Bash script restores the version of the exercise repository, including
# its parameters, scripts, package, and subpackages for local execution,
# meaning execution on a local machine rather than a remote server. This
# restoration does not pull code from repositories on GitHub; rather, it copies
# code from a local, more permanent storage location to a local, more temporary
# location for convenient process execution.

# It is unnecessary to copy repositories to the working directory for
# execution; however, it is necessary for the parent directory of each Python
# package to have the appropriate name. In order to store Python code within
# the "package" subdirectories of repositories, it is convenient to copy and
# rename the package directories.

# It is possible for the paths specific to a Python virtual environment to
# become corrupted. It might be necessary to examine the paths in the
# environment's "activate" script. In particular, any change to the names of
# directories in the file system's path to the virtual environment will disrupt
# the environment since the "activate" script uses absolute paths. The files
# for the relevant scripts have names "activate", "activate.csh", and
# "activate.fish". If correcting the paths in these activation scripts is not
# sufficient, then it might be necessary to re-create the environment.
# https://stackoverflow.com/questions/55740123/pip-virtualenv-reset-the-path-after-reactivating

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private_source=$(<"$path_directory_paths/path_parameters_private_exercise.txt")
path_directory_parameters_private="$path_directory_dock/in_parameters_private"
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")
path_directory_repository_exercise=$(<"$path_directory_paths/path_directory_repository_exercise.txt")
path_directory_package="$path_directory_process/package"
path_directory_package_partner_source="$path_directory_repository_partner/package"
path_directory_package_exercise_source="$path_directory_repository_exercise/package"
path_directory_package_partner_product="$path_directory_package/partner"
path_directory_package_exercise_product="$path_directory_package/exercise"

# Initialize directories.
if [ -d $path_directory_parameters ] || [ -d $path_directory_package ] || [ -d $path_directory_parameters_private ] ; then
  # Remove previous versions of code from temporary location for execution.
  rm -rf $path_directory_parameters
  rm -rf $path_directory_parameters_private
  rm -rf $path_directory_package
fi

if [ ! -d $path_directory_process ] || [ ! -d $path_directory_package ] || [ ! -d $path_directory_dock ] || [ ! -d $path_directory_parameters ] || [ ! -d $path_directory_parameters_private ] ; then
  # Directory or directories do not already exist.
  # Create directories.
  mkdir -p $path_directory_process
  mkdir -p $path_directory_package
  mkdir -p $path_directory_dock
  mkdir -p $path_directory_parameters
  mkdir -p $path_directory_parameters_private
fi

################################################################################
# Execute procedure.

##########
# Organize parameters.
cp -r "$path_directory_repository_partner/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/partner"
cp -r "$path_directory_repository_exercise/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/exercise"
cp -r $path_directory_parameters_private_source $path_directory_dock
mv "${path_directory_dock}/parameters" "${path_directory_dock}/exercise"
mv "${path_directory_dock}/exercise" $path_directory_parameters_private

##########
# Organize Python packages.
# Package: "partner"
# Hierarchy: child, lower-level subpackage
# Scripts remain within original repository's structure.
# Python code transfers to a subpackage child directory within the parent
# directory of the main package.
if true; then
  cp -r "$path_directory_package_partner_source" "$path_directory_package"
  mv "$path_directory_package/package" "$path_directory_package_partner_product"
fi
# Package: "exercise"
# Hierarchy: child, lower-level subpackage
# Scripts remain within original repository's structure.
# Python code transfers to a subpackage child directory within the parent
# directory of the main package.
if true; then
  cp -r "$path_directory_package_exercise_source" "$path_directory_package"
  mv "$path_directory_package/package" "$path_directory_package_exercise_product"
fi

##########
# Initialize directory permission.
chmod -R 0777 $path_directory_process

###############################################################################
# Report.
echo "----------"
echo "project: exercise"
echo "script: restore_repository_local.sh"
echo $0 # Print full file path to script.
echo "done"
echo "----------"

###############################################################################
# End.
