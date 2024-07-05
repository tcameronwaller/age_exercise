#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 5 July 2024
# Date, last execution or modification: 5 July 2024
# Review: TCW; 5 July 2024
###############################################################################
# Note

# This Bash script restores the version of the exercise repository, including
# its parameters, scripts, package, and subpackages for remote execution,
# meaning execution on a remote server rather than a local machine. This
# restoration pulls code from repositories on GitHub and organizes this code in
# a more temporary location for convenient process execution.

# It is unnecessary to copy repositories to the working directory for
# execution; however, it is necessary for the parent directory of each Python
# package to have the appropriate name. In order to store Python code within
# the "package" subdirectories of repositories, it is convenient to copy and
# rename the package directories.

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_parent_project=$(<"./paths/endocrinology/parent_tcameronwaller.txt")
path_directory_process="${path_directory_parent_project}/process"
path_directory_dock="${path_directory_process}/dock"
#path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_repository_partner="${path_directory_process}/partner"
path_directory_repository_exercise="${path_directory_process}/exercise"
path_directory_package="$path_directory_process/package"
path_directory_package_partner_source="${path_directory_repository_partner}/package"
path_directory_package_exercise_source="${path_directory_repository_exercise}/package"
path_directory_package_partner_product="$path_directory_package/partner"
path_directory_package_exercise_product="$path_directory_package/exercise"

# Initialize directories.
if [ -d $path_directory_parameters ] || [ -d $path_directory_repository_partner ] || [ -d $path_directory_repository_exercise ] || [ -d $path_directory_package ]; then
  # Remove previous versions of code from temporary location for execution.
  rm -rf $path_directory_parameters
  rm -rf $path_directory_repository_partner
  rm -rf $path_directory_repository_exercise
  rm -rf $path_directory_package
fi

if [ ! -d $path_directory_process ] || [ ! -d $path_directory_package ] || [ ! -d $path_directory_dock ] || [ ! -d $path_directory_parameters ] ; then
  # Directory or directories do not already exist.
  # Create directories.
  mkdir -p $path_directory_process
  mkdir -p $path_directory_package
  mkdir -p $path_directory_dock
  mkdir -p $path_directory_parameters
fi

################################################################################
# Execute procedure.

##########
# Access and organize current versions of repositories.
# Repository: partner
# Scripts remain within original repository's structure.
# Parameters transfer to a common parent directory.
# Python code transfers to a common parent directory.
cd $path_directory_process
wget https://github.com/tcameronwaller/partner/archive/main.zip
unzip main.zip
rm main.zip
mv partner-main $path_directory_repository_partner
# Repository: exercise
# Scripts remain within original repository's structure.
# Parameters transfer to a common parent directory.
# Python code transfers to a common parent directory.
wget https://github.com/tcameronwaller/exercise/archive/main.zip
unzip main.zip
rm main.zip
mv exercise-main $path_directory_repository_exercise

##########
# Organize parameters.
cp -r "$path_directory_repository_partner/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/partner"
cp -r "$path_directory_repository_exercise/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/exercise"

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

################################################################################
# Report.
echo "----------"
echo "Script complete:"
echo $0 # Print full file path to script.
echo "restore_repository_remote.sh"
echo "----------"

################################################################################
# End.
