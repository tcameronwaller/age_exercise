#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 9 April 2025
# Date, last execution or modification: 9 April 2025
# Review: 9 April 2025
###############################################################################
# Note



###############################################################################
# Organize arguments.


################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"
path_directory_package_partner="$path_directory_package/partner"

path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_demonstration="$path_directory_dock/in_demonstration"
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private="$path_directory_dock/in_parameters_private"

path_directory_source_vandongen="${path_directory_demonstration}/partner/15081686_vandongen_2004"
path_directory_source_waller="${path_directory_data}/regression/age_exercise"
path_directory_product_vandongen="${path_directory_dock}/out_regression/compare_anova_regression/15081686_vandongen_2004"
path_directory_product_waller="${path_directory_dock}/out_regression/compare_anova_regression/0_waller_2025"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.
name_file_table_data_vandongen="table_data.tsv"
name_file_table_data_waller="table_sample.tsv"

# Scripts.
path_file_script_source="${path_directory_scripts}/age_exercise/record_analyses/2025-04-09_compare_anova_regression/compare_anova_regression_placebo.py"
path_file_script_product="${path_directory_package}/compare_anova_regression_placebo.py"

# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product

# Executable handles.
path_environment_main="$path_directory_tools/python/environments/main"
echo $path_environment_main

# Initialize directory.
rm -r $path_directory_product_vandongen # caution
mkdir -p $path_directory_product_vandongen
rm -r $path_directory_product_waller # caution
mkdir -p $path_directory_product_waller
#mkdir -p $path_directory_temporary

# Initialize file.

###############################################################################
# Organize parameters.

# Parameters.
threads=6
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error



###############################################################################
# Activate Python virtual environment.

# Activate Python virtual environment.
source "${path_environment_main}/bin/activate"

# Set paths for local packages and modules.
export OLD_PYTHONPATH="$PYTHONPATH"
export PYTHONPATH="$PYTHONPATH:$path_directory_package"

# Regulate concurrent or parallel process threads on node cores.
# Force Python program (especially SciPy) not to use all available cores on a
# cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "Python virtual environment: main"
  echo "path to Python installation:"
  which python3
  echo "Python path variable:"
  echo $PYTHONPATH
  sleep 1s
  echo "----------"
fi



###############################################################################
# Execute procedure.

# Execute program process in Python.
python3 $path_file_script_product \
$name_file_table_data_vandongen \
$path_directory_source_vandongen \
$path_directory_product_vandongen \
$name_file_table_data_waller \
$path_directory_source_waller \
$path_directory_product_waller \
$path_directory_dock \
$report



###############################################################################
# Deactivate Python virtual environment.

# Restore paths.
export PYTHONPATH="$OLD_PYTHONPATH"

# Deactivate Python virtual environment.
deactivate
#which python3

# Remove directory of temporary, intermediate files.
#rm -r $path_directory_temporary

###############################################################################
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: demonstration_compare_anova_regression_placebo.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  echo "----------"
  echo "----------"
fi

###############################################################################
# End.
