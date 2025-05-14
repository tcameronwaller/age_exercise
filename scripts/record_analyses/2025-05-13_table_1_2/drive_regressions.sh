#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 28 March 2025
# Date, last execution or modification: 17 April 2025
# Review: 17 April 2025
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

#path_directory_source="${path_directory_dock}/${demonstration}/source"
path_directory_product="${path_directory_dock}/out_regression/age_exercise/table_1_results"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.
#path_file_table_parameters="${path_directory_demonstration}/partner/table_regression_parameters.tsv"
#path_file_table_parameters="${path_directory_parameters_private}/age_exercise/regression/table_parameters_regression.tsv"

#path_file_table_parameters="${path_directory_parameters_private}/age_exercise/regression/table_parameters_regression_automatic_age.tsv"
#path_file_table_parameters="${path_directory_parameters_private}/age_exercise/regression/table_parameters_regression_automatic_omega3.tsv"
path_file_table_parameters="${path_directory_parameters_private}/age_exercise/regression/table_parameters_regression_automatic_omega3_mix.tsv"

#path_file_table_parameters="${path_directory_parameters_private}/age_exercise/regression/table_parameters_regression_automatic_omega3_female.tsv"
#path_file_table_parameters="${path_directory_parameters_private}/age_exercise/regression/table_parameters_regression_automatic_omega3_male.tsv"

# Scripts.
path_file_script_source="${path_directory_scripts}/partner/python/drive_regressions_from_table_parameters.py"
path_file_script_product="${path_directory_package}/drive_regressions_from_table_parameters.py"

# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product

# Executable handles.
path_environment_main="$path_directory_tools/python/environments/main"
echo $path_environment_main

# Initialize directory.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
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
$path_file_table_parameters \
$path_directory_product \
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
  echo "script: template_drive_regressions.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  echo "Convert identifiers or names of genes by query to MyGene.info"
  echo "path to file for table of parameters: " $path_file_table_parameters
  echo "path to dock directory: " $path_directory_dock
  echo "----------"
  echo "----------"
  echo "----------"
fi

###############################################################################
# End.
