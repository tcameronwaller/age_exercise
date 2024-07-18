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



###############################################################################
# Organize arguments.



###############################################################################
# Organize paths.



# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tool=$(<"$path_directory_paths/path_directory_tool.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_package="$path_directory_process/package"
path_directory_package_partner="$path_directory_package/partner"
path_directory_package_exercise="$path_directory_package/exercise"
path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_parameters="$path_directory_dock/in_parameters"
#path_directory_source="${path_directory_dock}/source"
#path_directory_product="${path_directory_dock}/product"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.

# Scripts.

# Executable handles.
path_environment_main="$path_directory_tool/python/environments/main"
echo $path_environment_main
# Initialize directory.
#rm -r $path_directory_product # caution
#mkdir -p $path_directory_product
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
#echo "Python path variable before update"
#echo $PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$path_directory_package
export PYTHONPATH=$PYTHONPATH:$path_directory_package_partner
export PYTHONPATH=$PYTHONPATH:$path_directory_package_exercise
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
python3 $path_directory_package_exercise/interface.py \
main \
--exercise_transcriptomics_organization \
--path_directory_dock $path_directory_dock

###############################################################################
# Deactivate Python virtual environment.

# Deactivate Python virtual environment.
deactivate
#which python3



###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "project: exercise"
  echo "technology: transcriptomics"
  echo "procedure: 2_organization"
  echo "script: 1_filter_fill_scale_normalize.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove directory of temporary, intermediate files.
#rm -r $path_directory_temporary



###############################################################################
# End.
