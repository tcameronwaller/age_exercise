#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 May 2025
# Date, last execution or modification: 16 May 2025
# Review: 16 May 2025
###############################################################################
# Note

# TODO: TCW; 9 May 2025
# pass a new parameter for the name of the file for the plot chart


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


path_directory_source="${path_directory_parameters_private}/age_exercise/regression"
path_directory_product="${path_directory_dock}/out_regression/age_exercise/figure_2"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.

path_file_table_data="${path_directory_source}/table_results_regression_placebo_omega3.tsv"

# Scripts.
path_file_script_source="${path_directory_scripts}/partner/python/drive_plot_dot_forest_from_table_data.py"
path_file_script_product="${path_directory_package}/drive_plot_dot_forest_from_table_data.py"

# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product

# Executable handles.
path_environment_main="$path_directory_tools/python/environments/main"
echo $path_environment_main

# Initialize directory.
#rm -r $path_directory_product # caution
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

# For "title", "legend_series_primary", and "legend_series_secondary", replace
# character "#" with white space.

# Format of parameters for names of columns.
# name_product: name_source

# Format of parameters for names of features.
# name_source: name_product

title="Effect#Placebo#Omega-3" # cannot accommodate white space
feature="feature:feature_response"
#feature="feature:name_combination"
#feature="feature:name"
#features="none"
features="age,body_mass_index,body_fat_percent,body_skeletal_muscle_index,activity_moderate_vigorous,activity_steps,oxygen_consumption,temperature,heart_rate,pressure_blood_systolic,pressure_blood_diastolic,glucose,insulin,insulin_sensitivity,homa_insulin_resist,alanine_transaminase,aspartate_transaminase,red_blood_cells,hemoglobin,hematocrit,mean_corpuscular_volume,rbc_distribution_width,platelets,prothrombin_time,blood_clot_inr,omega3_eicosapentaenoate,omega3_docosahexaenoate,triglyceride,cholesterol,lipoprotein_hdl,lipoprotein_nonhdl,lipoprotein_ldl,adipocyte_lipid_content,cd68_adipose_percent,cd14_adipose_percent,cd206_adipose_percent,p16_adipose_percent,c_react_protein,erythrocyte_sedimentation_rate,white_blood_cells,neutrophils,lymphocytes,monocytes,eosinophils,basophils,thyroid_stimulate_hormone,mitochondrial_respiration_maximum"
translation_features="none"
legend_series_primary="omega-3"
legend_series_secondary="placebo" # or "none"
legend_series_tertiary="none" # or "none"
values_intervals_primary="value_primary:predictor_0_parameter;interval_low_primary:predictor_0_interval_95;interval_high_primary:predictor_0_interval_95_copy"
#values_intervals_secondary="none"
values_intervals_secondary="value_secondary:predictor_1_parameter;interval_low_secondary:predictor_1_interval_95;interval_high_secondary:predictor_1_interval_95_copy"
values_intervals_tertiary="none"
#values_intervals_tertiary="value_tertiary:predictor_3_parameter;interval_low_tertiary:predictor_3_interval_95;interval_high_tertiary:predictor_3_interval_95_copy"

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
$path_file_table_data \
$title \
$feature \
$features \
$translation_features \
$legend_series_primary \
$legend_series_secondary \
$legend_series_tertiary \
$values_intervals_primary \
$values_intervals_secondary \
$values_intervals_tertiary \
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
  echo "script: template_drive_plot_dot_forest.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  echo "path to file for table of data: " $path_file_table_data
  echo "path to dock directory: " $path_directory_dock
  echo "----------"
  echo "----------"
  echo "----------"
fi

###############################################################################
# End.
