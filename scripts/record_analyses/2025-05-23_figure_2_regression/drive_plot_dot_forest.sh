#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 May 2025
# Date, last execution or modification: 23 May 2025
# Review: 23 May 2025
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

#path_file_table_data="${path_directory_source}/table_results_regression_placebo_omega3_ols_no_scale.tsv"
#path_file_table_data="${path_directory_source}/table_results_regression_placebo_omega3_ols_yes_scale.tsv"
#path_file_table_data="${path_directory_source}/table_results_regression_placebo_omega3_mix_no_scale.tsv"
path_file_table_data="${path_directory_source}/table_results_regression_placebo_omega3_mix_yes_scale.tsv"

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

#type="ols"
type="mix"

# For "title", "legend_series_primary", and "legend_series_secondary", replace
# character "#" with white space.

# Format of parameters for names of columns.
# name_product: name_source

# Format of parameters for names of features.
# name_source: name_product

##########
# Ordinary Least Squares (OLS) Linear Regression
if [ "$type" == "ols" ]; then
  title="Placebo#Omega-3;#OLS" # cannot accommodate white space
  feature="feature:feature_response"
  #feature="feature:name_combination"
  #feature="feature:name"
  #features="none"
  # all features from regression
  features="age_log,body_mass_index_log,body_fat_percent_log,body_skeletal_muscle_index_log,activity_moderate_vigorous_log,activity_steps_log,oxygen_consumption_log,temperature_log,heart_rate_log,pressure_blood_systolic_log,pressure_blood_diastolic_log,red_blood_cells_log,hemoglobin_log,hematocrit_log,mean_corpuscular_volume_log,rbc_distribution_width_log,erythrocyte_sedimentation_rate_log,platelets_log,prothrombin_time_log,blood_clot_inr_log,white_blood_cells_log,neutrophils_log,lymphocytes_log,monocytes_log,eosinophils_log,basophils_log,glucose_log,insulin_log,insulin_sensitivity_log,homa_insulin_resist_log,alanine_transaminase_log,aspartate_transaminase_log,ratio_ast_alt_log,thyroid_stimulate_hormone_log,c_react_protein_log,omega3_eicosapentaenoate_log,omega3_docosahexaenoate_log,triglyceride_log,cholesterol_log,lipoprotein_hdl_log,lipoprotein_nonhdl_log,lipoprotein_ldl_log,adipocyte_lipid_content_log,cd68_adipose_percent_log,cd14_adipose_percent_log,cd206_adipose_percent_log,p16_adipose_percent_log,mitochondrial_respiration_maximum_log"
  # selection of features from regression
  #features=""
  translation_features="none"
  legend_series_primary="omega-3"
  legend_series_secondary="placebo" # or "none"
  legend_series_tertiary="none" # or "none"
  values_intervals_primary="value_primary:predictor_1_parameter;interval_low_primary:predictor_1_interval_95;interval_high_primary:predictor_1_interval_95_copy"
  #values_intervals_secondary="none"
  values_intervals_secondary="value_secondary:predictor_2_parameter;interval_low_secondary:predictor_2_interval_95;interval_high_secondary:predictor_2_interval_95_copy"
  values_intervals_tertiary="none"
  #values_intervals_tertiary="value_tertiary:predictor_3_parameter;interval_low_tertiary:predictor_3_interval_95;interval_high_tertiary:predictor_3_interval_95_copy"
fi

##########
# Mixed Effects Linear Regression
if [ "$type" == "mix" ]; then
  title="Placebo#Omega-3;#Mixed#Effects" # cannot accommodate white space
  feature="feature:feature_response"
  #feature="feature:name_combination"
  #feature="feature:name"
  #features="none"
  # all features from regression
  features="age_log,body_mass_index_log,body_fat_percent_log,body_skeletal_muscle_index_log,activity_moderate_vigorous_log,activity_steps_log,oxygen_consumption_log,temperature_log,heart_rate_log,pressure_blood_systolic_log,pressure_blood_diastolic_log,red_blood_cells_log,hemoglobin_log,hematocrit_log,mean_corpuscular_volume_log,rbc_distribution_width_log,erythrocyte_sedimentation_rate_log,platelets_log,prothrombin_time_log,blood_clot_inr_log,white_blood_cells_log,neutrophils_log,lymphocytes_log,monocytes_log,eosinophils_log,basophils_log,glucose_log,insulin_log,insulin_sensitivity_log,homa_insulin_resist_log,alanine_transaminase_log,aspartate_transaminase_log,ratio_ast_alt_log,thyroid_stimulate_hormone_log,c_react_protein_log,omega3_eicosapentaenoate_log,omega3_docosahexaenoate_log,triglyceride_log,cholesterol_log,lipoprotein_hdl_log,lipoprotein_nonhdl_log,lipoprotein_ldl_log,adipocyte_lipid_content_log,cd68_adipose_percent_log,cd14_adipose_percent_log,cd206_adipose_percent_log,p16_adipose_percent_log,mitochondrial_respiration_maximum_log"
  # selection of features from regression
  #features=""
  translation_features="none"
  legend_series_primary="omega-3"
  legend_series_secondary="placebo" # or "none"
  legend_series_tertiary="none" # or "none"
  values_intervals_primary="value_primary:predictor_0_parameter;interval_low_primary:predictor_0_interval_95;interval_high_primary:predictor_0_interval_95_copy"
  #values_intervals_secondary="none"
  values_intervals_secondary="value_secondary:predictor_1_parameter;interval_low_secondary:predictor_1_interval_95;interval_high_secondary:predictor_1_interval_95_copy"
  values_intervals_tertiary="none"
  #values_intervals_tertiary="value_tertiary:predictor_3_parameter;interval_low_tertiary:predictor_3_interval_95;interval_high_tertiary:predictor_3_interval_95_copy"
fi

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
