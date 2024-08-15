#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 15 August 2024
# Review: TCW; 15 August 2024
###############################################################################
# Note



###############################################################################
# Organize arguments.



###############################################################################
# Organize paths.

# Project.
project_main="exercise"

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tool=$(<"$path_directory_paths/path_directory_tool.txt")
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_dock="$path_directory_process/dock"
path_directory_parameters_private="$path_directory_dock/in_parameters_private/exercise/transcriptomics"

# Files.
path_file_table_parameter="$path_directory_parameters_private/table_set_differential_expression.tsv"

# Scripts.
path_script_deseq2="${path_directory_repository_partner}/scripts/r/analyze_rnaseq_deseq2.R"

# Executable handles.
path_execution_r="${path_directory_tool}/r/r-4.4.1/bin/Rscript"

# Initialize directory.

# Initialize file.

###############################################################################
# Organize parameters.

# Parameters.
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error
threads="5"
report="true"

###############################################################################
# Drive procedure across instances with parameters from table.

##########
# Iteratively read lines from file, split fields within each line by space,
# tab, or new-line delimiters, extract parameters, and execute procedure with
# those parameters.
input=$path_file_table_parameter
while IFS=$' \t\n' read -r -a array
do
  # Extract values from individual columns within table's current row.
  raw_inclusion="${array[0]}"
  raw_name_set="${array[1]}"
  raw_tissue="${array[2]}"
  raw_cohort_selection="${array[3]}"
  raw_factor_availability="${array[4]}"
  raw_formula_text="${array[5]}"
  raw_condition="${array[6]}"
  raw_levels_condition="${array[7]}"
  raw_supplement_1="${array[8]}"
  raw_levels_supplement_1="${array[9]}"
  raw_supplement_2="${array[10]}"
  raw_levels_supplement_2="${array[11]}"
  raw_subject="${array[12]}"
  raw_threshold_significance="${array[13]}"
  raw_note="${array[14]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "--------------------------------------------------"
    echo "field 0, inclusion: ${raw_inclusion}"
    echo "field 1, name_set: ${raw_name_set}"
    echo "field 2, tissue: ${raw_tissue}"
    echo "field 3, cohort_selection: ${raw_cohort_selection}"
    echo "field 4, factor_availability: ${raw_factor_availability}"
    echo "field 5, formula_text: ${raw_formula_text}"
    echo "field 6, condition: ${raw_condition}"
    echo "field 7, levels_condition: ${raw_levels_condition}"
    echo "field 8, supplement_1: ${raw_supplement_1}"
    echo "field 9, levels_supplement_1: ${raw_levels_supplement_1}"
    echo "field 10, supplement_2: ${raw_supplement_2}"
    echo "field 11, levels_supplement_2: ${raw_levels_supplement_2}"
    echo "field 12, subject: ${raw_subject}"
    echo "field 13, threshold_significance: ${raw_threshold_significance}"
    echo "field 14, note: ${raw_note}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  #  && [ "$raw_name_set" == "muscle_exercise-0hr_sex" ]
  if [ "$raw_inclusion" == "1" ]; then

    ##########
    # Organize paths for current instance.
    # Directories.
    path_directory_source="$path_directory_dock/out_exercise/transcriptomics/organize_signal/${raw_tissue}/${raw_name_set}/data"
    path_directory_product="$path_directory_dock/out_exercise/transcriptomics/deseq2/${raw_tissue}/${raw_name_set}"
    # Files.
    path_file_source_table_sample="${path_directory_source}/table_sample.tsv"
    path_file_source_table_gene="${path_directory_source}/table_gene.tsv"
    path_file_source_table_signal="${path_directory_source}/table_signal.tsv"
    path_file_product_table="${path_directory_product}/table_result_deseq2.tsv"
    # Initialize directory.
    rm -r $path_directory_product # caution
    mkdir -p $path_directory_product

    # Report.
    if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
      echo "----------"
      echo "path_file_source_table_sample: $path_file_source_table_sample"
      echo "path_file_source_table_gene: $path_file_source_table_gene"
      echo "path_file_source_table_signal: $path_file_source_table_signal"
      echo "path_file_product_table: $path_file_product_table"
      echo "--------------------------------------------------"
    fi

    ##########
    # Execute program script in R.
    if true; then
      $path_execution_r $path_script_deseq2 \
      $path_file_source_table_sample \
      $path_file_source_table_gene \
      $path_file_source_table_signal \
      $path_file_product_table \
      $raw_formula_text \
      $raw_condition \
      $raw_levels_condition \
      $raw_supplement_1 \
      $raw_levels_supplement_1 \
      $raw_supplement_2 \
      $raw_levels_supplement_2 \
      $raw_subject \
      $raw_threshold_significance \
      $threads \
      $report
    fi

  fi
done < "${input}"



###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "project: ${project_main}"
  echo "technology: transcriptomics"
  echo "procedure: 3_analysis"
  echo "script: 1_analyze_rnaseq_deseq2.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi



###############################################################################
# End.
