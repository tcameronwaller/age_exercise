#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 11 November 2024
# Review: TCW; 11 November 2024
###############################################################################
# Note



###############################################################################
# Organize arguments.



###############################################################################
# Organize paths.

# Project.
project_main="age_exercise"

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_dock="$path_directory_process/dock"
path_directory_parameters_private="$path_directory_dock/in_parameters_private/age_exercise/transcriptomics"
path_directory_product_parent="$path_directory_dock/out_age_exercise/transcriptomics/deseq2"

# Files.
path_file_table_parameter="$path_directory_parameters_private/table_differential_expressions_genes_samples.tsv"

# Scripts.
path_script_deseq2="${path_directory_repository_partner}/scripts/r/analyze_rnaseq_deseq2.R"

# Executable handles.
path_execution_r="${path_directory_tools}/r/r-4.4.1/bin/Rscript"

# Initialize directory.
rm -r $path_directory_product_parent # caution
mkdir -p $path_directory_product_parent


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
  raw_tissue="${array[1]}"
  raw_sort="${array[2]}"
  raw_group="${array[3]}"
  raw_instance="${array[4]}"
  raw_selection_samples_primary="${array[5]}"
  raw_selection_samples_secondary="${array[6]}"
  raw_continuity_scale="${array[7]}"
  raw_selection_genes="${array[8]}"
  raw_formula_text="${array[9]}"
  raw_condition="${array[10]}"
  raw_levels_condition="${array[11]}"
  raw_supplement_1="${array[12]}"
  raw_levels_supplement_1="${array[13]}"
  raw_supplement_2="${array[14]}"
  raw_levels_supplement_2="${array[15]}"
  raw_supplement_3="${array[16]}"
  raw_levels_supplement_3="${array[17]}"
  raw_subject="${array[18]}"
  raw_results_contrast="${array[19]}"
  raw_results_name="${array[20]}"
  raw_threshold_significance="${array[21]}"
  raw_name_set_gene_emphasis="${array[22]}"
  raw_name_set_gene_exclusion="${array[23]}"
  raw_review="${array[24]}"
  raw_note="${array[25]}"
  name_instance="${raw_tissue}_${raw_sort}_${raw_instance}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "--------------------------------------------------"
    echo "raw fields"
    echo "--------------------------------------------------"
    echo "field 0, inclusion: ${raw_inclusion}"
    echo "field 1, tissue: ${raw_tissue}"
    echo "field 2, sort: ${raw_sort}"
    echo "field 3, group: ${raw_group}"
    echo "field 4, instance: ${raw_instance}"
    echo "field 5, selection_samples_primary: ${raw_selection_samples_primary}"
    echo "field 6, selection_samples_secondary: ${raw_selection_samples_secondary}"
    echo "field 7, continuity_scale: ${raw_continuity_scale}"
    echo "field 8, selection_genes: ${raw_selection_genes}"
    echo "field 9, formula_text: ${raw_formula_text}"
    echo "field 10, condition: ${raw_condition}"
    echo "field 11, levels_condition: ${raw_levels_condition}"
    echo "field 12, supplement_1: ${raw_supplement_1}"
    echo "field 13, levels_supplement_1: ${raw_levels_supplement_1}"
    echo "field 14, supplement_2: ${raw_supplement_2}"
    echo "field 15, levels_supplement_2: ${raw_levels_supplement_2}"
    echo "field 16, supplement_3: ${raw_supplement_3}"
    echo "field 17, levels_supplement_3: ${raw_levels_supplement_3}"
    echo "field 18, subject: ${raw_subject}"
    echo "field 19, results_contrast: ${raw_results_contrast}"
    echo "field 20, results_name: ${raw_results_name}"
    echo "field 21, threshold_significance: ${raw_threshold_significance}"
    echo "field 22, name_set_gene_emphasis: ${raw_name_set_gene_emphasis}"
    echo "field 23, name_set_gene_exclusion: ${raw_name_set_gene_exclusion}"
    echo "field 24, review: ${raw_review}"
    echo "field 25, note: ${raw_note}"
    echo "----------"
    echo "derivation fields"
    echo "----------"
    echo "name_instance: ${name_instance}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  #  && [ "$raw_name_instance" == "muscle_exercise-0hr_sex" ]
  if [ "$raw_inclusion" == "1" ]; then

    ##########
    # Organize paths for current instance.
    # Directories.
    path_directory_source="$path_directory_dock/out_age_exercise/transcriptomics/organize_signal/parts/${raw_tissue}/${name_instance}/data"
    path_directory_product="$path_directory_product_parent/${raw_tissue}/${raw_group}"
    # Files.
    path_file_source_table_sample="${path_directory_source}/table_sample.tsv"
    path_file_source_table_gene="${path_directory_source}/table_gene.tsv"
    path_file_source_table_signal="${path_directory_source}/table_signal.tsv"
    path_file_product_table="${path_directory_product}/table_result_deseq2_${name_instance}.tsv"
    # Initialize directory.
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
      $raw_supplement_3 \
      $raw_levels_supplement_3 \
      $raw_subject \
      $raw_results_contrast \
      $raw_results_name \
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
