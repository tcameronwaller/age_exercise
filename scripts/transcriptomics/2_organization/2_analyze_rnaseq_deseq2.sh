#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 25 July 2024
# Review: TCW; 25 July 2024
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
# tissue: muscle
#path_directory_source="$path_directory_dock/out_exercise/transcriptomics/organization/muscle/data"
#path_directory_product="$path_directory_dock/out_exercise/transcriptomics/deseq2/muscle"
# tissue: adipose
path_directory_source="$path_directory_dock/out_exercise/transcriptomics/organization/adipose/data"
path_directory_product="$path_directory_dock/out_exercise/transcriptomics/deseq2/adipose"

# Files.
path_file_source_table_sample="${path_directory_source}/table_sample.tsv"
path_file_source_table_signal="${path_directory_source}/table_signal.tsv"
path_file_product_table="${path_directory_product}/table_result.tsv"

# Scripts.
path_script_deseq2="${path_directory_repository_partner}/scripts/r/analyze_rnaseq_deseq2.R"

# Executable handles.
path_execution_r="${path_directory_tool}/r/r-4.4.1/bin/Rscript"

# Initialize directory.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product

# Initialize file.

###############################################################################
# Organize parameters.

# Parameters.
threads=4
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error


###############################################################################
# Execute procedure.

# Execute program script in R.
$path_execution_r $path_script_deseq2 \
$path_file_source_table_sample \
$path_file_source_table_signal \
$path_file_product_table \
$threads \
$report



###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "project: ${project_main}"
  echo "technology: transcriptomics"
  echo "procedure: 2_organization"
  echo "script: 2_analyze_rnaseq_deseq2.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove directory of temporary, intermediate files.
#rm -r $path_directory_temporary



###############################################################################
# End.
