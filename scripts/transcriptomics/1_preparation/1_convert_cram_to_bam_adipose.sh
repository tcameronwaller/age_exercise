#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 5 July 2024
# Date, last execution or modification: 10 July 2024
# Review: TCW; 10 July 2024
###############################################################################
# Note

# View header of file in BAM or CRAM format to discern details of its prior
# processing.
# samtools head --headers 100 /path/to/input/file.bam
# samtools view --header-only --output /path/to/output/header/file.txt /path/to/input/file.bam

###############################################################################


###############################################################################
# Organize arguments.

###############################################################################
# Organize paths.

# Directories.
cd ~
path_directory_parent_project=$(<"./paths/endocrinology/parent_tcameronwaller.txt")
path_directory_reference="${path_directory_parent_project}/reference"
path_directory_tool="${path_directory_parent_project}/tool"
path_directory_process="${path_directory_parent_project}/process"
path_directory_dock="${path_directory_process}/dock"

#path_directory_source=$("./paths/endocrinology/transcriptomics_adipose.txt")
#path_directory_source=$("./paths/endocrinology/transcriptomics_muscle.txt")
path_directory_source="${path_directory_dock}/test_lanza_rnaseq_adipose_2024/raw"
#path_directory_source="${path_directory_dock}/test_lanza_rnaseq_muscle_2022/raw"
path_directory_product="${path_directory_dock}/test_lanza_rnaseq_adipose_2024/bam"
path_directory_parallel="${path_directory_product}/parallel"

# Executable handles.
path_execution_samtools="${path_directory_tool}/samtools-1.20/bin/samtools"

# Scripts.
path_script_convert_cram_to_bam_1="${path_directory_process}/partner/scripts/samtools/convert_cram_to_bam_slurm_1.sh"
path_script_convert_cram_to_bam_2="${path_directory_process}/partner/scripts/samtools/convert_cram_to_bam_slurm_2.sh"
path_script_convert_cram_to_bam_3="${path_directory_process}/partner/scripts/samtools/convert_cram_to_bam.sh"

# Files.
#path_file_reference_genome=$(<"./paths/community/reference_alignment_human_genome_grch38.txt")
path_file_reference_genome="${path_directory_reference}/human_genome/grch38_bsi_alignment/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
path_file_reference_genome_index="${path_directory_reference}/human_genome/grch38_bsi_alignment/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
path_file_parallel_instances="${path_directory_parallel}/instances_parallel.txt"

# Initialize directory.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_parallel
# Initialize file.
rm $path_file_parallel_instances # caution

###############################################################################
# Organize parameters.

threads=2
report="true"

#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

###############################################################################
# Execute procedure.

##########
# Simple execution.
if false; then
  /usr/bin/bash $path_script_convert_cram_to_bam_3 \
  $path_file_source \
  $path_file_product \
  $path_file_reference_genome \
  $path_file_reference_genome_index \
  $threads \
  $report \
  $path_execution_samtools
fi

##########
# Parallel batch of job instances
# Test with three instances.
if false; then
  # Define explicit instances.
  # Organize information within multi-dimensional array.
  instances_parallel=()
  path_file_source="${path_directory_source}/AAK959-AAT-B.FC22K7H7LT3_L1_IAACACTGTTA-GCAAGTCTCA.cram"
  path_file_product="${path_directory_product}/AAK959-AAT-B.FC22K7H7LT3_L1_IAACACTGTTA-GCAAGTCTCA.bam"
  instances_parallel+=("${path_file_source};${path_file_product}")
  path_file_source="${path_directory_source}/AAK959-AAT-B.FC22K7H7LT3_L2_IAACACTGTTA-GCAAGTCTCA.cram"
  path_file_product="${path_directory_product}/AAK959-AAT-B.FC22K7H7LT3_L2_IAACACTGTTA-GCAAGTCTCA.bam"
  instances_parallel+=("${path_file_source};${path_file_product}")
  path_file_source="${path_directory_source}/AAK959-AAT-B.FC22K7H7LT3_L3_IAACACTGTTA-GCAAGTCTCA.cram"
  path_file_product="${path_directory_product}/AAK959-AAT-B.FC22K7H7LT3_L3_IAACACTGTTA-GCAAGTCTCA.bam"
  instances_parallel+=("${path_file_source};${path_file_product}")
  # Write to file parameters for job instances.
  for instance in "${instances_parallel[@]}"; do
    echo $instance >> $path_file_parallel_instances
  done
  # Call script to submit parallel batch of job instances.
  /usr/bin/bash $path_script_convert_cram_to_bam_1 \
  $path_file_parallel_instances \
  $path_directory_parallel \
  $path_file_reference_genome \
  $path_file_reference_genome_index \
  $threads \
  $report \
  $path_script_convert_cram_to_bam_2 \
  $path_script_convert_cram_to_bam_3 \
  $path_execution_samtools
fi

##########
# Parallel batch of job instances
if true; then
  # Define explicit instances.
  # Collect files from source directory.
  #cd $path_directory_source
  # Bash version 4.4 introduced the "-d" option for "readarray".
  #readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
  paths_file_source=()
  while IFS= read -r -d $'\0'; do
    paths_file_source+=("$REPLY")
  done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.cram" -print0)
  count_paths_file_source=${#paths_file_source[@]}
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "script:"
    echo $0 # Print full file path to script.
    echo "1_convert_cram_to_bam_adipose.sh"
    echo "----------"
    echo "count of source files: " $count_paths_file_source
    echo "----------"
  fi
  # Organize information within multi-dimensional array.
  instances_parallel=()
  for path_file_source in "${paths_file_source[@]}"; do
    # Extract base name of file.
    name_base_file="$(basename $path_file_source .cram)"
    path_file_product="${path_directory_product}/${name_base_file}.bam" # hopefully unique
    instance="${path_file_source};${path_file_product}"
    instances_parallel+=$instance
  done
  # Write to file parameters for job instances.
  for instance in "${instances_parallel[@]}"; do
    echo $instance >> $path_file_parallel_instances
  done
  # Call script to submit parallel batch of job instances.
  /usr/bin/bash $path_script_convert_cram_to_bam_1 \
  $path_file_parallel_instances \
  $path_directory_parallel \
  $path_file_reference_genome \
  $path_file_reference_genome_index \
  $threads \
  $report \
  $path_script_convert_cram_to_bam_2 \
  $path_script_convert_cram_to_bam_3 \
  $path_execution_samtools
fi



##########
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Convert genomic or transcriptomic sequence data from CRAM to BAM file format."
  echo "----------"
fi


###############################################################################
# End.
