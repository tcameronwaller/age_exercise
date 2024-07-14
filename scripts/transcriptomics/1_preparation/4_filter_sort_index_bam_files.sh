#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 July 2024
# Date, last execution or modification: 12 July 2024
# Review: TCW; 12 July 2024
###############################################################################
# Note

# View header of file in BAM or CRAM format to discern details of its prior
# processing.
# samtools head --headers 100 /path/to/input/file.bam
# samtools view --header-only --output /path/to/output/header/file.txt /path/to/input/file.bam

# In next-generation sequencing, each sample was analyzed simultaneously on
# 8 separate flow cells, such that there are 8 different files in BAM format
# for each sample. It is necessary to merge these 8 files for each sample.

# Adipose
# Slurm batch job: 10501635
# - instances: 154
# - date: 12 July 2024

# Muscle
# Slurm batch job: 10522113
# - instances: 185
# - date: 13 July 2024


###############################################################################


###############################################################################
# Organize arguments.

###############################################################################
# Organize paths.

# Directories.
cd ~
path_directory_parent_project=$(<"./paths/endocrinology/parent_tcameronwaller.txt")
#path_directory_reference="${path_directory_parent_project}/reference"
path_directory_tool="${path_directory_parent_project}/tool"
path_directory_process="${path_directory_parent_project}/process"
path_directory_dock="${path_directory_process}/dock"

path_directory_source_adipose="${path_directory_dock}/consolidation_adipose_2024-05-31"
path_directory_source_muscle="${path_directory_dock}/consolidation_muscle_2022-07-13"
path_directory_product_adipose="${path_directory_dock}/consolidation_adipose_2024-05-31/filter_sort_index"
path_directory_product_muscle="${path_directory_dock}/consolidation_muscle_2022-07-13/filter_sort_index"
path_directory_parallel_adipose="${path_directory_product_adipose}/parallel"
path_directory_parallel_muscle="${path_directory_product_muscle}/parallel"

# Files.
path_file_parallel_instances_adipose="${path_directory_parallel_adipose}/instances_parallel.txt"
path_file_parallel_instances_muscle="${path_directory_parallel_muscle}/instances_parallel.txt"

# Scripts.
path_script_sort_index_bam_1="${path_directory_process}/partner/scripts/samtools/filter_sort_index_bam_slurm_1.sh"
path_script_sort_index_bam_2="${path_directory_process}/partner/scripts/samtools/filter_sort_index_bam_slurm_2.sh"
path_script_sort_index_bam_3="${path_directory_process}/partner/scripts/samtools/filter_sort_index_bam.sh"

# Executable handles.
path_execution_samtools="${path_directory_tool}/samtools-1.20/bin/samtools"

# Initialize directory.
#rm -r $path_directory_product_adipose # caution
#mkdir -p $path_directory_product_adipose
#mkdir -p $path_directory_parallel_adipose
rm -r $path_directory_product_muscle # caution
mkdir -p $path_directory_product_muscle
mkdir -p $path_directory_parallel_muscle
# Initialize file.
#rm $path_file_parallel_instances_adipose # caution
rm $path_file_parallel_instances_muscle # caution

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
# Adipose
if false; then
  # Define explicit instances for parallel batch of jobs.
  # Collect paths to files from source directory.
  #cd $path_directory_source
  # Bash version 4.4 introduced the "-d" option for "readarray".
  #readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
  paths_file_source=()
  while IFS= read -r -d $'\0'; do
    paths_file_source+=("$REPLY")
  done < <(find $path_directory_source_adipose -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
  count_paths_file_source_adipose=${#paths_file_source[@]}
  # Organize information within multi-dimensional array.
  instances_parallel=()
  for path_file_source in "${paths_file_source[@]}"; do
    # Extract base name of file.
    name_base_file="$(basename $path_file_source .bam)"
    path_file_product="${path_directory_product_adipose}/${name_base_file}.bam" # hopefully unique
    instance="${path_file_source};${path_file_product}"
    instances_parallel+=($instance)
  done
  # Write to file parameters for job instances.
  for instance in "${instances_parallel[@]}"; do
    echo $instance >> $path_file_parallel_instances_adipose
  done
  # Call script to submit parallel batch of job instances.
  /usr/bin/bash $path_script_sort_index_bam_1 \
  $path_file_parallel_instances_adipose \
  $path_directory_parallel_adipose \
  $threads \
  $report \
  $path_script_sort_index_bam_2 \
  $path_script_sort_index_bam_3 \
  $path_execution_samtools
fi

##########
# Muscle
if true; then
  # Define explicit instances for parallel batch of jobs.
  # Collect paths to files from source directory.
  #cd $path_directory_source
  # Bash version 4.4 introduced the "-d" option for "readarray".
  #readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
  paths_file_source=()
  while IFS= read -r -d $'\0'; do
    paths_file_source+=("$REPLY")
  done < <(find $path_directory_source_muscle -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
  count_paths_file_source_muscle=${#paths_file_source[@]}
  # Organize information within multi-dimensional array.
  instances_parallel=()
  for path_file_source in "${paths_file_source[@]}"; do
    # Extract base name of file.
    name_base_file="$(basename $path_file_source .bam)"
    path_file_product="${path_directory_product_muscle}/${name_base_file}.bam" # hopefully unique
    instance="${path_file_source};${path_file_product}"
    instances_parallel+=($instance)
  done
  # Write to file parameters for job instances.
  for instance in "${instances_parallel[@]}"; do
    echo $instance >> $path_file_parallel_instances_muscle
  done
  # Call script to submit parallel batch of job instances.
  /usr/bin/bash $path_script_sort_index_bam_1 \
  $path_file_parallel_instances_muscle \
  $path_directory_parallel_muscle \
  $threads \
  $report \
  $path_script_sort_index_bam_2 \
  $path_script_sort_index_bam_3 \
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
  echo "4_filter_sort_index_bam_files.sh"
  echo "Filter reads, sort reads by name, and create index for file in BAM format."
  echo "----------"
  echo "Adipose"
  #echo "count of source files: " $count_paths_file_source_adipose
  echo "----------"
  echo "Muscle"
  echo "count of source files: " $count_paths_file_source_muscle
  echo "----------"
fi

###############################################################################
# End.
