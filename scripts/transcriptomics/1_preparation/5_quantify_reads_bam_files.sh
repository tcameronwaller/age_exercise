#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 13 July 2024
# Date, last execution or modification: 13 July 2024
# Review: TCW; 13 July 2024
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

path_directory_source_adipose="${path_directory_dock}/consolidation_adipose_2024-05-31/filter_sort_index"
path_directory_source_muscle="${path_directory_dock}/consolidation_muscle_2022-07-13/filter_sort_index"
path_directory_product_adipose="${path_directory_dock}/consolidation_adipose_2024-05-31/quantification"
path_directory_product_muscle="${path_directory_dock}/consolidation_muscle_2022-07-13/quantification"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product_adipose}/temporary_${stamp_date}" # hopefully unique
path_directory_parallel="${path_directory_product_adipose}/parallel"

# Files.
path_file_annotation_gtf_gzip="${path_directory_reference}/human_genome/gencode/grch38/annotation/gencode.v46.primary_assembly.annotation.gtf.gz"
#path_file_annotation_gtf="${path_directory_temporary}/gencode.v46.primary_assembly.annotation.gtf"
path_file_product_adipose="${path_directory_product_adipose}/test_quantification.tsv"
path_file_product_muscle="${path_directory_product_muscle}/test_quantification.tsv"
path_file_parallel_instances="${path_directory_parallel_adipose}/instances_parallel.txt"

# Scripts.
path_script_quantify_1="${path_directory_process}/partner/scripts/htseq/quantify_rna_reads_slurm_1.sh"
path_script_quantify_2="${path_directory_process}/partner/scripts/htseq/quantify_rna_reads_slurm_2.sh"
path_script_quantify_3="${path_directory_process}/partner/scripts/htseq/quantify_rna_reads.sh"

# Executable handles.
path_environment_htseq="${path_directory_tool}/python/environments/htseq"

# Initialize directory.
rm -r $path_directory_product_adipose # caution
mkdir -p $path_directory_product_adipose
rm -r $path_directory_product_muscle # caution
mkdir -p $path_directory_product_muscle
#mkdir -p $path_directory_temporary
mkdir -p $path_directory_parallel
# Initialize file.
rm $path_file_parallel_instances # caution

###############################################################################
# Organize parameters.

threads=10
report="true"

#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

###############################################################################
# Execute procedure.


##########
if true; then
  # Define explicit instances for parallel batch of jobs.
  # Organize information within multi-dimensional array.
  instances_parallel=()
  instance="${path_directory_source_adipose};${path_file_product_adipose}"
  instances_parallel+=($instance)
  instance="${path_directory_source_muscle};${path_file_product_muscle}"
  instances_parallel+=($instance)
  # Write to file parameters for job instances.
  for instance in "${instances_parallel[@]}"; do
    echo $instance >> $path_file_parallel_instances
  done
  # Call script to submit parallel batch of job instances.
  /usr/bin/bash $path_script_quantify_1 \
  $path_file_parallel_instances \
  $path_directory_parallel \
  $path_file_annotation_gtf_gzip \
  $threads \
  $report \
  $path_script_quantify_2 \
  $path_script_quantify_3 \
  $path_environment_htseq
fi


##########
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "5_quantify_reads_bam_files.sh"
  echo "Quantify reads allocatable to specific genomic features."
  echo "----------"
  echo "Adipose"
  #echo "count of source files: " $count_paths_file_source_adipose
  echo "----------"
  echo "Muscle"
  #echo "count of source files: " $count_paths_file_source_muscle
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
#rm -r $path_directory_temporary

###############################################################################
# End.
