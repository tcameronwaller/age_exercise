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

# In next-generation sequencing, each sample was analyzed simultaneously on
# 8 separate flow cells, such that there are 8 different files in BAM format
# for each sample. It is necessary to merge these 8 files for each sample.

# Slurm batch job: ___
# - instances: ___
# - date: ___

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
path_directory_source="${path_directory_dock}/20240531_LH00386_0066_A22K7H7LT3/bam"
#path_directory_source=$("./paths/endocrinology/transcriptomics_muscle.txt")
#path_directory_source="${path_directory_dock}/test_lanza_rnaseq_adipose_2024/raw"
#path_directory_source="${path_directory_dock}/test_lanza_rnaseq_muscle_2022/raw"
path_directory_product="${path_directory_dock}/place_holder_temporary"
#path_directory_product="${path_directory_dock}/test_lanza_rnaseq_adipose_2024/bam"
stamp_date=$(date +%Y-%m-%d)
path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique
path_directory_parallel="${path_directory_product}/parallel"
path_directory_parallel_sample_files="${path_directory_product}/parallel/sample_files"

# Executable handles.
path_execution_samtools="${path_directory_tool}/samtools-1.20/bin/samtools"

# Scripts.
path_script_convert_cram_to_bam_1="${path_directory_process}/partner/scripts/samtools/convert_cram_to_bam_slurm_1.sh"
path_script_convert_cram_to_bam_2="${path_directory_process}/partner/scripts/samtools/convert_cram_to_bam_slurm_2.sh"
path_script_convert_cram_to_bam_3="${path_directory_process}/partner/scripts/samtools/convert_cram_to_bam.sh"

# Files.
#path_file_reference_genome=$(<"./paths/community/reference_alignment_human_genome_grch38.txt")
#path_file_reference_genome="${path_directory_reference}/human_genome/grch38_bsi_alignment/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#path_file_reference_genome_index="${path_directory_reference}/human_genome/grch38_bsi_alignment/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
path_file_temporary_1="${path_directory_temporary}/temporary_${stamp_date}_1.txt"
path_file_parallel_instances="${path_directory_parallel}/instances_parallel.txt"

# Initialize directory.
#rm -r $path_directory_product # caution
rm -r $path_directory_temporary # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
mkdir -p $path_directory_parallel
mkdir -p $path_directory_parallel_sample_files
# Initialize file.
rm $path_file_temporary_1 # caution
#rm $path_file_parallel_instances # caution

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
# 1. Determine unique identifiers of samples

# Collect paths to files from source directory.
#cd $path_directory_source
# Bash version 4.4 introduced the "-d" option for "readarray".
#readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
paths_file_source=()
while IFS= read -r -d $'\0'; do
  paths_file_source+=("$REPLY")
done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
count_paths_file_source=${#paths_file_source[@]}

# Extract base names of files from source directory.
names_file_base=()
for path_file_source in "${paths_file_source[@]}"; do
  # Extract base name of file.
  name_base_file="$(basename $path_file_source .bam)"
  names_file_base+=($name_base_file)
done
count_names_file_base=${#names_file_base[@]}

# Write array to temporary file.
for item in "${names_file_base[@]}"; do
  # Replace multi-character delimiter with single-character delimiter.
  echo $item | sed 's/_L/;/g' >> $path_file_temporary_1
done

# Extract identifiers of samples from base names of files.
identifiers_sample=()
#input=$names_file_base
input=$path_file_temporary_1
while IFS=$'\n' read -r -a array_lines; do
  for line in "${array_lines}"; do
    # Separate segments within current line.
    IFS=$';' read -r -a array_segments <<< "${line}"
    # Select identifier of sample from segments of current file's base name.
    identifier_sample="${array_segments[0]}"
    identifiers_sample+=($identifier_sample)
  done
done < "${input}"
#done <<< "${input}"
count_identifiers_sample=${#identifiers_sample[@]}

# Select unique identifiers of samples.
declare -A temporary_unique # initialize an associative array
for identifier_sample in "${identifiers_sample[@]}"; do
  temporary_unique[$identifier_sample]=0 # assign a key-value pair
done
# Organize array from keys of associative array.
#identifiers_sample_unique=${!temporary_unique[@]} # transfers a space-delimited list of keys
IFS=$' ' read -r -a identifiers_sample_unique <<< "${!temporary_unique[@]}"
count_identifiers_sample_unique=${#identifiers_sample_unique[@]}

##########
# 2. Define explicit instances for parallel batch of jobs.
# Organize information within multi-dimensional array.
instances_parallel=()
# Determine paths to files for each unique sample identifier.
for identifier_sample in "${identifiers_sample_unique[@]}"; do
  # Collect paths to files from source directory corresponding to current sample identifier.
  paths_file_sample=()
  while IFS= read -r -d $'\0'; do
    paths_file_sample+=("$REPLY")
  done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "${identifier_sample}*.bam" -print0)
  # Write paths to files to a parameter file.
  path_file_sample_files="${path_directory_parallel_sample_files}"
  for item in "${names_file_base[@]}"; do
    # Replace multi-character delimiter with single-character delimiter.
    echo $item | sed 's/_L/;/g' >> $path_file_temporary_1
  done



  # Extract base name of file.
  name_base_file="$(basename $path_file_source .cram)"
  path_file_product="${path_directory_product}/${name_base_file}.bam" # hopefully unique
  # Define instance.
  instance="${identifier_sample};${path_file_source};${path_file_product}"
  instances_parallel+=($instance)
done


##########
# 3. Parallel batch of job instances

if true; then
  # Define explicit instances.
  # Collect paths to files from source directory.
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
    instances_parallel+=($instance)
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
  echo "2_merge_files_adipose.sh"
  echo "Merge files for identical sample analyzed in parallel flow cells."
  echo "----------"
  echo "count of source files: " $count_paths_file_source
  echo "count of base names of source files: " $count_names_file_base
  echo "example of base name: " "${names_file_base[0]}"
  echo "count of sample identifiers: " $count_identifiers_sample
  echo "example of sample identifier: " "${identifiers_sample[0]}"
  echo "count of unique sample identifiers: " $count_identifiers_sample_unique
  echo "example of unique sample identifier: " "${identifiers_sample_unique[0]}"
  echo "example of unique sample identifier: " "${identifiers_sample_unique[1]}"
  echo "example of unique sample identifier: " "${identifiers_sample_unique[2]}"
  echo "example of unique sample identifier: " "${identifiers_sample_unique[3]}"
  echo "example of unique sample identifier: " "${identifiers_sample_unique[-1]}"
  echo "----------"
fi


##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

###############################################################################
# End.
