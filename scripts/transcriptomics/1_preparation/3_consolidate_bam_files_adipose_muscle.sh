#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 July 2024
# Date, last execution or modification: 11 July 2024
# Review: TCW; 11 July 2024
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
#path_directory_reference="${path_directory_parent_project}/reference"
path_directory_tool="${path_directory_parent_project}/tool"
path_directory_process="${path_directory_parent_project}/process"
path_directory_dock="${path_directory_process}/dock"

path_directory_source_adipose="${path_directory_dock}/20240531_LH00386_0066_A22K7H7LT3/merge_test"
path_directory_source_muscle="${path_directory_dock}/220713_A00938_0293_AHVFVFDSX3/test"

path_directory_product_adipose="${path_directory_dock}/consolidation_adipose_2024-05-31"
path_directory_product_muscle="${path_directory_dock}/consolidation_muscle_2022-07-13"

stamp_date=$(date +%Y-%m-%d)
path_directory_temporary_adipose="${path_directory_product_adipose}/temporary_${stamp_date}_adipose" # hopefully unique
path_directory_temporary_muscle="${path_directory_product_muscle}/temporary_${stamp_date}_muscle" # hopefully unique

# Executable handles.

# Scripts.

# Files.

# Initialize directory.
rm -r $path_directory_product_adipose # caution
rm -r $path_directory_product_muscle # caution
#rm -r $path_directory_temporary_adipose # caution
#rm -r $path_directory_temporary_muscle # caution
mkdir -p $path_directory_product_adipose
mkdir -p $path_directory_product_muscle
mkdir -p $path_directory_temporary_adipose
mkdir -p $path_directory_temporary_muscle
# Initialize file.

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
# 1. Copy files in BAM format to temporary directory.
cp $path_directory_source_adipose/*.bam $path_directory_temporary_adipose
cp $path_directory_source_muscle/*.bam $path_directory_temporary_muscle

##########
# 2. Modify the names of files.

# Adipose
# Collect paths to files from source directory.
#cd $path_directory_source
# Bash version 4.4 introduced the "-d" option for "readarray".
#readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
paths_file_source=()
while IFS= read -r -d $'\0'; do
  paths_file_source+=("$REPLY")
done < <(find $path_directory_temporary_adipose -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
count_paths_file_source_adipose=${#paths_file_source[@]}
# Modify the names of files.
for path_file_source in "${paths_file_source[@]}"; do
  # Extract base name of file.
  name_base_file="$(basename $path_file_source .bam)"
  # Determine new path and name for file.
  name_base_file_new=$name_base_file
  path_file_product="${path_directory_product_adipose}/${name_base_file_new}.bam"
  # Move file.
  mv $path_file_source $path_file_product
done

# Muscle
# Collect paths to files from source directory.
#cd $path_directory_source
# Bash version 4.4 introduced the "-d" option for "readarray".
#readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
paths_file_source=()
while IFS= read -r -d $'\0'; do
  paths_file_source+=("$REPLY")
done < <(find $path_directory_temporary_muscle -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
count_paths_file_source_muscle=${#paths_file_source[@]}
# Modify the names of files.
for path_file_source in "${paths_file_source[@]}"; do
  # Extract base name of file.
  name_base_file="$(basename $path_file_source .bam)"
  # Replace multi-character delimiter with single-character delimiter.
  echo $name_base_file | sed 's/_L/;/g' > $name_base_file_simple
  # Separate segments within current line.
  IFS=$';' read -r -a array_segments <<< "${name_base_file_simple}"
  # Select identifier of sample from segments of current file's base name.
  name_base_file_new="${array_segments[0]}"
  # Determine new path and name for file.
  path_file_product="${path_directory_product_muscle}/${name_base_file_new}.bam"
  # Move file.
  mv $path_file_source $path_file_product
done

##########
# Report.

if false; then
  if [ "$report" == "true" ]; then
    echo "----------"
    echo "----------"
    echo "----------"
    echo "Script:"
    echo $0 # Print full file path to script.
    echo "3_consolidate_bam_files_adipose_muscle.sh"
    echo "Consolidate files for adipose and muscle tissues."
    echo "----------"
    echo "Adipose"
    echo "count of source files: " $count_paths_file_source_adipose
    echo "----------"
    echo "Muscle"
    echo "count of source files: " $count_paths_file_source_muscle
    echo "----------"
  fi
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary_adipose
rm -r $path_directory_temporary_muscle

###############################################################################
# End.
