#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 13 July 2024
# Date, last execution or modification: 14 July 2024
# Review: TCW; 14 July 2024
###############################################################################
# Note

# View header of file in BAM or CRAM format to discern details of its prior
# processing.
# samtools head --headers 100 /path/to/input/file.bam
# samtools view --header-only --output /path/to/output/header/file.txt /path/to/input/file.bam

# Slurm batch job: 10545340
# - instances: 2
# - date: 14 July 2024
# threads=16
###SBATCH --nodes=1                            # count of cluster nodes (CPUs)
###SBATCH --ntasks-per-node=16                 # count of CPU cores or threads on node
###SBATCH --mem=10G                            # memory per node (per CPU)
###SBATCH --time=0-48:00:00                    # time allocation request (days-hours:minutes:seconds)


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

path_directory_source="${path_directory_dock}/consolidation_adipose_2024-05-31/filter_sort_index" # 14 July 2024
path_directory_product="${path_directory_dock}/test_test_test"

# Files.

# Scripts.
path_script_test="${path_directory_process}/partner/scripts/htseq/test_variable_expansion.sh"

# Executable handles.
path_environment_htseq="${path_directory_tool}/python/environments/htseq"

# Initialize directory.
#rm -r $path_directory_product # caution
#mkdir -p $path_directory_product

###############################################################################
# Organize parameters.

threads=16
report="true"

#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

###############################################################################
# Execute procedure.


paths_file_source=()
while IFS= read -r -d $'\0'; do
  paths_file_source+=("$REPLY")
done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
IFS=";"; paths_file_source_expansion="${paths_file_source[*]}"
echo $paths_file_source_expansion

# Divide the array of files into chunks of 30 or fewer elements.
count_chunk=10
for ((index=0; index < ${#paths_file_source[@]}; index+=count_chunk)); do
  chunk=("${paths_file_source[@]:index:count_chunk}")
  IFS=";"; chunk_expansion="${chunk[*]}"
  echo $chunk_expansion
done

###############################################################################
# End.
