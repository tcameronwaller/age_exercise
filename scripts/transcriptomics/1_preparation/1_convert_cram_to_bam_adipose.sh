#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 5 July 2024
# Date, last execution or modification: 5 July 2024
# Review: TCW; 5 July 2024
###############################################################################
# Note


################################################################################


################################################################################
# Organize arguments.

################################################################################
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
path_directory_source="${path_directory_dock}/dock/test_lanza_rnaseq_adipose_2024/raw"
#path_directory_source="${path_directory_dock}/dock/test_lanza_rnaseq_muscle_2022/raw"
path_directory_product="${path_directory_dock}/dock/test_lanza_rnaseq_adipose_2024/raw/bam"

# Executable handles.
path_execution_samtools="${path_directory_tool}/samtools-1.20/bin/samtools"

# Scripts.
path_script_convert_cram_to_bam="${path_directory_process}/partner/scripts/convert_cram_to_bam.sh"

# Files.
path_file_reference_genome=$(<"./paths/community/reference_alignment_human_genome_grch38.txt")
path_file_source="${path_directory_source}/AAK959-AAT-B.FC22K7H7LT3_L1_IAACACTGTTA-GCAAGTCTCA.cram"
path_file_product="${path_directory_product}/AAK959-AAT-B.FC22K7H7LT3_L1_IAACACTGTTA-GCAAGTCTCA.bam"

# Initialize directory.
#rm $path_directory_product
mkdir -p $path_directory_product

# Remove any previous version of the product file.
#rm $path_file_product

###############################################################################
# Organize parameters.

# Report.
report="true"
set -x # enable print commands to standard error
#set +x # disable print commands to standard error
set -v # enable print input to standard error
#set +v # disable print input to standard error

###############################################################################
# Execute procedure.

/usr/bin/bash $path_script_convert_cram_to_bam \
$path_file_source \
$name_file_product \
$name_file_reference_genome \
$path_execution_samtools \
$report

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
  echo "path to source file: " $path_file_source
  echo "path to product file: " $path_file_product
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
#rm -r $path_directory_product_temporary

###############################################################################
# End.
