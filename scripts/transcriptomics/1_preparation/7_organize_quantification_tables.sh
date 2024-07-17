#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 July 2024
# Date, last execution or modification: 17 July 2024
# Review: TCW; 17 July 2024
###############################################################################
# Note

# HTSeq includes auxiliary information about genes in the first few columns of
# the quantification report; however, these columns do not have header names.
# This script manages the insertion of names for these auxiliary columns.

# This script also manages the translation of identifiers for samples for the
# sake of increasing deidentification anonymity.



###############################################################################
# Organize arguments.



###############################################################################
# Organize paths.



# Directories.
cd ~
path_directory_parent_project=$(<"./paths/endocrinology/parent_tcameronwaller.txt")
#path_directory_reference="${path_directory_parent_project}/reference"
#path_directory_tool="${path_directory_parent_project}/tool"
path_directory_process="${path_directory_parent_project}/process"
path_directory_dock="${path_directory_process}/dock"

path_directory_source="${path_directory_dock}/quantification_2024-07-14"
path_directory_product="${path_directory_source}/organization"
stamp_date=$(date +%Y-%m-%d)
path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.
path_file_source_adipose="${path_directory_source}/quantification_rna_reads_gene_adipose.tsv"
path_file_source_muscle="${path_directory_source}/quantification_rna_reads_gene_muscle.tsv"
path_file_temporary_adipose="${path_directory_temporary}/temporary_${stamp_date}_adipose.tsv"
path_file_temporary_muscle="${path_directory_temporary}/temporary_${stamp_date}_muscle.tsv"
path_file_translations_adipose="${path_directory_product}/table_sample_identifiers_adipose.tsv"
path_file_translations_muscle="${path_directory_product}/table_sample_identifiers_muscle.tsv"
path_file_product_adipose="${path_directory_product}/quantification_rna_reads_gene_adipose.tsv"
path_file_product_muscle="${path_directory_product}/quantification_rna_reads_gene_muscle.tsv"

# Scripts.
path_script_name_auxiliary="${path_directory_process}/partner/scripts/htseq/name_table_first_auxiliary_columns.sh"
path_script_translate_samples="${path_directory_process}/partner/scripts/htseq/translate_table_column_sample_identifiers.sh"

# Executable handles.

# Initialize directory.
#rm -r $path_directory_product # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary

# Initialize file.
rm $path_file_temporary_adipose
rm $path_file_temporary_muscle
rm $path_file_product_adipose
rm $path_file_product_muscle


###############################################################################
# Organize parameters.

threads=1
report="true"

#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

###############################################################################
# Execute procedure.

##########
# Organize information in tables.

# Insert header names for columns of auxiliary information about genes.
/usr/bin/bash $path_script_name_auxiliary \
$path_file_source_adipose \
$path_file_temporary_adipose \
$report
/usr/bin/bash $path_script_name_auxiliary \
$path_file_source_muscle \
$path_file_temporary_muscle \
$report

# Translate identifiers for samples.
/usr/bin/bash $path_script_translate_samples \
$path_file_translations_adipose \
$path_file_temporary_adipose \
$path_file_product_adipose \
$report
/usr/bin/bash $path_script_translate_samples \
$path_file_translations_muscle \
$path_file_temporary_muscle \
$path_file_product_muscle \
$report

##########
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "7_organize_quantification_tables.sh"
  echo "----------"
fi


##########
# Remove directory of temporary, intermediate files.
rm -r $path_directory_temporary


###############################################################################
# End.
