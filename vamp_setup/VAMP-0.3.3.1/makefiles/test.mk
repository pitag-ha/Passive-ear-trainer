###############################
# Test VAMP Makefile 
#  
# author Lance Parsons  <lparsons@princeton.edu>
################################


# Analysis specific configuration
###################################

# Makefile directory
VAMP_HOME:=/Users/lparsons/Documents/workspace/vamp

# Fastq input file
READ1_FASTQ:=/Volumes/Fantom\ HD/test_data/vamp_test_data/drosophila-santomea-sto4-read1.fastq
READ2_FASTQ:=/Volumes/Fantom\ HD/test_data/vamp_test_data/drosophila-santomea-sto4-read2.fastq

# Output file location and name
OUTPUT_DIR:=/Volumes/Fantom\ HD/test_data/vamp_test_output
OUTPUT_BASE:=vamp_test


######################################################################

# Preprocessing Configuration
CUTADAPT_ADAPTERS_3PRIME:=AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
CUTADAPT_ADAPTERS_ANYWHERE:=
CUTADAPT_ERROR_RATE:=0.1 
CUTADAPT_MIN_OVERLAP_LENGTH:=5
CUTADAPT_MIN_LENGTH:=15
CUTADAPT_QUALITY_CUTOFF:=30

HOST_REFERENCE_FASTA:=/Users/lparsons/Documents/projects/sequencing/genomes/phiX/phi_plus_SNPs.fa

BOWTIE_THREADS:=4

#####################################################################
MAKEFILE_DIR:=${VAMP_HOME}/makefiles
OUTPUT_FILEBASE:=${OUTPUT_DIR}/${OUTPUT_BASE}

# Import global configuration
include ${MAKEFILE_DIR}/config.mk

# Import Proprocessing targets
include ${MAKEFILE_DIR}/preprocess.mk


# Execution targets
all: preprocess

.PHONY: clean-secondary
clean-secondary: clean-secondary-preprocess

.PHONY: clean-all
clean-all:
	rm -rf $(OUTPUT_FILEBASE)*
