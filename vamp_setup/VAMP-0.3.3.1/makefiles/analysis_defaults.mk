###############################
# Test VAMP Makefile 
#  
# author Lance Parsons  <lparsons@princeton.edu>
################################


# Analysis specific configuration
###################################

# Makefile directory
VAMP_HOME?=${vamp_home_dir}
MAKEFILE_DIR:=${VAMP_HOME}/makefiles

# Fastq input file
READ1_FASTQ:=${forward_read}
READ2_FASTQ:=${reverse_read}

# Output file location and name
OUTPUT_DIR:=${output_dir}
OUTPUT_BASE:=${output_base}
OUTPUT_FILEBASE:=${OUTPUT_DIR}/${OUTPUT_BASE}



# Preprocessing Configuration
######################################################################
CUTADAPT_ADAPTERS_3PRIME:=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
CUTADAPT_ADAPTERS_ANYWHERE:=
CUTADAPT_ERROR_RATE:=0.05 
CUTADAPT_MIN_OVERLAP_LENGTH:=10
CUTADAPT_MIN_LENGTH:=15
CUTADAPT_QUALITY_CUTOFF:=30

TQS_THRESH:=20
TQS_CONSEC:=20

${set_host_reference_fasta}
${set_host_reference_bowtie_index}

# Paired Sequence Match Configuration
#####################################
# Old Illumina Read IDs use: "/"
# New Illumina Read IDs use: " "
# Others might use: ":"
ID_SPLIT_CHAR:=" "

BOWTIE_THREADS:=1
BOWTIE_OPTIONS:=--seedmms 3 --maqerr 120 --maxins 700 -p ${BOWTIE_THREADS}

# Analysis modules
#####################################################################

# Import global configuration
include ${MAKEFILE_DIR}/config.mk

# Import Proprocessing targets
include ${MAKEFILE_DIR}/preprocess.mk


# Execution targets
.DEFAULT_GOAL := all

all: preprocess

.PHONY: clean-secondary
clean-secondary: clean-secondary-preprocess

.PHONY: clean-all
clean-all: clean-all-preprocess
