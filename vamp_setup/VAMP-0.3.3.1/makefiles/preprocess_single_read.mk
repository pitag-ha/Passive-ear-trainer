##########################
# Preprocess FASTQ File
# Steps for single read
#   Quality and Adatper trimming
#   Artifact filtering 
#
# author Lance Parsons <lparsons@princeton.edu>
##########################

# Setup read targets
READ$(READ): READ := $(READ)
READ$(READ)%: READ := $(READ)
CUTADAPT_FASTQ_FILE := ${CUTADAPT_FILEBASE}_read${READ}.fastq.gz
CUTADAPT_LOG_FILE := ${CUTADAPT_FILEBASE}_read${READ}.log
CUTADAPT_FASTQC_FILE := ${CUTADAPT_FILEBASE}_read${READ}_fastqc.zip
ARTIFACT_FASTQ_FILE := ${ARTIFACT_FILTERED_FILEBASE}_read${READ}.fastq
ARTIFACT_LOG_FILE := ${ARTIFACT_FILTERED_FILEBASE}_read${READ}.log
ARTIFACT_FASTQC_FILE := ${ARTIFACT_FILTERED_FILEBASE}_read${READ}_fastqc.zip
READ$(READ): ${FASTQC_FILE} ${CUTADAPT_FASTQC_FILE} ${ARTIFACT_FASTQC_FILE}


# Initial FastQC
################
#FASTQ_BASENAME:=$(basename $(notdir ${INPUT_FASTQ}))
#FASTQC_FILE:=${OUTPUT_FILEBASE}_read${READ}_fastqc.zip
#FASTQC_DIR:=${OUTPUT_FILEBASE}_read${READ}_fastqc
#${FASTQC_FILE} : ${INPUT_FASTQ}
#	${FASTQC} "$<" -o ${OUTPUT_DIR}
#	mv "${OUTPUT_DIR}/${FASTQ_BASENAME}_fastqc" "${OUTPUT_FILEBASE}_read${READ}_fastqc"
#	mv "${OUTPUT_DIR}/${FASTQ_BASENAME}_fastqc.zip" "$@"
#
#TEST: ${FASTQC_FILE}

# Cutadapt Quality and Adapter Trimming
#######################################
# Cuadapt default parameters
CUTADAPT_ERROR_RATE?=0.1 
CUTADAPT_MIN_OVERLAP_LENGTH?=5
CUTADAPT_MIN_LENGTH?=15
CUTADAPT_QUALITY_CUTOFF?=30

# Setup Cutadapt Adapters
CUTADAPT_ADAPTERS_3PRIME?=
CUTADAPT_ADAPTERS_ANYWHERE?=
ADAPTER_STRING:=$(foreach adapter, ${CUTADAPT_ADAPTERS_3PRIME},-a ${adapter} )
ifneq (${CUTADAPT_ADAPTERS_ANYWHERE},)
	ADAPTER_STRING+=$(foreach adapter, ${CUTADAPT_ADAPTERS_ANYWHERE},-b ${adapter} )
endif

.INTERMEDIATE: ${CUTADAPT_FASTQ_FILE}
${CUTADAPT_FASTQ_FILE} : ${INPUT_FASTQ} | ${OUTPUT_DIR}
	${CUTADAPT} -f fastq ${ADAPTER_STRING} -e ${CUTADAPT_ERROR_RATE} -O ${CUTADAPT_MIN_OVERLAP_LENGTH} -m ${CUTADAPT_MIN_LENGTH} -q ${CUTADAPT_QUALITY_CUTOFF} -o "$@" "$<" 1>${CUTADAPT_FILEBASE}_read${READ}.log

# FASTQC report of filtered file
.SECONDARY : ${CUTADAPT_FASTQC_FILE}
${CUTADAPT_FASTQC_FILE} : ${CUTADAPT_FASTQ_FILE}
	${FASTQC} "$<"

	
# Remove low complexity sequences
.INTERMEDIATE: ${ARTIFACT_FASTQ_FILE}
${ARTIFACT_FASTQ_FILE} : ${CUTADAPT_FASTQ_FILE}
	gunzip -c "$<" | ${FASTX_ARTIFACT_FILTER} -Q33 -v -o "$@" - 1>${ARTIFACT_FILTERED_FILEBASE}_read${READ}.log

# FASTQC report of clipped file
.SECONDARY : ${ARTIFACT_FASTQC_FILE}
${ARTIFACT_FASTQC_FILE} : ${ARTIFACT_FASTQ_FILE}
	${FASTQC} "$<"
