##########################
# Preprocess FASTQ File 
#
# author Lance Parsons <lparsons@princeton.edu>
##########################

.DELETE_ON_ERROR:

# Filenames
CUTADAPT_FILEBASE:=${OUTPUT_FILEBASE}_cutadapt
ARTIFACT_FILTERED_FILEBASE:=${CUTADAPT_FILEBASE}_artifact_filtered
HOST_FILTERED_FILEBASE:=${ARTIFACT_FILTERED_FILEBASE}_host_filtered

HOST_REFERENCE_BOWTIE_INDEX?=$(basename ${HOST_REFERENCE_FASTA})
HOST_REFERENCE_BOWTIE_INDEX_FILES = ${HOST_REFERENCE_BOWTIE_INDEX}.1.ebwt ${HOST_REFERENCE_BOWTIE_INDEX}.2.ebwt ${HOST_REFERENCE_BOWTIE_INDEX}.3.ebwt ${HOST_REFERENCE_BOWTIE_INDEX}.4.ebwt ${HOST_REFERENCE_BOWTIE_INDEX}.rev.1.ebwt ${HOST_REFERENCE_BOWTIE_INDEX}.rev.2.ebwt

PREPROCESSING_SUMMARY:=${OUTPUT_FILEBASE}_preprocessing_summary.txt

# Bowtie Configuration
######################
BOWTIE_THREADS?=1
BOWTIE_OPTIONS?=-p ${BOWTIE_THREADS}

# TQSFastq Configuration
########################
#TQS_THRESH?=20
#TQS_CONSEC?=20

# Paired Sequence Match Configuration
#####################################
ID_SPLIT_CHAR?=" "


# Input Files
READ2_FASTQ?=
FASTQ_FILES:= ${READ1_FASTQ} ${READ2_FASTQ}
ifneq (${READ2_FASTQ},)
	READS:=1 2
else
	READS:=1
endif
CURRENT_READS:=${READS}

# Remove files
.PHONY: clean clean-secondary-preprocess clean-host-bowtie-index clean-all-preprocess

clean:
	rm -f ${CUTADAPT_FILEBASE}_read?.fastq.gz
	rm -f ${ARTIFACT_FILTERED_FILEBASE}_read?.fastq
	rm -f ${ARTIFACT_FILTERED_FILEBASE}_read?_paired.fastq
	rm -f ${OUTPUT_FILEBASE}_bowtie_mapped.sam
	rm -f ${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq

clean-secondary-preprocess: 
	rm -f ${CUTADAPT_FILEBASE}_read?.fastq.gz
	rm -f ${ARTIFACT_FILTERED_FILEBASE}_read?.fastq

clean-host-bowtie-index:
	rm -f ${HOST_REFERENCE_BOWTIE_INDEX_FILES}

clean-all-preprocess:
	rm -rf ${CUTADAPT_FILEBASE}*
	rm -rf ${ARTIFACT_FILTERED_FILEBASE}*
	rm -rf ${OUTPUT_FILEBASE}_bowtie_mapped.*
	rm -rf ${HOST_FILTERED_FILEBASE}*
	rm -rf ${PREPROCESSING_SUMMARY}

# Output Directory
${OUTPUT_DIR}:
	mkdir ${OUTPUT_DIR}


# READ 1
########

# Adapter and quality trim, filter artifacts
READ:=1
INPUT_FASTQ:=${READ1_FASTQ}
include ${MAKEFILE_DIR}/preprocess_single_read.mk

# Initial FastQC
READ1_FASTQ_BASENAME:=$(basename $(notdir ${READ1_FASTQ}))
ifeq "$(suffix ${READ1_FASTQ})" ".gz"
READ1_FASTQ_BASENAME:=$(basename $(basename $(notdir ${READ1_FASTQ})))
endif
READ1_FASTQC_FILE:=${OUTPUT_FILEBASE}_read${READ}_fastqc.zip
READ1_FASTQC_DIR:=${OUTPUT_FILEBASE}_read${READ}_fastqc

${READ1_FASTQC_FILE} : ${READ1_FASTQ}
	${FASTQC} "$<" -o ${OUTPUT_DIR}
	-mv "${OUTPUT_DIR}/${READ1_FASTQ_BASENAME}_fastqc" "${READ1_FASTQC_DIR}"
	-mv "${OUTPUT_DIR}/${READ1_FASTQ_BASENAME}_fastqc.zip" "$@"

${PREPROCESSING_SUMMARY}: READ1 ${READ1_FASTQC_FILE}

# READ 2
########
ifneq (${READ2_FASTQ},)

# Adapter and quality trim, filter artifacts
READ:=2
INPUT_FASTQ:=${READ2_FASTQ}
include ${MAKEFILE_DIR}/preprocess_single_read.mk

# Initial FastQC
READ2_FASTQ_BASENAME:=$(basename $(notdir ${READ2_FASTQ}))
ifeq "$(suffix ${READ2_FASTQ})" ".gz"
READ2_FASTQ_BASENAME:=$(basename $(basename $(notdir ${READ2_FASTQ})))
endif
READ2_FASTQC_FILE:=${OUTPUT_FILEBASE}_read${READ}_fastqc.zip
READ2_FASTQC_DIR:=${OUTPUT_FILEBASE}_read${READ}_fastqc
${READ2_FASTQC_FILE} : ${READ2_FASTQ} 
	${FASTQC} "$<" -o ${OUTPUT_DIR}
	-mv "${OUTPUT_DIR}/${READ2_FASTQ_BASENAME}_fastqc" "${READ2_FASTQC_DIR}"
	-mv "${OUTPUT_DIR}/${READ2_FASTQ_BASENAME}_fastqc.zip" "$@"

${PREPROCESSING_SUMMARY}: READ2 ${READ2_FASTQC_FILE}

# Combined filtered reads
.INTERMEDIATE: ${ARTIFACT_FILTERED_FILEBASE}_read1_paired.fastq ${ARTIFACT_FILTERED_FILEBASE}_read2_paired.fastq ${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq

${ARTIFACT_FILTERED_FILEBASE}_read1_paired.fastq: ${ARTIFACT_FILTERED_FILEBASE}_read1.fastq ${ARTIFACT_FILTERED_FILEBASE}_read2.fastq
	${PAIRED_SEQUENCE_MATCH} ${ARTIFACT_FILTERED_FILEBASE}_read1.fastq ${ARTIFACT_FILTERED_FILEBASE}_read2.fastq -p ${ARTIFACT_FILTERED_FILEBASE}_read1_paired.fastq -p ${ARTIFACT_FILTERED_FILEBASE}_read2_paired.fastq -s ${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq --id-split=${ID_SPLIT_CHAR} --index-in-memory

${ARTIFACT_FILTERED_FILEBASE}_read2_paired.fastq: ${ARTIFACT_FILTERED_FILEBASE}_read1_paired.fastq

${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq: ${ARTIFACT_FILTERED_FILEBASE}_read2_paired.fastq 

${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq.gz: ${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq
	gzip -f ${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq 

${PREPROCESSING_SUMMARY}: ${ARTIFACT_FILTERED_FILEBASE}_unpaired.fastq.gz

endif


# Remove Host Contaminants
##########################

# Build Bowtie Index
#${HOST_REFERENCE_BOWTIE_INDEX} : ${HOST_REFERENCE_BOWTIE_INDEX}.1.ebwt
${HOST_REFERENCE_BOWTIE_INDEX}.1.ebwt: ${HOST_REFERENCE_FASTA}
	${BOWTIE_BUILD} "$<" ${HOST_REFERENCE_BOWTIE_INDEX}

ifeq (${READ2_FASTQ},)
# Single end reads
# Map Reads to Host Reference
${HOST_FILTERED_FILEBASE}_1.fastq.gz: ${ARTIFACT_FILTERED_FILEBASE}_read1.fastq ${HOST_REFERENCE_BOWTIE_INDEX}.1.ebwt
	gunzip -c "$<" | ${BOWTIE} -S ${BOWTIE_OPTIONS} --un ${HOST_FILTERED_FILEBASE}_1.fastq ${HOST_REFERENCE_BOWTIE_INDEX} - 1> ${OUTPUT_FILEBASE}_bowtie_mapped.sam 2> ${OUTPUT_FILEBASE}_bowtie_mapped.log
	gzip ${HOST_FILTERED_FILEBASE}.fastq
# FASTQC report of host filtered file
${HOST_FILTERED_FILEBASE}_1_fastqc.zip: ${HOST_FILTERED_FILEBASE}_1.fastq.gz
	${FASTQC} "$<"
# Add host filtering step to prereqs
${PREPROCESSING_SUMMARY}: ${HOST_FILTERED_FILEBASE}_1_fastqc.zip
else
# Paired end reads
${HOST_FILTERED_FILEBASE}_1.fastq.gz: ${ARTIFACT_FILTERED_FILEBASE}_read1_paired.fastq ${HOST_REFERENCE_BOWTIE_INDEX}.1.ebwt
	${BOWTIE} -S ${BOWTIE_OPTIONS} --un ${HOST_FILTERED_FILEBASE}.fastq ${HOST_REFERENCE_BOWTIE_INDEX} -1 ${ARTIFACT_FILTERED_FILEBASE}_read1_paired.fastq -2 ${ARTIFACT_FILTERED_FILEBASE}_read2_paired.fastq 1> /dev/null 2> ${OUTPUT_FILEBASE}_bowtie_mapped.log
	gzip -f ${HOST_FILTERED_FILEBASE}_1.fastq
	gzip -f ${HOST_FILTERED_FILEBASE}_2.fastq
${HOST_FILTERED_FILEBASE}_2.fastq.gz: ${HOST_FILTERED_FILEBASE}_1.fastq.gz

# FASTQC report of host filtered file
${HOST_FILTERED_FILEBASE}_1_fastqc.zip: ${HOST_FILTERED_FILEBASE}_1.fastq.gz
	${FASTQC} "$<"
${HOST_FILTERED_FILEBASE}_2_fastqc.zip: ${HOST_FILTERED_FILEBASE}_2.fastq.gz
	${FASTQC} "$<"
# Add host filtering step to prereqs
${PREPROCESSING_SUMMARY}: ${HOST_FILTERED_FILEBASE}_1_fastqc.zip ${HOST_FILTERED_FILEBASE}_2_fastqc.zip
endif

# Summary
#########
${PREPROCESSING_SUMMARY}:
	${VAMP_HOME}/bin/preprocess.py ${OUTPUT_DIR} ${OUTPUT_FILEBASE} > "$@"
	
preprocess: ${PREPROCESSING_SUMMARY}
