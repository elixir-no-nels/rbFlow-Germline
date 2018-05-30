#!/bin/bash

bwa mem -t 4 \
	  -R "@RG\tID:ID\tSM:SM\tLB:LB\tPL:ILLUMINA\tPU:NotDefined" \
	  -M /home/oskar/01-workspace/01-data/refdata/hg19/human_g1k_v37_decoy.fasta \
	  /home/oskar/01-workspace/01-data/fastq/test_R1.fastq \
	  /home/oskar/01-workspace/01-data/fastq/test_R2.fastq \
    \
		|    \
    \
	  java -jar /home/oskar/01-workspace/00-temp/wdl_pipeline/tools/gatk4/gatk-package-4.beta.6-local.jar \
	  MergeSamFiles -I /dev/stdin -O mergesam.bam --USE_THREADING FALSE \
    \
		&&   \
    \
	  java -jar /home/oskar/01-workspace/00-temp/wdl_pipeline/tools/gatk4/gatk-package-4.beta.6-local.jar \
	  MarkDuplicates --input mergesam.bam \
	  -O markdup.bam \
	  --VALIDATION_STRINGENCY LENIENT \
	  --METRICS_FILE markdup.metrics \
	  --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 200000 \
	  --CREATE_INDEX true     \
    \
		&&   \
		\
  	  java -jar /home/oskar/01-workspace/00-temp/wdl_pipeline/tools/gatk4/gatk-package-4.beta.6-local.jar \
	  BaseRecalibrator \
	  --reference /home/oskar/01-workspace/01-data/refdata/hg19/human_g1k_v37_decoy.fasta \
	  --input markdup.bam \
	  -O baserecal.grp \
	  --knownSites /home/oskar/01-workspace/01-data/refdata/hg19/dbsnp.vcf
