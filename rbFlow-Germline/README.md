# Inputs

Readgroups will be build following information provided in the input file.
The merging/sorting step will merge all file from a same sample, resulting a single bam file for each sample.

Files descriptor :

```bash
Flowcell_ID	SAMPLE_ID	LIB_ID	LANE	File_R1	File_R2
```

# Classic Germline calling pipeline workflow based on the GATK best practice

## Preprocessing

### Mapping

```bash
bwa mem -R "<READGROUPS>" -M <REFERENCE> <R1.fastq> <R2.fastq> | samtools view -bS - > <OUTPUT.BAM>
```

### Merging/Sorting

```bash
GATK4 MergeSamFiles --INPUTS=<ALL_BAM_FROM_SAMPLE> --OUTPUT=<BAMFILE_1> --USE_THREADING=true --CREATE_INDEX=true --MAX_RECORDS_IN_RAM=1000000 --SORT_ORDER=coordinate --VALIDATION_STRINGENCY=LENIENT
```

### MarkDuplicates

```bash
GATK4 MarkDuplicates --INPUT=<BAMFILE_1> --OUTPUT=<BAMFILE_2> --METRICS_FILE=<METRICS> --CREATE_INDEX=true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200000 --VALIDATION_STRINGENCY=LENIENT
```

### Recalibration (spark)

```bash
GATK4 BaseRecalibratorSpark -R <REFERENCES.2bit> -I <BAMFILE_2> -O <GROUPS> --known-sites <1000genomes> --known-sites <mills> --known-sites <dbsnp> -L <INTERVAL>

GATK4 ApplyBQSRSpark -R <REFERENCES.2bit> -I <BAMFILE_2> -O <BAMFILE_3> -bqsr <GROUPS>
```

## Calling variants

### Haplotype caller (spark)

```bash
GATK4 HaplotypeCallerSpark -R <REFERENCES.2bit> -I <BAMFILE_3> --emit-ref-confidence GVCF" -O <G.VCF_1>
```

### GenotypeGVCF

```bash
GATK4 GenotypeGVCFs -R <REFERENCES> --variant <G.VCF_1> -O <G.VCF_2> -L <INTERVALS>
```

### Variants recalibration SNP

```bash
GATK4 VariantRecalibrator -R <REFERENCE> -V <G.VCF_2> -mode SNP --tranches-file <TRANCHES> --output <G.VCF_3> --resource v1000G,known=false,training=true,truth=false,prior=10.0:1000genomes --resource omni,known=false,training=true,truth=true,prior=12.0:omni --resource dbsnp,known=true,training=false,truth=false,prior=2.0:dbsnp --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap -an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -L <INTERVALS>

GATK4  ApplyVQSR -R <REFERENCES> -V <G.VCF_3> --output <G.VCF_5> -mode SNP --tranches-file <TRANCHES> --recal_file <RECALL> --ts_filter_level 99.6
```

### Variants recalibration INDELS

```bash
GATK4  VariantRecalibrator -R <REFERENCE> -V <G.VCF_2> -mode INDEL --tranches-file <TRANCHES> --output <G.VCF_3> --resource mills,known=false,training=true,truth=true,prior=12.0:mills --resource dbsnp,known=true,training=false,truth=false,prior=2.0:dbsnp -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -mG 4 -L INTERVAL

GATK4 -R <REFERENCE> -V <G.VCF_3> --output <G.VCF_4> -mode INDEL --tranches-file <TRANCHES> --recal_file <RECALL> --ts_filter_level 95.0
```

# Germline calling pipeline workflow based on the GATK best practice using ReadsPipelineSpark All in one

## Preprocessing

### Mapping

```bash
bwa mem -R "<READGROUPS>" -M <REFERENCE> <R1.fastq> <R2.fastq> | samtools view -bS - > <OUTPUT.BAM>
```

### Merging/Sorting

```bash
GATK4 MergeSamFiles --INPUTS=<ALL_BAM_FROM_SAMPLE> --OUTPUT=<BAMFILE_1> --USE_THREADING=true --CREATE_INDEX=true --MAX_RECORDS_IN_RAM=1000000 --SORT_ORDER=coordinate --VALIDATION_STRINGENCY=LENIENT'
```

### ReadsPipelineSpark

```bash
GATK4 ReadsPipelineSpark -I <BAMFILE_1> -O <G.VCF_1> --outputBam <BAMFILE_2> --emit-ref-confidence GVCF --known-sites 1000genomes --known-sites mills --known-sites dbsnp -L <INTERVAL>
```

## Calling variants

### GenotypeGVCF

```bash
GATK4 GenotypeGVCFs -R <REFERENCES> --variant <G.VCF_1> -O <G.VCF_2> -L <INTERVALS>
```

### Variants recalibration SNP

```bash
GATK4 VariantRecalibrator -R <REFERENCE> -V <G.VCF_2> -mode SNP --tranches-file <TRANCHES> --output <G.VCF_3> --resource v1000G,known=false,training=true,truth=false,prior=10.0:1000genomes --resource omni,known=false,training=true,truth=true,prior=12.0:omni --resource dbsnp,known=true,training=false,truth=false,prior=2.0:dbsnp --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap -an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -L <INTERVALS>

GATK4  ApplyVQSR -R <REFERENCES> -V <G.VCF_3> --output <G.VCF_5> -mode SNP --tranches-file <TRANCHES> --recal_file <RECALL> --ts_filter_level 99.6
```

### Variants recalibration INDELS

```bash
GATK4  VariantRecalibrator -R <REFERENCE> -V <G.VCF_2> -mode INDEL --tranches-file <TRANCHES> --output <G.VCF_3> --resource mills,known=false,training=true,truth=true,prior=12.0:mills --resource dbsnp,known=true,training=false,truth=false,prior=2.0:dbsnp -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -mG 4 -L <INTERVAL>

GATK4 -R <REFERENCE> -V <G.VCF_3> --output <G.VCF_4> -mode INDEL --tranches-file <TRANCHES> --recal_file <RECALL> --ts_filter_level 95.0
```

# Spark parallel processing

The parallel processing can be activated on Spark complient tools by using a local or distributed Spark cluster.
The local Spark doesn't require to start a specific Spark master server.

```bash
--spark-master local[<PARALLEL_TASKS>]
```