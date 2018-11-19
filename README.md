## Instructions for how to use the rbFlow-Germline slurm workflow on TSD
### Short instructions  
Create a tsv file, add information about the samples with either the example file [here](https://github.com/elixir-no-nels/rbFlow-Germline/blob/master/Inputs/input.tsv) or copy/paste the example below into a new file:  
```bash
Flowcell_ID	SAMPLE_ID	LIB_ID	LANE	File_R1	File_R2
``` 
Fill the columns with appropriate information, then save the file and name it something good. Keep in mind that the columns must be tab separated. Don't use the space bar.  
Assuming you already have the input files ready, and that the output directory exists, you can now start the workflow.  
```bash 
./configure -i /input/data/directory -t /path/to/tsv-file -o /path/to/output-directory -r hg38
```  
This will use hg38 reference files, you can also use b37 reference files.

### Long instructions  
Let's begin with a thought experiment to understand how to supply the workflow with correct inputs:  
Your input data in this thought experiment is located in `/tsd/p11/username/input-data`, this directory has two files and one folder that also contains two files. The first flag that we can set based on this information is the `-i` flag, the `-i` flag takes the input file directory as argument, in this case it will look like this `./configure -i /tsd/p11/username/input-data [...]`.  
The next step is to prepare the tsv file.

**Preparing the tsv file**  
The tsv file is used to tell the workflow which files belong to the same sample. If you want to analyse many samples at once you simply put all information per sample and file in the tsv file. The sample information is also used to create proper read groups, this will help you identify the files since this information is put in the header of the bam files. This information is also used by the workflow to name the output files.  
Either use the sample file [here](https://github.com/elixir-no-nels/rbFlow-Germline/blob/master/Inputs/input.tsv) or copy/paste the example below into a new file:  
```bash
Flowcell_ID	SAMPLE_ID	LIB_ID	LANE	File_R1	File_R2
```  
If you don't have information about flowcell id, sample id, library id or lane number, you still need to put something there, in that case you will need to come up with something that makes sense to you. What matters most to get things to run is that you write the file names, or directory and file names, correctly. 

Let's continue the thought experiment from above and assume the two files that are located in the `/tsd/p11/username/input-data` directory are called `human_adenoma_R1.fastq.gz` and `human_adenoma_R2.fastq.gz`, simply put the file names in the last two columns like so: 
```bash
Flowcell_ID	SAMPLE_ID	LIB_ID	LANE	human_adenoma_R1.fastq.gz	human_adenoma_R2.fastq.gz
```
And let's also assume that the directory that is placed in the `/tsd/p11/username/input-data` directory is called `breast_cancer` and that the files are called `ductal_carcinoma_R1.fastq.gz` and `ductal_carcinoma_R2.fastq.gz`, the resulting columns in the tsv file would look like so:
```bash
Flowcell_ID	SAMPLE_ID	LIB_ID	LANE	human_adenoma_R1.fastq.gz	human_adenoma_R2.fastq.gz
Flowcell_ID	SAMPLE_ID	LIB_ID	LANE	breast_cancer/ductal_carcinoma_R1.fastq.gz	breast_cancer/ductal_carcinoma_R2.fastq.gz
```
Now that the tsv file has been written correctly you can add `-t /path/to/tsv/file.tsv` in the start command like so: 
```bash
./configure -i /tsd/p11/username/input-data -t /path/to/tsv/file.tsv [...]
```

Now let's move on the output directory.  

**Create the output directory**  
The output directory cannot be created by the workflow, therefore you need to create it manually before you start the workflow. Simply run:
```bash
mkdir /path/to/output/directory/on/the/shared/fileserver
```
and then add `-o /path/to/output/directory` to the command line like so: 
```bash
./configure -i /tsd/p11/username/input-data -t /path/to/tsv/file.tsv -o /path/to/output/directory/on/the/shared/fileserver [...]
```

**Selecting reference file version**  
You have a choice of two reference file versions, either the `b37` decoy version, or the `hg38` version.
This paper has some interesting comparisons of the b37/GRCh37 and hg38/GRCh38 reference files: https://arxiv.org/pdf/1404.0929.pdf

If you don't know which one to choose you should probably use hg38, it's genereally more complete compared to b37 according to the article above.

The flag is `-r`, so the resulting command line so far looks like so:  
```bash
./configure -i /tsd/p11/username/input-data -t /path/to/tsv/file.tsv -o /path/to/output/directory/on/the/shared/fileserver -r hg38
```

**Optional interval file**  
If you don't know what an interval file is you probably don't need one. But if you want to use one you need to be certain that the interval file is compatible with the default reference files.  
The hg38 reference files have chromosome names like this: `chr1`  
For a complete list of contigs in the hg38 reference fasta file you can check out this .dict file: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict  
In case it asks for password there is none.

The b37 reference files have chromosome names like this: `1`, i.e no `chr` before the chromosome number.  
For a complete list of contigs in the b37 reference fasta file you can check out this .dict file: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz  
In case it asks for password there is none.

There's a default hg38 interval list that you can use, there's not much documentation about it, from what I can gather it's tailored for the provided hg38 fasta file and seeks to exclude e.g centromeric regions since they don't add any useful information. Here's some more information: https://software.broadinstitute.org/gatk/documentation/article?id=11009  
Here's the most relevant part from that link:  
> **Whole genomes (WGS)**  
	For whole genome sequence, the intervals lists don’t depend on the prep (since in principle you captured the “whole genome”) so instead it depends on what regions of the genome you want to blacklist (e.g. centromeric regions that waste your time for nothing) and how the reference genome build enables you to cut up regions (separated by Ns) for scatter-gather parallelizing.  
	We make our WGS interval lists available, and the good news is that, as long as you're using the same genome reference build as us, you can use them with your own data even if it comes from somewhere else -- assuming you agree with our decisions about which regions to blacklist! Which you can examine by looking at the intervals themselves. However, we don't currently have documentation on their provenance, sorry -- baby steps.

If you want to use an interval list the flag is `-l` and you need provide a complete file path to it, for the default interval list it would look like this: `-l /put/the/actual/absolute/path/here/oskar.list`

And that's it! You should be able to run the workflow now, if you have any additional questions you can contact [Kjell Petersen](mailto:kjell.petersen@uib.no) or [Oskar Vidarsson](mailto:oskar.vidarsson@uib.no).  
The rest of the documentation covers technical nitty gritty things for systems developers. 

## Setting up the rbFlow-Germline slurm workflow from scratch

`git clone https://github.com/elixir-no-nels/rbFlow-Germline/`  
`cd rbFlow-Germline`



# rbFlow-Germline slurm start script explanation
The following text is intended to explain to a developer how the workflow is started.

## Short introduction
The `./configure` script creates variables that are exported to the `runOnNode.sbatch` script. These variables are in turn used to e.g copy files and set the output directory. Read the detailed description below for a more thorough explanation. 

### Detailed line by line description of `./configure` and `runOnNode.sbatch`
1. User runs `./configure`  
Mandatory flags are:  
	* -i for input file directory.  
	* -t for tab separated file with fastq info.  
	* -o for the output directory, directory must exist.  
	* -r for choosing reference file version, b37 or hg38.  
	Optional flag  
	* -l for interval file.  
	Other  
	* -h to print help message.

2. The `./configure` script checks that all required flags have been used.
3. The `./configure` script checks that all file paths are absolute.
4. When every flag has been added and the file paths are verified to be absolute the sbatch command is constructed.
5. The sbatch command is currently as such:  
`sbatch --export=INPUTS=$INPUTS,TSV=$TSV,OUTPUTDIR=$OUTPUTDIR,REFERENCE=$REFERENCE,INTERVAL=$INTERVAL runOnNode.sbatch`.  
The `runOnNode.sbatch` script receives the variables and here's a description of what they do:
	* INPUTS is used to copy the input files in the INPUTS directory to the `/scratch/rbFlowGermline/Inputs` directory.
	* TSV and REFERENCE puts the path to the tsv file and reference file choices in the start command for rbFlow, i.e `ruby wf_elixir_GermlineCaller_csv_singularity.rb -i $TSV -o Outputs -r $REFERENCE -l Outputs/Log.log -e singularity -t 16 -m 55 `.
	* INTERVAL is also added to the start command if it has been set, i.e `ruby wf_elixir_GermlineCaller_csv_singularity.rb -i $TSV -o Outputs -r $REFERENCE -l Outputs/Log.log -e singularity -t 16 -m 55 -b $INTERVAL`. See the explanation for the runOnNode.sbatch script for information about the if, then, else statement that runs the workflow either with or without an interval file.
	* Then the `runOnNode.sbatch` script is started.
6. The `runOnNode.sbatch` script does the following
	* First there are hardcoded settings for slurm to set the following variables.
		* job name: rbFlow-Germline.
		* account name: p172.
		* time: 0-50 hours.
		* mem-per-cpu: 3750.
		* cpus-per-task: 16.
		* This reserves a 60GB, 20 core execute node but only 16 cores are actually used.
		* HaplotypeCaller risks crashing if more cores are used, so this is as good as it gets.
	* Two modules are loaded: ruby (default version) and singularity (version 2.5.0).
	* The rbFlowGermline directory path is put in a variable.
	* This variable is used to create a basename.
	* This variable is also used to copy the rbFlowGermline directory to the scratch disk.
	* The input file directory is then copied to the input file directory in the rbFlowGermline directory on the scratch disk.
	* `ls /scratch/rbflow/Inputs` is executed only for troubleshooting purposes, this line can be deleted.
	* The last step before the pipeline is actually started is an if statement that checks if the INTERVAL variable was set and executes the appropriate command based on the result.
	* The last step is `cp -r /scratch/rbflow/outputs/* $OUTPUTDIR` and is executed no matter the exit status of the workflow, so any and all contents of the output directory are transfered back. 

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