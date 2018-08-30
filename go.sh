#!/bin/bash

echo "Clean (Some/All) previous results"
#rm -rf Outputs/01*
#rm -rf Outputs/02*
#rm -rf Outputs/03*
#rm -rf Outputs/04*
#rm -rf Outputs/05*
#rm -rf Outputs/06*
#rm -rf Outputs/07*

echo "Clean previous logs"
rm     Outputs/Log.log

echo "Run the Workflow (Classic)"
## Uncomment to run a basic test
#ruby wf_elixir_GermlineCaller_csv_singularity.rb -i Inputs/input.tsv -o Outputs -r test -l Outputs/Log.log -e singularity -t 16 -m 55

## Uncomment to run with na12878 test input files
ruby wf_elixir_GermlineCaller_csv_singularity.rb -i Inputs/input_na12878.tsv -o Outputs -r test -l Outputs/Log.log -e singularity -t 16 -m55

## Uncomment both lines to run with the bed file
#ruby wf_elixir_GermlineCaller_csv_singularity.rb -i Inputs/input.tsv -o Outputs -r test -l Outputs/Log.log -e singularity -t 16 -m 55 \
#-b Inputs/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed

## Don't uncomment this because it is untested ;)
#echo "Run the Workflow (Based on ReadsPipelineSpark)"
#ruby wf_elixir_GermlineCaller_csv_GATK-spark-singularity.rb -i Inputs/input.tsv -o Outputs -r ghis_test -l Outputs/Log.log
