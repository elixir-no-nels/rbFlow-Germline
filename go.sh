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
ruby wf_elixir_GermlineCaller_csv_singularity.rb -i Inputs/input_na12878.tsv -o Outputs -r test -l Outputs/Log.log -e singularity -t 16 -m 16

#echo "Run the Workflow (Based on ReadsPipelineSpark)"
#ruby wf_elixir_GermlineCaller_csv_GATK-spark-singularity.rb -i Inputs/input.tsv -o Outputs -r ghis_test -l Outputs/Log.log
