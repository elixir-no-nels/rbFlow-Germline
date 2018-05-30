#!/usr/bin/env ruby
################################################################################

################################################################################


Process.setproctitle('workflow') # Change the name of the process
require 'optparse'

$LOAD_PATH.unshift("#{File.dirname(__FILE__)}/")
require 'rbFlowLite/rbFlowLite.rb'
require 'rbFlowLite/setup.rb'
require 'rbFlowLite/toolbox_elixir_GermlineCaller_singularity.rb'

# Help to debug
DEBUGINFO=false
if DEBUGINFO
  puts "\n\n!!!!------ LOAD AWESOME_PRINT ------!!!!\n\n"
  require 'awesome_print'
end

## Arguments Parser

# Default values
options = {}
options[:inputs_csv]         = false
options[:ouputs_folder]      = false
options[:genome_version]     = false
options[:intervals_file]     = false
options[:wf_log]             = false
options[:set_threads]        = false
options[:set_ram]            = false
options[:engine]             = false

# Help message
def usage(options)
  puts 'Application Usage : '
  puts "-i    --inputs_csv         Directory containing inputs                    Default Value: #{options[:inputs_csv]}"
  puts "-o    --ouputs_folder      Directory where to write outputs               Default Value: #{options[:ouputs_folder]}"
  puts "-r    --genome_version     Directory containing references                Default Value: #{options[:genome_version]}"
  puts "-b    --intervals_file     Intervals file                                 Default Value: #{options[:intervals_file]}"
  puts "-l    --workflow_log       logs file for the workflow (stdout by default) Default Value: #{options[:wf_log]}"
  puts "-t    --set_threads        setup a custom number of threads to use        Default Value: #{options[:set_threads]}"
  puts "-m    --set_ram            setup a custom amount of memory for GATK       Default Value: #{options[:set_ram]}"
  puts "-e    --set_engine         setup a container execution engine             Default Value: #{options[:set_engine]}"
  puts "\n\n"
  puts 'TSV file format for the input samples list : text file format with unix return style (LF), tab separated fields, use # to comment a line :'
  puts 'flowcell_id    sample_id    library_id    lane    path/file_R1    path/file_R2'
  puts 'flowcell_id    sample_id    library_id    lane    path/file_R1    path/file_R2'
  puts '...'
  puts "\n\n"
  exit
  exit
end

# Args
OptionParser.new do |opts|
  opts.banner = 'Usage: example.rb [options]'

  opts.on('-i', '--inputs_csv folder', String, "Directory containing inputs: #{options[:inputs_csv]}") do |arg|
    options[:inputs_csv] = arg
  end

  opts.on('-o', '--ouputs_folder folder', String, "Directory where to write outputs: #{options[:ouputs_folder]}") do |arg|
    options[:ouputs_folder] = arg
  end

  opts.on('-r', '--genome_version name', String, "Name of the reference genome to use: #{options[:genome_version]}") do |arg|
    options[:genome_version] = arg
  end

  opts.on('-b', '--intervals_file file', String, "intervals file : #{options[:intervals_file]}") do |arg|
    options[:intervals_file] = arg
  end

  opts.on('-l', '--log file', String, "Log file (stdout if false) : #{options[:wf_log]}") do |arg|
    options[:log] = arg
  end

  opts.on('-t', '--set_threads threads', String, "setup a custom number of threads to use : #{options[:set_threads]}") do |arg|
    options[:set_threads] = arg
  end

  opts.on('-m', '--set_ram ram', String, "setup a custom amount of memory for GATK : #{options[:set_ram]}") do |arg|
    options[:set_ram] = arg
  end

  opts.on('-e', '--set_engine container_engine', String, "setup a container execution engine : #{options[:engine]}") do |arg|
    options[:engine] = arg
  end

  opts.on('-h', '--help', 'Help') do |arg|
    usage(options) if arg
  end

end.parse!


## ---- Setup the Workflow ----

# if one of this ptions is false print help and exit
error_list = []
error_list.push 'provide an input csv file'             if not options[:inputs_csv]
error_list.push 'provide an output'                     if not options[:ouputs_folder]
error_list.push 'provide a reference'                   if not options[:genome_version]
if error_list.size > 0
  spacer  = "\n" + ' ' * 48 + '- '
  message = error_list.join(spacer).to_s
  Logger.error message: "Error missing Args :#{spacer}#{message}", logfile: false
  puts
  usage options
end

# ---- Variables ----
inputs_csv         = File.expand_path options[:inputs_csv]
ouputs_folder      = File.expand_path options[:ouputs_folder]
genome_version     = options[:genome_version]
container_repo     = false
container_repo     = File.expand_path options[:singularity_repo] if options[:singularity_repo]
intervals_file     = false
intervals_file     = File.expand_path options[:intervals_file]   if options[:intervals_file]

# Setup Dataset and engine
references         = set_reference(genome_version: genome_version, intervals_file: intervals_file)
engine             = options[:engine]      if options[:engine]      # shell, docker , singularity...
container_setup    = set_engine(engine: engine)
config             = set_resources(type: :laptop)
config[:cpu]       = options[:set_threads] if options[:set_threads] # overide default threads value
config[:ram]       = options[:set_ram]     if options[:set_ram]     # overide default ram value
# ---- vars ----
#@ProjectDir        = '/tsd/'
#@Tmp               = '/tmp'
wf_log             = File.expand_path options[:log] if options[:log]


# ---- Create path for logs ----
FileUtils.mkdir_p File.dirname options[:log] if options[:log]

# ---- Print options ----
Logger.info message: '================  Germline variant calling Workflow Starting ================', logfile: wf_log
Logger.info message: "Parameters\tinput_csv           \t#{options[:inputs_csv]}",     logfile: wf_log
Logger.info message: "Parameters\toutput_folder       \t#{options[:ouputs_folder]}/", logfile: wf_log
Logger.info message: "Parameters\treference version   \t#{options[:genome_version]}", logfile: wf_log
Logger.info message: "Parameters\tintervals_file      \t#{intervals_file}",           logfile: wf_log
Logger.info message: "Parameters\tresource type       \t#{config}",                   logfile: wf_log
Logger.info message: "Parameters\tLog                 \t#{:wf_log}",                  logfile: wf_log
Logger.info message: "Parameters\tEngine              \t#{engine}",                   logfile: wf_log
Logger.info message: "Parameters\tcontainers config   \t#{container_setup}",          logfile: wf_log
Logger.info message: '==============================================================================', logfile: wf_log


# ---- Start the Workflow ----

## ---- Preprocessing ----

# Mapping
input  = inputs_csv
output = ouputs_folder + '/01_mapping'
mapping_from_list(input_csv: input, output_dir: output, reference: references, config: config, container: container_setup, wf_log: wf_log)

# Mergin/Sorting
input  = ouputs_folder + '/01_mapping'
output = ouputs_folder + '/02_sorting'
merge_sort(input_dir: input, output_dir: output, config: config, group_by_samples: true ,container: container_setup, wf_log: wf_log)

# GATK4 ReadsPipelineSpark **Beta**
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.6/org_broadinstitute_hellbender_tools_spark_pipelines_ReadsPipelineSpark.php
input  = ouputs_folder + '/02_sorting'
output = ouputs_folder + '/03_ReadsPipelineSpark'
gatk_ReadsPipelineSpark(input_dir: input, output_dir: output, reference: references, config: config, container: container_setup, wf_log: wf_log)

# GenotypeGVCF
input  = ouputs_folder + '/03_ReadsPipelineSpark'
output = ouputs_folder + '/04_genotype_GVCF'
gatk_genotype_GVCF(input_dir: input, output_dir: output, reference: references, config: config, container: container_setup, wf_log: wf_log)

# Variants recalibration SNP
input  = ouputs_folder + '/04_genotype_GVCF'
output = ouputs_folder + '/05_variant_recalibrator_applied_snp'
gatk_variant_recalibrator_snp(input_dir: input, output_dir: output, reference: references, config: config, container: container_setup, wf_log: wf_log)

# Variants recallibration INDELS
input  = ouputs_folder + '/04_genotype_GVCF'
output = ouputs_folder + '/05_variant_recalibrator_applied_indel'
gatk_variant_recalibrator_indel(input_dir: input, output_dir: output, reference: references, config: config, container: container_setup, wf_log: wf_log)
