require 'pathname'
require 'fileutils'

# GATK setup (DRY)
GATK_JAR = '/opt/conda/share/gatk4-4.0.8.1-0/gatk-package-4.0.8.1-local.jar'
BWA_IMG = 'Tools.simg'
GATK_IMG = 'Tools.simg'

## BWA_MEM Mapping with CSV input list
def mapping_from_list(input_csv: '', output_dir: '', add_suffix: '', reference: {}, config: {}, container: {container: false}, wf_log: false)
  cpu_num = config[:cpu]
  FileUtils.mkdir_p(output_dir)
  Logger.info message: "bwa mem mapping from #{input_csv} list", logfile: wf_log
  File.open(input_csv, "r") do |f|
    f.each do |line|
      line.chomp!
      # Some checks
      # check the number of columns
      split_char        = "\t"
      line_size         = line.split(split_char).size
      expected_field    = 6
      next if line_size == 0    # skip empty lines
      next if line[0]   == '#'  # skip line starting by "#""
      if line_size != expected_field
        Logger.error message: "Error on #{input_csv} parsing", logfile: wf_log
        Logger.error message: "a line contain #{line_size} field, #{expected_field} are expected :", logfile: wf_log
        Logger.error message: "#{line}", logfile: wf_log
        exit()
      end
      # Check if there is non allowed characters (no spaces,quotes or backslash)
      if not line.match(/[\"\'\\]/).nil?
        Logger.error message: "Error on #{input_csv} parsing", logfile: wf_log
        Logger.error message: "a line contain a non allowed symbol (quote, backslash) :", logfile: wf_log
        Logger.error message: "#{line}", logfile: wf_log
        exit()
      end
      # Build ReadGroups (based on : http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)
      # FLOWCELL_ID,RGSM,RGLB,Lane,File_R1,File_R2
      # RGID : the read group ID, must be unique flowcell_barcode.lane is a good ID
      # RGSM : the sample, to identify a sample splited on several files/sequencing Run
      # RGLB : Library ID, used for mark/remove duplicate.
      # RGPU : {flowcell_barcode}.{lane}.{sample_barcode}
      flowcell_id, sample, library_id, lane, file_name_R1, file_name_R2 = line.split(split_char)
      rgid = flowcell_id + '.' + lane
      rgpu = flowcell_id + '.' + lane + '.' + sample
      rgsm = sample
      rglb = library_id
      readgroup_string = "-R \"@RG\\tID:#{rgid}\\tSM:#{rgsm}\\tLB:#{rglb}\\tPL:ILLUMINA\\tPU:#{rgpu}\""
      binary = 'bwa'
      # Mapping
      base_name    = "#{sample}_#{lane}_#{flowcell_id}"
      output_file  = base_name + add_suffix + '.bam'
      skip_test    = File.exist?("#{output_dir}/#{output_file}")
      step         = Runner.new task_name: "BWA_MEM_map_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
      full_name_R1 = File.expand_path file_name_R1
      full_name_R2 = File.expand_path file_name_R2
      step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{BWA_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
      step.cmd_line << binary
      step.cmd_line << 'mem'
      step.cmd_line << "-t #{cpu_num}"
      step.cmd_line << readgroup_string
      step.cmd_line << '-M'
      step.cmd_line << "#{reference[:genome_fas]}"
      step.cmd_line << "#{full_name_R1}"
      step.cmd_line << "#{full_name_R2}"
      step.cmd_line << "| samtools view -bS -"
      step.cmd_line << "> #{output_dir}/#{output_file}"
      step.run skip: skip_test, debug: false
    end
  end
end

## GATK4 merge all input for each sample in a sorted Bam file
def merge_sort(input_dir: '', input_filter: '*', output_dir: '', add_suffix: '', reference: {}, group_by_samples: false, config: {}, container: {container: false}, wf_log: false)
  ram            = config[:ram]
  task_name      = "GATK_MergeSortSam"
  group_spliter  = '_' # Split file names on this character
  group_on_index = 0   # Sample are on the first par of the splited file name
  tmp_dir        = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK merge and sort', logfile: wf_log
  file_list   = get_file_names path: input_dir, filter: "*#{input_filter}*.bam", full_path: true
  # Get a list of samples
  groups_list = ['*']
  if group_by_samples
    Logger.info message: " Regroup file by Samples, sample name on index #{group_on_index} after spliting file names with #{group_spliter}", logfile: wf_log
    groups_list = []
    file_list.each do |file|
      filename = file.split('/')[-1] # Get the file name only
      groupname = filename.split(group_spliter)[group_on_index]
      groups_list.push groupname if not groups_list.include? groupname
    end
    Logger.info message: " Group by samples: #{groups_list.size} samples found", logfile: wf_log
  end
  error_list.push 'Cannot find sample tags to group files' if groups_list.size == 0
  # For each sample merge all files
  groups_list.each do |sample_tag|
    inputs_list = []
    if sample_tag == '*'
      inputs_list = file_list
    else
      file_list.each do |file|
        inputs_list.push file if file.include? sample_tag
      end
    end
    Logger.info message: " Sample: #{sample_tag} : #{inputs_list.size} files found", logfile: wf_log
    bin_cmd = []
    bin_cmd << 'java'
    bin_cmd << "-Xmx#{ram}G"
    bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
    bin_cmd << "-jar #{GATK_JAR}"
    binary  = bin_cmd.join(' ').to_s
    inputs  = []
    for file in inputs_list
      inputs << "--INPUT #{file}"
    end
    base_name   = sample_tag
    base_name   = 'mapped_sequences' if sample_tag == '*'
    output_file = base_name + add_suffix + '.bam'
    skip_test   = File.exist?("#{output_dir}/#{output_file}")
    step        = Runner.new task_name: "#{task_name}_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << 'MergeSamFiles'
    step.cmd_line << inputs.join(' ').to_s
    step.cmd_line << "--OUTPUT=#{output_dir}/#{output_file}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << '--USE_THREADING=true'
    step.cmd_line << '--CREATE_INDEX=true'
    step.cmd_line << '--MAX_RECORDS_IN_RAM=1000000'
    step.cmd_line << '--SORT_ORDER=coordinate'
    step.cmd_line << '--VALIDATION_STRINGENCY=LENIENT'
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end


## GATK MarkDuplicate
def mark_duplicates( input_dir: '', input_filter: '*', output_dir: '', config: {}, container: {container: false}, wf_log: false)
  ram     = config[:ram]
  tmp_dir = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK Mark Duplicate', logfile: wf_log
  file_list = get_file_names path: input_dir, filter: "*#{input_filter}*.bam", full_path: true, wf_log: wf_log
  bin_cmd = []
  bin_cmd << 'java'
  bin_cmd << "-Xmx#{ram}G"
  bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
  bin_cmd << "-jar #{GATK_JAR}"
  binary  = bin_cmd.join(' ').to_s
  for file in file_list
    base_name   = File.basename file, '.bam'
    output_file = base_name + '_dedup' + '.bam'
    metric_file = base_name + '_dedup' + '.metric'
    skip_test   = File.exist?("#{output_dir}/#{output_file}")
    step        = Runner.new task_name: "GATK_MarkDuplicate_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << 'MarkDuplicates'
    step.cmd_line << "--INPUT=#{file}"
    step.cmd_line << "--OUTPUT=#{output_dir}/#{output_file}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << "--METRICS_FILE=#{output_dir}/#{metric_file}"
    step.cmd_line << '--CREATE_INDEX=true'
    step.cmd_line << '--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200000'
    step.cmd_line << '--VALIDATION_STRINGENCY=LENIENT'
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end


## GATK BaseRecalibrator
def gatk_base_recalibrator(input_dir: '', input_filter: '*', output_dir: '', reference: {}, config: {}, container: {container: false}, wf_log: false)
  cpu_num = config[:cpu]
  ram     = config[:ram]
  tmp_dir = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK base recalibration', logfile: wf_log
  file_list = get_file_names path: input_dir, filter: "*#{input_filter}*.bam", full_path: true, wf_log: wf_log
  bin_cmd = []
  bin_cmd << 'java'
  bin_cmd << "-Xmx#{ram}G"
  bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
  bin_cmd << "-jar #{GATK_JAR}"
  binary  = bin_cmd.join(' ').to_s
  interval_opt = ''
  if reference[:intervals] != false
    interval_opt  = "--intervals #{reference[:intervals]}"
  end
  for file in file_list
    # Recalibration
    base_name  = File.basename file, '.bam'
    group_file = base_name + '.group'
    skip_test  = File.exist?("#{output_dir}/#{group_file}")
    step       = Runner.new task_name: "GATK_Recalib_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    if config[:spark]
      step.cmd_line << 'BaseRecalibratorSpark'
      step.cmd_line << "--spark-master local[#{cpu_num}]"
      step.cmd_line << "-R #{reference[:genome_2bit]}"
    else
      step.cmd_line << 'BaseRecalibrator'
      step.cmd_line << "-R #{reference[:genome_fas]}"
    end
    step.cmd_line << "-I #{file}"
    step.cmd_line << "-O #{output_dir}/#{group_file}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << "--known-sites #{reference[:kgenomes]} "
    step.cmd_line << "--known-sites #{reference[:mills]} "
    step.cmd_line << "--known-sites #{reference[:dbsnp]} "
    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
    #ApplyBQSR
    output_file = base_name + '.bam'
    skip_test   = File.exist?("#{output_dir}/#{output_file}")
    step        = Runner.new task_name: "GATK_ApplyBQSR_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    if config[:spark]
      step.cmd_line << 'ApplyBQSRSpark'
      step.cmd_line << "--spark-master local[#{cpu_num}]"
      step.cmd_line << "-R #{reference[:genome_2bit]}"
    else
      step.cmd_line << 'ApplyBQSR'
      step.cmd_line << "-R #{reference[:genome_fas]}"
    end
    step.cmd_line << "-I #{file}"
    step.cmd_line << "-O #{output_dir}/#{output_file}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << "-bqsr #{output_dir}/#{group_file}"
    step.cmd_line << interval_opt
    #step.cmd_line << '--createOutputBamIndex true'
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end


## GATK_HaplotypeCaller
def gatk_haplotype_caller(input_dir: '', input_filter: '*', output_dir: '', reference: {}, config: {}, container: {container: false}, wf_log: false)
  cpu_num = config[:cpu]
  ram     = config[:ram]
  tmp_dir = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK Haplotype Caller', logfile: wf_log
  file_list = get_file_names path: input_dir, filter: '*.bam', full_path: true, wf_log: wf_log
  bin_cmd = []
  bin_cmd << 'java'
  bin_cmd << "-Xmx#{ram}G"
  bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
  bin_cmd << "-jar #{GATK_JAR}"
  binary  = bin_cmd.join(' ').to_s
  interval_opt = ''
  if reference[:intervals] != false
    interval_opt  = "--intervals #{reference[:intervals]}"
  end
  for file in file_list
    base_name = File.basename file, '.bam'
    output    = base_name + '_haplotype' + '.g.vcf'
    skip_test = File.exist?("#{output_dir}/#{output}")
    step      = Runner.new task_name: "GATK_HaplotypeCaller_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    if config[:spark]
      step.cmd_line << 'HaplotypeCallerSpark'
      step.cmd_line << "--spark-master local[#{cpu_num = config[:cpu]}]"
      step.cmd_line << "-R #{reference[:genome_2bit]}"
    else
      step.cmd_line << 'HaplotypeCaller'
      step.cmd_line << "-R #{reference[:genome_fas]}"
    end
    step.cmd_line << "-I #{file}"
    step.cmd_line << "--emit-ref-confidence  GVCF"
    #step.cmd_line << "-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation" # Allele specific
    step.cmd_line << "-O #{output_dir}/#{output}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end


## GATK_GenotypeGVCFs
def gatk_genotype_GVCF(input_dir: '', input_filter: '*', output_dir: '', reference: {}, config: {}, container: {container: false}, wf_log: false)
  cpu_num = config[:cpu]
  ram     = config[:ram]
  tmp_dir = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK GenotypeGVCFs', logfile: wf_log
  file_list = get_file_names path: input_dir, filter: '*.g.vcf', full_path: true, wf_log: wf_log
  bin_cmd = []
  bin_cmd << 'java'
  bin_cmd << "-Xmx#{ram}G"
  bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
  bin_cmd << "-jar #{GATK_JAR}"
  binary  = bin_cmd.join(' ').to_s
  interval_opt = ''
  if reference[:intervals] != false
    interval_opt  = "--intervals #{reference[:intervals]}"
  end
  for file in file_list
    base_name = File.basename file, '.g.vcf'
    output    = base_name + "_genotype" + '.g.vcf'
    skip_test = File.exist?("#{output_dir}/#{output}")
    step      = Runner.new task_name: 'GATK_GenotypeGVCFs', wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << 'GenotypeGVCFs'
    step.cmd_line << "-R #{reference[:genome_fas]}"
    step.cmd_line << "--variant #{file}"
    step.cmd_line << "-O #{output_dir}/#{output}"
    #step.cmd_line << "-G Standard -G AS_Standard" # Allele specific
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end


## GATK_VariantRecalibrator INDEL Mode
def gatk_variant_recalibrator_indel(input_dir: '', input_filter: '*', output_dir: '', reference: {}, config: {}, container: {container: false}, wf_log: false)
  ram     = config[:ram]
  tmp_dir = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK VariantRecalibrator INDEL / ApplyRecalibration', logfile: wf_log
  file_list = get_file_names path: input_dir, filter: '*.g.vcf', full_path: true, wf_log: wf_log
  bin_cmd = []
  bin_cmd << 'java'
  bin_cmd << "-Xmx#{ram}G"
  bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
  bin_cmd << "-jar #{GATK_JAR}"
  binary  = bin_cmd.join(' ').to_s
  interval_opt = ''
  if reference[:intervals] != false
    interval_opt  = "--intervals #{reference[:intervals]}"
  end
  for file in file_list
    base_name   = File.basename file, '.g.vcf'
    tranche_out = base_name + '_indel_tranche' + '.tranches'
    recall_out  = base_name + '_indel_recall'  + '.recall'
    skip_test   = File.exist?("#{output_dir}/#{tranche_out}")
    skip_test   = File.exist?("#{output_dir}/#{recall_out}")
    step        = Runner.new task_name: "GATK_VariantRecalibrator_INDEL_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << ' VariantRecalibrator'
    step.cmd_line << "-R #{reference[:genome_fas]}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << "-V #{file}"
    step.cmd_line << '-mode INDEL'
    step.cmd_line << "--tranches-file #{output_dir}/#{tranche_out}"
    step.cmd_line << "--output        #{output_dir}/#{recall_out}"
	  step.cmd_line << "--resource mills,known=false,training=true,truth=true,prior=12.0:#{reference[:mills]} "
	  step.cmd_line << "--resource dbsnp,known=true,training=false,truth=false,prior=2.0:#{reference[:dbsnp]}"
	  step.cmd_line << '-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum '
	  step.cmd_line << '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0'
	  step.cmd_line << '-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5'
	  step.cmd_line << '-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0'
    step.cmd_line << '--max-gaussians 4'
    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
    #ApplyRecalibration
    output_file = base_name + '.g.vcf'
    skip_test   = File.exist?("#{output_dir}/#{output_file}")
    step        = Runner.new task_name: "GATK_ApplyVQSR_INDEL_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << 'ApplyVQSR'
    step.cmd_line << "-R #{reference[:genome_fas]}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << "-V #{file}"
    step.cmd_line << "--output #{output_dir}/#{output_file}"
    step.cmd_line << '-mode INDEL'
    step.cmd_line << "--tranches-file #{output_dir}/#{tranche_out}"
    step.cmd_line << "--recal-file    #{output_dir}/#{recall_out}"
    step.cmd_line << '--ts-filter-level 95.0 '

    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end

## GATK_VariantRecalibrator SNP Mode
def gatk_variant_recalibrator_snp(input_dir: '', input_filter: '*', output_dir: '', reference: {}, config: {}, container: {container: false}, wf_log: false)
  ram     = config[:ram]
  tmp_dir = output_dir + '/tmp'
  FileUtils.mkdir_p(output_dir)
  FileUtils.mkdir_p(tmp_dir)
  Logger.info message: 'GATK VariantRecalibrator SNP / ApplyRecalibration', logfile: wf_log
  file_list = get_file_names path: input_dir, filter: '*.g.vcf', full_path: true, wf_log: wf_log
  bin_cmd = []
  bin_cmd << 'java'
  bin_cmd << "-Xmx#{ram}G"
  bin_cmd << "-Djava.io.tempdir=#{tmp_dir}"
  bin_cmd << "-jar #{GATK_JAR}"
  binary  = bin_cmd.join(' ').to_s
  interval_opt = ''
  if reference[:intervals] != false
    interval_opt  = "--intervals #{reference[:intervals]}"
  end
  for file in file_list
    base_name   = File.basename file, '.g.vcf'
    tranche_out = base_name + '_snp_tranche' + '.tranches'
    recall_out  = base_name + '_snp_recall'  + '.recall'
    skip_test   = File.exist?("#{output_dir}/#{tranche_out}")
    skip_test   = File.exist?("#{output_dir}/#{recall_out}")
    step        = Runner.new task_name: "GATK_VariantRecalibrator_SNP_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << 'VariantRecalibrator'
    step.cmd_line << "-R #{reference[:genome_fas]}"
    step.cmd_line << "-V #{file}"
    step.cmd_line << '-mode SNP'
    step.cmd_line << "--tranches-file #{output_dir}/#{tranche_out}"
    step.cmd_line << "--output        #{output_dir}/#{recall_out}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
	  step.cmd_line << "--resource v1000G,known=false,training=true,truth=false,prior=10.0:#{reference[:kgenomes]} "
	  step.cmd_line << "--resource omni,known=false,training=true,truth=true,prior=12.0:#{reference[:omni]}"
	  step.cmd_line << "--resource dbsnp,known=true,training=false,truth=false,prior=2.0:#{reference[:dbsnp]}"
	  step.cmd_line << "--resource hapmap,known=false,training=true,truth=true,prior=15.0:#{reference[:hapmap]}"
	  step.cmd_line << '-an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR'
	  step.cmd_line << '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6'
	  step.cmd_line << '-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0'
	  step.cmd_line << '-tranche 97.0 -tranche 90.0'
    step.cmd_line << '--max-gaussians 4'
    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
    #ApplyRecalibration
    output_file = base_name + '.g.vcf'
    skip_test   = File.exist?("#{output_dir}/#{output_file}")
    step        = Runner.new task_name: "GATK_ApplyVQSR_SNP_#{base_name}", wf_log: wf_log, task_log_dir: output_dir
    step.in_container(engine: container[:type], imgURL: "#{container[:repo]}#{GATK_IMG}", additional_options: container[:optional_args], mountList: container[:mountList]) if container[:container]
    step.cmd_line << binary
    step.cmd_line << 'ApplyVQSR'
    step.cmd_line << "-R #{reference[:genome_fas]}"
    step.cmd_line << "--TMP_DIR #{tmp_dir}"
    step.cmd_line << "-V #{file}"
    step.cmd_line << "--output #{output_dir}/#{output_file}"
    step.cmd_line << '-mode SNP'
    step.cmd_line << "--tranches-file #{output_dir}/#{tranche_out}"
    step.cmd_line << "--recal-file    #{output_dir}/#{recall_out}"
    step.cmd_line << '--ts-filter-level 99.6'
    step.cmd_line << interval_opt
    step.run skip: skip_test, debug: false
  end
  FileUtils.rm_rf(tmp_dir)
end
