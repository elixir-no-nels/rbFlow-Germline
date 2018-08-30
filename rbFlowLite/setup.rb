# ---- Engines ----
# Structure containing information needed to run a specific engine
def set_engine(engine: false)
  container = {}
  if engine == 'singularity'
    container = {
      container:  true,
      type:      :singularity,
      repo:      '/media/01-workspace/04-pipelines/rbFlow-Germline/Singularity/',
      mountList: ['/media'],
      optional_args: []
    }
  elsif engine == 'singularity_TSD_VM'
    container = {
      container:  true,
      type:      :singularity,
      repo:      '/cluster/projects/p172/ghis/singularity-tsd/Singularity/',
      mountList: ['/tsd', '/net', '/cluster'],
      optional_args: ['-c']
    }
  elsif engine == 'singularity_colossus'
    container = {
      container:  true,
      type:      :singularity,
      repo:      '/cluster/projects/p172/ghis/singularity-tsd/Singularity/',
      mountList: ['/tsd', '/net', '/cluster', '/work'],
      optional_args: ['-c', '--home /cluster/projects/p172/ghis/singularity-tsd/:/srv']
    }
  elsif engine == 'docker'
    container = {
      container:  true,
      type:      :docker,
      repo:      '',
      mountList: [''],
      optional_args: []
    }
  else
    container = {
      container:  false
    }
  end
  return container
end




# ---- references ----
# Structure containing information on references
def set_reference(genome_version: '', intervals_file: false)
  ref = {}
  if genome_version == 'b37'
    ref = {
      name:        'b37',
      genome_fas:  '/XXXXXX/human.fas',
      genome_2bit: '/XXXXXX/human.fas.2bit',
      kgenomes:    '/XXXXXX/1000g.vcf',
      mills:       '/XXXXXX/Mills_and_1000G_gold_standard.indels.b37.vcf',
      dbsnp:       '/XXXXXX/dbsnp.vcf',
      omni:        '/XXXXXX/omni.vcf',
      hapmap:      '/XXXXXX/hapmap.vcf',
      intervals:   intervals_file
    }
  elsif genome_version == 'b38'
    ref = {
      name:        'b38',
      genome_fas:  '/XXXXXX/human.fas',
      genome_2bit: '/XXXXXX/human.fas.2bit',
      kgenomes:    '/XXXXXX/1000g.vcf',
      mills:       '/XXXXXX/Mills_and_1000G_gold_standard.indels.b37.vcf',
      dbsnp:       '/XXXXXX/dbsnp.vcf',
      omni:        '/XXXXXX/omni.vcf',
      hapmap:      '/XXXXXX/hapmap.vcf',
      intervals:   intervals_file
    }
  elsif genome_version == 'test'
    ref = {
      name:        'hg38_local_test',
      genome_fas:  '/media/01-workspace/01-data/refdata/hg38/Homo_sapiens_assembly38.fasta',
      genome_2bit: '/media/01-workspace/01-data/refdata/hg38/Homo_sapiens_assembly38.fasta.2bit',
      kgenomes:    '/media/01-workspace/01-data/refdata/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
      mills:       '/media/01-workspace/01-data/refdata/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
      dbsnp:       '/media/01-workspace/01-data/refdata/hg38/dbsnp_146.hg38.vcf.gz',
      omni:        '/media/01-workspace/01-data/refdata/hg38/1000G_omni2.5.hg38.vcf.gz',
      hapmap:      '/media/01-workspace/01-data/refdata/hg38/hapmap_3.3.hg38.vcf.gz',
      intervals:   intervals_file
    }
  elsif genome_version == 'tsd_test'
    ref = {
      name:        'hg38_local_test',
      genome_fas:  '/media/01-workspace/01-data/refdata/hg38/Homo_sapiens_assembly38.fasta',
      genome_2bit: '/media/01-workspace/01-data/refdata/hg38/Homo_sapiens_assembly38.fasta.2bit',
      kgenomes:    '/media/01-workspace/01-data/refdata/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
      mills:       '/media/01-workspace/01-data/refdata/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
      dbsnp:       '/media/01-workspace/01-data/refdata/hg38/dbsnp_146.hg38.vcf.gz',
      omni:        '/media/01-workspace/01-data/refdata/hg38/1000G_omni2.5.hg38.vcf.gz',
      hapmap:      '/media/01-workspace/01-data/refdata/hg38/hapmap_3.3.hg38.vcf.gz',
      intervals:   intervals_file
    }
  else
    puts '--genome_version must be b37 or b38'
  end
  return ref
end


# ---- Resources to use ----
def set_resources(type: '')
  config = {}
  if type == :laptop
    config = {
      spark: true,
      cpu: '8',
      ram: '24'
    }
  elsif type == :small
    config = {
      spark: true,
      cpu: '12',
      ram: '24'
    }
  elsif type == :medium
    config = {
      spark: true,
      cpu: '16',
      ram: '32'
    }
  elsif type == :large
    config = {
      spark: true,
      cpu: '32',
      ram: '64'
    }
  else
    config = {
      spark: false,
      cpu: '1',
      ram: '16'
    }
  end
  return config
end
