################################################################################

################################################################################
require 'open3'
require 'pathname'
require 'fileutils'
require 'time'

# Looging with level and nice colors
class Logger
  def self.info(task: 'Workflow', message: '', logfile: false)
    tag = green('Info')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  def self.warning(task: 'Workflow', message: '', logfile: false)
    tag = yellow('Warning')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  def self.execute(task: 'Workflow', message: '', logfile: false)
    tag = blue('Execute')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  def self.finished(task: 'Workflow', message: '', logfile: false)
    tag = cyan('Finished')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  def self.skip(task: 'Workflow', message: '', logfile: false)
    tag = cyan('Skipped')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  def self.error(task: 'Workflow', message: '', logfile: false)
    tag = red('Error')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  def self.debug(task: 'Workflow', message: '', logfile: false)
    tag = magenta('Debug')
    printlog(task: task, message: message, logfile: logfile, log_tag: tag)
  end

  # private

  # Add color methods to the String Class
  def self.black(m);          "\033[30m#{m}\033[0m" end
  def self.red(m);            "\033[31m#{m}\033[0m" end
  def self.green(m);          "\033[32m#{m}\033[0m" end
  def self.yellow(m);         "\033[33m#{m}\033[0m" end
  def self.blue(m);           "\033[34m#{m}\033[0m" end
  def self.magenta(m);        "\033[35m#{m}\033[0m" end
  def self.cyan(m);           "\033[36m#{m}\033[0m" end
  def self.gray(m);           "\033[37m#{m}\033[0m" end
  def self.bg_black(m);       "\033[40m#{m}\033[0m" end
  def self.bg_red(m);         "\033[41m#{m}\033[0m" end
  def self.bg_green(m);       "\033[42m#{m}\033[0m" end
  def self.bg_yellow(m);      "\033[43m#{m}\033[0m" end
  def self.bg_blue(m);        "\033[44m#{m}\033[0m" end
  def self.bg_magenta(m);     "\033[45m#{m}\033[0m" end
  def self.bg_cyan(m);        "\033[46m#{m}\033[0m" end
  def self.bg_gray(m);        "\033[47m#{m}\033[0m" end
  def self.bold(m);           "\033[1m#{m}\033[22m" end
  def self.reverse_color(m);  "\033[7m#{m}\033[27m" end

  def self.printlog(task:, message:, logfile: false, log_tag: 'No Tag info')
    task       = bold(task.ljust(15, ' '))
    log_tag    = log_tag.ljust(15, ' ')
    logmessage = "#{Time.now.strftime("%d_%m_%Y-%H_%M_%S")}\t#{log_tag}\t#{task}  \t#{message}"
    puts logmessage
    save logmessage, logfile
  end

  def self.save(logmessage, logfile)
    # Write Log File
    if logfile
      dest = File.open(logfile, 'a') # a = append
      dest.puts logmessage
      dest.close
    end
  end
end

# The class to run commands
class Runner
  attr_accessor :cmd_line
  attr_reader   :captured_stdout, :captured_stderr, :pid, :exit_status

  def initialize(task_name: '', wf_log: false, task_log_dir: false)
    @task_name       = task_name
    @wf_log          = wf_log
    @task_log_dir    = task_log_dir
    @cmd_line        = []
    @docker_cmd      = false
    @singularity_cmd = false
    @pid             = false # pid of the started process.
    @exit_status     = false
    @captured_stdout = ''
    @captured_stderr = ''
    # Check log dir for the tasks
    if @task_log_dir
      @task_log_dir = File.expand_path @task_log_dir
      @task_log_dir = @task_log_dir + "/"
      checkTaskLogDir @task_log_dir
    end
  end

  def checkTaskLogDir(path)
    if File.exist? @task_log_dir
      if not File.directory? @task_log_dir
        Logger.error task: @task_name, message: "Cannot use #{@task_log_dir} as Log directory for the current task", logfile: @wf_log
        exit
      end
    else
      Logger.info task: @task_name, message: "create log directory for this task : #{@task_log_dir}", logfile: @wf_log
      Dir.mkdir @task_log_dir
    end
  end

  # generic function to run the command in a container
  def in_container( engine: false, imgURL: false, additional_options: [''],  mountList: [''])
    if engine == :docker
      in_docker( imgURL: imgURL, additional_options: additional_options, mountList: mountList)
    elsif engine == :singularity
      in_singularity(imgURL: imgURL, additional_options: additional_options, mountList: mountList)
    end
  end

  # wrap the command inside Docker
  def in_docker(imgURL: false, additional_options: [''],  mountList: [''])
    if @singularity_cmd
      Logger.warning task: @task_name, message: 'singularity found, you cannot define a docker', logfile: @wf_log
      return
    end
    mount_list_param = self.mount_opts(engine: :docker, mountList: mountList)
    options_params = self.additional_opts(engine: :docker, additional_options: additional_options)
    base_opts = '--rm -ti'
    opts_list = base_opts + ' ' + options_params
    @docker_cmd = "docker run -u=#{guid}:#{ggid} #{opts_list} #{mount_list_param} --name=#{@task_name} #{imgURL} sh -c "
  end

  # wrap the command inside singularity
  def in_singularity(imgURL: false, additional_options: [''], mountList: [])
    if @docker_cmd
      Logger.warning task: @task_name, message: 'docker found, you cannot define a singularity', logfile: @wf_log
      return
    end
    mount_list_param = self.mount_opts(engine: :singularity, mountList: mountList)
    options_params = self.additional_opts(engine: :docker, additional_options: additional_options)
    base_opts = ''
    opts_list = base_opts + ' ' + options_params
    @singularity_cmd = "singularity exec #{mount_list_param} #{opts_list} #{imgURL} sh -c "
  end

 # build mount list argument
 def mount_opts(engine: nil, mountList: nil)
  mount_list_param = ''
  mountList.each do |mount_point|
    next if mount_point == ''
    #if not File.exist? mount_point and not Dir.exist? mount_point
    #  mount_point = File.expand_path mount_point
    #else
    #  Logger.error task: @task_name, message: "cannot bind #{mount_point} to #{engine}", logfile: @wf_log
    #  exit(1)
    #end
    if engine == :singularity
      mount_list_param = mount_list_param + ' -B ' + mount_point
    elsif engine == :docker
      mount_list_param = mount_list_param + ' -v=' + mount_point + ':' + mount_point
    end
  end
  return mount_list_param
 end

 # build other optional args
 def additional_opts(engine: nil, additional_options: nil)
  optional_param = ''
  additional_options.each do |option|
    if engine == :singularity
      next if option.include? '-u'
      next if option.include? '-v'
      next if option.include? '--name'
      optional_param = optional_param + ' ' + option
    elsif engine == :docker
      optional_param = optional_param + ' ' + option
    end
  end
  return optional_param
 end

  # Execute the command
  def run(debug: false, skip: false)
    self.build!
    self.compress_spaces!
    start_time = Time.new
    end_time   = 'none'
    Logger.execute task: @task_name, message: @cmd_line, logfile: @wf_log
    if (not debug) and (not skip)
      stdout_log  = ''
      stderr_log  = ''
      exit_status = ''
      # Execute Here
      stdout_log, stderr_log, exit_status = Open3.capture3(@cmd_line)
      end_time = Time.new
      # Write Log from the app
      stdout_file = "#{@task_log_dir}Log_#{@task_name}_#{start_time.strftime("%d_%m_%Y-%H_%M_%S")}_stdout.log"
      Logger.info task: @task_name, message: "stdout log file : #{stdout_file}", logfile: @wf_log
      stdout_file = File.open(stdout_file, 'a') # a = append
      stdout_file.puts stdout_log
      stdout_file.close
      stderr_file = "#{@task_log_dir}Log_#{@task_name}_#{start_time.strftime("%d_%m_%Y-%H_%M_%S")}_stderr.log"
      Logger.info task: @task_name, message: "stderr log file : #{stderr_file}", logfile: @wf_log
      stderr_file = File.open(stderr_file, 'a') # a = append
      stderr_file.puts stderr_log
      stderr_file.close
      # # Display infos
      #Logger.finished task: @task_name, message: "stdout : \n  #{stdout_log}",  logfile: @wf_log
      Logger.finished task: @task_name, message: "exit level : #{exit_status}", logfile: @wf_log
      if exit_status != 0
        Logger.warning task: @task_name, message: "exit level is not nul", logfile: @wf_log
      end
      if stderr_log != ''
        Logger.info  task: @task_name, message: "stderr is not empty, you should check the stderr log file",  logfile: @wf_log
        # Logger.warning  task: @task_name, message: "stderr : \n  #{stderr_log}",  logfile: @wf_log
      end
      Logger.finished task: @task_name, message: "Start\t#{start_time.to_i}\tEnd\t#{end_time.to_i}\tExecution time\t#{end_time - start_time}", logfile: @wf_log
    elsif debug
      Logger.skip task: @task_name, message: 'Debug mode', logfile: @wf_log
    else
      Logger.skip task: @task_name, message: 'Results already exists', logfile: @wf_log
    end
  end

  # command line Array to String
  def build!
    @cmd_line = @cmd_line.join ' ' if @cmd_line.class == Array
    if @docker_cmd
      @cmd_line = @docker_cmd + "' " + @cmd_line + " '"
    end
    if @singularity_cmd
      @cmd_line = @singularity_cmd + "' " + @cmd_line + " '"
    end
  end

  # remove multiple spaces
  def compress_spaces!
    @cmd_line.gsub!(/\s+/, ' ')
  end

  def guid
    return `id -u`.chomp!
  end

  def ggid
    return `id -g`.chomp!
  end
end


# Return an array of file name (without path) from a path and a shell regular expression
def get_file_names(path: '', filter: '', full_path: false, wf_log: nil)
  files_list = Dir.glob(path + '/' + filter)
  name_list  = []
  files_list.each do |file|
    file = File.basename(file) if not full_path
    name_list.push file
  end
  if name_list.empty?
    Logger.error message: "No file found #{path}/#{filter}", logfile: wf_log
    exit
  end
  return name_list
end


# Return an array of Array of 2 file names (without path) from a path and a shell regular expression
# Files are paired by Tag
def get_file_pair_names(path: '', filter: '', tag1: '', tag2: '', full_path: false, wf_log: nil)
  files_list = Dir.glob(path + '/' + filter)
  files_F    = []
  files_R    = []
  file_pairs = []
  files_list.each do |file|
    file = File.basename(file) if not full_path
    files_F.push file if file.include? tag1
    files_R.push file if file.include? tag2
  end
  files_F.sort!
  files_R.sort!
  if files_F.size != files_R.size
    Logger.warning message: "found #{files_F.size} files with tag #{tag1} and #{files_R.size} files with tag #{tag2}", logfile: wf_log
  end
  for i in 0..(files_F.size - 1) do
    file_pairs.push [files_F[i],files_R[i]]
  end
  if file_pairs.empty?
    Logger.error message: "No file found #{path}/#{filter} with tag #{tag1} and #{tag2}", logfile: wf_log
    exit
  end
  return file_pairs
end

# Remove the last extention ex :  aa.bb.cc -> aa.bb
def remove_extention(file)
  file_tmp = file.split '.'
  file_tmp.pop
  return file_tmp.join  '.'
end
