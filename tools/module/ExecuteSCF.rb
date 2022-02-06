#!/usr/bin/env ruby
#
# this program is a part of KANSAI16F suite
#   written by Hiroki Funashima, 2016 Kobe University
#

#
# execution of LAPW in KANSAI package class
#
class ExecuteSCF
  require 'fileutils'
  include FileUtils
  require 'PrintIO'
  include PrintIO
  def initialize(configure)
    #
    # initialize variables
    #
    @kansai = configure[:kansai]
    @executable = configure[:executable]
    @index = configure[:index].to_s
    @filelist = configure[:linkfiles]
    @dir = configure[:dir]
    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    @previous_files = configure[:previous_files]
    @tmpfiles = configure[:tmpfiles]
    @outcheck = configure[:outcheck]
    @exec_overwrite = configure[:exec_overwrite]
    #
    # logfile
    #
    logfile = @prefix + @index + '.log'
    @success = true
    log_add(logfile)
    go
  end

  def go
    write_common_log(:start)
    print_date
    check_old_data
    execution
    filecopy if @success
    delete_tmpfiles
    log_finalize
    write_common_log(:finish)
  end

  def write_common_log(type)
    logfile = ENV['HOME'] + '/calc.log'
    target = @prefix + @index.to_s
    open(logfile, 'a') do |fout|
      case(type)
      when :start
        fout.print '++ '
        fout.print " #{@executable} has been started: "
      when :finish
        fout.print '-- '
        fout.print " #{@executable} has been finished:"
      else
        fout.print '-- '
        fout.print " #{@executable} has been stopped: "
      end
      fout.print " #{Time.now}  Hostname:#{`hostname`.chomp}  Workspace:#{ENV['PWD']}  Target:#{target}\n"
    end
  end

  def print_date
    @fout.print "=== Execution SCF calculation package  ===\n"
    @fout.print " Date:#{Time.now}\n"
    @fout.print " Hostname:#{`hostname`}"
    @fout.print "\n"
  end

  def check_outfile(infile, filename)
    if File.exist?(filename)
      @fout.print "========== outfile is exist ===========\n"
      @fout.print "your inputfile is #{infile}\n"
      @fout.print "output file:#{filename} is exist.\n" 
      @fout.print "Maybe, your calculation has already done.\n"
      log_finalize
      write_common_log('')
      exit
    end
  end

  def execution
    infile = @dir + '/' + @prefix + @index + @suffix
    unless File.exist?(infile)
      @fout.print "================ Error ===============\n"
      @fout.print "file:#{infile} is not found.\n"
      @fout.print "======================================\n"
      log_finalize
      write_common_log('')
      exit
    end
    outfile = @dir + '/' + @prefix + @index + '.out'
    check_outfile(infile, outfile) if @outcheck
    errfile = @prefix + @index + '.err'
    @fout.print "execute scf calculation\n"
    if @exec_overwrite
      @fout.print " -> copy executable #{@executable} from #{@kansai}\n"
      origexec = @kansai + '/' + @executable
      if File.exist?(origexec)
        FileUtils.cp(origexec, '.')
        `chmod u+x #{@executable}`
      else
        @fout.print "===== Error =====\n"
        @fout.print "executable:#{origexec} is not found.\n"
        @fout.print "please `make' in #{@kansai}\n"
        log_finalize
        write_common_log('')
        exit
      end
    end
    @fout.print " file check for #{@executable} ..."
    if File.exist?("./#{@executable}")
      @fout.print "ok\n"
      FileUtils.rm(errfile) if File.exist?(errfile)
      t1 = Time.now
      `./#{@executable} < #{infile} > #{outfile} 2> #{errfile}`
      t2 = Time.now
      if File.size?(errfile)
        @fout.print "******** Warning ********\n"
        @fout.print "Error file was generated.\n"
        @fout.print "Check: #{errfile}\n"
        @fout.print "*************************\n"
        @success = false
      else
        FileUtils.rm(errfile)
      end
      @fout.print "#{@executable} has done (time:#{t2 - t1} sec).\n\n"
    else
      @fout.print "failed\n"
      log_finalize
      write_common_log('')
      exit
    end
  end

  def check_old_data
    @fout.print "check old data...\n"
    check = true
    @previous_files.each do |dev_num, contents|
      original_file = 'fort.' + dev_num.to_s
      @fout.print "-> #{original_file}  :(#{contents})..."
      if File.exist?(original_file)
        @fout.print "ok\n"
      else
        @fout.print "not found.\n"
        check = false
      end
    end
    unless check
      @fout.print "===== Error: previous data is not found.  =====\n"
      log_finalize
      write_common_log('')
      exit
    end
    @fout.print "\n"
  end

  def filecopy
    @fout.print "copy calculation results...\n"
    @filelist.each do |dev_num, filename|
      original_file = 'fort.' + dev_num.to_s
      newfile = @dir + '/' + @prefix + @index + filename + @suffix
      if File.exist?(original_file)
        FileUtils.cp(original_file, newfile)
        @fout.print " ->#{newfile} was generated from #{original_file}\n"
      end
    end
    @fout.print "\n"
  end

  def delete_tmpfiles
    @fout.print "remove symbolic link to temporary files...\n"
    @tmpfiles.each do |fn|
      filename = 'fort.' + fn.to_s
      FileUtils.rm(filename)
      @fout.print " ->file:#{filename} was deleted.\n"
    end
  end
end
