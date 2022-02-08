#!/usr/bin/env ruby
#
# this program is a part of KANSAI16F suite
#   written by Hiroki Funashima, 2016 Kobe University
#
# execution of LAPW in KANSAI package class
#
class ExecuteLAPW
  require 'fileutils'
  def initialize(configure)
    
    #
    # initialize variables
    #
    @executable = configure[:executable]
    @index = configure[:index].to_s
    @filelist = configure[:linkfiles]
    @dir = configure[:dir]
    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    @previous_files = configure[:previous_files]
    @tmpfiles = configure[:tmpfiles]
    #
    # logfile
    #
    logfile = @prefix + @index + '.log'
    @fout = open(logfile, 'a')
    go
  end

  def go
    check_old_data
    execution
    filecopy
    delete_tmpfiles
    finalize
  end

  def finalize
    @fout.close
  end

  def execution
    infile = @dir + '/' + @prefix + @index + @suffix
    unless File.exist?(infile)
      @fout.print "================ Error ===============\n"
      @fout.print "file:#{infile} is not found.\n"
      @fout.print "======================================\n"
      finalize
      exit
    end
    outfile = @dir + '/' + @prefix + @index + '.out'
    errfile = @prefix + @index + '.err'
    @fout.print "execute scf calculation\n"
    @fout.print " file check for #{@executable} ..."
    if File.exist?("./#{@executable}")
      @fout.print "ok\n"
      FileUtils.rm(errfile) if File.exist?(errfile)
      t1 = Time.now
      `./#{@executable} < #{infile} > #{outfile} 2> #{errfile}`
      t2 = Time.now
      FileUtils.rm(errfile) unless File.size?(errfile)
      @fout.print "#{@executable} has done(#{t2- t1} sec).\n\n"
    else
      @fout.print "failed\n"
      finalize
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
      finalize
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
    @fout.print "remove temporary files...\n"
    @tmpfiles.each do |fn|
      filename = 'fort.' + fn.to_s
      if File.exist?(filename)
        FileUtils.rm(filename)
        @fout.print " ->file:#{filename} was deleted.\n"
      end
    end
  end
end
