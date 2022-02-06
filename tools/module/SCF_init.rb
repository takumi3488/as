#!/usr/bin/env ruby
#
#

#
# initialize data 
#
class SCF_init
  require 'fileutils'
  def initialize(configure)
    @data_prefix = configure[:data_prefix]
    @filelist = configure[:filelist]
    @initdir = configure[:initdir]
    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    @dir = configure[:dir]
    @results_list = configure[:results_list]
    main
  end

  def main
    print "initialize data files.\n"
    clean_fort_files
    clean_result_files
    clean_out_files
    clean_log_files
    clean_err_files
    clean_old_datafiles
    copy_previous_series
  end

  def clean_fort_files
    list = Dir.glob("fort.*")
    if list.size > 0
      print "fort files ...\n"
      list.each do |filename|
        print " -> #{filename} was deleted.\n"
        FileUtils.rm(filename)
      end
    end
  end

  def clean_log_files
    list = Dir.glob("#{@prefix}*.log")
    if list.size > 0
      print "clean log files ...\n"
      list.each do |filename|
        print " -> #{filename} was deleted.\n"
        FileUtils.rm(filename)
      end
    end
  end

  def clean_err_files
    list = Dir.glob("#{@prefix}*.err")
    if list.size > 0
      print "clean log files ...\n"
      list.each do |filename|
        print " -> #{filename} was deleted.\n"
        FileUtils.rm(filename)
      end
    end
  end

  def clean_result_files
    @results_list.each do |result_name|
      list = Dir.glob("#{@dir}/#{@prefix}*#{result_name}#{@suffix}")
      if list.size > 0
        print "clean #{result_name} files ...\n"
        list.each do |filename|
          print " -> #{filename} was deleted.\n"
          FileUtils.rm(filename)
        end
      end
    end
  end

  def clean_out_files
    list = Dir.glob("#{@dir}/#{@prefix}*.out")
    if list.size > 0
      print "clean out files ...\n"
      list.each do |filename|
        print " -> #{filename} was deleted.\n"
        FileUtils.rm(filename)
      end
    end
  end

  def clean_old_datafiles
    delete_file_list = Dir.glob("#{@data_prefix}.*")
    if delete_file_list.size > 0
      delete_file_list.each do |filename|
        print " -> #{filename} was deleted.\n"
        FileUtils.rm(filename)
      end
    end
  end

  def copy_previous_series
    @filelist.each do |fn|
      filename = @data_prefix + fn.to_s
      original = @initdir + '/' + filename
      if File.exist?(original)
        print " -> copy #{filename} from #{@initdir}.\n"
        FileUtils.cp(original, '.')
      end
    end
  end
end
