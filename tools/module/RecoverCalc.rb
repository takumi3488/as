#!/usr/bin/env ruby
#
# Recover previous calculation Results
#   written by Hiroki Funashima, Kobe University, 2016
#
class RecoverCalc
  require 'fileutils'
  include FileUtils
  require 'MyTmpDir'
  include MyTmpDir
  require 'CheckTmpFiles'

  def initialize(configure)
    @tmpdir = my_tmp_dir
    @prefix = configure[:prefix]
    @index  = configure[:index].to_s
    @suffix = configure[:suffix]
    @linkfiles = configure[:linkfiles]
    @tmpfiles = configure[:tmpfiles]
    @dir = configure[:dir]
    delete_links
    delete_tmpfiles
    recover
  end

  def recover
    print "Recover calc: #{@prefix}#{@index} series...\n"
    @linkfiles.each do |dev_num, filename|
      original_file = @dir + '/' + @prefix + @index + filename + @suffix
      newfile = 'fort.' + dev_num.to_s
      if File.exist?(original_file)
        FileUtils.cp(original_file, newfile)
        print " -> #{newfile} was generated from #{original_file}\n"
      end
    end
  end

  def delete_links
    print "remove symbolic link to temporary files...\n"
    @tmpfiles.each do |fn|
      filename = 'fort.' + fn.to_s
      if File.exist?(filename)
        FileUtils.rm(filename)
        print " ->file:#{filename} was deleted.\n"
      end
    end
  end

  def delete_tmpfiles
    print "remove temporary file in #{@tmpdir}\n"
    @tmpfiles.each do |fn|
      filename = "#{@tmpdir}/fort.#{fn}"
      print " -> remove temporary file: fort.#{fn}\n"
      FileUtils.rm(filename) if File.exist?(filename)
    end
  end
end
