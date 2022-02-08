#!/usr/bin/env ruby
#
# this program is a part of KANSAI16F suite
#   written by Hiroki Funashima, 2016 Kobe University
#
# check temp directory and link class
#
class CheckTmpFiles
  require 'fileutils'
  include FileUtils
  require 'MyTmpDir'
  include MyTmpDir
  require 'PrintIO'
  include PrintIO
  def initialize(configure)
    @filelist = configure[:tmpfiles]
    @tmpdir = my_tmp_dir

    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    @index = configure[:index]
    #
    # logfile
    #
    logfile = @prefix + @index.to_s + '.log'
    @success = true
    log_add(logfile)
    header
  end

  def header
    @fout.print "=== Generate Link package  ===\n"
  end

  def linking
    generate_my_tmp_dir unless Dir.exist?(@tmpdir)
    @fout.print "temorary directory:#{@tmpdir}\n"
    @filelist.each do |fn|
      filename = "fort.#{fn}"
      `ln -sf #{@tmpdir}/#{filename} .`
      @fout.print " -> temporary file #{filename} was linked.\n"
    end
    log_finalize
  end

  def deleting
    logfile = @prefix + @index + '.log'
    log_add(logfile)
    @fout.print "remove temporary file in #{@tmpdir}\n"
    @filelist.each do |fn|
      filename = "#{@tmpdir}/fort.#{fn}"
      @fout.print " -> remove temporary file: fort.#{fn}\n"
      FileUtils.rm(filename) if File.exist?(filename)
    end
    log_finalize
  end
end
