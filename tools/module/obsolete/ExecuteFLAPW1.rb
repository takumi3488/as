#!/usr/bin/env ruby
#
# this program is a part of KANSAI16F suite
#   written by Hiroki Funashima, 2016 Kobe University
#
# execution of FLAPW in KANSAI package class
#
class ExecuteFLAPW1 < ExecuteLAPW
  require 'fileutils'
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
    @fout.print "execute flapw initial calculation\n"
    @fout.print " file check for #{@executable} ..."
    if File.exist?("./#{@executable}")
      @fout.print "ok\n"
      FileUtils.rm(errfile) unless File.exist?(errfile)
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
end
