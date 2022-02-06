#!/usr/bin/env ruby
#
# this program is a part of KANSAI16F suite
#   written by Hiroki Funashima, 2016 Kobe University
#

#
# display Envrionmental Variables
#   written by Hiroki Funashima, 2016 Kobe University
#
class DisplayEnvs
  require 'fileutils'
  include FileUtils
  require 'PrintIO'
  include PrintIO
  def initialize(configure)
    @configure = configure
    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    @index = configure[:index].to_s
    #
    # logfile
    #
    logfile = @prefix + @index + '.log'
    @success = true
    log_initialize(logfile)
    grand_start
    print_date
    display_envs
    log_finalize
  end

  def display_envs
    @fout.print "-> KANSAI Package ... #{@configure[:kansai]}\n"
    @fout.print "-> Executable ... #{@configure[:executable]}\n"
    @fout.print "-> overwrite executable ... "
    if @configure[:exec_overwrite]
      @fout.print "yes\n"
    else
      @fout.print "no\n"
    end
    @fout.print "-> Dir is #{@configure[:dir]}\n"
    @fout.print '-> check previous output file ...'
    if @outcheck
      @fout.print "yes\n"
    else
      @fout.print "no\n"
    end
    @fout.print "-> Series prefix ... #{@configure[:prefix]}\n"
    @fout.print "-> index of execution ... #{@index}\n"
    @fout.print "-> Series suffix ... #{@configure[:suffix]}\n"
  end

  def print_date
    @fout.print "=== Package Environmental Variables  ===\n"
    @fout.print "\n"
  end
end
