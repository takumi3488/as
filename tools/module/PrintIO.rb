#!/usr/bin/env ruby

#
# submodule PrintIO
#
module PrintIO
  def log_initialize(logfile)
    @fout = open(logfile, 'w')
  end

  def log_add(logfile)
    @fout = open(logfile, 'a')
  end

  def grand_start
    @fout.print "++++++++++++++ KANSAI CALCULATION PACKAGE ++++++++++++++++++\n"
    @fout.print "\n"
  end

  def log_finalize
    @fout.print "==========================================\n"
    @fout.print "\n"
    @fout.close
  end
end
