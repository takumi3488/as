#!/usr/bin/env ruby
#
# this program is a part of KANSAI16F suite
#   written by Hiroki Funashima, 2016 Kobe University
#
# check class Self-build or not.
#   if there is not Fortran compiler, self-build is not done.
#
class MakeCheck
  require 'fileutils'
  include FileUtils
  require 'PrintIO'
  include PrintIO

  def initialize(configure)
    @configure = configure
    @kansaidir = @configure[:kansai]
    @pwd = ENV['PWD']

    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    #
    # logfile
    #
    logfile = @prefix + @index + '.log'
    log_add(logfile)
  end

  #
  # check fortran compiler
  #
  def check_fc(compiler)
    if `which #{compiler}`.empty?
      false
    else
      true
    end
  end

  def makeKANSAI
    executable = @configure[:executable]
    compiler = @configure[:compiler]
    if check_fc(compiler)
      if File.exist?("#{@pwd}/prmtsp.f")
        @fout.print "i found prmtsp.f\n"
        if File.exist?('Makefile')
          @fout.print "i found Makefile here.\n"
          `make tapwmn`
        else
          remote_makefile = @kansaidir + '/Makefile'
          if File.exist?(remote_makefile)
            @fout.print "i will build executable using #{@kansaidir}/Makefile\n"
            `make tapwmn -f #{@kansaidir}/Makefile`
          else
            copy_executable(executable)
          end
        end
      else
        copy_executable unless File.exist?(executable)
      end
    else
      copy_executable(executable)
    end
  end

  def copy_executable(executable)
    original = @kansaidir + "/#{executable}"
    if File.exist?(original)
      @fout.print "copy executable from #{@kansai}\n"
      FileUtils.cp(original, '.')
    else
      @fout.print "executable tapwmn is not found in #{@kansai}\n"
      log_finalize
      exit
    end
  end
end
