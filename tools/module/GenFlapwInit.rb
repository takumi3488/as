#!/usr/bin/env ruby
#
# Data convert lapw input file to flapw initial file
#   written by Hiroki Funashima, 2 May 2017 in Kobe
#

#
# class for fw1.dat from sl?.dat
#
require 'fileutils'
class GenFlapwInit
  def initialize(filename)
    @filename = filename
    @prefix = filename
    @original = filename
    @newfile = 'fw1.dat'
    @header = ''
    @body = ''
    @footer = ''
    @info_level = 1
    @cut_off_flapw_sw = 6
    @efg_param = 'EFG1'
    @extlist = %w(dst.dat eld.dat esd.dat chg.dat fchg.dat pot.dat fpot.dat )
    parse_prefix
    @filename = @prefix + '.dat'
    set_header
    set_cutoff
    write_fw1
    copy_fort_files
  end

  def parse_prefix
    @extlist.each do |ext|
      @prefix = @prefix.split(ext)[0]
    end
    ext = File.extname(@filename)
    @prefix = @prefix.split(ext)[0] unless ext.empty?
  end

  def write_fw1
    print "generate fw1.dat from #{@filename}\n"
    fout = open(@newfile, 'w')
    fout.print @header
    fout.print "FULL POTENTIAL\n"
    fout.print " #{@cut_off_flapw_pw} #{@cut_off_flapw_sw}    # cutoff for spw and ssw (Full Potential)\n"
    fout.print " #{@info_level}          # infomation level\n"
    fout.print "#{@efg_param}        # Electric Field Gradient. EFG1:on EFG0:off\n"
    fout.print "LAPW\n"
    fout.print " #{@cut_off_lapw_pw} #{@cut_off_lapw_sw}     # cutoff for spw and ssw (MTA)\n"
    fout.print @body
    fout.print @footer
    fout.close
  end

  def set_header
    open(@filename, 'r') do |fin|
      fin.each_line do |linedata|
        @header << linedata if fin.lineno <= 2
      end
    end
  end

  def set_cutoff
    iter_flag = false
    iter_count = 0
    open(@filename, 'r') do |fin|
      fin.each_line do |linedata|
        if fin.lineno == 4
          @cut_off_lapw_pw = linedata.split[0].to_f
          @cut_off_flapw_pw = @cut_off_lapw_pw * 2.0
          @cut_off_lapw_sw = linedata.split[1].to_i
        elsif fin.lineno > 4
          if iter_flag
            iter_count += 1
            if iter_count == 1
              @footer << " 0      # mixing method(0:simple iteration. 1:simple mixing. 2:broyden method(you need fort.36 and 37).\n"
            else
              jpr, pmix, amix = linedata.chomp.split
              if iter_count == 2
                if pmix.to_f > 0
                  @footer << " #{jpr}  #{pmix}  #{amix}\n"
                else
                  @footer << " #{jpr}  #{pmix.split('-')[1]}  #{amix}\n"
                end
              elsif iter_count == 3
                if pmix.to_f > 0
                  @footer << " #{jpr} -#{pmix}  #{amix}\n"
                else
                  @footer << " #{jpr}  #{pmix}  #{amix}\n"
                end
              else
                if pmix.to_f > 0
                  @footer << " #{jpr}  #{pmix}  #{amix}\n"
                else
                  @footer << " #{jpr} #{pmix}  #{amix}\n"
                end
              end
            end
          else iter_flag
            iter_flag = true if linedata.chomp.strip =~ /^ITER/i
            @body << linedata
          end
        end
      end
    end
  end

  def copy_fort_files
    path = 'wk'
    FileUtils.mkdir_p(path) unless FileTest.exist?(path)
    lapwdir = File.dirname(@original)
    [1, 30, 35].each do |i|
      next if lapwdir.nil?
      original = lapwdir + '/wk/fort.' + i.to_s
      newfile = 'wk/fort.' + i.to_s 
      if File.exist?(newfile)
        next if FileUtils.cmp(original, newfile)
      end
      if File.exist?(original)
        print " file copy from #{original} to #{newfile}\n"
        FileUtils.copy(original, newfile)
      end
    end
  end
end
