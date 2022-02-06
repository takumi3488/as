#!/usr/bin/env ruby

#
# rewrite parse levst
#   written by H.Funashima, 2017 in Kobe
#
require 'KSUtil'
include KSUtil
include Math
class ParseLevsts
  attr_reader :levstsdata
  attr_reader :niter, :nspin, :nka
  def initialize(filename)
    filename = parse_prefix(filename) + '.out' unless filename.nil?
    if File.exist?(filename)
      add_file(filename)
    elsif @niter.nil?
      print "============ Error ===========\n"
      print "file:#{filename} is not found.\n"
      exit
    end
  end

  def add_file(filename)
    filename = parse_prefix(filename) + '.out' unless filename.nil?
    if File.exist?(filename)
      @filename = filename
      parse
    end
  end

  def parse
    linecount = 0
    body = false
    energy = nil
    @niter = 0 if @niter.nil?
    iter = @niter - 1
    isp = nil
    ia = nil
    l = nil
    @levstsdata = [] if @levstsdata.nil?
    open(@filename, 'r') do |fin|
      fin.each_line do |line|
        linedata = line.chomp.strip
        if linedata =~ /^IND0/
          body = true
          linecount = 0
          isp, ia, l = parse_header(linedata)
          iter += 1 if [isp, ia, l] == [1, 1, 0]
          @levstsdata[iter] = [] if @levstsdata[iter].nil?
          @levstsdata[iter][isp - 1] = [] if @levstsdata[iter][isp - 1].nil?
          @levstsdata[iter][isp - 1][ia - 1] =
            [] if @levstsdata[iter][isp - 1][ia - 1].nil?
          @levstsdata[iter][isp - 1][ia - 1][l] =
            {} if @levstsdata[iter][isp - 1][ia - 1][l].nil?
        elsif body
          if linedata.empty?
            body = false
          else
            linecount += 1
            case linecount
            when 1 # energy
              energy = float_array(linedata)
              @levstsdata[iter][isp - 1][ia - 1][l][:energy] = energy
            when 2 # phi0
              @levstsdata[iter][isp - 1][ia - 1][l][:phi0] =
                float_array(linedata)
              @levstsdata[iter][isp - 1][ia - 1][l][:e0] =
                solve_zero(linedata, energy)
              @levstsdata[iter][isp - 1][ia - 1][l][:grad0] =
                monotonic_func?(linedata)
            when 3 # phi1
              @levstsdata[iter][isp - 1][ia - 1][l][:phi1] =
                float_array(linedata)
              @levstsdata[iter][isp - 1][ia - 1][l][:e1] =
                solve_zero(linedata, energy)
              @levstsdata[iter][isp - 1][ia - 1][l][:grad1] =
                monotonic_func?(linedata)
              @levstsdata[iter][isp - 1][ia - 1][l][:parity] =
                even_or_odd?(l, linedata)
            end
          end
        end
      end
    end
    @nspin = isp
    @nka = ia
    @niter = iter + 1
  end

  def float_array(linedata)
    linedata.split.map(&:to_f)
  end

  def parse_header(linedata)
    l = keyword_value('L', linedata)
    ia = keyword_value('IA', linedata)
    isp = keyword_value('ISP', linedata)
    [isp, ia, l]
  end

  def even_or_odd?(l, linedata)
    grad = monotonic_func?(linedata)
    return nil if grad.nil?

    #
    # even: 0
    # odd:  1
    #
    if l.even?
      if grad < 0
        1
      else
        0
      end
    else
      if grad > 0
        1
      else
        0
      end
    end
  end

  def keyword_value(keyword, linedata)
    linedata.upcase.split(keyword + '=')[1].split[0].to_i
  end

  def monotonic_func?(linedata)
    data = linedata.split.map(&:to_f)
    check = true
    0.upto(data.size - 2) do |i|
      if (data[-1] - data[0]) * (data[i + 1] - data[i]) < 0
        check = nil
        break
      end
    end

    if check
      check =
        if (data[-1] - data[0]) > 0
          1
        else
          -1
        end
    end
    check
  end

  def solve_zero(linedata, energy)
    return nil unless monotonic_func?(linedata)
    data = linedata.split.map(&:to_f)
    return nil if data[0] * data[-1] > 0.0
    ans = nil
    0.upto(data.size - 2) do |i|
      next if (data[i + 1] * data[i]) > 0
      x = [energy[i], energy[i + 1]]
      y = [data[i], data[i + 1]]
      ans = linear_solver(x, y)
    end
    ans
  end

  def linear_solver(x, y)
    x[1] - y[1] * (x[1] - x[0]) / (y[1] - y[0])
  end

  def calc_erange(iter, ispin, ika, l, type)
    #
    # type = 0 ... phi0
    #      = 1 ... phi1
    #
    etype =
      if type == 0
        :e0
      else
        :e1
      end
    return nil unless @levstsdata[iter - 1][ispin - 1][ika - 1][l][etype]
    solution = @levstsdata[iter - 1][ispin - 1][ika - 1][l][etype]
    emin = solution - @levstsdata[iter - 1][ispin - 1][ika - 1][l][:energy][0]
    emax = @levstsdata[iter - 1][ispin - 1][ika - 1][l][:energy][-1] - solution
    [emin, emax]
  end

  def standard_deviation(ispin, ika, l, type, nwin)
    #
    # nwin = 0   ... 1 window calc
    #      = 1   ... 2 window calc and semi-core state
    #      = 2   ... 2 window calc and valence state
    #

    #
    # type = 0 ... phi0
    #      = 1 ... phi1
    #
    etype =
      if type == 0
        :e0
      else
        :e1
      end
    sum = 0.0
    sample_size = 0
    solution = 0.0
    #
    # calc average
    #
    1.upto(@niter) do |iter|
      next if skip_iter?(iter, nwin)
      solution = @levstsdata[iter - 1][ispin - 1][ika - 1][l][etype]
      unless solution.nil?
        sum += solution
        sample_size += 1
      end
    end
    return nil if sample_size == 0

    average = sum / sample_size
    variance = 0.0
    #
    # calc variance
    #
    1.upto(@niter) do |iter|
      next if skip_iter?(iter, nwin)
      solution = @levstsdata[iter - 1][ispin - 1][ika - 1][l][etype]
      variance += (solution - average)** 2 unless solution.nil?
    end
    variance /= sample_size
    sqrt(variance)
  end

  def skip_iter?(iter, nwin)
    #
    # nwin = 0   ... 1 window calc
    #      = 1   ... 2 window calc and semi-core state
    #      = 2   ... 2 window calc and valence state
    #
    return false if nwin == 0
    return false if nwin == 1 && iter.odd?
    return false if nwin == 2 && iter.even?
    true
  end
end
