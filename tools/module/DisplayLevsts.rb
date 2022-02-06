#!/usr/bin/env ruby
#
# analysis tool for results for levsts
#   ver. 2.0
#     written by HF in Kobe 14, April 2017
#
#

require 'ParseLevsts'
require 'ParseOrbital'
require 'CheckNwindows'
require 'IOUtils'
include IOUtils

class DisplayLevsts
  def initialize(filename)
    @levsts_obj = ParseLevsts.new(filename)
    @orbital_obj = ParseOrbital.new
    @nwin = CheckNwindows.new(filename).nwin
    @eps = 0.2
  end

  def add_outfile(filename)
    @levsts_obj.add_file(filename)
  end

  def display
    obj = @levsts_obj
    data = obj.levstsdata
    1.upto(@levsts_obj.nka) do |ika|
      1.upto(@levsts_obj.nspin) do |ispin|
        print "=== IA = #{ika}"
        print "  ISPIN = #{ispin} :"
        print " #{@orbital_obj.atom_name[ika - 1]}\n"
        print format('ITER')
        max_l = @orbital_obj.find_max_l(ika)
        0.upto(max_l) do |l|
          print format('                 phi0(L=%1s)', l)
          print format('                          phi1(L=%1s)', l)
        end
        print "\n"
        1.upto(@levsts_obj.niter) do |iter|
          if @nwin.nil?
            print format(' %3s', iter)
          elsif @nwin == 1
            print format(' %3s', iter)
          else
            print format(' %3s', (iter + 1) / 2)
          end
          display_energy_info(max_l, data, iter, ispin, ika)
        end
        print "----\n"
        display_deviation(ispin, ika)
        print "\n"
      end
    end
  end

  def display_deviation(ispin, ika)
    max_l = @orbital_obj.find_max_l(ika)
    if @nwin.nil? || (!@nwin.nil? && @nwin == 1)
      nwin = 0
    else
      nwin = 1
    end
    0.upto(nwin) do |iwin|
      print 'Dev:'
      0.upto(max_l) do |l|
        0.upto(1) do |type|
          if @nwin.nil? || (!@nwin.nil? && @nwin == 1)
            dev_win1(ispin, ika, l, type)
          else
            dev_win2(ispin, ika, l, type, iwin)
          end
        end
      end
      print "\n"
    end
  end

  def dev_win1(ispin, ika, l, type)
    nwin = 0
    deviation = @levsts_obj.standard_deviation(ispin, ika, l, type, nwin)
    nspace =
      if type == 0
        16
      else
        25
      end
    nspace.times { print ' ' }
    if deviation.nil?
      print '        NA'
    else
      print format('%10.4e', deviation)
    end
  end

  def dev_win2(ispin, ika, l, type, iwin)
    deviation = @levsts_obj.standard_deviation(ispin, ika, l, type, iwin+1)
    nspace =
      if type == 0
        16
      else
        25
      end
    nspace.times { print ' ' }
    if deviation.nil?
      print '        NA'
    else
      print format('%10.4e', deviation)
    end
  end

  def display_energy_info(max_l, data, iter, ispin, ika)
    0.upto(max_l) do |l|
      0.upto(1) do |i|
        display_erange(data, iter, ispin, ika, l, i)
        display_grad(data, iter, ispin, ika, l, i)
        if i == 0
          print '    '
        else
          display_parity(data, iter, ispin, ika, l, i)
        end
      end
    end
    print "\n"
  end

  def display_erange(data, iter, ispin, ika, l, i)
    ekind = [:e0, :e1]
    zero_point = data[iter - 1][ispin - 1][ika - 1][l][ekind[i]]

    if zero_point # have solution or not.
      erange = @levsts_obj.calc_erange(iter, ispin, ika, l, i)
      print format('  %7.4f', zero_point)
      print ':('
      erange_print(erange[0], @eps)
      print '->'
      erange_print(erange[1], @eps)
      print ')'
    else
      print format('      NA:(             )')
    end
  end

  def display_grad(data, iter, ispin, ika, l, i)
    gkind = [:grad0, :grad1]
    if grad = data[iter - 1][ispin - 1][ika - 1][l][gkind[i]].nil?
      print ':x '
    else
      if data[iter - 1][ispin - 1][ika - 1][l][gkind[i]] > 0
        print ':+ '
      else
        print ':- '
      end
    end
  end

  def erange_print(energy, eps)
    if energy < eps
      print_green format('%5.2f', energy)
    else
      print format('%5.2f', energy)
    end
  end

  def display_parity(data, iter, ispin, ika, l, i)
    parity = data[iter - 1][ispin - 1][ika - 1][l][:parity]
    if parity
      if @orbital_obj.valence?(ika, l).nil?
        print '   ' if i == 1
      else
        print_parity_sign(iter, ika, l, parity)
      end
    else
      print ':NA'
    end
  end

  def define_icore(iter, ika, l, parity)
    return 0 if @nwin.nil?
    icore = 0
    if @nwin == 2 && @orbital_obj.valence?(ika, l).size > 1
      icore =
        if iter.odd?
          0
        else
          1
        end
    end
    icore
  end

  def print_parity_sign(iter, ika, l, parity)
    icore = define_icore(iter, ika, l, parity)
    if @orbital_obj.valence?(ika, l)[icore]
      sign = define_sign(parity)
      if parity.even? && @orbital_obj.valence?(ika, l)[icore].even?
        print_blue format('(%1s)', sign)
      elsif parity.odd? && @orbital_obj.valence?(ika, l)[icore].odd?
        print_blue format('(%1s)', sign)
      else
        print_red format('(%1s)', sign)
      end
    else
      print_red format('(%1s)', sign)
    end
  end

  def define_sign(parity)
    if parity.even?
      'e'
    else
      'o'
    end
  end
end
