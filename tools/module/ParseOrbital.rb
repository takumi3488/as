#!/usr/bin/env ruby
#

#
# clas file for parsing fort.1
#   toparse fort.1 for parse levst
#    to get infomation about valence/core and something.
#      written by HF in Kobe 2017
#
class ParseOrbital
  attr_reader :orbital, :nka
  attr_reader :atom_name
  def initialize
    if File.exist?('fort.1')
      @filename = 'fort.1'
    elsif File.exist?('wk/fort.1')
      @filename = 'wk/fort.1'
    else
      print "file: fort.1 is not found.\n"
      exit
    end
    @orbital = []
    parse
  end

  def parse
    ika = 0
    @atom_name = []
    linecounter = nil
    header = true
    open(@filename, 'r') do |fin|
      fin.each_line do |line|
        next if fin.lineno == 1
        if line.chomp.strip =~ /^[A-Z,a-z]/
          header = false
          linecounter = 0
          @atom_name.push(line.chomp.split[0])
          ika += 1
          @orbital[ika - 1] = [] if @orbital[ika - 1].nil?
        else
          unless header
            linecounter += 1
            if linecounter == 1
              @orbital[ika - 1] = line.chomp.split.map(&:to_i) if linecounter == 1
            elsif linecounter == 2
              if line.chomp.strip =~ /./
                line.chomp.split.map(&:to_i).each { |x| @orbital[ika - 1].push(x) }
              end
            end
          end
        end
      end
    end
    @nka = ika
  end

  def valence?(ika, ll)
    vals = []
    @orbital[ika - 1].each do |orb|
      n = orb / 100
      l = (orb / 10) % 10
      val =
        if orb % 10 > 0
          true
        else
          false
        end
      vals.push(n) if ll == l && val
    end
    vals = nil if vals.empty?
    vals
  end

  def find_max_l(ika)
    max_l = 0
    0.upto(7) { |l| max_l = l if valence?(ika, l) && max_l < l }
    max_l
  end
end
