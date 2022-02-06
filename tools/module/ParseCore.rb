#!/usr/bin/env ruby
#
#
#

#
# class for ChecCore
#
class ParseCore
  def initialize(outfile)
    if File.exist?(outfile)
      @filename = outfile
      @ika = 0
    else
      print "===== Error(CheckCore) =====\n"
      print "file:#{outfile} is not found.\n"
      print "============================\n"
      exit
    end
  end

  def parse_core_state
    header = true
    open(@filename, 'r') do |fin|
      fin.each_line do |line|
        if line =~ /EIGENENERGY\sOF\sCORE\sSTATE/i
          header = false
          next
        end
        next if header
        header = true if line.chomp.strip.empty?
        next if header
        print "\n" if line.chomp.split[0] == '100'
        parse_energy_state(line)
      end
    end
    self
  end

  def parse_core_electrons
    header = true
    open(@filename, 'r') do |fin|
      fin.each_line do |line|
        if line =~ /INSIDE/
          header = false
          next
        end
        next if header
        header = true if line.chomp.strip.empty?
        next if header
        display_core_electron(line)
      end
    end
    self
  end

  def display_core_electron(line)
    if line =~ /ATOM/
      data = line.chomp.split('ATOM=')[1].strip.split
      data.shift
      @ika += 1
    elsif line =~ /ELCOM/
      data = line.chomp.split('ELCOM,ELCO')[1].strip.split
      @ika += 1
    else
      return
    end
    data.map! { |x| x.tr('D', 'e') }.map!(&:to_f)
    case data.size
    when 2 # nonmagnetic case
      print format(" ia = %2d # of elec = %8.4f inside MT. %8.6e elec is lost.\n", @ika, data[0], data[0].round(0) - data[0])
    when 3 # ferromagnetic case
      print format(" ia = %2d # of elec(upspin)    = %8.4f inside MT. %8.6e elec is lost.\n", @ika, data[0], data[0].round(0) - data[0])
      print format("          # of elec(down-spin) = %8.4f inside MT. %8.6e elec is lost.\n", data[1], data[1].round(1) - data[1])
      print format("          # of elec(total)     = %8.4f inside MT. %8.6e elec is lost.\n", data[2], data[2].round(2) - data[2])
    end
  end

  def parse_energy_state(line)
    elevel, energy = line.chomp.split
    elevel = elevel.to_i
    energy = energy.to_f
    n = elevel / 100
    l = (elevel / 10) % 10
    lname = %w(s p d f)
    print format("  %1d%1s: %15.6f(Ry.)\n", n, lname[l], energy)
  end
end

