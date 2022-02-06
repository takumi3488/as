#!/usr/bin/env ruby
#
# Configure for Kansai Directory and file name-rule
#    written by Hiroki Funashima, 2017 in Kobe
#
class ParseInit2LapwConfig
  attr_reader :configure
  def initialize(configfile)
    if File.exist?(configfile)
      @configfile = configfile
      parse
    else
      print "file:#{configfile} is not found.\n"
      print "default values are used.\n"
      @configfile = 'prp.dat'
      parse
    end
  end

  def parse
    pwd = ENV['PWD'].split("/")[-1]
    token = '='
    prpdir = def_value('prpdir', token, pwd ).to_s
    prpin = def_value('prpin', token, 'prp.dat').to_s
    prpout = def_value('prpout', token, 'prp.out').to_s
    #
    initdir = def_value('initdir', token, pwd).to_s
    initin = def_value('initin', token, 'init.dat').to_s
    initout = def_value('initout', token, 'init.out').to_s
    #
    info_level = def_value('info_level', token, 1).to_i
    rmt_const = def_value('rmt_const', token, 1.0).to_f
    band_ratio = def_value('band_ratio', token, 0.75).to_f
    emin = def_value('emin', token, 0.0).to_f
    emax = def_value('emax', token, 1.0).to_f
    no_valence_erange = [emin, emax]
    vxc = def_value('vxc', token, 'GL').to_s
    mixing_ratio = def_value('mixing_ratio', token, 0.1).to_f
    @configure = {
      prpdir: prpdir,
      prpin: prpin,
      prpout: prpout,
      initdir: initdir,
      initin: initin,
      initout: initout,
      info_level: info_level,
      rmt_const: rmt_const,
      band_ratio: band_ratio,
      no_valence_erange: no_valence_erange,
      vxc: vxc,
      mixing: [mixing_ratio, mixing_ratio]
    }
  end

  #
  # define value from @configfile
  #
  def def_value(keyword, token, default_value = nil)
    value = nil
    open(@configfile, 'r') do |fin|
      fin.each_line do |linedata|
        if linedata.strip !~ /^#/ && linedata.chomp.strip !~ /^$/
          key, val = line_to_key_and_val(linedata, token)
          value = val if key.upcase == keyword.upcase
        end
      end
    end
    if value.nil?
      check_default_value(default_value)
    else
      value
    end
  end

  def check_default_value(default_value)
    if default_value.nil?
      print "==================== Error(Ising_IO) ====================\n"
      print "#{value} has not been defined in #{@configfile}.\n"
      print "=========================================================\n"
      exit
    else
      default_value unless default_value.to_s.strip.empty?
    end
  end


  def line_to_key_and_val(linedata, token)
    keyword = linedata.chomp.split(token)[0].strip.downcase
    begin
      value   = linedata.chomp.split(token)[1].split('#')[0].split('!')[0].strip
    rescue
      value = nil
    end
    [keyword, value]
  end
end
