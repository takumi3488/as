#!/usr/bin/env ruby

#
# generate sl1.dat
#
require 'fileutils'
class GenLapwInit
  def initialize(configure)
    @configure = configure
    @initobj = ParseInit.new(@configure)
    @nka = @initobj.nka
    @valence = @initobj.valence
    @valence_level = @initobj.valence_level
    @atom_name = @initobj.atom_name
    @no_valence_erange = @configure[:no_valence_erange]
    header
    info_val
    write_sl
    footer
    copy_fort_files
  end


  def header
    print "generate sl1.dat from prp and init.\n"
    print "\n"
  end

  def footer
    print "\n"
    print "end of generation sl1.dat\n"
  end
  #
  # for debug function
  #
  def info_val
    print "total number of valence electron = #{@initobj.total_val}\n"
    print "\n"
    1.upto(@nka) do |ika|
      print "#{@atom_name[ika - 1]} "
      print "num of valence = #{@initobj.nval[ika - 1]}\n"
      0.upto(7) do |l|
        next unless @valence[ika - 1][l]
        print "  l=#{l}  "
        print 'valence'
        if occupy?(ika, l)
          print ' and full-occupied     '
        else
          print ' and partially-occupied'
        end
        print " Energy = #{@valence_level[ika - 1][l]} (Ry.)"
        erange = energy_range(@valence_level[ika - 1][l], ika, l)
        print format(' :range[ %4.1f-> %4.1f ]', erange[0], erange[1])
        print "\n"
      end
    end
  end

  def write_sl
    filename = 'sl1.dat'
    @fout = open(filename, 'w')
    write_magtype
    write_rmt
    write_nband
    1.upto(@nspin) do |ispin|
      1.upto(@nka) { |ika| write_levst(ika, ispin) }
    end
    write_vxc
    write_kpoint
    write_occ
    write_iteration
    @fout.close
  end

  def write_magtype
    if @initobj.magtype == :nonmag
      @fout.print "NONMAGNETIC\n"
      @nspin = 1
    else
      @fout.print "MAGNETIC\n"
      @nspin = 2
    end
  end

  def write_rmt
    const = @configure[:rmt_const]
    @fout.print ' '
    @initobj.rmt.each { |rmt| @fout.print "#{rmt}  " }
    @fout.print '      # radius for muffin-tin for '
    @initobj.atom_name.each_with_index do |aname,ind|
      if @initobj.atom_name.size - 1 ==  ind
        if @initobj.atom_name.size > 1
          @fout.print "and #{aname}"
        else
          @fout.print "#{aname}"
        end
      else
        if @initobj.atom_name.size > 3 && @initobj.atom_name.size - 1 - ind > 3 then
          @fout.print "#{aname}, "
        else
          @fout.print "#{aname} "
        end
      end
    end 
    @fout.print "\n"
    @fout.print "LAPW\n"
    min_rmt = @initobj.rmt.map(&:to_f).min
    @fout.print format(' %5.2f 7', const / min_rmt)
    @fout.print "    # cutoff for spw and ssw\n"
  end

  def write_nband
    info_level = @configure[:info_level]
    band_ratio = @configure[:band_ratio]
    nband = (@initobj.total_val * band_ratio).to_i
    @fout.print " #{nband}         # number of band\n"
    @fout.print " #{info_level}          # infomation level\n"
  end

  def write_levst(ika, ispin)
    @fout.print " #{ispin} #{ika}"
    0.upto(7) { |l| @fout.print "     #{l}" }
    @fout.print " # energy level for logarithm derivative for #{@initobj.atom_name[ika - 1]}"
    if @nspin == 2
      if ispin == 1
        @fout.print "(up-spin)"
      else
        @fout.print "(down-spin)"
      end
    end
    @fout.print "\n"
    0.upto(1) do |i|
      @fout.print '    '
      0.upto(7) do |l|
        erange = set_erange(ika, l)
        @fout.print format('  %4.1f', erange[i])
      end
      @fout.print "\n"
    end
  end

  def write_vxc
    @fout.print "LSDF DATA\n"
    @fout.print "#{@configure[:vxc]}"
    @fout.print "         # exchange-correlation potential\n"
  end

  def write_kpoint
    @fout.print "K POINT\n"
    @initobj.kmesh.each { |k| @fout.print " #{k}" }
    @fout.print "      # k-point mesh(nkx, nky, nkz)\n"
  end

  def write_occ
    #
    # to calculalate fermi energy
    #
    @fout.print " 1          # method to calculate occupation(1: tetrahedron / 0:gaussian )\n"
    @fout.print " 0.10       # smearing parameter for gaussian method.\n"
  end

  def write_iteration
    mixing = @configure[:mixing]
    @fout.print "ITERATION DATA;\n"
    @fout.print " 0              # mixing method(0:simple iteration. 1:simple mixing. 2:broyden method(you need fort.36 and 37).\n"
    @fout.print " 3  #{mixing[0]}  #{mixing[1]}    # info level,  ratio of mixing for new potential and not use.\n"
    @fout.print " 1 -1.0  -1.0   # end of calculation\n"
  end

  def set_erange(ika, l)
    erange =
      if @valence[ika - 1][l]
        energy_range(@valence_level[ika - 1][l], ika, l)
      else
        @no_valence_erange
      end
    erange
  end

  def occupy?(ika, l)
    if @initobj.filled[ika - 1][l]
      true
    else
      false
    end
  end

  def energy_range(energy, ika, l)
    ewidth =
      if l < 2
        0.5
      else
        0.6
      end
    de = 0.5

    case l
    when 0
      if occupy?(ika, l)
        de = -0.2
      else
        de = -0.1
      end
    when 1
      if occupy?(ika, l)
        de = -0.2
      else
        de = -0.1
      end
    when 2
      if occupy?(ika, l)
        de = 0.6
      else
        de = 0.7
      end
    when 3
      if occupy?(ika, l)
        de = 0.6
      else
        de = 0.7
      end
    else
      de = 0.5
    end
    if energy.nil?
      print "ika = #{ika}\n"
      print "l = #{l}\n"
      exit
    end
    [(energy + de).round(1) - ewidth, (energy + de).round(1) + ewidth]
  end

  def copy_fort_files
    path = 'wk'
    FileUtils.mkdir_p(path) unless FileTest.exist?(path)
    prpdir = @configure[:prpdir]
    unless prpdir.nil?
      original = prpdir + '/fort.1'
      newfile = 'wk/fort.1'
      if File.exist?(original)
        print " file copy from #{original} to #{newfile}\n"
        FileUtils.copy(original, newfile) if File.exist?(original)
      end
    end
    [30, 35].each do |i|
      initdir = @configure[:initdir]
      next if initdir.nil?
      original = initdir + '/fort.' + i.to_s
      newfile = 'wk/fort.' + i.to_s 
      if File.exist?(original)
        print " file copy from #{original} to #{newfile}\n"
        FileUtils.copy(original, newfile)
      end
    end
  end
end

#
# parse data
#   from prp and init
#     made by Hiroki Funashima in Kobe
#
class ParseInit
  attr_reader :magtype
  attr_reader :rmt
  attr_reader :nka, :atom_name, :nval
  attr_reader :total_val, :valence, :valence_level, :filled
  attr_reader :kmesh
  def initialize(configure)
    @configure = configure # hash
    init_params
    read_in
    calc_kmesh
  end

  def calc_kmesh
    @kmesh = []
    if @il > 0
      mina = [@lattice_const[0], @lattice_const[1], @lattice_const[2]].min
      0.upto(2) do |i| 
        a =  @lattice_const[i] 
        @kmesh[i] = ((mina / a) + 0.5).to_i * 4 
      end
    else
      mina = [@lattice_const[0], @lattice_const[1], @lattice_const[2]].min
      maxa = [@lattice_const[0], @lattice_const[1], @lattice_const[2]].max
      if maxa / mina > 2.0
        if @lattice_const[0] < @lattice_const[2]
          @kmesh[0] = ( (@lattice_const[2] / @lattice_const[0]) + 0.5 ).to_i * 3
          @kmesh[1] = @kmesh[0]
          @kmesh[2] = 2
        else
          @kmesh[0] = 6
          @kmesh[1] = 6
          @kmesh[2] = ( (@lattice_const[0] / @lattice_const[2]) + 0.5 ).to_i * 2
        end
      else
        @kmesh = [6, 6, 4]
      end
    end
  end

  def init_params
    @orbital = []
    @atom_name = []
    @nelec = {}
    @nval = []
  end

  def calc_total_val
    @total_val = 0
    1.upto(@nka) { |ika| @total_val += @atom_array[ika - 1] * @nval[ika - 1] }
  end

  def read_in
    parse_initin
    parse_prpin
    parse_init_energy_level
    check_n_for_atom
    check_total_electron
    check_valence
    calc_total_val
  end

  def check_valence
    @valence = []
    @filled = []
    @valence_level = []
    1.upto(@nka) do |ika| 
      @valence[ika - 1], @filled[ika - 1], @valence_level[ika - 1] = treat_valence?(ika)
    end
  end

  def treat_valence?(ika)
    @nval[ika - 1] = 0 if @nval[ika - 1].nil?
    maxl = 7
    state = []
    filled = []
    valence_level = []
    filled_electron = [2.0, 6.0, 10.0, 14.0]
    eps = 0.001
    maxl.times { state.push(nil) }
    maxl.times { filled.push(false) }
    @orbital[ika - 1].each_with_index do |orb, ind|
      l = (orb / 10) % 10
      treat_elec = orb % 10
      if treat_elec > 0
        state[l] = true
#        print "orbital ..."
#        p @orbital[ika - 1][ind]
#        print "energy ..."
#        p @energy_level[ika - 1][ind]
        valence_level[l] = @energy_level[ika - 1][ind]
        @nval[ika - 1] += @total_nelec[ika - 1][ind]
        if (filled_electron[l] - @total_nelec[ika - 1][ind]).abs < eps
          filled[l] = true
        end
      else
        state[l] = false
      end
    end
    [state, filled, valence_level]
  end

  def check_total_electron
    @total_nelec = []
    1.upto(@nka) { |ika| count_elec(ika) }
  end

  def count_elec(ika)
    norbital = @orbital[ika - 1].size
    @total_nelec[ika - 1] = [] if @total_nelec[ika].nil?
    0.upto(norbital - 1) do |i|
      total =
        case @magtype
        when :nonmag then @nelec[ika - 1][:up][i].to_f
        when :ferro then @nelec[ika - 1][:up][i].to_f + @nelec[ika - 1][:down][i].to_f
        end
      @total_nelec[ika - 1][i] = total
    end
  end

  def filecheck(filename)
    unless File.exist?(filename)
      print "file:#{filename} is not found.\n"
      exit
    end
  end

  def check_n_for_atom
    atomcount = []
    linecount = 0
    header = true
    filename = @configure[:prpdir] + '/' + @configure[:prpin]
    filecheck(filename)
    open(filename, 'r') do |fin|
      fin.each_line do |linedata|
        unless header
          linecount += 1
          if linecount == 1
            @lattice_const = linedata.chomp.split.map(&:to_f).reject { |x| x.nil? }
          elsif linecount == 5
            atomcount = linedata.chomp.split.map(&:to_i)
          end
        end
        header = false if linedata.chomp.strip.downcase =~ /^latt/
      end
    end
    @atom_array = []
    1.upto(@nka) { |ika| @atom_array[ika - 1] = atomcount[ 2*ika - 1] - atomcount[ 2*ika - 2]  + 1 }
  end

  def parse_init_energy_level
    @energy_level = []
    filename = @configure[:initdir] + '/' + @configure[:initout]
    filecheck(filename)
    header = true
    ika = 1
    open(filename, 'r') do |fin|
      fin.each_line do |line|
        linedata = line.chomp.strip
        unless header 
          if linedata.empty? || linedata =~ /^0TOTAL\sCPU\sTIME/i
            ika += 1
            break if @nka < ika
          else
            @energy_level[ika - 1] = [] if @energy_level[ika - 1].nil?
            iorbital = linedata.split[0].to_i
            energy_level = linedata.split[3].to_f
            @energy_level[ika - 1][iorbital - 1] = energy_level
          end
        end
        if linedata.empty?
          header = true unless header
        elsif linedata =~ /orbit\s+electrons/i
          header = false if header
        end
      end
    end
  end

  def parse_magtype(linedata)
    case linedata.chomp.strip.upcase
    when /^NONM/ then @magtype = :nonmag
    when /^MAGN/ then @magtype = :ferro
    else
      print "unknown magnetic type:#{linedata}"
      exit
    end
  end

  def parse_initin
    filename = @configure[:initdir] + '/' + @configure[:initin]
    filecheck(filename)
    open(filename, 'r') do |fin|
      fin.each_line do |linedata|
        case fin.lineno
        when 1
          parse_magtype(linedata)
        when 2
          @nka, @na = linedata.chomp.split.map(&:to_i)
        when 3
          @rmt = linedata.chomp.split
        when 4
          @vxc = linedata.chomp.split[0]
        end
      end
    end
  end

  def parse_prpin
    header = true
    linecount = 0
    filename = @configure[:prpdir] + '/' + @configure[:prpin]
    filecheck(filename)
    open(filename, 'r') do |fin|
      fin.each_line do |line|
        @il = line.chomp.split[0].to_i if fin.lineno == 3
        unless header
          linecount += 1
          if linecount <= @nka * 5
            linedata = line.chomp.strip
            parse_atominfo(linedata, linecount)
          end
        end
        header = false if line.chomp.strip.downcase =~ /^atomic/
      end
    end
  end

  def parse_atominfo(linedata, linecount)
    atom_index = (linecount - 1) / 5
    case (linecount - 1) % 5
    when 0
      @atom_name[atom_index] = linedata.split[0]
    when 2
      @orbital[atom_index] = linedata.split.map(&:to_i)
    when 3
      @nelec[atom_index] = {}
      @nelec[atom_index][:up] = linedata.split
    when 4
      @nelec[atom_index][:down] = linedata.split
    end
  end
end
