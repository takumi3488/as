#!/usr/bin/env ruby
#
# class get inputfile
#
class GetInputFile
  attr_reader :filename, :index
  def initialize(configure)
    @prefix = configure[:prefix]
    @suffix = configure[:suffix]
    @dir = configure[:dir]
    @dir = '..'
    inputfile
  end

  def check_header(header)
    case header
    when /^MAGN/, /^NONM/, /^SPIN/
      scfinput = true
    else
      scfinput = false
    end
    scfinput
  end

  def inputfile
    search = @dir + '/' + @prefix + '*' + @suffix
    @filename = nil
    `ls  #{search}`.split("\n").each do |target|
      header = `head -1 #{target}`.chomp
      @filename = target if check_header(header)
    end
    if @filename.nil?
      print "no inputfile\n"
      exit
    end
    @index = File.basename(@filename, @suffix).split(@prefix)[1]
  end

end
